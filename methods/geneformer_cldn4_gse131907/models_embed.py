"""Cell embeddings for the CLDN4-high vs low contrast.

Geneformer follows ctheodoris/Geneformer tokenization (total-count scale to
10_000, divide by the pretrained nonzero median, rank, truncate) and the
official mean-pool of the second-to-last hidden state.

scGPT follows wanglab/scGPT-human ``embed_data``: within-cell quantile bins
(51), ``<cls>`` at position 0 with pad value -2, post-norm flash-style
transformer, cell vector = L2-normalized final hidden state at ``<cls>``.
Flash-attn is not installed; the same Wqkv weights run with eager attention.
"""

from __future__ import annotations

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F


def tokenize_geneformer(
    counts: np.ndarray,
    ensembl_ids: list[str],
    token_of: dict,
    median_of: dict,
    n_counts: np.ndarray,
    *,
    holdout_ensembl: str | None = None,
    max_len: int = 2048,
    cls_id: int | None = None,
    eos_id: int | None = None,
) -> tuple[list[np.ndarray], dict]:
    """Rank-value encode each cell. counts is cells x genes, raw UMI."""
    keep_j: list[int] = []
    toks: list[int] = []
    meds: list[float] = []
    for j, ens in enumerate(ensembl_ids):
        if holdout_ensembl is not None and ens == holdout_ensembl:
            continue
        tok = token_of.get(ens)
        med = median_of.get(ens)
        if tok is None or med is None:
            continue
        med_f = float(med)
        if not np.isfinite(med_f) or med_f <= 0:
            continue
        keep_j.append(j)
        toks.append(int(tok))
        meds.append(med_f)
    token_arr = np.asarray(toks, dtype=np.int64)
    med_arr = np.asarray(meds, dtype=np.float64)
    sub = counts[:, keep_j]
    special = cls_id is not None and eos_id is not None
    room = max_len - 2 if special else max_len
    seqs: list[np.ndarray] = []
    n_empty = 0
    for i in range(sub.shape[0]):
        row = sub[i]
        nz = np.flatnonzero(row)
        if nz.size == 0 or n_counts[i] <= 0:
            n_empty += 1
            seqs.append(np.zeros(0, dtype=np.int64))
            continue
        scaled = row[nz] / float(n_counts[i]) * 10000.0 / med_arr[nz]
        order = np.argsort(-scaled, kind="mergesort")
        chosen = token_arr[nz[order]][:room]
        if special:
            chosen = np.concatenate(
                [
                    np.asarray([cls_id], dtype=np.int64),
                    chosen,
                    np.asarray([eos_id], dtype=np.int64),
                ]
            )
        seqs.append(chosen)
    info = {
        "n_genes_in_vocab": int(len(keep_j)),
        "n_empty_cells": int(n_empty),
        "holdout": holdout_ensembl,
        "special_tokens": bool(special),
        "max_len": int(max_len),
    }
    return seqs, info


def embed_geneformer(
    model,
    seqs: list[np.ndarray],
    *,
    layer_index: int,
    special_tokens: bool,
    batch_size: int = 8,
    vocab_size: int | None = None,
    collect_tokens: bool = False,
) -> tuple[np.ndarray, np.ndarray | None, np.ndarray | None]:
    """Mean-pool ``hidden_states[layer_index]`` over real gene tokens."""
    hidden = int(model.config.hidden_size)
    embs = np.full((len(seqs), hidden), np.nan, dtype=np.float32)
    gene_sum = None
    gene_n = None
    if collect_tokens:
        if vocab_size is None:
            raise ValueError("vocab_size is required when collect_tokens is set")
        gene_sum = np.zeros((vocab_size, hidden), dtype=np.float64)
        gene_n = np.zeros(vocab_size, dtype=np.int32)
    model.eval()
    with torch.inference_mode():
        for start in range(0, len(seqs), batch_size):
            batch = seqs[start : start + batch_size]
            lens = [int(s.size) for s in batch]
            if max(lens) == 0:
                continue
            maxlen = max(lens)
            ids = torch.zeros(len(batch), maxlen, dtype=torch.long)
            attn = torch.zeros(len(batch), maxlen, dtype=torch.long)
            for i, s in enumerate(batch):
                if s.size == 0:
                    continue
                ids[i, : s.size] = torch.from_numpy(s)
                attn[i, : s.size] = 1
            out = model(
                input_ids=ids,
                attention_mask=attn,
                output_hidden_states=True,
            )
            hs = out.hidden_states[layer_index]
            for i, length in enumerate(lens):
                if length <= 0:
                    continue
                if special_tokens:
                    if length < 3:
                        continue
                    sl = hs[i, 1 : length - 1]
                    tok = batch[i][1:-1]
                else:
                    sl = hs[i, :length]
                    tok = batch[i]
                embs[start + i] = sl.mean(dim=0).float().cpu().numpy()
                if collect_tokens and tok.size:
                    np.add.at(gene_sum, tok, sl.double().cpu().numpy())
                    np.add.at(gene_n, tok, 1)
            if start % (batch_size * 20) == 0:
                print(f"  geneformer {start + len(batch)}/{len(seqs)}", flush=True)
    return embs, gene_sum, gene_n


def _bin_row(row: np.ndarray, n_bins: int, rng: np.random.Generator) -> np.ndarray:
    """scGPT ``binning`` + ``_digitize(side='both')`` with an explicit RNG."""
    out = np.zeros(row.shape[0], dtype=np.float32)
    nz = np.flatnonzero(row)
    if nz.size == 0:
        return out
    vals = row[nz].astype(np.float64)
    bins = np.quantile(vals, np.linspace(0, 1, n_bins - 1))
    left = np.digitize(vals, bins)
    right = np.digitize(vals, bins, right=True)
    u = rng.random(vals.size)
    digits = np.ceil(u * (right - left) + left).astype(np.int64)
    out[nz] = digits.astype(np.float32)
    return out


def tokenize_scgpt(
    counts: np.ndarray,
    gene_ids: np.ndarray,
    *,
    holdout_gene_id: int | None,
    cls_id: int,
    pad_id: int,
    pad_value: float = -2.0,
    max_length: int = 1200,
    n_bins: int = 51,
    seed: int = 1,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Return (gene_id, binned_value) per cell, CLS prepended, length <= max_length."""
    seqs: list[tuple[np.ndarray, np.ndarray]] = []
    room = max_length - 1
    for i in range(counts.shape[0]):
        row = counts[i]
        nz = np.flatnonzero(row)
        if holdout_gene_id is not None and nz.size:
            nz = nz[gene_ids[nz] != holdout_gene_id]
        g = gene_ids[nz]
        v = row[nz].astype(np.float64)
        # Drop genes that are not in the vocabulary (id < 0).
        ok = g >= 0
        g = g[ok]
        v = v[ok]
        rng_np = np.random.default_rng(seed + i)
        binned = _bin_row(v, n_bins, rng_np) if v.size else np.zeros(0, dtype=np.float32)
        if g.size > room:
            gen = torch.Generator()
            gen.manual_seed(seed + i)
            take = torch.randperm(int(g.size), generator=gen)[:room].numpy()
            # Stable order so the permutation is the only randomness.
            take.sort()
            g = g[take]
            binned = binned[take]
        genes = np.concatenate([np.asarray([cls_id], dtype=np.int64), g.astype(np.int64)])
        vals = np.concatenate(
            [np.asarray([pad_value], dtype=np.float32), binned.astype(np.float32)]
        )
        seqs.append((genes, vals))
    return seqs


class _Block(nn.Module):
    def __init__(self, d_model: int, nhead: int, d_hid: int):
        super().__init__()
        self.nhead = nhead
        self.Wqkv = nn.Linear(d_model, 3 * d_model)
        self.out_proj = nn.Linear(d_model, d_model)
        self.linear1 = nn.Linear(d_model, d_hid)
        self.linear2 = nn.Linear(d_hid, d_model)
        self.norm1 = nn.LayerNorm(d_model)
        self.norm2 = nn.LayerNorm(d_model)

    def forward(self, src: torch.Tensor, pad_mask: torch.Tensor) -> torch.Tensor:
        bsz, length, d_model = src.shape
        qkv = self.Wqkv(src)
        q, k, v = qkv.chunk(3, dim=-1)
        head = d_model // self.nhead
        def split(t: torch.Tensor) -> torch.Tensor:
            return t.view(bsz, length, self.nhead, head).transpose(1, 2)

        q, k, v = split(q), split(k), split(v)
        scores = torch.matmul(q, k.transpose(-2, -1)) / (head ** 0.5)
        scores = scores.masked_fill(pad_mask[:, None, None, :], float("-inf"))
        weights = torch.softmax(scores, dim=-1)
        weights = torch.nan_to_num(weights, nan=0.0)
        mixed = torch.matmul(weights, v).transpose(1, 2).contiguous().view(bsz, length, d_model)
        src = self.norm1(src + self.out_proj(mixed))
        src = self.norm2(src + self.linear2(F.relu(self.linear1(src))))
        return src


class ScGPTEncoder(nn.Module):
    """Eager post-norm encoder matching wanglab/scGPT-human weights."""

    def __init__(self, n_token: int, d_model: int, nhead: int, d_hid: int, nlayers: int, pad_id: int):
        super().__init__()
        self.gene_emb = nn.Embedding(n_token, d_model, padding_idx=pad_id)
        self.gene_norm = nn.LayerNorm(d_model)
        self.val_l1 = nn.Linear(1, d_model)
        self.val_l2 = nn.Linear(d_model, d_model)
        self.val_norm = nn.LayerNorm(d_model)
        self.layers = nn.ModuleList(_Block(d_model, nhead, d_hid) for _ in range(nlayers))
        self.max_value = 512

    def encode(self, genes: torch.Tensor, values: torch.Tensor, pad_mask: torch.Tensor) -> torch.Tensor:
        x = self.gene_norm(self.gene_emb(genes))
        v = values.unsqueeze(-1).clamp(max=self.max_value)
        v = self.val_norm(self.val_l2(F.relu(self.val_l1(v))))
        h = x + v
        for layer in self.layers:
            h = layer(h, pad_mask)
        return h


def load_scgpt(weight_path: str, pad_id: int, n_token: int = 60697) -> ScGPTEncoder:
    sd = torch.load(weight_path, map_location="cpu", weights_only=False)
    model = ScGPTEncoder(n_token, 512, 8, 512, 12, pad_id)
    mapped = {
        "gene_emb.weight": sd["encoder.embedding.weight"],
        "gene_norm.weight": sd["encoder.enc_norm.weight"],
        "gene_norm.bias": sd["encoder.enc_norm.bias"],
        "val_l1.weight": sd["value_encoder.linear1.weight"],
        "val_l1.bias": sd["value_encoder.linear1.bias"],
        "val_l2.weight": sd["value_encoder.linear2.weight"],
        "val_l2.bias": sd["value_encoder.linear2.bias"],
        "val_norm.weight": sd["value_encoder.norm.weight"],
        "val_norm.bias": sd["value_encoder.norm.bias"],
    }
    for i in range(12):
        p = f"transformer_encoder.layers.{i}"
        mapped[f"layers.{i}.Wqkv.weight"] = sd[f"{p}.self_attn.Wqkv.weight"]
        mapped[f"layers.{i}.Wqkv.bias"] = sd[f"{p}.self_attn.Wqkv.bias"]
        mapped[f"layers.{i}.out_proj.weight"] = sd[f"{p}.self_attn.out_proj.weight"]
        mapped[f"layers.{i}.out_proj.bias"] = sd[f"{p}.self_attn.out_proj.bias"]
        mapped[f"layers.{i}.linear1.weight"] = sd[f"{p}.linear1.weight"]
        mapped[f"layers.{i}.linear1.bias"] = sd[f"{p}.linear1.bias"]
        mapped[f"layers.{i}.linear2.weight"] = sd[f"{p}.linear2.weight"]
        mapped[f"layers.{i}.linear2.bias"] = sd[f"{p}.linear2.bias"]
        mapped[f"layers.{i}.norm1.weight"] = sd[f"{p}.norm1.weight"]
        mapped[f"layers.{i}.norm1.bias"] = sd[f"{p}.norm1.bias"]
        mapped[f"layers.{i}.norm2.weight"] = sd[f"{p}.norm2.weight"]
        mapped[f"layers.{i}.norm2.bias"] = sd[f"{p}.norm2.bias"]
    missing, unexpected = model.load_state_dict(mapped, strict=False)
    if missing:
        raise RuntimeError(f"scGPT load missing keys: {missing}")
    model.eval()
    return model


def embed_scgpt(
    model: ScGPTEncoder,
    seqs: list[tuple[np.ndarray, np.ndarray]],
    pad_id: int,
    *,
    batch_size: int = 4,
) -> np.ndarray:
    """L2-normalized CLS embedding, matching ``get_batch_cell_embeddings``."""
    embs = np.full((len(seqs), 512), np.nan, dtype=np.float32)
    model.eval()
    with torch.inference_mode():
        for start in range(0, len(seqs), batch_size):
            batch = seqs[start : start + batch_size]
            lens = [int(g.size) for g, _ in batch]
            maxlen = max(lens) if lens else 1
            genes = torch.full((len(batch), maxlen), pad_id, dtype=torch.long)
            values = torch.full((len(batch), maxlen), -2.0, dtype=torch.float32)
            for i, (g, v) in enumerate(batch):
                genes[i, : g.size] = torch.from_numpy(g)
                values[i, : v.size] = torch.from_numpy(v)
            pad_mask = genes.eq(pad_id)
            h = model.encode(genes, values, pad_mask)
            cls = h[:, 0, :]
            cls = cls / cls.norm(dim=1, keepdim=True).clamp(min=1e-8)
            embs[start : start + len(batch)] = cls.float().cpu().numpy()
            if start % (batch_size * 20) == 0:
                print(f"  scgpt {start + len(batch)}/{len(seqs)}", flush=True)
    return embs
