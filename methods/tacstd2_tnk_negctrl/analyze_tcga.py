#!/usr/bin/env python3
"""TCGA: TACSTD2 vs a T/NK score, partialled on CLDN4 or on negative controls.

Primary tumors only (sample type 01). UCSC Xena GDC STAR log2(TPM+1).
Replicate aliquots of the same patient are averaged.

Primary cohorts are the locked keratin-funnel set:
LUAD, BRCA, CESC, KIRC, STAD, BLCA, PAAD.
LUSC is reported and is included only in a labeled sensitivity meta.

Primary outcome is the mean of CD3D, CD3E, CD3G, CD8A, NKG7, GNLY, and KLRD1
(the concordant-4 T/NK markers plus CD3G). CD3-only and CD8A are sensitivities.
No purity term and no extra keratin term: KRT19 is itself a negative control.
"""

from __future__ import annotations

import gzip
import hashlib
import os
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from stats import bh_fdr, dl_meta, partial_spearman, spearman, spearman_p

ROOT = Path(__file__).resolve().parent
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"
CACHE = Path(os.environ.get("TCGA_NEGCTRL_CACHE", "/tmp/tcga_negctrl"))
GDC = "https://gdc-hub.s3.us-east-1.amazonaws.com/download"
PROBEMAP_URL = f"{GDC}/gencode.v36.annotation.gtf.gene.probemap"

FUNNEL = ["LUAD", "BRCA", "CESC", "KIRC", "STAD", "BLCA", "PAAD"]
EXTRA = ["LUSC"]
COHORTS = FUNNEL + EXTRA
GENES = [
    "TACSTD2", "CLDN3", "CLDN4", "CLDN7", "EPCAM", "MUC1", "KRT19",
    "CD3D", "CD3E", "CD3G", "CD8A", "NKG7", "GNLY", "KLRD1",
]
CONDITIONERS = ["CLDN4", "CLDN3", "CLDN7", "EPCAM", "MUC1", "KRT19"]
CONTROLS = ["CLDN3", "CLDN7", "EPCAM", "MUC1", "KRT19"]
TNK_GENES = ["CD3D", "CD3E", "CD3G", "CD8A", "NKG7", "GNLY", "KLRD1"]
CD3_GENES = ["CD3D", "CD3E", "CD3G"]
N_PERM = 10000
SEED = 1


def say(msg: str) -> None:
    print(msg, flush=True)


def download(url: str, dest: Path) -> None:
    if dest.exists() and dest.stat().st_size > 0:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": "tcga-tacstd2-negctrl/1.0"})
    with urllib.request.urlopen(req, timeout=180) as resp, tmp.open("wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(dest)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_probemap(symbols: list[str]) -> dict[str, str]:
    path = CACHE / "gencode.v36.annotation.gtf.gene.probemap"
    download(PROBEMAP_URL, path)
    pm = pd.read_csv(path, sep="\t")
    mapping = {}
    for symbol in symbols:
        sub = pm[pm["gene"] == symbol].copy()
        if sub.empty:
            raise SystemExit(f"probemap missing {symbol}")
        sub = sub[sub["chrom"].astype(str).str.match(r"^chr([0-9]+|X|Y)$")]
        if sub.empty:
            raise SystemExit(f"no primary-chrom id for {symbol}")
        sub["span"] = sub["chromEnd"] - sub["chromStart"]
        sub = sub.sort_values("span", ascending=False)
        mapping[symbol] = str(sub.iloc[0]["id"])
    if len(set(mapping.values())) != len(mapping):
        raise SystemExit(f"Ensembl id collision: {mapping}")
    return mapping


def patient_id(barcode: str) -> str | None:
    parts = barcode.replace(".", "-").split("-")
    if len(parts) < 4 or not parts[3].startswith("01"):
        return None
    return "-".join(parts[:3])


def extract_cohort(cohort: str, id_to_gene: dict[str, str]) -> pd.DataFrame:
    cache_path = CACHE / f"{cohort}.primary01.tsv.gz"
    if cache_path.exists() and cache_path.stat().st_size > 0:
        say(f"[cache] {cohort}")
        return pd.read_csv(cache_path, sep="\t")
    url = f"{GDC}/TCGA-{cohort}.star_tpm.tsv.gz"
    want = set(id_to_gene)
    say(f"[get] {cohort} {url}")
    req = urllib.request.Request(url, headers={"User-Agent": "tcga-tacstd2-negctrl/1.0"})
    collected: dict[str, np.ndarray] = {}
    with urllib.request.urlopen(req, timeout=300) as resp:
        with gzip.GzipFile(fileobj=resp) as gz:
            header = gz.readline().decode("utf-8").rstrip("\n").split("\t")
            samples = header[1:]
            keep_idx = []
            patients = []
            for i, sample in enumerate(samples):
                pid = patient_id(sample)
                if pid is None:
                    continue
                keep_idx.append(i)
                patients.append(pid)
            if len(keep_idx) < 15:
                raise RuntimeError(f"{cohort}: only {len(keep_idx)} primary-tumor columns")
            keep_idx_arr = np.asarray(keep_idx, dtype=int)
            n_seen = 0
            for raw in gz:
                line = raw.decode("utf-8")
                gid, rest = line.split("\t", 1)
                n_seen += 1
                if n_seen % 10000 == 0:
                    say(f"  {cohort} scanned {n_seen} genes, kept {len(collected)}")
                if gid not in want:
                    continue
                parts = rest.rstrip("\n").split("\t")
                values = np.array([float(x) if x else np.nan for x in parts], dtype=float)
                if values.size != len(samples):
                    raise RuntimeError(f"{cohort} {gid} width {values.size} != {len(samples)}")
                collected[id_to_gene[gid]] = values[keep_idx_arr]
                say(f"  hit {id_to_gene[gid]}")
                if len(collected) == len(want):
                    break
    missing = [sym for sym in set(id_to_gene.values()) if sym not in collected]
    if missing:
        raise RuntimeError(f"{cohort} missing genes: {missing}")
    frame = pd.DataFrame(collected)
    frame.insert(0, "patient", patients)
    collapsed = frame.groupby("patient", as_index=False).mean(numeric_only=True)
    collapsed.insert(0, "cohort", cohort)
    collapsed.to_csv(cache_path, sep="\t", index=False, compression="gzip")
    say(f"  {cohort} patients {collapsed.shape[0]} sha {sha256_file(cache_path)[:12]}")
    return collapsed


def add_scores(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["tnk"] = out[TNK_GENES].mean(axis=1)
    out["cd3"] = out[CD3_GENES].mean(axis=1)
    out["cd8a"] = out["CD8A"]
    return out


def cohort_rho(x, y, z) -> tuple[float, float, int]:
    mask = np.isfinite(x) & np.isfinite(y)
    if z is not None:
        mask &= np.isfinite(z)
    n = int(mask.sum())
    if z is None:
        rho = spearman(x[mask], y[mask])
        return rho, spearman_p(rho, n, 0), n
    rho = partial_spearman(x[mask], y[mask], [z[mask]])
    return rho, spearman_p(rho, n, 1), n


def label_swap_p(frames, outcome, control, rng, n_perm: int) -> dict:
    def meta_from(use_swap: bool, swaps: list[np.ndarray] | None):
        ra, na, rb, nb = [], [], [], []
        for i, frame in enumerate(frames):
            x = frame["TACSTD2"].to_numpy(dtype=float)
            y = frame[outcome].to_numpy(dtype=float)
            a = frame["CLDN4"].to_numpy(dtype=float)
            b = frame[control].to_numpy(dtype=float)
            if use_swap:
                sw = swaps[i]
                aa = np.where(sw, b, a)
                bb = np.where(sw, a, b)
            else:
                aa, bb = a, b
            rho_a, _, n = cohort_rho(x, y, aa)
            rho_b, _, _ = cohort_rho(x, y, bb)
            ra.append(rho_a)
            rb.append(rho_b)
            na.append(n)
            nb.append(n)
        ma = dl_meta(ra, na, 1)
        mb = dl_meta(rb, nb, 1)
        return ma["rho"] - mb["rho"], ma["rho"], mb["rho"]

    obs, rho_a, rho_b = meta_from(False, None)
    count = 0
    for _ in range(n_perm):
        swaps = [rng.random(len(frame)) < 0.5 for frame in frames]
        delta, _, _ = meta_from(True, swaps)
        if np.isfinite(delta) and abs(delta) >= abs(obs) - 1e-15:
            count += 1
    return {
        "delta": obs,
        "partial_cldn4": rho_a,
        "partial_control": rho_b,
        "p": (count + 1) / (n_perm + 1),
    }


def shuffle_p(frames, outcome, gene, rng, n_perm: int) -> float:
    def meta_partial(permuted: bool):
        rhos, ns = [], []
        r0, n0 = [], []
        for frame in frames:
            x = frame["TACSTD2"].to_numpy(dtype=float)
            y = frame[outcome].to_numpy(dtype=float)
            z = frame[gene].to_numpy(dtype=float)
            rho0, _, n = cohort_rho(x, y, None)
            if permuted:
                z = z.copy()
                rng.shuffle(z)
            rho, _, n = cohort_rho(x, y, z)
            rhos.append(rho)
            ns.append(n)
            r0.append(rho0)
            n0.append(n)
        return dl_meta(rhos, ns, 1)["rho"] - dl_meta(r0, n0, 0)["rho"]

    obs = meta_partial(False)
    count = 0
    for _ in range(n_perm):
        att = meta_partial(True)
        if np.isfinite(att) and abs(att) >= abs(obs) - 1e-15:
            count += 1
    return (count + 1) / (n_perm + 1)


def run_block(frames: list[pd.DataFrame], block: str, outcome: str, rng, do_perm: bool) -> tuple[list[dict], list[dict], list[dict]]:
    say(f"block {block} outcome {outcome} n_cohorts={len(frames)}")
    cohort_rows = []
    meta_rows = []
    # unadjusted
    rhos, ns, ps = [], [], []
    for frame in frames:
        rho, p, n = cohort_rho(frame["TACSTD2"].to_numpy(float), frame[outcome].to_numpy(float), None)
        rhos.append(rho)
        ns.append(n)
        ps.append(p)
        cohort_rows.append(
            {
                "block": block,
                "outcome": outcome,
                "cohort": frame["cohort"].iloc[0],
                "conditioner": "none",
                "n": n,
                "rho_unadj": rho,
                "p_unadj": p,
                "rho_partial": rho,
                "p_partial": p,
                "attenuation": 0.0,
                "percent_attenuated": 0.0,
            }
        )
    meta_u = dl_meta(rhos, ns, 0)
    meta_rows.append(
        {
            "block": block,
            "outcome": outcome,
            "conditioner": "none",
            "rho_unadj": meta_u["rho"],
            "p_unadj": meta_u["p"],
            "rho_partial": meta_u["rho"],
            "p_partial": meta_u["p"],
            "I2": meta_u["I2"],
            "ci_lo": meta_u["ci_lo"],
            "ci_hi": meta_u["ci_hi"],
            "attenuation": 0.0,
            "percent_attenuated": 0.0,
            "perm_p_attenuation": float("nan"),
            "k": meta_u["k"],
            "N": meta_u["N"],
        }
    )
    unadj_by = {r["cohort"]: r["rho_unadj"] for r in cohort_rows}
    for gene in CONDITIONERS:
        say(f"  {gene}")
        rhos, ns = [], []
        for frame in frames:
            x = frame["TACSTD2"].to_numpy(float)
            y = frame[outcome].to_numpy(float)
            z = frame[gene].to_numpy(float)
            rho0, p0, n = cohort_rho(x, y, None)
            rho, p, n = cohort_rho(x, y, z)
            rhos.append(rho)
            ns.append(n)
            cohort_rows.append(
                {
                    "block": block,
                    "outcome": outcome,
                    "cohort": frame["cohort"].iloc[0],
                    "conditioner": gene,
                    "n": n,
                    "rho_unadj": rho0,
                    "p_unadj": p0,
                    "rho_partial": rho,
                    "p_partial": p,
                    "attenuation": rho - rho0,
                    "percent_attenuated": float("nan") if abs(rho0) < 1e-8 else 100.0 * (rho0 - rho) / rho0,
                }
            )
        meta_p = dl_meta(rhos, ns, 1)
        att = meta_p["rho"] - meta_u["rho"]
        pct = float("nan") if abs(meta_u["rho"]) < 1e-8 else 100.0 * (meta_u["rho"] - meta_p["rho"]) / meta_u["rho"]
        # Permute only the primary outcome. Sensitivities keep the point estimate.
        if do_perm:
            pp = shuffle_p(frames, outcome, gene, rng, N_PERM)
        else:
            pp = float("nan")
        meta_rows.append(
            {
                "block": block,
                "outcome": outcome,
                "conditioner": gene,
                "rho_unadj": meta_u["rho"],
                "p_unadj": meta_u["p"],
                "rho_partial": meta_p["rho"],
                "p_partial": meta_p["p"],
                "I2": meta_p["I2"],
                "ci_lo": meta_p["ci_lo"],
                "ci_hi": meta_p["ci_hi"],
                "attenuation": att,
                "percent_attenuated": pct,
                "perm_p_attenuation": pp,
                "k": meta_p["k"],
                "N": meta_p["N"],
            }
        )
    head = []
    if do_perm:
        for control in CONTROLS:
            say(f"  swap CLDN4 vs {control}")
            swap = label_swap_p(frames, outcome, control, rng, N_PERM)
            head.append(
                {
                    "block": block,
                    "outcome": outcome,
                    "index": "CLDN4",
                    "control": control,
                    "meta_partial_CLDN4": swap["partial_cldn4"],
                    "meta_partial_control": swap["partial_control"],
                    "delta_partial_CLDN4_minus_control": swap["delta"],
                    "perm_p_two_sided": swap["p"],
                    "n_perm": N_PERM,
                    "N": int(sum(len(f) for f in frames)),
                }
            )
        q = bh_fdr([r["perm_p_two_sided"] for r in head])
        for row, qq in zip(head, q):
            row["q_bh"] = qq
    _ = unadj_by
    return cohort_rows, meta_rows, head


def collinearity(frames: list[pd.DataFrame], block: str) -> list[dict]:
    rows = []
    for gene in CONDITIONERS:
        rhos, ns = [], []
        for frame in frames:
            x = frame["TACSTD2"].to_numpy(float)
            z = frame[gene].to_numpy(float)
            rho, p, n = cohort_rho(x, z, None)
            rhos.append(rho)
            ns.append(n)
            rows.append(
                {
                    "block": block,
                    "level": "cohort",
                    "cohort": frame["cohort"].iloc[0],
                    "gene": gene,
                    "rho_with_TACSTD2": rho,
                    "p": p,
                    "n": n,
                }
            )
        meta = dl_meta(rhos, ns, 0)
        rows.append(
            {
                "block": block,
                "level": "meta",
                "cohort": "DL",
                "gene": gene,
                "rho_with_TACSTD2": meta["rho"],
                "p": meta["p"],
                "n": meta["N"],
                "I2": meta["I2"],
                "ci_lo": meta["ci_lo"],
                "ci_hi": meta["ci_hi"],
            }
        )
    return rows


def plot(meta_rows: list[dict], cohort_rows: list[dict]) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    sub = [r for r in meta_rows if r["block"] == "funnel7" and r["outcome"] == "tnk"]
    by = {r["conditioner"]: r for r in sub}
    order = ["none"] + CONDITIONERS
    labels = ["unadjusted"] + CONDITIONERS
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.4))
    ax = axes[0]
    ypos = np.arange(len(order))[::-1]
    rhos = [by[g]["rho_partial"] for g in order]
    lo = [by[g]["rho_partial"] - by[g]["ci_lo"] for g in order]
    hi = [by[g]["ci_hi"] - by[g]["rho_partial"] for g in order]
    colors = ["#4d4d4d"] + ["#b2182b" if g == "CLDN4" else "#2166ac" for g in CONDITIONERS]
    ax.errorbar(rhos, ypos, xerr=[lo, hi], fmt="none", ecolor="#666666", elinewidth=1, capsize=2)
    ax.scatter(rhos, ypos, c=colors, s=36, zorder=3)
    ax.axvline(0, color="#888888", lw=0.6)
    ax.set_yticks(ypos)
    ax.set_yticklabels(labels)
    ax.set_xlabel("DL Spearman ρ  TACSTD2 vs T/NK score")
    ax.set_title("TCGA funnel-7 partial ρ")
    ax = axes[1]
    genes = CONDITIONERS
    ypos = np.arange(len(genes))[::-1]
    att = [by[g]["attenuation"] for g in genes]
    colors = ["#b2182b" if g == "CLDN4" else "#2166ac" for g in genes]
    ax.axvline(0, color="#888888", lw=0.6)
    ax.scatter(att, ypos, c=colors, s=36, zorder=3)
    for i, g in enumerate(genes):
        xs = [
            r["attenuation"]
            for r in cohort_rows
            if r["block"] == "funnel7" and r["outcome"] == "tnk" and r["conditioner"] == g
        ]
        ax.scatter(xs, np.full(len(xs), ypos[i]), s=12, c="#999999", zorder=2)
    ax.set_yticks(ypos)
    ax.set_yticklabels(genes)
    ax.set_xlabel("Change in ρ  (partial − unadjusted)")
    ax.set_title("TCGA funnel-7 shift")
    fig.tight_layout()
    fig.savefig(FIG / "tcga_attenuation.png", dpi=160)
    fig.savefig(FIG / "tcga_attenuation.pdf")
    plt.close(fig)


def main() -> None:
    CACHE.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    symbol_to_id = load_probemap(GENES)
    id_to_gene = {v: k for k, v in symbol_to_id.items()}
    pd.DataFrame(
        [{"symbol": s, "ensembl": symbol_to_id[s]} for s in GENES]
    ).to_csv(TAB / "tcga_probemap.tsv", sep="\t", index=False)
    frames = []
    for cohort in COHORTS:
        frame = add_scores(extract_cohort(cohort, id_to_gene))
        frames.append(frame)
        say(f"{cohort} n={len(frame)}")
    patient = pd.concat(frames, ignore_index=True)
    keep = ["cohort", "patient", "TACSTD2", *CONDITIONERS, *TNK_GENES, "tnk", "cd3", "cd8a"]
    # CONDITIONERS includes genes also listed; dict-order unique
    cols = []
    for c in keep:
        if c not in cols and c in patient.columns:
            cols.append(c)
    patient[cols].to_csv(TAB / "tcga_patient_scores.tsv.gz", sep="\t", index=False, compression="gzip")
    by_cohort = {f["cohort"].iloc[0]: f for f in frames}
    rng = np.random.default_rng(SEED)
    blocks = {
        "funnel7": [by_cohort[c] for c in FUNNEL],
        "funnel7_plus_lusc": [by_cohort[c] for c in COHORTS],
        "LUAD": [by_cohort["LUAD"]],
        "LUSC": [by_cohort["LUSC"]],
    }
    all_cohort, all_meta, all_head, all_coli = [], [], [], []
    for block, group in blocks.items():
        # Single-cohort blocks still get DL-of-one via dl_meta's n=1 path.
        for outcome in ("tnk", "cd3", "cd8a"):
            do_perm = block == "funnel7" and outcome == "tnk"
            c_rows, m_rows, h_rows = run_block(group, block, outcome, rng, do_perm)
            all_cohort.extend(c_rows)
            all_meta.extend(m_rows)
            all_head.extend(h_rows)
        all_coli.extend(collinearity(group if block == "funnel7" else group, block))
    pd.DataFrame(all_cohort).to_csv(TAB / "tcga_cohort_attenuation.tsv", sep="\t", index=False)
    pd.DataFrame(all_meta).to_csv(TAB / "tcga_meta_attenuation.tsv", sep="\t", index=False)
    pd.DataFrame(all_head).to_csv(TAB / "tcga_head_to_head.tsv", sep="\t", index=False)
    pd.DataFrame(all_coli).to_csv(TAB / "tcga_collinearity.tsv", sep="\t", index=False)
    plot(all_meta, all_cohort)
    say("TCGA done")


if __name__ == "__main__":
    main()
