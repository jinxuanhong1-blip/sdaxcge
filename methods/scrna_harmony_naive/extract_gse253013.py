#!/usr/bin/env python3
"""Stream-extract a gene panel + cell metadata from the GSE253013 RDS.

The GEO object is a double-gzipped XDR RDS (~256k cells). A full Seurat/CDS
load does not fit in 16 GB RAM, so this walker skips large INT/REAL payloads
to disk and keeps only requested gene rows plus character metadata.
"""

from __future__ import annotations

import argparse
import gzip
import json
import struct
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

from config import EXTRACTED, FILES, GENE_PANEL

TYPE_NAMES = {
    0: "NIL",
    1: "SYM",
    2: "LIST",
    3: "CLO",
    4: "ENV",
    5: "PROM",
    6: "LANG",
    7: "SPECIAL",
    8: "BUILTIN",
    9: "CHAR",
    10: "LGL",
    13: "INT",
    14: "REAL",
    15: "CPLX",
    16: "STR",
    17: "DOT",
    18: "ANY",
    19: "VEC",
    20: "EXPR",
    21: "BCODE",
    22: "EXTPTR",
    23: "WEAKREF",
    24: "RAW",
    25: "S4",
    238: "ALTREP",
    239: "ATTRLIST",
    240: "ATTRLANG",
    241: "BASEENV",
    242: "EMPTYENV",
    243: "BCREPREF",
    244: "BCREPDEF",
    249: "NAMESPACE",
    251: "MISSINGARG",
    253: "GLOBALENV",
    254: "NILVALUE",
    255: "REF",
}

BYTECODE_SPECIAL = {21, 243, 244, 6, 2, 240, 239}
LARGE = 50_000


def bits(data: int, start: int, stop: int) -> int:
    mask = ((1 << (stop - start)) - 1) << start
    return (data & mask) >> start


class Info:
    __slots__ = ("typ", "obj", "attr", "tag", "gp", "ref")

    def __init__(self, info_int: int):
        typ = bits(info_int, 0, 8)
        self.typ = typ
        if typ in (254, 255):
            self.obj = False
            self.attr = False
            self.tag = False
            self.gp = 0
        else:
            self.obj = bool(bits(info_int, 8, 9))
            self.attr = bool(bits(info_int, 9, 10))
            self.tag = bool(bits(info_int, 10, 11))
            self.gp = bits(info_int, 12, 28)
        self.ref = bits(info_int, 8, 32) if typ == 255 else 0


class Stream:
    def __init__(self, fh):
        self.fh = fh
        self.n = 0

    def read(self, n: int) -> bytes:
        buf = b""
        while len(buf) < n:
            chunk = self.fh.read(n - len(buf))
            if not chunk:
                raise EOFError(f"EOF at {self.n}, wanted {n} more")
            buf += chunk
        self.n += n
        return buf

    def skip(self, n: int) -> None:
        left = n
        while left:
            chunk = self.fh.read(min(left, 8 * 1024 * 1024))
            if not chunk:
                raise EOFError(f"EOF while skipping {n} at {self.n}")
            left -= len(chunk)
            self.n += len(chunk)

    def i32(self) -> int:
        return struct.unpack(">i", self.read(4))[0]


class Extractor:
    def __init__(self, stream: Stream, tmp: Path, genes: list[str], log):
        self.s = stream
        self.tmp = tmp
        self.want = set(genes)
        self.log = log
        self.refs: list = [None]
        self.arrays: list[dict] = []
        self.str_vecs: list[dict] = []
        self.small_int: list[dict] = []
        self.df_cols: list[dict] = []
        self.classes: list[str] = []
        self.gene_hits: dict[str, list[str]] = defaultdict(list)
        self.depth = 0
        self.path: list[str] = []
        self.n_obj = 0

    def emit(self, msg: str) -> None:
        self.log.write(msg + "\n")
        self.log.flush()
        print(msg, flush=True)

    def parse_info(self, info_int: int | None = None) -> Info:
        if info_int is None:
            info_int = self.s.i32()
        return Info(info_int)

    def parse_obj(self, info_int: int | None = None):
        self.n_obj += 1
        info = self.parse_info(info_int)
        tag = None
        attributes = None
        value = None
        tag_read = False
        attr_read = False
        add_ref = False
        typ = info.typ

        if typ == 0:  # NIL
            value = None
        elif typ == 1:  # SYM
            value = self.parse_obj()
            add_ref = True
        elif typ in (2, 3, 5, 6, 17, 240):  # LIST/CLO/PROM/LANG/DOT/ATTRLANG
            if typ == 240:
                info.attr = True
                typ = 6
            if info.attr:
                attributes = self.parse_obj()
                attr_read = True
            if info.tag:
                tag = self.parse_obj()
                tag_read = True
                tname = self._tag_name(tag)
                if tname:
                    self.path.append(tname)
            car = self.parse_obj()
            cdr = self.parse_obj()
            if info.tag and tag is not None and self.path and self._tag_name(tag) == self.path[-1]:
                self.path.pop()
            # Drop pairlist nodes after walking so the Seurat tree is not retained.
            value = None
        elif typ == 4:  # ENV
            idx = len(self.refs)
            self.refs.append({"type": "ENV"})
            _locked = bool(self.s.i32())
            self.parse_obj()
            self.parse_obj()
            self.parse_obj()
            attributes = self.parse_obj()
            attr_read = True
            value = None
            self.refs[idx] = {"type": "ENV"}
            result = {"type": "ENV", "value": None, "attr": None, "tag": tag}
            if add_ref:
                pass
            return result
        elif typ in (7, 8):  # SPECIAL/BUILTIN
            length = self.s.i32()
            value = self.s.read(length) if length > 0 else b""
        elif typ == 9:  # CHAR
            length = self.s.i32()
            if length >= 0:
                value = self.s.read(length)
            elif length == -1:
                value = None
            else:
                raise NotImplementedError(f"CHAR length {length}")
        elif typ in (10, 13):  # LGL / INT
            value = self._int_array(signed=True)
        elif typ == 14:
            value = self._real_array()
        elif typ == 15:
            length = self.s.i32()
            if length > LARGE:
                self.s.skip(length * 16)
                value = {"skipped_cplx": length}
            else:
                value = np.frombuffer(self.s.read(length * 16), dtype=">c16").copy()
        elif typ in (16, 19, 20):  # STR / VEC / EXPR
            value = self._vector(typ)
        elif typ == 21:  # BCODE
            value = self._bytecode()
            tag_read = True
        elif typ == 22:  # EXTPTR
            idx = len(self.refs)
            self.refs.append({"_extptr": True})
            protected = self.parse_obj()
            ext_tag = self.parse_obj()
            value = (protected, ext_tag)
            self.refs[idx] = value
        elif typ == 23:  # WEAKREF
            idx = len(self.refs)
            self.refs.append({"_weak": True})
            self.parse_obj()
            self.parse_obj()
            self.parse_obj()
            value = None
            self.refs[idx] = value
        elif typ == 24:  # RAW
            length = self.s.i32()
            if length > LARGE:
                self.s.skip(length)
                value = {"skipped_raw": length}
            else:
                value = self.s.read(length)
        elif typ == 25:  # S4
            value = None
        elif typ == 238:  # ALTREP
            altrep_info = self.parse_obj()
            altrep_state = self.parse_obj()
            altrep_attr = self.parse_obj()
            value = self._expand_altrep(altrep_info, altrep_state)
            if isinstance(altrep_attr, dict) and altrep_attr.get("type") not in (
                "NILVALUE",
                "NIL",
            ):
                attributes = altrep_attr
                attr_read = True
                info.attr = True
        elif typ in (241, 242, 251, 253, 254):
            value = None
        elif typ == 243:  # BCREPREF
            _ = self.s.i32()
            value = None
        elif typ == 244:  # BCREPDEF
            _ = self.s.i32()
            inner = self.s.i32()
            return self.parse_obj(inner)
        elif typ == 249:  # NAMESPACE
            assert self.s.i32() == 0
            value = self._vector(19)
            add_ref = True
        elif typ == 255:  # REF
            ref = info.ref
            value = {"ref": ref}
        else:
            raise NotImplementedError(f"type {typ} ({TYPE_NAMES.get(typ)}) at {self.s.n}")

        if info.attr and not attr_read:
            attributes = self.parse_obj()

        result = {
            "type": TYPE_NAMES.get(typ, str(typ)),
            "value": value,
            "attr": attributes,
            "tag": tag,
        }
        if add_ref:
            self.refs.append(result)
        self._ingest(result)
        return result

    def _tag_name(self, tag) -> str | None:
        if not isinstance(tag, dict):
            return None
        if tag.get("type") == "SYM":
            ch = tag.get("value")
            if isinstance(ch, dict) and isinstance(ch.get("value"), (bytes, bytearray)):
                return ch["value"].decode("utf-8", "replace")
        if tag.get("type") == "REF":
            ref = tag.get("value", {}).get("ref")
            if ref and ref < len(self.refs):
                return self._tag_name(self.refs[ref])
        return None

    def _class_names(self, attr) -> list[str]:
        names: list[str] = []
        node = attr
        seen = 0
        while isinstance(node, dict) and node.get("type") in ("LIST", "LANG") and seen < 20:
            seen += 1
            tname = self._tag_name(node.get("tag"))
            car = node.get("value")
            if isinstance(car, tuple) and len(car) == 2:
                val, cdr = car
                if tname == "class":
                    names.extend(self._as_str_list(val))
                node = cdr
            else:
                break
        if isinstance(attr, dict) and attr.get("type") == "STR":
            names.extend(self._as_str_list(attr))
        return names

    def _as_str_list(self, obj) -> list[str]:
        if not isinstance(obj, dict):
            return []
        if obj.get("type") == "STR" and isinstance(obj.get("value"), list):
            out = []
            for el in obj["value"]:
                if isinstance(el, dict) and isinstance(el.get("value"), (bytes, bytearray)):
                    out.append(el["value"].decode("utf-8", "replace"))
                elif isinstance(el, (bytes, bytearray)):
                    out.append(el.decode("utf-8", "replace"))
                elif isinstance(el, str):
                    out.append(el)
            return out
        if obj.get("type") == "CHAR" and isinstance(obj.get("value"), (bytes, bytearray)):
            return [obj["value"].decode("utf-8", "replace")]
        return []

    def _ingest(self, result: dict) -> None:
        classes = self._class_names(result.get("attr"))
        if classes:
            key = "/".join(classes)
            if key not in self.classes:
                self.classes.append(key)
                self.emit(f"class {key} path={'/'.join(self.path)}")

    def _vector(self, typ: int):
        length = self.s.i32()
        if typ == 16:  # STR
            keep = length <= 400_000
            hits = []
            decoded: list[str] = []
            for i in range(length):
                el = self.parse_obj()
                s = ""
                if isinstance(el, dict) and isinstance(el.get("value"), (bytes, bytearray)):
                    s = el["value"].decode("utf-8", "replace")
                if keep:
                    decoded.append(s)
                if s in self.want:
                    hits.append((s, i))
            rec = {
                "length": length,
                "path": "/".join(self.path),
                "n_hits": len(hits),
                "hits": hits[:50],
            }
            if keep:
                rec["strings"] = decoded
                self.str_vecs.append(rec)
                for s, i in hits:
                    self.gene_hits[s].append(f"{rec['path']}#{i}")
                if rec["path"].endswith("class") or rec["path"] == "class":
                    key = "/".join(decoded)
                    if key and key not in self.classes:
                        self.classes.append(key)
                        self.emit(f"class {key} path={rec['path']}")
            if hits:
                self.emit(f"STR n={length} hits={[h[0] for h in hits]} path={rec['path']}")
            return {"type": "STR", "length": length}
        # VEC / EXPR
        items = []
        keep = length <= 200
        for i in range(length):
            self.path.append(f"[{i}]")
            el = self.parse_obj()
            self.path.pop()
            if keep:
                items.append(el)
        return items if keep else {"vec_len": length}

    def _int_array(self, signed: bool = True):
        length = self.s.i32()
        nbytes = length * 4
        path = "/".join(self.path)
        if length > LARGE:
            dest = self.tmp / f"arr_{len(self.arrays):04d}_int32.npy"
            self._dump_numeric(dest, length, 4, ">i4", np.int32)
            rec = {"id": len(self.arrays), "dtype": "int32", "length": length, "path": path, "file": str(dest)}
            self.arrays.append(rec)
            self.emit(f"INT n={length} -> {dest.name} path={path}")
            return rec
        raw = self.s.read(nbytes)
        arr = np.frombuffer(raw, dtype=">i4").astype(np.int32, copy=True)
        rec = {"dtype": "int32", "length": length, "path": path, "values": arr}
        if length <= 8 or path:
            self.small_int.append({"path": path, "length": length, "head": arr[:8].tolist()})
        return rec

    def _real_array(self):
        length = self.s.i32()
        path = "/".join(self.path)
        if length > LARGE:
            dest = self.tmp / f"arr_{len(self.arrays):04d}_f64.npy"
            self._dump_numeric(dest, length, 8, ">f8", np.float64)
            rec = {"id": len(self.arrays), "dtype": "float64", "length": length, "path": path, "file": str(dest)}
            self.arrays.append(rec)
            self.emit(f"REAL n={length} -> {dest.name} path={path}")
            return rec
        raw = self.s.read(length * 8)
        arr = np.frombuffer(raw, dtype=">f8").astype(np.float64, copy=True)
        return {"dtype": "float64", "length": length, "path": path, "values": arr}

    def _dump_numeric(self, dest: Path, length: int, item: int, be_dtype: str, out_dtype):
        dest.parent.mkdir(parents=True, exist_ok=True)
        mm = np.lib.format.open_memmap(dest, mode="w+", dtype=out_dtype, shape=(length,))
        left = length
        off = 0
        chunk_n = 1_000_000
        while left:
            n = min(left, chunk_n)
            raw = self.s.read(n * item)
            mm[off : off + n] = np.frombuffer(raw, dtype=be_dtype)
            off += n
            left -= n
        mm.flush()
        del mm

    def _bytecode(self):
        n_rep = self.s.i32()
        self.parse_obj()
        n_const = self.s.i32()
        for _ in range(n_const):
            inner = self.s.i32()
            self.parse_obj(inner)
        return {"bcode": n_rep}

    def _expand_altrep(self, info, state):
        name = self._altrep_name(info)
        if name in (b"compact_intseq", b"compact_realseq"):
            # state is REAL vector [n, start, step] typically
            vals = None
            if isinstance(state, dict) and isinstance(state.get("value"), dict):
                vals = state["value"].get("values")
            if vals is None and isinstance(state, dict) and "values" in state:
                vals = state["values"]
            return {"altrep": name.decode(), "state": "compact"}
        return {"altrep": name.decode() if isinstance(name, (bytes, bytearray)) else str(name)}

    def _altrep_name(self, info) -> bytes:
        # info is LIST: CAR = class SYM
        if isinstance(info, dict) and info.get("type") == "LIST":
            car = info.get("value")
            if isinstance(car, tuple):
                return self._sym_bytes(car[0]) or b"unknown"
        return b"unknown"

    def _sym_bytes(self, obj) -> bytes | None:
        if not isinstance(obj, dict):
            return None
        if obj.get("type") == "REF":
            ref = obj.get("value", {}).get("ref")
            if ref and ref < len(self.refs):
                return self._sym_bytes(self.refs[ref])
        if obj.get("type") == "SYM":
            ch = obj.get("value")
            if isinstance(ch, dict) and isinstance(ch.get("value"), (bytes, bytearray)):
                return bytes(ch["value"])
        if obj.get("type") == "CHAR" and isinstance(obj.get("value"), (bytes, bytearray)):
            return bytes(obj["value"])
        return None


def pairlist_to_items(node) -> list[tuple[str | None, object]]:
    items = []
    seen = 0
    while isinstance(node, dict) and node.get("type") in ("LIST", "LANG") and seen < 500:
        seen += 1
        val = node.get("value")
        if not (isinstance(val, tuple) and len(val) == 2):
            break
        car, cdr = val
        name = None
        tag = node.get("tag")
        if isinstance(tag, dict):
            if tag.get("type") == "SYM":
                ch = tag.get("value")
                if isinstance(ch, dict) and isinstance(ch.get("value"), (bytes, bytearray)):
                    name = ch["value"].decode("utf-8", "replace")
        items.append((name, car))
        node = cdr
        if isinstance(node, dict) and node.get("type") in ("NILVALUE", "NIL"):
            break
    return items


def find_gene_vector(str_vecs: list[dict], panel: list[str]) -> dict | None:
    panel_set = set(panel)
    best = None
    best_n = -1
    for rec in str_vecs:
        strings = rec.get("strings")
        if not strings:
            continue
        hits = panel_set.intersection(strings)
        if len(hits) > best_n:
            best_n = len(hits)
            best = rec
    return best


def find_barcode_vector(str_vecs: list[dict], n_cells_guess: int | None) -> dict | None:
    best = None
    best_score = -1
    for rec in str_vecs:
        strings = rec.get("strings")
        if not strings or len(strings) < 1000:
            continue
        sample = strings[:200]
        score = sum(1 for s in sample if "-" in s or s.endswith("-1") or "MRC" in s)
        if n_cells_guess and abs(len(strings) - n_cells_guess) < 10:
            score += 100
        if score > best_score:
            best_score = score
            best = rec
    return best


def extract_genes_from_csc(
    i_path: Path,
    p_path: Path,
    x_path: Path,
    gene_idx: dict[str, int],
    n_cells: int,
) -> dict[str, np.ndarray]:
    i = np.load(i_path, mmap_mode="r")
    p = np.load(p_path, mmap_mode="r")
    x = np.load(x_path, mmap_mode="r")
    want = {idx: name for name, idx in gene_idx.items()}
    out = {name: np.zeros(n_cells, dtype=np.float32) for name in gene_idx}
    # p is CSC column pointers; length should be n_cells+1
    if len(p) != n_cells + 1:
        # try n_cells from p
        n_cells = len(p) - 1
        out = {name: np.zeros(n_cells, dtype=np.float32) for name in gene_idx}
    for col in range(n_cells):
        a, b = int(p[col]), int(p[col + 1])
        if b <= a:
            continue
        rows = i[a:b]
        vals = x[a:b]
        for row, val in zip(rows.tolist(), vals.tolist()):
            name = want.get(int(row))
            if name is not None:
                out[name][col] = float(val)
        if col and col % 50000 == 0:
            print(f"  CSC col {col}/{n_cells}", flush=True)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--rds", type=Path, default=FILES["GSE253013_rds"])
    ap.add_argument("--outdir", type=Path, default=EXTRACTED / "GSE253013")
    args = ap.parse_args()
    out = args.outdir
    tmp = out / "tmp_arrays"
    out.mkdir(parents=True, exist_ok=True)
    tmp.mkdir(parents=True, exist_ok=True)

    log_path = out / "extract_log.txt"
    log = log_path.open("w")
    print(f"opening {args.rds}", flush=True)
    with gzip.open(args.rds, "rb") as f1, gzip.open(f1, "rb") as f2:
        s = Stream(f2)
        magic = s.read(2)
        if magic != b"X\n":
            raise SystemExit(f"not XDR RDS, magic={magic!r}")
        ver = s.i32()
        rver = s.i32()
        minver = s.i32()
        log.write(f"RDS format={ver} r={rver:#x} min={minver:#x}\n")
        if ver >= 3:
            enc_len = s.i32()
            enc = s.read(enc_len)
            log.write(f"encoding={enc!r}\n")
        ex = Extractor(s, tmp, GENE_PANEL, log)
        try:
            root = ex.parse_obj()
        except Exception as exc:
            log.write(f"PARSE_FAIL at byte {s.n} objs={ex.n_obj}: {exc}\n")
            log.close()
            print(f"PARSE_FAIL at byte {s.n}: {exc}", file=sys.stderr)
            raise

    log.write(f"parsed bytes={s.n} objects={ex.n_obj} arrays={len(ex.arrays)} str_vecs={len(ex.str_vecs)}\n")
    log.write("classes: " + ", ".join(ex.classes) + "\n")
    log.write("gene_hits: " + json.dumps(ex.gene_hits) + "\n")
    log.close()

    (out / "arrays.json").write_text(json.dumps(ex.arrays, indent=2) + "\n")
    (out / "small_int.json").write_text(json.dumps(ex.small_int, indent=2) + "\n")
    (out / "classes.json").write_text(json.dumps(ex.classes, indent=2) + "\n")
    (out / "gene_hits.json").write_text(json.dumps(ex.gene_hits, indent=2) + "\n")

    # Persist string vectors of interest (genes, barcodes, metadata-like).
    str_index = []
    for i, rec in enumerate(ex.str_vecs):
        strings = rec.get("strings")
        if not strings:
            continue
        fname = out / f"str_{i:03d}_n{len(strings)}.txt"
        # only write if it looks useful
        useful = bool(rec.get("hits")) or any(
            k in (rec.get("path") or "").lower()
            for k in ("ident", "sample", "patient", "cell", "type", "garnett", "tissue", "group")
        )
        n_unique = len(set(strings[:1000]))
        if useful or n_unique < min(200, len(strings)):
            fname.write_text("\n".join(strings) + "\n")
            rec_out = {k: rec[k] for k in rec if k != "strings"}
            rec_out["file"] = str(fname)
            rec_out["n_unique_head"] = n_unique
            str_index.append(rec_out)
        else:
            rec_out = {k: rec[k] for k in rec if k != "strings"}
            rec_out["n_unique_head"] = n_unique
            str_index.append(rec_out)
    (out / "str_index.json").write_text(json.dumps(str_index, indent=2) + "\n")

    gene_vec = find_gene_vector(ex.str_vecs, GENE_PANEL)
    gene_present = []
    gene_idx = {}
    if gene_vec and gene_vec.get("strings"):
        for i, g in enumerate(gene_vec["strings"]):
            if g in ex.want:
                gene_idx[g] = i
                gene_present.append(g)
    (out / "gene_index.json").write_text(
        json.dumps(
            {
                "n_genes_in_vector": gene_vec["length"] if gene_vec else None,
                "present": gene_present,
                "absent": [g for g in GENE_PANEL if g not in gene_idx],
                "index": gene_idx,
            },
            indent=2,
        )
        + "\n"
    )

    # Infer sparse CSC: look for Dim-like small int[2] and matching p/i/x lengths.
    dims = [r for r in ex.small_int if r["length"] == 2]
    (out / "extract_summary.json").write_text(
        json.dumps(
            {
                "n_arrays": len(ex.arrays),
                "n_str_vecs": len(ex.str_vecs),
                "classes": ex.classes,
                "gene_present": gene_present,
                "both_targets_present": ("TACSTD2" in gene_idx and "CLDN4" in gene_idx),
                "either_target_present": ("TACSTD2" in gene_idx or "CLDN4" in gene_idx),
                "dim_candidates": dims,
                "array_lengths": [{"id": a["id"], "dtype": a["dtype"], "length": a["length"]} for a in ex.arrays],
            },
            indent=2,
        )
        + "\n"
    )
    print("extract walk complete; gene_present=", gene_present)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
