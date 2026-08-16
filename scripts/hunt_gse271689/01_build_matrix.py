"""Build a gene x AOI count matrix for GSE271689 from the deposited GeoMx DCC files.

Inputs (downloaded by 00_fetch_data.sh):
  data/dcc/*.dcc.gz                 raw DCC files from GSE271689_RAW.tar
  data/Hs_R_NGS_WTA_v1.0.pkc        WTA v1.0 probe kit configuration (RTS_ID -> target)
  data/GSE271689_series_matrix.txt.gz   GEO sample annotation (spot id, segment, treatment)

Outputs (data/processed/):
  counts_raw.tsv.gz     target x AOI deduplicated read counts
  aoi_annotation.tsv    per-AOI annotation + sequencing QC
  q3_log2.tsv.gz        Q3-normalised log2 expression (targets x AOIs passing QC)
"""

import gzip
import json
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(DATA, "processed")
os.makedirs(OUT, exist_ok=True)


def parse_dcc(path):
    """Return (attributes dict, {RTS_ID: count})."""
    attrs, counts = {}, {}
    section = None
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("<") and line.endswith(">"):
                section = line.strip("<>")
                continue
            if section == "Code_Summary":
                rts, val = line.split(",")
                counts[rts] = int(val)
            elif "," in line:
                k, v = line.split(",", 1)
                attrs[k] = v.strip('"')
    return attrs, counts


def load_pkc(path):
    pkc = json.load(open(path))
    rts2target, negative_rts = {}, set()
    for target in pkc["Targets"]:
        name = target["DisplayName"]
        is_neg = target["CodeClass"].startswith("Negative")
        for probe in target["Probes"]:
            rts2target[probe["RTS_ID"]] = name
            if is_neg:
                negative_rts.add(probe["RTS_ID"])
    return rts2target, negative_rts


def load_geo_annotation(path):
    fields = defaultdict(list)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if not line.startswith("!Sample_"):
                continue
            parts = line.rstrip("\n").split("\t")
            fields[parts[0]].append([p.strip('"') for p in parts[1:]])
    titles = fields["!Sample_title"][0]
    gsms = fields["!Sample_geo_accession"][0]
    chars = fields["!Sample_characteristics_ch1"]

    def char(idx, prefix):
        return [c.split(": ", 1)[1] if c.startswith(prefix) else np.nan for c in chars[idx]]

    ann = pd.DataFrame(
        {
            "gsm": gsms,
            "title": titles,
            "spot_id": char(0, "spotid"),
            "segment": char(2, "cell type"),
            "treatment": char(3, "treatment"),
        }
    )
    dsp = ann["title"].str.extract(r"(DSP-(\d+)-([A-Z])-([A-H]\d+))")
    ann["dsp_id"] = dsp[0]
    ann["plate"] = dsp[1]
    ann["well"] = dsp[3]
    return ann


def main():
    rts2target, negative_rts = load_pkc(os.path.join(DATA, "Hs_R_NGS_WTA_v1.0.pkc"))
    ann = load_geo_annotation(os.path.join(DATA, "GSE271689_series_matrix.txt.gz"))

    dcc_dir = os.path.join(DATA, "dcc")
    files = sorted(f for f in os.listdir(dcc_dir) if f.endswith(".dcc.gz"))
    print(f"parsing {len(files)} DCC files")

    counts_by_aoi, qc_rows = {}, []
    for fname in files:
        gsm = fname.split("_")[0]
        attrs, counts = parse_dcc(os.path.join(dcc_dir, fname))
        target_counts = defaultdict(int)
        unmapped = 0
        neg_probe_counts = {}
        for rts, val in counts.items():
            target = rts2target.get(rts)
            if target is None:
                unmapped += val
                continue
            if rts in negative_rts:
                neg_probe_counts[rts] = val
            else:
                target_counts[target] += val
        counts_by_aoi[gsm] = pd.Series(target_counts, dtype="int64")

        raw = float(attrs.get("Raw", np.nan))
        aligned = float(attrs.get("Aligned", np.nan))
        dedup = float(sum(counts.values()))
        neg_vals = np.array(list(neg_probe_counts.values()), dtype=float)
        qc_rows.append(
            {
                "gsm": gsm,
                "dsp_id": attrs.get("ID"),
                "raw_reads": raw,
                "trimmed": float(attrs.get("Trimmed", np.nan)),
                "stitched": float(attrs.get("Stitched", np.nan)),
                "aligned": aligned,
                "dedup_reads": dedup,
                "pct_aligned": 100 * aligned / raw if raw else np.nan,
                "pct_trimmed": 100 * float(attrs.get("Trimmed", np.nan)) / raw if raw else np.nan,
                "pct_stitched": 100 * float(attrs.get("Stitched", np.nan)) / raw if raw else np.nan,
                "seq_saturation": 100 * (1 - dedup / aligned) if aligned else np.nan,
                "unmapped_rts_counts": unmapped,
                "n_neg_probes": len(neg_probe_counts),
                # geometric mean of the 139 WTA negative probes
                "neg_geomean": float(np.exp(np.mean(np.log(neg_vals + 1))) - 1) if len(neg_vals) else np.nan,
                "neg_geosd": float(np.exp(np.std(np.log(neg_vals + 1)))) if len(neg_vals) else np.nan,
            }
        )

    counts = pd.DataFrame(counts_by_aoi).fillna(0).astype("int64")
    qc = pd.DataFrame(qc_rows)
    ann = ann.merge(qc, on="gsm", how="left")
    ann["n_targets_detected"] = (counts > 0).sum(axis=0).reindex(ann["gsm"]).values
    ann["total_counts"] = counts.sum(axis=0).reindex(ann["gsm"]).values
    # limit of quantitation, the GeoMx default: negative geomean * geoSD^2
    ann["loq"] = np.maximum(ann["neg_geomean"] * ann["neg_geosd"] ** 2, 2.0)

    # AOI QC: NanoString defaults for raw reads / alignment; NTC and failed wells drop out
    ann["qc_pass"] = (
        (ann["raw_reads"] >= 1000)
        & (ann["pct_aligned"] >= 75)
        & (ann["pct_trimmed"] >= 80)
        & (ann["pct_stitched"] >= 80)
        & (ann["seq_saturation"] >= 50)
        & (ann["segment"].isin(["CK", "CD45", "CD68"]))
        & (ann["total_counts"] >= 10000)
    )
    print(ann["qc_pass"].value_counts())

    keep = ann.loc[ann["qc_pass"], "gsm"].tolist()
    mat = counts[keep]
    q3 = mat.apply(lambda col: np.percentile(col[col > 0], 75), axis=0)
    q3norm = mat.divide(q3, axis=1) * q3.mean()
    log2q3 = np.log2(q3norm + 1)

    counts.to_csv(os.path.join(OUT, "counts_raw.tsv.gz"), sep="\t")
    ann.assign(q3=ann["gsm"].map(q3)).to_csv(os.path.join(OUT, "aoi_annotation.tsv"), sep="\t", index=False)
    log2q3.to_csv(os.path.join(OUT, "q3_log2.tsv.gz"), sep="\t")
    print("targets:", counts.shape[0], "AOIs:", counts.shape[1], "QC-pass AOIs:", len(keep))


if __name__ == "__main__":
    main()
