#!/usr/bin/env python3
"""Verify + inventory the extra mouse-lung ICI RNA candidates found in the hunt.

Writes notes/fable_mouse_ici/hunt_candidates.json.
"""
import json, re, time, urllib.request
from pathlib import Path

CANDIDATES = [
    "GSE114300",  # PD-1 inhibition early-stage lung, n=62
    "GSE114601",  # BET + PD-1 GEMM NSCLC, n=8
    "GSE169196",  # anti-PD-1 resistance lung models, n=48
    "GSE190264",  # MEK + chemoIO LLC1, n=24
    "GSE157880",  # radiation + ICB bulk, n=16
    "GSE157883",  # radiation + ICB, n=22
    "GSE194166",  # AXL + PD-1 LKB1 NSCLC, n=10
    "GSE277610",  # endothelial PD-L1 NSCLC, n=17
    "GSE309199",  # HDAC + IO SCLC, n=12
    "GSE309192",  # HDAC + IO SCLC, n=6
    "GSE184000",  # whole lung anti-PD-1 irAE, n=12
    "GSE129298",  # SCLC scRNA exp2, n=4
    "GSE315010",  # KRAS(ON) + ICB, n=44 mixed
    "GSE303940",  # XIAP RNA-seq, n=2
]
NOTES = Path(__file__).resolve().parents[2] / "notes" / "fable_mouse_ici"
UA = {"User-Agent": "Mozilla/5.0 (fable-mouse-ici hunt)"}
TWO_GB = 2 * 1024 ** 3


def fetch(url, tries=4, binary=False):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=60) as r:
                data = r.read()
                return data if binary else data.decode("utf-8", "replace")
        except Exception as e:
            last = e
            time.sleep(2 ** i)
    raise last


def head_size(url):
    try:
        req = urllib.request.Request(url, headers=UA, method="HEAD")
        with urllib.request.urlopen(req, timeout=40) as r:
            cl = r.headers.get("Content-Length")
            return int(cl) if cl else None
    except Exception:
        return None


def parse_geo_text(txt):
    fields = {}
    for line in txt.splitlines():
        m = re.match(r"^!(\w+?)_(\w+) = (.*)$", line)
        if not m:
            continue
        fields.setdefault(m.group(2), []).append(m.group(3))
    return fields


def main():
    out = []
    for acc in CANDIDATES:
        print("[hunt]", acc)
        rec = {"accession": acc}
        url = (f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"
               f"acc={acc}&targ=self&form=text&view=quick")
        try:
            txt = fetch(url)
            f = parse_geo_text(txt)
            rec["exists"] = True
            rec["title"] = " | ".join(f.get("title", []))
            rec["organism"] = f.get("sample_organism", []) or f.get("platform_organism", [])
            rec["type"] = f.get("type", [])
            rec["summary"] = " ".join(f.get("summary", []))[:600]
            rec["overall_design"] = " ".join(f.get("overall_design", []))[:600]
            rec["n_samples"] = len(f.get("sample_id", []))
            rec["supplementary_files"] = f.get("supplementary_file", [])
        except Exception as e:
            rec["exists"] = False
            rec["error"] = str(e)
            out.append(rec)
            continue

        nnn = acc[:-3] + "nnn"
        base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{nnn}/{acc}/suppl/"
        rec["files"] = []
        try:
            html = fetch(base)
            for m in re.finditer(r'<a href="([^"?/][^"]*)">', html):
                name = m.group(1)
                if name.endswith("/") or "hhs.gov" in name:
                    continue
                urlf = base + name
                sz = head_size(urlf)
                rec["files"].append({"name": name, "url": urlf, "size_bytes": sz,
                                     "skip_gt_2gb": bool(sz and sz > TWO_GB)})
                time.sleep(0.15)
        except Exception as e:
            rec["suppl_error"] = str(e)
        out.append(rec)
        time.sleep(0.3)

    (NOTES / "hunt_candidates.json").write_text(json.dumps(out, indent=2))
    print("\n==== HUNT SUMMARY ====")
    for r in out:
        print("=" * 70)
        print(r["accession"], "| exists", r.get("exists"), "| n=", r.get("n_samples"),
              "| org", r.get("organism"), "| type", r.get("type"))
        print("  title:", (r.get("title") or "")[:130])
        for f in r.get("files") or []:
            sz = f["size_bytes"]
            human = f"{sz/1024**2:.1f} MB" if sz else "?"
            flag = " SKIP>2GB" if f.get("skip_gt_2gb") else ""
            print(f"   {human:>10}  {f['name']}{flag}")


if __name__ == "__main__":
    main()
