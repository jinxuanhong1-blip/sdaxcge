#!/usr/bin/env python3
"""Extra public metadata: family XML peek, SuperSeries, file sizes, E-MTAB file JSON."""
from __future__ import annotations

import json
import tarfile
import time
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "metadata"
UA = "Mozilla/5.0 leftover-scrna-2426/1.0"


def get(url: str, dest: Path, retries: int = 4) -> None:
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=120) as r:
                dest.write_bytes(r.read())
            print(f"OK {dest.name} {dest.stat().st_size}")
            return
        except Exception as e:
            last = e
            time.sleep(2 ** (i + 1))
            print(f"retry {url} {e}")
    raise RuntimeError(last)


def head_size(url: str) -> str:
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": UA})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            return r.headers.get("Content-Length", "?")
    except Exception as e:
        return f"ERR:{e}"


def parse_family(acc: str) -> None:
    tgz = OUT / f"{acc}_family.xml.tgz"
    if not tgz.exists() or tgz.stat().st_size < 100:
        return
    with tarfile.open(tgz, "r:gz") as tf:
        names = tf.getnames()
        xml_name = [n for n in names if n.endswith(".xml")][0]
        f = tf.extractfile(xml_name)
        root = ET.parse(f).getroot()
    ns = ""
    if root.tag.startswith("{"):
        ns = root.tag.split("}")[0] + "}"
    rows = []
    for sample in root.iter(f"{ns}Sample"):
        title = ""
        gsm = sample.attrib.get("iid", "")
        chars = {}
        for child in sample:
            tag = child.tag.replace(ns, "")
            if tag == "Title":
                title = (child.text or "").strip()
            if tag == "Channel":
                for ch in child:
                    cht = ch.tag.replace(ns, "")
                    if cht == "Characteristics":
                        key = ch.attrib.get("tag", "char")
                        chars[key] = (ch.text or "").strip()
        rows.append({"gsm": gsm, "title": title, **chars})
    dest = OUT / f"{acc}_samples.json"
    dest.write_text(json.dumps(rows, indent=2))
    print(f"samples {acc}: {len(rows)} -> {dest.name}")


def main() -> None:
    for acc in [
        "GSE274584",
        "GSE274588",
        "GSE274595",
        "GSE302113",
        "GSE267108",
        "GSE176021",
        "GSE176022",
        "GSE186446",
        "GSE337519",
        "GSE308745",
    ]:
        parse_family(acc)

    extras = {
        "GSE185206.soft.txt": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185206&targ=self&form=text&view=brief",
        "GSE274596.soft.txt": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274596&targ=self&form=text&view=brief",
        "GSE173351.soft.txt": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173351&targ=self&form=text&view=brief",
        "GSE176021_CD3_annotations.rds.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176021/suppl/GSE176021_CD3_annotations.rds.gz",
        "GSE176021_CD8_annotations.rds.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176021/suppl/GSE176021_CD8_annotations.rds.gz",
        "GSE186446_fcount_aggr.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186446/suppl/GSE186446_fcount_aggr.txt.gz",
        "GSE274584_countTable.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274584/suppl/GSE274584_countTable.txt.gz",
        "E-MTAB-13526.files.json": "https://www.ebi.ac.uk/biostudies/files/E-MTAB-13526",
    }
    for name, url in extras.items():
        dest = OUT / name
        try:
            get(url, dest)
        except Exception as e:
            print(f"FAIL {name}: {e}")
        time.sleep(0.3)

    sizes = {
        "GSE267108_processed_data.tar.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE267nnn/GSE267108/suppl/GSE267108_processed_data.tar.gz",
        "GSE274595_RAW.tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274595/suppl/GSE274595_RAW.tar",
        "GSE337519_RAW.tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE337nnn/GSE337519/suppl/GSE337519_RAW.tar",
        "GSE186446_fcount_aggr.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186446/suppl/GSE186446_fcount_aggr.txt.gz",
        "tumour_h5ad": "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/135/E-MTAB-13526/Files/10X_Lung_Tumour_Annotated_v2.h5ad",
        "bg_h5ad": "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/135/E-MTAB-13526/Files/10X_Lung_Healthy_Background_Annotated_v2.h5ad",
    }
    rec = {k: head_size(v) for k, v in sizes.items()}
    (OUT / "file_sizes.json").write_text(json.dumps(rec, indent=2))
    print(json.dumps(rec, indent=2))


if __name__ == "__main__":
    main()
