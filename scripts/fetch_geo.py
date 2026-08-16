#!/usr/bin/env python3
"""Download public GEO supplementary files for GSE190265 / GSE190266."""
import os
import urllib.request

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn"

FILES = {
    "GSE190265": [
        "GSE190265_TPM_France3.csv.gz",
        "GSE190265_samples_info_France3.csv.gz",
        "GSE190265_series_matrix.txt.gz",
    ],
    "GSE190266": [
        "GSE190266_TPM_France4.csv.gz",
        "GSE190266_series_matrix.txt.gz",
    ],
}


def url_for(acc, fname):
    if fname.endswith("_series_matrix.txt.gz"):
        return f"{BASE}/{acc}/matrix/{fname}"
    return f"{BASE}/{acc}/suppl/{fname}"


def main():
    for acc, files in FILES.items():
        dest = os.path.join(REPO, "data", acc)
        os.makedirs(dest, exist_ok=True)
        for fname in files:
            path = os.path.join(dest, fname)
            if os.path.exists(path) and os.path.getsize(path) > 0:
                print(f"exists {acc}/{fname}")
                continue
            u = url_for(acc, fname)
            print(f"fetch {u}")
            urllib.request.urlretrieve(u, path)
            print(f"  -> {path} ({os.path.getsize(path)} bytes)")


if __name__ == "__main__":
    main()
