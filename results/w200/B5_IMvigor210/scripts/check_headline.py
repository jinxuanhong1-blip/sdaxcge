#!/usr/bin/env python3
"""Sanity-check primary 2x2 arithmetic against headline_OR_n_p.tsv."""
import csv
from pathlib import Path

root = Path(__file__).resolve().parents[1]
two = list(csv.DictReader((root / "tables" / "orr_2x2_primary.tsv").open(), delimiter="\t"))
by = {r["CLDN4"]: r for r in two}
low_r, low_n = int(by["low"]["CR_PR"]), int(by["low"]["SD_PD"]) + int(by["low"]["CR_PR"])
high_r, high_n = int(by["high"]["CR_PR"]), int(by["high"]["SD_PD"]) + int(by["high"]["CR_PR"])
assert low_n == 149 and high_n == 149, (low_n, high_n)
assert low_r == 29 and high_r == 39, (low_r, high_r)
assert low_n + high_n == 298
# Unconditional MLE OR (same as logistic). Fisher reports the conditional MLE 1.465.
or_uncond = (39 / 110) / (29 / 120)
assert abs(or_uncond - 1.467085) < 1e-5, or_uncond
head = {r["item"]: r["value"] for r in csv.DictReader((root / "tables" / "headline_OR_n_p.tsv").open(), delimiter="\t")}
assert head["primary_ORR_n_eval"] == "298"
assert head["primary_ORR_n_high_resp_over_n_high"] == "39/149"
assert head["primary_ORR_n_low_resp_over_n_low"] == "29/149"
assert head["primary_ORR_fisher_OR_high_vs_low"] == "1.465"
assert head["primary_OS_n"] == "347"
assert head["primary_OS_HR_high_vs_low"] == "0.890"
print("headline 2x2 checks OK; uncond OR", round(or_uncond, 3), "Fisher headline 1.465; n=298")
