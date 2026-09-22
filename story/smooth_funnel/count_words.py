#!/usr/bin/env python3
"""Word-count helper for parallel EN/ZH elevator pitches (~120 each)."""
from __future__ import annotations

import re
from pathlib import Path

TEXT = Path(__file__).with_name("ELEVATOR_PITCH.md").read_text(encoding="utf-8")


def extract(lang: str) -> str:
    if lang == "en":
        m = re.search(
            r"## English[^\n]*\n\n(.+?)\n\n\*\*Word count:\*\*",
            TEXT,
            re.S,
        )
    else:
        m = re.search(
            r"## 中文[^\n]*\n\n(.+?)\n\n\*\*词数",
            TEXT,
            re.S,
        )
    if not m:
        raise SystemExit(f"block not found: {lang}")
    # strip markdown bold markers for counting
    return re.sub(r"\*\*([^*]+)\*\*", r"\1", m.group(1)).strip()


def tokens(s: str) -> list[str]:
    """Whitespace words; keep Latin/CJK/num tokens, drop lone punctuation."""
    out = []
    for w in s.split():
        w = w.strip("，。：；、（）()《》\"“”‘’,.;:!?")
        if not w or w in {",", ".", ";", ":", "—", "-", "/", "、"}:
            continue
        out.append(w)
    return out


def main() -> None:
    en = extract("en")
    zh = extract("zh")
    ew, zw = tokens(en), tokens(zh)
    print(f"EN words: {len(ew)}")
    print(f"ZH words: {len(zw)}")
    if not (118 <= len(ew) <= 122):
        raise SystemExit(f"EN out of 120±2: {len(ew)}")
    if not (118 <= len(zw) <= 122):
        raise SystemExit(f"ZH out of 120±2: {len(zw)}")
    # required mentions
    blob = en + "\n" + zh
    need = [
        ("SKB264", "SKB264"),
        ("ICI", "ICI"),
        ("TROP2", "TROP2"),
        ("junction", "junction"),
        ("CLDN4", "CLDN4"),
        ("GSEA", "GSEA"),
        ("#1 or 第 1 / TJ#1 honesty", "TJ#1"),
    ]
    for label, key in need:
        if key not in blob and key.replace("#", "") not in blob:
            # allow Chinese 第1 without hash
            if key == "TJ#1" and ("TJ#1" in blob or "第 1" in zh or "第1" in blob):
                continue
            raise SystemExit(f"missing required mention: {label}")
    if "screen and pin" not in en.lower() and "screen and pin CLDN4" not in en:
        if "screen and pin" not in en:
            raise SystemExit("EN missing screen and pin")
    if "筛选" not in zh or "钉住" not in zh:
        raise SystemExit("ZH missing 筛选/钉住")
    print("OK — checklist + word band")


if __name__ == "__main__":
    main()
