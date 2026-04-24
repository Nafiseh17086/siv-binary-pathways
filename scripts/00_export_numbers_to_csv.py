#!/usr/bin/env python3
"""
Convert an Apple Numbers workbook to CSV.

Exports *every* table in *every* sheet so you can inspect the full workbook,
then copies the main pathway-enrichment table to the expected CSV path.

Usage:
    python 00_export_numbers_to_csv.py <input.numbers> <output.csv>
"""
import argparse
import csv
import os
import re
import sys
from pathlib import Path


def _safe(s: str) -> str:
    s = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s)).strip("_")
    return s[:80] or "table"


def export_all(numbers_path: Path, out_dir: Path) -> list[tuple[str, str, int, int, str]]:
    try:
        from numbers_parser import Document
    except ImportError:
        sys.exit(
            "numbers-parser is not installed. "
            "Run: pip install numbers-parser"
        )

    doc = Document(str(numbers_path))
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = []
    for sheet in doc.sheets:
        for table in sheet.tables:
            fname = f"{_safe(sheet.name)}__{_safe(table.name)}.csv"
            path = out_dir / fname
            with path.open("w", newline="") as fh:
                w = csv.writer(fh)
                for row in table.rows():
                    w.writerow([c.value if c.value is not None else "" for c in row])
            manifest.append((sheet.name, table.name, table.num_rows, table.num_cols, fname))
    return manifest


MAIN_TABLE_RE = re.compile(r"AllSets.*Pathways.*Enrich", re.IGNORECASE)


def pick_main_csv(manifest, out_dir: Path) -> Path:
    """Pick the largest table whose name looks like the main enrichment table."""
    candidates = [m for m in manifest if MAIN_TABLE_RE.search(m[1])]
    if not candidates:
        # fall back to the single largest table by row count
        candidates = sorted(manifest, key=lambda m: -m[2])[:1]
    best = max(candidates, key=lambda m: m[2])
    return out_dir / best[-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input", help="Input .numbers file")
    ap.add_argument("output", help="Output .csv path for the main enrichment table")
    ap.add_argument("--all-dir", default=None,
                    help="Directory to dump every sheet/table as its own CSV "
                         "(default: sibling folder named '<input>_csv')")
    args = ap.parse_args()

    in_path = Path(args.input)
    out_path = Path(args.output)
    all_dir = Path(args.all_dir) if args.all_dir else in_path.with_suffix("").with_name(in_path.stem + "_csv")

    if not in_path.exists():
        sys.exit(f"Input not found: {in_path}")

    print(f"[00] Reading  : {in_path}")
    manifest = export_all(in_path, all_dir)
    print(f"[00] Exported {len(manifest)} tables to {all_dir}/")
    for sheet, table, rows, cols, fname in manifest:
        print(f"       {rows:>6}×{cols:<3}  {sheet:<45}  {table:<50}  -> {fname}")

    main_csv = pick_main_csv(manifest, all_dir)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_bytes(main_csv.read_bytes())
    print(f"[00] Main enrichment table -> {out_path}")


if __name__ == "__main__":
    main()
