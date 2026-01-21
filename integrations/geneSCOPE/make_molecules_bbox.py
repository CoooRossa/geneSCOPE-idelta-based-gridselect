#!/usr/bin/env python3
"""
make_molecules_bbox.py

Stdlib-only helper to build a bbox-prefiltered molecules TSV for the external
`idelta-gridselect` tool by streaming a Xenium transcripts table.

Default Xenium column mapping:
  gene: feature_name
  x:    x_location
  y:    y_location

Output schema (TSV):
  gene    x_um    y_um
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import os
import sys
from pathlib import Path
from typing import Iterable, Tuple


DEFAULT_EXCLUDE_PREFIXES = (
    "Unassigned",
    "NegControl",
    "Background",
    "DeprecatedCodeword",
    "SystemControl",
    "Negative",
    "BlankCodeword",
    "Blank",
)


def _iter_non_comment_lines(path: Path) -> Iterable[str]:
    with path.open("r", newline="") as f:
        for line in f:
            if not line.strip():
                continue
            if line.lstrip().startswith("#"):
                continue
            yield line


def _sniff_csv_dialect(sample: str) -> csv.Dialect:
    try:
        return csv.Sniffer().sniff(sample, delimiters=[",", "\t"])
    except csv.Error:
        class _D(csv.Dialect):
            delimiter = ","
            quotechar = '"'
            doublequote = True
            skipinitialspace = False
            lineterminator = "\n"
            quoting = csv.QUOTE_MINIMAL

        return _D()


def read_roi_bbox(roi_csv: Path) -> Tuple[float, float, float, float]:
    lines = list(_iter_non_comment_lines(roi_csv))
    if not lines:
        raise RuntimeError(f"ROI file had no data rows after stripping comments: {roi_csv}")

    sample = "".join(lines[:50])
    dialect = _sniff_csv_dialect(sample)
    reader = csv.DictReader(io.StringIO("".join(lines)), dialect=dialect)
    if not reader.fieldnames:
        raise RuntimeError(f"ROI file missing header row: {roi_csv}")

    field_map = {name.strip().lower(): name for name in reader.fieldnames}
    x_key = field_map.get("x_um") or field_map.get("x")
    y_key = field_map.get("y_um") or field_map.get("y")
    if x_key is None or y_key is None:
        raise RuntimeError(
            f"ROI header must include x_um/y_um (or x/y); got: {reader.fieldnames}"
        )

    min_x = float("inf")
    max_x = float("-inf")
    min_y = float("inf")
    max_y = float("-inf")
    n = 0
    for row in reader:
        try:
            x = float(row[x_key])
            y = float(row[y_key])
        except Exception as e:
            raise RuntimeError(f"Invalid ROI vertex row: {row} ({e})") from e
        min_x = min(min_x, x)
        max_x = max(max_x, x)
        min_y = min(min_y, y)
        max_y = max(max_y, y)
        n += 1

    if n < 3 or not all(map(lambda v: v == v and abs(v) != float("inf"), [min_x, max_x, min_y, max_y])):
        raise RuntimeError(f"ROI bbox invalid (need >=3 vertices): {roi_csv}")

    return min_x, max_x, min_y, max_y


def open_maybe_gzip_text(path: Path):
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", newline="")
    return path.open("r", newline="")


def main() -> int:
    p = argparse.ArgumentParser(
        description="Stream Xenium transcripts.csv(.gz) and emit bbox-prefiltered molecules TSV (gene,x_um,y_um)."
    )
    p.add_argument("--roi_csv", required=True, help="ROI vertices CSV/TSV with x_um,y_um (comments allowed).")
    p.add_argument("--transcripts_gz", required=True, help="Xenium transcripts.csv.gz (or transcripts.csv).")
    p.add_argument("--out_tsv", required=True, help="Output TSV path (gene, x_um, y_um).")
    p.add_argument("--gene_col", default="feature_name", help="Gene column name in transcripts.")
    p.add_argument("--x_col", default="x_location", help="X column name in transcripts (µm).")
    p.add_argument("--y_col", default="y_location", help="Y column name in transcripts (µm).")
    p.add_argument(
        "--exclude_prefixes",
        default=",".join(DEFAULT_EXCLUDE_PREFIXES),
        help="Comma-separated list of feature_name prefixes to exclude as non-gene labels.",
    )
    args = p.parse_args()

    roi_csv = Path(os.path.expanduser(args.roi_csv)).resolve()
    transcripts = Path(os.path.expanduser(args.transcripts_gz)).resolve()
    out_tsv = Path(os.path.expanduser(args.out_tsv)).resolve()

    min_x, max_x, min_y, max_y = read_roi_bbox(roi_csv)
    print(
        f"[make_molecules_bbox] ROI bbox: x=[{min_x}, {max_x}] y=[{min_y}, {max_y}]",
        file=sys.stderr,
    )

    exclude_prefixes = [p.strip() for p in str(args.exclude_prefixes).split(",")]
    exclude_prefixes = [p for p in exclude_prefixes if p]
    if exclude_prefixes:
        print(
            "[make_molecules_bbox] exclude_prefixes=" + ",".join(exclude_prefixes),
            file=sys.stderr,
        )

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    n_in = 0
    n_out = 0
    n_excl_prefix = 0

    with open_maybe_gzip_text(transcripts) as fin, out_tsv.open("w", newline="") as fout:
        reader = csv.DictReader(fin)
        if reader.fieldnames is None:
            raise RuntimeError(f"transcripts file missing header: {transcripts}")

        for col in (args.gene_col, args.x_col, args.y_col):
            if col not in reader.fieldnames:
                raise RuntimeError(f"Missing transcripts column '{col}'. Available: {reader.fieldnames}")

        w = csv.writer(fout, delimiter="\t", lineterminator="\n")
        w.writerow(["gene", "x_um", "y_um"])

        for row in reader:
            n_in += 1
            try:
                x = float(row[args.x_col])
                y = float(row[args.y_col])
            except Exception:
                continue
            if x < min_x or x > max_x or y < min_y or y > max_y:
                continue
            gene = row[args.gene_col]
            if not gene:
                continue
            if exclude_prefixes and any(gene.startswith(pref) for pref in exclude_prefixes):
                n_excl_prefix += 1
                continue
            w.writerow([gene, x, y])
            n_out += 1

    print(
        f"[make_molecules_bbox] wrote {n_out} / {n_in} transcripts rows to {out_tsv} (excluded_by_prefix={n_excl_prefix})",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

