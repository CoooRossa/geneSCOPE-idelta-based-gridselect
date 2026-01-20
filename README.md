# `idelta-gridselect`

Dedicated, environment-isolated CLI for **recommending an optimal grid size** (µm) using per-gene Morisita’s Iδ curves inside a polygon ROI.

This tool is intended to be consumed by a thin downstream wrapper (e.g. in geneSCOPE). It is **not** marketed as a general “Iδ curve generator”, though it computes curves internally.

## Public command

- `recommend-grid`: outputs a recommended grid width and a minimal diagnostics bundle.

## Build

```bash
cargo build --release
./target/release/idelta-gridselect --help
```

## Example

```bash
./target/release/idelta-gridselect recommend-grid \
  --molecules path/to/molecules.csv \
  --roi path/to/roi.geojson \
  --out out/gridselect_run_001 \
  --w0-um 1 \
  --tile-um 250 \
  --prescreen-fraction 0.25 \
  --prescreen-reps 25 \
  --seed 1 \
  --knee-window-um 20 \
  --threads 8
```

## Inputs

### Molecules (`--molecules`)

Accepted formats:
- CSV (`.csv`)
- TSV (`.tsv` / `.txt`)
- Parquet (`.parquet` / `.pq`)

Required columns:
- `gene` (string)
- `x_um` (float; micrometers) or `x`
- `y_um` (float; micrometers) or `y`

### ROI polygon (`--roi`)

Accepted formats:
- GeoJSON Polygon (`.geojson` / `.json`)
- WKT Polygon (`.wkt` / `.txt`): `POLYGON((x y, x y, ...))`
- Vertices table (`.csv` / `.tsv`) with header `x_um,y_um` (or `x,y`)

## Outputs (one `--out` directory, stable bundle)

- `final_curves.tsv.gz`
- `avg_curve.tsv`
- `knee_summary.tsv`

### `final_curves.tsv.gz` schema (TSV.GZ)

The downstream wrapper should treat these outputs as source of truth.

Columns:
- `gene`, `width_um`, `idelta_raw`, `idelta_smooth`, `status`
- `N_total`, `Q_bins`, `n_bins_ge2`, `frac_bins_ge2`, `max_bin_count`, `sum_n_n1`

### `knee_summary.tsv` schema (TSV)

This file contains:
- one summary row (`gene="__SUMMARY__"`) with overall statistics and the recommended grid size
- one row per gene with `gene, knee_um, gene_status` (other columns are `NA`)

## Degenerate “0/1 occupancy” mitigation

At very small widths many genes have bin counts dominated by `{0,1}`, leading to `sum_n_n1 = 0` and misleadingly small Iδ.

For each gene and width, this tool records:
- `n_bins_ge2`, `frac_bins_ge2`, `sum_n_n1`, `max_bin_count`

A point is marked `status=degenerate_01` if either:
- `n_bins_ge2 < B_min` (default `--b-min 50`), or
- `sum_n_n1 == 0 AND N_total >= N_min` (default `--n-min 100`)

Degenerate points are still written for transparency, but `idelta_raw` is set to `NA` and they are excluded from smoothing/knee estimation.

A gene is excluded from knee estimation if it has fewer than `--min-points` non-degenerate widths.

The final recommendation is derived only from informative genes and includes a global sanity check:
at the recommended width, at least `--min-informative-frac-ge2` of informative genes must have `n_bins_ge2 >= B_min`; otherwise the tool deterministically falls back to a larger evaluated width and records the fallback in `knee_summary.tsv`.

## Reproducibility

- `--seed` controls prescreen tile sampling RNG; identical inputs + seed yield deterministic outputs after sorting.
- `--threads` controls rayon worker threads; output ordering is always sorted by `(gene, width_um)`.

## Knee algorithm

Fixed algorithm (no per-gene tuning):
- x = `log(width_um)`
- y = smoothed Iδ
- normalize x and y to `[0,1]`
- knee = point of maximum perpendicular distance to the chord connecting the endpoints

## Containers

- Dockerfile: `Dockerfile`
- Apptainer: `apptainer.def`

