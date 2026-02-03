#!/usr/bin/env sh
set -eu

usage() {
  cat <<'EOF' >&2
run_idelta_gridselect.sh

One-step runner for local `idelta-gridselect recommend-grid` that:
  - normalizes ROI vertices
  - generates bbox-prefiltered molecules TSV
  - runs the idelta-gridselect binary from this repo (or builds it if configured)
  - writes outputs under <out>/tool and logs under <out>/logs
  - prints: RECOMMENDED_GRID_UM=<value>

Usage:
  integrations/geneSCOPE/run_idelta_gridselect.sh \
    --xenium_dir <outs_dir> | --transcripts <transcripts.csv.gz|csv|parquet> \
    --roi_csv <roi_vertices_csv> \
    --out <out_dir> \
    [--knee-search-min-um <int>] [--descending-delta-um <int>] \
    [--max-width-um <int>] \
    [--threads <int>] \
    [--skip_molecules] [--skip_tool] \
    [--repo_pin <path>]

Notes:
  - If --repo_pin is provided, BIN_REL_PATH and BUILD_IF_MISSING are read from it.
  - This script never clones or pulls git; repo presence is handled outside (e.g., by geneSCOPE R wrappers).
EOF
}

err() {
  echo "ERROR: $*" >&2
  exit 1
}

abspath() {
  python3 -c 'import os,sys; print(os.path.abspath(os.path.expanduser(sys.argv[1])))' "$1"
}

XENIUM_DIR=""
TRANSCRIPTS=""
ROI_CSV=""
OUTDIR=""
THREADS="32"
KNEE_SEARCH_MIN_UM="5"
DESCENDING_DELTA_UM="2"
MAX_WIDTH_UM=""
SKIP_MOLECULES="0"
SKIP_TOOL="0"
REPO_PIN=""

while [ $# -gt 0 ]; do
  case "$1" in
    --xenium_dir) XENIUM_DIR="$2"; shift 2 ;;
    --transcripts) TRANSCRIPTS="$2"; shift 2 ;;
    --roi_csv) ROI_CSV="$2"; shift 2 ;;
    --out) OUTDIR="$2"; shift 2 ;;
    --knee-search-min-um) KNEE_SEARCH_MIN_UM="$2"; shift 2 ;;
    --descending-delta-um) DESCENDING_DELTA_UM="$2"; shift 2 ;;
    --max-width-um) MAX_WIDTH_UM="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --skip_molecules) SKIP_MOLECULES="1"; shift 1 ;;
    --skip_tool) SKIP_TOOL="1"; shift 1 ;;
    --repo_pin) REPO_PIN="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) err "Unknown argument: $1 (use --help)" ;;
  esac
done

[ -n "$ROI_CSV" ] || err "--roi_csv is required"
[ -n "$OUTDIR" ] || err "--out is required"

command -v python3 >/dev/null 2>&1 || err "python3 is required"

ROI_CSV="$(abspath "$ROI_CSV")"
OUTDIR="$(abspath "$OUTDIR")"
[ -f "$ROI_CSV" ] || err "roi_csv not found: $ROI_CSV"

if [ -n "$XENIUM_DIR" ] && [ -n "$TRANSCRIPTS" ]; then
  err "Provide only one of: --xenium_dir or --transcripts"
fi

if [ -n "$XENIUM_DIR" ]; then
  XENIUM_DIR="$(abspath "$XENIUM_DIR")"
  [ -d "$XENIUM_DIR" ] || err "xenium_dir not found: $XENIUM_DIR"
  if [ -f "$XENIUM_DIR/transcripts.csv.gz" ]; then
    TRANSCRIPTS="$XENIUM_DIR/transcripts.csv.gz"
  elif [ -f "$XENIUM_DIR/transcripts.csv" ]; then
    TRANSCRIPTS="$XENIUM_DIR/transcripts.csv"
  elif [ -f "$XENIUM_DIR/transcripts.parquet" ]; then
    TRANSCRIPTS="$XENIUM_DIR/transcripts.parquet"
  elif [ -f "$XENIUM_DIR/transcripts.pq" ]; then
    TRANSCRIPTS="$XENIUM_DIR/transcripts.pq"
  else
    err "Could not find transcripts.csv.gz (or transcripts.csv / transcripts.parquet) under: $XENIUM_DIR"
  fi
fi

[ -n "$TRANSCRIPTS" ] || err "Provide either --xenium_dir or --transcripts"
TRANSCRIPTS="$(abspath "$TRANSCRIPTS")"
[ -f "$TRANSCRIPTS" ] || err "transcripts not found: $TRANSCRIPTS"

SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
REPO_ROOT="$(CDPATH= cd -- "$SCRIPT_DIR/../.." && pwd)"

# Defaults; may be overridden by --repo_pin.
BIN_REL_PATH="target/release/idelta-gridselect"
BUILD_IF_MISSING="0"

if [ -n "$REPO_PIN" ]; then
  REPO_PIN="$(abspath "$REPO_PIN")"
  [ -f "$REPO_PIN" ] || err "repo pin file not found: $REPO_PIN"
  while IFS= read -r line || [ -n "$line" ]; do
    case "$line" in
      "" ) continue ;;
      \#* ) continue ;;
    esac
    key="${line%%=*}"
    val="${line#*=}"
    case "$key" in
      BIN_REL_PATH) BIN_REL_PATH="$val" ;;
      BUILD_IF_MISSING) BUILD_IF_MISSING="$val" ;;
    esac
  done <"$REPO_PIN"
fi

BIN_PATH="$REPO_ROOT/$BIN_REL_PATH"
if [ ! -x "$BIN_PATH" ]; then
  if [ "$BUILD_IF_MISSING" = "1" ]; then
    command -v cargo >/dev/null 2>&1 || err "Binary missing and cargo not found; install Rust or build manually: $REPO_ROOT"
    (cd "$REPO_ROOT" && cargo build --release) >/dev/null 2>&1 || err "cargo build --release failed in $REPO_ROOT"
  else
    err "Binary not found: $BIN_PATH. Build it in the repo (cargo build --release) or set BUILD_IF_MISSING=1 in repo pin."
  fi
fi
[ -x "$BIN_PATH" ] || err "Binary still not executable: $BIN_PATH"

STAGING_DIR="$OUTDIR/staging"
LOGS_DIR="$OUTDIR/logs"
TOOL_OUT_DIR="$OUTDIR/tool"
mkdir -p "$STAGING_DIR" "$LOGS_DIR" "$TOOL_OUT_DIR"

ROI_NORM="$STAGING_DIR/roi_vertices.csv"
MOLS_BBOX="$STAGING_DIR/molecules_bbox.tsv"
MOLS_LOG="$LOGS_DIR/make_molecules_bbox.log"
TOOL_LOG="$LOGS_DIR/idelta_gridselect.log"
KNEE_SUMMARY="$TOOL_OUT_DIR/knee_summary.tsv"

# Optional max-width cap via generated widths files (passed to tool).
WIDTHS_PRESCREEN=""
WIDTHS_ANCHOR=""
if [ -n "$MAX_WIDTH_UM" ]; then
  WIDTHS_PRESCREEN="$STAGING_DIR/widths_prescreen_v2_max${MAX_WIDTH_UM}.txt"
  WIDTHS_ANCHOR="$STAGING_DIR/widths_anchor_v2_max${MAX_WIDTH_UM}.txt"

  if [ ! -s "$WIDTHS_PRESCREEN" ]; then
    python3 - "$WIDTHS_PRESCREEN" "$MAX_WIDTH_UM" <<'PY'
import sys
path=sys.argv[1]; mx=int(float(sys.argv[2]))
base=[1,2,4,8,16,24,32,40,50,60,80,100,120,150,200,250,300]
out=[w for w in base if w<=mx]
if len(out)<3:
    raise SystemExit("max-width-um too small: prescreen widths would have <3 entries")
with open(path,"w") as f:
    f.write("# Generated by idelta-gridselect integrations/geneSCOPE/run_idelta_gridselect.sh\n")
    f.write(f"# base=widths_prescreen_v2 max_width_um={mx}\n")
    for w in out:
        f.write(str(w)+"\n")
PY
  fi
  if [ ! -s "$WIDTHS_ANCHOR" ]; then
    python3 - "$WIDTHS_ANCHOR" "$MAX_WIDTH_UM" <<'PY'
import sys
path=sys.argv[1]; mx=int(float(sys.argv[2]))
base=list(range(1,41))+[45,50,55,60,65,70,80,90,100,110,120,130,140,150,175,200,225,250,275,300]
out=[w for w in base if w<=mx]
if len(out)<6:
    raise SystemExit("max-width-um too small: anchor widths would have <6 entries (min_points=6)")
with open(path,"w") as f:
    f.write("# Generated by idelta-gridselect integrations/geneSCOPE/run_idelta_gridselect.sh\n")
    f.write(f"# base=widths_anchor_v2 max_width_um={mx}\n")
    for w in out:
        f.write(str(w)+"\n")
PY
  fi
fi

# 1) Normalize ROI vertices (strip # comments; normalize header X,Y -> x_um,y_um).
if [ ! -s "$ROI_NORM" ] || [ "$ROI_NORM" -ot "$ROI_CSV" ]; then
  python3 - "$ROI_CSV" "$ROI_NORM" <<'PY'
import csv, io, sys
from pathlib import Path

inp = Path(sys.argv[1])
out = Path(sys.argv[2])

lines = []
with inp.open("r", newline="") as f:
    for line in f:
        if not line.strip():
            continue
        if line.lstrip().startswith("#"):
            continue
        lines.append(line)
if not lines:
    raise SystemExit(f"ROI file had no data after stripping comments: {inp}")

sample = "".join(lines[:50])
try:
    dialect = csv.Sniffer().sniff(sample, delimiters=[",", "\t"])
except csv.Error:
    dialect = csv.excel

r = csv.reader(io.StringIO("".join(lines)), dialect=dialect)
hdr = next(r)
hdr_low = [h.strip().lower() for h in hdr]
if len(hdr_low) < 2:
    raise SystemExit(f"ROI header must have >=2 columns; got: {hdr}")

if hdr_low[0] not in ("x_um", "x") or hdr_low[1] not in ("y_um", "y"):
    raise SystemExit(f"ROI header must start with X,Y or x_um,y_um; got: {hdr}")

out.parent.mkdir(parents=True, exist_ok=True)
with out.open("w", newline="") as f:
    w = csv.writer(f, delimiter=",", lineterminator="\n")
    w.writerow(["x_um", "y_um"])
    n = 0
    for row in r:
        if not row or len(row) < 2:
            continue
        x = float(row[0])
        y = float(row[1])
        w.writerow([x, y])
        n += 1
if n < 3:
    raise SystemExit(f"ROI needs >=3 vertices; got {n}: {inp}")
PY
fi

# 2) Generate bbox-prefiltered molecules (tool-side helper; supports CSV/TSV(.gz) and Parquet transcripts).

if [ "$SKIP_MOLECULES" = "1" ]; then
  [ -s "$MOLS_BBOX" ] || err "--skip_molecules was set but missing: $MOLS_BBOX"
else
  if [ ! -s "$MOLS_BBOX" ] || [ "$MOLS_BBOX" -ot "$TRANSCRIPTS" ] || [ "$MOLS_BBOX" -ot "$ROI_NORM" ]; then
    rm -f "$MOLS_LOG" 2>/dev/null || true
    "$BIN_PATH" make-molecules-bbox \
      --roi-csv "$ROI_NORM" \
      --transcripts "$TRANSCRIPTS" \
      --out-tsv "$MOLS_BBOX" >"$MOLS_LOG" 2>&1 || {
        tail -n 50 "$MOLS_LOG" >&2 || true
        err "molecules bbox extraction failed; see $MOLS_LOG"
      }
  fi
fi

# 3) Run tool (or skip if up-to-date).
if [ "$SKIP_TOOL" = "1" ]; then
  [ -s "$KNEE_SUMMARY" ] || err "--skip_tool was set but missing: $KNEE_SUMMARY"
else
  tool_uptodate=0
  if [ -s "$KNEE_SUMMARY" ] && [ "$KNEE_SUMMARY" -nt "$MOLS_BBOX" ] && [ "$KNEE_SUMMARY" -nt "$ROI_NORM" ]; then
    tool_uptodate=1
  fi
  if [ "$tool_uptodate" -eq 1 ] && [ -n "$MAX_WIDTH_UM" ] && ( [ "$KNEE_SUMMARY" -ot "$WIDTHS_PRESCREEN" ] || [ "$KNEE_SUMMARY" -ot "$WIDTHS_ANCHOR" ] ); then
    tool_uptodate=0
  fi
  if [ "$tool_uptodate" -eq 0 ]; then
    rm -f "$TOOL_LOG" 2>/dev/null || true
    if [ -n "$MAX_WIDTH_UM" ]; then
      "$BIN_PATH" recommend-grid \
        --molecules "$MOLS_BBOX" \
        --roi "$ROI_NORM" \
        --out "$TOOL_OUT_DIR" \
        --w0-um 1 \
        --tile-um 250 \
        --prescreen-fraction 0.25 \
        --prescreen-reps 25 \
        --seed 1 \
        --knee-window-um 20 \
        --knee-search-min-um "$KNEE_SEARCH_MIN_UM" \
        --descending-delta-um "$DESCENDING_DELTA_UM" \
        --threads "$THREADS" \
        --b-min 50 \
        --n-min 100 \
        --min-points 6 \
        --min-informative-frac-ge2 0.6 \
        --widths-prescreen "$WIDTHS_PRESCREEN" \
        --widths-anchor "$WIDTHS_ANCHOR" \
        >"$TOOL_LOG" 2>&1 || {
          tail -n 50 "$TOOL_LOG" >&2 || true
          err "external tool failed; see $TOOL_LOG"
        }
    else
      "$BIN_PATH" recommend-grid \
        --molecules "$MOLS_BBOX" \
        --roi "$ROI_NORM" \
        --out "$TOOL_OUT_DIR" \
        --w0-um 1 \
        --tile-um 250 \
        --prescreen-fraction 0.25 \
        --prescreen-reps 25 \
        --seed 1 \
        --knee-window-um 20 \
        --knee-search-min-um "$KNEE_SEARCH_MIN_UM" \
        --descending-delta-um "$DESCENDING_DELTA_UM" \
        --threads "$THREADS" \
        --b-min 50 \
        --n-min 100 \
        --min-points 6 \
        --min-informative-frac-ge2 0.6 \
        >"$TOOL_LOG" 2>&1 || {
          tail -n 50 "$TOOL_LOG" >&2 || true
          err "external tool failed; see $TOOL_LOG"
        }
    fi
  fi
fi

# 4) Emit RECOMMENDED_GRID_UM=<value> by parsing knee_summary.tsv.
[ -f "$KNEE_SUMMARY" ] || err "missing tool output: $KNEE_SUMMARY"
REC="$(python3 - "$KNEE_SUMMARY" <<'PY'
import csv, sys
path=sys.argv[1]
with open(path, newline="") as f:
    r=csv.DictReader(f, delimiter="\t")
    for row in r:
        if row.get("gene") == "__SUMMARY__":
            v=row.get("recommended_grid_um","")
            if v is None or v == "" or v.lower() == "na":
                raise SystemExit(1)
            print(v)
            raise SystemExit(0)
raise SystemExit(1)
PY
)" || err "could not parse recommended_grid_um from $KNEE_SUMMARY"

echo "RECOMMENDED_GRID_UM=$REC"

exit 0
