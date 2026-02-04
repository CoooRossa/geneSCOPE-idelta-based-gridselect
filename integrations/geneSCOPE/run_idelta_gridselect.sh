#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
One-click grid-size recommendation (Xenium) using the external `idelta-gridselect` tool.

Usage:
  scripts/run_gridselect_xenium.sh \
    --xenium_dir <outs_dir> \
    --roi_csv <roi_vertices_csv> \
    --out <outdir> \
    [--tool-bin <path> | --docker-image <image> | --apptainer-sif <path>] \
    [--auto-find] [--search-root <dir> ...] \
    [--knee-search-min-um <int>] [--descending-delta-um <int>] \
    [--max-width-um <int>] \
    [--threads <n>] \
    [--skip-tool] [--skip-plots] \
    [--strict <0|1>]

External tool resolution (in order):
  1) --tool-bin
  2) IDELTA_GRIDSELECT_BIN
  3) --docker-image
  4) IDELTA_GRIDSELECT_DOCKER_IMAGE
  5) --apptainer-sif
  6) IDELTA_GRIDSELECT_APPTAINER_SIF
  7) --auto-find (search common local build paths)

Outputs under --out:
  staging/roi_vertices.csv
  staging/molecules_bbox.tsv
  tool/knee_summary.tsv (+ other tool outputs)
  plots/*.png
  logs/idelta_gridselect.log
  logs/idelta_gridselect.resolve.txt
  logs/make_molecules_bbox.log
  logs/wrapper.log
  logs/best_grid_summary.txt
EOF
}

err() { echo "ERROR: $*" >&2; exit 1; }

abspath() {
  python3 - "$1" <<'PY'
import os,sys
print(os.path.abspath(os.path.expanduser(sys.argv[1])))
PY
}

XENIUM_DIR=""
ROI_CSV=""
OUTDIR=""
THREADS=""
KNEE_SEARCH_MIN_UM="5"
DESCENDING_DELTA_UM="2"
MAX_WIDTH_UM=""
SKIP_TOOL=0
SKIP_PLOTS=0
STRICT=1
TOOL_BIN_CLI=""
DOCKER_IMAGE_CLI=""
APPTAINER_SIF_CLI=""
AUTO_FIND=0
SEARCH_ROOTS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --xenium_dir) XENIUM_DIR="$2"; shift 2 ;;
    --roi_csv) ROI_CSV="$2"; shift 2 ;;
    --out) OUTDIR="$2"; shift 2 ;;
    --knee-search-min-um) KNEE_SEARCH_MIN_UM="$2"; shift 2 ;;
    --descending-delta-um) DESCENDING_DELTA_UM="$2"; shift 2 ;;
    --max-width-um) MAX_WIDTH_UM="$2"; shift 2 ;;
    --tool-bin) TOOL_BIN_CLI="$2"; shift 2 ;;
    --docker-image) DOCKER_IMAGE_CLI="$2"; shift 2 ;;
    --apptainer-sif) APPTAINER_SIF_CLI="$2"; shift 2 ;;
    --auto-find) AUTO_FIND=1; shift ;;
    --search-root) SEARCH_ROOTS+=("$2"); shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --skip-tool) SKIP_TOOL=1; shift ;;
    --skip-plots) SKIP_PLOTS=1; shift ;;
    --strict) STRICT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) err "Unknown argument: $1 (use --help)" ;;
  esac
done

[[ -n "$XENIUM_DIR" ]] || err "--xenium_dir is required"
[[ -n "$ROI_CSV" ]] || err "--roi_csv is required"
[[ -n "$OUTDIR" ]] || err "--out is required"

XENIUM_DIR="$(abspath "$XENIUM_DIR")"
ROI_CSV="$(abspath "$ROI_CSV")"
OUTDIR="$(abspath "$OUTDIR")"

[[ -d "$XENIUM_DIR" ]] || err "xenium_dir not found: $XENIUM_DIR"
[[ -f "$ROI_CSV" ]] || err "roi_csv not found: $ROI_CSV"

TRANSCRIPTS=""
TRANSCRIPTS_KIND=""
if [[ -f "$XENIUM_DIR/transcripts.csv.gz" ]]; then
  TRANSCRIPTS="$XENIUM_DIR/transcripts.csv.gz"
  TRANSCRIPTS_KIND="csv"
elif [[ -f "$XENIUM_DIR/transcripts.csv" ]]; then
  TRANSCRIPTS="$XENIUM_DIR/transcripts.csv"
  TRANSCRIPTS_KIND="csv"
elif [[ -f "$XENIUM_DIR/transcripts.parquet" ]]; then
  TRANSCRIPTS="$XENIUM_DIR/transcripts.parquet"
  TRANSCRIPTS_KIND="parquet"
elif [[ -f "$XENIUM_DIR/transcripts.pq" ]]; then
  TRANSCRIPTS="$XENIUM_DIR/transcripts.pq"
  TRANSCRIPTS_KIND="parquet"
elif [[ -f "$XENIUM_DIR/transcripts.zarr.zip" ]]; then
  err "Found transcripts.zarr.zip under xenium_dir, but Zarr is not supported. Please export transcripts.parquet or transcripts.csv.gz: $XENIUM_DIR"
elif [[ -d "$XENIUM_DIR/transcripts.zarr" ]]; then
  err "Found transcripts.zarr under xenium_dir, but Zarr is not supported. Please export transcripts.parquet or transcripts.csv.gz: $XENIUM_DIR"
else
  err "Could not find transcripts.csv.gz (or transcripts.csv / transcripts.parquet) under: $XENIUM_DIR"
fi

if [[ -z "$THREADS" ]]; then
  if command -v nproc >/dev/null 2>&1; then
    THREADS="$(nproc)"
  elif command -v sysctl >/dev/null 2>&1; then
    THREADS="$(sysctl -n hw.ncpu)"
  else
    THREADS="1"
  fi
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
PIN_FILE="$REPO_ROOT/docs/idelta_gridselect_tool_version.txt"

STAGING_DIR="$OUTDIR/staging"
LOGS_DIR="$OUTDIR/logs"
TOOL_OUT_DIR="$OUTDIR/tool"
PLOTS_OUT_DIR="$OUTDIR/plots"

WRAPPER_LOG="$LOGS_DIR/wrapper.log"
RESOLVE_FILE="$LOGS_DIR/idelta_gridselect.resolve.txt"

ROI_NORM="$STAGING_DIR/roi_vertices.csv"
MOLS_BBOX="$STAGING_DIR/molecules_bbox.tsv"
TOOL_SUMMARY="$TOOL_OUT_DIR/knee_summary.tsv"
MOLS_LOG="$LOGS_DIR/make_molecules_bbox.log"

mkdir -p "$STAGING_DIR" "$LOGS_DIR" "$TOOL_OUT_DIR" "$PLOTS_OUT_DIR"

wlog() {
  echo "[gridselect] $(date '+%Y-%m-%d %H:%M:%S') $*" | tee -a "$WRAPPER_LOG"
}

T0="$(date +%s)"
wlog "outdir=$OUTDIR"
wlog "xenium_dir=$XENIUM_DIR"
wlog "transcripts=$TRANSCRIPTS"
wlog "roi_csv=$ROI_CSV"
wlog "threads=$THREADS"
wlog "knee_search_min_um=$KNEE_SEARCH_MIN_UM"
wlog "descending_delta_um=$DESCENDING_DELTA_UM"
wlog "max_width_um=${MAX_WIDTH_UM:-}"
wlog "pin_file=$PIN_FILE"

if [[ -f "$PIN_FILE" ]]; then
  wlog "tool_pin_begin"
  cat "$PIN_FILE" | tee -a "$WRAPPER_LOG"
  wlog "tool_pin_end"
else
  wlog "WARNING: missing pin file: $PIN_FILE"
fi

## 0) Resolve external tool early (fail-fast).
MODE=""
TOOL_CMD=()
RESOLVED_TOOL_BIN=""
RESOLVED_SOURCE=""

resolve_bin() {
  local source="$1"
  local bin="$2"
  [[ -n "$bin" ]] || return 1
  bin="$(abspath "$bin")"
  if [[ -x "$bin" ]]; then
    MODE="bin"
    TOOL_CMD=("$bin")
    RESOLVED_TOOL_BIN="$bin"
    RESOLVED_SOURCE="$source"
    return 0
  fi
  return 1
}

autofind_bin() {
  local candidates=()
  candidates+=("$PWD/target/release/idelta-gridselect")
  candidates+=("$PWD/../idelta-gridselect/target/release/idelta-gridselect")
  candidates+=("$REPO_ROOT/../idelta-gridselect/target/release/idelta-gridselect")

  local c=""
  for c in "${candidates[@]}"; do
    if [[ -x "$c" ]]; then
      echo "$c"
      return 0
    fi
  done

  local roots=()
  if [[ "${#SEARCH_ROOTS[@]}" -gt 0 ]]; then
    roots+=("${SEARCH_ROOTS[@]}")
  fi
  if [[ -d "$HOME/Documents/Code" ]]; then
    roots+=("$HOME/Documents/Code")
  fi

  local root=""
  local maxdepth="6"
  for root in "${roots[@]}"; do
    root="$(abspath "$root")"
    [[ -d "$root" ]] || continue
    local found=""
    found="$(find "$root" -maxdepth "$maxdepth" -type f -perm -111 -path "*/idelta-gridselect/target/release/idelta-gridselect" 2>/dev/null | head -n 1 || true)"
    if [[ -n "$found" && -x "$found" ]]; then
      echo "$found"
      return 0
    fi
  done

  return 1
}

if [[ -n "$TOOL_BIN_CLI" ]]; then
  resolve_bin "--tool-bin" "$TOOL_BIN_CLI" || err "--tool-bin was provided but not executable: $TOOL_BIN_CLI"
elif [[ -n "${IDELTA_GRIDSELECT_BIN:-}" ]]; then
  resolve_bin "IDELTA_GRIDSELECT_BIN" "${IDELTA_GRIDSELECT_BIN}" || err "IDELTA_GRIDSELECT_BIN is set but not executable: ${IDELTA_GRIDSELECT_BIN}"
fi

if [[ -z "$MODE" ]]; then
  DOCKER_IMAGE=""
  if [[ -n "$DOCKER_IMAGE_CLI" ]]; then
    DOCKER_IMAGE="$DOCKER_IMAGE_CLI"
    RESOLVED_SOURCE="--docker-image"
  elif [[ -n "${IDELTA_GRIDSELECT_DOCKER_IMAGE:-}" ]]; then
    DOCKER_IMAGE="${IDELTA_GRIDSELECT_DOCKER_IMAGE}"
    RESOLVED_SOURCE="IDELTA_GRIDSELECT_DOCKER_IMAGE"
  fi
  if [[ -n "$DOCKER_IMAGE" ]]; then
    MODE="docker"
    command -v docker >/dev/null 2>&1 || err "docker not found but docker execution was selected"
    TOOL_CMD=(docker run --rm -u "$(id -u)":"$(id -g)" -v "$OUTDIR":/out -v "$XENIUM_DIR":/xenium:ro -w /out --entrypoint idelta-gridselect "$DOCKER_IMAGE")
  fi
fi

if [[ -z "$MODE" ]]; then
  SIF=""
  if [[ -n "$APPTAINER_SIF_CLI" ]]; then
    SIF="$APPTAINER_SIF_CLI"
    RESOLVED_SOURCE="--apptainer-sif"
  elif [[ -n "${IDELTA_GRIDSELECT_APPTAINER_SIF:-}" ]]; then
    SIF="${IDELTA_GRIDSELECT_APPTAINER_SIF}"
    RESOLVED_SOURCE="IDELTA_GRIDSELECT_APPTAINER_SIF"
  fi
  if [[ -n "$SIF" ]]; then
    SIF="$(abspath "$SIF")"
    [[ -f "$SIF" ]] || err "Apptainer SIF not found: $SIF"
    MODE="apptainer"
    if command -v apptainer >/dev/null 2>&1; then
      TOOL_CMD=(apptainer run --bind "$OUTDIR":/out --bind "$XENIUM_DIR":/xenium "$SIF")
    elif command -v singularity >/dev/null 2>&1; then
      TOOL_CMD=(singularity run --bind "$OUTDIR":/out --bind "$XENIUM_DIR":/xenium "$SIF")
    else
      err "apptainer/singularity not found but apptainer execution was selected"
    fi
  fi
fi

if [[ -z "$MODE" && "$AUTO_FIND" -eq 1 ]]; then
  FOUND_BIN="$(autofind_bin || true)"
  if [[ -n "$FOUND_BIN" ]]; then
    resolve_bin "--auto-find" "$FOUND_BIN" || true
  fi
fi

if [[ -z "$MODE" ]]; then
  err "No external tool configured. Example: scripts/run_gridselect_xenium.sh --xenium_dir ... --roi_csv ... --out ... --tool-bin /path/to/idelta-gridselect"
fi

wlog "execution_mode=$MODE"

if [[ "$MODE" == "bin" && -n "$RESOLVED_TOOL_BIN" ]]; then
  echo "RESOLVED_TOOL_BIN=$RESOLVED_TOOL_BIN"
fi

{
  echo "resolved_at=$(date '+%Y-%m-%d %H:%M:%S')"
  echo "mode=$MODE"
  echo "source=$RESOLVED_SOURCE"
  if [[ "$MODE" == "bin" ]]; then
    echo "tool_bin=$RESOLVED_TOOL_BIN"
  elif [[ "$MODE" == "docker" ]]; then
    echo "docker_cmd=${TOOL_CMD[*]}"
  elif [[ "$MODE" == "apptainer" ]]; then
    echo "apptainer_cmd=${TOOL_CMD[*]}"
  fi
  echo "auto_find=$AUTO_FIND"
  if [[ "${#SEARCH_ROOTS[@]}" -gt 0 ]]; then
    echo "search_roots=${SEARCH_ROOTS[*]}"
  fi
} >"$RESOLVE_FILE"

TOOL_LOG="$LOGS_DIR/idelta_gridselect.log"

if [[ "$MODE" == "bin" ]]; then
  TOOL_MOLECULES="$MOLS_BBOX"
  TOOL_ROI="$ROI_NORM"
  TOOL_OUT="$TOOL_OUT_DIR"
  TOOL_TRANSCRIPTS="$TRANSCRIPTS"
else
  TOOL_MOLECULES="/out/staging/molecules_bbox.tsv"
  TOOL_ROI="/out/staging/roi_vertices.csv"
  TOOL_OUT="/out/tool"
  TOOL_TRANSCRIPTS="/xenium/$(basename "$TRANSCRIPTS")"
fi

wlog "step=resolve_tool seconds=$(( $(date +%s) - T0 ))"

## 1) Normalize ROI vertices -> staging/roi_vertices.csv
T1="$(date +%s)"
if [[ ! -s "$ROI_NORM" || "$ROI_NORM" -ot "$ROI_CSV" ]]; then
  echo "[gridselect] Normalizing ROI vertices -> $ROI_NORM"
  python3 - "$ROI_CSV" "$ROI_NORM" <<'PY'
import csv, io, os, sys
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

def _find_idx(keys):
    for k in keys:
        try:
            return hdr_low.index(k)
        except ValueError:
            pass
    return None

idx_x = _find_idx(["x_um", "x"])
idx_y = _find_idx(["y_um", "y"])
if idx_x is None or idx_y is None:
    raise SystemExit(f"ROI header must include X,Y or x_um,y_um; got: {hdr}")

out.parent.mkdir(parents=True, exist_ok=True)
with out.open("w", newline="") as f:
    w = csv.writer(f, delimiter=",", lineterminator="\n")
    w.writerow(["x_um", "y_um"])
    n = 0
    for row in r:
        if not row or len(row) <= max(idx_x, idx_y):
            continue
        x = float(row[idx_x])
        y = float(row[idx_y])
        w.writerow([x, y])
        n += 1
if n < 3:
    raise SystemExit(f"ROI needs >=3 vertices; got {n}: {inp}")
PY
else
  echo "[gridselect] ROI staging is up-to-date: $ROI_NORM"
fi

wlog "step=normalize_roi seconds=$(( $(date +%s) - T1 ))"

## 2) Stream transcripts and emit bbox-prefiltered molecules TSV.
T2="$(date +%s)"
PY_HELPER="$REPO_ROOT/scripts/gridselect_make_molecules_bbox.py"
MOLS_UPTODATE=0
if [[ -s "$MOLS_BBOX" && "$MOLS_BBOX" -nt "$TRANSCRIPTS" && "$MOLS_BBOX" -nt "$ROI_NORM" ]]; then
  if [[ "$TRANSCRIPTS_KIND" == "csv" && "$MOLS_BBOX" -ot "$PY_HELPER" ]]; then
    MOLS_UPTODATE=0
  else
    MOLS_UPTODATE=1
  fi
fi

if [[ "$MOLS_UPTODATE" -eq 1 ]]; then
  echo "[gridselect] Molecules staging is up-to-date: $MOLS_BBOX"
else
  echo "[gridselect] Building bbox-prefiltered molecules -> $MOLS_BBOX"
  rm -f "$MOLS_LOG" 2>/dev/null || true

  if [[ "$TRANSCRIPTS_KIND" == "parquet" ]]; then
    {
      echo "[gridselect] cmd: ${TOOL_CMD[*]} make-molecules-bbox --roi-csv \"$TOOL_ROI\" --transcripts \"$TOOL_TRANSCRIPTS\" --out-tsv \"$TOOL_MOLECULES\""
      ( "${TOOL_CMD[@]}" make-molecules-bbox --roi-csv "$TOOL_ROI" --transcripts "$TOOL_TRANSCRIPTS" --out-tsv "$TOOL_MOLECULES" )
    } >"$MOLS_LOG" 2>&1 || {
      echo "[gridselect] make-molecules-bbox failed. See: $MOLS_LOG" >&2
      tail -n 50 "$MOLS_LOG" >&2 || true
      err "molecules bbox extraction failed; update idelta-gridselect to a version that supports make-molecules-bbox"
    }
  else
    python3 "$PY_HELPER" \
        --roi_csv "$ROI_NORM" \
        --transcripts_gz "$TRANSCRIPTS" \
        --out_tsv "$MOLS_BBOX" 2>&1 | tee "$MOLS_LOG"
  fi
fi

wlog "step=make_molecules_bbox seconds=$(( $(date +%s) - T2 ))"

## 2b) Optional max-width cap via generated widths files (passed to tool).
WIDTHS_PRESCREEN=""
WIDTHS_ANCHOR=""
TOOL_WIDTHS_PRESCREEN=""
TOOL_WIDTHS_ANCHOR=""
if [[ -n "${MAX_WIDTH_UM:-}" ]]; then
  WIDTHS_PRESCREEN="$STAGING_DIR/widths_prescreen_v2_max${MAX_WIDTH_UM}.txt"
  WIDTHS_ANCHOR="$STAGING_DIR/widths_anchor_v2_max${MAX_WIDTH_UM}.txt"
  if [[ "$MODE" == "bin" ]]; then
    TOOL_WIDTHS_PRESCREEN="$WIDTHS_PRESCREEN"
    TOOL_WIDTHS_ANCHOR="$WIDTHS_ANCHOR"
  else
    TOOL_WIDTHS_PRESCREEN="/out/staging/$(basename "$WIDTHS_PRESCREEN")"
    TOOL_WIDTHS_ANCHOR="/out/staging/$(basename "$WIDTHS_ANCHOR")"
  fi
  if [[ ! -s "$WIDTHS_PRESCREEN" ]]; then
    cat >"$WIDTHS_PRESCREEN" <<EOF
# Generated by scripts/run_gridselect_xenium.sh
# base=widths_prescreen_v2 max_width_um=${MAX_WIDTH_UM}
1
2
4
8
16
24
32
40
50
60
80
100
120
150
200
250
300
EOF
    python3 - "$WIDTHS_PRESCREEN" "$MAX_WIDTH_UM" <<'PY'
import sys
path=sys.argv[1]; mx=int(float(sys.argv[2]))
out=[]
for line in open(path):
    s=line.strip()
    if not s or s.startswith("#"): 
        continue
    out.append(int(s))
out=[w for w in out if w<=mx]
if len(out)<3:
    raise SystemExit("max-width-um too small: prescreen widths would have <3 entries")
with open(path,"w") as f:
    f.write("# Generated by scripts/run_gridselect_xenium.sh\n")
    f.write(f"# base=widths_prescreen_v2 max_width_um={mx}\n")
    for w in out:
        f.write(str(w)+"\n")
PY
  fi
  if [[ ! -s "$WIDTHS_ANCHOR" ]]; then
    python3 - "$WIDTHS_ANCHOR" "$MAX_WIDTH_UM" <<'PY'
import sys
path=sys.argv[1]; mx=int(float(sys.argv[2]))
base=list(range(1,41))+[45,50,55,60,65,70,80,90,100,110,120,130,140,150,175,200,225,250,275,300]
out=[w for w in base if w<=mx]
if len(out)<6:
    raise SystemExit("max-width-um too small: anchor widths would have <6 entries (min_points=6)")
with open(path,"w") as f:
    f.write("# Generated by scripts/run_gridselect_xenium.sh\n")
    f.write(f"# base=widths_anchor_v2 max_width_um={mx}\n")
    for w in out:
        f.write(str(w)+"\n")
PY
  fi
fi

T3="$(date +%s)"
if [[ "$SKIP_TOOL" -eq 1 ]]; then
  [[ -s "$TOOL_SUMMARY" ]] || err "--skip-tool was set but missing: $TOOL_SUMMARY"
  echo "[gridselect] Skipping tool (--skip-tool)."
else
  NEED_RERUN_SCHEMA=0
  if [[ -s "$TOOL_SUMMARY" ]]; then
    HDR="$(head -n 1 "$TOOL_SUMMARY" || true)"
    if [[ "$HDR" != *"knee_of_mean_curve_um"* ]]; then
      NEED_RERUN_SCHEMA=1
      echo "[gridselect] Detected legacy knee_summary.tsv schema; forcing tool rerun to refresh outputs."
    fi
  fi

  NEED_RERUN_PARAMS=0
  if [[ "$NEED_RERUN_SCHEMA" -eq 0 && -s "$TOOL_SUMMARY" ]]; then
    python3 - "$TOOL_SUMMARY" "$KNEE_SEARCH_MIN_UM" "$DESCENDING_DELTA_UM" <<'PY' || NEED_RERUN_PARAMS=1
import csv, sys
path = sys.argv[1]
want_min = int(float(sys.argv[2]))
want_delta = int(float(sys.argv[3]))
with open(path, newline="") as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        if row.get("gene") == "__SUMMARY__":
            try:
                got_min = int(float(row.get("knee_search_min_um", "")))
                got_delta = int(float(row.get("descending_delta_um", "")))
            except Exception:
                raise SystemExit(1)
            raise SystemExit(0 if (got_min == want_min and got_delta == want_delta) else 1)
raise SystemExit(1)
PY
    if [[ "$NEED_RERUN_PARAMS" -eq 1 ]]; then
      echo "[gridselect] Detected changed knee-search params; forcing tool rerun."
    fi
  fi

  NEED_RERUN_MAXWIDTH=0
  if [[ -n "${MAX_WIDTH_UM:-}" && "$NEED_RERUN_SCHEMA" -eq 0 && -s "$TOOL_SUMMARY" ]]; then
    python3 - "$TOOL_SUMMARY" "$MAX_WIDTH_UM" <<'PY' || NEED_RERUN_MAXWIDTH=1
import csv, sys
path=sys.argv[1]; mx=int(float(sys.argv[2]))
with open(path,newline="") as f:
    r=csv.DictReader(f,delimiter="\\t")
    for row in r:
        if row.get("gene")=="__SUMMARY__":
            wav=row.get("widths_anchor_version","")
            want=f"file:widths_anchor_v2_max{mx}.txt"
            raise SystemExit(0 if wav==want else 1)
raise SystemExit(1)
PY
    if [[ "$NEED_RERUN_MAXWIDTH" -eq 1 ]]; then
      echo "[gridselect] max-width-um was requested; forcing tool rerun with --widths-* override."
    fi
  fi

  TOOL_UPTODATE=0
  if [[ "$NEED_RERUN_SCHEMA" -eq 0 && "$NEED_RERUN_PARAMS" -eq 0 && "$NEED_RERUN_MAXWIDTH" -eq 0 && -s "$TOOL_SUMMARY" && "$TOOL_SUMMARY" -nt "$MOLS_BBOX" && "$TOOL_SUMMARY" -nt "$ROI_NORM" ]]; then
    TOOL_UPTODATE=1
    if [[ -n "${MAX_WIDTH_UM:-}" && ( "$TOOL_SUMMARY" -ot "$WIDTHS_PRESCREEN" || "$TOOL_SUMMARY" -ot "$WIDTHS_ANCHOR" ) ]]; then
      TOOL_UPTODATE=0
      echo "[gridselect] Width override files are newer than tool summary; forcing tool rerun."
    fi
  fi

  if [[ "$TOOL_UPTODATE" -eq 1 ]]; then
    echo "[gridselect] Tool output up-to-date; skipping tool run."
  else
    echo "[gridselect] Running external tool (log: $TOOL_LOG)"
    # Fixed defaults (deterministic); adjust by editing this script if needed.
    TOOL_ARGS=(
      recommend-grid
      --molecules "$TOOL_MOLECULES"
      --roi "$TOOL_ROI"
      --out "$TOOL_OUT"
      --w0-um 1
      --tile-um 250
      --prescreen-fraction 0.25
      --prescreen-reps 25
      --seed 1
      --knee-window-um 20
      --knee-search-min-um "$KNEE_SEARCH_MIN_UM"
      --descending-delta-um "$DESCENDING_DELTA_UM"
      --threads "$THREADS"
      --b-min 50
      --n-min 100
      --min-points 6
      --min-informative-frac-ge2 0.6
    )
    if [[ -n "${MAX_WIDTH_UM:-}" ]]; then
      TOOL_ARGS+=( --widths-prescreen "$TOOL_WIDTHS_PRESCREEN" --widths-anchor "$TOOL_WIDTHS_ANCHOR" )
    fi
    # Run and capture stdout/stderr to log.
    {
      echo "[gridselect] cmd: ${TOOL_CMD[*]} ${TOOL_ARGS[*]}"
      ( "${TOOL_CMD[@]}" "${TOOL_ARGS[@]}" )
    } >"$TOOL_LOG" 2>&1 || {
      echo "[gridselect] External tool failed. See: $TOOL_LOG" >&2
      tail -n 50 "$TOOL_LOG" >&2 || true
      exit 1
    }
  fi
fi

wlog "step=run_tool seconds=$(( $(date +%s) - T3 ))"

## 4) R plotting + summary (best_grid.r style)
PLOT1="$PLOTS_OUT_DIR/scope_iDelta_mean_curve_knee_stable.png"
PLOT2="$PLOTS_OUT_DIR/scope_iDelta_all_genes_annotated_knee_stable.png"
PLOT3="$PLOTS_OUT_DIR/scope_geneWiseKnee_UIK_violin.png"
PLOT4="$PLOTS_OUT_DIR/scope_geneWiseKnee_UIK_histogram.png"

BEST_GRID_STYLE="1"
if [[ "$SKIP_PLOTS" -eq 1 ]]; then
  [[ -s "$PLOT1" && -s "$PLOT2" && -s "$PLOT3" && -s "$PLOT4" ]] || err "--skip-plots was set but plots are missing under: $PLOTS_OUT_DIR"
  echo "[gridselect] Skipping plot generation (--skip-plots)."
  BEST_GRID_STYLE="0"
elif [[ -s "$PLOT1" && -s "$PLOT2" && -s "$PLOT3" && -s "$PLOT4" && "$PLOT1" -nt "$TOOL_SUMMARY" && "$PLOT2" -nt "$TOOL_SUMMARY" && "$PLOT3" -nt "$TOOL_SUMMARY" && "$PLOT4" -nt "$TOOL_SUMMARY" ]]; then
  echo "[gridselect] Plots up-to-date; skipping plot regeneration."
  BEST_GRID_STYLE="0"
fi

SUMMARY_LOG="$LOGS_DIR/best_grid_summary.txt"
T4="$(date +%s)"
wlog "step=summary begin best_grid_style=$BEST_GRID_STYLE"
echo "[gridselect] Generating summary (best_grid_style=$BEST_GRID_STYLE) (log: $SUMMARY_LOG)"
Rscript "$REPO_ROOT/inst/scripts/recommend_grid_external.R" \
  --out "$OUTDIR" \
  --log_out "$SUMMARY_LOG" \
  --strict "$STRICT" \
  --best_grid_style "$BEST_GRID_STYLE"
T4_DONE="$(date +%s)"
{
  echo "[gridselect] $(date '+%Y-%m-%d %H:%M:%S') step=summary seconds=$((T4_DONE - T4))"
  echo "[gridselect] $(date '+%Y-%m-%d %H:%M:%S') total_seconds=$((T4_DONE - T0))"
} >>"$WRAPPER_LOG"
