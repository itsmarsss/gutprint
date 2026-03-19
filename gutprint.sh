#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python}"

usage() {
  cat <<EOF
Usage:
  ./gutprint.sh [options]

Runs the gutprint pipeline in one command:
  1) Optionally fetch/refresh marker files
  2) Build index stats
  3) Generate report outputs

Options:
  --fastq PATH            Input FASTQ file (default: data/SRR23930995_1.fastq)
  --min-quality N         Minimum average PHRED quality (default: 20)
  --k N                   k-mer size (default: 15)
  --min-score FLOAT       Minimum normalized score to report (default: 0.001)
  --top-n N               Number of species in terminal report (default: 10)
  --chart-top-n N         Number of species in chart (default: 15)
  --results-tsv PATH      TSV output (default: results/results.tsv)
  --chart-output PATH     HTML chart output (default: results.html)
  --no-chart              Skip chart generation

  --animals-output PATH   Animals marker TSV (default: markers/markers_animals.tsv)
  --plants-output PATH    Plants marker TSV (default: markers/markers_plants.tsv)
  --refresh-markers       Force marker re-fetch from NCBI before report
  --email EMAIL           Email for NCBI fetch (default: \$NCBI_EMAIL or your@email.com)
  --max-results N         Max NCBI IDs per query during fetch (default: 1)
  --delay-seconds FLOAT   Delay between fetch requests (default: 0.4)

  -h, --help              Show this help

Examples:
  ./gutprint.sh --fastq data/SRR23930995_1.fastq
  ./gutprint.sh --refresh-markers --email you@example.com --fastq data/SRR23930994_1.fastq
EOF
}

if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  echo "Error: Python interpreter '$PYTHON_BIN' not found."
  echo "Set a different interpreter with: PYTHON_BIN=python3 ./gutprint.sh ..."
  exit 1
fi

fastq="data/SRR23930995_1.fastq"
min_quality="20"
kmer="15"
min_score="0.001"
top_n="10"
chart_top_n="15"
results_tsv="results/results.tsv"
chart_output="results.html"
no_chart="0"

animals_output="markers/markers_animals.tsv"
plants_output="markers/markers_plants.tsv"
refresh_markers="0"
email="${NCBI_EMAIL:-your@email.com}"
max_results="1"
delay_seconds="0.4"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --fastq) fastq="$2"; shift 2 ;;
    --min-quality) min_quality="$2"; shift 2 ;;
    --k) kmer="$2"; shift 2 ;;
    --min-score) min_score="$2"; shift 2 ;;
    --top-n) top_n="$2"; shift 2 ;;
    --chart-top-n) chart_top_n="$2"; shift 2 ;;
    --results-tsv) results_tsv="$2"; shift 2 ;;
    --chart-output) chart_output="$2"; shift 2 ;;
    --no-chart) no_chart="1"; shift ;;
    --animals-output) animals_output="$2"; shift 2 ;;
    --plants-output) plants_output="$2"; shift 2 ;;
    --refresh-markers) refresh_markers="1"; shift ;;
    --email) email="$2"; shift 2 ;;
    --max-results) max_results="$2"; shift 2 ;;
    --delay-seconds) delay_seconds="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *)
      echo "Unknown option: $1"
      usage
      exit 1
      ;;
  esac
done

mkdir -p "$ROOT_DIR/markers" "$ROOT_DIR/results"

if [[ "$refresh_markers" == "1" || ! -f "$ROOT_DIR/$animals_output" || ! -f "$ROOT_DIR/$plants_output" ]]; then
  echo "Step 1/3: fetching markers..."
  "$PYTHON_BIN" "$ROOT_DIR/fetch_markers.py" \
    --animals-output "$animals_output" \
    --plants-output "$plants_output" \
    --email "$email" \
    --max-results "$max_results" \
    --delay-seconds "$delay_seconds"
else
  echo "Step 1/3: using existing marker files."
fi

echo "Step 2/3: building index stats..."
"$PYTHON_BIN" "$ROOT_DIR/build_index.py" \
  --marker-files "$animals_output" "$plants_output" \
  -k "$kmer"

echo "Step 3/3: generating report..."
report_cmd=(
  "$PYTHON_BIN" "$ROOT_DIR/report.py"
  --fastq "$fastq"
  --min-quality "$min_quality"
  --marker-files "$animals_output" "$plants_output"
  -k "$kmer"
  --min-score "$min_score"
  --top-n "$top_n"
  --chart-top-n "$chart_top_n"
  --results-tsv "$results_tsv"
  --chart-output "$chart_output"
)
if [[ "$no_chart" == "1" ]]; then
  report_cmd+=(--no-chart)
fi
"${report_cmd[@]}"

echo "Done."
