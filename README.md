# gutprint

gutprint is a lightweight Python pipeline and web app for estimating likely dietary species from stool metagenomic FASTQ reads using marker k-mer matching.

## What This Project Does

The pipeline:

1. Parses FASTQ reads and filters by minimum average PHRED quality.
2. Builds a species marker k-mer index from TSV marker files.
3. Scores species by weighted k-mer overlap across reads.
4. Normalizes and ranks species confidence scores.
5. Exports ranked results as TSV and optional chart HTML.

This repository includes:

- CLI scripts for each stage.
- A one-command pipeline runner script.
- A Flask + HTMX web interface that can use local FASTQ, uploaded FASTQ, or NCBI SRA run accessions.

## Repository Layout

- `parse_fastq.py`: FASTQ parsing and quality filtering.
- `build_index.py`: marker loading and k-mer index construction.
- `score.py`: weighted scoring, normalization, and ranking.
- `report.py`: report printing, TSV output, and optional chart output.
- `fetch_markers.py`: NCBI marker retrieval for configured species lists.
- `gutprint.sh`: one-command pipeline execution script.
- `app.py`: Flask web application.

## Requirements

- Python 3.10+ recommended.
- pip for dependency installation.
- SRA Toolkit (`fasterq-dump` or `fastq-dump`) if using NCBI SRA FASTQ download.
- Optional: `plotly` for chart export.

Install web app dependencies:

```bash
pip install -r requirements.txt
```

Optional chart dependency:

```bash
pip install plotly
```

## Quick Start CLI

Run full pipeline in one command:

```bash
./gutprint.sh --fastq data/SRR23930995_1.fastq
```

See all options:

```bash
./gutprint.sh --help
```

Example with custom outputs:

```bash
./gutprint.sh \
  --fastq data/SRR23930994_1.fastq \
  --min-quality 20 \
  --k 15 \
  --min-score 0.001 \
  --top-n 10 \
  --results-tsv results/results.tsv \
  --chart-output results/results.html
```

## Web App

Start server:

```bash
python app.py
```

Open:

- `http://127.0.0.1:5000`

Web UI FASTQ source precedence:

1. Uploaded FASTQ file
2. NCBI SRA accession (example `SRR23930994`)
3. Selected existing local FASTQ

## NCBI SRA Download Notes

The web app can download FASTQ from SRA runs when an accession is provided.

- Valid pattern: `SRR` followed by digits.
- Read limit is controlled by `NCBI max reads` in the form.
- The app uses `fasterq-dump` first and falls back to `fastq-dump`.

You can also inspect run IDs from a BioProject using E-utilities:

```bash
esearch -db sra -query PRJNA947193 | efetch -format runinfo | head
```

## Markers

Expected marker files:

- `markers/markers_animals.tsv`
- `markers/markers_plants.tsv`

Regenerate marker files:

```bash
python fetch_markers.py --email you@example.com
```

## Output

- Ranked TSV outputs are written under `results/`.
- Chart HTML output is generated when `plotly` is installed.
- Web app runs write timestamped result files.

## Limitations

- This is a research and development prototype.
- Output scores are heuristic and depend on marker quality, read depth, and parameter settings.
- This project is not a clinical diagnostic tool.

## References and Credit

This project is inspired by and builds on the MEDI method:

Diener, C., Holscher, H. D., Filek, K., Corbin, K. D., Moissl-Eichinger, C., & Gibbons, S. M. (2025). Metagenomic estimation of dietary intake from human stool. Nature Metabolism, 7(3), 617-630. https://doi.org/10.1038/s42255-025-01220-1

PubMed record:

- https://pubmed.ncbi.nlm.nih.gov/39966520/
