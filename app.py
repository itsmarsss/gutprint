from __future__ import annotations

import os
import re
import shutil
import subprocess
import time
import uuid
from pathlib import Path

from flask import Flask, render_template, request

from build_index import build_index
from parse_fastq import parse_fastq
from report import write_tsv
from score import normalize_scores, rank_species, score_reads

ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
MARKERS_DIR = ROOT / "markers"
RESULTS_DIR = ROOT / "results"
UPLOAD_DIR = DATA_DIR / "uploads"
DOWNLOAD_DIR = DATA_DIR / "downloads"

DEFAULT_FASTQ = "data/SRR23930995_1.fastq"
DEFAULT_MARKER_FILES = [
    "markers/markers_animals.tsv",
    "markers/markers_plants.tsv",
]

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 2 * 1024 * 1024 * 1024  # 2GB upload limit


def _list_fastq_files() -> list[str]:
    fastq_files: list[str] = []
    for ext in ("*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"):
        for p in sorted(DATA_DIR.rglob(ext)):
            fastq_files.append(str(p.relative_to(ROOT)))
    if DEFAULT_FASTQ not in fastq_files and (ROOT / DEFAULT_FASTQ).exists():
        fastq_files.insert(0, DEFAULT_FASTQ)
    return sorted(set(fastq_files))


def _safe_uploaded_name(filename: str) -> str:
    base = os.path.basename(filename).replace(" ", "_")
    return "".join(ch for ch in base if ch.isalnum() or ch in ("-", "_", ".", "+")) or "upload.fastq"


def _download_sra_fastq(accession: str, max_reads: int) -> tuple[Path, str]:
    accession = accession.strip().upper()
    if not re.fullmatch(r"SRR\d+", accession):
        raise ValueError("NCBI accession must look like SRR12345678.")

    DOWNLOAD_DIR.mkdir(parents=True, exist_ok=True)
    r1 = DOWNLOAD_DIR / f"{accession}_1.fastq"
    single = DOWNLOAD_DIR / f"{accession}.fastq"
    if r1.exists():
        return r1, str(r1.relative_to(ROOT))
    if single.exists():
        return single, str(single.relative_to(ROOT))

    tool = shutil.which("fasterq-dump") or shutil.which("fastq-dump")
    if not tool:
        raise RuntimeError("SRA toolkit missing. Install `fasterq-dump` or `fastq-dump`.")

    cmd = [tool, "--split-files"]
    if max_reads > 0:
        cmd.extend(["-X", str(max_reads)])
    cmd.extend([accession, "-O", str(DOWNLOAD_DIR)])
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        stderr = (proc.stderr or "").strip()
        stdout = (proc.stdout or "").strip()
        details = stderr or stdout or "unknown error"
        raise RuntimeError(f"Failed to download {accession}: {details}")

    if r1.exists():
        return r1, str(r1.relative_to(ROOT))
    if single.exists():
        return single, str(single.relative_to(ROOT))
    raise RuntimeError(f"Download completed but FASTQ output not found for {accession}.")


def _selected_fastq_path() -> tuple[Path, str]:
    chosen = request.form.get("existing_fastq", "").strip()
    ncbi_accession = request.form.get("ncbi_accession", "").strip()
    max_reads = int(request.form.get("ncbi_max_reads", "500000"))
    uploaded = request.files.get("fastq_upload")

    if uploaded and uploaded.filename:
        UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
        safe_name = _safe_uploaded_name(uploaded.filename)
        dest = UPLOAD_DIR / f"{int(time.time())}_{uuid.uuid4().hex[:8]}_{safe_name}"
        uploaded.save(dest)
        return dest, str(dest.relative_to(ROOT))

    if ncbi_accession:
        return _download_sra_fastq(ncbi_accession, max_reads=max_reads)

    if not chosen:
        chosen = DEFAULT_FASTQ

    path = ROOT / chosen
    if not path.exists():
        raise FileNotFoundError(f"FASTQ file not found: {chosen}")
    return path, chosen


def _maybe_make_chart(ranked: list[tuple[str, float]], top_n: int, out_path: Path) -> bool:
    try:
        import plotly.graph_objects as go
    except ImportError:
        return False

    top = ranked[:top_n]
    species = [s for s, _ in top]
    scores = [sc for _, sc in top]
    labels = [" ".join(s.split()[:2]) for s in species]

    fig = go.Figure(
        go.Bar(
            x=scores,
            y=labels,
            orientation="h",
            marker_color="#2a7f62",
        )
    )
    fig.update_layout(
        title="gutprint detected species",
        xaxis_title="confidence score",
        yaxis_title="species",
        yaxis=dict(autorange="reversed"),
        height=500,
        margin=dict(l=160, r=40, t=60, b=60),
    )
    fig.write_html(str(out_path))
    return True


@app.get("/")
def index():
    return render_template(
        "index.html",
        fastq_files=_list_fastq_files(),
        defaults={
            "min_quality": 20,
            "k": 15,
            "min_score": 0.001,
            "top_n": 10,
            "chart_top_n": 15,
            "ncbi_max_reads": 500000,
        },
    )


@app.post("/run")
def run_pipeline():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    try:
        fastq_path, fastq_display = _selected_fastq_path()
        min_quality = int(request.form.get("min_quality", "20"))
        k = int(request.form.get("k", "15"))
        min_score = float(request.form.get("min_score", "0.001"))
        top_n = int(request.form.get("top_n", "10"))
        chart_top_n = int(request.form.get("chart_top_n", "15"))

        marker_files = [str(ROOT / p) for p in DEFAULT_MARKER_FILES]
        missing_markers = [p for p in marker_files if not Path(p).exists()]
        if missing_markers:
            raise FileNotFoundError(f"Missing marker files: {', '.join(missing_markers)}")

        reads = parse_fastq(str(fastq_path), min_quality=min_quality)
        index_data = build_index(marker_files, k=k)
        raw = score_reads(reads, index_data, k=k)
        norm = normalize_scores(raw, index_data, k=k)
        ranked = rank_species(norm, min_score=min_score)

        stamp = f"{int(time.time())}_{uuid.uuid4().hex[:6]}"
        tsv_path = RESULTS_DIR / f"web_results_{stamp}.tsv"
        chart_path = RESULTS_DIR / f"web_chart_{stamp}.html"

        write_tsv(ranked, output_file=str(tsv_path))
        chart_built = _maybe_make_chart(ranked, top_n=chart_top_n, out_path=chart_path)

        top_rows = ranked[:top_n]
        match_reads = sum(
            1
            for read in reads
            for i in range(len(read) - k + 1)
            if read[i : i + k] in index_data
        )
        match_rate = (match_reads / len(reads)) if reads else 0.0

        return render_template(
            "_results.html",
            ok=True,
            error="",
            fastq=fastq_display,
            reads_count=len(reads),
            species_count=len(ranked),
            match_reads=match_reads,
            match_rate=match_rate,
            top_rows=top_rows,
            tsv_rel=str(tsv_path.relative_to(ROOT)),
            chart_rel=str(chart_path.relative_to(ROOT)) if chart_built else "",
            chart_built=chart_built,
        )
    except Exception as exc:
        return render_template(
            "_results.html",
            ok=False,
            error=str(exc),
            fastq="",
            reads_count=0,
            species_count=0,
            match_reads=0,
            match_rate=0.0,
            top_rows=[],
            tsv_rel="",
            chart_rel="",
            chart_built=False,
        ), 400


@app.get("/download")
def download_file():
    rel = request.args.get("path", "")
    if not rel:
        return "Missing path", 400
    target = (ROOT / rel).resolve()
    if ROOT.resolve() not in target.parents and target != ROOT.resolve():
        return "Invalid path", 400
    if not target.exists():
        return "File not found", 404

    # Avoid adding extra Flask dependency import in global scope.
    from flask import send_file

    return send_file(target, as_attachment=True)


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5000, debug=True)
