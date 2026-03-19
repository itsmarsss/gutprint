from __future__ import annotations

import os
import re
import shutil
import subprocess
import threading
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
RUNS: dict[str, dict] = {}
RUNS_LOCK = threading.Lock()


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


def _resolve_fastq_path(params: dict) -> tuple[Path, str]:
    upload_rel = params.get("upload_rel", "").strip()
    chosen = params.get("existing_fastq", "").strip()
    ncbi_accession = params.get("ncbi_accession", "").strip()
    max_reads = int(params.get("ncbi_max_reads", 500000))

    if upload_rel:
        path = ROOT / upload_rel
        if not path.exists():
            raise FileNotFoundError(f"Uploaded FASTQ not found: {upload_rel}")
        return path, upload_rel

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


def _run_init(run_id: str) -> None:
    with RUNS_LOCK:
        RUNS[run_id] = {
            "status": "running",
            "step": "Queued",
            "logs": [],
            "result": None,
            "error": "",
        }


def _run_log(run_id: str, message: str) -> None:
    line = f"[{time.strftime('%H:%M:%S')}] {message}"
    with RUNS_LOCK:
        if run_id in RUNS:
            RUNS[run_id]["logs"].append(line)


def _run_update(run_id: str, **kwargs) -> None:
    with RUNS_LOCK:
        if run_id in RUNS:
            RUNS[run_id].update(kwargs)


def _run_snapshot(run_id: str) -> dict | None:
    with RUNS_LOCK:
        run = RUNS.get(run_id)
        if not run:
            return None
        return {
            "status": run["status"],
            "step": run["step"],
            "logs": list(run["logs"]),
            "result": run["result"],
            "error": run["error"],
        }


def _execute_pipeline_run(run_id: str, params: dict) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    try:
        _run_update(run_id, step="Resolving FASTQ source")
        _run_log(run_id, "Selecting FASTQ source.")
        fastq_path, fastq_display = _resolve_fastq_path(params)
        _run_log(run_id, f"Using FASTQ: {fastq_display}")

        min_quality = int(params.get("min_quality", 20))
        k = int(params.get("k", 15))
        min_score = float(params.get("min_score", 0.001))
        top_n = int(params.get("top_n", 10))
        chart_top_n = int(params.get("chart_top_n", 15))

        marker_files = [str(ROOT / p) for p in DEFAULT_MARKER_FILES]
        missing_markers = [p for p in marker_files if not Path(p).exists()]
        if missing_markers:
            raise FileNotFoundError(f"Missing marker files: {', '.join(missing_markers)}")

        _run_update(run_id, step="Parsing FASTQ reads")
        _run_log(run_id, f"Parsing FASTQ with min_quality={min_quality}.")
        reads = parse_fastq(str(fastq_path), min_quality=min_quality)
        _run_log(run_id, f"Reads kept: {len(reads)}")

        _run_update(run_id, step="Building marker index")
        _run_log(run_id, f"Building k-mer index with k={k}.")
        index_data = build_index(marker_files, k=k)
        _run_log(run_id, f"Indexed k-mers: {len(index_data)}")

        _run_update(run_id, step="Scoring species")
        _run_log(run_id, "Scoring reads against marker index.")
        raw = score_reads(reads, index_data, k=k)
        norm = normalize_scores(raw, index_data, k=k)
        ranked = rank_species(norm, min_score=min_score)
        _run_log(run_id, f"Species above threshold: {len(ranked)}")

        _run_update(run_id, step="Writing outputs")
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
        _run_log(run_id, "Run completed successfully.")
        _run_update(
            run_id,
            status="done",
            step="Completed",
            result={
                "fastq": fastq_display,
                "reads_count": len(reads),
                "species_count": len(ranked),
                "match_reads": match_reads,
                "match_rate": match_rate,
                "top_rows": top_rows,
                "tsv_rel": str(tsv_path.relative_to(ROOT)),
                "chart_rel": str(chart_path.relative_to(ROOT)) if chart_built else "",
                "chart_built": chart_built,
            },
        )
    except Exception as exc:
        _run_log(run_id, f"Run failed: {exc}")
        _run_update(run_id, status="failed", step="Failed", error=str(exc))


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


@app.get("/health")
def health():
    return {"status": "ok"}, 200


@app.post("/start-run")
def start_run():
    uploaded = request.files.get("fastq_upload")
    upload_rel = ""
    if uploaded and uploaded.filename:
        UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
        safe_name = _safe_uploaded_name(uploaded.filename)
        dest = UPLOAD_DIR / f"{int(time.time())}_{uuid.uuid4().hex[:8]}_{safe_name}"
        uploaded.save(dest)
        upload_rel = str(dest.relative_to(ROOT))

    params = {
        "existing_fastq": request.form.get("existing_fastq", "").strip(),
        "ncbi_accession": request.form.get("ncbi_accession", "").strip(),
        "ncbi_max_reads": request.form.get("ncbi_max_reads", "500000"),
        "min_quality": request.form.get("min_quality", "20"),
        "k": request.form.get("k", "15"),
        "min_score": request.form.get("min_score", "0.001"),
        "top_n": request.form.get("top_n", "10"),
        "chart_top_n": request.form.get("chart_top_n", "15"),
        "upload_rel": upload_rel,
    }

    run_id = uuid.uuid4().hex
    _run_init(run_id)
    thread = threading.Thread(target=_execute_pipeline_run, args=(run_id, params), daemon=True)
    thread.start()
    return render_template("_run_status.html", run_id=run_id, run=_run_snapshot(run_id))


@app.get("/run-status/<run_id>")
def run_status(run_id: str):
    run = _run_snapshot(run_id)
    if not run:
        return render_template(
            "_run_status.html",
            run_id=run_id,
            run={"status": "failed", "step": "Missing", "logs": [], "result": None, "error": "Run not found."},
        ), 404
    return render_template("_run_status.html", run_id=run_id, run=run)


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
