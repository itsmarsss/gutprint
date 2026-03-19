import csv

def write_tsv(ranked: list[tuple[str, float]], output_file: str = "results/results.tsv"):
    """Write ranked results to a TSV file."""
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["species", "score"])
        for species, score in ranked:
            writer.writerow([species, f"{score:.6f}"])
    print(f"Results saved to {output_file}")


def print_report(ranked: list[tuple[str, float]], top_n: int = 10):
    """Print a clean summary report to the terminal."""
    print("\n" + "=" * 50)
    print("  gutprint — dietary DNA report")
    print("=" * 50)
    print(f"  {'Species':<28} {'Score':>8}  {'Bar'}")
    print("-" * 50)

    max_score = ranked[0][1] if ranked else 1
    for species, score in ranked[:top_n]:
        bar_len = int((score / max_score) * 20)
        bar = "█" * bar_len
        print(f"  {species:<28} {score:>8.4f}  {bar}")

    print("=" * 50)
    print(f"  {len(ranked)} species detected above threshold")


def plot_chart(ranked: list[tuple[str, float]], top_n: int = 15, output_file: str = "results.html"):
    """Plot a bar chart of species scores using plotly."""
    try:
        import plotly.graph_objects as go
    except ImportError:
        print("plotly not installed — run: pip install plotly")
        return

    top = ranked[:top_n]
    species = [s for s, _ in top]
    scores  = [sc for _, sc in top]

    # shorten species names to just genus + species
    labels = [" ".join(s.split()[:2]) for s in species]

    fig = go.Figure(go.Bar(
        x=scores,
        y=labels,
        orientation="h",
        marker_color="steelblue",
    ))

    fig.update_layout(
        title="gutprint — detected food species",
        xaxis_title="confidence score",
        yaxis_title="species",
        yaxis=dict(autorange="reversed"),
        height=500,
        margin=dict(l=160, r=40, t=60, b=60),
    )

    fig.write_html(output_file)
    print(f"Chart saved to {output_file}")


if __name__ == "__main__":
    from parse_fastq import parse_fastq
    from build_index import build_index
    from score import score_reads, normalize_scores, rank_species
    import os

    print("Loading reads...")
    reads = parse_fastq("data/SRR23930995_1.fastq", min_quality=20)

    print("Building index...")
    index = build_index([
        "markers/markers_animals.tsv",
        "markers/markers_plants.tsv",
    ], k=15)

    print("Scoring...")
    raw    = score_reads(reads, index, k=15)
    norm   = normalize_scores(raw, index, k=15)
    ranked = rank_species(norm)

    os.makedirs("results", exist_ok=True)

    print_report(ranked)
    write_tsv(ranked)
    plot_chart(ranked)