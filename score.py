import argparse


def score_reads(reads: list[str], index: dict, k: int = 15) -> dict[str, float]:
    """
    Score each species by how many k-mers from the reads match its markers.
    
    Uses weighted scoring — k-mers shared across many species contribute less.
    Returns a dict mapping species → confidence score.
    """
    scores = {}

    for read in reads:
        # extract every k-mer from this read
        for i in range(len(read) - k + 1):
            kmer = read[i:i+k]

            if kmer not in index:
                continue

            matched_species = index[kmer]
            weight = 1 / len(matched_species)  # 1 if unique, 0.5 if shared by 2, etc.

            for species in matched_species:
                if species not in scores:
                    scores[species] = 0.0
                scores[species] += weight

    return scores

def normalize_scores(scores: dict[str, float], index: dict, k: int = 15) -> dict[str, float]:
    """
    Normalize raw scores by the total possible k-mers per species.
    This makes scores comparable across species with different marker lengths.
    
    Score of 1.0 = every possible k-mer for that species was found in the reads.
    Score of 0.01 = 1% of k-mers found — weak signal.
    """
    # count how many k-mers in the index belong to each species
    possible_kmers = {}
    for kmer, species_set in index.items():
        for species in species_set:
            if species not in possible_kmers:
                possible_kmers[species] = 0
            possible_kmers[species] += 1

    normalized = {}
    for species, raw_score in scores.items():
        total = possible_kmers.get(species, 1)
        normalized[species] = raw_score / total

    return normalized

def rank_species(normalized_scores: dict[str, float], min_score: float = 0.001) -> list[tuple[str, float]]:
    """
    Sort species by score descending, filter out noise below min_score.
    Returns list of (species, score) tuples.
    """
    ranked = [
        (species, score)
        for species, score in normalized_scores.items()
        if score >= min_score
    ]
    return sorted(ranked, key=lambda x: x[1], reverse=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Score species from FASTQ reads.")
    parser.add_argument(
        "--fastq",
        default="data/SRR23930994_1.fastq",
        help="Path to input FASTQ file.",
    )
    parser.add_argument(
        "--min-quality",
        type=int,
        default=20,
        help="Minimum average PHRED quality to keep a read.",
    )
    parser.add_argument(
        "--marker-files",
        nargs="+",
        default=[
            "markers/markers_animals.tsv",
            "markers/markers_plants.tsv",
        ],
        help="One or more marker TSV files.",
    )
    parser.add_argument("-k", type=int, default=15, help="k-mer size.")
    parser.add_argument(
        "--min-score",
        type=float,
        default=0.001,
        help="Minimum normalized score to include in ranked output.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    from parse_fastq import parse_fastq
    from build_index import build_index

    args = parse_args()

    print("Loading reads...")
    reads = parse_fastq(args.fastq, min_quality=args.min_quality)

    print("Building index...")
    index = build_index(args.marker_files, k=args.k)

    print("Scoring...")
    raw_scores  = score_reads(reads, index, k=args.k)
    norm_scores = normalize_scores(raw_scores, index, k=args.k)
    ranked      = rank_species(norm_scores, min_score=args.min_score)

    print("\nResults:")
    print(f"{'Species':<30} {'Score':>10}")
    print("-" * 42)
    for species, score in ranked:
        print(f"{species:<30} {score:>10.6f}")

    # how many reads matched anything at all?
    total_reads = len(reads)
    matched_reads = sum(
        1 for read in reads
        for i in range(len(read) - args.k + 1)
        if read[i:i+args.k] in index
    )
    print(f"\nReads with at least one k-mer match: ~{matched_reads}")
    print(f"Match rate: {matched_reads/total_reads:.2%}" if total_reads else "Match rate: n/a")
