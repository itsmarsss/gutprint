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

if __name__ == "__main__":
    from parse_fastq import parse_fastq
    from build_index import build_index

    print("Loading reads...")
    reads = parse_fastq("data/SRR23930994_1.fastq", min_quality=20)

    print("Building index...")
    index = build_index([
        "markers/markers_animals.tsv",
        "markers/markers_plants.tsv",
    ], k=15)

    print("Scoring...")
    raw_scores  = score_reads(reads, index, k=15)
    norm_scores = normalize_scores(raw_scores, index, k=15)
    ranked      = rank_species(norm_scores)

    print("\nResults:")
    print(f"{'Species':<30} {'Score':>10}")
    print("-" * 42)
    for species, score in ranked:
        print(f"{species:<30} {score:>10.6f}")

    # how many reads matched anything at all?
    total_reads = len(reads)
    matched_reads = sum(
        1 for read in reads
        for i in range(len(read) - 15 + 1)
        if read[i:i+15] in index
    )
    print(f"\nReads with at least one k-mer match: ~{matched_reads}")
    print(f"Match rate: {matched_reads/total_reads:.2%}")