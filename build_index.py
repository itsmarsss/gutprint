import argparse


def load_markers(filepath: str) -> list[dict]:
    """
    Load marks TSV into a list of dicts.
    Each dict has keys: species, marker_type, header, sequence.
    """
    markers = []
    with open(filepath, "r") as f:
        header = f.readline() # skip header row
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 4: # invalid line
                continue
            markers.append({
                "species": parts[0],
                "marker_type": parts[1],
                "header": parts[2],
                "sequence": parts[3]
            })
    return markers

def build_index(marker_files: list[str], k: int = 15) -> dict:
    """
    Build a k-mer index from one or more index TSV files.

    Returns a dict mapping each k-mer to the set of species that contain it:
        { "<sequence>": {"<species1>", "<species2>" }, ... }
    """
    index = {}

    for filepath in marker_files:
        markers = load_markers(filepath)
        for marker in markers:
            species = marker["species"]
            sequence = marker["sequence"]
            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i+k]
                if kmer not in index:
                    index[kmer] = set()
                index[kmer].add(species)
    return index

def index_stats(index: dict) -> None:
    """Print some basic stats about the index."""
    total_kmers = len(index)
    unique_kmers = sum(1 for species in index.values() if len(species) == 1)
    shared_kmers = total_kmers - unique_kmers

    print(f"Total k-mers in index:    {total_kmers:,}")
    print(f"Species-unique k-mers:    {unique_kmers:,}")
    print(f"Shared k-mers:            {shared_kmers:,}")
    print(f"Unique k-mer ratio:       {unique_kmers/total_kmers:.1%}")

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build and inspect a marker k-mer index.")
    parser.add_argument(
        "--marker-files",
        nargs="+",
        default=[
            "markers/markers_animals.tsv",
            "markers/markers_plants.tsv",
        ],
        help="One or more marker TSV files.",
    )
    parser.add_argument(
        "-k",
        type=int,
        default=15,
        help="k-mer size to index.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    index = build_index(args.marker_files, k=args.k)
    index_stats(index)
