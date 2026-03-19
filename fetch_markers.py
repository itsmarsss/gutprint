import argparse
import os
import urllib.request
import urllib.parse
import time
import json

BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
DEFAULT_EMAIL = "your@email.com"

# Each marker type has its own NCBI search terms
MARKER_QUERIES = {
    "CO1":  "{species}[Organism] COX1[Gene] mitochondrion",
    "ITS":  "{species}[Organism] internal transcribed spacer",
    "16S":  "{species}[Organism] 16S ribosomal RNA",
    "18S":  "{species}[Organism] 18S ribosomal RNA",
}

def search_ncbi(query: str, email: str, max_results: int = 1) -> list[str]:
    params = urllib.parse.urlencode({
        "db": "nucleotide",
        "term": query,
        "retmax": max_results,
        "retmode": "json",
        "email": email,
    })
    url = BASE_URL + "esearch.fcgi?" + params
    with urllib.request.urlopen(url) as response:
        data = json.loads(response.read())
    return data["esearchresult"]["idlist"]

def fetch_sequences(ids: list[str], email: str) -> str:
    params = urllib.parse.urlencode({
        "db": "nucleotide",
        "id": ",".join(ids),
        "rettype": "fasta",
        "retmode": "text",
        "email": email,
    })
    url = BASE_URL + "efetch.fcgi?" + params
    with urllib.request.urlopen(url) as response:
        return response.read().decode("utf-8")

def parse_fasta(fasta_text: str) -> list[tuple[str, str]]:
    sequences = []
    current_header = None
    current_seq = []
    for line in fasta_text.strip().splitlines():
        if line.startswith(">"):
            if current_header:
                sequences.append((current_header, "".join(current_seq)))
            current_header = line[1:]
            current_seq = []
        else:
            current_seq.append(line.strip().upper())
    if current_header:
        sequences.append((current_header, "".join(current_seq)))
    return sequences

def fetch_markers(
    species_list: list[str],
    marker_types: list[str],
    output_file: str = "markers/markers.tsv",
    max_results: int = 1,
    email: str = DEFAULT_EMAIL,
    delay_seconds: float = 0.4,
):
    """
    Fetch markers for each species x marker type combination.
    Falls back gracefully if a species has no result for a given marker type.
    """
    if invalid := [m for m in marker_types if m not in MARKER_QUERIES]:
        raise ValueError(f"Unknown marker types: {invalid}. Choose from {list(MARKER_QUERIES)}")

    with open(output_file, "w") as out:
        out.write("species\tmarker_type\theader\tsequence\n")

        for species in species_list:
            for marker_type in marker_types:
                print(f"Fetching {marker_type} for {species}...")
                query = MARKER_QUERIES[marker_type].format(species=species)
                ids = search_ncbi(query, email=email, max_results=max_results)

                if not ids:
                    print(f"  No {marker_type} found for {species}, skipping")
                    continue

                fasta = fetch_sequences(ids, email=email)
                parsed = parse_fasta(fasta)

                for header, seq in parsed:
                    out.write(f"{species}\t{marker_type}\t{header}\t{seq}\n")

                time.sleep(delay_seconds)

    print(f"Done. Saved to {output_file}")


def default_animals() -> list[str]:
    return [
        "Salmo salar",          # salmon
        "Gallus gallus",        # chicken
        "Bos taurus",           # beef
        "Sus scrofa",           # pork
        "Thunnus albacares",    # yellowfin tuna
        "Oncorhynchus mykiss",  # rainbow trout
        "Gadus morhua",         # cod
        "Ovis aries",           # lamb
        "Meleagris gallopavo",  # turkey
        "Penaeus vannamei",     # shrimp
    ]


def default_plants() -> list[str]:
    return [
        "Solanum lycopersicum",  # tomato
        "Triticum aestivum",     # wheat
        "Oryza sativa",          # rice
        "Zea mays",              # corn
        "Musa acuminata",        # banana
        "Glycine max",           # soybean
        "Solanum tuberosum",     # potato
        "Brassica oleracea",     # broccoli/cabbage
        "Daucus carota",         # carrot
        "Capsicum annuum",       # pepper
        "Lactuca sativa",        # lettuce
        "Allium cepa",           # onion
        "Lycopersicon esculentum", # tomato variant
        "Arachis hypogaea",      # peanut
        "Citrullus lanatus",     # watermelon
    ]


def _load_species_file(path: str) -> list[str]:
    species = []
    with open(path, "r") as f:
        for line in f:
            name = line.strip()
            if not name or name.startswith("#"):
                continue
            species.append(name)
    return species


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fetch marker sequences from NCBI.")
    parser.add_argument(
        "--species",
        action="append",
        default=[],
        help="Species name to include (repeat for multiple values).",
    )
    parser.add_argument(
        "--species-file",
        action="append",
        default=[],
        help="Path to a file containing one species name per line.",
    )
    parser.add_argument(
        "--marker-types",
        nargs="+",
        default=["CO1"],
        choices=sorted(MARKER_QUERIES),
        help="Marker types to query.",
    )
    parser.add_argument(
        "--output-file",
        default="markers/markers.tsv",
        help="Output TSV file for custom species mode.",
    )
    parser.add_argument(
        "--max-results",
        type=int,
        default=1,
        help="Maximum NCBI IDs to fetch per species x marker query.",
    )
    parser.add_argument(
        "--delay-seconds",
        type=float,
        default=0.4,
        help="Delay between requests in seconds.",
    )
    parser.add_argument(
        "--email",
        default=os.environ.get("NCBI_EMAIL", DEFAULT_EMAIL),
        help="Email value sent to NCBI E-utilities.",
    )
    parser.add_argument(
        "--animals-output",
        default="markers/markers_animals.tsv",
        help="Output file for default animals preset mode.",
    )
    parser.add_argument(
        "--plants-output",
        default="markers/markers_plants.tsv",
        help="Output file for default plants preset mode.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    os.makedirs("markers", exist_ok=True)

    species = list(args.species)
    for species_file in args.species_file:
        species.extend(_load_species_file(species_file))

    if species:
        fetch_markers(
            species,
            marker_types=args.marker_types,
            output_file=args.output_file,
            max_results=args.max_results,
            email=args.email,
            delay_seconds=args.delay_seconds,
        )
    else:
        fetch_markers(
            default_animals(),
            marker_types=["CO1"],
            output_file=args.animals_output,
            max_results=args.max_results,
            email=args.email,
            delay_seconds=args.delay_seconds,
        )
        fetch_markers(
            default_plants(),
            marker_types=["ITS"],
            output_file=args.plants_output,
            max_results=args.max_results,
            email=args.email,
            delay_seconds=args.delay_seconds,
        )
