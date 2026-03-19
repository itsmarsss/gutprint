import urllib.request
import urllib.parse
import time
import json

BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
EMAIL = "your@email.com"

# Each marker type has its own NCBI search terms
MARKER_QUERIES = {
    "CO1":  "{species}[Organism] COX1[Gene] mitochondrion",
    "ITS":  "{species}[Organism] internal transcribed spacer",
    "16S":  "{species}[Organism] 16S ribosomal RNA",
    "18S":  "{species}[Organism] 18S ribosomal RNA",
}

def search_ncbi(query: str, max_results: int = 1) -> list[str]:
    params = urllib.parse.urlencode({
        "db": "nucleotide",
        "term": query,
        "retmax": max_results,
        "retmode": "json",
        "email": EMAIL,
    })
    url = BASE_URL + "esearch.fcgi?" + params
    with urllib.request.urlopen(url) as response:
        data = json.loads(response.read())
    return data["esearchresult"]["idlist"]

def fetch_sequences(ids: list[str]) -> str:
    params = urllib.parse.urlencode({
        "db": "nucleotide",
        "id": ",".join(ids),
        "rettype": "fasta",
        "retmode": "text",
        "email": EMAIL,
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
    output_file: str = "markers/markers.tsv"
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
                ids = search_ncbi(query, max_results=1)

                if not ids:
                    print(f"  No {marker_type} found for {species}, skipping")
                    continue

                fasta = fetch_sequences(ids)
                parsed = parse_fasta(fasta)

                for header, seq in parsed:
                    out.write(f"{species}\t{marker_type}\t{header}\t{seq}\n")

                time.sleep(0.4)

    print(f"Done. Saved to {output_file}")


if __name__ == "__main__":
    import os

    # Animals → CO1 is best
    # Plants → ITS is best
    # Mix and match per species
    animals = [
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

    plants = [
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

    animals_path = "markers/markers_animals.tsv"
    plants_path = "markers/markers_plants.tsv"

    os.makedirs("markers", exist_ok=True)

    fetch_markers(animals, marker_types=["CO1"], output_file=animals_path)
    fetch_markers(plants,  marker_types=["ITS"], output_file=plants_path)
