def parse_fastq(filepath: str, min_quality: int = 20) -> list[str]:
    """
    Parse a FASTQ file and return a list of high-quality read sequences.
    Filters out reads where the average quality score is below min_quality.
    """
    reads = []
    with open(filepath, "r") as f:
        while True:
            header  = f.readline().strip()  # line 1: @read_id
            seq     = f.readline().strip()  # line 2: sequence
            plus    = f.readline().strip()  # line 3: + (ignored)
            quality = f.readline().strip()  # line 4: quality scores (! to ~ ASCII)

            if not header:  # end of file
                break

            if average_quality(quality) >= min_quality:
                reads.append(seq.upper())

    print(f"Loaded {len(reads)} reads from {filepath}")
    return reads

def average_quality(quality_str: str) -> float:
    """
    Convert ASCII quality string to average Phred quality score.
    Phred score = ASCII value - 33.

    Character  →  ASCII  →  minus 33  →  Phred score  →  accuracy
    '!'        →  33     →  0         →  0            →  50%
    '5'        →  53     →  20        →  20           →  99%
    'I'        →  73     →  40        →  40           →  99.99%

    A score of 20 means 99% base call accuracy.
    A score of 30 means 99.9% base call accuracy.
    """
    scores = [ord(char) - 33 for char in quality_str]
    return sum(scores) / len(scores)

if __name__ == "__main__":
    reads = parse_fastq("data/SRR23930995_1.fastq", min_quality=20)
    print(f"Total reads loaded: {len(reads)}")
    print(f"First read: {reads[0]}")
    print(f"Read length: {len(reads[0])}bp")
