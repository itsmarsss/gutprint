## First look up which SRR numbers are in this project

```bash
esearch -db sra -query PRJNA947193 | efetch -format runinfo | head
```

## Then pull just 500k reads from one sample

```bash
fastq-dump --split-files -X 500000 SRR_NUMBER_HERE -O data/
```

## Example:

```bash
esearch -db sra -query PRJNA947193 | efetch -format runinfo | head

fastq-dump --split-files -X 500000 SRR23930994 -O
```

> Run from Experiement `SRX19741078`
