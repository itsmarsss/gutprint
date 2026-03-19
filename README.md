## First look up which SRR numbers are in this project

```bash
esearch -db sra -query PRJNA947193 | efetch -format runinfo | head
```

## Then pull just 50k reads from one sample

```bash
fastq-dump --split-files -X 50000 SRR_NUMBER_HERE -O data/
```

## Example:

```bash
esearch -db sra -query PRJNA947193 | efetch -format runinfo | head

fastq-dump --split-files -X 50000 SRR23930995 -O
```

> Run from Experiement `SRX19741077`
