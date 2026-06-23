# OmicsFlow — Complete Usage Guide

This guide walks you through a complete RNA-seq analysis from raw FASTQ files to differential expression results, step by step.

---

## Prerequisites

You only need two tools installed on your machine:

1. **Docker** — [docs.docker.com/get-docker](https://docs.docker.com/get-docker/)
2. **Nextflow** (optional, for automated pipeline) — [nextflow.io](https://www.nextflow.io/docs/latest/install.html)

Pull the OmicsFlow Docker image:
```bash
docker pull smill/omicsflow:1.0.0
```

---

## Step 0 — Organize your data

```
my_experiment/
├── data/
│   ├── raw/
│   │   ├── ctrl_rep1_R1.fastq.gz
│   │   ├── ctrl_rep1_R2.fastq.gz
│   │   ├── ctrl_rep2_R1.fastq.gz
│   │   ├── ctrl_rep2_R2.fastq.gz
│   │   ├── treat_rep1_R1.fastq.gz
│   │   ├── treat_rep1_R2.fastq.gz
│   │   ├── treat_rep2_R1.fastq.gz
│   │   └── treat_rep2_R2.fastq.gz
│   ├── reference/
│   │   ├── genome.fa          ← reference genome FASTA
│   │   └── genes.gtf          ← gene annotation GTF
│   └── samplesheet.csv
└── results/                   ← pipeline output goes here
```

---

## Step 1 — Create your samplesheet

Create `samplesheet.csv` with your samples:

```csv
sample,fastq_1,fastq_2,strandedness,condition
ctrl_rep1,data/raw/ctrl_rep1_R1.fastq.gz,data/raw/ctrl_rep1_R2.fastq.gz,reverse,control
ctrl_rep2,data/raw/ctrl_rep2_R1.fastq.gz,data/raw/ctrl_rep2_R2.fastq.gz,reverse,control
treat_rep1,data/raw/treat_rep1_R1.fastq.gz,data/raw/treat_rep1_R2.fastq.gz,reverse,treatment
treat_rep2,data/raw/treat_rep2_R1.fastq.gz,data/raw/treat_rep2_R2.fastq.gz,reverse,treatment
```

**Columns:**
- `sample` — unique sample name
- `fastq_1` — path to R1 FASTQ file
- `fastq_2` — path to R2 FASTQ file (leave empty for single-end)
- `strandedness` — `reverse` (TruSeq), `forward`, or `unstranded`
- `condition` — experimental group (used by DESeq2)

---

## Step 2 — Build reference indexes (one-time)

### 2a. STAR genome index

```bash
docker run --rm -v $(pwd):/data smill/omicsflow:1.0.0 \
  bash -c "mkdir -p /data/star_index && STAR \
  --runMode genomeGenerate \
  --genomeDir /data/star_index \
  --genomeFastaFiles /data/reference/genome.fa \
  --sjdbGTFfile /data/reference/genes.gtf \
  --runThreadN 8"
```

> ⚠️ Build once, reuse forever. For human GRCh38: ~45 min, ~30 GB disk.

### 2b. Salmon transcriptome index

```bash
# First, extract transcript sequences from genome + GTF
docker run --rm -v $(pwd):/data smill/omicsflow:1.0.0 \
  bash -c "export LC_ALL=C && salmon index \
  -t /data/reference/transcriptome.fa \
  -i /data/salmon_index \
  --threads 8"
```

### 2c. Generate tx2gene.csv

Required for gene-level summarization in DESeq2:

```bash
docker run --rm -v $(pwd):/data smill/omicsflow:1.0.0 \
  python3 scripts/generate_tx2gene.py \
  --gtf /data/reference/genes.gtf \
  --output /data/tx2gene.csv
```

---

## Step 3 — Quality Control

Check raw read quality:

```bash
docker run --rm -v $(pwd):/data smill/omicsflow:1.0.0 \
  bash -c "fastqc /data/raw/*.fastq.gz --outdir /data/results/qc --threads 4"
```

Open `results/qc/*.html` to inspect quality scores, adapter content, and GC content.

---

## Step 4 — Adapter Trimming

```bash
docker run --rm -v $(pwd):/data smill/omicsflow:1.0.0 \
  bash -c "trim_galore --paired --cores 4 \
  /data/raw/ctrl_rep1_R1.fastq.gz /data/raw/ctrl_rep1_R2.fastq.gz \
  -o /data/trimmed"
```

Repeat for each sample, or use the automated Nextflow pipeline (Step 8).

---

## Step 5 — Alignment with STAR

```bash
docker run --rm -v $(pwd):/data smill/omicsflow:1.0.0 \
  bash -c "mkdir -p /data/aligned/ctrl_rep1 && STAR \
  --runMode alignReads \
  --genomeDir /data/star_index \
  --readFilesIn /data/trimmed/ctrl_rep1_R1_val_1.fq.gz \
                /data/trimmed/ctrl_rep1_R2_val_2.fq.gz \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix /data/aligned/ctrl_rep1/ \
  --runThreadN 4"
```

Check alignment stats:
```bash
docker run --rm -v $(pwd):/data smill/omicsflow:1.0.0 \
  bash -c "cat /data/aligned/ctrl_rep1/Log.final.out"
```

---

## Step 6 — Quantification with Salmon

```bash
docker run --rm -v $(pwd):/data smill/omicsflow:1.0.0 \
  bash -c "export LC_ALL=C && salmon quant \
  --index /data/salmon_index \
  --libType A \
  -1 /data/trimmed/ctrl_rep1_R1_val_1.fq.gz \
  -2 /data/trimmed/ctrl_rep1_R2_val_2.fq.gz \
  --output /data/counts/ctrl_rep1 \
  --threads 4 --validateMappings"
```

---

## Step 7 — Aggregated QC Report (MultiQC)

```bash
docker run --rm -v $(pwd):/data smill/omicsflow:1.0.0 \
  bash -c "export LC_ALL=C.UTF-8 && multiqc /data --outdir /data/results/multiqc"
```

Open `results/multiqc/multiqc_report.html` for a summary of all samples.

---

## Step 8 — Full automated pipeline (Nextflow)

Run everything in one command:

```bash
nextflow run Millimono/OmicsFlow \
  --input samplesheet.csv \
  --star_index data/star_index \
  --salmon_index data/salmon_index \
  --tx2gene data/tx2gene.csv \
  --outdir results/ \
  -profile docker
```

Or clone locally first:
```bash
git clone https://github.com/Millimono/OmicsFlow.git
cd OmicsFlow

nextflow run rnaseq.nf \
  --input data/test/samplesheet.csv \
  --star_index /data/star_index \
  --salmon_index /data/salmon_index \
  --tx2gene /data/tx2gene.csv \
  --outdir results/ \
  -profile docker
```

---

## Step 9 — Interpret your results

After the pipeline completes, check `results/deseq2/`:

| File | Description |
|---|---|
| `deseq2_results.csv` | All genes with LFC, p-value, padj |
| `normalized_counts.csv` | VST-normalized expression matrix |
| `volcano_plot.pdf` | Up/down regulated genes visualization |
| `heatmap_top50.pdf` | Top 50 DE genes heatmap |
| `pca_plot.pdf` | Sample clustering — verify replicate quality |

**Key columns in `deseq2_results.csv`:**
- `log2FoldChange` — positive = up in treatment, negative = down
- `padj` — adjusted p-value (use this, not `pvalue`)
- `gene` — gene ID

Significant DE genes: `padj < 0.05` AND `|log2FoldChange| > 1`

---

## Troubleshooting

| Problem | Solution |
|---|---|
| STAR index error | Check `--genomeSAindexNbases` — use 14 for human, 7 for small genomes |
| Salmon locale error | Add `export LC_ALL=C` before salmon commands |
| DESeq2 condition error | Check `condition` column exists in samplesheet |
| Low alignment rate (<70%) | Check genome matches your species; verify strandedness |
| tx2gene missing | Run `generate_tx2gene.py` script — see Step 2c |
| Docker permission error | Add `--rm` flag; check volume mount path |

---

## Getting help

- 📧 millimono64.sm@gmail.com
- 🔗 [GitHub Issues](https://github.com/Millimono/OmicsFlow/issues)
- 📚 [Full documentation](https://millimono.github.io/OmicsFlow)
