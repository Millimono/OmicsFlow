# 🧬 OmicsFlow

> **A modular, containerized NGS pipeline for RNA-seq, long-read, and metagenomic analysis**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A522.10-brightgreen)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/docker-ready-blue)](https://hub.docker.com/r/smill/omicsflow)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![CI](https://github.com/Millimono/OmicsFlow/actions/workflows/ci.yml/badge.svg)](https://github.com/Millimono/OmicsFlow/actions)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20677900.svg)](https://doi.org/10.5281/zenodo.20677900)

---

## 📋 Overview

**OmicsFlow** is a production-ready bioinformatics pipeline built with [Nextflow](https://www.nextflow.io/) and Docker, designed for reproducible multi-omics data analysis. It supports three major sequencing technologies:

| Workflow | Technology | Key tools | Status |
|---|---|---|---|
| `rnaseq.nf` | Illumina short reads | FastQC · STAR · Salmon · DESeq2 | ✅ Stable |
| `longread.nf` | Oxford Nanopore (ONT) | NanoStat · Minimap2 · Samtools | 🚧 In development |
| `metagenomics.nf` | Illumina / ONT | Kraken2 · Bracken | 🚧 In development |

All workflows are fully containerized via Docker and can run locally, on HPC clusters (SLURM/PBS), or in the cloud (AWS Batch).

---

## 📸 Pipeline Results

### Pipeline Architecture
![Pipeline Architecture](docs/images/pipeline_architecture.png)

### Differential Expression — Volcano Plot
![Volcano Plot](docs/images/volcano_plot.png)

### Top 50 DE Genes — Heatmap
![Heatmap](docs/images/heatmap_top50.png)

### Sample Clustering — PCA
![PCA Plot](docs/images/pca_plot.png)

---

## 📊 Validation Metrics

Benchmarked on nf-core test dataset (S. cerevisiae, GSE110004, 4 samples × 50,000 reads):

| Metric | Value |
|---|---|
| Input reads per sample | 50,000 |
| Reads passing QC | 99.5% |
| Adapter contamination (auto-detected & removed) | 40.3% |
| Uniquely mapped reads (STAR) | 81.8% – 84.6% |
| Properly paired reads | 100% |
| Mismatch rate | 0.9% |
| Pipeline execution time (4 samples, 4 CPUs) | ~8 min |
| Docker image size | 4.63 GB |

---

## 🚀 Quick Start — Run in 3 steps

### Step 1 — Pull the Docker image
```bash
docker pull smill/omicsflow:1.0.0
```

### Step 2 — Prepare your samplesheet
Create a CSV file with your samples:
```csv
sample,fastq_1,fastq_2,strandedness,condition
ctrl_rep1,/data/ctrl_rep1_R1.fastq.gz,/data/ctrl_rep1_R2.fastq.gz,reverse,control
ctrl_rep2,/data/ctrl_rep2_R1.fastq.gz,/data/ctrl_rep2_R2.fastq.gz,reverse,control
treat_rep1,/data/treat_rep1_R1.fastq.gz,/data/treat_rep1_R2.fastq.gz,reverse,treatment
treat_rep2,/data/treat_rep2_R1.fastq.gz,/data/treat_rep2_R2.fastq.gz,reverse,treatment
```

> **Strandedness:** `reverse` for most Illumina TruSeq kits · `forward` for some stranded protocols · `unstranded` if unsure

### Step 3 — Run the pipeline
```bash
nextflow run Millimono/OmicsFlow \
  --input samplesheet.csv \
  --star_index /data/star_index \
  --salmon_index /data/salmon_index \
  --outdir results/ \
  -profile docker
```

> **Windows users:** replace `$(pwd)` with `%cd%` in CMD

> **Need to build reference indexes first?** Use the one-command script:
> ```bash
> docker run --rm -v $(pwd):/data smill/omicsflow:1.0.0 \
>   bash /data/scripts/prepare_references.sh \
>   --genome /data/genome.fa \
>   --gtf    /data/genes.gtf \
>   --outdir /data/references/
> # → builds star_index/ + salmon_index/ + tx2gene.csv automatically
> ```
> Then run the pipeline with `--star_index references/star_index --salmon_index references/salmon_index --tx2gene references/tx2gene.csv`

---

## 📦 What you need before starting

OmicsFlow is flexible — you can use the full pipeline or individual tools depending on your needs.

### Requirements per use case

| Use case | What you need |
|---|---|
| Quality control only | FASTQ files |
| Trimming only | FASTQ files |
| Alignment (STAR) | FASTQ files + STAR index |
| Quantification (Salmon) | FASTQ files + Salmon index |
| Full RNA-seq (FASTQ → DE genes) | FASTQ + STAR index + Salmon index + tx2gene.csv |
| Statistics on existing BAM | BAM file |
| Python / R analysis | Your own data + scripts |

### A. Reference genome & STAR index

If you already have a STAR index on your server — point to it with `--star_index`. No need to rebuild.

If you need to build one:
```bash
# 1. Download reference genome (example: human GRCh38)
wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# 2. Download gene annotation
wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
gunzip Homo_sapiens.GRCh38.109.gtf.gz

# 3. Build STAR index using OmicsFlow Docker
docker run --rm -v $(pwd):/data smill/omicsflow:1.0.0 \
  bash -c "mkdir -p /data/star_index && STAR --runMode genomeGenerate \
  --genomeDir /data/star_index \
  --genomeFastaFiles /data/genome/GRCh38.fa \
  --sjdbGTFfile /data/genome/GRCh38.gtf \
  --runThreadN 8"
```
> ⚠️ Build the STAR index **once per genome**, then reuse it for all experiments.
> Any STAR-compatible index works — it does not need to be built with this Docker image.

### B. Salmon index

```bash
# 1. Download transcriptome FASTA
wget https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

# 2. Build Salmon index
docker run --rm -v $(pwd):/data smill/omicsflow:1.0.0 \
  bash -c "export LC_ALL=C && salmon index \
  -t /data/transcriptome.fa.gz \
  -i /data/salmon_index \
  --threads 8"
```

### C. tx2gene.csv (required for gene-level DESeq2)

This file maps transcript IDs to gene IDs — required for Salmon → DESeq2 gene-level summarization.

**Generate automatically from your GTF:**
```bash
docker run --rm -v $(pwd):/data smill/omicsflow:1.0.0 \
  python3 /data/scripts/generate_tx2gene.py \
  --gtf /data/genes.gtf \
  --output /data/tx2gene.csv
```

Or use the one-liner:
```bash
docker run --rm -v $(pwd):/data smill/omicsflow:1.0.0 \
  bash -c "python3 -c \"
import re
with open('/data/genes.gtf') as f, open('/data/tx2gene.csv','w') as out:
    out.write('tx_id,gene_id\n')
    for line in f:
        if '\ttranscript\t' in line:
            tx = re.search('transcript_id \\\"([^\\\"]+)\\\"', line)
            gn = re.search('gene_id \\\"([^\\\"]+)\\\"', line)
            if tx and gn:
                out.write(f'{tx.group(1)},{gn.group(1)}\n')
\""
```

Then pass it to the pipeline:
```bash
nextflow run Millimono/OmicsFlow \
  --input samplesheet.csv \
  --star_index /data/star_index \
  --salmon_index /data/salmon_index \
  --tx2gene /data/tx2gene.csv \
  --outdir results/ \
  -profile docker
```

### D. What you do NOT need to install

Everything is already inside the Docker image `smill/omicsflow:1.0.0`:

| Tool | Without OmicsFlow | With OmicsFlow |
|---|---|---|
| FastQC | Manual install | ✅ Included |
| Trim Galore | Manual install | ✅ Included |
| STAR | Compile from source | ✅ Included |
| Salmon | Manual install | ✅ Included |
| Samtools | Compile from source | ✅ Included |
| DESeq2 | R + Bioconductor setup | ✅ Included |
| MultiQC | pip install | ✅ Included |
| BioPython | pip install | ✅ Included |
| numpy / pandas / matplotlib | pip install | ✅ Included |
| NanoStat / NanoPlot | pip install | ✅ Included |
| Kraken2 | Manual install | ✅ Included |
| Minimap2 | Manual install | ✅ Included |

---

## 🐳 Run with Docker only (no Nextflow required)

```bash
# Pull the image
docker pull smill/omicsflow:1.0.0

# Step 1 — Quality control
docker run --rm -v $(pwd)/data:/data smill/omicsflow:1.0.0 \
  bash -c "fastqc /data/*.fastq.gz --outdir /data/qc"

# Step 2 — Adapter trimming
docker run --rm -v $(pwd)/data:/data smill/omicsflow:1.0.0 \
  bash -c "trim_galore --paired --cores 4 \
  /data/sample_R1.fastq.gz /data/sample_R2.fastq.gz -o /data/trimmed"

# Step 3 — Alignment
docker run --rm -v $(pwd)/data:/data smill/omicsflow:1.0.0 \
  bash -c "STAR --runMode alignReads \
  --genomeDir /data/star_index \
  --readFilesIn /data/trimmed/sample_R1_val_1.fq.gz /data/trimmed/sample_R2_val_2.fq.gz \
  --readFilesCommand zcat \
  --outSAMtype BAM SortedByCoordinate \
  --outFileNamePrefix /data/aligned/sample. \
  --runThreadN 4"

# Step 4 — Quantification
docker run --rm -v $(pwd)/data:/data smill/omicsflow:1.0.0 \
  bash -c "export LC_ALL=C && salmon quant \
  --index /data/salmon_index \
  --libType A \
  -1 /data/trimmed/sample_R1_val_1.fq.gz \
  -2 /data/trimmed/sample_R2_val_2.fq.gz \
  --output /data/counts/sample \
  --threads 4 --validateMappings"

# Step 5 — BAM statistics
docker run --rm -v $(pwd)/data:/data smill/omicsflow:1.0.0 \
  bash -c "samtools flagstat /data/aligned/sample.Aligned.sortedByCoord.out.bam"

# Step 6 — Aggregated QC report
docker run --rm -v $(pwd)/data:/data smill/omicsflow:1.0.0 \
  bash -c "export LC_ALL=C.UTF-8 && multiqc /data --outdir /data/multiqc"

# Interactive R session (DESeq2, ggplot2...)
docker run --rm -it -v $(pwd)/data:/data smill/omicsflow:1.0.0 R

# Interactive Python session (biopython, pandas, matplotlib...)
docker run --rm -it -v $(pwd)/data:/data smill/omicsflow:1.0.0 python3
```

> **Windows users:** replace `$(pwd)` with `%cd%`

---

## 🚀 Run with Nextflow (recommended for production)

### Prerequisites
- [Nextflow](https://www.nextflow.io/docs/latest/install.html) ≥ 22.10
- [Docker](https://docs.docker.com/get-docker/) or Singularity
- Java 17+

### Full RNA-seq pipeline
```bash
# Clone the repository
git clone https://github.com/Millimono/OmicsFlow.git
cd OmicsFlow

# Run with all options
nextflow run rnaseq.nf \
  --input data/test/samplesheet.csv \
  --star_index /data/star_index \
  --salmon_index /data/salmon_index \
  --tx2gene /data/tx2gene.csv \
  --outdir results/ \
  --fdr_cutoff 0.05 \
  --lfc_threshold 1.0 \
  -profile docker
```

### Skip DESeq2 (alignment + quantification only)
```bash
nextflow run rnaseq.nf \
  --input samplesheet.csv \
  --star_index /data/star_index \
  --salmon_index /data/salmon_index \
  --skip_deseq2 \
  --outdir results/ \
  -profile docker
```

### Run on HPC cluster (SLURM)
```bash
nextflow run rnaseq.nf \
  --input samplesheet.csv \
  --star_index /data/star_index \
  --salmon_index /data/salmon_index \
  --tx2gene /data/tx2gene.csv \
  --outdir results/ \
  -profile cluster
```

### All parameters

| Parameter | Default | Description |
|---|---|---|
| `--input` | required | Path to samplesheet CSV |
| `--star_index` | null | Path to STAR genome index |
| `--salmon_index` | null | Path to Salmon transcriptome index |
| `--tx2gene` | null | Path to tx2gene.csv (transcript→gene) |
| `--outdir` | results | Output directory |
| `--fdr_cutoff` | 0.05 | DESeq2 FDR threshold |
| `--lfc_threshold` | 1.0 | DESeq2 log2 fold-change threshold |
| `--min_counts` | 10 | Minimum counts filter for DESeq2 |
| `--skip_trimming` | false | Skip Trim Galore step |
| `--skip_fastqc` | false | Skip FastQC step |
| `--skip_deseq2` | false | Skip DESeq2 step |
| `--skip_multiqc` | false | Skip MultiQC step |

---

## 🗂️ Project Structure

```
OmicsFlow/
├── rnaseq.nf                    # ✅ RNA-seq pipeline (stable)
├── workflows/
│   ├── longread.nf              # 🚧 Nanopore pipeline (in development)
│   └── metagenomics.nf          # 🚧 Metagenomic pipeline (in development)
├── modules/
│   ├── qc/                      # FastQC, MultiQC, Trim Galore
│   ├── alignment/               # STAR, Samtools, Minimap2
│   └── quantification/          # Salmon, DESeq2
├── analysis/
│   ├── deseq2.R                 # Standalone DESeq2 script
│   ├── plots.py                 # Visualization scripts
│   └── report.Rmd               # HTML report template
├── scripts/
│   ├── prepare_references.sh    # Build STAR + Salmon index + tx2gene in one command
│   └── generate_tx2gene.py      # Generate tx2gene.csv from GTF
├── containers/
│   └── Dockerfile               # Reproducible environment
├── data/test/
│   └── samplesheet.csv          # Example samplesheet
├── docs/
│   ├── usage.md                 # Full usage guide
│   └── images/                  # Pipeline figures
├── .github/workflows/ci.yml     # GitHub Actions CI/CD
└── nextflow.config              # Execution profiles
```

---

## 📈 Results & Outputs

```
results/
├── qc/
│   ├── fastqc/                  # Per-sample FastQC HTML reports
│   └── multiqc_report.html      # Aggregated QC report
├── trimmed/
│   └── logs/                    # Trim Galore reports
├── aligned/
│   ├── sample.Aligned.sortedByCoord.out.bam
│   └── sample.Log.final.out     # Mapping statistics
├── counts/
│   └── salmon/quant.sf          # Transcript-level counts
├── deseq2/
│   ├── deseq2_results.csv       # DE genes table
│   ├── normalized_counts.csv    # VST-normalized counts
│   ├── volcano_plot.pdf         # Volcano plot
│   ├── heatmap_top50.pdf        # Top 50 DE genes
│   └── pca_plot.pdf             # Sample clustering PCA
└── pipeline_info/
    ├── execution_report.html
    └── execution_timeline.html
```

---

## 📊 Workflows in Detail

### 1. RNA-seq Pipeline (`rnaseq.nf`) ✅ Stable

```
Input FASTQ
    │
    ▼
[FastQC] ──────────────────────> QC report
    │
    ▼
[Trim Galore] ──> Trimmed reads
    │
    ├──────────────────────────────────────────┐
    ▼                                          ▼
[STAR] ──> BAM ──> [Samtools]          [Salmon] ──> quant.sf
           (QC + IGV visualization)         │
                                            ▼
                                    [tximport] ──> [DESeq2]
                                                       │
                                                       ▼
                                            volcano · heatmap · PCA
    │
    ▼
[MultiQC] ──> Aggregated HTML report
```

### 2. Long-read Pipeline (`longread.nf`) 🚧 In development

Tools already available in Docker — Nextflow modules coming soon.

```bash
# Use individually now:
docker run --rm -v $(pwd):/data smill/omicsflow:1.0.0 \
  bash -c "NanoStat --fastq /data/sample.fastq.gz"

docker run --rm -v $(pwd):/data smill/omicsflow:1.0.0 \
  bash -c "minimap2 -ax splice /data/reference.fa /data/sample.fastq.gz \
  | samtools sort -o /data/aligned.bam"
```

### 3. Metagenomic Pipeline (`metagenomics.nf`) 🚧 In development

Kraken2 already available in Docker — Nextflow modules coming soon.

```bash
# Use individually now:
docker run --rm -v $(pwd):/data smill/omicsflow:1.0.0 \
  bash -c "kraken2 --db /data/kraken2_db --paired \
  /data/R1.fastq.gz /data/R2.fastq.gz \
  --output /data/kraken2_output.txt \
  --report /data/kraken2_report.txt"
```

---

## 🛠️ Technical Stack

| Category | Tools | Versions |
|---|---|---|
| **Pipeline orchestration** | Nextflow DSL2 | ≥ 22.10 |
| **Containerization** | Docker · Singularity | 28.x |
| **QC** | FastQC · MultiQC · NanoStat · NanoPlot | 0.12.1 · 1.35 · 1.6.0 |
| **Alignment** | STAR · Minimap2 | 2.7.11b · 2.31 |
| **Quantification** | Salmon | 1.12.0 |
| **Variant calling** | Samtools · BCFtools | 1.23.1 |
| **Metagenomics** | Kraken2 | 2.1.3 |
| **Statistical analysis** | DESeq2 · tximport · R | R 4.5.2 |
| **Visualization** | ggplot2 · matplotlib · seaborn | — |
| **Languages** | Python · R · Bash · C · C++ | Python 3.x |
| **CI/CD** | GitHub Actions | — |
| **Documentation** | GitHub Pages | — |

---

## 🧪 Test Data

| Dataset | Source | Size | Used for |
|---|---|---|---|
| GSE110004 / SRR6357070-71 (4 samples) | nf-core test datasets | ~8 MB | RNA-seq validation |
| S. cerevisiae R64-1-1 genome | nf-core test datasets | ~230 KB | Reference genome |
| S. cerevisiae gene annotation | nf-core test datasets | ~200 KB | Gene annotation |

---

## ⚙️ Configuration Profiles

```groovy
profiles {
    docker  { docker.enabled = true; process.executor = 'local' }
    cluster { process.executor = 'slurm'; singularity.enabled = true }
    cloud   { process.executor = 'awsbatch'; aws.region = 'ca-central-1' }
    test    { params.input = "${projectDir}/data/test/samplesheet.csv" }
    stub    { /* CI testing without real tools */ }
}
```

---

## 🔗 Link to Research

- **MalariaScan** — AI detection of malaria. Prix Jean-Marc Léger, UdeM 2025.
- **HAtt-CNN** — Adaptive visual attention for CNN interpretability. *(Under review 2026)*
- **EpitopeNet** — Prototype learning for mammography classification. *(Under review 2026)*  
  → [github.com/Millimono/EpitopeNet](https://github.com/Millimono/EpitopeNet)

---

## 📚 Documentation

Full documentation: **[millimono.github.io/OmicsFlow](https://millimono.github.io/OmicsFlow)**

---

## 📄 Citation

```
Millimono, S. (2026). OmicsFlow: A modular containerized NGS pipeline
for reproducible multi-omics analysis (v1.0.2). Zenodo.
https://doi.org/10.5281/zenodo.20677900
```

---

## 👤 Author

**Sory Millimono** — PhD Candidate in AI · Bioinformatician  
Mohammed V University – ENSIAS · Université de Montréal

- 📧 millimono64.sm@gmail.com
- 🔗 [LinkedIn](https://linkedin.com/in/sory-millimono-ai-searcher-820314162)
- 🎓 [Google Scholar](https://scholar.google.com/citations?user=5M-zcxYAAAAJ) — h-index 1 · 24 citations
- 🔬 [ORCID: 0009-0005-1960-9136](https://orcid.org/0009-0005-1960-9136)

---

## 📜 License

MIT License — see [LICENSE](LICENSE) for details.
