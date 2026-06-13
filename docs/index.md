# 🧬 OmicsFlow Documentation

> A modular, containerized NGS pipeline for RNA-seq, long-read, and metagenomic analysis

**GitHub:** [github.com/Millimono/OmicsFlow](https://github.com/Millimono/OmicsFlow)  
**Docker Hub:** [hub.docker.com/r/smill/omicsflow](https://hub.docker.com/r/smill/omicsflow)  
**Author:** Sory Millimono — PhD Candidate in AI & Bioinformatics

---

## Quick Start

```bash
# Pull the Docker image
docker pull smill/omicsflow:1.0.0

# Run FastQC on your data
docker run --rm -v $(pwd)/data:/data smill/omicsflow:1.0.0 \
  bash -c "fastqc /data/*.fastq.gz --outdir /data/qc"
```

---

## Contents

- [Installation](installation.md)
- [Usage & parameters](usage.md)
- [Output description](outputs.md)
- [Adding new modules](contributing.md)

---

## Available Tools

All tools are pre-installed in the Docker image `smill/omicsflow:1.0.0`:

| Tool | Version | Use |
|---|---|---|
| FastQC | 0.12.1 | Quality control |
| Trim Galore | latest | Adapter trimming |
| STAR | 2.7.11b | RNA-seq alignment |
| Salmon | 1.12.0 | Quantification |
| Samtools | 1.23.1 | BAM manipulation |
| MultiQC | 1.35 | Aggregated QC report |
| DESeq2 | latest | Differential expression |
| NanoStat | 1.6.0 | Nanopore QC |
| Minimap2 | 2.31 | Long-read alignment |
| Kraken2 | 2.1.3 | Metagenomics |
| BioPython | latest | Python bioinformatics |
| R 4.5.2 | 4.5.2 | Statistical analysis |

---

## Contact

📧 millimono64.sm@gmail.com  
🔗 [LinkedIn](https://linkedin.com/in/sory-millimono-ai-searcher-820314162)  
🎓 [Google Scholar](https://scholar.google.com/citations?user=5M-zcxYAAAAJ)  
🔬 [ORCID: 0009-0005-1960-9136](https://orcid.org/0009-0005-1960-9136)
