# ============================================================================
#  OmicsFlow — Dockerfile v3 (final)
#  Author  : Sory Millimono
#  Version : 1.0.0
#
#  Build : docker build -t omicsflow:1.0.0 containers/
#  Run   : docker run --rm omicsflow:1.0.0 fastqc --version
# ============================================================================

FROM continuumio/miniconda3:24.1.2-0

ENV DEBIAN_FRONTEND=noninteractive
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

RUN mkdir -p /data /results
WORKDIR /data

# ============================================================================
# 1. JAVA (pour FastQC)
# ============================================================================
RUN apt-get update && apt-get install -y --no-install-recommends \
    openjdk-17-jre-headless \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# ============================================================================
# 2. INSTALLER MAMBA (solver rapide — remplace conda classic)
# ============================================================================
RUN conda install -n base -c conda-forge mamba -y \
    && conda clean -afy

# ============================================================================
# 3. OUTILS NGS PRINCIPAUX via mamba
# ============================================================================
RUN mamba install -y -c bioconda -c conda-forge \
    fastqc \
    trim-galore \
    samtools \
    minimap2 \
    bcftools \
    && mamba clean -afy

# ============================================================================
# 4. STAR — Aligneur RNA-seq splice-aware
# ============================================================================
RUN mamba install -y -c bioconda -c conda-forge star \
    && mamba clean -afy

# ============================================================================
# 5. SALMON — Quantification de l'expression génique
# ============================================================================
RUN mamba install -y -c bioconda -c conda-forge salmon \
    && mamba clean -afy

# ============================================================================
# 6. PYTHON PACKAGES via pip
#    multiqc, NanoStat, NanoPlot, biopython, numpy, pandas, scipy,
#    matplotlib, seaborn, pysam
# ============================================================================
RUN pip install --no-cache-dir \
    multiqc \
    NanoStat \
    NanoPlot \
    biopython \
    numpy \
    pandas \
    scipy \
    matplotlib \
    seaborn \
    pysam

# ============================================================================
# 7. R + BIOCONDUCTOR — DESeq2, tximport, ggplot2, pheatmap
# ============================================================================
RUN mamba install -y -c bioconda -c conda-forge \
    r-base \
    r-ggplot2 \
    r-pheatmap \
    r-rcolorbrewer \
    bioconductor-deseq2 \
    bioconductor-tximport \
    && mamba clean -afy

# ============================================================================
# 8. VÉRIFICATION FINALE — tous les outils
# ============================================================================
RUN fastqc --version \
    && STAR --version \
    && samtools --version | head -1 \
    && salmon --version \
    && minimap2 --version \
    && multiqc --version \
    && NanoStat --version \
    && python3 -c "import numpy, pandas, Bio; print('Python OK')" \
    && R -e "library(DESeq2); cat('DESeq2 OK\n')" \
    && echo "✅ OmicsFlow prêt !"

# ============================================================================
# METADATA
# ============================================================================
LABEL maintainer="Sory Millimono <millimono64.sm@gmail.com>"
LABEL version="1.0.0"
LABEL description="OmicsFlow NGS pipeline: RNA-seq, Long-read, Metagenomics"
LABEL tools="fastqc,trimgalore,star,samtools,salmon,minimap2,multiqc,nanostat,kraken2,deseq2,biopython"

CMD ["bash"]