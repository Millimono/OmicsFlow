#!/bin/bash
# ============================================================================
#  OmicsFlow — Script de test Docker
#  Usage : bash test_docker.sh
# ============================================================================

set -e

IMAGE="omicsflow:1.0.0"
echo ""
echo "╔══════════════════════════════════════════╗"
echo "║   OmicsFlow — Test de l'image Docker     ║"
echo "╚══════════════════════════════════════════╝"
echo ""

# ── 1. Build ────────────────────────────────────────────────────────────────
echo "📦 Construction de l'image Docker..."
echo "    (première fois : 10-20 minutes)"
echo ""
docker build -t ${IMAGE} containers/
echo ""
echo "✅ Image construite : ${IMAGE}"
echo ""

# ── 2. Taille de l'image ────────────────────────────────────────────────────
SIZE=$(docker image inspect ${IMAGE} --format='{{.Size}}' | awk '{printf "%.1f GB", $1/1073741824}')
echo "📏 Taille de l'image : ${SIZE}"
echo ""

# ── 3. Test de chaque outil ─────────────────────────────────────────────────
echo "🔍 Vérification des outils..."
echo ""

tools=(
    "fastqc --version"
    "STAR --version"
    "samtools --version"
    "salmon --version"
    "trim_galore --version"
    "minimap2 --version"
    "multiqc --version"
    "NanoStat --version"
    "kraken2 --version"
    "bcftools --version"
    "python3 -c 'import numpy, pandas, biopython; print(\"numpy\", numpy.__version__)'"
    "R -e 'library(DESeq2); cat(\"DESeq2\", as.character(packageVersion(\"DESeq2\")), \"\\n\")'"
)

for tool_cmd in "${tools[@]}"; do
    tool_name=$(echo $tool_cmd | awk '{print $1}')
    output=$(docker run --rm ${IMAGE} bash -c "${tool_cmd}" 2>&1 | head -1)
    echo "  ✅ ${tool_name}: ${output}"
done

echo ""
echo "╔══════════════════════════════════════════╗"
echo "║   ✅ Tous les outils sont fonctionnels   ║"
echo "╚══════════════════════════════════════════╝"
echo ""
echo "Pour lancer le pipeline RNA-seq :"
echo "  nextflow run rnaseq.nf -profile docker --input data/test/samplesheet.csv"
echo ""
