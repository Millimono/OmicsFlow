#!/bin/bash
# ============================================================================
#  OmicsFlow — prepare_references.sh
#  Prepare all reference files needed for the RNA-seq pipeline:
#    1. STAR genome index
#    2. Salmon transcriptome index
#    3. tx2gene.csv (transcript → gene mapping)
#
#  Usage:
#    bash scripts/prepare_references.sh \
#      --genome /path/to/genome.fa \
#      --gtf    /path/to/genes.gtf \
#      --outdir /path/to/references/
#
#  Via Docker (recommended):
#    docker run --rm -v $(pwd):/data smill/omicsflow:1.0.0 \
#      bash /data/scripts/prepare_references.sh \
#        --genome /data/genome.fa \
#        --gtf    /data/genes.gtf \
#        --outdir /data/references/
#
#  Requirements: smill/omicsflow:1.0.0 Docker image
#  Author: Sory Millimono — github.com/Millimono/OmicsFlow
# ============================================================================

set -euo pipefail

# ── COLORS ───────────────────────────────────────────────────────────────────
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# ── DEFAULTS ─────────────────────────────────────────────────────────────────
GENOME=""
GTF=""
OUTDIR="references"
THREADS=8
GENOME_SA_BASES=14   # 14 for human/mouse, 11 for chr22, 7 for yeast
SKIP_STAR=false
SKIP_SALMON=false
SKIP_TX2GENE=false

# ── USAGE ────────────────────────────────────────────────────────────────────
usage() {
    echo ""
    echo "Usage: bash prepare_references.sh [OPTIONS]"
    echo ""
    echo "Required:"
    echo "  --genome FILE     Reference genome FASTA (.fa or .fa.gz)"
    echo "  --gtf    FILE     Gene annotation GTF (.gtf or .gtf.gz)"
    echo ""
    echo "Optional:"
    echo "  --outdir DIR      Output directory (default: references/)"
    echo "  --threads INT     Number of threads (default: 8)"
    echo "  --sa-bases INT    STAR genomeSAindexNbases (default: 14)"
    echo "                    14=human/mouse, 11=chr22 only, 7=yeast"
    echo "  --skip-star       Skip STAR index generation"
    echo "  --skip-salmon     Skip Salmon index generation"
    echo "  --skip-tx2gene    Skip tx2gene.csv generation"
    echo "  -h, --help        Show this help"
    echo ""
    echo "Examples:"
    echo "  # Human GRCh38"
    echo "  bash prepare_references.sh \\"
    echo "    --genome Homo_sapiens.GRCh38.dna.fa \\"
    echo "    --gtf    Homo_sapiens.GRCh38.109.gtf \\"
    echo "    --outdir references/ --threads 8 --sa-bases 14"
    echo ""
    echo "  # Mouse GRCm39"
    echo "  bash prepare_references.sh \\"
    echo "    --genome Mus_musculus.GRCm39.dna.fa \\"
    echo "    --gtf    Mus_musculus.GRCm39.109.gtf \\"
    echo "    --outdir references/ --threads 8 --sa-bases 14"
    echo ""
    echo "  # S. cerevisiae (small genome)"
    echo "  bash prepare_references.sh \\"
    echo "    --genome Saccharomyces_cerevisiae.R64-1-1.dna.fa \\"
    echo "    --gtf    Saccharomyces_cerevisiae.R64-1-1.gtf \\"
    echo "    --outdir references/ --sa-bases 7"
    echo ""
}

# ── PARSE ARGUMENTS ───────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case $1 in
        --genome)     GENOME="$2";       shift 2 ;;
        --gtf)        GTF="$2";          shift 2 ;;
        --outdir)     OUTDIR="$2";       shift 2 ;;
        --threads)    THREADS="$2";      shift 2 ;;
        --sa-bases)   GENOME_SA_BASES="$2"; shift 2 ;;
        --skip-star)  SKIP_STAR=true;    shift ;;
        --skip-salmon) SKIP_SALMON=true; shift ;;
        --skip-tx2gene) SKIP_TX2GENE=true; shift ;;
        -h|--help)    usage; exit 0 ;;
        *) echo "Unknown option: $1"; usage; exit 1 ;;
    esac
done

# ── VALIDATE ──────────────────────────────────────────────────────────────────
if [[ -z "$GENOME" || -z "$GTF" ]]; then
    echo -e "${RED}ERROR: --genome and --gtf are required.${NC}"
    usage
    exit 1
fi

if [[ ! -f "$GENOME" ]]; then
    echo -e "${RED}ERROR: Genome file not found: $GENOME${NC}"
    exit 1
fi

if [[ ! -f "$GTF" ]]; then
    echo -e "${RED}ERROR: GTF file not found: $GTF${NC}"
    exit 1
fi

# ── SETUP ─────────────────────────────────────────────────────────────────────
mkdir -p "$OUTDIR"
STAR_INDEX="$OUTDIR/star_index"
SALMON_INDEX="$OUTDIR/salmon_index"
TRANSCRIPTOME="$OUTDIR/transcriptome.fa"
TX2GENE="$OUTDIR/tx2gene.csv"
LOG="$OUTDIR/prepare_references.log"

echo "" | tee "$LOG"
echo -e "${BLUE}╔══════════════════════════════════════════════════╗${NC}" | tee -a "$LOG"
echo -e "${BLUE}║   OmicsFlow — Reference Preparation              ║${NC}" | tee -a "$LOG"
echo -e "${BLUE}╠══════════════════════════════════════════════════╣${NC}" | tee -a "$LOG"
echo -e "${BLUE}║  Genome  : $GENOME${NC}" | tee -a "$LOG"
echo -e "${BLUE}║  GTF     : $GTF${NC}" | tee -a "$LOG"
echo -e "${BLUE}║  Output  : $OUTDIR${NC}" | tee -a "$LOG"
echo -e "${BLUE}║  Threads : $THREADS${NC}" | tee -a "$LOG"
echo -e "${BLUE}╚══════════════════════════════════════════════════╝${NC}" | tee -a "$LOG"
echo "" | tee -a "$LOG"

# Decompress if needed
GENOME_UNZIPPED="$GENOME"
if [[ "$GENOME" == *.gz ]]; then
    echo -e "${YELLOW}Decompressing genome...${NC}" | tee -a "$LOG"
    GENOME_UNZIPPED="${GENOME%.gz}"
    gunzip -c "$GENOME" > "$GENOME_UNZIPPED"
fi

GTF_UNZIPPED="$GTF"
if [[ "$GTF" == *.gz ]]; then
    echo -e "${YELLOW}Decompressing GTF...${NC}" | tee -a "$LOG"
    GTF_UNZIPPED="${GTF%.gz}"
    gunzip -c "$GTF" > "$GTF_UNZIPPED"
fi

# ── STEP 1: STAR INDEX ───────────────────────────────────────────────────────
if [[ "$SKIP_STAR" == false ]]; then
    echo -e "${YELLOW}[1/3] Building STAR genome index...${NC}" | tee -a "$LOG"
    echo -e "      This may take 10-45 minutes depending on genome size." | tee -a "$LOG"
    mkdir -p "$STAR_INDEX"

    STAR \
        --runMode genomeGenerate \
        --genomeDir "$STAR_INDEX" \
        --genomeFastaFiles "$GENOME_UNZIPPED" \
        --sjdbGTFfile "$GTF_UNZIPPED" \
        --genomeSAindexNbases "$GENOME_SA_BASES" \
        --runThreadN "$THREADS" \
        2>&1 | tee -a "$LOG"

    echo -e "${GREEN}✅ STAR index built: $STAR_INDEX${NC}" | tee -a "$LOG"
else
    echo -e "${YELLOW}⏭️  Skipping STAR index (--skip-star)${NC}" | tee -a "$LOG"
fi

# ── STEP 2: TRANSCRIPTOME FASTA + SALMON INDEX ───────────────────────────────
if [[ "$SKIP_SALMON" == false ]]; then
    echo "" | tee -a "$LOG"
    echo -e "${YELLOW}[2/3] Extracting transcriptome and building Salmon index...${NC}" | tee -a "$LOG"

    # Extract transcript sequences from genome + GTF using Python
    python3 - << PYEOF 2>&1 | tee -a "$LOG"
import re
import sys

genome_file = "$GENOME_UNZIPPED"
gtf_file    = "$GTF_UNZIPPED"
out_file    = "$TRANSCRIPTOME"

print(f"Loading genome from {genome_file}...")
genome = {}
current_chrom = None
current_seq   = []

with open(genome_file) as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if current_chrom:
                genome[current_chrom] = "".join(current_seq)
            current_chrom = line[1:].split()[0]
            current_seq   = []
        else:
            current_seq.append(line)
    if current_chrom:
        genome[current_chrom] = "".join(current_seq)

print(f"Loaded {len(genome)} chromosomes/contigs")

print(f"Extracting transcripts from {gtf_file}...")
transcripts = {}

with open(gtf_file) as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.strip().split("\t")
        if len(fields) < 9 or fields[2] != "exon":
            continue

        chrom  = fields[0]
        start  = int(fields[3]) - 1  # GTF is 1-based
        end    = int(fields[4])
        strand = fields[6]
        attrs  = fields[8]

        tx_match   = re.search(r'transcript_id "([^"]+)"', attrs)
        gene_match = re.search(r'gene_id "([^"]+)"', attrs)

        if not tx_match or not gene_match:
            continue

        tx_id   = tx_match.group(1)
        gene_id = gene_match.group(1)

        if chrom not in genome:
            continue

        if tx_id not in transcripts:
            transcripts[tx_id] = {"gene": gene_id, "strand": strand,
                                   "chrom": chrom, "exons": []}
        transcripts[tx_id]["exons"].append((start, end))

print(f"Found {len(transcripts)} transcripts")

# Write FASTA
complement = str.maketrans("ACGTacgt", "TGCAtgca")
n_written = 0

with open(out_file, "w") as out:
    for tx_id, info in transcripts.items():
        chrom  = info["chrom"]
        strand = info["strand"]
        exons  = sorted(info["exons"])
        seq    = "".join(genome[chrom][s:e] for s, e in exons)

        if strand == "-":
            seq = seq.translate(complement)[::-1]

        if seq:
            out.write(f">{tx_id}\n{seq}\n")
            n_written += 1

print(f"✅ Transcriptome written: {n_written} sequences → {out_file}")
PYEOF

    # Build Salmon index
    echo -e "${YELLOW}Building Salmon index...${NC}" | tee -a "$LOG"
    export LC_ALL=C
    salmon index \
        -t "$TRANSCRIPTOME" \
        -i "$SALMON_INDEX" \
        --threads "$THREADS" \
        2>&1 | tee -a "$LOG"

    echo -e "${GREEN}✅ Salmon index built: $SALMON_INDEX${NC}" | tee -a "$LOG"
else
    echo -e "${YELLOW}⏭️  Skipping Salmon index (--skip-salmon)${NC}" | tee -a "$LOG"
fi

# ── STEP 3: TX2GENE ──────────────────────────────────────────────────────────
if [[ "$SKIP_TX2GENE" == false ]]; then
    echo "" | tee -a "$LOG"
    echo -e "${YELLOW}[3/3] Generating tx2gene.csv...${NC}" | tee -a "$LOG"

    python3 - << PYEOF 2>&1 | tee -a "$LOG"
import re

gtf_file = "$GTF_UNZIPPED"
out_file = "$TX2GENE"
n = 0

with open(gtf_file) as f, open(out_file, "w") as out:
    out.write("tx_id,gene_id\n")
    for line in f:
        if line.startswith("#") or "\ttranscript\t" not in line:
            continue
        tx   = re.search(r'transcript_id "([^"]+)"', line)
        gene = re.search(r'gene_id "([^"]+)"', line)
        if tx and gene:
            out.write(f"{tx.group(1)},{gene.group(1)}\n")
            n += 1

print(f"✅ tx2gene.csv written: {n} transcripts → {out_file}")
PYEOF

    echo -e "${GREEN}✅ tx2gene.csv: $TX2GENE${NC}" | tee -a "$LOG"
else
    echo -e "${YELLOW}⏭️  Skipping tx2gene (--skip-tx2gene)${NC}" | tee -a "$LOG"
fi

# ── SUMMARY ───────────────────────────────────────────────────────────────────
echo "" | tee -a "$LOG"
echo -e "${GREEN}╔══════════════════════════════════════════════════╗${NC}" | tee -a "$LOG"
echo -e "${GREEN}║   ✅ Reference preparation complete!             ║${NC}" | tee -a "$LOG"
echo -e "${GREEN}╠══════════════════════════════════════════════════╣${NC}" | tee -a "$LOG"
echo -e "${GREEN}║  STAR index  : $STAR_INDEX${NC}" | tee -a "$LOG"
echo -e "${GREEN}║  Salmon index: $SALMON_INDEX${NC}" | tee -a "$LOG"
echo -e "${GREEN}║  tx2gene     : $TX2GENE${NC}" | tee -a "$LOG"
echo -e "${GREEN}╠══════════════════════════════════════════════════╣${NC}" | tee -a "$LOG"
echo -e "${GREEN}║  Next step — run the pipeline:                   ║${NC}" | tee -a "$LOG"
echo -e "${GREEN}║                                                  ║${NC}" | tee -a "$LOG"
echo -e "${GREEN}║  nextflow run Millimono/OmicsFlow \\              ║${NC}" | tee -a "$LOG"
echo -e "${GREEN}║    --input samplesheet.csv \\                     ║${NC}" | tee -a "$LOG"
echo -e "${GREEN}║    --star_index $STAR_INDEX \\  ║${NC}" | tee -a "$LOG"
echo -e "${GREEN}║    --salmon_index $SALMON_INDEX \\║${NC}" | tee -a "$LOG"
echo -e "${GREEN}║    --tx2gene $TX2GENE \\          ║${NC}" | tee -a "$LOG"
echo -e "${GREEN}║    --outdir results/ -profile docker             ║${NC}" | tee -a "$LOG"
echo -e "${GREEN}╚══════════════════════════════════════════════════╝${NC}" | tee -a "$LOG"
echo "" | tee -a "$LOG"
echo -e "Full log saved to: $LOG" | tee -a "$LOG"
