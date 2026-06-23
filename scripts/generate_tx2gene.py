#!/usr/bin/env python3
"""
OmicsFlow — generate_tx2gene.py
Generate tx2gene.csv from a GTF annotation file.

Usage:
    python3 generate_tx2gene.py --gtf genes.gtf --output tx2gene.csv

Via Docker:
    docker run --rm -v $(pwd):/data smill/omicsflow:1.0.0 \
      python3 /data/scripts/generate_tx2gene.py \
      --gtf /data/genes.gtf --output /data/tx2gene.csv
"""

import re
import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate tx2gene.csv from GTF for tximport/DESeq2"
    )
    parser.add_argument("--gtf",    required=True, help="Path to GTF annotation file")
    parser.add_argument("--output", required=True, help="Output tx2gene.csv path")
    parser.add_argument("--feature", default="transcript",
                        help="GTF feature type to parse (default: transcript)")
    return parser.parse_args()


def extract_attribute(line, key):
    """Extract attribute value from GTF line."""
    match = re.search(rf'{key} "([^"]+)"', line)
    return match.group(1) if match else None


def generate_tx2gene(gtf_path, output_path, feature="transcript"):
    n_written = 0
    n_skipped = 0

    print(f"Reading GTF: {gtf_path}")
    print(f"Feature type: {feature}")

    with open(gtf_path, "r") as gtf, open(output_path, "w") as out:
        out.write("tx_id,gene_id\n")

        for line in gtf:
            # Skip comment lines
            if line.startswith("#"):
                continue

            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            # Only process requested feature type
            if fields[2] != feature:
                continue

            tx_id   = extract_attribute(fields[8], "transcript_id")
            gene_id = extract_attribute(fields[8], "gene_id")

            if tx_id and gene_id:
                out.write(f"{tx_id},{gene_id}\n")
                n_written += 1
            else:
                n_skipped += 1

    print(f"✅ tx2gene.csv generated: {n_written} transcripts written")
    if n_skipped > 0:
        print(f"⚠️  {n_skipped} lines skipped (missing transcript_id or gene_id)")
    print(f"Output: {output_path}")


if __name__ == "__main__":
    args = parse_args()
    try:
        generate_tx2gene(args.gtf, args.output, args.feature)
    except FileNotFoundError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
