#!/bin/bash

#SBATCH --job-name=gblocks
#SBATCH --partition=small
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

# Trim codon alignments with Gblocks using the
# parameters applied in the final analysis.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

IN="$ROOT/output/codon_alignments"
OUT="$ROOT/output/gblocks_alignments"

GBLOCKS="${GBLOCKS:-$(command -v Gblocks || true)}"

if [ -z "$GBLOCKS" ]; then
    echo "ERROR: Gblocks not found" >&2
    exit 1
fi

mkdir -p "$OUT"

for fasta in "$IN"/*.codon.full.fasta
do
    gene=$(basename "$fasta" .codon.full.fasta)
    work="$OUT/$gene"

    mkdir -p "$work"
    cp "$fasta" "$work/${gene}.codon.full.fasta"

    nseq=$(grep -c '^>' "$fasta")
    half=$(((nseq + 1) / 2))
    b1=$((half + 2))
    b2=$((half + 4))

    (
        cd "$work"

        "$GBLOCKS" \
            "${gene}.codon.full.fasta" \
            -t=c \
            -b1="$b1" \
            -b2="$b2" \
            -b3=8 \
            -b4=15 \
            -b5=a
    )

    if [ ! -s "$work/${gene}.codon.full.fasta-gb" ]; then
        echo "ERROR: Gblocks output missing for $gene" >&2
        exit 1
    fi
done