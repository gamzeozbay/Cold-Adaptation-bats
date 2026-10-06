#!/bin/bash

#SBATCH --job-name=pal2nal
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G

# Convert protein alignments and matching CDS sequences
# into codon alignments with PAL2NAL.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

PROT="$ROOT/output/protein_alignments"
CDS="$ROOT/output/cds_by_gene"
OUT="$ROOT/output/codon_alignments"

PAL2NAL="${PAL2NAL:-$(command -v pal2nal.pl || true)}"

if [ -z "$PAL2NAL" ]; then
    echo "ERROR: pal2nal.pl not found" >&2
    exit 1
fi

mkdir -p "$OUT"

for aln in "$PROT"/*.protein.aln.faa
do
    gene=$(basename "$aln" .protein.aln.faa)
    cds="$CDS/${gene}.cds.fasta"

    if [ ! -s "$cds" ]; then
        echo "ERROR: missing CDS file for $gene" >&2
        exit 1
    fi

    if [[ "$gene" == "ATP6" || "$gene" == "ND3" || "$gene" == "ND4" ]]; then
        codontable=2
    else
        codontable=1
    fi

    perl "$PAL2NAL" \
        "$aln" \
        "$cds" \
        -output fasta \
        -codontable "$codontable" \
        > "$OUT/${gene}.codon.full.fasta"
done