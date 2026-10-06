#!/bin/bash

# Mask sequence-specific erroneous regions in Gblocks-trimmed
# nucleotide alignments using TAPER.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

IN="$ROOT/output/gblocks_alignments"
OUT="$ROOT/output/taper_alignments"

TAPER_DIR="${TAPER_DIR:?Please set TAPER_DIR to the TAPER installation directory}"
JULIA="${JULIA:-$(command -v julia || true)}"

if [ -z "$JULIA" ]; then
    echo "ERROR: Julia not found" >&2
    exit 1
fi

if [ ! -s "$TAPER_DIR/correction_multi.jl" ]; then
    echo "ERROR: correction_multi.jl not found" >&2
    exit 1
fi

mkdir -p "$OUT"

for fasta in "$IN"/*/*.codon.full.fasta-gb
do
    gene=$(basename "$fasta" .codon.full.fasta-gb)

    "$JULIA" \
        "$TAPER_DIR/correction_multi.jl" \
        -m N \
        -a N \
        "$fasta" \
        > "$OUT/${gene}.codon.gblocks.taper.fasta"
done