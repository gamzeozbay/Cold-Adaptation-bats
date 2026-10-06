#!/bin/bash

#SBATCH --job-name=mafft_prot
#SBATCH --partition=small
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

# Align final selected protein sequences with MAFFT.
#
# Input:
#   output/proteins_by_gene/*.faa
#
# Output:
#   output/protein_alignments/*.protein.aln.faa

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

IN="$ROOT/output/proteins_by_gene"
OUT="$ROOT/output/protein_alignments"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

command -v mafft >/dev/null || {
    echo "ERROR: mafft not found" >&2
    exit 1
}

mkdir -p "$OUT"

for faa in "$IN"/*.faa
do
    gene=$(basename "$faa" .faa)
    aln="$OUT/${gene}.protein.aln.faa"

    mafft \
        --auto \
        --thread "$THREADS" \
        "$faa" \
        > "$aln"
done