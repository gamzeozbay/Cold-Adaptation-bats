#!/bin/bash

#SBATCH --job-name=iqtree
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

# Infer nuclear gene trees from TAPER-masked codon alignments.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

IN="$ROOT/output/taper_alignments"
OUT="$ROOT/output/gene_trees"

THREADS="${SLURM_CPUS_PER_TASK:-4}"

IQTREE="${IQTREE:-$(command -v iqtree2 || true)}"

if [ -z "$IQTREE" ]; then
    echo "ERROR: iqtree2 not found" >&2
    exit 1
fi

mkdir -p "$OUT"

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

for aln in "$IN"/*.codon.gblocks.taper.fasta
do
    gene=$(basename "$aln" .codon.gblocks.taper.fasta)

    if [[ "$gene" == "ATP6" || "$gene" == "ND3" || "$gene" == "ND4" ]]; then
        continue
    fi

    species_aln="$WORK/${gene}.species.fasta"

    awk '
        /^>/ {
            split(substr($0, 2), a, "|")
            print ">" a[2]
            next
        }
        { print }
    ' "$aln" > "$species_aln"

    "$IQTREE" \
        -s "$species_aln" \
        -pre "$OUT/$gene" \
        -m MFP \
        -B 1000 \
        -alrt 1000 \
        -nt "$THREADS"
done