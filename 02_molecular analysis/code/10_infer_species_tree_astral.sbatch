#!/bin/bash

#SBATCH --job-name=astral
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

# Infer the nuclear species tree from IQ-TREE gene trees using ASTRAL.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

TREE_DIR="$ROOT/output/gene_trees"
OUT="$ROOT/output/species_tree"

ASTRAL_JAR="${ASTRAL_JAR:?Please set ASTRAL_JAR to astral.5.7.8.jar}"

if [ ! -s "$ASTRAL_JAR" ]; then
    echo "ERROR: ASTRAL jar not found" >&2
    exit 1
fi

mkdir -p "$OUT"

cat "$TREE_DIR"/*.treefile \
    > "$OUT/nuclear_gene_trees.tre"

java -jar "$ASTRAL_JAR" \
    -i "$OUT/nuclear_gene_trees.tre" \
    -o "$OUT/astral_species_tree.tre"