#!/bin/bash

#SBATCH --job-name=rbh_full
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

# Validate selected bat proteins by reciprocal BLASTP searches
# against the complete Myotis myotis proteome.
#
# Inputs:
#   output/proteins_by_gene/*.faa
#   config/gene_keep_patterns.txt
#   local RefSeq mammalian BLAST database
#
# Outputs:
#   output/rbh_full/rbh_full_hits.tsv
#   output/rbh_full/rbh_full_summary.tsv
#   output/rbh_full/myotis_proteome_size.tsv
#
# Before running:
#   export BLAST_DB=/path/to/refseq_mammal_taxdb

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

PROT="$ROOT/output/proteins_by_gene"
PATTERN_FILE="$ROOT/config/gene_keep_patterns.txt"
OUT="$ROOT/output/rbh_full"

DB="${BLAST_DB:?Please set BLAST_DB to the RefSeq mammalian BLAST database prefix}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "$OUT"

# Use a temporary directory for files that do not need to be retained.
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# Combine all final selected protein sequences.
cat "$PROT"/*.faa > "$WORK/final_proteins.faa"

# Extract the complete Myotis myotis proteome from the local RefSeq database.
blastdbcmd \
    -db "$DB" \
    -taxids 51298 \
    -outfmt "%f" \
    > "$WORK/myotis_proteome.faa"

N=$(grep -c '^>' "$WORK/myotis_proteome.faa" || true)

if [[ "$N" -lt 1000 ]]; then
    echo "ERROR: unexpectedly few Myotis myotis proteins: $N" >&2
    exit 1
fi

# Record the number of proteins in the extracted M. myotis proteome.
printf "taxid\tn_proteins\n51298\t%s\n" "$N" \
    > "$OUT/myotis_proteome_size.tsv"

# Build a temporary BLAST database from the M. myotis proteome.
makeblastdb \
    -in "$WORK/myotis_proteome.faa" \
    -dbtype prot \
    -out "$WORK/myotis_db" \
    >/dev/null

# Search all selected proteins against the complete M. myotis proteome.
blastp \
    -query "$WORK/final_proteins.faa" \
    -db "$WORK/myotis_db" \
    -outfmt "6 qseqid sseqid pident bitscore evalue stitle" \
    -max_target_seqs 5 \
    -evalue 1e-5 \
    -num_threads "$THREADS" \
    > "$OUT/rbh_full_hits.tsv"

# Identify the best M. myotis hit and compare its annotation
# with the expected target-gene annotation.
python3 - "$ROOT" "$OUT" <<'PY'

import sys
import csv
from pathlib import Path

ROOT = Path(sys.argv[1])
OUT = Path(sys.argv[2])

PATTERN_FILE = ROOT / "config" / "gene_keep_patterns.txt"
HITS_FILE = OUT / "rbh_full_hits.tsv"
SUMMARY_FILE = OUT / "rbh_full_summary.tsv"

MITOCHONDRIAL_GENES = {"ATP6", "ND3", "ND4"}

patterns = {}

with PATTERN_FILE.open() as handle:
    for line in handle:
        line = line.strip()

        if "|||" not in line:
            continue

        gene, pattern_string = line.split("|||", 1)

        patterns[gene] = [
            pattern.strip().lower()
            for pattern in pattern_string.split(";;")
            if pattern.strip()
        ]

best_hits = {}

with HITS_FILE.open() as handle:
    for line in handle:
        columns = line.rstrip("\n").split("\t")

        if len(columns) < 6:
            continue

        query = columns[0]
        subject = columns[1]
        bitscore = float(columns[3])
        title = columns[5]

        if query not in best_hits or bitscore > best_hits[query][1]:
            best_hits[query] = (
                subject,
                bitscore,
                title,
            )

rows = []

for query, (subject, bitscore, title) in sorted(best_hits.items()):
    gene = query.split("|")[0]

    # The reported RBH result refers to the nuclear candidate-gene set.
    if gene in MITOCHONDRIAL_GENES:
        continue

    expected_patterns = patterns.get(gene, [])

    correct_gene = any(
        pattern in title.lower()
        for pattern in expected_patterns
    )

    status = "RBH_OK" if correct_gene else "RBH_MISMATCH"

    rows.append([
        gene,
        query,
        subject,
        status,
        title,
    ])

with SUMMARY_FILE.open("w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")

    writer.writerow([
        "gene",
        "query",
        "best_myotis_hit",
        "status",
        "hit_title",
    ])

    writer.writerows(rows)

n_total = len(rows)
n_ok = sum(row[3] == "RBH_OK" for row in rows)
n_mismatch = n_total - n_ok

print(f"Nuclear queries: {n_total}")
print(f"RBH validated: {n_ok}")
print(f"RBH mismatches: {n_mismatch}")

PY