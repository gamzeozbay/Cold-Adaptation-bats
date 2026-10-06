#!/bin/bash

# Retrieve protein and coding sequences for the final selected accessions.
#
# Input:
#   input/selected_by_gene/*.selected.tsv
#
# Outputs:
#   output/proteins_by_gene/*.faa
#   output/cds_by_gene/*.cds.fasta
#
# Required:
#   BLAST+ (blastdbcmd)
#   NCBI EDirect (efetch)
#
# Before running:
#   export BLAST_DB=/path/to/refseq_mammal_taxdb

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

IN="$ROOT/input/selected_by_gene"
PROTEIN_OUT="$ROOT/output/proteins_by_gene"
CDS_OUT="$ROOT/output/cds_by_gene"

DB="${BLAST_DB:?Please set BLAST_DB to the RefSeq mammalian BLAST database prefix}"

command -v blastdbcmd >/dev/null || {
    echo "ERROR: blastdbcmd not found" >&2
    exit 1
}

command -v efetch >/dev/null || {
    echo "ERROR: efetch not found" >&2
    exit 1
}

mkdir -p "$PROTEIN_OUT" "$CDS_OUT"

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

for hits in "$IN"/*.selected.tsv
do
    gene=$(basename "$hits" .selected.tsv)

    accfile="$WORK/${gene}.accessions.txt"
    raw_proteins="$WORK/${gene}.raw.faa"

    protein_file="$PROTEIN_OUT/${gene}.faa"
    cds_file="$CDS_OUT/${gene}.cds.fasta"

    tail -n +2 "$hits" | cut -f3 > "$accfile"

    # Retrieve the selected proteins from the local RefSeq database.
    blastdbcmd \
        -db "$DB" \
        -entry_batch "$accfile" \
        -out "$raw_proteins" \
        -outfmt "%f"

    # Rename protein FASTA headers as gene|species|accession.
    python3 - "$gene" "$hits" "$raw_proteins" "$protein_file" <<'PY'

import sys
import re

gene, selection_path, raw_path, output_path = sys.argv[1:]


def norm_acc(accession):
    accession = accession.split("|")[-1]
    return re.sub(r"\.\d+$", "", accession)


selected = []

with open(selection_path) as handle:
    next(handle)

    for line in handle:
        if not line.strip():
            continue

        gene_name, species, accession = line.rstrip("\n").split("\t")[:3]

        species = species.replace(" ", "_")

        selected.append(
            (
                norm_acc(accession),
                accession,
                f"{gene_name}|{species}|{accession}",
            )
        )


sequences = {}
current = None
buffer = []

with open(raw_path) as handle:
    for line in handle:
        line = line.rstrip("\n")

        if line.startswith(">"):
            if current is not None:
                sequences[current] = "".join(buffer)

            current = norm_acc(line[1:].split()[0])
            buffer = []

        else:
            buffer.append(line.strip())

if current is not None:
    sequences[current] = "".join(buffer)


with open(output_path, "w") as out:
    for accession_norm, accession, header in selected:

        if accession_norm not in sequences:
            raise RuntimeError(
                f"Protein not found: {gene} {accession}"
            )

        sequence = sequences[accession_norm]

        out.write(f">{header}\n")

        for i in range(0, len(sequence), 60):
            out.write(sequence[i:i + 60] + "\n")

PY

    : > "$cds_file"

    # Retrieve the CDS corresponding to each selected protein accession.
    tail -n +2 "$hits" |
    while IFS=$'\t' read -r gene_name species accession
    do
        [ -n "$accession" ] || continue

        species_clean=$(echo "$species" | tr ' ' '_')
        sequence=""

        for attempt in 1 2 3 4 5
        do
            sequence=$(
                efetch \
                    -db protein \
                    -id "$accession" \
                    -format fasta_cds_na \
                    2>/dev/null |
                awk '
                    /^>/ {
                        if (n++) exit
                        next
                    }
                    {
                        printf "%s", $0
                    }
                '
            )

            if [ -n "$sequence" ]; then
                break
            fi

            sleep $((attempt * 3))
        done

        if [ -z "$sequence" ]; then
            echo "ERROR: CDS not found: $gene_name $species_clean $accession" >&2
            exit 1
        fi

        echo ">${gene_name}|${species_clean}|${accession}" >> "$cds_file"

        echo "$sequence" |
        fold -w 60 >> "$cds_file"

        sleep 0.35

    done

done