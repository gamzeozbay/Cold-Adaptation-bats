#!/bin/bash
#SBATCH --job-name=blast_chiro_clean
#SBATCH --account=project_2004342
#SBATCH --partition=small
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --logs/blast_chiro_%A_%a.out
#SBATCH --error=/logs/blast_chiro_%A_%a.err

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

QUERY_LIST="$ROOT/config/query_list.txt"
OUT_DIR="$ROOT/output/blast"
LOG_DIR="$ROOT/logs"

mkdir -p "$OUT_DIR" "$LOG_DIR"
query=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$QUERY_LIST")
base=$(basename "$query" .fasta)
out="$OUT_DIR/${base}.tsv"
tmp="${out}.tmp.${SLURM_JOB_ID}.${SLURM_ARRAY_TASK_ID}"

#running blastp

blastp \
  -query "$query" \
  -db "$DB" \
  -taxids 9397 \
  -out "$tmp" \
  -outfmt "6 qseqid sacc pident length qlen slen qcovs evalue bitscore stitle staxids" \
  -evalue 1e-10 \
  -max_target_seqs 50000 \
  -num_threads 8


