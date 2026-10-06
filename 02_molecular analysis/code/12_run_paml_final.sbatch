#!/bin/bash

#SBATCH --job-name=paml
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=48G

# Run all prepared CODEML analyses in parallel.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

JOBLIST="$ROOT/output/paml/job_lists/all_paml_jobs.tsv"
MAX_PARALLEL="${SLURM_CPUS_PER_TASK:-24}"

CODEML="${CODEML:-$(command -v codeml || true)}"

if [ -z "$CODEML" ]; then
    echo "ERROR: codeml not found" >&2
    exit 1
fi

if [ ! -s "$JOBLIST" ]; then
    echo "ERROR: PAML job list not found" >&2
    exit 1
fi

N=$(( $(wc -l < "$JOBLIST") - 1 ))

if [ "$N" -le 0 ]; then
    echo "ERROR: PAML job list is empty" >&2
    exit 1
fi

run_one() {
    idx="$1"

    line=$(sed -n "$((idx + 1))p" "$JOBLIST")
    IFS=$'\t' read -r gene family test run_dir ctl <<< "$line"

    run_dir="$ROOT/$run_dir"
    ctl="$ROOT/$ctl"

    if [ ! -s "$ctl" ]; then
        echo "ERROR: missing CTL for $gene $family $test" >&2
        return 1
    fi

    mkdir -p "$run_dir"
    cp "$ctl" "$run_dir/codeml.ctl"

    (
        cd "$run_dir"
        "$CODEML" codeml.ctl > codeml.stdout 2> codeml.stderr
    )
}

export -f run_one
export ROOT JOBLIST CODEML

seq 1 "$N" | xargs -n 1 -P "$MAX_PARALLEL" bash -c 'run_one "$1"' _
