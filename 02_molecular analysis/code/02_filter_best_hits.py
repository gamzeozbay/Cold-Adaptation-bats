#!/usr/bin/env python3
"""
Filter BLASTP results and select protein accesion per species and gene.
Inputs:
- input/blast_results/*_mmyo.tsv
- config/gene_keep_patterns.txt
- config/rescue_accessions.tsv

Outputs:
- output/best_hits_by_gene/
- output/best_accessions_by_gene/
- output/best_hits_all_genes.tsv
"""

from pathlib import Path
from collections import defaultdict
import csv
import re
import math


# Resolve paths relative to this script:
# 02_molecular_analysis/
# ├── code/02_filter_best_hits.py
# ├── config/
# ├── input/blast_results/
# └── output/
SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent

BLAST_DIR = ROOT / "input" / "blast_results"
PATTERN_FILE = ROOT / "config" / "gene_keep_patterns.txt"
RESCUE_FILE = ROOT / "config" / "rescue_accessions.tsv"

OUT_DIR = ROOT / "output"
BEST_DIR = OUT_DIR / "best_hits_by_gene"
ACC_DIR = OUT_DIR / "best_accessions_by_gene"


# Sequence-quality thresholds used for orthologue filtering.
MIN_PIDENT = 70.0
MIN_QCOVS = 70.0
MIN_RATIO = 0.70
MAX_RATIO = 1.30


FIELDS = [
    "gene", "species", "sacc", "pident", "qcovs", "qlen", "slen", "ratio",
    "bitscore", "evalue", "predicted", "like", "partial", "low_quality",
    "manual_rescue", "failed_filters", "matched_pattern", "stitle", "staxids"
]


# Recreate only this step's outputs without removing later pipeline results.
for directory in (BEST_DIR, ACC_DIR):
    if directory.exists():
        for path in directory.iterdir():
            if path.is_file():
                path.unlink()
    else:
        directory.mkdir(parents=True, exist_ok=True)

summary_file = OUT_DIR / "best_hits_all_genes.tsv"
if summary_file.exists():
    summary_file.unlink()

OUT_DIR.mkdir(parents=True, exist_ok=True)
BEST_DIR.mkdir(parents=True, exist_ok=True)
ACC_DIR.mkdir(parents=True, exist_ok=True)


def load_patterns():
    patterns = {}

    for line in PATTERN_FILE.read_text().splitlines():
        line = line.strip()
        if not line:
            continue

        gene, pats = line.split("|||", 1)
        patterns[gene] = [
            p.strip().lower()
            for p in pats.split(";;")
            if p.strip()
        ]

    return patterns


def load_rescues():
    rescues = set()

    if not RESCUE_FILE.exists():
        return rescues

    with RESCUE_FILE.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rescues.add((row["gene"], row["accession"]))

    return rescues


def species_from_title(title):
    # RefSeq protein titles typically end with the species name in brackets.
    matches = re.findall(r"\[([^\[\]]+)\]", title)
    return matches[-1] if matches else "UNKNOWN_SPECIES"


def fnum(value):
    try:
        return float(value)
    except Exception:
        return math.nan


def title_match(title, patterns):
    title_lower = title.lower()

    for pattern in patterns:
        if pattern in title_lower:
            return pattern

    return ""


def is_like_title(gene, title):
    title_lower = title.lower()

    # For CREB3L1 and CREB3L4, "3-like protein 1" and "3-like protein 4"
    # are part of the expected RefSeq protein names and are not treated
    # as generic "-like" annotations.
    if gene in {"CREB3L1", "CREB3L4"}:
        official = {
            "CREB3L1": "cyclic amp-responsive element-binding protein 3-like protein 1",
            "CREB3L4": "cyclic amp-responsive element-binding protein 3-like protein 4",
        }[gene]
        title_lower = title_lower.replace(official, "")

    return ("-like" in title_lower) or (" like" in title_lower)


def fail_reasons(row):
    reasons = []

    if not row["matched_pattern"]:
        reasons.append("annotation_not_matched")

    if row["partial"]:
        reasons.append("partial")

    if row["low_quality"]:
        reasons.append("low_quality")

    if math.isnan(row["pident"]) or row["pident"] < MIN_PIDENT:
        reasons.append("pident_lt_70")

    if math.isnan(row["qcovs"]) or row["qcovs"] < MIN_QCOVS:
        reasons.append("qcovs_lt_70")

    if (
        math.isnan(row["ratio"])
        or row["ratio"] < MIN_RATIO
        or row["ratio"] > MAX_RATIO
    ):
        reasons.append("ratio_outside_0.70_1.30")

    return reasons


def parse_gene(gene, patterns, rescues):
    path = BLAST_DIR / f"{gene}_mmyo.tsv"
    rows = []

    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")

            # Ignore malformed BLAST rows that do not match the expected
            # 11-column outfmt used by the BLASTP script.
            if len(parts) != 11:
                continue

            (
                qseqid,
                sacc,
                pident,
                length,
                qlen,
                slen,
                qcovs,
                evalue,
                bitscore,
                stitle,
                staxids,
            ) = parts

            qlen_f = fnum(qlen)
            slen_f = fnum(slen)

            ratio = math.nan
            if not math.isnan(qlen_f) and qlen_f > 0 and not math.isnan(slen_f):
                ratio = slen_f / qlen_f

            title_lower = stitle.lower()

            row = {
                "gene": gene,
                "species": species_from_title(stitle),
                "sacc": sacc,
                "pident": fnum(pident),
                "qcovs": fnum(qcovs),
                "qlen": qlen_f,
                "slen": slen_f,
                "ratio": ratio,
                "bitscore": fnum(bitscore),
                "evalue": evalue,
                "predicted": "predicted:" in title_lower,
                "like": is_like_title(gene, stitle),
                "partial": "partial" in title_lower,
                "low_quality": "low quality" in title_lower,
                "manual_rescue": (gene, sacc) in rescues,
                "matched_pattern": title_match(stitle, patterns[gene]),
                "stitle": stitle,
                "staxids": staxids,
            }

            failed = fail_reasons(row)
            row["failed_filters"] = ";".join(failed) if failed else "PASS"

            # Retain manually reviewed rescue accessions despite failed
            # automated filters.
            if not failed or row["manual_rescue"]:
                rows.append(row)

    return rows


def select_key(row):
    # Select the longest isoform for each species.
    # Ties are resolved by bitscore, query coverage, identity, and
    # preference for a non-"like" annotation.
    non_like = 0 if row["like"] else 1

    return (
        row["slen"],
        row["bitscore"],
        row["qcovs"],
        row["pident"],
        non_like,
    )


def write_rows(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=FIELDS,
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()

        for row in rows:
            output_row = row.copy()

            for key in ["pident", "qcovs", "qlen", "slen", "ratio", "bitscore"]:
                value = output_row[key]
                if isinstance(value, float):
                    output_row[key] = (
                        "" if math.isnan(value) else f"{value:.4f}"
                    )

            for key in [
                "predicted",
                "like",
                "partial",
                "low_quality",
                "manual_rescue",
            ]:
                output_row[key] = "YES" if output_row[key] else "NO"

            writer.writerow(output_row)


patterns = load_patterns()
rescues = load_rescues()

all_best = []

for gene in sorted(patterns):
    rows = parse_gene(gene, patterns, rescues)

    # Group passing hits by species and retain one accession per species.
    by_species = defaultdict(list)
    for row in rows:
        by_species[row["species"]].append(row)

    best = []
    for species, hits in by_species.items():
        best.append(sorted(hits, key=select_key, reverse=True)[0])

    best = sorted(best, key=lambda row: row["species"])
    all_best.extend(best)

    write_rows(BEST_DIR / f"{gene}_best_hits.tsv", best)

    with (ACC_DIR / f"{gene}.accessions.txt").open("w") as handle:
        for row in best:
            handle.write(row["sacc"] + "\n")


write_rows(OUT_DIR / "best_hits_all_genes.tsv", all_best)

print(f"Processed genes: {len(patterns)}")
print(f"Selected protein accessions: {len(all_best)}")
print(f"Outputs written to: {OUT_DIR}")