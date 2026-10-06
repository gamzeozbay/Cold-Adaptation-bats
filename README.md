# Cold adaptation in mammals: mechanisms and genomic insights from bats: code and data

This repository contains the code and data package used for the systematic-review synthesis and the candidate-gene molecular evolutionary analyses reported in the manuscript.

## Repository structure

```text
Manuscript_codes/
├── README.md
├── supplementary_tables.xlsx
│
├── 01_systematic_review/
│   ├── code/
│   ├── input/
│   └── output/
│
└── 02_molecular_analysis/
    ├── code/
    ├── config/
    ├── input/
    └── output/
```

## 1. Systematic review

The systematic-review analysis starts from the final curated dataset of 137 included studies.

Screening and eligibility decisions were completed before this analysis stage and are documented in the manuscript and supplementary material.

## 2. Molecular analysis

The molecular workflow analyzes candidate genes associated with thermogenesis, metabolism, mitochondrial function, and cold adaptation in bats.

The final nuclear codon-model dataset contains **33 genes**.

Three mitochondrial genes (**ATP6, ND3, and ND4**) were processed during sequence recovery and alignment but were excluded from the nuclear PAML analyses.

### Workflow

```text
query proteins
    ↓
BLASTP against Chiroptera
    ↓
orthologue filtering and accession selection
    ↓
manual final QC / locked accession corrections
    ↓
protein and CDS retrieval
    ↓
reciprocal best-hit validation
    ↓
MAFFT protein alignment
    ↓
PAL2NAL codon alignment
    ↓
Gblocks codon trimming
    ↓
TAPER sequence-specific masking
    ↓
IQ-TREE gene trees
    ↓
ASTRAL species tree
    ↓
PAML codon models
    ↓
likelihood-ratio tests, FDR correction and BEB parsing
    ↓
final figures
```

### Molecular input files

```text
02_molecular_analysis/input/
├── query_proteins/
├── blast_results/
└── selected_by_gene/
```

The large local RefSeq mammalian BLAST database is not distributed with this repository.

To rerun sequence retrieval, set the `BLAST_DB` environment variable to the corresponding local BLAST database prefix.

### Molecular configuration files

```text
02_molecular_analysis/config/
├── gene_keep_patterns.txt
├── rescue_accessions.tsv
└── final_qc_actions.tsv
```

`gene_keep_patterns.txt` contains gene-specific annotation patterns used during orthologue filtering.

`rescue_accessions.tsv` contains predefined accession-level rescue rules used during filtering.

`final_qc_actions.tsv` records the locked manual accession replacements and taxon removals applied before construction of the final alignments.

The scripts use relative paths within `02_molecular_analysis/` wherever possible. External executables or databases are supplied through the local software environment or environment variables.

## PAML output

The distributed PAML directory contains the analysis inputs, control files, labelled trees, raw CODEML outputs, and job definitions used in the final analysis.

```text
02_molecular_analysis/output/paml/
├── alignments/
├── trees/
├── ctl/
├── results/
├── job_lists/
└── tree_label_report.tsv
```

Parsed result tables can be regenerated with:

```text
13_parse_paml_results.py
```

The parser writes:

```text
paml_model_parameters.tsv
paml_master_results.tsv
paml_beb_sites.tsv
```

## Other distributed molecular outputs

```text
02_molecular_analysis/output/
├── gblocks_alignments/
├── gene_trees/
├── species_tree/
└── paml/
```

The package does not include every historical intermediate file.

Protein FASTA files, CDS FASTA files, protein alignments, and pre-Gblocks codon alignments can be regenerated from the locked accession tables with the supplied scripts.

## Main software

The final workflow used the following principal software:

- BLAST+ 2.17
- NCBI EDirect 13.4
- MAFFT 7.505
- PAL2NAL 14
- Gblocks 0.91b
- TAPER 1.02
- IQ-TREE2 2.3.6
- ASTRAL 5.7.8
- PAML codeml 4.9
- R
- Python 3

Software-specific paths are not hard-coded in the reviewer-facing package.

Where needed, executables or resources can be provided through the user's environment, including variables such as `BLAST_DB`, `TAPER_DIR`, and `ASTRAL_JAR`.

## Reproducibility note

This repository is intended to reproduce the final manuscript analysis rather than preserve the complete development history of the project.

The authoritative starting points are:

- the curated systematic-review dataset,
- the raw BLASTP tables,
- the locked final accession tables, and
- the final configuration/QC decision files.

Superseded pipeline versions, exploratory scripts, temporary cluster logs, sensitivity-test directories that are not part of the reported results, and obsolete result tables were intentionally excluded.

For exact biological interpretation, sample counts, supplementary table definitions, and statistical reporting, refer to the manuscript and accompanying supplementary material.
