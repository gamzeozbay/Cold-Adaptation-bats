#!/usr/bin/env python3

# Prepare nuclear PAML alignments, foreground trees, control files and job list.

from pathlib import Path
import os

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent

ALIGNMENT_DIR = ROOT / "output" / "taper_alignments"
SPECIES_TREE = ROOT / "output" / "species_tree" / "astral_species_tree.tre"
PAML_DIR = ROOT / "output" / "paml"

ALIGN_DIR = PAML_DIR / "alignments"
TREE_SITE_DIR = PAML_DIR / "trees" / "site"
TREE_LABEL_DIR = PAML_DIR / "trees" / "labelled"
CTL_DIR = PAML_DIR / "ctl"
RESULT_DIR = PAML_DIR / "results"
JOB_DIR = PAML_DIR / "job_lists"

for directory in [
    ALIGN_DIR,
    TREE_SITE_DIR,
    TREE_LABEL_DIR,
    CTL_DIR,
    RESULT_DIR,
    JOB_DIR,
]:
    directory.mkdir(parents=True, exist_ok=True)

MITO = {"ATP6", "ND3", "ND4"}

MYOTIS7 = {
    "Myotis_lucifugus",
    "Myotis_yumanensis",
    "Myotis_brandtii",
    "Myotis_davidii",
    "Myotis_daubentonii",
    "Myotis_myotis",
    "Myotis_nattereri",
}

EPTESICUS2 = {
    "Cnephaeus_nilssonii",
    "Eptesicus_fuscus",
}

VESPERTILIONIDAE = MYOTIS7 | {
    "Pipistrellus_kuhlii",
    "Cnephaeus_nilssonii",
    "Eptesicus_fuscus",
    "Plecotus_auritus",
}

HYPOTHESES = {
    "nilssonii_terminal": {
        "mode": "terminal",
        "taxa": {"Cnephaeus_nilssonii"},
    },
    "ept_stem": {
        "mode": "stem",
        "taxa": EPTESICUS2,
    },
    "ept_terminals": {
        "mode": "terminal",
        "taxa": EPTESICUS2,
    },
    "myotis_stem": {
        "mode": "stem",
        "taxa": MYOTIS7,
    },
    "myotis_all7_terminals": {
        "mode": "terminal",
        "taxa": MYOTIS7,
    },
    "myotis_wholeclade": {
        "mode": "stem_internal",
        "taxa": MYOTIS7,
    },
    "vespertilionidae_stem": {
        "mode": "stem",
        "taxa": VESPERTILIONIDAE,
    },
}


class Node:
    def __init__(self, name=None, children=None):
        self.name = name
        self.children = children or []
        self.mark = False

    def is_leaf(self):
        return not self.children

    def leaves(self):
        if self.is_leaf():
            return {self.name}

        result = set()
        for child in self.children:
            result |= child.leaves()
        return result

    def to_newick(self):
        label = "#1" if self.mark else ""

        if self.is_leaf():
            return f"{self.name}{label}"

        children = ",".join(child.to_newick() for child in self.children)
        return f"({children}){label}"


def parse_newick(text):
    text = "".join(text.split()).rstrip(";")

    def read_token(index):
        start = index
        while index < len(text) and text[index] not in ",():;":
            index += 1
        return text[start:index], index

    def skip_branch_length(index):
        if index < len(text) and text[index] == ":":
            index += 1
            while index < len(text) and text[index] not in ",();":
                index += 1
        return index

    def parse_subtree(index):
        if text[index] == "(":
            index += 1
            children = []

            while True:
                child, index = parse_subtree(index)
                children.append(child)

                if index < len(text) and text[index] == ",":
                    index += 1
                    continue

                if index < len(text) and text[index] == ")":
                    index += 1
                    break

            _, index = read_token(index)
            index = skip_branch_length(index)
            return Node(children=children), index

        name, index = read_token(index)
        index = skip_branch_length(index)
        return Node(name=name), index

    root, _ = parse_subtree(0)
    return root


def prune(node, keep):
    if node.is_leaf():
        return Node(name=node.name) if node.name in keep else None

    children = []

    for child in node.children:
        pruned = prune(child, keep)
        if pruned is not None:
            children.append(pruned)

    if not children:
        return None

    if len(children) == 1:
        return children[0]

    return Node(children=children)


def mark_terminals(node, targets):
    if node.is_leaf():
        if node.name in targets:
            node.mark = True
            return 1
        return 0

    return sum(mark_terminals(child, targets) for child in node.children)


def find_mrca(node, targets):
    if not targets <= node.leaves():
        return None

    for child in node.children:
        candidate = find_mrca(child, targets)
        if candidate is not None:
            return candidate

    return node


def exact_clade(node, targets):
    candidate = find_mrca(node, targets)

    if candidate is None:
        return None

    if candidate.leaves() != targets:
        return None

    return candidate


def mark_stem(node, targets):
    candidate = exact_clade(node, targets)

    if candidate is None:
        return "NOT_MONOPHYLETIC_AFTER_PRUNING", 0

    candidate.mark = True
    return "STEM_LABELLED", 1


def mark_internal_descendants(node):
    count = 0

    for child in node.children:
        if child.is_leaf():
            continue

        child.mark = True
        count += 1
        count += mark_internal_descendants(child)

    return count


def mark_stem_and_internal(node, targets):
    candidate = exact_clade(node, targets)

    if candidate is None:
        return "NOT_MONOPHYLETIC_AFTER_PRUNING", 0

    candidate.mark = True
    count = 1 + mark_internal_descendants(candidate)

    return "STEM_INTERNAL_LABELLED", count


def read_fasta(path):
    records = {}
    header = None
    sequence = []

    with path.open() as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(sequence).upper()

                header = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)

    if header is not None:
        records[header] = "".join(sequence).upper()

    return records


def species_records(path):
    raw = read_fasta(path)

    if not raw:
        raise RuntimeError(f"Empty alignment: {path}")

    records = {}

    for header, sequence in raw.items():
        parts = header.split("|")

        if len(parts) < 2:
            raise RuntimeError(
                f"Expected gene|species|accession FASTA header: {header}"
            )

        species = parts[1]

        if species in records:
            raise RuntimeError(
                f"Duplicate species in {path.name}: {species}"
            )

        records[species] = sequence

    lengths = {len(sequence) for sequence in records.values()}

    if len(lengths) != 1:
        raise RuntimeError(f"Unequal sequence lengths in {path}")

    length = lengths.pop()

    if length % 3 != 0:
        raise RuntimeError(
            f"Alignment length is not divisible by three: {path}"
        )

    return records


def write_paml_alignment(records, path):
    length = len(next(iter(records.values())))

    with path.open("w") as handle:
        handle.write(f"{len(records)} {length}\n")

        for species, sequence in records.items():
            handle.write(f"{species}\n{sequence}\n")


def relative_path(path, start):
    return os.path.relpath(path, start=start)


def ctl_text(
    seqfile,
    treefile,
    outfile,
    model,
    nssites,
    fix_omega,
    omega,
):
    return f"""seqfile = {seqfile}
treefile = {treefile}
outfile = {outfile}

noisy = 0
verbose = 1
runmode = 0

seqtype = 1
CodonFreq = 2
clock = 0
aaDist = 0

model = {model}
NSsites = {nssites}

icode = 0
Mgene = 0

fix_kappa = 0
kappa = 2

fix_omega = {fix_omega}
omega = {omega}

fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 4

getSE = 0
RateAncestor = 0
Small_Diff = 0.5e-6
cleandata = 1
fix_blength = 0
method = 0
"""


if not SPECIES_TREE.is_file():
    raise RuntimeError(f"Species tree not found: {SPECIES_TREE}")

species_root = parse_newick(SPECIES_TREE.read_text())
all_species = species_root.leaves()

site_models = {
    "M0": {
        "model": 0,
        "NSsites": 0,
        "fix_omega": 0,
        "omega": 0.4,
    },
    "M1a": {
        "model": 0,
        "NSsites": 1,
        "fix_omega": 0,
        "omega": 0.4,
    },
    "M2a": {
        "model": 0,
        "NSsites": 2,
        "fix_omega": 0,
        "omega": 1.5,
    },
    "M7": {
        "model": 0,
        "NSsites": 7,
        "fix_omega": 0,
        "omega": 0.4,
    },
    "M8": {
        "model": 0,
        "NSsites": 8,
        "fix_omega": 0,
        "omega": 1.5,
    },
    "M8a": {
        "model": 0,
        "NSsites": 8,
        "fix_omega": 1,
        "omega": 1.0,
    },
}

jobs = []
tree_report = PAML_DIR / "tree_label_report.tsv"

with tree_report.open("w") as report:
    report.write(
        "gene\thypothesis\tmode\tstatus\t"
        "n_alignment_taxa\tn_requested\tn_present\t"
        "foreground_present\tforeground_missing\tn_labelled\n"
    )

    for alignment in sorted(
        ALIGNMENT_DIR.glob("*.codon.gblocks.taper.fasta")
    ):
        gene = alignment.name.replace(
            ".codon.gblocks.taper.fasta",
            "",
        )

        if gene in MITO:
            continue

        records = species_records(alignment)
        taxa = set(records)

        missing_from_tree = taxa - all_species

        if missing_from_tree:
            raise RuntimeError(
                f"{gene}: taxa missing from species tree: "
                f"{sorted(missing_from_tree)}"
            )

        seqfile = ALIGN_DIR / f"{gene}.phy"
        write_paml_alignment(records, seqfile)

        pruned = prune(species_root, taxa)

        if pruned is None:
            raise RuntimeError(f"Empty pruned tree for {gene}")

        site_tree_dir = TREE_SITE_DIR / gene
        site_tree_dir.mkdir(parents=True, exist_ok=True)

        site_tree = (
            site_tree_dir /
            f"{gene}.pruned.unlabelled.tree"
        )
        site_tree.write_text(pruned.to_newick() + ";\n")

        for model_name, parameters in site_models.items():
            run_dir = RESULT_DIR / gene / "site_models" / model_name
            ctl_dir = CTL_DIR / gene / "site_models"

            run_dir.mkdir(parents=True, exist_ok=True)
            ctl_dir.mkdir(parents=True, exist_ok=True)

            ctl = ctl_dir / f"{gene}.{model_name}.ctl"
            outfile = f"{gene}.{model_name}.out"

            ctl.write_text(
                ctl_text(
                    seqfile=relative_path(seqfile, run_dir),
                    treefile=relative_path(site_tree, run_dir),
                    outfile=outfile,
                    model=parameters["model"],
                    nssites=parameters["NSsites"],
                    fix_omega=parameters["fix_omega"],
                    omega=parameters["omega"],
                )
            )

            jobs.append(
                (
                    gene,
                    "site",
                    model_name,
                    relative_path(run_dir, ROOT),
                    relative_path(ctl, ROOT),
                )
            )

        for hypothesis, config in HYPOTHESES.items():
            requested = set(config["taxa"])
            present = requested & taxa
            missing = requested - taxa

            if not present:
                report.write(
                    f"{gene}\t{hypothesis}\t{config['mode']}\t"
                    f"SKIPPED_NO_FOREGROUND\t{len(taxa)}\t"
                    f"{len(requested)}\t0\tNA\t"
                    f"{','.join(sorted(missing))}\t0\n"
                )
                continue

            tree = prune(species_root, taxa)

            if config["mode"] == "terminal":
                n_labelled = mark_terminals(tree, present)
                status = (
                    "TERMINAL_LABELLED"
                    if n_labelled
                    else "SKIPPED_ZERO_LABELS"
                )

            elif config["mode"] == "stem":
                status, n_labelled = mark_stem(tree, present)

            elif config["mode"] == "stem_internal":
                status, n_labelled = mark_stem_and_internal(
                    tree,
                    present,
                )

            else:
                raise RuntimeError(
                    f"Unknown foreground mode: {config['mode']}"
                )

            if (
                status.startswith("SKIPPED")
                or status == "NOT_MONOPHYLETIC_AFTER_PRUNING"
            ):
                report.write(
                    f"{gene}\t{hypothesis}\t{config['mode']}\t"
                    f"{status}\t{len(taxa)}\t{len(requested)}\t"
                    f"{len(present)}\t"
                    f"{','.join(sorted(present))}\t"
                    f"{','.join(sorted(missing))}\t0\n"
                )
                continue

            labelled_dir = TREE_LABEL_DIR / gene / hypothesis
            labelled_dir.mkdir(parents=True, exist_ok=True)

            labelled_tree = (
                labelled_dir /
                f"{gene}.{hypothesis}.labelled.tree"
            )
            labelled_tree.write_text(tree.to_newick() + ";\n")

            report.write(
                f"{gene}\t{hypothesis}\t{config['mode']}\t"
                f"{status}\t{len(taxa)}\t{len(requested)}\t"
                f"{len(present)}\t"
                f"{','.join(sorted(present))}\t"
                f"{','.join(sorted(missing))}\t"
                f"{n_labelled}\n"
            )

            for model_name, fix_omega, omega in [
                ("BS_null", 1, 1.0),
                ("BS_alt", 0, 1.5),
            ]:
                run_dir = (
                    RESULT_DIR /
                    gene /
                    "branchsite" /
                    hypothesis /
                    model_name
                )
                ctl_dir = (
                    CTL_DIR /
                    gene /
                    "branchsite" /
                    hypothesis
                )

                run_dir.mkdir(parents=True, exist_ok=True)
                ctl_dir.mkdir(parents=True, exist_ok=True)

                ctl = (
                    ctl_dir /
                    f"{gene}.{hypothesis}.{model_name}.ctl"
                )
                outfile = (
                    f"{gene}.{hypothesis}.{model_name}.out"
                )

                ctl.write_text(
                    ctl_text(
                        seqfile=relative_path(seqfile, run_dir),
                        treefile=relative_path(
                            labelled_tree,
                            run_dir,
                        ),
                        outfile=outfile,
                        model=2,
                        nssites=2,
                        fix_omega=fix_omega,
                        omega=omega,
                    )
                )

                jobs.append(
                    (
                        gene,
                        "branchsite",
                        f"{hypothesis}_{model_name}",
                        relative_path(run_dir, ROOT),
                        relative_path(ctl, ROOT),
                    )
                )

            run_dir = (
                RESULT_DIR /
                gene /
                "branch_model" /
                hypothesis /
                "two_ratio"
            )
            ctl_dir = (
                CTL_DIR /
                gene /
                "branch_model" /
                hypothesis
            )

            run_dir.mkdir(parents=True, exist_ok=True)
            ctl_dir.mkdir(parents=True, exist_ok=True)

            ctl = (
                ctl_dir /
                f"{gene}.{hypothesis}.two_ratio.ctl"
            )
            outfile = f"{gene}.{hypothesis}.two_ratio.out"

            ctl.write_text(
                ctl_text(
                    seqfile=relative_path(seqfile, run_dir),
                    treefile=relative_path(
                        labelled_tree,
                        run_dir,
                    ),
                    outfile=outfile,
                    model=2,
                    nssites=0,
                    fix_omega=0,
                    omega=0.4,
                )
            )

            jobs.append(
                (
                    gene,
                    "branch",
                    f"{hypothesis}_two_ratio",
                    relative_path(run_dir, ROOT),
                    relative_path(ctl, ROOT),
                )
            )

job_file = JOB_DIR / "all_paml_jobs.tsv"

with job_file.open("w") as handle:
    handle.write(
        "gene\tfamily\ttest\trun_dir\tctl\n"
    )

    for row in jobs:
        handle.write("\t".join(row) + "\n")