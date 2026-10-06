#!/usr/bin/env python3

# Parse completed CODEML outputs and build final PAML result tables.

from pathlib import Path
import csv
import math
import re

from statsmodels.stats.multitest import multipletests

SCRIPT_DIR = Path(__file__).resolve().parent
ROOT = SCRIPT_DIR.parent

PAML_DIR = ROOT / "output" / "paml"
JOBLIST = PAML_DIR / "job_lists" / "all_paml_jobs.tsv"

OUT_MODELS = PAML_DIR / "paml_model_parameters.tsv"
OUT_RESULTS = PAML_DIR / "paml_master_results.tsv"
OUT_BEB = PAML_DIR / "paml_beb_sites.tsv"

FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"

LNL_RE = re.compile(
    r"lnL\(.*?np:\s*(\d+).*?\):\s*(" + FLOAT + r")"
)
KAPPA_RE = re.compile(
    r"kappa\s*\(ts/tv\)\s*=\s*(" + FLOAT + r")"
)
OMEGA_RE = re.compile(
    r"omega\s*\(dN/dS\)\s*=\s*(" + FLOAT + r")"
)
BEB_SITE_RE = re.compile(
    r"^\s*(\d+)\s+([A-Za-z*])\s+"
    r"(0?\.\d+|1\.0+)(\*{0,2})"
    r"(?:\s+(" + FLOAT + r")\s*\+\-\s*(" + FLOAT + r"))?",
    flags=re.MULTILINE,
)


def read_ctl_outfile(ctl):
    with ctl.open() as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith("outfile") and "=" in line:
                return line.split("=", 1)[1].strip()

    raise RuntimeError(f"outfile was not found in {ctl}")


def extract_section(text, start_pattern, end_patterns):
    match = re.search(start_pattern, text, flags=re.IGNORECASE)

    if not match:
        return ""

    start = match.start()
    end = len(text)

    for pattern in end_patterns:
        next_match = re.search(
            pattern,
            text[match.end():],
            flags=re.IGNORECASE,
        )
        if next_match:
            candidate = match.end() + next_match.start()
            end = min(end, candidate)

    return text[start:end]


def parse_beb(section):
    sites = []

    for match in BEB_SITE_RE.finditer(section):
        probability = float(match.group(3))

        sites.append({
            "site": int(match.group(1)),
            "aa": match.group(2),
            "probability": probability,
            "stars": match.group(4) or "",
            "postmean_omega": match.group(5) or "NA",
            "se_omega": match.group(6) or "NA",
        })

    return sites


def parse_output(path):
    result = {
        "status": "MISSING",
        "lnL": None,
        "np": None,
        "kappa": None,
        "omega_main": None,
        "p_values": "",
        "w_values": "",
        "w_max": None,
        "background_w": "",
        "background_w_max": None,
        "foreground_w": "",
        "foreground_w_max": None,
        "branch_w": "",
        "branch_w_max": None,
    }

    if not path.is_file() or path.stat().st_size == 0:
        return result, []

    text = path.read_text(errors="replace")

    lnl_matches = LNL_RE.findall(text)

    if not lnl_matches:
        result["status"] = "PARTIAL"
        return result, []

    np_value, lnl_value = lnl_matches[-1]
    result["np"] = int(np_value)
    result["lnL"] = float(lnl_value)
    result["status"] = "DONE"

    kappa_matches = KAPPA_RE.findall(text)
    if kappa_matches:
        result["kappa"] = float(kappa_matches[-1])

    omega_matches = OMEGA_RE.findall(text)
    if omega_matches:
        result["omega_main"] = float(omega_matches[-1])

    p_values = []
    w_values = []

    for line in text.splitlines():
        if re.match(r"^\s*p\s*:", line):
            p_values = re.findall(FLOAT, line)

        elif re.match(r"^\s*w\s*:", line):
            w_values = re.findall(FLOAT, line)

        elif re.match(r"^\s*background\s+w", line, flags=re.IGNORECASE):
            values = re.findall(FLOAT, line)
            result["background_w"] = ",".join(values)
            if values:
                result["background_w_max"] = max(map(float, values))

        elif re.match(r"^\s*foreground\s+w", line, flags=re.IGNORECASE):
            values = re.findall(FLOAT, line)
            result["foreground_w"] = ",".join(values)
            if values:
                result["foreground_w_max"] = max(map(float, values))

        elif "w (dN/dS) for branches" in line:
            values = re.findall(FLOAT, line)
            result["branch_w"] = ",".join(values)
            if values:
                result["branch_w_max"] = max(map(float, values))

    result["p_values"] = ",".join(p_values)
    result["w_values"] = ",".join(w_values)

    if w_values:
        result["w_max"] = max(map(float, w_values))

    beb_section = extract_section(
        text,
        r"Bayes Empirical Bayes|BEB",
        [
            r"Naive Empirical Bayes",
            r"NEB",
            r"Time used",
            r"The grid",
        ],
    )

    return result, parse_beb(beb_section)


def chi2_sf_df1(value):
    if value <= 0:
        return 1.0
    return math.erfc(math.sqrt(value / 2.0))


def chi2_sf_df2(value):
    if value <= 0:
        return 1.0
    return math.exp(-value / 2.0)


def lrt_pvalue(lrt, distribution):
    if distribution == "chi2_1":
        return chi2_sf_df1(lrt)

    if distribution == "chi2_2":
        return chi2_sf_df2(lrt)

    if distribution == "mix_1":
        if lrt <= 0:
            return 1.0
        return 0.5 * chi2_sf_df1(lrt)

    raise RuntimeError(f"Unknown null distribution: {distribution}")


def format_value(value):
    if value is None:
        return "NA"

    if isinstance(value, float):
        return f"{value:.12g}"

    return str(value)


if not JOBLIST.is_file():
    raise RuntimeError(f"Job list not found: {JOBLIST}")

jobs = []

with JOBLIST.open() as handle:
    reader = csv.DictReader(handle, delimiter="\t")

    required = {"gene", "family", "test", "run_dir", "ctl"}

    if set(reader.fieldnames or []) != required:
        raise RuntimeError(
            "Unexpected job-list columns. Expected: "
            "gene, family, test, run_dir, ctl"
        )

    for row in reader:
        run_dir = ROOT / row["run_dir"]
        ctl = ROOT / row["ctl"]

        if not ctl.is_file():
            raise RuntimeError(f"CTL not found: {ctl}")

        outfile_name = read_ctl_outfile(ctl)
        outfile = run_dir / outfile_name

        parsed, beb_sites = parse_output(outfile)

        info = {
            "gene": row["gene"],
            "family": row["family"],
            "test": row["test"],
            "run_dir": row["run_dir"],
            "ctl": row["ctl"],
            "outfile": str(outfile.relative_to(ROOT)),
            **parsed,
        }

        jobs.append((info, beb_sites))

runs = {
    (info["gene"], info["family"], info["test"]): info
    for info, _ in jobs
}

model_fields = [
    "gene",
    "family",
    "test",
    "status",
    "lnL",
    "np",
    "kappa",
    "omega_main",
    "p_values",
    "w_values",
    "w_max",
    "background_w",
    "background_w_max",
    "foreground_w",
    "foreground_w_max",
    "branch_w",
    "branch_w_max",
    "outfile",
    "ctl",
]

with OUT_MODELS.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=model_fields,
        delimiter="\t",
        extrasaction="ignore",
    )
    writer.writeheader()

    for info, _ in jobs:
        writer.writerow({
            key: format_value(info.get(key))
            for key in model_fields
        })

beb_rows = []

for info, sites in jobs:
    for site in sites:
        beb_rows.append({
            "gene": info["gene"],
            "family": info["family"],
            "test": info["test"],
            "site": site["site"],
            "aa": site["aa"],
            "probability": site["probability"],
            "stars": site["stars"],
            "postmean_omega": site["postmean_omega"],
            "se_omega": site["se_omega"],
            "outfile": info["outfile"],
        })

with OUT_BEB.open("w", newline="") as handle:
    fields = [
        "gene",
        "family",
        "test",
        "site",
        "aa",
        "probability",
        "stars",
        "postmean_omega",
        "se_omega",
        "outfile",
    ]

    writer = csv.DictWriter(
        handle,
        fieldnames=fields,
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(beb_rows)

results = []


def add_lrt(
    gene,
    comparison,
    test_family,
    null_key,
    alternative_key,
    distribution,
    fdr_group,
):
    null = runs.get(null_key)
    alternative = runs.get(alternative_key)

    if null is None or alternative is None:
        return

    status = (
        "DONE"
        if null["status"] == "DONE" and alternative["status"] == "DONE"
        else "INCOMPLETE"
    )

    if status == "DONE":
        raw_lrt = 2.0 * (alternative["lnL"] - null["lnL"])
        lrt = max(0.0, raw_lrt)
        p_value = lrt_pvalue(lrt, distribution)
    else:
        raw_lrt = None
        lrt = None
        p_value = None

    results.append({
        "gene": gene,
        "comparison": comparison,
        "test_family": test_family,
        "null_model": null_key[2],
        "alternative_model": alternative_key[2],
        "status": status,
        "lnL_null": null["lnL"],
        "lnL_alt": alternative["lnL"],
        "np_null": null["np"],
        "np_alt": alternative["np"],
        "raw_LRT": raw_lrt,
        "LRT": lrt,
        "null_distribution": distribution,
        "p_value": p_value,
        "fdr_group": fdr_group,
        "q_value": None,
        "omega_background": alternative["background_w_max"],
        "omega_foreground": alternative["foreground_w_max"],
        "omega_alt_max": alternative["w_max"],
        "null_outfile": null["outfile"],
        "alternative_outfile": alternative["outfile"],
    })


genes = sorted({key[0] for key in runs})

for gene in genes:
    add_lrt(
        gene,
        "M1a_vs_M2a",
        "site",
        (gene, "site", "M1a"),
        (gene, "site", "M2a"),
        "chi2_2",
        "site_M1a_vs_M2a",
    )

    add_lrt(
        gene,
        "M7_vs_M8",
        "site",
        (gene, "site", "M7"),
        (gene, "site", "M8"),
        "chi2_2",
        "site_M7_vs_M8",
    )

    add_lrt(
        gene,
        "M8a_vs_M8",
        "site",
        (gene, "site", "M8a"),
        (gene, "site", "M8"),
        "mix_1",
        "site_M8a_vs_M8",
    )

    branchsite_alternatives = sorted(
        test
        for g, family, test in runs
        if (
            g == gene
            and family == "branchsite"
            and test.endswith("_BS_alt")
        )
    )

    for alternative_test in branchsite_alternatives:
        hypothesis = alternative_test.removesuffix("_BS_alt")

        add_lrt(
            gene,
            hypothesis,
            "branchsite",
            (
                gene,
                "branchsite",
                f"{hypothesis}_BS_null",
            ),
            (
                gene,
                "branchsite",
                alternative_test,
            ),
            "mix_1",
            "branchsite",
        )

    branch_alternatives = sorted(
        test
        for g, family, test in runs
        if (
            g == gene
            and family == "branch"
            and test.endswith("_two_ratio")
        )
    )

    for alternative_test in branch_alternatives:
        hypothesis = alternative_test.removesuffix("_two_ratio")

        add_lrt(
            gene,
            hypothesis,
            "branch",
            (gene, "site", "M0"),
            (gene, "branch", alternative_test),
            "chi2_1",
            "branch",
        )

groups = sorted({row["fdr_group"] for row in results})

for group in groups:
    indices = [
        index
        for index, row in enumerate(results)
        if (
            row["fdr_group"] == group
            and row["status"] == "DONE"
            and row["p_value"] is not None
        )
    ]

    if not indices:
        continue

    p_values = [results[index]["p_value"] for index in indices]

    q_values = multipletests(
        p_values,
        method="fdr_bh",
    )[1]

    for index, q_value in zip(indices, q_values):
        results[index]["q_value"] = float(q_value)

result_fields = [
    "gene",
    "comparison",
    "test_family",
    "null_model",
    "alternative_model",
    "status",
    "lnL_null",
    "lnL_alt",
    "np_null",
    "np_alt",
    "raw_LRT",
    "LRT",
    "null_distribution",
    "p_value",
    "fdr_group",
    "q_value",
    "omega_background",
    "omega_foreground",
    "omega_alt_max",
    "null_outfile",
    "alternative_outfile",
]

with OUT_RESULTS.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=result_fields,
        delimiter="\t",
    )
    writer.writeheader()

    for row in results:
        writer.writerow({
            key: format_value(row.get(key))
            for key in result_fields
        })