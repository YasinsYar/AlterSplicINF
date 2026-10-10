#!/usr/bin/env python3
"""Aggregate the six-sample methods branch and generate reports from the current run."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from junction_strand_check import inspect_junctions, junction_policy_signature


MARKERS = {
    "IFIT1": "ENSG00000185745",
    "ISG15": "ENSG00000187608",
    "MX1": "ENSG00000157601",
    "OAS1": "ENSG00000089127",
    "RSAD2": "ENSG00000134321",
    "STAT1": "ENSG00000115415",
}

FRESH_RESOURCE_COLUMNS = [
    "stage",
    "started_at",
    "finished_at",
    "elapsed_seconds",
    "free_before_gb",
    "minimum_free_gb",
    "free_after_gb",
    "transient_disk_delta_gb",
    "output_bytes",
    "exit_code",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_fresh_resource_metrics(metrics_dir: Path) -> pd.DataFrame:
    """Read the successful production stages recorded in one repeat namespace."""
    current = metrics_dir / "fresh_end_to_end.tsv"
    resources = pd.read_csv(current, sep="\t", encoding="utf-8-sig")
    if not set(FRESH_RESOURCE_COLUMNS).issubset(resources.columns):
        raise RuntimeError(f"Malformed repeat metric table: {current}")
    resources = resources[FRESH_RESOURCE_COLUMNS].copy()
    numeric_columns = [
        "elapsed_seconds",
        "free_before_gb",
        "minimum_free_gb",
        "free_after_gb",
        "transient_disk_delta_gb",
        "output_bytes",
        "exit_code",
    ]
    for column in numeric_columns:
        resources[column] = pd.to_numeric(resources[column], errors="coerce")
    resources = resources.loc[resources["exit_code"].eq(0)].copy()
    resources.sort_values("started_at", inplace=True)
    resources.drop_duplicates(subset=["stage"], keep="last", inplace=True)
    resources.sort_values("started_at", inplace=True)
    resources["elapsed_seconds"] = resources["elapsed_seconds"].round(1)
    disk_columns = ["free_before_gb", "minimum_free_gb", "free_after_gb", "transient_disk_delta_gb"]
    resources[disk_columns] = resources[disk_columns].round(3)
    resources[["output_bytes", "exit_code"]] = resources[["output_bytes", "exit_code"]].astype("int64")
    resources.reset_index(drop=True, inplace=True)
    return resources


def parse_star(path: Path) -> dict[str, float]:
    values: dict[str, float] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "|" not in line:
            continue
        key, raw = (part.strip() for part in line.split("|", 1))
        raw = raw.rstrip("%").strip()
        try:
            values[key] = float(raw)
        except ValueError:
            continue
    return values


def parse_featurecounts(path: Path) -> tuple[int, int]:
    table = pd.read_csv(path, sep="\t")
    count_col = table.columns[-1]
    counts = dict(zip(table["Status"], table[count_col], strict=False))
    assigned = int(counts.get("Assigned", 0))
    return assigned, int(table[count_col].sum())


def parse_rseqc(path: Path) -> tuple[float, float, float]:
    fractions: list[float] = []
    failed = math.nan
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = re.search(r"Fraction of reads failed to determine:\s*([0-9.]+)", line)
        if match:
            failed = float(match.group(1))
        match = re.search(r"Fraction of reads explained by .*:\s*([0-9.]+)", line)
        if match:
            fractions.append(float(match.group(1)))
    first = fractions[0] if fractions else math.nan
    second = fractions[1] if len(fractions) > 1 else math.nan
    dominant = max(first, second) if fractions else math.nan
    return failed, second, dominant


def read_gene_counts(path: Path) -> dict[str, int]:
    table = pd.read_csv(path, sep="\t", comment="#")
    count_col = table.columns[-1]
    counts: dict[str, int] = {}
    for gene_id, value in zip(table["Geneid"], table[count_col], strict=False):
        normalized = re.sub(r"\.\d+(?=_PAR_Y$|$)", "", str(gene_id))
        if normalized in counts:
            raise RuntimeError(f"Duplicate normalized gene ID: {normalized}")
        counts[normalized] = int(value)
    return counts


def count_nonempty_lines(path: Path) -> int:
    with path.open(encoding="utf-8", errors="replace") as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("track"))


def human_duration(seconds: float) -> str:
    hours = seconds / 3600
    return f"{hours:.2f} h" if hours >= 1 else f"{seconds / 60:.1f} min"


def safe_run(args: list[str]) -> str:
    try:
        return subprocess.run(args, check=True, capture_output=True, text=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def resolve_data_path(data_root: Path, configured: str) -> Path:
    path = Path(configured)
    return path.resolve() if path.is_absolute() else (data_root / path).resolve()


def redact_data_root(text: str, data_root: Path) -> str:
    variants = {
        str(data_root),
        str(data_root).replace("\\", "/"),
        str(data_root).replace("/", "\\"),
    }
    variants.update(value.replace("\\", "\\\\") for value in tuple(variants))
    for value in sorted(variants, key=len, reverse=True):
        text = text.replace(value, "<DATA_ROOT>")
    text = re.sub(r'[A-Za-z]:[/\\]+[^"\r\n]*[/\\]+config\.yaml', "<CODE_ROOT>/config.yaml", text, flags=re.I)
    return text


def validate_portable_release(branch: Path, config: dict) -> dict[str, object]:
    required = [
        "config.yaml", "samples.tsv", "source_metadata.tsv", "multiqc_config.yaml",
        "requirements.lock", "scripts/load_config.py", "scripts/common.ps1",
        "scripts/build_report.py", "scripts/junction_strand_check.py",
        "scripts/hash_fastq_content.py", "scripts/gtf_to_rseqc_bed.py",
        "scripts/run_pipeline.ps1", "scripts/prepare_reference.ps1",
        "scripts/process_samples.ps1", "scripts/run_sample.ps1", "scripts/run_repeat.ps1",
    ]
    missing = [name for name in required if not (branch / name).is_file()]
    absolute_config_paths = [
        name for name, value in config["paths"].items()
        if name != "code_outputs" and Path(str(value)).is_absolute()
    ]
    mutable_images = [name for name, image in config["containers"].items() if "@sha256:" not in image]
    result = {
        "branch_root": str(branch),
        "required_files_present": not missing,
        "missing_files": missing,
        "configured_data_paths_relative": not absolute_config_paths,
        "absolute_config_paths": absolute_config_paths,
        "containers_pinned_by_digest": not mutable_images,
        "mutable_container_entries": mutable_images,
        "analysis_release": config["analysis_release"],
    }
    if missing or absolute_config_paths or mutable_images:
        raise RuntimeError(f"Release portability check failed: {json.dumps(result)}")
    return result


def save_figure(fig: plt.Figure, stem: Path) -> None:
    fig.savefig(stem.with_suffix(".png"), dpi=240, bbox_inches="tight", facecolor="white")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def workflow_figure(path: Path, config: dict) -> None:
    fig, ax = plt.subplots(figsize=(12, 3.2))
    ax.axis("off")
    labels = [
        "Frozen metadata\n6 libraries",
        "Final reference\nGENCODE v50",
        f"QC + STAR\n{config['library']['star_twopass_mode']} 2-pass",
        f"Counts + junctions\n-s {config['library']['featurecounts_s']}; PE RF / SE XS",
        "MultiQC +\nreporting",
        "Generated technical\noutputs",
    ]
    colors = ["#dbeafe", "#ede9fe", "#dcfce7", "#fef3c7", "#fee2e2", "#e0f2fe"]
    xs = np.linspace(0.08, 0.92, len(labels))
    for idx, (x, label, color) in enumerate(zip(xs, labels, colors, strict=True)):
        ax.text(x, 0.68, label, ha="center", va="center", fontsize=10,
                bbox={"boxstyle": "round,pad=0.7", "facecolor": color, "edgecolor": "#334155"})
        if idx < len(labels) - 1:
            ax.annotate("", xy=(xs[idx + 1] - 0.065, 0.68), xytext=(x + 0.065, 0.68),
                        arrowprops={"arrowstyle": "->", "color": "#475569", "lw": 1.6})
    ax.text(0.67, 0.20, "Fresh repeats +\nstrand benchmark", ha="center", va="center", fontsize=10,
            bbox={"boxstyle": "round,pad=0.7", "facecolor": "#f1f5f9", "edgecolor": "#334155"})
    ax.annotate("", xy=(0.67, 0.33), xytext=(xs[3], 0.54),
                arrowprops={"arrowstyle": "->", "color": "#475569", "lw": 1.6})
    ax.annotate("", xy=(xs[5], 0.54), xytext=(0.78, 0.20),
                arrowprops={"arrowstyle": "->", "color": "#475569", "lw": 1.6})
    ax.set_title("Resource-aware RNA-seq preprocessing workflow", fontsize=14, weight="bold", pad=12)
    save_figure(fig, path)


def alignment_figure(qc: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 4.8))
    x = np.arange(len(qc))
    width = 0.36
    ax.bar(x - width / 2, qc["unique_mapping_percent"], width, label="Unique mapping", color="#2563eb")
    ax.bar(x + width / 2, qc["featurecounts_assigned_percent_of_star_input"], width, label="Gene-assigned input reads/fragments", color="#0d9488")
    ax.set_xticks(x, qc["sample_id"], rotation=25, ha="right")
    ax.set_ylabel("Percent")
    ax.set_ylim(0, 100)
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, ncol=2)
    ax.set_title("Alignment and stranded gene-assignment performance")
    save_figure(fig, path)


def junction_figure(qc: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    x = np.arange(len(qc))
    panels = [
        (qc["regtools_junction_records"] / 1000, "RegTools junction records (thousands)", "#7c3aed"),
        (qc["rseqc_reverse_among_classifiable"] * 100, "Reverse among classifiable reads (%)", "#ea580c"),
        (qc["rseqc_failed_fraction"] * 100, "Unclassifiable reads (%)", "#64748b"),
    ]
    for index, (values, ylabel, color) in enumerate(panels):
        axes[index].bar(x, values, color=color)
        axes[index].set_xticks(x, qc["sample_id"], rotation=35, ha="right")
        axes[index].set_ylabel(ylabel)
        if index > 0:
            axes[index].set_ylim(0, 100)
        axes[index].grid(axis="y", alpha=0.25)
    fig.suptitle("Extracted junctions and library-orientation diagnostics")
    fig.tight_layout()
    save_figure(fig, path)


def resource_figure(
    resources: pd.DataFrame,
    path: Path,
    fresh_resources: pd.DataFrame | None = None,
    se_fresh_resources: pd.DataFrame | None = None,
) -> None:
    compute_stages = {"fastqc", "star_2pass", "samtools_index", "samtools_flagstat", "featurecounts_s2", "regtools_rf", "regtools_xs", "rseqc_infer_experiment"}
    compute = resources[resources["stage"].isin(compute_stages)]
    per_sample = compute.groupby("sample_id", sort=False).agg(
        wall_hours=("elapsed_seconds", lambda x: x.sum() / 3600),
        peak_memory_gb=("peak_memory_mb", lambda x: x.max() / 1024),
    ).reset_index()
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.5))
    specs = [("wall_hours", "Compute stage time (h)", "#0284c7"), ("peak_memory_gb", "Observed peak RAM (GiB)", "#16a34a")]
    x = np.arange(len(per_sample))
    for ax, (column, label, color) in zip(axes[:2], specs, strict=True):
        ax.bar(x, per_sample[column], color=color)
        ax.set_xticks(x, per_sample["sample_id"], rotation=45, ha="right", fontsize=8)
        ax.set_ylabel(label)
        ax.grid(axis="y", alpha=0.25)
    disk_labels: list[str] = []
    disk_values: list[float] = []
    if fresh_resources is not None and not fresh_resources.empty:
        disk_labels.append("5INA PE100")
        disk_values.append(max(0.0, float(fresh_resources.iloc[0]["free_before_gb"]) - float(fresh_resources["minimum_free_gb"].min())))
    if se_fresh_resources is not None and not se_fresh_resources.empty:
        disk_labels.append("TIGK SE75")
        disk_values.append(max(0.0, float(se_fresh_resources.iloc[0]["free_before_gb"]) - float(se_fresh_resources["minimum_free_gb"].min())))
    if disk_values:
        axes[2].bar(disk_labels, disk_values, color=["#ca8a04", "#9333ea"][:len(disk_values)])
        axes[2].set_ylabel("Observed transient disk use (GiB)")
        axes[2].set_title("Fresh-from-SRA repeats: transient disk use")
        axes[2].grid(axis="y", alpha=0.25)
    else:
        axes[2].axis("off")
        axes[2].text(0.5, 0.5, "Transient disk monitoring\npending fresh repeat", ha="center", va="center")
    fig.suptitle("Recorded resource use for the protocol and fresh repeats")
    fig.tight_layout()
    save_figure(fig, path)


def marker_figure(markers: pd.DataFrame, samples: pd.DataFrame, path: Path) -> None:
    control = markers.loc[markers["condition"] == "control", ["dataset", "gene", "cpm"]].rename(columns={"cpm": "control_cpm"})
    treated = markers.loc[markers["condition"].isin(["IFN-I", "IFN-III"])].merge(control, on=["dataset", "gene"], how="left")
    treated["contrast"] = treated["dataset"] + " " + treated["condition"] + "/control"
    treated["log2_ratio"] = np.log2((treated["cpm"] + 0.5) / (treated["control_cpm"] + 0.5))
    contrast_order = ["GSE190899 IFN-I/control", "GSE190899 IFN-III/control", "GSE184456 IFN-I/control", "GSE184456 IFN-III/control"]
    matrix = treated.pivot(index="gene", columns="contrast", values="log2_ratio").reindex(index=MARKERS, columns=contrast_order)
    fig, ax = plt.subplots(figsize=(9.5, 4.5))
    limit = max(1.0, float(np.nanmax(np.abs(matrix.values))))
    image = ax.imshow(matrix.values, aspect="auto", cmap="coolwarm", vmin=-limit, vmax=limit)
    ax.set_xticks(np.arange(matrix.shape[1]), matrix.columns, rotation=30, ha="right")
    ax.set_yticks(np.arange(matrix.shape[0]), matrix.index)
    ax.set_title("Within-system descriptive interferon-response ratios (n=1 per condition)")
    bar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.03)
    bar.set_label("log2 treated/control CPM ratio")
    fig.tight_layout()
    save_figure(fig, path)




def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--branch-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--data-root", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    branch = args.branch_root.resolve()
    config_path = args.config.resolve() if args.config else branch / "config.yaml"
    config = yaml.safe_load(config_path.read_text(encoding="utf-8-sig"))
    configured_data_root = Path(config["paths"]["data_root"])
    data_root = (
        args.data_root.resolve()
        if args.data_root
        else configured_data_root.resolve() if configured_data_root.is_absolute()
        else (branch / configured_data_root).resolve()
    )
    portability = validate_portable_release(branch, config)
    if args.dry_run:
        print(json.dumps({**portability, "data_root": str(data_root), "dry_run": True}, ensure_ascii=False))
        return
    outputs = branch / "outputs"
    tables_dir = outputs / "tables"
    figures_dir = outputs / "figures"
    manifest_dir = outputs / "manifests"
    multiqc_release_dir = outputs / "multiqc"
    for directory in (tables_dir, figures_dir, manifest_dir, multiqc_release_dir):
        directory.mkdir(parents=True, exist_ok=True)

    samples = pd.read_csv(branch / "samples.tsv", sep="\t")
    if len(samples) != 6 or samples["sample_id"].nunique() != 6:
        raise RuntimeError("The methods manifest must contain exactly six unique samples")

    result_root = resolve_data_path(data_root, config["paths"]["data_results"])
    qc_rows: list[dict[str, object]] = []
    marker_rows: list[dict[str, object]] = []
    for sample in samples.itertuples(index=False):
        star = parse_star(result_root / "star" / sample.sample_id / f"{sample.sample_id}.Log.final.out")
        assigned, fc_total = parse_featurecounts(result_root / "counts" / f"{sample.sample_id}.featureCounts.txt.summary")
        failed, reverse_like, dominant = parse_rseqc(result_root / "rseqc" / f"{sample.sample_id}.infer_experiment.txt")
        classifiable = 1 - failed if not math.isnan(failed) else math.nan
        reverse_among_classifiable = reverse_like / classifiable if classifiable and classifiable > 0 else math.nan
        input_reads = int(star["Number of input reads"])
        total_splices = int(star.get("Number of splices: Total", 0))
        annotated_splices = int(star.get("Number of splices: Annotated (sjdb)", 0))
        gene_counts = read_gene_counts(result_root / "counts" / f"{sample.sample_id}.featureCounts.txt")
        library_total = sum(gene_counts.values())
        junction_path = result_root / "junctions" / f"{sample.sample_id}.junc"
        junction_check = inspect_junctions(junction_path)
        policy = json.loads(junction_path.with_suffix(".junc.policy.json").read_text(encoding="utf-8-sig"))
        if policy.get("signature") != junction_policy_signature(config, sample.layout) or policy.get("junction_sha256") != junction_check["junction_sha256"]:
            raise RuntimeError(f"Junction policy/artifact mismatch for {sample.sample_id}")
        if sample.layout == "SINGLE":
            if (
                policy.get("signature") != junction_policy_signature(config, "SINGLE")
                or policy.get("strand") != config["library"]["regtools_strand_se"]
                or policy.get("junction_sha256") != junction_check["junction_sha256"]
            ):
                raise RuntimeError(f"Stale single-end junction policy for {sample.sample_id}")
            if not all(junction_check["records_by_strand"][strand] for strand in ("+", "-")):
                raise RuntimeError(f"Single-end strand check failed for {sample.sample_id}")
        qc_rows.append({
            "sample_id": sample.sample_id,
            "dataset": sample.dataset,
            "condition": sample.condition,
            "layout": sample.layout,
            "input_reads": input_reads,
            "unique_mapping_percent": star["Uniquely mapped reads %"],
            "featurecounts_assigned": assigned,
            "featurecounts_total_classified": fc_total,
            "featurecounts_summary_assigned_percent": 100 * assigned / fc_total if fc_total else math.nan,
            "featurecounts_assigned_percent_of_star_input": 100 * assigned / input_reads if input_reads else math.nan,
            "star_total_splices": total_splices,
            "star_sjdb_matching_splice_observations": annotated_splices,
            "star_splice_observations_matching_sjdb_percent": 100 * annotated_splices / total_splices if total_splices else math.nan,
            "regtools_junction_records": junction_check["junction_records"],
            "regtools_plus_records": junction_check["records_by_strand"]["+"],
            "regtools_minus_records": junction_check["records_by_strand"]["-"],
            "regtools_unknown_records": junction_check["records_by_strand"]["?"] + junction_check["records_by_strand"]["."],
            "rseqc_failed_fraction": failed,
            "rseqc_reverse_like_fraction": reverse_like,
            "rseqc_dominant_fraction": dominant,
            "rseqc_classifiable_fraction": classifiable,
            "rseqc_reverse_among_classifiable": reverse_among_classifiable,
            "rseqc_reverse_among_classifiable_percent": 100 * reverse_among_classifiable,
            "rseqc_unclassifiable_percent": 100 * failed,
        })
        for gene, gene_id in MARKERS.items():
            count = gene_counts.get(gene_id, 0)
            cpm = count / library_total * 1_000_000 if library_total else 0
            marker_rows.append({
                "sample_id": sample.sample_id,
                "dataset": sample.dataset,
                "condition": sample.condition,
                "gene": gene,
                "gene_id": gene_id,
                "count": count,
                "cpm": cpm,
            })

    qc = pd.DataFrame(qc_rows)
    markers = pd.DataFrame(marker_rows)
    metric_frames = [pd.read_csv(path, sep="\t", encoding="utf-8-sig") for path in sorted((result_root / "metrics" / "pipeline").glob("*.steps.tsv"))]
    resources = pd.concat(metric_frames, ignore_index=True)
    numeric_cols = ["elapsed_seconds", "cpu_limit", "memory_limit_gb", "peak_memory_mb", "peak_cpu_percent", "output_bytes", "exit_code"]
    for column in numeric_cols:
        resources[column] = pd.to_numeric(resources[column], errors="coerce").fillna(0)

    # Cache reuse and failed attempts must not inflate stage time.
    resources = resources.loc[resources["exit_code"].eq(0) & resources["state"].eq("executed")].copy()
    resources.sort_values("started_at", inplace=True)
    resources.drop_duplicates(["sample_id", "stage"], keep="last", inplace=True)
    for sample in samples.itertuples(index=False):
        stages = set(resources.loc[resources.sample_id.eq(sample.sample_id), "stage"])
        expected_junction_stage = "regtools_rf" if sample.layout == "PAIRED" else "regtools_xs"
        expected = {"fastqc", "star_2pass", "samtools_index", "samtools_flagstat", "featurecounts_s2", expected_junction_stage, "rseqc_infer_experiment"}
        if not expected.issubset(stages) or (sample.layout == "SINGLE" and "regtools_rf" in stages):
            raise RuntimeError(f"Resource measurements do not match the protocol for {sample.sample_id}")

    repro_root = data_root / "results" / config["reprocheck"]["result_namespace"]
    sample_id = config["reprocheck"]["sample_id"]
    repeat_id = f"{sample_id}_fresh_repeat"
    repro_required = [
        repro_root / "REPROCHECK_COMPLETE.json",
        repro_root / "star" / f"{repeat_id}.Log.final.out",
        repro_root / "counts" / f"{repeat_id}.featureCounts.txt",
        repro_root / "counts" / f"{repeat_id}.featureCounts.txt.summary",
        repro_root / "junctions" / f"{repeat_id}.junc",
        repro_root / "rseqc" / f"{repeat_id}.infer_experiment.txt",
        repro_root / "input_manifest.json",
        repro_root / "metrics" / "fresh_end_to_end.tsv",
    ]
    repro_available = all(path.exists() for path in repro_required)
    repro_summary: dict[str, object] | None = None
    fresh_resources: pd.DataFrame | None = None
    strand_sensitivity: pd.DataFrame | None = None
    if repro_available:
        primary_star = parse_star(result_root / "star" / sample_id / f"{sample_id}.Log.final.out")
        repeat_star = parse_star(repro_root / "star" / f"{repeat_id}.Log.final.out")
        primary_assigned, _ = parse_featurecounts(result_root / "counts" / f"{sample_id}.featureCounts.txt.summary")
        repeat_assigned, _ = parse_featurecounts(repro_root / "counts" / f"{repeat_id}.featureCounts.txt.summary")
        primary_counts = read_gene_counts(result_root / "counts" / f"{sample_id}.featureCounts.txt")
        repeat_counts = read_gene_counts(repro_root / "counts" / f"{repeat_id}.featureCounts.txt")
        gene_ids = sorted(set(primary_counts) | set(repeat_counts))
        primary_vector = np.array([primary_counts.get(gene, 0) for gene in gene_ids], dtype=float)
        repeat_vector = np.array([repeat_counts.get(gene, 0) for gene in gene_ids], dtype=float)
        count_correlation = float(np.corrcoef(primary_vector, repeat_vector)[0, 1])
        primary_failed, primary_reverse, _ = parse_rseqc(result_root / "rseqc" / f"{sample_id}.infer_experiment.txt")
        repeat_failed, repeat_reverse, _ = parse_rseqc(repro_root / "rseqc" / f"{repeat_id}.infer_experiment.txt")
        repeat_manifest = json.loads((repro_root / "input_manifest.json").read_text(encoding="utf-8-sig"))
        primary_manifest = json.loads((result_root / "input_checksums" / f"{sample_id}.json").read_text(encoding="utf-8-sig"))
        fastq_content_equal = all(
            primary_manifest[mate]["decompressed_content_sha256"] == repeat_manifest[mate]["decompressed_content_sha256"]
            for mate in ("fastq", "fastq_r2")
        )
        pe_complete = json.loads((repro_root / "REPROCHECK_COMPLETE.json").read_text(encoding="utf-8-sig"))
        if not pe_complete.get("ready") or pe_complete.get("junction_policy_signature") != junction_policy_signature(config, "PAIRED"):
            raise RuntimeError("Paired repeat policy differs from the current protocol")
        primary_input = int(primary_star["Number of input reads"])
        repeat_input = int(repeat_star["Number of input reads"])
        repro_summary = {
            "sample_id": sample_id,
            "repeat_type": "fresh-from-SRA computational repeatability",
            "decompressed_fastq_content_equal": fastq_content_equal,
            "input_reads_equal": primary_input == repeat_input,
            "unique_mapping_primary_percent": primary_star["Uniquely mapped reads %"],
            "unique_mapping_repeat_percent": repeat_star["Uniquely mapped reads %"],
            "unique_mapping_delta_pp": repeat_star["Uniquely mapped reads %"] - primary_star["Uniquely mapped reads %"],
            "assignment_primary_percent_of_star_input": 100 * primary_assigned / primary_input,
            "assignment_repeat_percent_of_star_input": 100 * repeat_assigned / repeat_input,
            "assignment_delta_pp": 100 * repeat_assigned / repeat_input - 100 * primary_assigned / primary_input,
            "gene_count_pearson_r": count_correlation,
            "gene_counts_exact": bool(np.array_equal(primary_vector, repeat_vector)),
            "junction_file_exact": sha256(result_root / "junctions" / f"{sample_id}.junc") == sha256(repro_root / "junctions" / f"{repeat_id}.junc"),
            "rseqc_reverse_among_classifiable_primary": primary_reverse / (1 - primary_failed),
            "rseqc_reverse_among_classifiable_repeat": repeat_reverse / (1 - repeat_failed),
        }
        pd.DataFrame([repro_summary]).to_csv(tables_dir / "computational_repeatability.tsv", sep="\t", index=False, float_format="%.8f")
        fresh_resources = read_fresh_resource_metrics(repro_root / "metrics")
        fresh_resources["minimum_free_gb"] = fresh_resources[["minimum_free_gb", "free_after_gb"]].min(axis=1).round(3)
        fresh_resources["transient_disk_delta_gb"] = (
            fresh_resources["free_before_gb"] - fresh_resources["minimum_free_gb"]
        ).clip(lower=0).round(3)
        fresh_resources.to_csv(tables_dir / "fresh_pe100_resources.tsv", sep="\t", index=False)
        (tables_dir / "fresh_end_to_end_resources.tsv").unlink(missing_ok=True)
        sensitivity_rows = []
        for mode in config["reprocheck"]["strand_benchmark_modes"]:
            summary_path = repro_root / "strand_benchmark" / f"{repeat_id}.s{mode}.featureCounts.txt.summary"
            assigned, classified = parse_featurecounts(summary_path)
            sensitivity_rows.append({
                "sample_id": sample_id,
                "featurecounts_s": mode,
                "assigned": assigned,
                "star_input_reads_or_fragments": repeat_input,
                "assigned_percent_of_star_input": 100 * assigned / repeat_input,
                "assigned_percent_of_featurecounts_summary": 100 * assigned / classified,
            })
        strand_sensitivity = pd.DataFrame(sensitivity_rows)
        strand_sensitivity.to_csv(tables_dir / "strand_sensitivity.tsv", sep="\t", index=False, float_format="%.6f")
    else:
        for stale in ("computational_repeatability.tsv", "fresh_pe100_resources.tsv", "fresh_end_to_end_resources.tsv", "strand_sensitivity.tsv", "independent_repeat.tsv"):
            (tables_dir / stale).unlink(missing_ok=True)

    se_config = config["reprocheck_se"]
    se_repro_root = data_root / "results" / se_config["result_namespace"]
    se_sample_id = se_config["sample_id"]
    se_repeat_id = f"{se_sample_id}_fresh_repeat"
    se_required = [
        se_repro_root / "REPROCHECK_COMPLETE.json",
        se_repro_root / "star" / f"{se_repeat_id}.Log.final.out",
        se_repro_root / "counts" / f"{se_repeat_id}.featureCounts.txt",
        se_repro_root / "counts" / f"{se_repeat_id}.featureCounts.txt.summary",
        se_repro_root / "junctions" / f"{se_repeat_id}.junc",
        se_repro_root / "rseqc" / f"{se_repeat_id}.infer_experiment.txt",
        se_repro_root / "input_manifest.json",
        se_repro_root / "metrics" / "fresh_end_to_end.tsv",
    ]
    se_repro_available = all(path.exists() for path in se_required)
    se_repro_summary: dict[str, object] | None = None
    se_fresh_resources: pd.DataFrame | None = None
    if se_repro_available:
        se_complete = json.loads((se_repro_root / "REPROCHECK_COMPLETE.json").read_text(encoding="utf-8-sig"))
        if (
            se_complete.get("regtools_strand") != config["library"]["regtools_strand_se"]
            or se_complete.get("junction_policy_signature") != junction_policy_signature(config, "SINGLE")
        ):
            raise RuntimeError("The SE repeat uses a different junction policy")
        se_primary_star = parse_star(result_root / "star" / se_sample_id / f"{se_sample_id}.Log.final.out")
        se_repeat_star = parse_star(se_repro_root / "star" / f"{se_repeat_id}.Log.final.out")
        se_primary_assigned, _ = parse_featurecounts(result_root / "counts" / f"{se_sample_id}.featureCounts.txt.summary")
        se_repeat_assigned, _ = parse_featurecounts(se_repro_root / "counts" / f"{se_repeat_id}.featureCounts.txt.summary")
        se_primary_counts = read_gene_counts(result_root / "counts" / f"{se_sample_id}.featureCounts.txt")
        se_repeat_counts = read_gene_counts(se_repro_root / "counts" / f"{se_repeat_id}.featureCounts.txt")
        se_gene_ids = sorted(set(se_primary_counts) | set(se_repeat_counts))
        se_primary_vector = np.array([se_primary_counts.get(gene, 0) for gene in se_gene_ids], dtype=float)
        se_repeat_vector = np.array([se_repeat_counts.get(gene, 0) for gene in se_gene_ids], dtype=float)
        se_count_correlation = float(np.corrcoef(se_primary_vector, se_repeat_vector)[0, 1])
        se_primary_failed, se_primary_reverse, _ = parse_rseqc(result_root / "rseqc" / f"{se_sample_id}.infer_experiment.txt")
        se_repeat_failed, se_repeat_reverse, _ = parse_rseqc(se_repro_root / "rseqc" / f"{se_repeat_id}.infer_experiment.txt")
        se_primary_input = int(se_primary_star["Number of input reads"])
        se_repeat_input = int(se_repeat_star["Number of input reads"])
        se_input_manifest = json.loads((se_repro_root / "input_manifest.json").read_text(encoding="utf-8-sig"))
        primary_input_path = result_root / "input_checksums" / f"{se_sample_id}.json"
        se_primary_input_manifest = json.loads(primary_input_path.read_text(encoding="utf-8-sig"))
        se_input_content_equal = (
            se_primary_input_manifest["srr"] == se_input_manifest["srr"]
            and se_primary_input_manifest["fastq"]["decompressed_content_sha256"]
            == se_input_manifest["fastq"]["decompressed_content_sha256"]
        )
        (manifest_dir / "primary_se75_input_manifest.json").write_text(
            json.dumps(se_primary_input_manifest, indent=2) + "\n", encoding="utf-8"
        )
        se_repro_summary = {
            "sample_id": se_sample_id,
            "repeat_type": "fresh-from-SRA computational repeatability",
            "layout": se_input_manifest["layout"],
            "srr": se_input_manifest["srr"],
            "fresh_sra_sha256": se_input_manifest["sra_sha256"],
            "fresh_fastq_decompressed_sha256": se_input_manifest["fastq"]["decompressed_content_sha256"],
            "primary_fastq_decompressed_sha256": se_primary_input_manifest["fastq"]["decompressed_content_sha256"],
            "decompressed_fastq_content_equal": se_input_content_equal,
            "input_reads_equal": se_primary_input == se_repeat_input,
            "unique_mapping_primary_percent": se_primary_star["Uniquely mapped reads %"],
            "unique_mapping_repeat_percent": se_repeat_star["Uniquely mapped reads %"],
            "unique_mapping_delta_pp": se_repeat_star["Uniquely mapped reads %"] - se_primary_star["Uniquely mapped reads %"],
            "assignment_primary_percent_of_star_input": 100 * se_primary_assigned / se_primary_input,
            "assignment_repeat_percent_of_star_input": 100 * se_repeat_assigned / se_repeat_input,
            "assignment_delta_pp": 100 * se_repeat_assigned / se_repeat_input - 100 * se_primary_assigned / se_primary_input,
            "gene_count_pearson_r": se_count_correlation,
            "gene_counts_exact": bool(np.array_equal(se_primary_vector, se_repeat_vector)),
            "junction_file_exact": sha256(result_root / "junctions" / f"{se_sample_id}.junc") == sha256(se_repro_root / "junctions" / f"{se_repeat_id}.junc"),
            "rseqc_reverse_among_classifiable_primary": se_primary_reverse / (1 - se_primary_failed),
            "rseqc_reverse_among_classifiable_repeat": se_repeat_reverse / (1 - se_repeat_failed),
            "rseqc_unclassifiable_primary": se_primary_failed,
            "rseqc_unclassifiable_repeat": se_repeat_failed,
        }
        pd.DataFrame([se_repro_summary]).to_csv(
            tables_dir / "computational_repeatability_se75.tsv", sep="\t", index=False, float_format="%.8f"
        )
        se_fresh_resources = read_fresh_resource_metrics(se_repro_root / "metrics")
        se_fresh_resources["minimum_free_gb"] = se_fresh_resources[["minimum_free_gb", "free_after_gb"]].min(axis=1).round(3)
        se_fresh_resources["transient_disk_delta_gb"] = (
            se_fresh_resources["free_before_gb"] - se_fresh_resources["minimum_free_gb"]
        ).clip(lower=0).round(3)
        se_fresh_resources.to_csv(tables_dir / "fresh_se75_resources.tsv", sep="\t", index=False)
    else:
        for stale in ("computational_repeatability_se75.tsv", "fresh_se75_resources.tsv"):
            (tables_dir / stale).unlink(missing_ok=True)

    qc.to_csv(tables_dir / "technical_qc.tsv", sep="\t", index=False, float_format="%.4f")
    public_resources = resources.drop(columns=["peak_cpu_percent"], errors="ignore")
    public_resources.to_csv(tables_dir / "resource_metrics_by_step.tsv", sep="\t", index=False, float_format="%.3f")
    (tables_dir / "se_junction_refresh_metrics.tsv").unlink(missing_ok=True)
    markers.to_csv(tables_dir / "marker_sanity_check.tsv", sep="\t", index=False, float_format="%.4f")
    samples.to_csv(tables_dir / "sample_manifest.tsv", sep="\t", index=False)
    source_metadata = pd.read_csv(branch / "source_metadata.tsv", sep="\t")
    if set(source_metadata["sample_id"]) != set(samples["sample_id"]):
        raise RuntimeError("Source metadata must contain exactly the six frozen sample IDs")
    source_metadata.to_csv(tables_dir / "source_metadata.tsv", sep="\t", index=False)

    workflow_figure(figures_dir / "Figure1_workflow", config)
    alignment_figure(qc, figures_dir / "Figure2_alignment_assignment")
    junction_figure(qc, figures_dir / "Figure3_junctions_strand")
    resource_figure(resources, figures_dir / "Figure4_resources", fresh_resources, se_fresh_resources)
    marker_figure(markers, samples, figures_dir / "FigureS1_marker_sanity")
    # The descriptive marker panel is supplementary; remove its former figure name.
    for suffix in (".png", ".pdf"):
        (figures_dir / ("Figure5_marker_sanity" + suffix)).unlink(missing_ok=True)

    reference_manifest = json.loads((data_root / "reference" / "BRANCH_A_REFERENCE_READY.json").read_text(encoding="utf-8-sig"))
    multiqc_dir = result_root / "multiqc"
    multiqc_candidates = [multiqc_dir / "multiqc_report.html", *sorted(multiqc_dir.glob("*_multiqc_report.html"))]
    multiqc_report = next((path for path in multiqc_candidates if path.exists()), multiqc_dir / "multiqc_report.html")
    shutil.copy2(data_root / "reference" / "BRANCH_A_REFERENCE_READY.json", manifest_dir / "BRANCH_A_REFERENCE_READY.json")
    if multiqc_report.exists():
        shutil.copy2(multiqc_report, multiqc_release_dir / "multiqc_report.html")
    log_root = resolve_data_path(data_root, config["paths"]["data_logs"])
    for key, available, destination in (("reprocheck", repro_available, "pe100"), ("reprocheck_se", se_repro_available, "se75")):
        if available:
            source_root = data_root / "results" / config[key]["result_namespace"]
            for source, name in ((source_root / "input_manifest.json", f"repeat_{destination}_input.json"),
                                 (result_root / "input_checksums" / f"{config[key]['sample_id']}.json", f"primary_{destination}_input.json"),
                                 (log_root / config[key]["log_namespace"] / "commands.jsonl", f"repeat_{destination}_commands.jsonl")):
                content = source.read_text(encoding="utf-8-sig")
                role = "<PE_REPEAT>" if key == "reprocheck" else "<SE_REPEAT>"
                content = re.sub(r'(/data/(?:results|work)/)[^/\\\s"]+', lambda match: match[1] + role, content)
                if source.suffix == ".jsonl":
                    records = [json.loads(line) for line in content.splitlines() if line.strip()]
                    production_records = [record for record in records if record["stage"] not in ("featurecounts_s0", "featurecounts_s1")]
                    experiment_records = [record for record in records if record["stage"] in ("featurecounts_s0", "featurecounts_s1")]
                    content = "\n".join(json.dumps(record) for record in production_records) + "\n"
                    if experiment_records:
                        (manifest_dir / "strand_benchmark_commands.jsonl").write_text(redact_data_root("\n".join(json.dumps(record) for record in experiment_records) + "\n", data_root), encoding="utf-8")
                (manifest_dir / name).write_text(redact_data_root(content, data_root), encoding="utf-8")

    execution_records = []
    command_checks = []
    for sample in samples.itertuples(index=False):
        source = log_root / "pipeline" / "samples" / sample.sample_id / "commands.jsonl"
        records = [json.loads(line) for line in source.read_text(encoding="utf-8-sig").splitlines() if line.strip()]
        measurements = resources.loc[resources.sample_id.eq(sample.sample_id)]
        timestamps = dict(zip(measurements.stage, measurements.started_at))
        records = [record for record in records if record.get("time") == timestamps.get(record["stage"])]
        status = json.loads((source.parent / "status.json").read_text(encoding="utf-8-sig"))
        if status.get("phase") != "done":
            raise RuntimeError(f"Sample is incomplete: {sample.sample_id}")
        commands = {record["stage"]: record["command"] for record in records}
        junction_stage = "regtools_rf" if sample.layout == "PAIRED" else "regtools_xs"
        strand = config["library"]["regtools_strand"] if sample.layout == "PAIRED" else config["library"]["regtools_strand_se"]
        tokens = {
            "star_2pass": [f"/data/reference/{config['reference']['star_index']}", f"--twopassMode {config['library']['star_twopass_mode']}", f"--outSAMstrandField {config['library']['star_out_sam_strand_field']}"],
            "featurecounts_s2": [f"-a /data/reference/{config['reference']['gtf']}", "-s 2"],
            junction_stage: [f"-s {strand}", f"-a {config['library']['regtools_anchor']}", f"-m {config['library']['regtools_min_intron']}", f"-M {config['library']['regtools_max_intron']}"],
            "rseqc_infer_experiment": [f"/data/reference/{config['library']['rseqc_bed']}"],
        }
        command_checks.append(all(all(token in commands.get(stage, "") for token in required) for stage, required in tokens.items())
                              and not (sample.layout == "SINGLE" and "regtools_rf" in commands))
        execution_records.extend(records)
    executed_commands = "\n".join(json.dumps(record, ensure_ascii=False) for record in execution_records)
    (manifest_dir / "pipeline_commands.jsonl").write_text(redact_data_root(executed_commands, data_root) + "\n", encoding="utf-8")
    for stale in ("main_executed_commands.txt", "fresh_repeat_commands.jsonl", "fresh_se75_commands.jsonl",
                  "fresh_repeat_input_checksums.json", "fresh_se75_input_manifest.json", "primary_se75_input_manifest.json"):
        (manifest_dir / stale).unlink(missing_ok=True)
    forbidden = [result_root / name for name in ("leafcutter", "rmats", "deseq2", "enrichment", "differential_splicing")]
    checks = [
        ("frozen_sample_count", len(samples) == 6, f"{len(samples)} samples"),
        ("reference_ready", bool(reference_manifest.get("ready")), reference_manifest.get("star_index", "missing")),
        ("reference_hashes", reference_manifest.get("fasta_sha256") == config["reference"]["fasta_sha256"] and reference_manifest.get("gtf_sha256") == config["reference"]["gtf_sha256"], "FASTA/GTF hashes match the frozen config"),
        ("featurecounts_strand", config["library"]["featurecounts_s"] == 2, "-s 2"),
        ("regtools_strand", config["library"]["regtools_strand"] == "RF" and config["library"]["regtools_strand_se"] == "XS", "PE: RF; SE: XS from STAR intronMotif tags"),
        ("all_qc_rows", len(qc) == 6, f"{len(qc)} rows"),
        ("all_resource_rows", resources["sample_id"].nunique() == 6, f"{resources['sample_id'].nunique()} samples"),
        ("multiqc_report", multiqc_report.exists(), "multiqc/multiqc_report.html"),
        ("source_metadata", len(source_metadata) == 6, "six frozen GEO/SRA records"),
        ("immutable_container_images", all("@sha256:" in image for image in config["containers"].values()), "all images pinned by digest"),
        ("executed_commands_match_config", all(command_checks), "per-sample STAR/GTF/layout/strand/RSeQC commands match the protocol"),
        ("forbidden_outputs_absent", not any(path.exists() for path in forbidden), "LeafCutter/rMATS/DESeq2/enrichment absent"),
    ]
    checks.extend([("paired_repeat_complete", repro_available, "representative PE100"),
                   ("single_repeat_complete", se_repro_available, "representative SE75"),
                   ("strand_experiment_complete", strand_sensitivity is not None, "controlled featureCounts modes 0/1/2")])
    if repro_available and repro_summary is not None:
        repeat_passed = bool(
            repro_summary["decompressed_fastq_content_equal"]
            and repro_summary["input_reads_equal"]
            and abs(float(repro_summary["unique_mapping_delta_pp"])) < 0.01
            and abs(float(repro_summary["assignment_delta_pp"])) < 0.01
            and float(repro_summary["gene_count_pearson_r"]) > 0.99999
            and repro_summary["gene_counts_exact"]
            and repro_summary["junction_file_exact"]
            and abs(float(repro_summary["rseqc_reverse_among_classifiable_repeat"]) - float(repro_summary["rseqc_reverse_among_classifiable_primary"])) < 1e-8
        )
        checks.append(("fresh_computational_repeatability", repeat_passed, "5INA repeated from a fresh SRA acquisition within frozen tolerances"))
    if se_repro_available and se_repro_summary is not None:
        se_repeat_passed = bool(
            se_repro_summary["decompressed_fastq_content_equal"]
            and se_repro_summary["input_reads_equal"]
            and abs(float(se_repro_summary["unique_mapping_delta_pp"])) < 0.01
            and abs(float(se_repro_summary["assignment_delta_pp"])) < 0.01
            and float(se_repro_summary["gene_count_pearson_r"]) > 0.99999
            and se_repro_summary["gene_counts_exact"]
            and se_repro_summary["junction_file_exact"]
            and abs(float(se_repro_summary["rseqc_reverse_among_classifiable_repeat"]) - float(se_repro_summary["rseqc_reverse_among_classifiable_primary"])) < 1e-8
            and abs(float(se_repro_summary["rseqc_unclassifiable_repeat"]) - float(se_repro_summary["rseqc_unclassifiable_primary"])) < 1e-8
        )
        checks.append(("fresh_se75_computational_repeatability", se_repeat_passed, "TIGK_Unstim_1 repeated from a fresh SRA acquisition with exact counts, junctions, and RSeQC metrics"))
    check_frame = pd.DataFrame(checks, columns=["check", "passed", "detail"])
    check_frame.to_csv(tables_dir / "reproducibility_checks.tsv", sep="\t", index=False)
    if not check_frame["passed"].all():
        failed = check_frame.loc[~check_frame["passed"], "check"].tolist()
        raise RuntimeError(f"Reproducibility checks failed: {failed}")

    image_digests = {}
    for name, image in config["containers"].items():
        image_digests[name] = {
            "configured_reference": image,
            "repo_digests": safe_run(["docker", "image", "inspect", image, "--format", "{{json .RepoDigests}}"]),
        }
    environment = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "analysis_release": config["analysis_release"],
        "python_version": sys.version,
        "requirements_lock_sha256": sha256(branch / "requirements.lock"),
        "docker_server_version": safe_run(["docker", "version", "--format", "{{.Server.Version}}"]),
        "docker_cpus": safe_run(["docker", "info", "--format", "{{.NCPU}}"]),
        "docker_memory_bytes": safe_run(["docker", "info", "--format", "{{.MemTotal}}"]),
        "hardware": config["hardware"],
        "reference": reference_manifest,
        "container_images": image_digests,
    }
    (manifest_dir / "environment.json").write_text(json.dumps(environment, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    mapping_range = (qc["unique_mapping_percent"].min(), qc["unique_mapping_percent"].max())
    assignment_range = (qc["featurecounts_assigned_percent_of_star_input"].min(), qc["featurecounts_assigned_percent_of_star_input"].max())
    annotated_range = (qc["star_splice_observations_matching_sjdb_percent"].min(), qc["star_splice_observations_matching_sjdb_percent"].max())
    compute_stages = {"fastqc", "star_2pass", "samtools_index", "samtools_flagstat", "featurecounts_s2", "regtools_rf", "regtools_xs", "rseqc_infer_experiment"}
    total_wall = resources.loc[resources["stage"].isin(compute_stages), "elapsed_seconds"].sum()
    peak_ram_gb = resources["peak_memory_mb"].max() / 1024
    release = {
        "analysis_release": config["analysis_release"],
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "sample_count": len(samples),
        "all_checks_passed": True,
        "figures": 5,
        "main_figures": 4,
        "supplementary_figures": 1,
        "tables": len(list(tables_dir.glob("*.tsv"))),
        "fresh_computational_repeat_available": repro_available or se_repro_available,
        "fresh_pe100_repeat_available": repro_available,
        "fresh_se75_repeat_available": se_repro_available,
        "fresh_repeat_layouts": [layout for layout, available in (("PE100", repro_available), ("SE75", se_repro_available)) if available],
        "mapping_percent_range": [round(mapping_range[0], 4), round(mapping_range[1], 4)],
        "assignment_percent_range": [round(assignment_range[0], 4), round(assignment_range[1], 4)],
        "second_pass_sjdb_splice_observation_percent_range": [round(annotated_range[0], 4), round(annotated_range[1], 4)],
        "cumulative_compute_stage_seconds": round(total_wall, 1),
        "peak_observed_memory_mb": round(resources["peak_memory_mb"].max(), 1),
        "scope": "technical methods and computational repeatability only",
    }
    pe_benchmark_metrics = repro_root / "metrics" / "strand_benchmark.tsv"
    if not pe_benchmark_metrics.exists():
        raise RuntimeError("The full report requires resource measurements for the controlled strand experiment")
    if pe_benchmark_metrics.exists():
        benchmark_resources = pd.read_csv(pe_benchmark_metrics, sep="\t", encoding="utf-8-sig")
        benchmark_resources = benchmark_resources.loc[benchmark_resources.exit_code.eq(0)].drop_duplicates("stage", keep="last")
        benchmark_resources.to_csv(tables_dir / "strand_benchmark_resources.tsv", sep="\t", index=False)
    resource_summary = {
        "compute_hours": f"{total_wall / 3600:.2f}",
        "peak_gib": f"{peak_ram_gb:.2f}",
        "index_gib": f"{reference_manifest['index_size_gb']:.3f}",
        "index_hours": f"{reference_manifest['elapsed_seconds'] / 3600:.2f}",
        "pe_hours": f"{fresh_resources.elapsed_seconds.sum() / 3600:.2f}",
        "se_hours": f"{se_fresh_resources.elapsed_seconds.sum() / 3600:.2f}",
        "pe_disk_gib": f"{float(fresh_resources.iloc[0].free_before_gb) - fresh_resources.minimum_free_gb.min():.2f}",
        "se_disk_gib": f"{float(se_fresh_resources.iloc[0].free_before_gb) - se_fresh_resources.minimum_free_gb.min():.2f}",
        "pe_with_benchmark_hours": f"{(fresh_resources.elapsed_seconds.sum() + benchmark_resources.elapsed_seconds.sum()) / 3600:.2f}",
    }
    resource_summary_path = manifest_dir / "resource_summary.json"
    resource_summary_path.write_text(json.dumps(resource_summary, indent=2) + "\n", encoding="utf-8")
    release_path = outputs / "RELEASE.json"
    release_path.write_text(json.dumps(release, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (outputs / "BRANCH_A_RELEASE.json").unlink(missing_ok=True)
    print(json.dumps(release, ensure_ascii=False))


if __name__ == "__main__":
    main()
