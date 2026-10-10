#!/usr/bin/env python3
"""Check the strand field in BED12 junction artifacts before BAM cleanup."""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path


def junction_policy_signature(config: dict, layout: str) -> str:
    """Use the same ordered policy identity as the PowerShell sample runner."""
    library = config["library"]
    policy = {
        "layout": layout,
        "strand": library["regtools_strand_se"] if layout == "SINGLE" else library["regtools_strand"],
        "anchor": library["regtools_anchor"],
        "min_intron": library["regtools_min_intron"],
        "max_intron": library["regtools_max_intron"],
        "regtools_image": config["containers"]["regtools"],
        "star_image": config["containers"]["star"],
        "star_mode": library["star_twopass_mode"],
        "star_strand_field": library["star_out_sam_strand_field"],
        "fasta_sha256": config["reference"]["fasta_sha256"],
        "gtf_sha256": config["reference"]["gtf_sha256"],
    }
    return hashlib.sha256(json.dumps(policy, separators=(",", ":"), ensure_ascii=False).encode("utf-8")).hexdigest()


def inspect_junctions(path: Path) -> dict:
    rows = Counter()
    support = Counter()
    with path.open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("track"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 12 or fields[5] not in {"+", "-", "?", "."}:
                raise ValueError(f"Invalid BED12 junction at line {number}")
            if int(fields[9]) != 2 or int(fields[4]) < 1:
                raise ValueError(f"Invalid junction blocks/support at line {number}")
            rows[fields[5]] += 1
            support[fields[5]] += int(fields[4])
    if not rows:
        raise ValueError("Junction file is empty")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return {
        "junction_records": sum(rows.values()),
        "records_by_strand": {key: rows[key] for key in ("+", "-", "?", ".")},
        "alignment_support_by_strand": {key: support[key] for key in ("+", "-", "?", ".")},
        "junction_sha256": digest.hexdigest(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--junctions", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--require-both-strands", action="store_true")
    args = parser.parse_args()
    report = inspect_junctions(args.junctions)
    if args.require_both_strands and not all(report["records_by_strand"][key] for key in ("+", "-")):
        raise ValueError("The evaluated bulk SE libraries must contain junctions on both genomic strands")
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report))


if __name__ == "__main__":
    main()
