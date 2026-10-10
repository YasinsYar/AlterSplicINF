#!/usr/bin/env python3
"""Build a minimal transcript BED for RSeQC infer_experiment.py from a GTF."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--bed", type=Path, required=True)
    args = parser.parse_args()
    count = 0
    with args.gtf.open(encoding="utf-8") as source, args.bed.open("w", encoding="utf-8") as target:
        for line in source:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            transcript_id = "NA"
            for attribute in fields[8].split(";"):
                attribute = attribute.strip()
                if attribute.startswith("transcript_id"):
                    transcript_id = attribute.split(" ")[1].strip('"')
                    break
            target.write(f"{fields[0]}\t{int(fields[3]) - 1}\t{int(fields[4])}\t{transcript_id}\t0\t{fields[6]}\n")
            count += 1
    print(f"transcripts\t{count}")


if __name__ == "__main__":
    main()

