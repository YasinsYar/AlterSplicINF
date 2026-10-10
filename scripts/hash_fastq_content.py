#!/usr/bin/env python3
"""Hash compressed FASTQ files by compressed bytes and decompressed content."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path


def digest_stream(handle) -> str:
    digest = hashlib.sha256()
    for chunk in iter(lambda: handle.read(1024 * 1024), b""):
        digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+", type=Path)
    args = parser.parse_args()
    result = {}
    for path in args.files:
        with path.open("rb") as compressed:
            compressed_sha256 = digest_stream(compressed)
        with gzip.open(path, "rb") as decompressed:
            content_sha256 = digest_stream(decompressed)
        result[str(path)] = {
            "size_bytes": path.stat().st_size,
            "compressed_sha256": compressed_sha256,
            "decompressed_content_sha256": content_sha256,
        }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

