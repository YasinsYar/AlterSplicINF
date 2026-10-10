#!/usr/bin/env python3
"""Verify frozen release files without dependencies or scientific computation."""
import argparse
import hashlib
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.root.resolve()
    manifest = root / 'checksums.sha256'
    checked = 0
    for line in manifest.read_text(encoding='utf-8').splitlines():
        expected, relative = line.split('  ', 1)
        path = (root / relative).resolve()
        if not path.is_relative_to(root) or not path.is_file():
            parser.exit(1, f'Missing or unsafe release file: {relative}\n')
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            parser.exit(1, f'Checksum mismatch: {relative}\n')
        checked += 1
    print(f'SHA-256 verified: {checked} release files.')

if __name__ == '__main__':
    main()
