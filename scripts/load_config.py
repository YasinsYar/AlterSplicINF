#!/usr/bin/env python3
"""Print the branch YAML configuration as JSON for PowerShell runners."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    print(json.dumps(config, ensure_ascii=False))


if __name__ == "__main__":
    main()

