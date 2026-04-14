#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import pandas as pd

REGION_RE = re.compile(r"(?P<prefix>.+\.region)(?P<region_num>\d+)$")


def convert_region_id(region_id: str) -> str:
    match = REGION_RE.match(region_id)
    if not match:
        raise ValueError(f"Cannot convert region id: {region_id}")

    prefix = match.group("prefix")
    region_num = match.group("region_num")
    return f"{prefix}{region_num}.gbk_region_{int(region_num)}"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert BGC region IDs like CP150300.1.region010 to CP150300.1.region010.gb_region_10"
    )
    parser.add_argument("input_csv", type=Path, help="Path to input CSV file")
    parser.add_argument(
        "output_csv",
        type=Path,
        help="Path to write output CSV file",
    )
    parser.add_argument(
        "--input-column",
        default="bgc_id",
        help="Name of the input ID column (default: bgc_id)",
    )
    parser.add_argument(
        "--output-column",
        default="bgc_id_converted",
        help="Name of the new output column (default: bgc_id_converted)",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite output file if it exists",
    )

    args = parser.parse_args()

    if args.output_csv.exists() and not args.overwrite:
        raise FileExistsError(
            f"Output file already exists: {args.output_csv}. Use --overwrite to replace."
        )

    df = pd.read_csv(args.input_csv)

    if args.input_column not in df.columns:
        raise KeyError(f"Input column not found: {args.input_column}")

    df[args.output_column] = df[args.input_column].astype(str).apply(convert_region_id)
    df.to_csv(args.output_csv, index=False)


if __name__ == "__main__":
    main()
