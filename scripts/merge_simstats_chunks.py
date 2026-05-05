#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge chunked sim-stats CSV files.")
    parser.add_argument("chunk_dir", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--chrom-dir", default="chrom1")
    args = parser.parse_args()

    pattern = f"sim_stats_{args.chrom_dir}_*.csv"
    rows = []
    fieldnames = None
    for path in sorted(args.chunk_dir.glob(pattern)):
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            if fieldnames is None:
                fieldnames = reader.fieldnames
            elif reader.fieldnames != fieldnames:
                raise ValueError(f"Header mismatch in {path}")
            rows.extend(reader)

    if fieldnames is None:
        raise FileNotFoundError(f"No chunk files matched {args.chunk_dir / pattern}")

    rows.sort(key=lambda row: int(row["ID"]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {args.output}")


if __name__ == "__main__":
    main()
