#!/usr/bin/env python3
import argparse
from pathlib import Path


def parse_args():
    ap = argparse.ArgumentParser(
        description="Subset a GFF to one chromosome region, preserving original coordinates."
    )
    ap.add_argument("--gff", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--out", required=True)
    return ap.parse_args()


def record_overlaps_region(fields, region_start, region_end):
    feature_type = fields[2]
    if feature_type == "chromosome":
        return True

    start = int(fields[3])
    end = int(fields[4])
    return start <= region_end and end >= region_start


def main():
    args = parse_args()
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(args.gff, newline="") as src, open(out_path, "w", newline="") as dst:
        dst.write(f"##subset chrom={args.chrom} start={args.start} end={args.end}\n")
        for raw_line in src:
            if not raw_line.strip():
                continue
            if raw_line.startswith("#"):
                dst.write(raw_line)
                continue

            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[0] != args.chrom:
                continue
            if record_overlaps_region(fields, args.start, args.end):
                dst.write(raw_line)

    print(f"{args.chrom}:{args.start}-{args.end} -> {out_path}")


if __name__ == "__main__":
    main()
