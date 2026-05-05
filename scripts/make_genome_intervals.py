#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def parse_args():
    ap = argparse.ArgumentParser(
        description=(
            "Build local-coordinate gene/intergene interval CSVs from a GFF by "
            "taking the first N gene features at or after a chromosome position."
        )
    )
    ap.add_argument("--gff", required=True)
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start-pos", type=int, required=True)
    ap.add_argument("--gene-count", type=int, default=100)
    ap.add_argument("--gene-out", required=True)
    ap.add_argument("--intergene-out", required=True)
    return ap.parse_args()


def iter_gene_features(gff_path, chrom, start_pos):
    with open(gff_path, newline="") as handle:
        for raw_line in handle:
            if not raw_line or raw_line.startswith("#"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid, _, feature_type, start, stop = fields[:5]
            if seqid != chrom or feature_type != "gene":
                continue
            start = int(start)
            stop = int(stop)
            if start >= start_pos:
                yield (start, stop)


def merge_intervals(intervals):
    merged = []
    for start, stop in intervals:
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, stop])
        else:
            merged[-1][1] = max(merged[-1][1], stop)
    return [(start, stop) for start, stop in merged]


def to_local(intervals, region_start):
    return [(start - region_start + 1, stop - region_start + 1) for start, stop in intervals]


def build_intergenic(merged_genes_local, region_length):
    intergenic = []
    if not merged_genes_local:
        return intergenic

    first_gene_start = merged_genes_local[0][0]
    if first_gene_start > 1:
        intergenic.append((1, first_gene_start - 1))

    for (_, prev_stop), (next_start, _) in zip(merged_genes_local, merged_genes_local[1:]):
        gap_start = prev_stop + 1
        gap_stop = next_start - 1
        if gap_start <= gap_stop:
            intergenic.append((gap_start, gap_stop))

    last_gene_stop = merged_genes_local[-1][1]
    if last_gene_stop < region_length:
        intergenic.append((last_gene_stop + 1, region_length))

    return intergenic


def write_intervals(path, intervals, label):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle)
        for start, stop in intervals:
            writer.writerow([start, stop, label])


def main():
    args = parse_args()

    raw_genes = []
    for start, stop in iter_gene_features(args.gff, args.chrom, args.start_pos):
        raw_genes.append((start, stop))
        if len(raw_genes) == args.gene_count + 1:
            break

    if len(raw_genes) < args.gene_count:
        raise SystemExit(
            f"Found only {len(raw_genes)} gene features on {args.chrom} at/after {args.start_pos}, "
            f"but {args.gene_count} were requested."
        )

    selected_genes = raw_genes[: args.gene_count]
    next_gene = raw_genes[args.gene_count] if len(raw_genes) > args.gene_count else None

    region_start = args.start_pos
    region_end = next_gene[0] - 1 if next_gene is not None else max(stop for _, stop in selected_genes)
    region_length = region_end - region_start + 1

    merged_genes = merge_intervals(selected_genes)
    merged_genes_local = to_local(merged_genes, region_start)
    intergenic_local = build_intergenic(merged_genes_local, region_length)

    gene_out = Path(args.gene_out)
    intergene_out = Path(args.intergene_out)
    gene_out.parent.mkdir(parents=True, exist_ok=True)
    intergene_out.parent.mkdir(parents=True, exist_ok=True)

    write_intervals(gene_out, merged_genes_local, "g2")
    write_intervals(intergene_out, intergenic_local, "g1")

    print(
        f"{args.chrom} start={args.start_pos} genes={len(selected_genes)} "
        f"merged_gene_segments={len(merged_genes_local)} intergenic_segments={len(intergenic_local)} "
        f"region_length={region_length}"
    )


if __name__ == "__main__":
    main()
