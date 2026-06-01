#!/usr/bin/env python3
"""Summarize missingness, callable denominators, and pi in all-sites VCF segments."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, median
from typing import Dict, Iterable, List, Sequence, TextIO, Tuple


DEFAULT_SEGMENT_SUMMARY = Path(
    "/group/gmonroegrp2/mlensink/ABC_data/real_data_segments/vcf_subsets_all_sites/"
    "segment_record_counts.tsv"
)
DEFAULT_OUT_DIR = Path("/group/gmonroegrp2/mlensink/1001data/invariantanalysis")
DEFAULT_THRESHOLDS = (0.50, 0.75, 0.80, 0.90, 0.95, 1.00)


@dataclass(frozen=True)
class Segment:
    name: str
    chrom: str
    start: int
    end: int
    full_length: int
    vcf_gz: Path


class ThresholdStats:
    def __init__(self) -> None:
        self.callable_sites = 0
        self.variant_records = 0
        self.biallelic_snp_records = 0
        self.polymorphic_biallelic_snps = 0
        self.pi_total = 0.0
        self.singletons = 0
        self.doubletons = 0
        self.low_frequency_variants = 0
        self.intermediate_frequency_variants = 0


class WindowStats:
    def __init__(self) -> None:
        self.n_sites = 0
        self.n_variant_records = 0
        self.callable_sites_q80 = 0
        self.missing_rate_sum = 0.0
        self.called_sum = 0
        self.pi_total_q80 = 0.0
        self.polymorphic_biallelic_snps_q80 = 0


def open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open("rt")


def open_maybe_gzip_write(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "wt", newline="")
    return path.open("w", newline="")


def parse_region(region: str) -> Tuple[str, int, int]:
    chrom, coords = region.split(":", 1)
    start_s, end_s = coords.split("-", 1)
    return chrom, int(start_s), int(end_s)


def load_segments(summary_tsv: Path) -> List[Segment]:
    segments: List[Segment] = []
    with summary_tsv.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            chrom, start, end = parse_region(row["region"])
            segments.append(
                Segment(
                    name=row["segment"],
                    chrom=chrom,
                    start=start,
                    end=end,
                    full_length=end - start + 1,
                    vcf_gz=Path(row["vcf_gz"]),
                )
            )
    return segments


def is_missing_gt(gt: str) -> bool:
    if gt in {"", ".", "./.", ".|."}:
        return True
    alleles = gt.replace("|", "/").split("/")
    return any(allele == "." or allele == "" for allele in alleles)


def gt_has_alt(gt: str) -> bool:
    alleles = gt.replace("|", "/").split("/")
    return any(allele not in {"0", "."} for allele in alleles)


def count_accession_genotypes(sample_fields: Sequence[str]) -> Tuple[int, int]:
    called = 0
    alt_carriers = 0
    for field in sample_fields:
        gt = field.split(":", 1)[0]
        if is_missing_gt(gt):
            continue
        called += 1
        if gt_has_alt(gt):
            alt_carriers += 1
    return called, alt_carriers


def is_biallelic_snp(ref: str, alt: str) -> bool:
    return len(ref) == 1 and len(alt) == 1 and alt not in {".", "*"} and "," not in alt


def quantile(values: Sequence[float], q: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    idx = (len(ordered) - 1) * q
    lo = math.floor(idx)
    hi = math.ceil(idx)
    if lo == hi:
        return ordered[lo]
    return ordered[lo] * (hi - idx) + ordered[hi] * (idx - lo)


def write_row(writer: csv.DictWriter, row: Dict[str, object]) -> None:
    writer.writerow(row)


def make_writer(path: Path, fieldnames: Sequence[str]) -> Tuple[TextIO, csv.DictWriter]:
    handle = open_maybe_gzip_write(path)
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    return handle, writer


def window_bounds(segment: Segment, pos: int, window_size: int) -> Tuple[int, int]:
    offset = pos - segment.start
    start = segment.start + (offset // window_size) * window_size
    end = min(start + window_size - 1, segment.end)
    return start, end


def summarize_segment(
    segment: Segment,
    out_dir: Path,
    thresholds: Sequence[float],
    window_size: int,
    per_site_writer: csv.DictWriter,
    per_sample_writer: csv.DictWriter,
    per_window_writer: csv.DictWriter,
    summary_writer: csv.DictWriter,
    denom_writer: csv.DictWriter,
    pi_writer: csv.DictWriter,
) -> None:
    site_missing_rates: List[float] = []
    site_called_counts: List[int] = []
    threshold_stats = {threshold: ThresholdStats() for threshold in thresholds}
    windows: Dict[Tuple[int, int], WindowStats] = {}

    n_samples = 0
    sample_ids: List[str] = []
    sample_missing: List[int] = []
    sample_called: List[int] = []

    n_sites = 0
    n_invariant_alt_dot = 0
    n_variant_alt_not_dot = 0
    n_sites_any_missing = 0
    n_biallelic_snp_records = 0
    n_polymorphic_biallelic_snps = 0

    with open_text(segment.vcf_gz) as handle:
        for raw_line in handle:
            if raw_line.startswith("##"):
                continue
            if raw_line.startswith("#CHROM"):
                fields = raw_line.rstrip("\n").split("\t")
                sample_ids = fields[9:]
                n_samples = len(sample_ids)
                sample_missing = [0] * n_samples
                sample_called = [0] * n_samples
                continue
            if not raw_line.strip():
                continue

            fields = raw_line.rstrip("\n").split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            samples = fields[9:]
            if n_samples == 0:
                raise ValueError(f"Missing #CHROM header before records in {segment.vcf_gz}")
            if len(samples) != n_samples:
                raise ValueError(f"Unexpected sample count at {segment.name}:{pos}")

            called, alt_carriers = count_accession_genotypes(samples)
            missing = n_samples - called
            call_fraction = called / n_samples if n_samples else 0.0
            missing_rate = missing / n_samples if n_samples else 0.0
            is_variant = alt != "."
            is_snp = is_biallelic_snp(ref, alt)
            is_poly_snp = is_snp and called >= 2 and 0 < alt_carriers < called

            n_sites += 1
            if is_variant:
                n_variant_alt_not_dot += 1
            else:
                n_invariant_alt_dot += 1
            if missing > 0:
                n_sites_any_missing += 1
            if is_snp:
                n_biallelic_snp_records += 1
            if is_poly_snp:
                n_polymorphic_biallelic_snps += 1

            site_missing_rates.append(missing_rate)
            site_called_counts.append(called)

            for idx, sample_field in enumerate(samples):
                gt = sample_field.split(":", 1)[0]
                if is_missing_gt(gt):
                    sample_missing[idx] += 1
                else:
                    sample_called[idx] += 1

            minor_allele_count = min(alt_carriers, called - alt_carriers) if called else 0
            allele_frequency = alt_carriers / called if called else 0.0
            pi_site = 0.0
            if is_poly_snp:
                pi_site = (alt_carriers * (called - alt_carriers)) / (called * (called - 1) / 2.0)

            write_row(
                per_site_writer,
                {
                    "subset_name": segment.name,
                    "chrom": chrom,
                    "position": pos,
                    "is_invariant_alt_dot": int(not is_variant),
                    "is_variant_alt_not_dot": int(is_variant),
                    "is_biallelic_snp": int(is_snp),
                    "is_polymorphic_biallelic_snp": int(is_poly_snp),
                    "n_called_individuals": called,
                    "n_missing_individuals": missing,
                    "call_fraction": f"{call_fraction:.8g}",
                    "missing_rate": f"{missing_rate:.8g}",
                    "alt_carrier_count": alt_carriers,
                    "ref_carrier_count": called - alt_carriers,
                    "minor_allele_count": minor_allele_count,
                    "allele_frequency_among_called": f"{allele_frequency:.8g}",
                    "pi_site": f"{pi_site:.12g}",
                },
            )

            win_start, win_end = window_bounds(segment, pos, window_size)
            win = windows.setdefault((win_start, win_end), WindowStats())
            win.n_sites += 1
            win.n_variant_records += int(is_variant)
            win.missing_rate_sum += missing_rate
            win.called_sum += called

            for threshold, stats in threshold_stats.items():
                if call_fraction < threshold:
                    continue
                stats.callable_sites += 1
                stats.variant_records += int(is_variant)
                stats.biallelic_snp_records += int(is_snp)
                if is_poly_snp:
                    stats.polymorphic_biallelic_snps += 1
                    stats.pi_total += pi_site
                    if minor_allele_count == 1:
                        stats.singletons += 1
                    if minor_allele_count == 2:
                        stats.doubletons += 1
                    if minor_allele_count <= 10:
                        stats.low_frequency_variants += 1
                    elif minor_allele_count <= max(10, called // 2):
                        stats.intermediate_frequency_variants += 1

                if math.isclose(threshold, 0.80):
                    win.callable_sites_q80 += 1
                    if is_poly_snp:
                        win.pi_total_q80 += pi_site
                        win.polymorphic_biallelic_snps_q80 += 1

    if n_sites == 0:
        raise ValueError(f"No VCF records found in {segment.vcf_gz}")

    for sample_id, missing, called in zip(sample_ids, sample_missing, sample_called):
        total = missing + called
        write_row(
            per_sample_writer,
            {
                "subset_name": segment.name,
                "individual_id": sample_id,
                "n_missing_genotypes": missing,
                "n_called_genotypes": called,
                "n_sites": total,
                "missing_rate": f"{(missing / total if total else 0.0):.8g}",
            },
        )

    for (win_start, win_end), win in sorted(windows.items()):
        write_row(
            per_window_writer,
            {
                "subset_name": segment.name,
                "chrom": segment.chrom,
                "window_start": win_start,
                "window_end": win_end,
                "window_length": win_end - win_start + 1,
                "n_sites": win.n_sites,
                "n_variant_records": win.n_variant_records,
                "callable_sites_q80": win.callable_sites_q80,
                "callable_fraction_q80": f"{win.callable_sites_q80 / win.n_sites:.8g}",
                "mean_site_missing_rate": f"{win.missing_rate_sum / win.n_sites:.8g}",
                "mean_called_individuals_per_site": f"{win.called_sum / win.n_sites:.8g}",
                "polymorphic_biallelic_snps_q80": win.polymorphic_biallelic_snps_q80,
                "pi_full_denominator_q80": f"{win.pi_total_q80 / (win_end - win_start + 1):.12g}",
                "pi_callable_denominator_q80": (
                    f"{win.pi_total_q80 / win.callable_sites_q80:.12g}"
                    if win.callable_sites_q80
                    else "0"
                ),
            },
        )

    mean_sample_missing = mean(
        missing / (missing + called) if (missing + called) else 0.0
        for missing, called in zip(sample_missing, sample_called)
    )
    max_sample_missing = max(
        missing / (missing + called) if (missing + called) else 0.0
        for missing, called in zip(sample_missing, sample_called)
    )

    write_row(
        summary_writer,
        {
            "subset_name": segment.name,
            "chrom": segment.chrom,
            "region_start": segment.start,
            "region_end": segment.end,
            "full_sequence_length": segment.full_length,
            "n_vcf_records": n_sites,
            "records_equal_full_length": int(n_sites == segment.full_length),
            "n_individuals": n_samples,
            "n_invariant_ALT_dot": n_invariant_alt_dot,
            "n_variant_ALT_not_dot": n_variant_alt_not_dot,
            "n_biallelic_snp_records": n_biallelic_snp_records,
            "n_polymorphic_biallelic_snps": n_polymorphic_biallelic_snps,
            "n_sites_with_any_missing_genotype": n_sites_any_missing,
            "mean_called_individuals_per_site": f"{mean(site_called_counts):.8g}",
            "median_called_individuals_per_site": f"{median(site_called_counts):.8g}",
            "fraction_sites_called_in_all_individuals": f"{sum(c == n_samples for c in site_called_counts) / n_sites:.8g}",
            "fraction_sites_called_in_at_least_90_percent_individuals": f"{sum(c / n_samples >= 0.90 for c in site_called_counts) / n_sites:.8g}",
            "mean_missing_rate_per_site": f"{mean(site_missing_rates):.8g}",
            "median_missing_rate_per_site": f"{median(site_missing_rates):.8g}",
            "p05_missing_rate_per_site": f"{quantile(site_missing_rates, 0.05):.8g}",
            "p95_missing_rate_per_site": f"{quantile(site_missing_rates, 0.95):.8g}",
            "mean_missing_rate_per_individual": f"{mean_sample_missing:.8g}",
            "max_missing_rate_per_individual": f"{max_sample_missing:.8g}",
        },
    )

    for threshold, stats in threshold_stats.items():
        callable_fraction = stats.callable_sites / segment.full_length if segment.full_length else 0.0
        write_row(
            denom_writer,
            {
                "subset_name": segment.name,
                "threshold": f"{threshold:.2f}",
                "full_sequence_length": segment.full_length,
                "callable_sites": stats.callable_sites,
                "callable_fraction": f"{callable_fraction:.8g}",
                "expected_full_denominator_pi_ratio": f"{callable_fraction:.8g}",
                "expected_downward_bias_fraction": f"{1.0 - callable_fraction:.8g}",
            },
        )
        write_row(
            pi_writer,
            {
                "subset_name": segment.name,
                "threshold": f"{threshold:.2f}",
                "n_sites_retained": stats.callable_sites,
                "n_sites_removed": segment.full_length - stats.callable_sites,
                "n_variant_records_retained": stats.variant_records,
                "fraction_variant_records_retained": (
                    f"{stats.variant_records / n_variant_alt_not_dot:.8g}"
                    if n_variant_alt_not_dot
                    else "0"
                ),
                "n_biallelic_snp_records_retained": stats.biallelic_snp_records,
                "n_polymorphic_biallelic_snps_retained": stats.polymorphic_biallelic_snps,
                "singleton_count": stats.singletons,
                "doubleton_count": stats.doubletons,
                "low_frequency_variant_count": stats.low_frequency_variants,
                "intermediate_frequency_variant_count": stats.intermediate_frequency_variants,
                "pairwise_difference_numerator": f"{stats.pi_total:.12g}",
                "full_sequence_length": segment.full_length,
                "callable_sites": stats.callable_sites,
                "pi_full_denominator": f"{stats.pi_total / segment.full_length:.12g}",
                "pi_callable_denominator": (
                    f"{stats.pi_total / stats.callable_sites:.12g}"
                    if stats.callable_sites
                    else "0"
                ),
                "relative_pi_full_vs_callable": (
                    f"{(stats.pi_total / segment.full_length) / (stats.pi_total / stats.callable_sites):.8g}"
                    if stats.callable_sites and stats.pi_total
                    else "0"
                ),
            },
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--segment-summary", type=Path, default=DEFAULT_SEGMENT_SUMMARY)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--window-size", type=int, default=10_000)
    parser.add_argument(
        "--thresholds",
        default=",".join(f"{x:.2f}" for x in DEFAULT_THRESHOLDS),
        help="Comma-separated call-fraction thresholds.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    thresholds = tuple(float(x) for x in args.thresholds.split(",") if x)
    if not thresholds:
        raise ValueError("At least one threshold is required")
    if any(threshold <= 0 or threshold > 1 for threshold in thresholds):
        raise ValueError("Thresholds must be in the interval (0, 1]")
    if args.window_size <= 0:
        raise ValueError("--window-size must be positive")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    segments = load_segments(args.segment_summary)

    summary_fields = [
        "subset_name",
        "chrom",
        "region_start",
        "region_end",
        "full_sequence_length",
        "n_vcf_records",
        "records_equal_full_length",
        "n_individuals",
        "n_invariant_ALT_dot",
        "n_variant_ALT_not_dot",
        "n_biallelic_snp_records",
        "n_polymorphic_biallelic_snps",
        "n_sites_with_any_missing_genotype",
        "mean_called_individuals_per_site",
        "median_called_individuals_per_site",
        "fraction_sites_called_in_all_individuals",
        "fraction_sites_called_in_at_least_90_percent_individuals",
        "mean_missing_rate_per_site",
        "median_missing_rate_per_site",
        "p05_missing_rate_per_site",
        "p95_missing_rate_per_site",
        "mean_missing_rate_per_individual",
        "max_missing_rate_per_individual",
    ]
    denom_fields = [
        "subset_name",
        "threshold",
        "full_sequence_length",
        "callable_sites",
        "callable_fraction",
        "expected_full_denominator_pi_ratio",
        "expected_downward_bias_fraction",
    ]
    pi_fields = [
        "subset_name",
        "threshold",
        "n_sites_retained",
        "n_sites_removed",
        "n_variant_records_retained",
        "fraction_variant_records_retained",
        "n_biallelic_snp_records_retained",
        "n_polymorphic_biallelic_snps_retained",
        "singleton_count",
        "doubleton_count",
        "low_frequency_variant_count",
        "intermediate_frequency_variant_count",
        "pairwise_difference_numerator",
        "full_sequence_length",
        "callable_sites",
        "pi_full_denominator",
        "pi_callable_denominator",
        "relative_pi_full_vs_callable",
    ]
    per_sample_fields = [
        "subset_name",
        "individual_id",
        "n_missing_genotypes",
        "n_called_genotypes",
        "n_sites",
        "missing_rate",
    ]
    per_site_fields = [
        "subset_name",
        "chrom",
        "position",
        "is_invariant_alt_dot",
        "is_variant_alt_not_dot",
        "is_biallelic_snp",
        "is_polymorphic_biallelic_snp",
        "n_called_individuals",
        "n_missing_individuals",
        "call_fraction",
        "missing_rate",
        "alt_carrier_count",
        "ref_carrier_count",
        "minor_allele_count",
        "allele_frequency_among_called",
        "pi_site",
    ]
    per_window_fields = [
        "subset_name",
        "chrom",
        "window_start",
        "window_end",
        "window_length",
        "n_sites",
        "n_variant_records",
        "callable_sites_q80",
        "callable_fraction_q80",
        "mean_site_missing_rate",
        "mean_called_individuals_per_site",
        "polymorphic_biallelic_snps_q80",
        "pi_full_denominator_q80",
        "pi_callable_denominator_q80",
    ]

    handles: List[TextIO] = []
    try:
        summary_handle, summary_writer = make_writer(args.out_dir / "segment_missingness_summary.tsv", summary_fields)
        denom_handle, denom_writer = make_writer(args.out_dir / "segment_callable_denominators.tsv", denom_fields)
        pi_handle, pi_writer = make_writer(args.out_dir / "segment_pi_by_call_threshold.tsv", pi_fields)
        sample_handle, sample_writer = make_writer(args.out_dir / "per_sample_missingness.tsv", per_sample_fields)
        site_handle, site_writer = make_writer(args.out_dir / "per_site_missingness.tsv.gz", per_site_fields)
        window_handle, window_writer = make_writer(args.out_dir / "per_window_missingness.tsv", per_window_fields)
        handles.extend([summary_handle, denom_handle, pi_handle, sample_handle, site_handle, window_handle])

        for segment in segments:
            print(f"Summarizing {segment.name}: {segment.vcf_gz}", flush=True)
            summarize_segment(
                segment=segment,
                out_dir=args.out_dir,
                thresholds=thresholds,
                window_size=args.window_size,
                per_site_writer=site_writer,
                per_sample_writer=sample_writer,
                per_window_writer=window_writer,
                summary_writer=summary_writer,
                denom_writer=denom_writer,
                pi_writer=pi_writer,
            )
    finally:
        for handle in handles:
            handle.close()

    manifest = args.out_dir / "analysis_manifest.txt"
    with manifest.open("w") as handle:
        handle.write(f"segment_summary={args.segment_summary}\n")
        handle.write(f"out_dir={args.out_dir}\n")
        handle.write(f"window_size={args.window_size}\n")
        handle.write(f"thresholds={','.join(f'{x:.2f}' for x in thresholds)}\n")
        handle.write("outputs=\n")
        for name in [
            "segment_missingness_summary.tsv",
            "segment_callable_denominators.tsv",
            "segment_pi_by_call_threshold.tsv",
            "per_sample_missingness.tsv",
            "per_site_missingness.tsv.gz",
            "per_window_missingness.tsv",
        ]:
            handle.write(f"  {args.out_dir / name}\n")

    print(f"Wrote invariant analysis outputs to {args.out_dir}", flush=True)


if __name__ == "__main__":
    main()
