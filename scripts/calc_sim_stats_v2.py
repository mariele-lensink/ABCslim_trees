#!/usr/bin/env python3
import csv
import gzip
import math
import os
import random
import shutil
import subprocess
import tempfile
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pyslim  # noqa: F401 - imported so SLiM metadata support is available
import tskit


SFS_BINS: List[Tuple[str, int, Optional[int]]] = [
    ("1", 1, 1),
    ("2", 2, 2),
    ("3_5", 3, 5),
    ("6_10", 6, 10),
    ("11_20", 11, 20),
    ("21_50", 21, 50),
    ("51_100", 51, 100),
    ("101_200", 101, 200),
    ("201_500", 201, 500),
    ("501_1000", 501, 1000),
    ("1001plus", 1001, None),
]

REGION_PREFIX = {
    "genic": "g",
    "intergenic": "i",
}

_WORKER_CONFIG = {}


def default_interval_paths(chrom_dir: str) -> Tuple[str, str]:
    mapping = {
        "chrom1": (
            "genomeinfo/gene100_chr1_10000000.csv",
            "genomeinfo/intergene100_chr1_10000000.csv",
        ),
        "chrom2": (
            "genomeinfo/gene100_chr2_18000000.csv",
            "genomeinfo/intergene100_chr2_18000000.csv",
        ),
        "chrom3": (
            "genomeinfo/gene100_chr3_22000000.csv",
            "genomeinfo/intergene100_chr3_22000000.csv",
        ),
        "chrom4": (
            "genomeinfo/gene100_chr4_2000000.csv",
            "genomeinfo/intergene100_chr4_2000000.csv",
        ),
        "chrom5": (
            "genomeinfo/gene100_chr5_1.csv",
            "genomeinfo/intergene100_chr5_1.csv",
        ),
    }
    return mapping.get(
        chrom_dir,
        ("genomeinfo/gene100_chr5_1.csv", "genomeinfo/intergene100_chr5_1.csv"),
    )


def env_int(name: str, default: int) -> int:
    return int(os.getenv(name, str(default)))


def env_float(name: str, default: float) -> float:
    return float(os.getenv(name, str(default)))


def env_optional_int(name: str) -> Optional[int]:
    value = os.getenv(name, "")
    if not value:
        return None
    parsed = int(value)
    if parsed < 1:
        raise ValueError(f"{name} must be >= 1 when set")
    return parsed


def env_nonnegative_int(name: str, default: int) -> int:
    parsed = int(os.getenv(name, str(default)))
    if parsed < 0:
        raise ValueError(f"{name} must be >= 0")
    return parsed


def env_bool(name: str, default: bool) -> bool:
    value = os.getenv(name, "")
    if not value:
        return default
    return value.strip().lower() in {"1", "true", "t", "yes", "y"}


def env_path(name: str, default: str) -> Path:
    return Path(os.getenv(name, default))


def env_optional_path(name: str) -> Optional[Path]:
    value = os.getenv(name, "")
    if not value:
        return None
    return Path(value)


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def read_intervals(path: Path, region_name: str) -> List[Tuple[int, int, str]]:
    intervals: List[Tuple[int, int, str]] = []
    with open_text(path) as handle:
        reader = csv.reader(handle)
        for row in reader:
            if not row:
                continue
            start = int(row[0])
            stop = int(row[1])
            if stop < start:
                raise ValueError(f"Invalid interval in {path}: {row}")
            intervals.append((start, stop, region_name))
    intervals.sort(key=lambda x: x[0])
    return intervals


def build_region_intervals(gene_csv: Path, intergene_csv: Path) -> List[Tuple[int, int, str]]:
    intervals = read_intervals(gene_csv, "genic") + read_intervals(intergene_csv, "intergenic")
    intervals.sort(key=lambda x: x[0])
    return intervals


def intervals_by_region(intervals: Iterable[Tuple[int, int, str]]) -> Dict[str, List[Tuple[int, int]]]:
    grouped = {"genic": [], "intergenic": []}
    for start, stop, region in intervals:
        grouped[region].append((start, stop))
    return grouped


def region_lengths(intervals: Iterable[Tuple[int, int, str]]) -> Dict[str, int]:
    lengths = {"genic": 0, "intergenic": 0}
    for start, stop, region in intervals:
        lengths[region] += stop - start + 1
    return lengths


class RegionLocator:
    def __init__(self, intervals: List[Tuple[int, int, str]]):
        self.intervals = intervals
        self.idx = 0

    def classify(self, pos: int) -> Optional[str]:
        while self.idx < len(self.intervals) and pos > self.intervals[self.idx][1]:
            self.idx += 1
        if self.idx >= len(self.intervals):
            return None
        start, stop, region = self.intervals[self.idx]
        if start <= pos <= stop:
            return region
        return None


def tajima_constants(n_haplotypes: int) -> Dict[str, float]:
    if n_haplotypes < 2:
        raise ValueError("Need at least two haplotypes to compute Tajima statistics")
    a1 = sum(1.0 / i for i in range(1, n_haplotypes))
    a2 = sum(1.0 / (i * i) for i in range(1, n_haplotypes))
    b1 = (n_haplotypes + 1.0) / (3.0 * (n_haplotypes - 1.0))
    b2 = 2.0 * (n_haplotypes * n_haplotypes + n_haplotypes + 3.0) / (
        9.0 * n_haplotypes * (n_haplotypes - 1.0)
    )
    c1 = b1 - (1.0 / a1)
    c2 = b2 - ((n_haplotypes + 2.0) / (a1 * n_haplotypes)) + (a2 / (a1 * a1))
    e1 = c1 / a1
    e2 = c2 / (a1 * a1 + a2)
    return {
        "a1": a1,
        "e1": e1,
        "e2": e2,
    }


def zero_region_stats(length: int) -> Dict[str, float]:
    stats = {
        "sites_total": length,
        "callable_sites": 0,
        "snps_copy": 0,
        "singletons_copy": 0,
        "doubletons_copy": 0,
        "pi": 0.0,
        "theta_w": 0.0,
        "td": 0.0,
        "pi_total": 0.0,
        "callable_records": 0,
        "low_call_records": 0,
        "no_alt_sites": 0,
        "variant_records": 0,
    }
    for label, _, _ in SFS_BINS:
        stats[f"sfs_copy_{label}"] = 0
        stats[f"sfs_copy_{label}_prop"] = 0.0
    return stats


def init_region_stats(lengths: Dict[str, int]) -> Dict[str, Dict[str, float]]:
    return {region: zero_region_stats(lengths[region]) for region in REGION_PREFIX}


def bin_label(allele_count: int) -> str:
    if allele_count == 1:
        return "1"
    if allele_count == 2:
        return "2"
    if allele_count <= 5:
        return "3_5"
    if allele_count <= 10:
        return "6_10"
    if allele_count <= 20:
        return "11_20"
    if allele_count <= 50:
        return "21_50"
    if allele_count <= 100:
        return "51_100"
    if allele_count <= 200:
        return "101_200"
    if allele_count <= 500:
        return "201_500"
    if allele_count <= 1000:
        return "501_1000"
    if allele_count >= 1001:
        return "1001plus"
    raise ValueError(f"Could not assign SFS bin for allele count {allele_count}")


def sfs_count(alt_count: int, haplotype_count: int, fold_sfs: bool) -> int:
    if fold_sfs:
        return min(alt_count, haplotype_count - alt_count)
    return alt_count


def is_biallelic_snp(ref: str, alt: str) -> bool:
    bases = {"A", "C", "G", "T"}
    return len(ref) == 1 and ref in bases and len(alt) == 1 and alt in bases


def sample_indices(total_count: int, sample_count: Optional[int], seed: int, file_id: int) -> Optional[List[int]]:
    if sample_count is None:
        return None
    if total_count < sample_count:
        raise ValueError(
            f"Cannot sample {sample_count} individuals from file ID {file_id}; only {total_count} are available"
        )
    rng = random.Random(seed + file_id * 1_000_003)
    return sorted(rng.sample(range(total_count), sample_count))


def sample_tree_individual_nodes(
    ts: tskit.TreeSequence, sample_count: Optional[int], seed: int, file_id: int
) -> List[List[int]]:
    nodes_by_individual: Dict[int, List[int]] = {}
    for node_id in ts.samples():
        individual_id = ts.node(node_id).individual
        if individual_id == tskit.NULL:
            continue
        nodes_by_individual.setdefault(individual_id, []).append(node_id)

    diploid_individuals = sorted(
        individual_id for individual_id, node_ids in nodes_by_individual.items() if len(node_ids) == 2
    )
    if sample_count is None:
        sampled_individuals = diploid_individuals
    else:
        if len(diploid_individuals) < sample_count:
            raise ValueError(
                f"Cannot sample {sample_count} diploid individuals from tree ID {file_id}; "
                f"only {len(diploid_individuals)} diploid individuals are available"
            )
        rng = random.Random(seed + file_id * 1_000_003)
        sampled_individuals = sorted(rng.sample(diploid_individuals, sample_count))

    return [sorted(nodes_by_individual[individual_id]) for individual_id in sampled_individuals]


def flatten_sampled_nodes(sampled_individual_nodes: List[List[int]]) -> np.ndarray:
    sampled_nodes: List[int] = []
    for node_ids in sampled_individual_nodes:
        sampled_nodes.extend(node_ids)
    return np.array(sampled_nodes, dtype=np.int32)


def count_allele_copies(genotypes: Sequence[str]) -> Tuple[int, int]:
    alt_count = 0
    haplotype_count = 0
    for field in genotypes:
        gt = field.split(":", 1)[0]
        if gt in {".", "./.", ".|."}:
            continue
        alleles = gt.replace("|", "/").split("/")
        called = [allele for allele in alleles if allele != "."]
        if not called:
            continue
        for allele in called:
            if allele not in {"0", "1"}:
                continue
            haplotype_count += 1
            if allele == "1":
                alt_count += 1
    return alt_count, haplotype_count


def finalize_region_stats(stats: Dict[str, float], n_haplotypes: int, denominator_key: str = "sites_total") -> None:
    snps = int(stats["snps_copy"])
    for label, _, _ in SFS_BINS:
        stats[f"sfs_copy_{label}_prop"] = stats[f"sfs_copy_{label}"] / snps if snps > 0 else 0.0

    denominator = int(stats[denominator_key])
    stats["pi"] = stats["pi_total"] / denominator if denominator > 0 else 0.0

    if snps > 0 and denominator > 0 and n_haplotypes >= 2:
        tajima = tajima_constants(n_haplotypes)
        theta_total = snps / tajima["a1"]
        stats["theta_w"] = theta_total / denominator
        denom = math.sqrt(tajima["e1"] * snps + tajima["e2"] * snps * (snps - 1))
        stats["td"] = (stats["pi_total"] - theta_total) / denom if denom > 0 else 0.0
    else:
        stats["theta_w"] = 0.0
        stats["td"] = 0.0


def result_row(file_id: int, stats: Dict[str, Dict[str, float]]) -> Dict[str, float]:
    row: Dict[str, float] = {"ID": file_id}
    for region, prefix in REGION_PREFIX.items():
        region_stats = stats[region]
        row[f"{prefix}_sites_total"] = int(region_stats["sites_total"])
        row[f"{prefix}_callable_sites"] = int(region_stats["callable_sites"])
        row[f"{prefix}_snps_copy"] = int(region_stats["snps_copy"])
        row[f"{prefix}_singletons_copy"] = int(region_stats["singletons_copy"])
        row[f"{prefix}_doubletons_copy"] = int(region_stats["doubletons_copy"])
        row[f"{prefix}_pi"] = region_stats["pi"]
        row[f"{prefix}_theta_w"] = region_stats["theta_w"]
        row[f"{prefix}_td"] = region_stats["td"]
        row[f"{prefix}_callable_records"] = int(region_stats["callable_records"])
        row[f"{prefix}_low_call_records"] = int(region_stats["low_call_records"])
        row[f"{prefix}_no_alt_sites"] = int(region_stats["no_alt_sites"])
        row[f"{prefix}_variant_records"] = int(region_stats["variant_records"])
        for label, _, _ in SFS_BINS:
            row[f"{prefix}_sfs_copy_{label}_prop"] = region_stats[f"sfs_copy_{label}_prop"]
    return row


def tskit_windows(intervals: Sequence[Tuple[int, int]], sequence_length: float) -> List[np.ndarray]:
    windows: List[np.ndarray] = []
    max_pos = float(sequence_length)
    for start, stop in intervals:
        left = max(0.0, float(start - 1))
        right = min(max_pos, float(stop))
        if right > left:
            windows.append(np.array([left, right], dtype=float))
    return windows


def as_window_vector(values) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.ndim == 0:
        return np.array([float(arr)])
    if arr.ndim == 1:
        return arr
    return np.squeeze(arr)


def sfs_array(values) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    arr = np.squeeze(arr)
    if arr.ndim == 1:
        return arr[np.newaxis, :]
    return arr


def region_tree_stats(
    ts: tskit.TreeSequence,
    sample_nodes: np.ndarray,
    windows: List[np.ndarray],
    n_haplotypes: int,
    length: int,
    fold_sfs: bool,
) -> Dict[str, float]:
    stats = zero_region_stats(length)
    stats["callable_sites"] = length
    if not windows:
        return stats

    pi_total = 0.0
    snps = 0
    sfs_counts = {label: 0 for label, _, _ in SFS_BINS}
    singletons = 0
    doubletons = 0

    for interval_window in windows:
        left = float(interval_window[0])
        right = float(interval_window[1])
        stat_windows = np.array(sorted({0.0, left, right, float(ts.sequence_length)}), dtype=float)
        target_idx = int(np.where(stat_windows == left)[0][0])
        diversity = ts.diversity(
            sample_sets=[sample_nodes],
            windows=stat_windows,
            mode="site",
            span_normalise=False,
        )
        _ = ts.segregating_sites(
            sample_sets=[sample_nodes],
            windows=stat_windows,
            mode="site",
            span_normalise=False,
        )
        afs = ts.allele_frequency_spectrum(
            sample_sets=[sample_nodes],
            windows=stat_windows,
            mode="site",
            polarised=not fold_sfs,
            span_normalise=False,
        )
        pi_total += float(as_window_vector(diversity)[target_idx])

        afs_rows = sfs_array(afs)[target_idx : target_idx + 1]
        for row in afs_rows:
            for allele_count, count in enumerate(row):
                count_int = int(round(float(count)))
                if count_int == 0 or allele_count == 0 or allele_count == n_haplotypes:
                    continue
                count_for_sfs = sfs_count(allele_count, n_haplotypes, fold_sfs)
                if count_for_sfs == 0:
                    continue
                label = bin_label(count_for_sfs)
                sfs_counts[label] += count_int
                if count_for_sfs == 1:
                    singletons += count_int
                if count_for_sfs == 2:
                    doubletons += count_int

    snps = sum(sfs_counts.values())
    stats["pi_total"] = pi_total
    stats["snps_copy"] = snps
    stats["singletons_copy"] = singletons
    stats["doubletons_copy"] = doubletons
    for label, count in sfs_counts.items():
        stats[f"sfs_copy_{label}"] = count

    finalize_region_stats(stats, n_haplotypes)
    return stats


def process_tree_file(file_id: int) -> Optional[Dict[str, float]]:
    config = _WORKER_CONFIG
    tree_path = config["derived_dir"] / f"{file_id}.trees"
    if not tree_path.exists():
        return None

    ts = tskit.load(tree_path)
    sampled_individual_nodes = sample_tree_individual_nodes(
        ts,
        config["sample_individuals"],
        config["sample_seed"],
        file_id,
    )
    if not sampled_individual_nodes:
        return None

    sample_nodes = flatten_sampled_nodes(sampled_individual_nodes)
    n_haplotypes = len(sample_nodes)
    if n_haplotypes < 2:
        return None

    stats: Dict[str, Dict[str, float]] = {}
    for region, intervals in config["intervals_by_region"].items():
        windows = tskit_windows(intervals, ts.sequence_length)
        stats[region] = region_tree_stats(
            ts,
            sample_nodes,
            windows,
            n_haplotypes,
            config["region_lengths"][region],
            config["fold_sfs"],
        )
    return result_row(file_id, stats)


def read_vcf_samples(vcf_path: Path) -> List[str]:
    with open_text(vcf_path) as handle:
        for raw_line in handle:
            if raw_line.startswith("#CHROM"):
                return raw_line.rstrip("\n").split("\t")[9:]
    raise ValueError(f"No #CHROM header found in {vcf_path}")


def write_keep_file(sample_names: Sequence[str], selected_indices: Optional[List[int]], path: Path) -> None:
    with path.open("w") as handle:
        if selected_indices is None:
            for sample in sample_names:
                handle.write(f"{sample}\n")
        else:
            for idx in selected_indices:
                handle.write(f"{sample_names[idx]}\n")


def vcf_arg(vcf_path: Path) -> str:
    return "--gzvcf" if vcf_path.suffix == ".gz" else "--vcf"


def run_vcftools(config: Dict[str, object], vcf_path: Path, keep_file: Path, prefix: Path) -> None:
    vcftools_bin = str(config["vcftools_bin"])
    if shutil.which(vcftools_bin) is None and not Path(vcftools_bin).exists():
        raise RuntimeError(
            f"VCF mode requires vcftools, but {vcftools_bin!r} was not found. "
            "Set VCFTOOLS_BIN or activate the abc_trees_vcf environment."
        )

    base_cmd = [
        vcftools_bin,
        vcf_arg(vcf_path),
        str(vcf_path),
        "--keep",
        str(keep_file),
    ]
    if config["vcf_chrom"] is not None:
        base_cmd.extend(["--chr", str(config["vcf_chrom"])])

    for option in ("--site-pi", "--missing-site", "--counts"):
        cmd = base_cmd + [option, "--out", str(prefix)]
        completed = subprocess.run(cmd, check=False, capture_output=True, text=True)
        if completed.returncode != 0:
            raise RuntimeError(
                f"vcftools failed for {option} with exit code {completed.returncode}\n"
                f"STDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
            )


def read_table_rows(path: Path) -> Iterable[Dict[str, str]]:
    if not path.exists():
        return
    with path.open() as handle:
        header_line = handle.readline()
        if not header_line:
            return
        header = header_line.strip().split()
        for raw_line in handle:
            parts = raw_line.strip().split()
            if len(parts) < len(header):
                continue
            yield dict(zip(header, parts))


def read_site_pi(path: Path) -> Dict[Tuple[str, int], float]:
    values: Dict[Tuple[str, int], float] = {}
    for row in read_table_rows(path) or []:
        chrom = row.get("CHROM") or row.get("CHR")
        pos = row.get("POS")
        pi = row.get("PI")
        if chrom is None or pos is None or pi is None:
            continue
        values[(chrom, int(pos))] = float(pi)
    return values


def read_site_missingness(path: Path) -> Dict[Tuple[str, int], float]:
    values: Dict[Tuple[str, int], float] = {}
    for row in read_table_rows(path) or []:
        chrom = row.get("CHR") or row.get("CHROM")
        pos = row.get("POS")
        f_miss = row.get("F_MISS")
        if chrom is None or pos is None or f_miss is None:
            continue
        values[(chrom, int(pos))] = float(f_miss)
    return values


def read_allele_counts(path: Path) -> Dict[Tuple[str, int], Tuple[int, Dict[str, int]]]:
    values: Dict[Tuple[str, int], Tuple[int, Dict[str, int]]] = {}
    if not path.exists():
        return values
    with path.open() as handle:
        header = handle.readline().strip().split()
        if not header:
            return values
        for raw_line in handle:
            row = raw_line.strip().split()
            if len(row) < 5:
                continue
            chrom = row[0]
            pos = int(row[1])
            n_chr = int(row[3])
            counts: Dict[str, int] = {}
            for entry in row[4:]:
                if ":" not in entry:
                    continue
                allele, count = entry.rsplit(":", 1)
                counts[allele] = int(count)
            values[(chrom, pos)] = (n_chr, counts)
    return values


def read_callable_mask(path: Path, vcf_chrom: Optional[str]) -> List[Tuple[int, int]]:
    intervals: List[Tuple[int, int]] = []
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.replace(",", "\t").split()
            if len(parts) >= 3:
                chrom, start_text, stop_text = parts[:3]
                if vcf_chrom is not None and chrom != vcf_chrom:
                    continue
                start = int(start_text) + 1
                stop = int(stop_text)
            elif len(parts) == 2:
                start = int(parts[0]) + 1
                stop = int(parts[1])
            else:
                continue
            if stop >= start:
                intervals.append((start, stop))
    intervals.sort(key=lambda x: x[0])
    return intervals


def overlap_length(a: Sequence[Tuple[int, int]], b: Sequence[Tuple[int, int]]) -> int:
    total = 0
    i = 0
    j = 0
    while i < len(a) and j < len(b):
        left = max(a[i][0], b[j][0])
        right = min(a[i][1], b[j][1])
        if right >= left:
            total += right - left + 1
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return total


def callable_lengths_by_region(
    mask_intervals: Sequence[Tuple[int, int]],
    grouped_intervals: Dict[str, List[Tuple[int, int]]],
) -> Dict[str, int]:
    return {
        region: overlap_length(grouped_intervals[region], mask_intervals)
        for region in REGION_PREFIX
    }


def process_vcf_file(file_id: int) -> Optional[Dict[str, float]]:
    config = _WORKER_CONFIG
    vcf_path = config["vcf_file"]
    if vcf_path is None:
        vcf_path = config["vcf_dir"] / f"{file_id}.vcf"
        if not vcf_path.exists():
            vcf_path = config["vcf_dir"] / f"{file_id}.vcf.gz"
        if not vcf_path.exists():
            return None

    stats = init_region_stats(config["region_lengths"])
    sample_names = read_vcf_samples(vcf_path)
    selected_sample_indices = sample_indices(
        len(sample_names),
        config["sample_individuals"],
        config["sample_seed"],
        file_id,
    )
    expected_individuals = (
        len(selected_sample_indices) if selected_sample_indices is not None else len(sample_names)
    )
    expected_haplotypes = expected_individuals * 2
    min_called_haplotypes = max(
        1,
        int(math.ceil(expected_haplotypes * config["min_vcf_call_fraction"])),
    )

    tmp_parent = Path(config["tmp_dir"])
    tmp_parent.mkdir(parents=True, exist_ok=True)
    tmp_context = None
    if config["keep_vcftools_tmp"]:
        tmp_dir = Path(tempfile.mkdtemp(prefix=f"vcftools_{file_id}_", dir=tmp_parent))
    else:
        tmp_context = tempfile.TemporaryDirectory(prefix=f"vcftools_{file_id}_", dir=tmp_parent)
        tmp_dir = Path(tmp_context.name)

    try:
        keep_file = tmp_dir / "samples.keep"
        write_keep_file(sample_names, selected_sample_indices, keep_file)
        prefix = tmp_dir / "vcftools"
        run_vcftools(config, vcf_path, keep_file, prefix)
        pi_by_site = read_site_pi(prefix.with_suffix(".sites.pi"))
        missing_by_site = read_site_missingness(prefix.with_suffix(".lmiss"))
        allele_counts_by_site = read_allele_counts(prefix.with_suffix(".frq.count"))

        locator = RegionLocator(config["intervals"])
        with open_text(vcf_path) as handle:
            for raw_line in handle:
                if raw_line.startswith("#"):
                    continue

                fields = raw_line.rstrip("\n").split("\t")
                chrom = fields[0]
                if config["vcf_chrom"] is not None and chrom != config["vcf_chrom"]:
                    continue
                pos = int(fields[1])
                local_pos = pos - config["vcf_pos_offset"]
                if local_pos < 1:
                    continue
                region = locator.classify(local_pos)
                if region is None:
                    continue

                site_key = (chrom, pos)
                vcftools_missing = missing_by_site.get(site_key)
                if vcftools_missing is not None:
                    called_haplotypes_estimate = int(round(expected_haplotypes * (1.0 - vcftools_missing)))
                    if config["require_complete_vcf_samples"]:
                        low_call = called_haplotypes_estimate != expected_haplotypes
                    else:
                        low_call = called_haplotypes_estimate < min_called_haplotypes
                    if low_call:
                        stats[region]["low_call_records"] += 1
                        continue

                counts_record = allele_counts_by_site.get(site_key)
                if counts_record is not None:
                    haplotype_count, allele_counts = counts_record
                    alt_count = allele_counts.get(fields[4], 0)
                else:
                    genotypes = fields[9:]
                    if selected_sample_indices is not None:
                        genotypes = [genotypes[idx] for idx in selected_sample_indices]
                    alt_count, haplotype_count = count_allele_copies(genotypes)

                if haplotype_count == 0:
                    stats[region]["low_call_records"] += 1
                    continue
                if config["require_complete_vcf_samples"]:
                    if haplotype_count != expected_haplotypes:
                        stats[region]["low_call_records"] += 1
                        continue
                elif haplotype_count < min_called_haplotypes:
                    stats[region]["low_call_records"] += 1
                    continue

                stats[region]["callable_records"] += 1
                if config["all_sites_vcf"]:
                    stats[region]["callable_sites"] += 1
                if fields[4] == ".":
                    stats[region]["no_alt_sites"] += 1
                    continue
                stats[region]["variant_records"] += 1

                if config["snp_only"] and not is_biallelic_snp(fields[3], fields[4]):
                    continue

                if alt_count == 0 or alt_count == haplotype_count:
                    continue

                count_for_sfs = sfs_count(alt_count, haplotype_count, config["fold_sfs"])
                if count_for_sfs == 0:
                    continue

                stats[region]["snps_copy"] += 1
                if count_for_sfs == 1:
                    stats[region]["singletons_copy"] += 1
                if count_for_sfs == 2:
                    stats[region]["doubletons_copy"] += 1

                site_pi = pi_by_site.get(site_key)
                if site_pi is None:
                    site_pi = (alt_count * (haplotype_count - alt_count)) / (
                        haplotype_count * (haplotype_count - 1) / 2.0
                    )
                stats[region]["pi_total"] += site_pi
                stats[region][f"sfs_copy_{bin_label(count_for_sfs)}"] += 1

    finally:
        if tmp_context is not None:
            tmp_context.cleanup()

    for region in REGION_PREFIX:
        if config["callable_lengths"] is not None:
            stats[region]["callable_sites"] = config["callable_lengths"][region]
            denominator_key = "callable_sites"
        elif config["all_sites_vcf"]:
            denominator_key = "callable_sites"
        else:
            denominator_key = "sites_total"
        finalize_region_stats(stats[region], expected_haplotypes, denominator_key=denominator_key)
    return result_row(file_id, stats)


def process_id(file_id: int) -> Optional[Dict[str, float]]:
    config = _WORKER_CONFIG
    if config["vcf_file"] is not None:
        return process_vcf_file(file_id)
    if (config["derived_dir"] / f"{file_id}.trees").exists():
        return process_tree_file(file_id)
    return process_vcf_file(file_id)


def worker_init(config: Dict[str, object]) -> None:
    global _WORKER_CONFIG
    _WORKER_CONFIG = config


def output_columns() -> List[str]:
    columns = ["ID"]
    for prefix in ("g", "i"):
        columns.extend(
            [
                f"{prefix}_sites_total",
                f"{prefix}_callable_sites",
                f"{prefix}_snps_copy",
                f"{prefix}_singletons_copy",
                f"{prefix}_doubletons_copy",
                f"{prefix}_pi",
                f"{prefix}_theta_w",
                f"{prefix}_td",
                f"{prefix}_callable_records",
                f"{prefix}_low_call_records",
                f"{prefix}_no_alt_sites",
                f"{prefix}_variant_records",
            ]
        )
        for label, _, _ in SFS_BINS:
            columns.append(f"{prefix}_sfs_copy_{label}_prop")
    return columns


def main() -> None:
    chrom_dir = os.getenv("CHROM_DIR", "chrom1")
    default_gene_csv, default_intergene_csv = default_interval_paths(chrom_dir)
    gene_csv = env_path("GENE_CSV", default_gene_csv)
    intergene_csv = env_path("INTERGENE_CSV", default_intergene_csv)
    base_out = env_path("BASE_OUT", "/home/mlensink/slimsimulations/ABCslim_trees/output")
    vcf_dir = env_path("VCF_DIR", str(base_out / "constant" / chrom_dir / "vcf"))
    vcf_file = env_optional_path("VCF_FILE")
    vcf_chrom = os.getenv("VCF_CHROM", "") or None
    vcf_pos_offset = env_nonnegative_int("VCF_POS_OFFSET", 0)
    derived_dir = env_path("DERIVED_DIR", str(base_out / "constant" / chrom_dir / "trees_derived"))
    start_id = env_int("START_ID", 1)
    end_id = env_int("END_ID", 50000)
    output_file = env_path(
        "SIM_STATS_OUT",
        f"sim_stats_v2_{chrom_dir}_{start_id}_{end_id}.csv",
    )
    n_cores = max(1, env_int("SIM_STATS_CORES", max(1, cpu_count() - 1)))
    sample_individuals = env_optional_int("SAMPLE_INDIVIDUALS")
    sample_seed = env_int("SAMPLE_SEED", 1729)
    fold_sfs = env_bool("FOLD_SFS", False)
    snp_only = env_bool("SNP_ONLY", False)
    require_complete_vcf_samples = env_bool("REQUIRE_COMPLETE_VCF_SAMPLES", True)
    min_vcf_call_fraction = env_float("MIN_VCF_CALL_FRACTION", 0.8)
    vcftools_bin = os.getenv("VCFTOOLS_BIN", "vcftools")
    tmp_dir = env_path("SIM_STATS_TMP_DIR", os.path.join(os.getenv("TMPDIR", "/tmp"), f"calc_sim_stats_v2_{os.getpid()}"))
    callable_mask = env_optional_path("CALLABLE_MASK")
    keep_vcftools_tmp = env_bool("KEEP_VCFTOOLS_TMP", False)
    all_sites_vcf = env_bool("ALL_SITES_VCF", False)

    if start_id < 1 or end_id < start_id:
        raise ValueError("Invalid START_ID/END_ID values")
    if min_vcf_call_fraction <= 0 or min_vcf_call_fraction > 1:
        raise ValueError("MIN_VCF_CALL_FRACTION must be in the interval (0, 1]")

    intervals = build_region_intervals(gene_csv, intergene_csv)
    grouped_intervals = intervals_by_region(intervals)
    callable_lengths = None
    if callable_mask is not None:
        mask_intervals = read_callable_mask(callable_mask, vcf_chrom)
        callable_lengths = callable_lengths_by_region(mask_intervals, grouped_intervals)
    config = {
        "intervals": intervals,
        "intervals_by_region": grouped_intervals,
        "region_lengths": region_lengths(intervals),
        "vcf_dir": vcf_dir,
        "vcf_file": vcf_file,
        "vcf_chrom": vcf_chrom,
        "vcf_pos_offset": vcf_pos_offset,
        "derived_dir": derived_dir,
        "sample_individuals": sample_individuals,
        "sample_seed": sample_seed,
        "fold_sfs": fold_sfs,
        "snp_only": snp_only,
        "require_complete_vcf_samples": require_complete_vcf_samples,
        "min_vcf_call_fraction": min_vcf_call_fraction,
        "vcftools_bin": vcftools_bin,
        "tmp_dir": tmp_dir,
        "keep_vcftools_tmp": keep_vcftools_tmp,
        "callable_lengths": callable_lengths,
        "all_sites_vcf": all_sites_vcf,
    }

    file_ids = [env_int("VCF_FILE_ID", 1)] if vcf_file is not None else list(range(start_id, end_id + 1))
    if n_cores == 1:
        worker_init(config)
        rows = [process_id(file_id) for file_id in file_ids]
    else:
        with Pool(processes=n_cores, initializer=worker_init, initargs=(config,)) as pool:
            rows = list(pool.imap(process_id, file_ids, chunksize=25))

    rows = [row for row in rows if row is not None]
    rows.sort(key=lambda row: int(row["ID"]))

    output_file.parent.mkdir(parents=True, exist_ok=True)
    columns = output_columns()
    with output_file.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"Wrote {len(rows)} rows to {output_file}")


if __name__ == "__main__":
    main()
