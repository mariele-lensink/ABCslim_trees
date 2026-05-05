#!/usr/bin/env python3
import csv
import gzip
import math
import os
import random
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
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
            intervals.append((start, stop, region_name))
    intervals.sort(key=lambda x: x[0])
    return intervals


def build_region_intervals(gene_csv: Path, intergene_csv: Path) -> List[Tuple[int, int, str]]:
    intervals = read_intervals(gene_csv, "genic") + read_intervals(intergene_csv, "intergenic")
    intervals.sort(key=lambda x: x[0])
    return intervals


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


def tajima_constants(n_accessions: int) -> Dict[str, float]:
    if n_accessions < 2:
        raise ValueError("Need at least two accessions to compute Tajima statistics")
    a1 = sum(1.0 / i for i in range(1, n_accessions))
    a2 = sum(1.0 / (i * i) for i in range(1, n_accessions))
    b1 = (n_accessions + 1.0) / (3.0 * (n_accessions - 1.0))
    b2 = 2.0 * (n_accessions * n_accessions + n_accessions + 3.0) / (9.0 * n_accessions * (n_accessions - 1.0))
    c1 = b1 - (1.0 / a1)
    c2 = b2 - ((n_accessions + 2.0) / (a1 * n_accessions)) + (a2 / (a1 * a1))
    e1 = c1 / a1
    e2 = c2 / (a1 * a1 + a2)
    return {
        "a1": a1,
        "e1": e1,
        "e2": e2,
    }


def zero_region_stats(length: int) -> Dict[str, float]:
    return {
        "sites_total": length,
        "snps": 0,
        "singletons": 0,
        "doubletons": 0,
        "pi": 0.0,
        "theta_w": 0.0,
        "td": 0.0,
        "pi_total": 0.0,
    }


def init_region_stats(lengths: Dict[str, int]) -> Dict[str, Dict[str, float]]:
    stats = {region: zero_region_stats(lengths[region]) for region in REGION_PREFIX}
    for region in stats:
        for label, _, _ in SFS_BINS:
            stats[region][f"sfs_{label}"] = 0
    return stats


def bin_label(derived_count: int) -> str:
    if derived_count == 1:
        return "1"
    if derived_count == 2:
        return "2"
    if derived_count <= 5:
        return "3_5"
    if derived_count <= 10:
        return "6_10"
    if derived_count <= 20:
        return "11_20"
    if derived_count <= 50:
        return "21_50"
    if derived_count <= 100:
        return "51_100"
    if derived_count <= 200:
        return "101_200"
    if derived_count <= 500:
        return "201_500"
    if derived_count <= 1000:
        return "501_1000"
    if derived_count >= 1001:
        return "1001plus"
    raise ValueError(f"Could not assign SFS bin for allele count {derived_count}")


def sfs_count(alt_count: int, accession_count: int) -> int:
    if _WORKER_CONFIG.get("fold_sfs", False):
        return min(alt_count, accession_count - alt_count)
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


def count_accessions(genotypes: List[str]) -> Tuple[int, int]:
    alt_count = 0
    accession_count = 0
    for field in genotypes:
        gt = field.split(":", 1)[0]
        if gt in {"0|0", "0/0"}:
            accession_count += 1
        elif gt in {"0|1", "1|0", "0/1", "1/0", "1|1", "1/1"}:
            accession_count += 1
            alt_count += 1
        else:
            digits = [char for char in gt if char in "01"]
            if not digits:
                continue
            accession_count += 1
            if "1" in digits:
                alt_count += 1
    return alt_count, accession_count


def finalize_region_stats(stats: Dict[str, float], tajima: Dict[str, float]) -> Dict[str, float]:
    snps = int(stats["snps"])
    if snps > 0:
        for label, _, _ in SFS_BINS:
            stats[f"sfs_{label}_prop"] = stats[f"sfs_{label}"] / snps
    else:
        for label, _, _ in SFS_BINS:
            stats[f"sfs_{label}_prop"] = 0.0

    length = int(stats["sites_total"])
    stats["pi"] = stats["pi_total"] / length if length > 0 else 0.0

    if snps > 0 and length > 0:
        theta_total = snps / tajima["a1"]
        stats["theta_w"] = theta_total / length
        denom = math.sqrt(tajima["e1"] * snps + tajima["e2"] * snps * (snps - 1))
        stats["td"] = (stats["pi_total"] - theta_total) / denom if denom > 0 else 0.0
    else:
        stats["theta_w"] = 0.0
        stats["td"] = 0.0

    return stats


def result_row(file_id: int, stats: Dict[str, Dict[str, float]]) -> Dict[str, float]:
    row: Dict[str, float] = {"ID": file_id}
    for region, prefix in REGION_PREFIX.items():
        region_stats = stats[region]
        row[f"{prefix}_sites_total"] = int(region_stats["sites_total"])
        row[f"{prefix}_snps"] = int(region_stats["snps"])
        row[f"{prefix}_singletons"] = int(region_stats["singletons"])
        row[f"{prefix}_doubletons"] = int(region_stats["doubletons"])
        row[f"{prefix}_pi"] = region_stats["pi"]
        row[f"{prefix}_theta_w"] = region_stats["theta_w"]
        row[f"{prefix}_td"] = region_stats["td"]
        for label, _, _ in SFS_BINS:
            row[f"{prefix}_sfs_{label}_prop"] = region_stats[f"sfs_{label}_prop"]
    return row


def process_vcf_file(file_id: int) -> Optional[Dict[str, float]]:
    config = _WORKER_CONFIG
    vcf_path = config["vcf_file"]
    if vcf_path is None:
        vcf_path = config["vcf_dir"] / f"{file_id}.vcf"
        if not vcf_path.exists():
            vcf_path = config["vcf_dir"] / f"{file_id}.vcf.gz"
        if not vcf_path.exists():
            return None

    locator = RegionLocator(config["intervals"])
    stats = init_region_stats(config["region_lengths"])
    expected_samples = None
    expected_accessions = None
    min_called_accessions = 1
    selected_sample_indices = None

    with open_text(vcf_path) as handle:
        for raw_line in handle:
            if raw_line.startswith("##"):
                continue
            if raw_line.startswith("#CHROM"):
                header_fields = raw_line.rstrip("\n").split("\t")
                expected_samples = len(header_fields) - 9
                selected_sample_indices = sample_indices(
                    expected_samples,
                    config["sample_individuals"],
                    config["sample_seed"],
                    file_id,
                )
                expected_accessions = (
                    len(selected_sample_indices) if selected_sample_indices is not None else expected_samples
                )
                min_called_accessions = max(
                    1,
                    int(math.ceil(expected_accessions * config["min_vcf_call_fraction"])),
                )
                continue

            fields = raw_line.rstrip("\n").split("\t")
            if config["vcf_chrom"] is not None and fields[0] != config["vcf_chrom"]:
                continue
            if config["snp_only"] and not is_biallelic_snp(fields[3], fields[4]):
                continue

            pos = int(fields[1])
            region = locator.classify(pos)
            if region is None:
                continue

            genotypes = fields[9:]
            if selected_sample_indices is not None:
                genotypes = [genotypes[idx] for idx in selected_sample_indices]

            alt_count, accession_count = count_accessions(genotypes)
            if accession_count == 0 or alt_count == 0 or alt_count == accession_count:
                continue

            if expected_accessions is not None:
                if config["require_complete_vcf_samples"]:
                    if accession_count != expected_accessions:
                        continue
                elif accession_count < min_called_accessions:
                    continue

            region_stats = stats[region]
            region_stats["snps"] += 1
            count_for_sfs = sfs_count(alt_count, accession_count)
            if count_for_sfs == 1:
                region_stats["singletons"] += 1
            if count_for_sfs == 2:
                region_stats["doubletons"] += 1

            pairwise_diff = (alt_count * (accession_count - alt_count)) / (
                accession_count * (accession_count - 1) / 2.0
            )
            region_stats["pi_total"] += pairwise_diff
            region_stats[f"sfs_{bin_label(count_for_sfs)}"] += 1

    if expected_accessions is None:
        return None

    tajima = tajima_constants(expected_accessions)
    for region in REGION_PREFIX:
        finalize_region_stats(stats[region], tajima)
    return result_row(file_id, stats)


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
    n_accessions = len(sampled_individual_nodes)
    if n_accessions < 2:
        return None

    sample_nodes = flatten_sampled_nodes(sampled_individual_nodes)
    stats = init_region_stats(config["region_lengths"])
    locator = RegionLocator(config["intervals"])

    for variant in ts.variants(samples=sample_nodes):
        pos = int(math.floor(variant.site.position)) + 1
        region = locator.classify(pos)
        if region is None:
            continue

        genotypes = variant.genotypes
        if genotypes.size != len(sample_nodes) or np.any(genotypes < 0):
            continue
        if int(genotypes.min()) != 0 or int(genotypes.max()) != 1:
            continue

        alt_count = 0
        offset = 0
        for node_ids in sampled_individual_nodes:
            width = len(node_ids)
            if np.any(genotypes[offset : offset + width] > 0):
                alt_count += 1
            offset += width

        if alt_count == 0 or alt_count == n_accessions:
            continue

        region_stats = stats[region]
        region_stats["snps"] += 1
        count_for_sfs = sfs_count(alt_count, n_accessions)
        if count_for_sfs == 1:
            region_stats["singletons"] += 1
        if count_for_sfs == 2:
            region_stats["doubletons"] += 1

        pairwise_diff = (alt_count * (n_accessions - alt_count)) / (
            n_accessions * (n_accessions - 1) / 2.0
        )
        region_stats["pi_total"] += pairwise_diff
        region_stats[f"sfs_{bin_label(count_for_sfs)}"] += 1

    tajima = tajima_constants(n_accessions)
    for region in REGION_PREFIX:
        finalize_region_stats(stats[region], tajima)
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
                f"{prefix}_snps",
                f"{prefix}_singletons",
                f"{prefix}_doubletons",
                f"{prefix}_pi",
                f"{prefix}_theta_w",
                f"{prefix}_td",
            ]
        )
        for label, _, _ in SFS_BINS:
            columns.append(f"{prefix}_sfs_{label}_prop")
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
    derived_dir = env_path("DERIVED_DIR", str(base_out / "constant" / chrom_dir / "trees_derived"))
    start_id = env_int("START_ID", 1)
    end_id = env_int("END_ID", 50000)
    output_file = env_path(
        "SIM_STATS_OUT",
        f"sim_stats_{chrom_dir}_{start_id}_{end_id}.csv",
    )
    n_cores = max(1, env_int("SIM_STATS_CORES", max(1, cpu_count() - 1)))
    sample_individuals = env_optional_int("SAMPLE_INDIVIDUALS")
    sample_seed = env_int("SAMPLE_SEED", 1729)
    fold_sfs = env_bool("FOLD_SFS", False)
    snp_only = env_bool("SNP_ONLY", False)
    require_complete_vcf_samples = env_bool("REQUIRE_COMPLETE_VCF_SAMPLES", True)
    min_vcf_call_fraction = env_float("MIN_VCF_CALL_FRACTION", 0.8)

    if start_id < 1 or end_id < start_id:
        raise ValueError("Invalid START_ID/END_ID values")
    if min_vcf_call_fraction <= 0 or min_vcf_call_fraction > 1:
        raise ValueError("MIN_VCF_CALL_FRACTION must be in the interval (0, 1]")

    intervals = build_region_intervals(gene_csv, intergene_csv)
    config = {
        "intervals": intervals,
        "region_lengths": region_lengths(intervals),
        "vcf_dir": vcf_dir,
        "vcf_file": vcf_file,
        "vcf_chrom": vcf_chrom,
        "derived_dir": derived_dir,
        "sample_individuals": sample_individuals,
        "sample_seed": sample_seed,
        "fold_sfs": fold_sfs,
        "snp_only": snp_only,
        "require_complete_vcf_samples": require_complete_vcf_samples,
        "min_vcf_call_fraction": min_vcf_call_fraction,
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
