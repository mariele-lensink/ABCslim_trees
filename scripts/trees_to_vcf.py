#!/usr/bin/env python3
import argparse
import os
import numpy as np
import pandas as pd
import pyslim
import msprime
import tskit

def read_intervals(csv_path):
    df = pd.read_csv(csv_path, header=None)
    starts = df.iloc[:, 0].astype(int).to_numpy()
    stops  = df.iloc[:, 1].astype(int).to_numpy()
    return list(zip(starts, stops))

def make_rate_map(L, gene_intervals, inter_intervals, gene_rate, inter_rate):
    # Build breakpoints from interval boundaries (half-open [start, stop+1))
    breaks = {0, L}
    for s, e in gene_intervals + inter_intervals:
        breaks.add(max(0, min(L, int(s))))
        breaks.add(max(0, min(L, int(e) + 1)))
    positions = np.array(sorted(breaks), dtype=float)

    # Helper: membership test using interval list (small; simple scan is fine)
    def in_any(intervals, x):
        for s, e in intervals:
            if s <= x <= e:
                return True
        return False

    rates = []
    for a, b in zip(positions[:-1], positions[1:]):
        mid = int((a + b) // 2)
        if in_any(gene_intervals, mid):
            rates.append(gene_rate)
        elif in_any(inter_intervals, mid):
            rates.append(inter_rate)
        else:
            rates.append(0.0)

    return msprime.RateMap(position=positions, rate=np.array(rates, dtype=float))

def keep_biallelic_only(ts):
    # Drop sites that are multi-allelic among samples (vcftools can be picky)
    to_drop = []
    for site in ts.sites():
        alleles = set()
        for mut in site.mutations:
            alleles.add(mut.derived_state)
        # include ancestral too
        alleles.add(site.ancestral_state)
        if len(alleles) > 2:
            to_drop.append(site.id)
    if len(to_drop) == 0:
        return ts
    return ts.delete_sites(to_drop)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trees", required=True)
    ap.add_argument("--out_vcf", required=True)
    ap.add_argument("--gene_csv", required=True)
    ap.add_argument("--intergene_csv", required=True)
    ap.add_argument("--L", type=int, required=True)

    ap.add_argument("--recapitate", action="store_true")
    ap.add_argument("--Ne", type=float, default=1135.0)
    ap.add_argument("--recomb", type=float, default=1.2e-5)

    ap.add_argument("--gmu", type=float, required=True)
    ap.add_argument("--imu", type=float, required=True)
    ap.add_argument("--gd", type=float, required=True)
    ap.add_argument("--id", dest="idp", type=float, required=True)

    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--biallelic_only", action="store_true")
    args = ap.parse_args()

    #neutral mutation rates are independent of eachother between genic and intergenic regions
    gene_intervals = read_intervals(args.gene_csv)
    inter_intervals = read_intervals(args.intergene_csv)

    ts = pyslim.load(args.trees)
#fill in the missing deep ancestry by recapitation
#For summaries like SFS and Tajima’s D, deep ancestry affects baseline diversity and the allele frequency spectrum.
    if args.recapitate:
        ts = pyslim.recapitate(
            ts,
            recombination_rate=args.recomb,
            Ne=args.Ne,
            random_seed=args.seed
        )

    # Overlay neutral mutations (keep existing deleterious ones)
    gmu_neu = args.gmu * (1.0 - args.gd)
    imu_neu = args.imu * (1.0 - args.idp)
    rm = make_rate_map(args.L, gene_intervals, inter_intervals, gmu_neu, imu_neu)
    ts = msprime.mutate(ts, rate=rm, keep=True, random_seed=args.seed)

    if args.biallelic_only:
        ts = keep_biallelic_only(ts)

    os.makedirs(os.path.dirname(args.out_vcf) or ".", exist_ok=True)
    with open(args.out_vcf, "w") as f:
        ts.write_vcf(f, ploidy=2)

    print(f"✅ wrote {args.out_vcf}")

if __name__ == "__main__":
    main()
