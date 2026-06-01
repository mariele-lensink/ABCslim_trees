# calc_sim_stats_v2.py Methods Notes

Date: 2026-05-26

Script: `scripts/calc_sim_stats_v2.py`

## Purpose

`calc_sim_stats_v2.py` is an additive replacement candidate for the original
`scripts/calc_sim_stats.py`. It leaves the original script unchanged and writes a
separate v2 summary-statistic CSV. The main changes are:

- tree-sequence simulations are summarized with `tskit`/`pyslim` instead of
  iterating through every variant manually;
- VCF inputs use a hybrid `vcftools` + Python path for per-site statistics and
  interval aggregation;
- SFS, singleton, doubleton, and SNP-count columns are explicitly allele-copy
  based rather than accession-carrier based;
- genic and intergenic regions are handled as two separate concatenated interval
  sets with their own denominators and statistics.

## Inputs and Processing Order

The script keeps the same environment-variable style as the original script.
Important inputs are:

- `CHROM_DIR`: chromosome label, for example `chrom5`.
- `GENE_CSV` and `INTERGENE_CSV`: interval files. If unset, chromosome-specific
  defaults are used.
- `DERIVED_DIR`: directory containing derived `.trees` files.
- `VCF_DIR` or `VCF_FILE`: VCF source.
- `START_ID` and `END_ID`: simulation ID range.
- `SAMPLE_INDIVIDUALS`: number of diploid individuals to sample.
- `SAMPLE_SEED`: seed for reproducible individual/sample selection.
- `FOLD_SFS`: if true, SFS counts use the minor allele count.
- `SIM_STATS_OUT`: output CSV path.

For each ID, the script uses this order:

1. If `VCF_FILE` is set, run VCF mode for `VCF_FILE_ID`.
2. Otherwise, if `${DERIVED_DIR}/${ID}.trees` exists, run tree mode.
3. Otherwise, look for `${VCF_DIR}/${ID}.vcf` or `${VCF_DIR}/${ID}.vcf.gz`.

## Genic and Intergenic Sections

The script reads both interval files and labels each interval as either `genic`
or `intergenic`. Coordinates in the CSVs are treated as 1-based inclusive
intervals: `start` through `stop`, including both endpoints.

The regions are not analyzed as individual genes or individual intergenic
segments. Instead, all genic intervals are treated as one concatenated genic
region, and all intergenic intervals are treated as one concatenated intergenic
region.

For each region:

- `g_*` columns summarize all genic intervals together.
- `i_*` columns summarize all intergenic intervals together.
- `*_sites_total` is the sum of interval lengths:

```text
sum(stop - start + 1)
```

This means `g_sites_total` and `i_sites_total` are interval-denominator totals,
not counts of observed variant records.

## Windows and Interval Handling in Tree Mode

In tree mode, the script loads a `.trees` file with `tskit.load`. It imports
`pyslim` so SLiM metadata support is available, but the current implementation
does not need to inspect SLiM metadata directly.

The script first samples diploid individuals:

- sample nodes are grouped by individual;
- only individuals with exactly two sample nodes are treated as diploid;
- if `SAMPLE_INDIVIDUALS` is set, the script reproducibly samples that many
  diploid individuals using `SAMPLE_SEED + ID * 1_000_003`;
- the selected diploid individuals are flattened to haploid sample nodes.

For tskit statistics, each 1-based inclusive interval is converted to a
0-based half-open tree-sequence window:

```text
CSV interval:       [start, stop]
tskit interval:     [start - 1, stop)
```

`tskit` requires window arrays to be an increasing list from 0 to the full tree
sequence length. Therefore, for each target interval the script creates:

```text
[0, interval_left, interval_right, ts.sequence_length]
```

and then extracts only the statistic bin corresponding to
`[interval_left, interval_right)`. Statistics are summed across all intervals in
the region.

For complete simulations, all bases in the interval set are treated as callable:

```text
callable_sites = sites_total
```

## Denominators

The script tracks two denominator concepts:

- `*_sites_total`: total interval length from the genic/intergenic CSVs.
- `*_callable_sites`: callable denominator, currently equal to `sites_total` in
  tree mode. In VCF mode this is only populated from `CALLABLE_MASK`; otherwise
  it remains 0 because a variant-only VCF cannot distinguish missing invariant
  sites from true callable invariant sites.

For tree-sequence simulations:

- `pi`, `theta_w`, and `td` use `sites_total` as the denominator.
- This is appropriate because simulated tree sequences are treated as complete
  over the target intervals.

For VCF inputs:

- if `CALLABLE_MASK` is provided, `pi` and `theta_w` use `callable_sites`;
- if `CALLABLE_MASK` is not provided, `pi` and `theta_w` fall back to
  `sites_total`;
- `callable_records` and `low_call_records` report how many VCF records in the
  interval passed or failed call-rate filters, but these are record counts, not
  base-pair denominators.

Important consequence: for variant-only observed VCFs without a callable mask,
the script still cannot infer callable invariant sites. The denominator is then
the same interval length used for simulations, but this should be interpreted
carefully for observed data.

## Statistics

### `*_snps_copy`

This is the number of segregating sites represented in the allele-copy SFS after
excluding invariant and fixed sites. In tree mode, it is computed as the sum of
the binned allele-frequency-spectrum counts. In VCF mode, it is incremented for
biallelic segregating records that pass filters.

This is not the same as the old `*_snps` column in `calc_sim_stats.py`, which was
accession-carrier based.

### `*_singletons_copy`

This is the number of segregating sites whose allele-copy count is 1 after SFS
folding rules are applied.

With `FOLD_SFS=1`, a singleton is a site where the minor allele is present in one
haploid copy among the sampled haplotypes.

### `*_doubletons_copy`

This is the number of segregating sites whose allele-copy count is 2 after SFS
folding rules are applied.

With `FOLD_SFS=1`, a doubleton is a site where the minor allele is present in two
haploid copies among the sampled haplotypes.

### `*_sfs_copy_<bin>_prop`

These are folded or unfolded SFS bin proportions, depending on `FOLD_SFS`.

The bins are:

- `1`
- `2`
- `3_5`
- `6_10`
- `11_20`
- `21_50`
- `51_100`
- `101_200`
- `201_500`
- `501_1000`
- `1001plus`

Each bin is divided by `*_snps_copy`. For a region with at least one SNP, the SFS
proportions should sum to 1.0. If a region has no SNPs, all SFS proportions are
0.

These values describe the shape of the allele-frequency distribution. More mass
in low-count bins means more rare alleles; more mass in high-count bins means
more intermediate or common alleles.

### `*_pi`

Nucleotide diversity, pi, measures mean pairwise genetic differences per base.

In tree mode:

- `ts.diversity(..., mode="site", span_normalise=False)` is called for each
  interval;
- the unnormalized diversity totals are summed across intervals;
- the sum is divided by the region denominator.

In VCF mode:

- `vcftools --site-pi` provides per-site pi values;
- Python sums those values across passing records in each region;
- the sum is divided by `callable_sites` if a callable mask is supplied,
  otherwise by `sites_total`.

Interpretation: higher pi means more average pairwise diversity in that region.
Comparing `g_pi` and `i_pi` asks whether diversity differs between genic and
intergenic sequence.

### `*_theta_w`

Watterson's theta estimates population mutation rate from the number of
segregating sites.

The script computes:

```text
theta_total = snps_copy / a1
theta_w = theta_total / denominator
```

where:

```text
a1 = sum(1 / i) for i = 1 to n_haplotypes - 1
```

Interpretation: higher `theta_w` means more segregating sites per base, adjusted
for sample size. Differences between pi and theta_w contribute to Tajima's D.

### `*_td`

Tajima's D compares pairwise diversity to Watterson's theta:

```text
D = (pi_total - theta_total) / sqrt(e1*S + e2*S*(S - 1))
```

where `S` is `snps_copy`, and `e1` and `e2` are the standard Tajima constants
for the number of sampled haplotypes.

Interpretation:

- negative values mean an excess of rare alleles relative to neutral equilibrium
  expectations;
- positive values mean an excess of intermediate-frequency alleles;
- values near zero mean pi and theta_w are more concordant.

The script uses a single haploid sample size for a region. For tree simulations
this is straightforward. For VCF data with missing genotypes, this is more
approximate because per-site called haplotype counts can vary.

### `*_callable_records`

VCF-only diagnostic. This counts VCF records in the region that pass call-rate
filters. It is 0 in tree mode.

This is not a base-pair denominator. It is useful for checking how many VCF
records contributed to VCF statistics.

### `*_low_call_records`

VCF-only diagnostic. This counts VCF records in the region that failed call-rate
filters. It is 0 in tree mode.

The filter uses:

- `REQUIRE_COMPLETE_VCF_SAMPLES=1`: require all selected haplotypes to be called;
- otherwise require at least `MIN_VCF_CALL_FRACTION` of selected haplotypes.

## VCF Mode Details

VCF mode is intended for observed or VCF-formatted data. It uses `vcftools` to
write temporary per-site files:

- `--site-pi`
- `--missing-site`
- `--counts`

Python then reads those outputs, classifies each VCF position as genic or
intergenic, applies call-rate and SNP filters, and aggregates the statistics.

If `SNP_ONLY=1`, only biallelic A/C/G/T SNPs are retained for SNP-based
statistics. Multiallelic sites and indels are excluded by this filter.

Temporary vcftools output is written below `SIM_STATS_TMP_DIR`. It is removed
unless `KEEP_VCFTOOLS_TMP=1`.

## Output Columns

For each region prefix `g` and `i`, the script writes:

- `*_sites_total`
- `*_callable_sites`
- `*_snps_copy`
- `*_singletons_copy`
- `*_doubletons_copy`
- `*_pi`
- `*_theta_w`
- `*_td`
- `*_callable_records`
- `*_low_call_records`
- `*_sfs_copy_<bin>_prop` for all SFS bins

The output rows are sorted by numeric `ID`.

## Chromosome 5 Smoke Test

A tree-mode test was run on chromosome 5 IDs 1-3 using:

```bash
CHROM_DIR=chrom5
DERIVED_DIR=/group/gmonroegrp2/mlensink/ABC_data/chrom5/trees_derived
START_ID=1
END_ID=3
SAMPLE_INDIVIDUALS=1135
FOLD_SFS=1
SIM_STATS_CORES=1
SIM_STATS_OUT=/tmp/sim_stats_v2_chrom5_1_3.csv
```

The script wrote 3 rows and 43 columns. SFS proportions summed to 1.0 for both
genic and intergenic regions for all tested IDs.

Key output preview:

```text
ID,g_snps_copy,g_pi,i_snps_copy,i_pi
1,1780,0.00059039361,2260,0.0012853547
2,3110,0.00012310286,8370,0.00077475817
3,177,1.6357227e-05,8363,0.0013394097
```

`vcftools` was not available in the tested `abc_trees` environment, so VCF mode
was not smoke-tested there.
