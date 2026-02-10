# ABCslim: Forward Simulations for Intragenomic Mutation Rate Heterogeneity

This repository contains a reproducible simulation and inference framework for studying **intragenomic mutation rate heterogeneity** using forward evolutionary simulations in **SLiM**, combined with **Approximate Bayesian Computation (ABC)** and **ABC Random Forests (ABC-RF)**.

The project is designed to disentangle the effects of **mutation rate variation**, **selection**, and **demography** on patterns of genetic variation—particularly differences between **genic** and **intergenic** regions—using summary statistics comparable to population-genomic data.

### Key Features
- Forward simulations in **SLiM** with **tree-sequence recording**
- Log-uniform parameter sampling over biologically realistic ranges
- Large-scale parallelization via **SLURM**
- Modular summary-statistic calculation
- ABC and ABC-RF–based parameter inference in R

---
## Initialize Environment
```conda env create -f env/abc_trees_env.yaml```

## Parameterization

Each simulation draws a unique parameter set from prior distributions describing mutation rate heterogeneity and selection differences between genic and intergenic regions.

### Parameters

| Parameter | Description |
|----------|-------------|
| `gmu` | Genic mutation rate |
| `imu` | Intergenic mutation rate |
| `gd` | Fraction of genic sites under selection |
| `id` | Fraction of intergenic sites under selection |
| `gdfe` | Mean deleterious effect size for genic sites |
| `idfe` | Mean deleterious effect size for intergenic sites |

### Prior Distributions
- Mutation rates (`gmu`, `imu`): **log-uniform**
- DFEs (`gdfe`, `idfe`): **log-uniform (negative)**
- Selection fractions (`gd`, `id`): **uniform**

Parameter tables are generated programmatically to ensure reproducibility and coverage of biologically realistic parameter space:

```Rscript scripts/generate_priors.R```

## Running simulations on the cluster (SLURM job arrays)

Simulations are run in parallel using a SLURM array job. Each array task reads **one row** from the priors CSV (one parameter set) and runs **three demography models** (constant, bottleneck, expansion). Within each array task, the three SLiM runs are launched in parallel (backgrounded), using **1 CPU each** (so request `--cpus-per-task=3`).

### 1) Basic submission

From the repo root:

```sbatch scripts/bgs_test_array.sbatch```

