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

Note: The SBATCH script also activates this environment internally using the same YAML-defined environment. 

### 1) Basic submission

From the repo root:

```sbatch scripts/bgs_test_array.sbatch```

Overriding parameters at submit time:

```sbatch --export=ALL,<FLAG>=<value> scripts/bgs_test_array.sbatch```

## Configurable SBATCH flags

The SLURM array script (`scripts/bgs_test_array.sbatch`) exposes user-configurable flags that can be overridden at submission time using:

### Core simulation flags

| Flag | Default | Role | One-line explanation |
|------|---------|------|----------------------|
| `PARAM_FILE` | `priors/priors_2.6.25.csv` | Fixed (input) | CSV file containing one row per simulation with the parameters to be inferred |
| `END` | `20000` | Fixed (experimental) | Total number of generations simulated (controls time to quasi-equilibrium) |
| `OFFSET` | `0` | Execution | Number of rows to skip in the priors file (used to resume or shard runs) |
| `SIMPLIFY_INTERVAL` | `200` | Fixed (technical) | How often tree sequences are simplified to reduce memory usage |
| `SELF` | `0.98` | Fixed (biological) | Selfing rate applied uniformly across all simulations |

---

### Inference target parameters (from priors file)

These parameters are **jointly inferred via ABC / ABC-RF**.  
Each array task reads exactly one row from `PARAM_FILE`.

| Parameter | Role | One-line explanation |
|----------|------|----------------------|
| `gmu` | Inference target | Mutation rate in genic regions |
| `imu` | Inference target | Mutation rate in intergenic regions |
| `gd` | Inference target | Fraction of genic sites under selection |
| `id` | Inference target | Fraction of intergenic sites under selection |
| `gdfe` | Inference target | Mean deleterious effect size for genic mutations |
| `idfe` | Inference target | Mean deleterious effect size for intergenic mutations |

---

### Bottleneck demography flags  
*(used only by `models/ABCtrees_bottleneck.slim`)*

| Flag | Default | Role | One-line explanation |
|------|---------|------|----------------------|
| `BOT_FRAC` | `0.2` | Fixed (experimental) | Fraction of population size retained during the bottleneck |
| `BOT_FRAC_TBOT` | `0.50` | Fixed (experimental) | Timing of bottleneck onset as a fraction of `END` |
| `BOT_FRAC_TREC` | `0.75` | Fixed (experimental) | Timing of population recovery as a fraction of `END` |

**Internally converted to:**
- `T_BOT = round(END × BOT_FRAC_TBOT)`
- `T_REC = round(END × BOT_FRAC_TREC)`

---

### Expansion demography flags  
*(used only by `models/ABCtrees_expansion.slim`)*

| Flag | Default | Role | One-line explanation |
|------|---------|------|----------------------|
| `EXP_MULT` | `3.0` | Fixed (experimental) | Population size multiplier applied at expansion |
| `EXP_FRAC_TEXP` | `0.75` | Fixed (experimental) | Timing of expansion as a fraction of `END` |

**Internally converted to:**
- `T_EXP = round(END × EXP_FRAC_TEXP)`


