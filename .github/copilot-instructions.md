# AI Coding Agent Instructions for ABCslim_trees

## Overview
This project simulates evolutionary scenarios using SLiM tree-sequence recording and processes the resulting data into VCF files. It includes Python scripts for post-simulation processing and SLiM scripts for defining evolutionary models.

### Key Components
1. **SLiM Scripts**:
   - `ABCtrees_bottleneck.slim`: Simulates a bottleneck and recovery scenario.
   - `ABCtrees_expansion.slim`: Simulates a population expansion scenario.
   - `ABCtrees.slim`: General-purpose simulation script.
2. **Python Script**:
   - `trees_to_vcf.py`: Converts `.trees` files from SLiM simulations into `.vcf` files, applying mutation models and filtering.
3. **Genome Information**:
   - `genomeinfo/gene100.csv` and `genomeinfo/intergene100.csv`: Define genomic intervals for genes and intergenic regions.

## Developer Workflows

### Running Simulations
- Use SLiM to execute simulation scripts. Example:
  ```bash
  slim -d ID=1 -d gmu=1e-8 -d imu=1e-9 -d gd=0.5 -d id=0.5 ABCtrees_bottleneck.slim
  ```
- Key parameters:
  - `ID`: Unique identifier for the simulation.
  - `gmu`, `imu`: Mutation rates for genic and intergenic regions.
  - `gd`, `id`: Proportions of deleterious mutations.

### Converting Trees to VCF
- Run `trees_to_vcf.py` to process `.trees` files:
  ```bash
  python3 trees_to_vcf.py --trees output.trees --out_vcf output.vcf \
      --gene_csv genomeinfo/gene100.csv --intergene_csv genomeinfo/intergene100.csv \
      --L 1000000 --gmu 1e-8 --imu 1e-9 --gd 0.5 --id 0.5
  ```
- Optional flags:
  - `--recapitate`: Adds deep ancestry to the tree sequence.
  - `--biallelic_only`: Filters out multi-allelic sites.

## Project-Specific Conventions
- **SLiM Scripts**:
  - Use `defineConstant` for configurable parameters.
  - Ensure `initializeTreeSeq` is called for tree-sequence recording.
  - Validate parameter ranges (e.g., `END <= 20000`).
- **Python Scripts**:
  - Use `pyslim` for loading `.trees` files and `msprime` for mutation overlays.
  - Define mutation rates separately for genic and intergenic regions.

## Integration Points
- **SLiM and Python**:
  - SLiM outputs `.trees` files, which are processed by `trees_to_vcf.py`.
- **Genome Information**:
  - `gene100.csv` and `intergene100.csv` are required for defining genomic intervals.

## Examples
### SLiM Bottleneck Simulation
```bash
slim -d ID=42 -d gmu=1e-8 -d imu=1e-9 -d gd=0.3 -d id=0.7 ABCtrees_bottleneck.slim
```

### Convert Trees to VCF
```bash
python3 trees_to_vcf.py --trees 42.trees --out_vcf 42.vcf \
    --gene_csv genomeinfo/gene100.csv --intergene_csv genomeinfo/intergene100.csv \
    --L 1000000 --gmu 1e-8 --imu 1e-9 --gd 0.3 --id 0.7
```

## Notes
- Ensure SLiM and Python dependencies (`pyslim`, `msprime`, `pandas`, `numpy`) are installed.
- Use `SIMPLIFY_INTERVAL` in SLiM scripts to manage memory usage during long simulations.
- Validate `.csv` files for genomic intervals to avoid runtime errors.