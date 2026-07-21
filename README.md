# simGO

`simGO` provides reproducible workflows for:

- simulating ancestry-specific genotypes from 1000 Genomes-style HAPGEN2 references;
- converting and quality-controlling HAPGEN2 output with PLINK 1.x; and
- preparing ancestry-specific train, tuning, and test datasets for PRS benchmarks.

Older experimental phenotype and omics prototypes are retained under
`inst/legacy/` for reference, but they are not loaded as package functions.

## Installation

```r
remotes::install_github("qsmei/simGO")
library(simGO)
```

HAPGEN2 and PLINK are external programs. If they are available on `PATH`, use
`hapgen2 = "hapgen2"` and `plink = "plink"`. Otherwise provide their absolute
executable paths.

## HAPGEN2 genotype simulation

```r
result <- simulate_1kg_hapgen2(
  # Root containing EUR/chr1.hap, EUR/chr1.legend, EAS/..., and genetic_map/.
  reference_path =
    "/restricted/projectnb/fhs_data/qsmei/Phase3_1000G_imputed/1000GP_Phase3/Hapmap3_1000G",

  # HAPGEN2 outputs are organized as rep1/EUR/, rep1/EAS/, and so on.
  output_path =
    "/restricted/projectnb/fhs_data/qsmei/Phase3_1000G_imputed/1000GP_Phase3/Hapmap3_1000G/simulation/hapgen2",

  ancestries = c("EUR", "EAS"),
  chr_set = 1,

  # One genotype replication is usually enough when causal effects and
  # phenotypes are independently replicated for every simulation setting.
  rep_range = 1,
  sample_size = c(EUR = 60000, EAS = 30000),
  n_cases = 0,

  hapgen2 = "hapgen2",
  plink = "plink",

  qc = TRUE,
  qc_maf = 0.01,   # Keep variants with MAF >= 1%.
  qc_hwe = 1e-6,   # Remove variants with HWE P < 1e-6.
  qc_geno = 0.05,  # Remove variants with >5% missing genotypes.
  qc_mind = 0.05,  # Remove individuals with >5% missing genotypes.
  qc_merge = FALSE, # Nothing needs merging for chromosome 1 only.

  check_files = TRUE,
  run = FALSE       # Write scripts for inspection; execute them later.
)

result$hapgen2_script
result$qc_script
```

Run the generated scripts on the cluster:

```bash
bash Hapgen2.sh
bash hapgen2_plink_qc.sh
```

For chromosomes 1–22, `qc_merge = TRUE` creates one PLINK dataset per ancestry
and replication. It combines chromosome columns, never EUR and EAS individuals.

## PRS sample design

```r
design <- prs_sample_design(
  total = c(EUR = 60000, EAS = 30000),
  training = c(EUR = 50000, EAS = 20000),
  tuning = 5000,
  testing = 5000
)

result <- prepare_prs_benchmark(
  genotype_prefixes = c(EUR = "geno/EUR", EAS = "geno/EAS"),
  phenotype_files = c(EUR = "phe/EUR.txt", EAS = "phe/EAS.txt"),
  trait_names = "Phe_Trait1",
  sample_design = design,
  output_dir = "prs_benchmark",
  seed = 2026,
  plink = "plink"
)
```
