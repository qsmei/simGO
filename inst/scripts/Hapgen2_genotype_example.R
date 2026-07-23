# Genotype-only HAPGEN2 example. No phenotype is generated.

library(simGO)

reference_path <- paste0(
  "/restricted/projectnb/fhs_data/qsmei/Phase3_1000G_imputed/",
  "1000GP_Phase3/Hapmap3_1000G"
)

result <- simulate_1kg_hapgen2(
  reference_path = reference_path,
  output_path = file.path(reference_path, "simulation", "hapgen2"),
  ancestries = c("EUR", "EAS"),
  chr_set = 1,
  rep_range = 1,
  sample_size = c(EUR = 60000, EAS = 30000),
  n_cases = 0,
  hapgen2 = "hapgen2",
  hap_prefix = "Hapmap3_1000GP_Phase3_chr",
  legend_path = file.path(
    reference_path,
    "legend"
  ),
  legend_prefix = "Hapmap3_1000GP_Phase3_chr",
  map_prefix = "genetic_map_chr",
  map_suffix = "_combined_b37.txt",
  qc = TRUE,
  plink = "plink",
  qc_maf = 0.01,   # Keep SNPs with minor-allele frequency >= 1%.
  qc_hwe = 1e-6,   # Remove SNPs with HWE P-value < 1e-6.
  qc_geno = 0.05,  # Remove SNPs with >5% missing genotypes.
  qc_mind = 0.05,  # Remove individuals with >5% missing genotypes.
  qc_merge = FALSE, # Nothing needs merging for chromosome 1 only.
  check_files = TRUE,
  run = FALSE
)

print(result)
cat("Run the generated scripts with:\n")
cat("  bash", result$hapgen2_script, "\n")
cat("  bash", result$qc_script, "\n")
