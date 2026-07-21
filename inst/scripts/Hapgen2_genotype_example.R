# Genotype-only HAPGEN2 example. No phenotype is generated.

library(simGO)

result <- simulate_1kg_hapgen2(
  reference_path = paste0(
    "/restricted/projectnb/fhs_data/qsmei/Phase3_1000G_imputed/",
    "1000GP_Phase3/Hapmap3_1000G"
  ),
  output_path = paste0(
    "/restricted/projectnb/fhs_data/qsmei/Phase3_1000G_imputed/",
    "1000GP_Phase3/Hapmap3_1000G/simulation/hapgen2"
  ),
  ancestries = c("EUR", "EAS"),
  chr_set = 1,
  rep_range = 1,
  sample_size = c(EUR = 60000, EAS = 30000),
  n_cases = 0,
  hapgen2 = "hapgen2",
  qc = TRUE,
  qc_output_file = "Hapgen2_plink_qc.sh",
  plink = "plink",
  qc_maf = 0.01,   # Keep SNPs with minor-allele frequency >= 1%.
  qc_hwe = 1e-6,   # Remove SNPs with HWE P-value < 1e-6.
  qc_geno = 0.05,  # Remove SNPs with >5% missing genotypes.
  qc_mind = 0.05,  # Remove individuals with >5% missing genotypes.
  qc_merge = FALSE, # Nothing needs merging for chromosome 1 only.
  output_file = "Hapgen2_genotype_simulation.sh",
  check_files = TRUE,
  run = FALSE
)

print(result)
cat("Run the generated scripts with:\n")
cat("  bash", result$hapgen2_script, "\n")
cat("  bash", result$qc_script, "\n")
