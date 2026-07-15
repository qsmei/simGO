# Genotype-only HAPGEN2 example. No phenotype is generated.

library(SimGO)

result <- simulate_1kg_hapgen2(
  reference_path = "/path/to/1KG_HAPGEN2_reference",
  output_path = "/path/to/simulated_genotypes",
  ancestries = c("EUR", "EAS"),
  chr_set = 1:22,
  rep_range = 1,
  sample_size = c(EUR = 60000, EAS = 30000),
  n_cases = 0,
  hapgen2 = "/path/to/hapgen2",
  genetic_map_path = "/path/to/1KG_HAPGEN2_reference/genetic_map",
  qc = TRUE,
  qc_output_file = "Hapgen2_plink_qc.sh",
  plink = "/path/to/plink",
  qc_maf = 0.01,
  qc_hwe = 1e-6,
  qc_geno = 0.05,
  qc_mind = 0.05,
  qc_merge = TRUE,
  output_file = "Hapgen2_genotype_simulation.sh",
  check_files = TRUE,
  run = FALSE
)

print(result)
cat("Run the generated scripts with:\n")
cat("  bash", result$hapgen2_script, "\n")
cat("  bash", result$qc_script, "\n")
