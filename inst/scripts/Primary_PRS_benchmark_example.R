# Primary EUR/EAS PRS benchmark: EUR 50k/5k/5k and EAS 20k/5k/5k.
# Replace the four input paths below with the paths on your cluster.

library(simGO)

genotype_prefixes <- c(
  EUR = "/path/to/genotypes/EUR",
  EAS = "/path/to/genotypes/EAS"
)

phenotype_files <- c(
  EUR = "/path/to/phenotypes/EUR_simGO_phe.txt",
  EAS = "/path/to/phenotypes/EAS_simGO_phe.txt"
)

result <- prepare_primary_prs_benchmark(
  genotype_prefixes = genotype_prefixes,
  phenotype_files = phenotype_files,
  trait_names = c("Phe_Trait1", "Phe_Trait2"),
  output_dir = "Primary_PRS_benchmark",
  seed = 2026,
  plink = "plink"
)

print(result$design)
print(result$manifest)
cat("Run the generated script with:\n  bash", result$script, "\n")
