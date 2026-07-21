library(simGO)

work <- tempfile("simgo-prs-")
dir.create(work)
prefix <- file.path(work, "EUR")

writeLines(paste("F", 1:6, 0, 0, 0, -9), paste0(prefix, ".fam"))
file.create(paste0(prefix, ".bed"), paste0(prefix, ".bim"))
phenotype <- file.path(work, "phenotype.txt")
writeLines(c("FID IID Trait1", paste("F", 1:6, rnorm(6))), phenotype)

design <- prs_sample_design(
  total = c(EUR = 6),
  training = 4,
  tuning = 2,
  testing = 0
)

manifest <- create_prs_benchmark_splits(
  genotype_prefixes = c(EUR = prefix),
  sample_design = design,
  output_dir = file.path(work, "benchmark"),
  seed = 2026
)
stopifnot(
  sum(manifest$n) == 6,
  manifest$n[manifest$split == "testing"] == 0
)

script <- write_prs_benchmark_plink_script(
  genotype_prefixes = c(EUR = prefix),
  phenotype_files = c(EUR = phenotype),
  trait_names = "Trait1",
  sample_design = design,
  output_dir = file.path(work, "benchmark")
)
stopifnot(file.exists(script), any(grepl("--linear", readLines(script), fixed = TRUE)))
