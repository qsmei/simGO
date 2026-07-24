library(simGO)

work <- tempfile("simgo-hapgen2-")
dir.create(work)

one_chr <- simulate_1kg_hapgen2(
  reference_path = file.path(work, "reference"),
  output_path = file.path(work, "genotypes"),
  ancestries = c("EUR", "EAS"),
  chr_set = 1,
  rep_range = 1,
  sample_size = c(EUR = 100, EAS = 50),
  legend_path = file.path(work, "legend"),
  map_prefix = "chr",
  map_suffix = ".map",
  qc = TRUE,
  scripts_path = file.path(work, "custom-scripts"),
  job_scheduler = "slurm",
  job_parameters = list(
    project = "test-project",
    memory = "16G",
    cpus = 2,
    walltime = "12:00:00",
    queue = "short",
    job_name = "test",
    log_path = file.path(work, "logs")
  ),
  check_files = FALSE,
  run = FALSE
)

stopifnot(
  identical(names(one_chr), c("output_path", "qc_output_path", "scripts_path")),
  file.exists(file.path(one_chr$scripts_path, "hapgen2_simGO.sh")),
  file.exists(file.path(one_chr$scripts_path, "hapgen2_simGO_QC.sh")),
  dir.exists(one_chr$output_path),
  identical(basename(one_chr$scripts_path), "custom-scripts"),
  dir.exists(file.path(work, "genotypes", "qc")),
  identical(
    normalizePath(one_chr$qc_output_path, mustWork = FALSE),
    normalizePath(file.path(work, "genotypes", "qc"), mustWork = FALSE)
  )
)

qc_lines <- readLines(file.path(one_chr$scripts_path, "hapgen2_simGO_QC.sh"))
stopifnot(
  any(grepl("QC_OUTPUT_PATH=", qc_lines, fixed = TRUE)),
  any(grepl("--maf.*0.01", qc_lines)),
  any(grepl("--geno.*0.05", qc_lines)),
  !any(grepl("sed -E", qc_lines, fixed = TRUE)),
  !any(grepl("--merge-list", qc_lines, fixed = TRUE))
)

hapgen2_lines <- readLines(file.path(one_chr$scripts_path, "hapgen2_simGO.sh"))
stopifnot(
  any(grepl("#SBATCH --job-name=test_hapgen2", hapgen2_lines, fixed = TRUE)),
  any(grepl("#SBATCH --account=test-project", hapgen2_lines, fixed = TRUE)),
  any(grepl("#SBATCH --mem=16G", hapgen2_lines, fixed = TRUE)),
  any(grepl("#SBATCH --cpus-per-task=2", hapgen2_lines, fixed = TRUE)),
  any(grepl("LEGEND_PATH=", hapgen2_lines, fixed = TRUE)),
  any(grepl("MAP_SUFFIX='.map'", hapgen2_lines, fixed = TRUE))
)

multi_chr_qc <- run_hapgen2_plink_qc(
  genotype_path = file.path(work, "genotypes"),
  ancestries = "EUR",
  chr_set = 1:2,
  rep_range = 1,
  scripts_path = file.path(work, "standalone-scripts"),
  run = FALSE
)
stopifnot(any(grepl("--merge-list", readLines(multi_chr_qc), fixed = TRUE)))

custom_work <- tempfile("simgo-custom-header-")
custom_result <- simulate_1kg_hapgen2(
  reference_path = file.path(custom_work, "reference"),
  output_path = file.path(custom_work, "output"),
  ancestries = "EUR",
  chr_set = 1,
  sample_size = 10,
  qc = TRUE,
  job_scheduler = "custom",
  job_parameters = "
#!/bin/bash
#$ -P sequencing
#$ -N simGO_{job_type}
#$ -j y
#$ -l mem_per_core=50G
#$ -pe omp 4
",
  check_files = FALSE
)
custom_hapgen2 <- readLines(file.path(custom_result$scripts_path, "hapgen2_simGO.sh"))
custom_qc <- readLines(file.path(custom_result$scripts_path, "hapgen2_simGO_QC.sh"))
stopifnot(
  sum(grepl("^#!", custom_hapgen2)) == 1L,
  any(custom_hapgen2 == "#$ -N simGO_hapgen2"),
  any(custom_qc == "#$ -N simGO_QC"),
  any(custom_hapgen2 == "#$ -l mem_per_core=50G")
)

invalid_maf <- try(
  run_hapgen2_plink_qc(
    genotype_path = work,
    maf = 0.8,
    scripts_path = file.path(work, "invalid-scripts")
  ),
  silent = TRUE
)
stopifnot(inherits(invalid_maf, "try-error"))
