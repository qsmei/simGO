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
  output_file = file.path(work, "hapgen2.sh"),
  qc_output_file = file.path(work, "qc.sh"),
  check_files = FALSE,
  run = FALSE
)

stopifnot(
  file.exists(one_chr$hapgen2_script),
  file.exists(one_chr$qc_script),
  identical(one_chr$executed, FALSE)
)

qc_lines <- readLines(one_chr$qc_script)
stopifnot(
  any(grepl("--maf.*0.01", qc_lines)),
  any(grepl("--geno.*0.05", qc_lines)),
  !any(grepl("sed -E", qc_lines, fixed = TRUE)),
  !any(grepl("--merge-list", qc_lines, fixed = TRUE))
)

hapgen2_lines <- readLines(one_chr$hapgen2_script)
stopifnot(
  any(grepl("LEGEND_PATH=", hapgen2_lines, fixed = TRUE)),
  any(grepl("MAP_SUFFIX='.map'", hapgen2_lines, fixed = TRUE))
)

multi_chr_qc <- run_hapgen2_plink_qc(
  genotype_path = file.path(work, "genotypes"),
  ancestries = "EUR",
  chr_set = 1:2,
  rep_range = 1,
  output_file = file.path(work, "qc-multi.sh"),
  run = FALSE
)
stopifnot(any(grepl("--merge-list", readLines(multi_chr_qc), fixed = TRUE)))

invalid_maf <- try(
  run_hapgen2_plink_qc(
    genotype_path = work,
    maf = 0.8,
    output_file = file.path(work, "invalid.sh")
  ),
  silent = TRUE
)
stopifnot(inherits(invalid_maf, "try-error"))
