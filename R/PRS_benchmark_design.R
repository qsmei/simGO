# Reproducible sample splitting and PLINK/GWAS preparation for PRS benchmarks.

simgo_validate_sample_design <- function(sample_design) {
  required <- c("ancestry", "total", "training", "tuning", "testing")
  missing_columns <- setdiff(required, names(sample_design))
  if (length(missing_columns) > 0) {
    stop("sample_design is missing: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }

  sample_design <- as.data.frame(sample_design, stringsAsFactors = FALSE)
  if (anyDuplicated(sample_design$ancestry)) {
    stop("Each ancestry must occur once in sample_design.", call. = FALSE)
  }
  size_columns <- c("total", "training", "tuning", "testing")
  sample_design[size_columns] <- lapply(sample_design[size_columns], as.integer)
  if (anyNA(sample_design[size_columns]) || any(unlist(sample_design[size_columns]) < 0)) {
    stop("All sample sizes must be non-negative integers.", call. = FALSE)
  }
  assigned <- sample_design$training + sample_design$tuning + sample_design$testing
  if (any(assigned != sample_design$total)) {
    bad <- sample_design$ancestry[assigned != sample_design$total]
    stop("training + tuning + testing must equal total for: ",
         paste(bad, collapse = ", "), call. = FALSE)
  }
  sample_design
}

simgo_primary_prs_design <- function() {
  prs_sample_design(
    total = c(EUR = 60000L, EAS = 30000L),
    training = c(EUR = 50000L, EAS = 20000L),
    tuning = 5000L,
    testing = 5000L
  )
}

prs_sample_design <- function(total,
                              training = NULL,
                              tuning = 5000L,
                              testing = tuning,
                              ancestry = names(total)) {
  # Construct a validated sample design for any number of ancestry groups.
  if (is.null(ancestry) || length(ancestry) == 0 || any(ancestry == "")) {
    stop("Supply ancestry names or provide a named total vector.", call. = FALSE)
  }
  ancestry <- as.character(ancestry)
  recycle_size <- function(x, name) {
    if (length(x) == 1L) {
      return(rep(as.integer(x), length(ancestry)))
    }
    if (!is.null(names(x))) {
      missing <- setdiff(ancestry, names(x))
      if (length(missing) > 0) {
        stop(name, " is missing: ", paste(missing, collapse = ", "), call. = FALSE)
      }
      return(as.integer(x[ancestry]))
    }
    if (length(x) != length(ancestry)) {
      stop(name, " must have length 1 or match ancestry.", call. = FALSE)
    }
    as.integer(x)
  }

  total <- recycle_size(total, "total")
  tuning <- recycle_size(tuning, "tuning")
  testing <- recycle_size(testing, "testing")
  if (is.null(training)) {
    training <- total - tuning - testing
  } else {
    training <- recycle_size(training, "training")
  }
  simgo_validate_sample_design(data.frame(
    ancestry = ancestry,
    total = total,
    training = training,
    tuning = tuning,
    testing = testing,
    stringsAsFactors = FALSE
  ))
}

simgo_read_fam_ids <- function(prefix) {
  fam_file <- paste0(prefix, ".fam")
  if (!file.exists(fam_file)) {
    stop("Missing PLINK FAM file: ", fam_file, call. = FALSE)
  }
  fam <- utils::read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(fam) < 2) {
    stop("FAM file must contain at least FID and IID: ", fam_file, call. = FALSE)
  }
  ids <- fam[, 1:2, drop = FALSE]
  names(ids) <- c("FID", "IID")
  if (anyDuplicated(paste(ids$FID, ids$IID, sep = "\r"))) {
    stop("Duplicated FID/IID pairs in: ", fam_file, call. = FALSE)
  }
  ids
}

create_prs_benchmark_splits <- function(genotype_prefixes,
                                        sample_design = simgo_primary_prs_design(),
                                        output_dir = "prs_benchmark",
                                        seed = 2026L) {
  sample_design <- simgo_validate_sample_design(sample_design)
  if (is.null(names(genotype_prefixes)) || any(names(genotype_prefixes) == "")) {
    stop("genotype_prefixes must be named by ancestry, e.g. c(EUR='...', EAS='...').",
         call. = FALSE)
  }
  missing_prefix <- setdiff(sample_design$ancestry, names(genotype_prefixes))
  if (length(missing_prefix) > 0) {
    stop("Missing genotype prefix for: ", paste(missing_prefix, collapse = ", "), call. = FALSE)
  }

  keep_dir <- file.path(output_dir, "keep")
  dir.create(keep_dir, recursive = TRUE, showWarnings = FALSE)
  manifest <- list()

  for (i in seq_len(nrow(sample_design))) {
    ancestry <- sample_design$ancestry[i]
    ids <- simgo_read_fam_ids(genotype_prefixes[[ancestry]])
    expected <- sample_design$total[i]
    if (nrow(ids) != expected) {
      stop(ancestry, " has ", nrow(ids), " individuals in its FAM file; expected ",
           expected, ".", call. = FALSE)
    }

    set.seed(as.integer(seed) + i - 1L)
    order_index <- sample.int(expected, expected, replace = FALSE)
    n_train <- sample_design$training[i]
    n_tune <- sample_design$tuning[i]
    train_index <- order_index[seq_len(n_train)]
    tune_index <- order_index[n_train + seq_len(n_tune)]
    test_index <- order_index[(n_train + n_tune + 1L):expected]
    split_indices <- list(training = train_index, tuning = tune_index, testing = test_index)

    for (split_name in names(split_indices)) {
      split_ids <- ids[sort(split_indices[[split_name]]), , drop = FALSE]
      keep_file <- file.path(keep_dir, paste0(ancestry, ".", split_name, ".keep"))
      utils::write.table(split_ids, keep_file, quote = FALSE, row.names = FALSE,
                         col.names = FALSE, sep = "\t")
      manifest[[length(manifest) + 1L]] <- data.frame(
        ancestry = ancestry,
        split = split_name,
        n = nrow(split_ids),
        keep_file = normalizePath(keep_file, mustWork = FALSE),
        stringsAsFactors = FALSE
      )
    }
  }

  manifest <- do.call(rbind, manifest)
  manifest_file <- file.path(output_dir, "sample_manifest.tsv")
  utils::write.table(manifest, manifest_file, quote = FALSE, row.names = FALSE, sep = "\t")
  invisible(manifest)
}

write_prs_benchmark_plink_script <- function(genotype_prefixes,
                                             phenotype_files,
                                             trait_names,
                                             sample_design = simgo_primary_prs_design(),
                                             output_dir = "prs_benchmark",
                                             plink = "plink",
                                             script_file = file.path(output_dir, "split_and_gwas.sh")) {
  sample_design <- simgo_validate_sample_design(sample_design)
  ancestries <- sample_design$ancestry
  if (!all(ancestries %in% names(genotype_prefixes))) {
    stop("genotype_prefixes must contain every ancestry in sample_design.", call. = FALSE)
  }
  if (!all(ancestries %in% names(phenotype_files))) {
    stop("phenotype_files must contain every ancestry in sample_design.", call. = FALSE)
  }
  if (length(trait_names) == 0 || anyNA(trait_names) || any(trait_names == "")) {
    stop("trait_names must contain at least one phenotype column name.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_abs <- normalizePath(output_dir, mustWork = FALSE)
  lines <- c("#!/usr/bin/env bash", "set -euo pipefail", "",
             paste0("PLINK=", shQuote(plink)),
             paste0("ROOT=", shQuote(output_abs)),
             "mkdir -p \"${ROOT}/geno\" \"${ROOT}/sumstats\"", "")

  for (ancestry in ancestries) {
    bfile <- normalizePath(genotype_prefixes[[ancestry]], mustWork = FALSE)
    pheno <- normalizePath(phenotype_files[[ancestry]], mustWork = FALSE)
    lines <- c(lines,
      paste0("echo ", shQuote(paste("Preparing", ancestry))),
      paste0("BFILE=", shQuote(bfile)),
      paste0("PHENO=", shQuote(pheno)))
    for (split_name in c("training", "tuning", "testing")) {
      keep <- file.path(output_abs, "keep", paste0(ancestry, ".", split_name, ".keep"))
      out <- paste0("${ROOT}/geno/", ancestry, ".", split_name)
      lines <- c(lines,
        paste0("\"${PLINK}\" --bfile \"${BFILE}\" --keep ", shQuote(keep),
               " --keep-allele-order --make-bed --out \"", out, "\""))
    }
    for (trait in trait_names) {
      out <- paste0("${ROOT}/sumstats/", ancestry, ".", trait)
      lines <- c(lines,
        paste0("\"${PLINK}\" --bfile \"${ROOT}/geno/", ancestry,
               ".training\" --pheno \"${PHENO}\" --pheno-name ", shQuote(trait),
               " --linear hide-covar --allow-no-sex --out \"", out, "\""))
    }
    lines <- c(lines, "")
  }

  writeLines(lines, script_file)
  Sys.chmod(script_file, mode = "0755")
  invisible(normalizePath(script_file, mustWork = FALSE))
}

prepare_prs_benchmark <- function(genotype_prefixes,
                                  phenotype_files,
                                  trait_names,
                                  sample_design,
                                  output_dir = "prs_benchmark",
                                  seed = 2026L,
                                  plink = "plink") {
  # General package entry point for ancestry- and trait-agnostic preparation.
  sample_design <- simgo_validate_sample_design(sample_design)
  manifest <- create_prs_benchmark_splits(
    genotype_prefixes = genotype_prefixes,
    sample_design = sample_design,
    output_dir = output_dir,
    seed = seed
  )
  script <- write_prs_benchmark_plink_script(
    genotype_prefixes = genotype_prefixes,
    phenotype_files = phenotype_files,
    trait_names = trait_names,
    sample_design = sample_design,
    output_dir = output_dir,
    plink = plink
  )
  invisible(list(design = sample_design, manifest = manifest, script = script))
}

prepare_primary_prs_benchmark <- function(genotype_prefixes,
                                          phenotype_files,
                                          trait_names,
                                          output_dir = "prs_benchmark",
                                          seed = 2026L,
                                          plink = "plink") {
  prepare_prs_benchmark(
    genotype_prefixes = genotype_prefixes,
    phenotype_files = phenotype_files,
    trait_names = trait_names,
    sample_design = simgo_primary_prs_design(),
    output_dir = output_dir,
    seed = seed,
    plink = plink
  )
}
