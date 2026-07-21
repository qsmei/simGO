# HAPGEN2 / 1000 Genomes simulation helpers.

.simgo_validate_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  x
}

.simgo_validate_labels <- function(x, name) {
  if (!is.character(x) || length(x) == 0L || anyNA(x) || any(x == "")) {
    stop(name, " must be a non-empty character vector.", call. = FALSE)
  }
  if (anyDuplicated(x)) {
    stop(name, " must not contain duplicates.", call. = FALSE)
  }
  if (any(!grepl("^[A-Za-z0-9_.-]+$", x))) {
    stop(name, " may contain only letters, numbers, '.', '_' and '-'.", call. = FALSE)
  }
  x
}

.simgo_validate_integer_vector <- function(x, name, minimum = 1L, maximum = Inf) {
  if (!is.numeric(x) || length(x) == 0L || anyNA(x) ||
      any(!is.finite(x)) || any(x != floor(x)) ||
      any(x < minimum) || any(x > maximum)) {
    stop(name, " must contain integers between ", minimum, " and ", maximum, ".",
         call. = FALSE)
  }
  as.integer(x)
}

.simgo_validate_qc_threshold <- function(x, name, maximum = 1) {
  if (is.null(x)) {
    return(NULL)
  }
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x < 0 || x > maximum) {
    stop(name, " must be NULL or a number between 0 and ", maximum, ".",
         call. = FALSE)
  }
  as.numeric(x)
}

.simgo_named_integer <- function(x, labels, name, minimum = 1L) {
  if (length(x) == 1L) {
    x <- stats::setNames(rep(x, length(labels)), labels)
  } else if (is.null(names(x))) {
    if (length(x) != length(labels)) {
      stop(name, " must have length 1 or match ancestries.", call. = FALSE)
    }
    names(x) <- labels
  }
  if (anyDuplicated(names(x))) {
    stop(name, " must not contain duplicated names.", call. = FALSE)
  }
  missing_labels <- setdiff(labels, names(x))
  if (length(missing_labels) > 0L) {
    stop(name, " is missing: ", paste(missing_labels, collapse = ", "), call. = FALSE)
  }
  values <- .simgo_validate_integer_vector(x[labels], name, minimum = minimum)
  stats::setNames(values, labels)
}

.simgo_prepare_script_path <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || path == "") {
    stop("Script output path must be one non-empty string.", call. = FALSE)
  }
  path <- path.expand(path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = FALSE)
}

.simgo_check_executable <- function(command, name) {
  if (!is.character(command) || length(command) != 1L || is.na(command) || command == "") {
    stop(name, " must be an executable name or path.", call. = FALSE)
  }
  available <- if (grepl("/", command, fixed = TRUE)) {
    file.exists(command) && file.access(command, mode = 1L) == 0L
  } else {
    nzchar(Sys.which(command))
  }
  if (!available) {
    stop(name, " executable was not found: ", command, call. = FALSE)
  }
  invisible(TRUE)
}

.simgo_run_script <- function(script_file, label) {
  status <- system2("bash", args = shQuote(script_file))
  if (!identical(status, 0L)) {
    stop(label, " failed with status ", status, ".", call. = FALSE)
  }
  invisible(TRUE)
}

#' Default HAPGEN2 effective population sizes.
hapgen2_ne_defaults <- function() {
  c(EUR = 11418L, EAS = 14269L, AMR = 11418L, SAS = 14269L, AFR = 17469L)
}

validate_hapgen2_ancestries <- function(ancestries, ne) {
  if (is.null(names(ne)) || anyNA(names(ne)) || any(names(ne) == "")) {
    stop("ne must be a named vector.", call. = FALSE)
  }
  missing_ne <- setdiff(ancestries, names(ne))
  if (length(missing_ne) > 0L) {
    stop("Missing effective population size (Ne) for: ",
         paste(missing_ne, collapse = ", "), call. = FALSE)
  }
  ne_values <- .simgo_validate_integer_vector(ne[ancestries], "ne", minimum = 1L)
  stats::setNames(ne_values, ancestries)
}

write_1kg_hapgen2_script <- function(reference_path,
                                     output_path,
                                     ancestries = c("EUR", "EAS"),
                                     chr_set = 1:22,
                                     rep_range = 1L,
                                     sample_size = 10000L,
                                     n_cases = 0L,
                                     hapgen2 = "hapgen2",
                                     genetic_map_path = file.path(reference_path, "genetic_map"),
                                     hap_prefix = "chr",
                                     legend_prefix = hap_prefix,
                                     map_prefix = "genetic_map_chr",
                                     map_suffix = "_combined_b37.txt",
                                     output_prefix = NULL,
                                     ne = hapgen2_ne_defaults(),
                                     no_haps_output = TRUE,
                                     output_file = "Hapgen2.sh",
                                     check_files = FALSE) {
  if (missing(reference_path) || missing(output_path)) {
    stop("reference_path and output_path must be provided.", call. = FALSE)
  }
  ancestries <- .simgo_validate_labels(ancestries, "ancestries")
  chr_set <- .simgo_validate_integer_vector(chr_set, "chr_set", 1L, 22L)
  rep_range <- .simgo_validate_integer_vector(rep_range, "rep_range", 1L)
  sample_size <- .simgo_named_integer(sample_size, ancestries, "sample_size", 1L)
  n_cases <- .simgo_validate_integer_vector(n_cases, "n_cases", 0L)[1L]
  ne <- validate_hapgen2_ancestries(ancestries, ne)
  no_haps_output <- .simgo_validate_flag(no_haps_output, "no_haps_output")
  check_files <- .simgo_validate_flag(check_files, "check_files")

  reference_path <- path.expand(reference_path)
  output_path <- path.expand(output_path)
  genetic_map_path <- path.expand(genetic_map_path)
  output_file <- .simgo_prepare_script_path(output_file)

  if (check_files) {
    required_dirs <- c(reference_path, genetic_map_path)
    missing_dirs <- required_dirs[!dir.exists(required_dirs)]
    if (length(missing_dirs) > 0L) {
      stop("Missing directory: ", paste(missing_dirs, collapse = ", "), call. = FALSE)
    }

    required_files <- unlist(lapply(chr_set, function(chr) {
      c(
        file.path(genetic_map_path, paste0(map_prefix, chr, map_suffix)),
        unlist(lapply(ancestries, function(ancestry) {
          c(
            file.path(reference_path, ancestry, paste0(hap_prefix, chr, ".hap")),
            file.path(reference_path, ancestry, paste0(legend_prefix, chr, ".legend"))
          )
        }), use.names = FALSE)
      )
    }), use.names = FALSE)
    missing_files <- required_files[!file.exists(required_files)]
    if (length(missing_files) > 0L) {
      shown <- utils::head(missing_files, 20L)
      stop(
        "Missing HAPGEN2 input file(s):\n",
        paste(shown, collapse = "\n"),
        if (length(missing_files) > length(shown)) "\n..." else "",
        call. = FALSE
      )
    }
  }

  prefix_value <- if (is.null(output_prefix) || identical(output_prefix, "")) {
    ""
  } else {
    if (length(output_prefix) != 1L) {
      stop("output_prefix must be one string.", call. = FALSE)
    }
    .simgo_validate_labels(output_prefix, "output_prefix")
    paste0(output_prefix, "_")
  }
  output_line <- if (no_haps_output) {
    "        -o \"${out_prefix}\" -no_haps_output"
  } else {
    "        -o \"${out_prefix}\""
  }
  sample_size_cases <- paste(
    sprintf("        %s) sample_n=%s ;;", shQuote(names(sample_size)), sample_size),
    collapse = "\n"
  )
  ne_cases <- paste(
    sprintf("        %s) Ne=%s ;;", shQuote(names(ne)), ne),
    collapse = "\n"
  )

  script <- c(
    "#!/usr/bin/env bash",
    "set -euo pipefail",
    "",
    paste0("REFERENCE_PATH=", shQuote(reference_path)),
    paste0("OUTPUT_PATH=", shQuote(output_path)),
    paste0("GENETIC_MAP_PATH=", shQuote(genetic_map_path)),
    paste0("HAPGEN2=", shQuote(hapgen2)),
    paste0("HAP_PREFIX=", shQuote(hap_prefix)),
    paste0("LEGEND_PREFIX=", shQuote(legend_prefix)),
    paste0("MAP_PREFIX=", shQuote(map_prefix)),
    paste0("MAP_SUFFIX=", shQuote(map_suffix)),
    paste0("OUTPUT_PREFIX=", shQuote(prefix_value)),
    paste0("N_CASES=", n_cases),
    paste0("ANCESTRIES=(", paste(shQuote(ancestries), collapse = " "), ")"),
    paste0("CHROMOSOMES=(", paste(chr_set, collapse = " "), ")"),
    paste0("REPLICATIONS=(", paste(rep_range, collapse = " "), ")"),
    "",
    "for i_rep in \"${REPLICATIONS[@]}\"; do",
    "  for i_chr in \"${CHROMOSOMES[@]}\"; do",
    "    map_file=\"${GENETIC_MAP_PATH}/${MAP_PREFIX}${i_chr}${MAP_SUFFIX}\"",
    "    for i_ancestry in \"${ANCESTRIES[@]}\"; do",
    "      case \"${i_ancestry}\" in",
    ne_cases,
    "        *) echo \"No Ne configured for ${i_ancestry}\" >&2; exit 1 ;;",
    "      esac",
    "      case \"${i_ancestry}\" in",
    sample_size_cases,
    "        *) echo \"No sample size configured for ${i_ancestry}\" >&2; exit 1 ;;",
    "      esac",
    "",
    "      hap_file=\"${REFERENCE_PATH}/${i_ancestry}/${HAP_PREFIX}${i_chr}.hap\"",
    "      legend_file=\"${REFERENCE_PATH}/${i_ancestry}/${LEGEND_PREFIX}${i_chr}.legend\"",
    "      out_dir=\"${OUTPUT_PATH}/rep${i_rep}/${i_ancestry}\"",
    "      out_prefix=\"${out_dir}/${OUTPUT_PREFIX}${i_ancestry}_chr${i_chr}\"",
    "      mkdir -p \"${out_dir}\"",
    "",
    "      dl_initial=$(awk 'NR==2{print $2}' \"${legend_file}\")",
    "      if [ -z \"${dl_initial}\" ]; then",
    "        echo \"Could not read a disease-locus position from ${legend_file}\" >&2",
    "        exit 1",
    "      fi",
    "      echo \"HAPGEN2 ancestry=${i_ancestry} chr=${i_chr} rep=${i_rep} n=${sample_n}\"",
    "      \"${HAPGEN2}\" -dl \"${dl_initial}\" 1 1 1 \\",
    "        -h \"${hap_file}\" \\",
    "        -m \"${map_file}\" \\",
    "        -l \"${legend_file}\" \\",
    "        -Ne \"${Ne}\" \\",
    "        -n \"${sample_n}\" \"${N_CASES}\" \\",
    output_line,
    "    done",
    "  done",
    "done",
    ""
  )

  writeLines(script, output_file)
  Sys.chmod(output_file, mode = "0755")
  invisible(output_file)
}

write_hapgen2_plink_qc_script <- function(genotype_path,
                                          ancestries = c("EUR", "EAS"),
                                          chr_set = 1:22,
                                          rep_range = 1L,
                                          prefix = NULL,
                                          plink = "plink",
                                          maf = 0.01,
                                          hwe = 1e-6,
                                          geno = 0.05,
                                          mind = 0.05,
                                          merge = length(chr_set) > 1L,
                                          output_file = "hapgen2_plink_qc.sh") {
  if (missing(genotype_path)) {
    stop("genotype_path must be provided.", call. = FALSE)
  }
  ancestries <- .simgo_validate_labels(ancestries, "ancestries")
  chr_set <- .simgo_validate_integer_vector(chr_set, "chr_set", 1L, 22L)
  if (!is.null(rep_range)) {
    rep_range <- .simgo_validate_integer_vector(rep_range, "rep_range", 1L)
  }
  maf <- .simgo_validate_qc_threshold(maf, "maf", maximum = 0.5)
  hwe <- .simgo_validate_qc_threshold(hwe, "hwe")
  geno <- .simgo_validate_qc_threshold(geno, "geno")
  mind <- .simgo_validate_qc_threshold(mind, "mind")
  merge <- .simgo_validate_flag(merge, "merge")
  output_file <- .simgo_prepare_script_path(output_file)

  prefix_value <- if (is.null(prefix) || identical(prefix, "")) {
    ""
  } else {
    if (length(prefix) != 1L) {
      stop("prefix must be one string.", call. = FALSE)
    }
    .simgo_validate_labels(prefix, "prefix")
    paste0(prefix, "_")
  }
  rep_dirs <- if (is.null(rep_range)) {
    "REPLICATION_DIRS=(\"\")"
  } else {
    paste0("REPLICATION_DIRS=(", paste(shQuote(paste0("rep", rep_range)), collapse = " "), ")")
  }
  qc_args <- c(
    if (!is.null(maf)) c("--maf", format(maf, scientific = FALSE, trim = TRUE)),
    if (!is.null(hwe)) c("--hwe", format(hwe, scientific = TRUE, trim = TRUE)),
    if (!is.null(geno)) c("--geno", format(geno, scientific = FALSE, trim = TRUE)),
    if (!is.null(mind)) c("--mind", format(mind, scientific = FALSE, trim = TRUE))
  )
  qc_arg_line <- paste(shQuote(qc_args), collapse = " ")

  append_to_merge_list <- if (merge) {
    "      echo \"${output_prefix}\" >> \"${merge_list}\""
  } else {
    character()
  }
  merge_block <- if (merge) {
    c(
      "    if [ ! -s \"${merge_list}\" ]; then",
      "      echo \"No chromosome datasets were available to merge in ${ancestry_dir}\" >&2",
      "      exit 1",
      "    fi",
      "    \"${PLINK}\" --merge-list \"${merge_list}\" --make-bed --allow-no-sex --out \"${PREFIX}${i_ancestry}\"",
      "    mv \"${merge_list}\" chr_file_list"
    )
  } else {
    character()
  }

  script <- c(
    "#!/usr/bin/env bash",
    "set -euo pipefail",
    "",
    paste0("GENOTYPE_PATH=", shQuote(path.expand(genotype_path))),
    paste0("PLINK=", shQuote(plink)),
    paste0("PREFIX=", shQuote(prefix_value)),
    paste0("QC_ARGS=(", qc_arg_line, ")"),
    rep_dirs,
    paste0("ANCESTRIES=(", paste(shQuote(ancestries), collapse = " "), ")"),
    paste0("CHROMOSOMES=(", paste(chr_set, collapse = " "), ")"),
    "",
    "for i_rep_dir in \"${REPLICATION_DIRS[@]}\"; do",
    "  for i_ancestry in \"${ANCESTRIES[@]}\"; do",
    "    if [ -n \"${i_rep_dir}\" ]; then",
    "      ancestry_dir=\"${GENOTYPE_PATH}/${i_rep_dir}/${i_ancestry}\"",
    "    else",
    "      ancestry_dir=\"${GENOTYPE_PATH}/${i_ancestry}\"",
    "    fi",
    "    if [ ! -d \"${ancestry_dir}\" ]; then",
    "      echo \"Missing ancestry directory: ${ancestry_dir}\" >&2",
    "      exit 1",
    "    fi",
    "    cd \"${ancestry_dir}\"",
    "    merge_list=.simgo_chr_file_list",
    "    rm -f \"${merge_list}\"",
    "",
    "    for chr in \"${CHROMOSOMES[@]}\"; do",
    "      input_prefix=\"${PREFIX}${i_ancestry}_chr${chr}.controls\"",
    "      output_prefix=\"${PREFIX}${i_ancestry}_chr${chr}\"",
    "      if [ ! -f \"${input_prefix}.gen\" ] || [ ! -f \"${input_prefix}.sample\" ]; then",
    "        echo \"Missing Oxford input: ${input_prefix}.gen or ${input_prefix}.sample\" >&2",
    "        exit 1",
    "      fi",
    "      echo \"PLINK QC rep=${i_rep_dir:-none} ancestry=${i_ancestry} chr=${chr}\"",
    "      \"${PLINK}\" --data \"${input_prefix}\" --oxford-single-chr \"${chr}\" \\",
    "        --make-bed --allow-no-sex \"${QC_ARGS[@]}\" --out \"${output_prefix}\"",
    append_to_merge_list,
    "    done",
    merge_block,
    "  done",
    "done",
    ""
  )

  writeLines(script, output_file)
  Sys.chmod(output_file, mode = "0755")
  invisible(output_file)
}

# Main entry point: generate scripts and optionally execute them.
simulate_1kg_hapgen2 <- function(reference_path,
                                 output_path,
                                 ancestries = c("EUR", "EAS"),
                                 chr_set = 1:22,
                                 rep_range = 1L,
                                 # Number of simulated controls per ancestry.
                                 sample_size = 10000L,
                                 # Keep zero for quantitative-trait simulations.
                                 n_cases = 0L,
                                 # Executable name, or an absolute HAPGEN2 path.
                                 hapgen2 = "hapgen2",
                                 genetic_map_path = file.path(reference_path, "genetic_map"),
                                 hap_prefix = "chr",
                                 legend_prefix = hap_prefix,
                                 map_prefix = "genetic_map_chr",
                                 map_suffix = "_combined_b37.txt",
                                 output_prefix = NULL,
                                 ne = hapgen2_ne_defaults(),
                                 no_haps_output = TRUE,
                                 output_file = "Hapgen2.sh",
                                 # Generate a PLINK conversion/QC script.
                                 qc = FALSE,
                                 qc_output_file = "hapgen2_plink_qc.sh",
                                 # Executable name, or an absolute PLINK 1.x path.
                                 plink = "plink",
                                 # Variant and sample QC thresholds.
                                 qc_maf = 0.01,
                                 qc_hwe = 1e-6,
                                 qc_geno = 0.05,
                                 qc_mind = 0.05,
                                 # Automatically merge only multi-chromosome runs.
                                 qc_merge = length(chr_set) > 1L,
                                 # Execute scripts now instead of only writing them.
                                 run = FALSE,
                                 # Validate all reference inputs before script creation.
                                 check_files = run) {
  qc <- .simgo_validate_flag(qc, "qc")
  run <- .simgo_validate_flag(run, "run")
  check_files <- .simgo_validate_flag(check_files, "check_files")

  if (run) {
    .simgo_check_executable(hapgen2, "hapgen2")
    if (qc) {
      .simgo_check_executable(plink, "plink")
    }
  }

  hapgen2_script <- write_1kg_hapgen2_script(
    reference_path = reference_path,
    output_path = output_path,
    ancestries = ancestries,
    chr_set = chr_set,
    rep_range = rep_range,
    sample_size = sample_size,
    n_cases = n_cases,
    hapgen2 = hapgen2,
    genetic_map_path = genetic_map_path,
    hap_prefix = hap_prefix,
    legend_prefix = legend_prefix,
    map_prefix = map_prefix,
    map_suffix = map_suffix,
    output_prefix = output_prefix,
    ne = ne,
    no_haps_output = no_haps_output,
    output_file = output_file,
    check_files = check_files
  )

  qc_script <- NULL
  if (qc) {
    qc_script <- write_hapgen2_plink_qc_script(
      genotype_path = output_path,
      ancestries = ancestries,
      chr_set = chr_set,
      rep_range = rep_range,
      prefix = output_prefix,
      plink = plink,
      maf = qc_maf,
      hwe = qc_hwe,
      geno = qc_geno,
      mind = qc_mind,
      merge = qc_merge,
      output_file = qc_output_file
    )
  }

  if (run) {
    .simgo_run_script(hapgen2_script, "HAPGEN2 simulation")
    if (!is.null(qc_script)) {
      .simgo_run_script(qc_script, "PLINK QC")
    }
  }

  invisible(list(
    hapgen2_script = hapgen2_script,
    qc_script = qc_script,
    executed = run
  ))
}

# Standalone QC entry point for previously generated HAPGEN2 genotypes.
run_hapgen2_plink_qc <- function(genotype_path,
                                 ancestries = c("EUR", "EAS"),
                                 chr_set = 1:22,
                                 rep_range = 1L,
                                 prefix = NULL,
                                 plink = "plink",
                                 maf = 0.01,
                                 hwe = 1e-6,
                                 geno = 0.05,
                                 mind = 0.05,
                                 merge = length(chr_set) > 1L,
                                 output_file = "hapgen2_plink_qc.sh",
                                 run = FALSE) {
  run <- .simgo_validate_flag(run, "run")
  if (run) {
    .simgo_check_executable(plink, "plink")
  }
  script_file <- write_hapgen2_plink_qc_script(
    genotype_path = genotype_path,
    ancestries = ancestries,
    chr_set = chr_set,
    rep_range = rep_range,
    prefix = prefix,
    plink = plink,
    maf = maf,
    hwe = hwe,
    geno = geno,
    mind = mind,
    merge = merge,
    output_file = output_file
  )
  if (run) {
    .simgo_run_script(script_file, "PLINK QC")
  }
  invisible(script_file)
}
