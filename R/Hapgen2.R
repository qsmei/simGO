# HAPGEN2 / 1000 Genomes simulation helpers for simGO.
# Updated: 2026-05-26
#
# This file intentionally contains function definitions only. Long-running
# cluster examples are kept under inst/scripts/ so the package can be loaded
# without launching simulations or depending on user-specific paths.

hapgen2_ne_defaults <- function() {
  # Effective population sizes used in the original simGO HAPGEN2 workflow.
  # Users can override these values by passing a named vector to `ne`.
  c(EUR = 11418L, EAS = 14269L, AMR = 11418L, SAS = 14269L, AFR = 17469L)
}

validate_hapgen2_ancestries <- function(ancestries, ne) {
  # Every requested ancestry must have a corresponding Ne value because HAPGEN2
  # requires this parameter for each simulation run.
  missing_ne <- setdiff(ancestries, names(ne))
  if (length(missing_ne) > 0) {
    stop(
      "Missing effective population size (Ne) for ancestry: ",
      paste(missing_ne, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

write_1kg_hapgen2_script <- function(reference_path,
                                     output_path,
                                     ancestries = c("EUR", "EAS"),
                                     chr_set = 1:22,
                                     rep_range = 1,
                                     sample_size = 10000,
                                     n_cases = 0,
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
  # Write a standalone bash script for HAPGEN2 genotype simulation.
  # The script loops over replication, chromosome, and ancestry, producing files
  # under output_path/rep<rep>/<ancestry>/ to avoid overwriting replicates.
  if (missing(reference_path) || missing(output_path)) {
    stop("Both reference_path and output_path must be provided.", call. = FALSE)
  }
  if (length(ancestries) == 0 || length(chr_set) == 0 || length(rep_range) == 0) {
    stop("ancestries, chr_set, and rep_range must not be empty.", call. = FALSE)
  }
  if (length(sample_size) == 1) {
    sample_size <- stats::setNames(rep(as.integer(sample_size), length(ancestries)), ancestries)
  } else if (is.null(names(sample_size))) {
    if (length(sample_size) != length(ancestries)) {
      stop(
        "When sample_size is unnamed, its length must be 1 or match the number of ancestries.",
        call. = FALSE
      )
    }
    names(sample_size) <- ancestries
  }
  if (anyNA(names(sample_size)) || any(names(sample_size) == "")) {
    stop("All sample_size values must be named when a vector is supplied.", call. = FALSE)
  }
  missing_sample_size <- setdiff(ancestries, names(sample_size))
  if (length(missing_sample_size) > 0) {
    stop(
      "Missing sample_size for ancestry: ",
      paste(missing_sample_size, collapse = ", "),
      call. = FALSE
    )
  }
  validate_hapgen2_ancestries(ancestries, ne)

  if (check_files) {
    # When the user requests execution, fail early if required reference files
    # are not present, instead of failing midway through a long cluster job.
    required_dirs <- c(reference_path, genetic_map_path)
    missing_dirs <- required_dirs[!dir.exists(required_dirs)]
    if (length(missing_dirs) > 0) {
      stop("Missing directory: ", paste(missing_dirs, collapse = ", "), call. = FALSE)
    }

    required_files <- character()
    for (i_chr in chr_set) {
      required_files <- c(
        required_files,
        file.path(genetic_map_path, paste0(map_prefix, i_chr, map_suffix))
      )
      for (i_ancestry in ancestries) {
        required_files <- c(
          required_files,
          file.path(reference_path, i_ancestry, paste0(hap_prefix, i_chr, ".hap")),
          file.path(reference_path, i_ancestry, paste0(legend_prefix, i_chr, ".legend"))
        )
      }
    }
    missing_files <- required_files[!file.exists(required_files)]
    if (length(missing_files) > 0) {
      stop(
        "Missing HAPGEN2 input file(s):\n",
        paste(utils::head(missing_files, 20), collapse = "\n"),
        if (length(missing_files) > 20) "\n..." else "",
        call. = FALSE
      )
    }
  }

  prefix_value <- if (is.null(output_prefix) || identical(output_prefix, "")) "" else paste0(output_prefix, "_")
  # HAPGEN2 can optionally suppress haplotype output; for simGO we generally
  # need only the simulated Oxford genotype files for PLINK conversion/QC.
  output_line <- if (isTRUE(no_haps_output)) {
    "        -o \"${out_prefix}\" -no_haps_output"
  } else {
    "        -o \"${out_prefix}\""
  }

  sample_size_cases <- paste(
    # Build bash case statements so each ancestry can have its own sample size.
    sprintf("        %s) sample_n=%s ;;", names(sample_size), as.integer(sample_size)),
    collapse = "\n"
  )
  ne_cases <- paste(
    # Build the matching ancestry-specific Ne case statement.
    sprintf("        %s) Ne=%s ;;", names(ne), as.integer(ne)),
    collapse = "\n"
  )

  script <- c(
    "#!/usr/bin/env bash",
    "set -euo pipefail",
    "",
    paste0("REFERENCE_PATH=", shQuote(path.expand(reference_path))),
    paste0("OUTPUT_PATH=", shQuote(path.expand(output_path))),
    paste0("GENETIC_MAP_PATH=", shQuote(path.expand(genetic_map_path))),
    paste0("HAPGEN2=", shQuote(hapgen2)),
    paste0("HAP_PREFIX=", shQuote(hap_prefix)),
    paste0("LEGEND_PREFIX=", shQuote(legend_prefix)),
    paste0("MAP_PREFIX=", shQuote(map_prefix)),
    paste0("MAP_SUFFIX=", shQuote(map_suffix)),
    paste0("OUTPUT_PREFIX=", shQuote(prefix_value)),
    paste0("N_CASES=", as.integer(n_cases)),
    paste0("ANCESTRIES=(", paste(shQuote(ancestries), collapse = " "), ")"),
    paste0("CHROMOSOMES=(", paste(as.integer(chr_set), collapse = " "), ")"),
    paste0("REPLICATIONS=(", paste(as.integer(rep_range), collapse = " "), ")"),
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
  invisible(normalizePath(output_file, mustWork = FALSE))
}

simulate_1kg_hapgen2 <- function(reference_path,
                                 output_path,
                                 ancestries = c("EUR", "EAS"),
                                 chr_set = 1:22,
                                 rep_range = 1,
                                 sample_size = 10000,
                                 n_cases = 0,
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
                                 qc = FALSE,
                                 qc_output_file = "hapgen2_plink_qc.sh",
                                 plink = "plink",
                                 qc_maf = 0.01,
                                 qc_hwe = 1e-6,
                                 qc_geno = NULL,
                                 qc_mind = NULL,
                                 qc_merge = TRUE,
                                 run = FALSE,
                                 check_files = run) {
  # High-level workflow function:
  # 1. write a HAPGEN2 simulation script;
  # 2. optionally write a PLINK QC script for the simulated genotypes;
  # 3. optionally run both scripts in sequence.
  script_file <- write_1kg_hapgen2_script(
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

  qc_script_file <- NULL
  if (isTRUE(qc)) {
    # The QC script is replication-aware because HAPGEN2 output is organized as
    # output_path/rep<rep>/<ancestry>/.
    qc_script_file <- write_hapgen2_plink_qc_script(
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

  if (isTRUE(run)) {
    # Keep execution simple and transparent: run HAPGEN2 first, then QC only if
    # the simulation succeeded and qc = TRUE.
    status <- system2("bash", script_file)
    if (!identical(status, 0L)) {
      stop("HAPGEN2 script failed with status ", status, call. = FALSE)
    }
    if (!is.null(qc_script_file)) {
      qc_status <- system2("bash", qc_script_file)
      if (!identical(qc_status, 0L)) {
        stop("PLINK QC script failed with status ", qc_status, call. = FALSE)
      }
    }
  }

  invisible(list(hapgen2_script = script_file, qc_script = qc_script_file))
}

generate_hapgen2_sh <- function(path,
                                ancestries,
                                rep_range,
                                chr_set,
                                EUR_Ne = 11418,
                                EAS_Ne = 14269,
                                AMR_Ne = 11418,
                                SAS_Ne = 14269,
                                AFR_Ne = 17469,
                                prefix = NULL,
                                sample_size = 12000,
                                output_path = NULL,
                                run = FALSE,
                                hap_prefix = "chr",
                                output_file = "Hapgen2.sh",
                                hapgen2 = "hapgen2",
                                genetic_map_path = file.path(path, "genetic_map")) {
  ne <- c(EUR = EUR_Ne, EAS = EAS_Ne, AMR = AMR_Ne, SAS = SAS_Ne, AFR = AFR_Ne)
  simulate_1kg_hapgen2(
    reference_path = path,
    output_path = output_path,
    ancestries = ancestries,
    chr_set = chr_set,
    rep_range = rep_range,
    sample_size = sample_size,
    n_cases = 0,
    hapgen2 = hapgen2,
    genetic_map_path = genetic_map_path,
    hap_prefix = hap_prefix,
    legend_prefix = hap_prefix,
    output_prefix = prefix,
    ne = ne,
    no_haps_output = TRUE,
    output_file = output_file,
    run = run
  )
}

write_hapgen2_plink_qc_script <- function(genotype_path,
                                          ancestries = c("EUR", "EAS"),
                                          chr_set = 1:22,
                                          rep_range = NULL,
                                          prefix = NULL,
                                          plink = "plink",
                                          maf = 0.01,
                                          hwe = 1e-6,
                                          geno = NULL,
                                          mind = NULL,
                                          merge = TRUE,
                                          output_file = "hapgen2_plink_qc.sh") {
  # Write a PLINK QC script for HAPGEN2 Oxford-format outputs.
  # Expected input naming:
  #   <prefix><ancestry>_chr<chr>.controls.gen/sample
  # Expected output naming:
  #   <prefix><ancestry>_chr<chr>.bed/bim/fam
  if (missing(genotype_path)) {
    stop("genotype_path must be provided.", call. = FALSE)
  }
  if (length(ancestries) == 0 || length(chr_set) == 0) {
    stop("ancestries and chr_set must not be empty.", call. = FALSE)
  }

  prefix_value <- if (is.null(prefix) || identical(prefix, "")) "" else paste0(prefix, "_")
  # If rep_range is supplied, QC each output_path/rep<rep>/<ancestry>/ folder.
  # If rep_range is NULL, QC the older flat output_path/<ancestry>/ layout.
  rep_dirs <- if (is.null(rep_range)) {
    "REPLICATION_DIRS=(\"\")"
  } else {
    paste0("REPLICATION_DIRS=(", paste(shQuote(paste0("rep", as.integer(rep_range))), collapse = " "), ")")
  }
  qc_args <- c(
    # Only include QC thresholds that the user explicitly keeps non-NULL.
    if (!is.null(maf)) c("--maf", as.character(maf)),
    if (!is.null(hwe)) c("--hwe", as.character(hwe)),
    if (!is.null(geno)) c("--geno", as.character(geno)),
    if (!is.null(mind)) c("--mind", as.character(mind))
  )
  qc_arg_line <- paste(shQuote(qc_args), collapse = " ")
  merge_block <- if (isTRUE(merge)) {
    # Merge per-chromosome PLINK files within each ancestry/replicate folder.
    c(
      "    if [ -s file_list ]; then",
      "      \"${PLINK}\" --merge-list file_list --make-bed --out \"${PREFIX}${i_ancestry}\" --allow-no-sex",
      "      mv file_list chr_file_list",
      "    fi"
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
    paste0("CHROMOSOMES=(", paste(as.integer(chr_set), collapse = " "), ")"),
    "",
    "for i_rep_dir in \"${REPLICATION_DIRS[@]}\"; do",
    "  for i_ancestry in \"${ANCESTRIES[@]}\"; do",
    "    if [ -n \"${i_rep_dir}\" ]; then",
    "      ancestry_dir=\"${GENOTYPE_PATH}/${i_rep_dir}/${i_ancestry}\"",
    "    else",
    "      ancestry_dir=\"${GENOTYPE_PATH}/${i_ancestry}\"",
    "    fi",
    "    cd \"${ancestry_dir}\"",
    "    rm -f file_list",
    "    for chr in \"${CHROMOSOMES[@]}\"; do",
    "      input_prefix=\"${PREFIX}${i_ancestry}_chr${chr}.controls\"",
    "      output_prefix=\"${PREFIX}${i_ancestry}_chr${chr}\"",
    "      gen_file=\"${input_prefix}.gen\"",
    "      if [ ! -f \"${gen_file}\" ]; then",
    "        echo \"Missing ${gen_file}; skipping\" >&2",
    "        continue",
    "      fi",
    "      echo \"PLINK QC rep=${i_rep_dir:-none} ancestry=${i_ancestry} chr=${chr}\"",
    "      sed -E \"s/^snp_[0-9]+ /${chr} /\" \"${gen_file}\" > \"${gen_file}.tmp\"",
    "      mv \"${gen_file}.tmp\" \"${gen_file}\"",
    "      \"${PLINK}\" --data \"${input_prefix}\" --oxford-single-chr \"${chr}\" --make-bed \"${QC_ARGS[@]}\" --out \"${output_prefix}\"",
    "      echo \"${output_prefix}\" >> file_list",
    "    done",
    merge_block,
    "  done",
    "done",
    ""
  )

  writeLines(script, output_file)
  invisible(normalizePath(output_file, mustWork = FALSE))
}

run_hapgen2_plink_qc <- function(genotype_path,
                                 ancestries = c("EUR", "EAS"),
                                 chr_set = 1:22,
                                 rep_range = NULL,
                                 prefix = NULL,
                                 plink = "plink",
                                 maf = 0.01,
                                 hwe = 1e-6,
                                 geno = NULL,
                                 mind = NULL,
                                 merge = TRUE,
                                 output_file = "hapgen2_plink_qc.sh",
                                 run = FALSE) {
  # Convenience wrapper: create the QC script and optionally run it.
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

  if (isTRUE(run)) {
    status <- system2("bash", script_file)
    if (!identical(status, 0L)) {
      stop("PLINK QC script failed with status ", status, call. = FALSE)
    }
  }

  invisible(script_file)
}

generate_qc_script <- function(path,
                               ancestries,
                               chromosomes,
                               rep_range = NULL,
                               prefix = NULL,
                               plink = "plink",
                               run = FALSE,
                               output_file = "script.sh",
                               maf = 0.01,
                               hwe = 1e-6,
                               geno = NULL,
                               mind = NULL,
                               merge = TRUE) {
  # Backward-compatible wrapper for older scripts that called
  # `generate_qc_script(path, ancestries, chromosomes, ...)`.
  run_hapgen2_plink_qc(
    genotype_path = path,
    ancestries = ancestries,
    chr_set = chromosomes,
    rep_range = rep_range,
    prefix = prefix,
    plink = plink,
    maf = maf,
    hwe = hwe,
    geno = geno,
    mind = mind,
    merge = merge,
    output_file = output_file,
    run = run
  )
}
