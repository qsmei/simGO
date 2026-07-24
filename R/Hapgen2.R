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

.simgo_create_directory <- function(path, name) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || path == "") {
    stop(name, " must be one non-empty directory path.", call. = FALSE)
  }
  path <- path.expand(path)
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(path)) {
    stop(
      "Unable to create ", name, ": ", path,
      ". Check the path and write permissions.",
      call. = FALSE
    )
  }
  normalizePath(path, mustWork = TRUE)
}

.simgo_job_value <- function(x, name) {
  if (is.null(x)) {
    return(NULL)
  }
  if (length(x) != 1L || is.na(x) || identical(as.character(x), "") ||
      grepl("[\r\n]", as.character(x))) {
    stop("job_parameters$", name, " must be one non-empty value.", call. = FALSE)
  }
  as.character(x)
}

.simgo_scheduler_header <- function(job_scheduler,
                                    job_parameters,
                                    job_type,
                                    output_path) {
  choices <- c("none", "slurm", "sge", "pbs", "custom")
  if (!is.character(job_scheduler) || length(job_scheduler) != 1L ||
      is.na(job_scheduler) || !(job_scheduler %in% choices)) {
    stop("job_scheduler must be one of: ", paste(choices, collapse = ", "),
         ".", call. = FALSE)
  }
  # Character input is always treated as a complete scheduler header. This
  # allows cluster-specific directives without requiring a special mode.
  if (is.character(job_parameters)) {
    if (length(job_parameters) == 0L || anyNA(job_parameters)) {
      stop("Character job_parameters must contain scheduler header text.",
           call. = FALSE)
    }
    header <- unlist(strsplit(job_parameters, "\r?\n"), use.names = FALSE)
    header <- header[nzchar(trimws(header))]
    # The generated script already contains one portable Bash shebang.
    header <- header[!grepl("^#!", trimws(header))]
    if (length(header) > 0L && any(!grepl("^#", trimws(header)))) {
      stop("Custom job header lines must begin with '#'.", call. = FALSE)
    }
    trimmed_header <- trimws(header)
    if (job_scheduler == "sge" && any(grepl("^#PBS", trimmed_header))) {
      stop("SGE was selected, but job_parameters contains PBS directives (#PBS).",
           call. = FALSE)
    }
    if (job_scheduler == "pbs" && any(grepl("^#\\$", trimmed_header))) {
      stop("PBS was selected, but job_parameters contains SGE directives (#$). Use job_scheduler = 'sge'.",
           call. = FALSE)
    }
    if (job_scheduler == "slurm" &&
        any(grepl("^#(PBS|\\$)", trimmed_header))) {
      stop("Slurm was selected, but job_parameters contains non-Slurm directives.",
           call. = FALSE)
    }
    # If the user does not specify log destinations, keep scheduler logs with
    # the analysis output instead of the R submission working directory.
    if (job_scheduler == "sge") {
      if (!any(grepl("^#\\$[[:space:]]+-o([[:space:]]|$)", trimmed_header))) {
        header <- c(header, paste0("#$ -o ", output_path))
      }
      merges_output <- any(grepl("^#\\$[[:space:]]+-j[[:space:]]+y", trimmed_header))
      if (!merges_output &&
          !any(grepl("^#\\$[[:space:]]+-e([[:space:]]|$)", trimmed_header))) {
        header <- c(header, paste0("#$ -e ", output_path))
      }
    } else if (job_scheduler == "slurm") {
      if (!any(grepl("^#SBATCH.*(--output|-o)(=|[[:space:]])", trimmed_header))) {
        header <- c(header, paste0("#SBATCH --output=", output_path, "/%x_%j.log"))
      }
      if (!any(grepl("^#SBATCH.*(--error|-e)(=|[[:space:]])", trimmed_header))) {
        header <- c(header, paste0("#SBATCH --error=", output_path, "/%x_%j.err"))
      }
    } else if (job_scheduler == "pbs") {
      if (!any(grepl("^#PBS[[:space:]]+-o([[:space:]]|$)", trimmed_header))) {
        header <- c(header, paste0("#PBS -o ", output_path))
      }
      merges_output <- any(grepl("^#PBS[[:space:]]+-j[[:space:]]+oe", trimmed_header))
      if (!merges_output &&
          !any(grepl("^#PBS[[:space:]]+-e([[:space:]]|$)", trimmed_header))) {
        header <- c(header, paste0("#PBS -e ", output_path))
      }
    }
    return(gsub("\\{job_type\\}", job_type, header))
  }

  if (job_scheduler == "custom") {
    stop("job_scheduler = 'custom' requires character job_parameters.",
         call. = FALSE)
  }

  if (!is.list(job_parameters)) {
    stop("job_parameters must be a named list or scheduler header text.", call. = FALSE)
  }
  if (length(job_parameters) > 0L &&
      (is.null(names(job_parameters)) || any(names(job_parameters) == ""))) {
    stop("job_parameters must be a named list.", call. = FALSE)
  }
  allowed <- c(
    "project", "memory", "cpus", "walltime", "queue", "job_name", "log_path"
  )
  unknown <- setdiff(names(job_parameters), allowed)
  if (length(unknown) > 0L) {
    stop("Unknown job_parameters: ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  if (job_scheduler == "none") {
    if (length(job_parameters) > 0L) {
      stop("job_parameters requires job_scheduler = 'slurm', 'sge', or 'pbs'.",
           call. = FALSE)
    }
    return(character())
  }

  project <- .simgo_job_value(job_parameters$project, "project")
  memory <- .simgo_job_value(job_parameters$memory, "memory")
  walltime <- .simgo_job_value(job_parameters$walltime, "walltime")
  queue <- .simgo_job_value(job_parameters$queue, "queue")
  base_name <- .simgo_job_value(job_parameters$job_name, "job_name")
  if (is.null(base_name)) {
    base_name <- "simGO"
  }
  job_name <- paste0(base_name, "_", job_type)

  cpus <- job_parameters$cpus
  if (!is.null(cpus)) {
    cpus <- .simgo_validate_integer_vector(cpus, "job_parameters$cpus", 1L)[1L]
  }

  log_path <- job_parameters$log_path
  if (is.null(log_path)) {
    log_path <- output_path
  } else {
    if (!grepl("^(/|~|[A-Za-z]:[/\\\\])", log_path)) {
      log_path <- file.path(output_path, log_path)
    }
  }
  log_path <- .simgo_create_directory(log_path, "job_parameters$log_path")

  if (job_scheduler == "slurm") {
    return(c(
      paste0("#SBATCH --job-name=", job_name),
      if (!is.null(project)) paste0("#SBATCH --account=", project),
      if (!is.null(memory)) paste0("#SBATCH --mem=", memory),
      if (!is.null(cpus)) paste0("#SBATCH --cpus-per-task=", cpus),
      if (!is.null(walltime)) paste0("#SBATCH --time=", walltime),
      if (!is.null(queue)) paste0("#SBATCH --partition=", queue),
      if (!is.null(log_path)) paste0("#SBATCH --output=", log_path, "/%x_%j.out"),
      if (!is.null(log_path)) paste0("#SBATCH --error=", log_path, "/%x_%j.err")
    ))
  }

  if (job_scheduler == "sge") {
    return(c(
      paste0("#$ -N ", job_name),
      if (!is.null(project)) paste0("#$ -P ", project),
      if (!is.null(memory)) paste0("#$ -l mem_per_core=", memory),
      if (!is.null(cpus)) paste0("#$ -pe omp ", cpus),
      if (!is.null(walltime)) paste0("#$ -l h_rt=", walltime),
      if (!is.null(queue)) paste0("#$ -q ", queue),
      "#$ -cwd",
      "#$ -j y",
      if (!is.null(log_path)) paste0("#$ -o ", log_path)
    ))
  }

  select_parts <- c(
    "select=1",
    if (!is.null(cpus)) paste0("ncpus=", cpus),
    if (!is.null(memory)) paste0("mem=", memory)
  )
  c(
    paste0("#PBS -N ", job_name),
    if (!is.null(project)) paste0("#PBS -A ", project),
    if (length(select_parts) > 1L) paste0("#PBS -l ", paste(select_parts, collapse = ":")),
    if (!is.null(walltime)) paste0("#PBS -l walltime=", walltime),
    if (!is.null(queue)) paste0("#PBS -q ", queue),
    "#PBS -j oe",
    if (!is.null(log_path)) paste0("#PBS -o ", log_path)
  )
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

.simgo_submit_script <- function(script_file,
                                 job_scheduler,
                                 label,
                                 dependency = NULL) {
  if (job_scheduler == "slurm") {
    command <- "sbatch"
    args <- c(
      "--parsable",
      if (!is.null(dependency)) paste0("--dependency=afterok:", dependency),
      shQuote(script_file)
    )
  } else if (job_scheduler == "sge") {
    command <- "qsub"
    args <- c(
      "-terse",
      if (!is.null(dependency)) c("-hold_jid", dependency),
      shQuote(script_file)
    )
  } else if (job_scheduler == "pbs") {
    command <- "qsub"
    args <- c(
      if (!is.null(dependency)) c("-W", paste0("depend=afterok:", dependency)),
      shQuote(script_file)
    )
  } else {
    stop("Automatic submission requires job_scheduler = 'slurm', 'sge', or 'pbs'.",
         call. = FALSE)
  }

  .simgo_check_executable(command, command)
  output <- system2(command, args = args, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }
  if (!identical(as.integer(status), 0L)) {
    stop(label, " submission failed:\n", paste(output, collapse = "\n"),
         call. = FALSE)
  }
  output <- trimws(output[nzchar(trimws(output))])
  if (length(output) == 0L) {
    stop(label, " submission did not return a job ID.", call. = FALSE)
  }
  job_id <- output[length(output)]
  if (job_scheduler == "slurm") {
    job_id <- sub(";.*$", "", job_id)
  } else if (job_scheduler == "sge") {
    # Normally qsub -terse returns only the ID, but support verbose SGE output.
    if (grepl("Your job[[:space:]]+[0-9]+", job_id)) {
      job_id <- sub(".*Your job[[:space:]]+([0-9]+).*", "\\1", job_id)
    } else {
      job_id <- sub("[[:space:]].*$", "", job_id)
    }
  } else if (job_scheduler == "pbs") {
    # PBS usually returns 12345.server; retain that complete identifier.
    job_id <- sub("[[:space:]].*$", "", job_id)
  }
  if (!grepl("^[A-Za-z0-9_.:-]+$", job_id)) {
    stop(label, " returned an invalid job ID: ", job_id, call. = FALSE)
  }
  message(label, " submitted as job ", job_id, ".")
  job_id
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
                                     legend_path = NULL,
                                     hap_prefix = "chr",
                                     legend_prefix = hap_prefix,
                                     map_prefix = "genetic_map_chr",
                                     map_suffix = "_combined_b37.txt",
                                     output_prefix = NULL,
                                     ne = hapgen2_ne_defaults(),
                                     no_haps_output = TRUE,
                                     scheduler_header = character(),
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
  if (!is.null(legend_path)) {
    if (!is.character(legend_path) || length(legend_path) != 1L ||
        is.na(legend_path) || legend_path == "") {
      stop("legend_path must be NULL or one non-empty directory path.", call. = FALSE)
    }
    legend_path <- path.expand(legend_path)
  }
  output_file <- .simgo_prepare_script_path(output_file)

  if (check_files) {
    required_dirs <- c(reference_path, genetic_map_path, legend_path)
    missing_dirs <- required_dirs[!dir.exists(required_dirs)]
    if (length(missing_dirs) > 0L) {
      stop("Missing directory: ", paste(missing_dirs, collapse = ", "), call. = FALSE)
    }

    required_files <- unlist(lapply(chr_set, function(chr) {
      hap_files <- file.path(
        reference_path,
        ancestries,
        paste0(hap_prefix, chr, ".hap")
      )
      legend_files <- if (is.null(legend_path)) {
        file.path(
          reference_path,
          ancestries,
          paste0(legend_prefix, chr, ".legend")
        )
      } else {
        file.path(legend_path, paste0(legend_prefix, chr, ".legend"))
      }
      c(
        file.path(genetic_map_path, paste0(map_prefix, chr, map_suffix)),
        hap_files,
        legend_files
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
    scheduler_header,
    "set -euo pipefail",
    "",
    paste0("REFERENCE_PATH=", shQuote(reference_path)),
    paste0("OUTPUT_PATH=", shQuote(output_path)),
    paste0("GENETIC_MAP_PATH=", shQuote(genetic_map_path)),
    paste0("LEGEND_PATH=", shQuote(if (is.null(legend_path)) "" else legend_path)),
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
    "      if [ -n \"${LEGEND_PATH}\" ]; then",
    "        legend_file=\"${LEGEND_PATH}/${LEGEND_PREFIX}${i_chr}.legend\"",
    "      else",
    "        legend_file=\"${REFERENCE_PATH}/${i_ancestry}/${LEGEND_PREFIX}${i_chr}.legend\"",
    "      fi",
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
                                          qc_output_path = file.path(genotype_path, "qc"),
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
                                          scheduler_header = character(),
                                          output_file = "hapgen2_plink_qc.sh") {
  if (missing(genotype_path)) {
    stop("genotype_path must be provided.", call. = FALSE)
  }
  if (!is.character(qc_output_path) || length(qc_output_path) != 1L ||
      is.na(qc_output_path) || qc_output_path == "") {
    stop("qc_output_path must be one non-empty directory path.", call. = FALSE)
  }
  genotype_path <- path.expand(genotype_path)
  qc_output_path <- .simgo_create_directory(qc_output_path, "qc_output_path")
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
      "      echo \"No chromosome datasets were available to merge in ${qc_dir}\" >&2",
      "      exit 1",
      "    fi",
      "    \"${PLINK}\" --merge-list \"${merge_list}\" --make-bed --allow-no-sex --out \"${qc_dir}/${PREFIX}${i_ancestry}\"",
      "    mv \"${merge_list}\" \"${qc_dir}/chr_file_list\""
    )
  } else {
    character()
  }

  script <- c(
    "#!/usr/bin/env bash",
    scheduler_header,
    "set -euo pipefail",
    "",
    paste0("GENOTYPE_PATH=", shQuote(genotype_path)),
    paste0("QC_OUTPUT_PATH=", shQuote(qc_output_path)),
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
    "      qc_dir=\"${QC_OUTPUT_PATH}/${i_rep_dir}/${i_ancestry}\"",
    "    else",
    "      ancestry_dir=\"${GENOTYPE_PATH}/${i_ancestry}\"",
    "      qc_dir=\"${QC_OUTPUT_PATH}/${i_ancestry}\"",
    "    fi",
    "    if [ ! -d \"${ancestry_dir}\" ]; then",
    "      echo \"Missing ancestry directory: ${ancestry_dir}\" >&2",
    "      exit 1",
    "    fi",
    "    mkdir -p \"${qc_dir}\"",
    "    merge_list=\"${qc_dir}/.simgo_chr_file_list\"",
    "    rm -f \"${merge_list}\"",
    "",
    "    for chr in \"${CHROMOSOMES[@]}\"; do",
    "      input_prefix=\"${ancestry_dir}/${PREFIX}${i_ancestry}_chr${chr}.controls\"",
    "      output_prefix=\"${qc_dir}/${PREFIX}${i_ancestry}_chr${chr}\"",
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
                                 # Optional shared directory containing chrN.legend files.
                                 legend_path = NULL,
                                 hap_prefix = "chr",
                                 legend_prefix = hap_prefix,
                                 map_prefix = "genetic_map_chr",
                                 map_suffix = "_combined_b37.txt",
                                 output_prefix = NULL,
                                 ne = hapgen2_ne_defaults(),
                                 no_haps_output = TRUE,
                                 # Directory for generated workflow scripts.
                                 scripts_path = file.path(output_path, "scripts"),
                                 # Generate a PLINK conversion/QC script.
                                 qc = FALSE,
                                 # Directory for QC-filtered PLINK files.
                                 qc_output_path = file.path(output_path, "qc"),
                                 # Cluster scheduler used for script directives.
                                 job_scheduler = "none",
                                 # Scheduler resources such as project, memory and CPUs.
                                 job_parameters = list(),
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

  output_path <- .simgo_create_directory(output_path, "output_path")

  hapgen2_scheduler_header <- .simgo_scheduler_header(
    job_scheduler = job_scheduler,
    job_parameters = job_parameters,
    job_type = "hapgen2",
    output_path = output_path
  )

  scripts_path <- .simgo_create_directory(scripts_path, "scripts_path")
  hapgen2_script_file <- file.path(scripts_path, "hapgen2_simGO.sh")

  if (qc) {
    qc_output_path <- .simgo_create_directory(qc_output_path, "qc_output_path")
    qc_script_file <- file.path(scripts_path, "hapgen2_simGO_QC.sh")
  }

  if (run) {
    if (job_scheduler == "none") {
      .simgo_check_executable(hapgen2, "hapgen2")
      if (qc) {
        .simgo_check_executable(plink, "plink")
      }
    } else if (job_scheduler == "custom") {
      stop("run = TRUE cannot infer a submission command for job_scheduler = 'custom'.",
           call. = FALSE)
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
    legend_path = legend_path,
    hap_prefix = hap_prefix,
    legend_prefix = legend_prefix,
    map_prefix = map_prefix,
    map_suffix = map_suffix,
    output_prefix = output_prefix,
    ne = ne,
    no_haps_output = no_haps_output,
    scheduler_header = hapgen2_scheduler_header,
    output_file = hapgen2_script_file,
    check_files = check_files
  )

  qc_script <- NULL
  if (qc) {
    qc_scheduler_header <- .simgo_scheduler_header(
      job_scheduler = job_scheduler,
      job_parameters = job_parameters,
      job_type = "QC",
      output_path = output_path
    )
    qc_script <- write_hapgen2_plink_qc_script(
      genotype_path = output_path,
      qc_output_path = qc_output_path,
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
      scheduler_header = qc_scheduler_header,
      output_file = qc_script_file
    )
  }

  if (run) {
    if (job_scheduler == "none") {
      .simgo_run_script(hapgen2_script, "HAPGEN2 simulation")
      if (!is.null(qc_script)) {
        .simgo_run_script(qc_script, "PLINK QC")
      }
    } else {
      hapgen2_job <- .simgo_submit_script(
        hapgen2_script, job_scheduler, "HAPGEN2 simulation"
      )
      if (!is.null(qc_script)) {
        .simgo_submit_script(
          qc_script,
          job_scheduler,
          "PLINK QC",
          dependency = hapgen2_job
        )
      }
    }
  }

  invisible(list(
    output_path = normalizePath(output_path, mustWork = FALSE),
    qc_output_path = if (qc) normalizePath(qc_output_path, mustWork = FALSE) else NULL,
    scripts_path = normalizePath(scripts_path, mustWork = FALSE)
  ))
}

# Standalone QC entry point for previously generated HAPGEN2 genotypes.
run_hapgen2_plink_qc <- function(genotype_path,
                                 qc_output_path = file.path(genotype_path, "qc"),
                                 scripts_path = file.path(genotype_path, "scripts"),
                                 job_scheduler = "none",
                                 job_parameters = list(),
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
                                 run = FALSE) {
  run <- .simgo_validate_flag(run, "run")
  scripts_path <- .simgo_create_directory(scripts_path, "scripts_path")
  output_file <- file.path(scripts_path, "hapgen2_simGO_QC.sh")
  scheduler_header <- .simgo_scheduler_header(
    job_scheduler = job_scheduler,
    job_parameters = job_parameters,
    job_type = "QC",
    output_path = genotype_path
  )
  if (run) {
    if (job_scheduler == "none") {
      .simgo_check_executable(plink, "plink")
    } else if (job_scheduler == "custom") {
      stop("run = TRUE cannot infer a submission command for job_scheduler = 'custom'.",
           call. = FALSE)
    }
  }
  script_file <- write_hapgen2_plink_qc_script(
    genotype_path = genotype_path,
    qc_output_path = qc_output_path,
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
    scheduler_header = scheduler_header,
    output_file = output_file
  )
  if (run) {
    if (job_scheduler == "none") {
      .simgo_run_script(script_file, "PLINK QC")
    } else {
      .simgo_submit_script(script_file, job_scheduler, "PLINK QC")
    }
  }
  invisible(script_file)
}
