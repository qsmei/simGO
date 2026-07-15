\name{simulate_1kg_hapgen2}
\alias{simulate_1kg_hapgen2}
\alias{write_1kg_hapgen2_script}
\alias{hapgen2_ne_defaults}
\alias{write_hapgen2_plink_qc_script}
\alias{run_hapgen2_plink_qc}
\title{Simulate Genotypes from 1000 Genomes Haplotype References}
\usage{
simulate_1kg_hapgen2(reference_path, output_path,
  ancestries = c("EUR", "EAS"), chr_set = 1:22, rep_range = 1,
  sample_size = 10000, n_cases = 0, hapgen2 = "hapgen2",
  genetic_map_path = file.path(reference_path, "genetic_map"),
  hap_prefix = "chr", legend_prefix = hap_prefix,
  map_prefix = "genetic_map_chr",
  map_suffix = "_combined_b37.txt", output_prefix = NULL,
  ne = hapgen2_ne_defaults(), no_haps_output = TRUE,
  output_file = "Hapgen2.sh", qc = FALSE,
  qc_output_file = "hapgen2_plink_qc.sh", plink = "plink",
  qc_maf = 0.01, qc_hwe = 1e-6, qc_geno = NULL,
  qc_mind = NULL, qc_merge = TRUE, run = FALSE,
  check_files = run)
}
\arguments{
  \item{reference_path}{Directory containing ancestry-specific HAPGEN2 haplotype and legend files.}
  \item{output_path}{Directory for simulated genotype files.}
  \item{ancestries}{Population labels matching subdirectories under \code{reference_path}.}
  \item{chr_set}{Chromosomes to simulate.}
  \item{rep_range}{Replication identifiers.}
  \item{sample_size}{A scalar or named ancestry-specific sample-size vector.}
  \item{n_cases}{Number of case samples; use zero for ordinary quantitative-trait genotype simulation.}
  \item{hapgen2}{HAPGEN2 executable name or path.}
  \item{genetic_map_path}{Directory containing genetic maps.}
  \item{hap_prefix,legend_prefix,map_prefix,map_suffix}{Reference-file naming components.}
  \item{output_prefix}{Optional prefix added to output filenames.}
  \item{ne}{Named effective population sizes.}
  \item{no_haps_output}{Suppress simulated haplotype output when \code{TRUE}.}
  \item{output_file}{Generated HAPGEN2 shell script.}
  \item{qc}{Generate a PLINK conversion and QC script.}
  \item{qc_output_file}{Generated PLINK script.}
  \item{plink}{PLINK executable name or path.}
  \item{qc_maf,qc_hwe,qc_geno,qc_mind}{Optional PLINK QC thresholds.}
  \item{qc_merge}{Merge chromosome-level PLINK files.}
  \item{run}{Execute generated scripts when \code{TRUE}.}
  \item{check_files}{Check all required reference files before execution.}
}
\value{Invisibly returns paths to the generated HAPGEN2 and optional PLINK scripts.}
\description{
Generates ancestry-, chromosome-, and replication-aware HAPGEN2 commands using
1000 Genomes-style haplotype references. The function simulates genotype data
only; it does not generate phenotypes.
}
\examples{
\dontrun{
simulate_1kg_hapgen2(
  reference_path = "/data/1kg_hapgen2",
  output_path = "/data/simulated_genotypes",
  ancestries = c("EUR", "EAS"),
  sample_size = c(EUR = 60000, EAS = 30000),
  qc = TRUE,
  run = FALSE
)
}
}
