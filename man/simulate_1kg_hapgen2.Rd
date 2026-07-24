\name{simulate_1kg_hapgen2}
\alias{simulate_1kg_hapgen2}
\alias{run_hapgen2_plink_qc}
\alias{hapgen2_ne_defaults}
\title{Prepare HAPGEN2 Genotype Simulation and PLINK QC Workflows}
\usage{
hapgen2_ne_defaults()

simulate_1kg_hapgen2(reference_path, output_path,
  ancestries = c("EUR", "EAS"), chr_set = 1:22,
  rep_range = 1L, sample_size = 10000L, n_cases = 0L,
  hapgen2 = "hapgen2",
  genetic_map_path = file.path(reference_path, "genetic_map"),
  legend_path = NULL,
  hap_prefix = "chr", legend_prefix = hap_prefix,
  map_prefix = "genetic_map_chr",
  map_suffix = "_combined_b37.txt", output_prefix = NULL,
  ne = hapgen2_ne_defaults(), no_haps_output = TRUE,
  scripts_path = file.path(output_path, "scripts"), qc = FALSE,
  qc_output_path = file.path(output_path, "qc"),
  job_scheduler = "none", job_parameters = list(), plink = "plink",
  qc_maf = 0.01, qc_hwe = 1e-6, qc_geno = 0.05,
  qc_mind = 0.05, qc_merge = length(chr_set) > 1L,
  run = FALSE, check_files = run)

run_hapgen2_plink_qc(genotype_path,
  qc_output_path = file.path(genotype_path, "qc"),
  scripts_path = file.path(genotype_path, "scripts"),
  job_scheduler = "none", job_parameters = list(),
  ancestries = c("EUR", "EAS"), chr_set = 1:22,
  rep_range = 1L, prefix = NULL, plink = "plink",
  maf = 0.01, hwe = 1e-6, geno = 0.05, mind = 0.05,
  merge = length(chr_set) > 1L, run = FALSE)
}
\arguments{
  \item{reference_path}{Root directory containing ancestry subdirectories with HAPGEN2 haplotype and legend files.}
  \item{output_path,genotype_path}{Directory for simulated genotypes, or an existing HAPGEN2 output directory for standalone QC.}
  \item{qc_output_path}{Directory for QC-filtered PLINK files. By default, a \code{qc} subdirectory is created under the genotype output directory.}
  \item{scripts_path}{Directory for generated workflow scripts. By default, a \code{scripts} subdirectory is created under the genotype output directory.}
  \item{job_scheduler}{Scheduler directive format: \code{"none"}, \code{"slurm"}, \code{"sge"}, or \code{"pbs"}.}
  \item{job_parameters}{Named list of scheduler resources. Supported fields are \code{project}, \code{memory}, \code{cpus}, \code{walltime}, \code{queue}, \code{job_name}, and \code{log_path}.}
  \item{ancestries}{Unique ancestry labels matching subdirectories under the reference and output roots.}
  \item{chr_set}{Integer chromosomes to simulate or process.}
  \item{rep_range}{Positive replication identifiers. Use one genotype replication when phenotypes will be replicated separately.}
  \item{sample_size}{A scalar or ancestry-named vector giving the number of simulated controls.}
  \item{n_cases}{Number of simulated cases. Keep zero for quantitative-trait simulations.}
  \item{hapgen2,plink}{Executable name available on \code{PATH}, or an absolute executable path. PLINK 1.x syntax is used.}
  \item{genetic_map_path}{Directory containing chromosome-specific genetic maps.}
  \item{legend_path}{Optional shared directory containing chromosome legend files. When \code{NULL}, each legend is read from its ancestry directory.}
  \item{hap_prefix,legend_prefix,map_prefix,map_suffix}{Components used to construct reference filenames.}
  \item{output_prefix,prefix}{Optional prefix added to generated genotype filenames.}
  \item{ne}{Named effective population sizes.}
  \item{no_haps_output}{Suppress HAPGEN2 haplotype output when \code{TRUE}.}
  \item{qc}{Generate a PLINK conversion and QC script after preparing the HAPGEN2 script.}
  \item{qc_maf,maf}{Minimum minor-allele frequency retained by PLINK.}
  \item{qc_hwe,hwe}{Remove variants with Hardy-Weinberg test P-values below this threshold.}
  \item{qc_geno,geno}{Maximum missing-genotype proportion allowed per variant.}
  \item{qc_mind,mind}{Maximum missing-genotype proportion allowed per individual.}
  \item{qc_merge,merge}{Merge chromosome-level PLINK files within each ancestry and replication. The default is \code{FALSE} for one chromosome and \code{TRUE} for multiple chromosomes. Ancestries are never merged.}
  \item{run}{Execute generated scripts when \code{TRUE}; otherwise only write them.}
  \item{check_files}{Validate all reference inputs before writing the simulation script.}
}
\value{
\code{hapgen2_ne_defaults()} returns a named integer vector.
\code{simulate_1kg_hapgen2()} invisibly returns only the genotype output,
QC output, and script directories.
\code{run_hapgen2_plink_qc()} invisibly returns the QC script path.
}
\description{
Creates ancestry-, chromosome-, and replication-aware HAPGEN2 commands and an
optional PLINK conversion/QC workflow. These functions simulate and process
genotypes only; they do not generate phenotypes.
}
\examples{
\dontrun{
result <- simulate_1kg_hapgen2(
  reference_path = "/data/Hapmap3_1000G",
  output_path = "/data/Hapmap3_1000G/simulation/hapgen2",
  ancestries = c("EUR", "EAS"),
  chr_set = 1,
  rep_range = 1,
  sample_size = c(EUR = 60000, EAS = 30000),
  qc = TRUE,
  qc_merge = FALSE,
  run = FALSE
)
}
}
