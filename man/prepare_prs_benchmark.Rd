\name{prepare_prs_benchmark}
\alias{prepare_prs_benchmark}
\alias{prepare_primary_prs_benchmark}
\alias{create_prs_benchmark_splits}
\alias{write_prs_benchmark_plink_script}
\title{Prepare Genotype Splits and GWAS Inputs for a PRS Benchmark}
\usage{
prepare_prs_benchmark(genotype_prefixes, phenotype_files, trait_names,
  sample_design, output_dir = "prs_benchmark", seed = 2026L,
  plink = "plink")

prepare_primary_prs_benchmark(genotype_prefixes, phenotype_files,
  trait_names, output_dir = "prs_benchmark", seed = 2026L,
  plink = "plink")
}
\arguments{
  \item{genotype_prefixes}{Named character vector of PLINK BED prefixes.}
  \item{phenotype_files}{Named character vector of phenotype files. Names must match the sample design.}
  \item{trait_names}{Phenotype column names to analyze.}
  \item{sample_design}{A data frame created by \code{prs_sample_design}.}
  \item{output_dir}{Directory for keep files, split PLINK files, summary statistics, and scripts.}
  \item{seed}{Integer random seed used for reproducible splitting.}
  \item{plink}{PLINK executable name or path.}
}
\value{Invisibly returns the design, sample manifest, and generated shell-script path.}
\description{
Creates reproducible, mutually exclusive training, tuning, and testing ID sets
for any number of populations. It also writes a shell script that creates PLINK
datasets and performs trait-specific GWAS in the training samples. External PRS
programs and LD resources are not installed or executed by this function.
}
\examples{
\dontrun{
design <- prs_sample_design(
  total = c(EUR = 60000, EAS = 30000),
  training = c(EUR = 50000, EAS = 20000),
  tuning = 5000,
  testing = 5000
)
prepare_prs_benchmark(
  genotype_prefixes = c(EUR = "geno/EUR", EAS = "geno/EAS"),
  phenotype_files = c(EUR = "phe/EUR.txt", EAS = "phe/EAS.txt"),
  trait_names = c("Phe_Trait1", "Phe_Trait2"),
  sample_design = design
)
}
}
