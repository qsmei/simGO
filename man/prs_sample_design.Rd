\name{prs_sample_design}
\alias{prs_sample_design}
\alias{simgo_primary_prs_design}
\title{Create a PRS Benchmark Sample Design}
\usage{
prs_sample_design(total, training = NULL, tuning = 5000L,
  testing = tuning, ancestry = names(total))

simgo_primary_prs_design()
}
\arguments{
  \item{total}{Named vector of total sample sizes, or a vector paired with \code{ancestry}.}
  \item{training}{Optional training sizes. If omitted, the remaining samples are assigned to training.}
  \item{tuning}{Tuning sizes; a scalar, ancestry-matched vector, or named vector.}
  \item{testing}{Testing sizes; by default equal to \code{tuning}.}
  \item{ancestry}{Ancestry or population labels.}
}
\value{A validated data frame containing total, training, tuning, and testing sizes.}
\description{
Constructs an ancestry-agnostic sample design for PRS simulation benchmarks.
The primary convenience preset returns EUR 50,000/5,000/5,000 and EAS
20,000/5,000/5,000 training/tuning/testing samples.
}
\examples{
prs_sample_design(
  total = c(EUR = 60000, EAS = 30000, AFR = 20000),
  tuning = 5000,
  testing = 5000
)
}
