\name{IncrementalAnalysis}
\alias{IncrementalAnalysis}
\title{Incremental model evaluation based on sharp variable selection}
\description{
Evaluates model performance as variables are added in order of importance,
optionally using stability information from \code{sharp::VariableSelection}.
}
\usage{
IncrementalAnalysis(
  xdata,
  ydata,
  n_predictors = NULL,
  stability = NULL,
  test_ratio = 0.3,
  k = 100,
  seed = 123,
  ncores = 4,
  verbose = TRUE,
  ci_level = 0.95,
  roc_points = 101,
  max_predictors = 40
)
}
\arguments{
\item{xdata}{Predictor matrix or data frame.}
\item{ydata}{Outcome vector or a \code{survival::Surv} object.}
\item{n_predictors}{Number of predictors to evaluate. Defaults to all, capped by \code{max_predictors}.}
\item{stability}{Optional output from \code{sharp::VariableSelection} to define predictor order.}
\item{test_ratio}{Proportion of data reserved for testing in each split.}
\item{k}{Number of resampling splits.}
\item{seed}{Random seed for reproducibility.}
\item{ncores}{Number of parallel workers for \code{future}.}
\item{verbose}{Logical; print progress.}
\item{ci_level}{Confidence level for empirical intervals.}
\item{roc_points}{Number of ROC grid points for binomial outcomes.}
\item{max_predictors}{Maximum number of predictors to consider.}
}
\value{
A list with a summary data frame, per-model details, and optional ROC curve
confidence bands for binomial outcomes.
}
\examples{
\dontrun{
set.seed(1)
X <- matrix(rnorm(200 * 10), nrow = 200)
colnames(X) <- paste0("x", 1:10)
y <- rbinom(200, 1, 0.5)
vs <- sharp::VariableSelection(xdata = X, ydata = y, family = "binomial", K = 20)
res <- IncrementalAnalysis(X, y, stability = vs, n_predictors = 5, k = 20, ncores = 1)
}
}
