\name{plotROC_from_ci}
\alias{plotROC_from_ci}
\title{Plot ROC curve with empirical confidence bands}
\description{
Plots an ROC curve from the empirical confidence band output of
\code{IncrementalAnalysis} for binomial outcomes.
}
\usage{
plotROC_from_ci(
  ci_df,
  title = "ROC",
  subtitle = NULL,
  ci_label = "95 percent CI",
  digits = 3,
  annotate_single = FALSE,
  annotate_pos = c(0.65, 0.15)
)
}
\arguments{
\item{ci_df}{Data frame containing \code{fpr}, \code{tpr_mean} and optional \code{tpr_lo}/\code{tpr_hi}.}
\item{title}{Plot title.}
\item{subtitle}{Plot subtitle.}
\item{ci_label}{Label used for confidence intervals.}
\item{digits}{Number of digits for AUC annotation.}
\item{annotate_single}{If \code{TRUE}, annotate a single curve directly on the plot.}
\item{annotate_pos}{x/y position for single-curve annotation.}
}
\value{
A ggplot object.
}
\examples{
\dontrun{
set.seed(1)
X <- matrix(rnorm(200 * 10), nrow = 200)
colnames(X) <- paste0("x", 1:10)
y <- rbinom(200, 1, 0.5)
vs <- sharp::VariableSelection(xdata = X, ydata = y, family = "binomial", K = 20)
inc <- IncrementalAnalysis(X, y, stability = vs, n_predictors = 5, k = 20, ncores = 1)
ci_best <- subset(inc$roc_ci, n_vars == which.max(inc$summary$metric_mean))
plotROC_from_ci(ci_best)
}
}
