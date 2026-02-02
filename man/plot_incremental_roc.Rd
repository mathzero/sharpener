\name{plot_incremental_roc}
\alias{plot_incremental_roc}
\title{Plot ROC curve for an incremental analysis}
\description{
Convenience wrapper that plots the ROC curve (with empirical bands) for a
chosen increment from \code{IncrementalAnalysis}.
}
\usage{
plot_incremental_roc(
  incremental_res,
  m = NULL,
  m_select = c("best_test", "best_train"),
  which = c("test", "train"),
  labels = NULL,
  label_predictors = FALSE,
  label_style = c("full", "added"),
  label_collapse = ", ",
  label_include_m = TRUE,
  title = NULL,
  subtitle = NULL,
  ...
)
}
\arguments{
\item{incremental_res}{Output from \code{IncrementalAnalysis}.}
\item{m}{Increment(s) to plot, or \code{"best_test"}/\code{"best_train"}.}
\item{m_select}{Selection rule used when \code{m} is \code{NULL}.}
\item{which}{Use test or train ROC bands.}
\item{labels}{Optional vector of legend labels; length 1 or length equal to \code{m}.}
\item{label_predictors}{If \code{TRUE}, include predictor names in legend labels.}
\item{label_style}{Use full predictor list or only the last added predictor for labels.}
\item{label_collapse}{Separator for predictor names when \code{label_style = "full"}.}
\item{label_include_m}{If \code{TRUE}, include \code{m} in auto-generated labels.}
\item{title}{Optional plot title.}
\item{subtitle}{Optional plot subtitle.}
\item{...}{Additional arguments passed to \code{plotROC_from_ci}.}
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
plot_incremental_roc(inc, m = "best_test")
}
}
