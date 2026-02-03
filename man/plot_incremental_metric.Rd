\name{plot_incremental_metric}
\alias{plot_incremental_metric}
\title{Plot incremental performance metrics}
\description{
Plots the incremental performance summary (mean and confidence interval) and
colors points by selection status.
}
\usage{
plot_incremental_metric(
  incremental_res,
  y_label = NULL,
  title = "Incremental performance",
  subtitle = NULL
)
}
\arguments{
\item{incremental_res}{Output from \code{IncrementalAnalysis}.}
\item{y_label}{Optional y-axis label. Defaults to the metric name stored in the object.}
\item{title}{Plot title.}
\item{subtitle}{Plot subtitle.}
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
plot_incremental_metric(inc)
}
}
