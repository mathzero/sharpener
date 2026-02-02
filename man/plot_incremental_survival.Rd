\name{plot_incremental_survival}
\alias{plot_incremental_survival}
\title{Plot predicted vs actual survival from an incremental analysis}
\description{
Convenience wrapper that plots predicted vs actual survival using the output
from \code{IncrementalAnalysis}.
}
\usage{
plot_incremental_survival(
  incremental_res,
  m = NULL,
  m_select = c("best_test", "best_train"),
  ...
)
}
\arguments{
\item{incremental_res}{Output from \code{IncrementalAnalysis}.}
\item{m}{Increment to plot; if \code{NULL}, selects by \code{m_select}.}
\item{m_select}{Selection rule for \code{m} when \code{m} is \code{NULL}.}
\item{...}{Additional arguments passed to \code{plot_predicted_vs_actual_survival}.}
}
\value{
A list containing a ggplot object and the data used to draw the curves.
}
\examples{
\dontrun{
set.seed(2)
X <- matrix(rnorm(200 * 10), nrow = 200)
colnames(X) <- paste0("x", 1:10)
linpred <- as.numeric(X %*% c(rep(0.5, 3), rep(0, 7)))
time_event <- rexp(200, rate = 0.1 * exp(linpred))
time_cens <- rexp(200, rate = 0.05)
time_obs <- pmin(time_event, time_cens)
event <- as.integer(time_event <= time_cens)
y <- survival::Surv(time_obs, event)
vs <- sharp::VariableSelection(xdata = X, ydata = y, family = "cox", K = 20)
inc <- IncrementalAnalysis(X, y, stability = vs, n_predictors = 5, k = 20, ncores = 1)
plot_incremental_survival(inc)
}
}
