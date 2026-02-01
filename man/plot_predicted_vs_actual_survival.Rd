\name{plot_predicted_vs_actual_survival}
\alias{plot_predicted_vs_actual_survival}
\title{Predicted vs actual survival curve}
\description{
Creates a predicted survival curve (with empirical uncertainty bands) and
compares it to the Kaplan--Meier curve using the output of \code{IncrementalAnalysis}.
}
\usage{
plot_predicted_vs_actual_survival(
  incremental_res,
  xdata,
  ydata,
  m = NULL,
  m_select = c("best_test", "best_train"),
  n_draws = 200,
  ci_level = NULL,
  times = NULL,
  curve_type = c("marginal", "mean_lp"),
  seed = 1,
  show_km_ci = FALSE
)
}
\arguments{
\item{incremental_res}{Output from \code{IncrementalAnalysis}.}
\item{xdata}{Predictor matrix or data frame used in the model.}
\item{ydata}{A \code{survival::Surv} object.}
\item{m}{Number of predictors to use; if \code{NULL}, selects by \code{m_select}.}
\item{m_select}{Selection rule for \code{m} when \code{m} is \code{NULL}.}
\item{n_draws}{Number of coefficient draws used to build prediction bands.}
\item{ci_level}{Confidence level for prediction bands.}
\item{times}{Optional vector of time points for evaluation.}
\item{curve_type}{\code{"marginal"} or \code{"mean_lp"} prediction type.}
\item{seed}{Random seed for draw selection.}
\item{show_km_ci}{Logical; include KM confidence bands.}
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
plot_predicted_vs_actual_survival(inc, X, y)
}
}
