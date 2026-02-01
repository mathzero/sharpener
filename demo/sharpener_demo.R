# Demo: sharpener + sharp variable selection pipeline

suppressPackageStartupMessages({
  library(sharp)
  library(survival)
})

set.seed(1)

# ---- 1) Binomial example -------------------------------------------------

n <- 200
p <- 20
x <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(x) <- paste0("x", seq_len(p))

beta <- c(rep(1, 5), rep(0, p - 5))
linpred <- as.numeric(x %*% beta)
prob <- 1 / (1 + exp(-linpred))
y_bin <- rbinom(n, size = 1, prob = prob)

vs_bin <- sharp::VariableSelection(
  xdata = x,
  ydata = y_bin,
  family = "binomial",
  K = 20,
  seed = 1,
  n_cores = 1,
  verbose = FALSE
)

inc_bin <- IncrementalAnalysis(
  xdata = x,
  ydata = y_bin,
  stability = vs_bin,
  n_predictors = 10,
  k = 20,
  test_ratio = 0.3,
  ncores = 1,
  seed = 1,
  verbose = TRUE
)

best_m_bin <- which.max(inc_bin$summary$metric_mean)
roc_df_best <- subset(inc_bin$roc_ci, n_vars == best_m_bin)
print(plotROC_from_ci(
  roc_df_best,
  title = "Incremental ROC (test)",
  subtitle = paste("m =", best_m_bin)
))

# ---- 2) Cox survival example --------------------------------------------

set.seed(2)

n2 <- 250
p2 <- 20
x2 <- matrix(rnorm(n2 * p2), nrow = n2, ncol = p2)
colnames(x2) <- paste0("x", seq_len(p2))

beta2 <- c(rep(0.7, 5), rep(0, p2 - 5))
linpred2 <- as.numeric(x2 %*% beta2)

base_rate <- 0.08
cens_rate <- 0.05

time_event <- rexp(n2, rate = base_rate * exp(linpred2))
time_cens <- rexp(n2, rate = cens_rate)

time_obs <- pmin(time_event, time_cens)
event <- as.integer(time_event <= time_cens)

y_surv <- survival::Surv(time_obs, event)

vs_surv <- sharp::VariableSelection(
  xdata = x2,
  ydata = y_surv,
  family = "cox",
  K = 20,
  seed = 2,
  n_cores = 1,
  verbose = FALSE
)

inc_surv <- IncrementalAnalysis(
  xdata = x2,
  ydata = y_surv,
  stability = vs_surv,
  n_predictors = 10,
  k = 20,
  test_ratio = 0.3,
  ncores = 1,
  seed = 2,
  verbose = TRUE
)

surv_plot <- plot_predicted_vs_actual_survival(
  inc_surv,
  x2,
  y_surv,
  m_select = "best_test"
)
print(surv_plot$plot)
