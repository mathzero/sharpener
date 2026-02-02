`%||%` <- function(a, b) if (!is.null(a)) a else b

.has_overreact <- requireNamespace("OverReact", quietly = TRUE)

mycols <- c("#000080", "#00bfff", "#fa8072", "#232333", "#7b68ee", "#006400", "#c71585")

theme_react_safe <- function(...) {
  if (.has_overreact) OverReact::theme_react(...) else ggplot2::theme_minimal()
}

scale_colour_imperial_safe <- function() {
  ggplot2::scale_colour_manual(values = mycols)
}

scale_fill_imperial_safe <- function() {
  ggplot2::scale_fill_manual(values = mycols)
}

.resolve_m <- function(incremental_res, m, m_select = c("best_test", "best_train")) {
  if (is.null(incremental_res$summary)) {
    stop("incremental_res must contain a summary table.")
  }
  m_select <- match.arg(m_select)
  if (is.null(m)) m <- m_select
  if (is.character(m)) {
    if (length(m) != 1L) stop("m must be length 1 when character.")
    m <- match.arg(m, c("best_test", "best_train"))
    if (m == "best_test") {
      return(as.integer(which.max(incremental_res$summary$metric_mean)))
    }
    return(as.integer(which.max(incremental_res$summary$metric_mean_train)))
  }
  m_vec <- as.integer(m)
  if (length(m_vec) == 0L) stop("m must be a non-empty numeric vector.")
  if (anyNA(m_vec)) stop("m contains NA values.")
  if (any(m_vec < 1L | m_vec > nrow(incremental_res$summary))) {
    stop("m is out of range for incremental_res$summary.")
  }
  m_vec <- m_vec[!duplicated(m_vec)]
  m_vec
}

.predictor_order <- function(incremental_res) {
  if (!is.null(incremental_res$predictors)) {
    return(as.character(incremental_res$predictors))
  }
  if (!is.null(incremental_res$summary) && "variable" %in% names(incremental_res$summary)) {
    labs <- as.character(incremental_res$summary$variable)
    return(gsub("^\\+", "", labs))
  }
  NULL
}

.labels_for_m <- function(m_vec,
                          incremental_res,
                          labels,
                          label_predictors,
                          label_style,
                          label_collapse,
                          label_include_m) {
  if (!is.null(labels)) {
    if (length(labels) == 1L && length(m_vec) > 1L) {
      labels <- rep(labels, length(m_vec))
    }
    if (length(labels) != length(m_vec)) {
      stop("labels must be length 1 or match length(m).")
    }
    return(as.character(labels))
  }

  summ_var <- NULL
  if (!is.null(incremental_res$summary) && "variable" %in% names(incremental_res$summary)) {
    summ_var <- as.character(incremental_res$summary$variable)
  }

  make_default <- function(m) {
    base <- NULL
    if (!is.null(summ_var) && m <= length(summ_var)) base <- summ_var[m]
    if (isTRUE(label_include_m)) {
      if (!is.null(base)) return(paste0("m=", m, " (", base, ")"))
      return(paste0("m=", m))
    }
    if (!is.null(base)) return(base)
    paste0("m=", m)
  }

  if (isTRUE(label_predictors)) {
    predictor_order <- .predictor_order(incremental_res)
    if (!is.null(predictor_order) && length(predictor_order) >= max(m_vec)) {
      return(vapply(m_vec, function(m) {
        preds <- predictor_order[seq_len(m)]
        label_part <- if (label_style == "added") tail(preds, 1) else paste(preds, collapse = label_collapse)
        if (isTRUE(label_include_m)) paste0("m=", m, ": ", label_part) else label_part
      }, character(1)))
    }
  }

  vapply(m_vec, make_default, character(1))
}


# IncrementalAnalysis -----------------------------------------------------

IncrementalAnalysis <- function(xdata, ydata, n_predictors = NULL, stability = NULL,
                                test_ratio = 0.3, k = 100, seed = 123, ncores = 4, verbose = TRUE,
                                ci_level = 0.95, roc_points = 101, max_predictors = 40,
                                store_data = TRUE) {
  if (!requireNamespace("future", quietly = TRUE)) stop("Package 'future' is required.")
  if (!requireNamespace("future.apply", quietly = TRUE)) stop("Package 'future.apply' is required.")

  set.seed(seed)

  xdf <- as.data.frame(xdata)
  if (is.null(n_predictors)) n_predictors <- min(ncol(xdf), max_predictors)
  if (!is.numeric(n_predictors) || length(n_predictors) != 1L || n_predictors < 1) {
    stop("n_predictors must be a positive integer.")
  }

  family <- NULL
  if (inherits(ydata, "Surv")) {
    family <- "cox"
  } else if (is.factor(ydata) || length(unique(na.omit(ydata))) == 2) {
    family <- "binomial"
  } else {
    family <- "gaussian"
  }

  metric <- if (family == "binomial") "AUC" else if (family == "cox") "C-index" else "R-squared"
  target_info <- NULL

  if (family == "binomial") {
    if (!requireNamespace("pROC", quietly = TRUE)) stop("Package 'pROC' is required for binomial outcomes.")

    target_vec <- ydata
    if (is.factor(target_vec)) {
      if (nlevels(target_vec) != 2) stop("Binary outcome factor must have exactly 2 levels.")
      pos_level <- levels(target_vec)[2]
      target_vec <- as.integer(target_vec == pos_level)
      target_info <- list(type = "factor", levels = levels(ydata), positive = pos_level)
    } else if (is.logical(target_vec)) {
      target_vec <- as.integer(target_vec)
      target_info <- list(type = "logical", levels = c(FALSE, TRUE), positive = TRUE)
    } else if (is.numeric(target_vec)) {
      u <- sort(unique(na.omit(as.numeric(target_vec))))
      if (length(u) != 2) stop("Binary outcome must have exactly two unique numeric values.")
      target_vec <- as.integer(as.numeric(target_vec) == u[2])
      target_info <- list(type = "numeric", levels = u, positive = u[2])
    } else stop("Unsupported ydata type for binomial.")

    if (anyNA(target_vec)) stop("ydata contains NA values; please impute or remove.")
    if (length(unique(target_vec)) < 2) stop("Binary outcome must contain both classes.")
  }

  if (family == "gaussian") {
    target_vec <- as.numeric(ydata)
    if (anyNA(target_vec)) stop("ydata contains NA values; please impute or remove.")
  }

  if (family == "cox") {
    if (!requireNamespace("survival", quietly = TRUE)) stop("Package 'survival' is required for Cox outcomes.")
    y_time  <- as.numeric(ydata[, 1])
    y_event <- as.integer(ydata[, 2])
    if (anyNA(y_time) || anyNA(y_event)) stop("Surv outcome contains NA values.")
    if (!all(y_event %in% c(0L, 1L))) stop("Surv event indicator must be coded 0/1.")
  }

  # Derive ordering and status from SHARP stability if provided
  if (!is.null(stability)) {
    selprop_raw <- sharp::SelectionProportions(stability)
    if (is.matrix(selprop_raw)) {
      if (ncol(selprop_raw) == 1) {
        selprop <- selprop_raw[, 1]
        names(selprop) <- rownames(selprop_raw)
      } else {
        selprop <- colMeans(selprop_raw, na.rm = TRUE)
        if (!is.null(colnames(selprop_raw))) names(selprop) <- colnames(selprop_raw)
      }
    } else {
      selprop <- as.numeric(selprop_raw)
      names(selprop) <- names(selprop_raw)
    }
    selprop <- selprop[!is.na(selprop)]
    selprop <- sort(selprop, decreasing = TRUE)

    selected_vars <- sharp::SelectedVariables(stability)
    selected_vars <- as.character(selected_vars)

    all_beta_names <- NULL
    if (!is.null(stability$Beta)) {
      dn <- dimnames(stability$Beta)
      if (length(dn) >= 2 && !is.null(dn[[2]])) all_beta_names <- dn[[2]]
    }
    if (is.null(all_beta_names)) all_beta_names <- colnames(xdf)

    unpenalised_vars <- setdiff(all_beta_names, names(selprop))

    penalised_order <- names(selprop)
    all_order <- c(unpenalised_vars, penalised_order)
    all_order <- all_order[all_order %in% colnames(xdf)]
    if (length(all_order) == 0) stop("No predictors found in stability object that match xdata.")

    all_predictors <- all_order[seq_len(min(length(all_order), n_predictors))]

    status_map <- setNames(rep("Not selected", length(all_order)), all_order)
    if (length(selected_vars)) status_map[intersect(selected_vars, names(status_map))] <- "Selected"
    if (length(unpenalised_vars)) status_map[intersect(unpenalised_vars, names(status_map))] <- "Unpenalised (forced)"
  } else {
    all_predictors <- colnames(xdf)[seq_len(min(n_predictors, ncol(xdf)))]
    status_map <- setNames(rep("Selected", length(all_predictors)), all_predictors)
  }

  if (family == "binomial") {
    data <- cbind.data.frame(target = target_vec, xdf)
  } else if (family == "gaussian") {
    data <- cbind.data.frame(target = target_vec, xdf)
  } else {
    data <- cbind.data.frame(time = y_time, event = y_event, xdf)
  }

  fpr_grid <- seq(0, 1, length.out = roc_points)
  spec_grid <- 1 - fpr_grid
  alpha <- 1 - ci_level
  q_lo <- alpha / 2
  q_hi <- 1 - alpha / 2

  quote_if_needed <- function(x) if (make.names(x) != x) paste0("`", x, "`") else x

  fill_na_linear <- function(v) {
    if (!anyNA(v)) return(v)
    if (all(is.na(v))) return(v)
    if (length(v) < 2L) return(v)
    good <- !is.na(v)
    if (sum(good) == 1L) {
      v[!good] <- v[good][1]
      return(v)
    }
    if (is.na(v[1])) v[1] <- v[which(good)[1]]
    if (is.na(v[length(v)])) v[length(v)] <- v[rev(which(good))[1]]
    na_idx <- which(is.na(v))
    if (length(na_idx) > 0) {
      v[!good] <- approx(x = which(good), y = v[good], xout = which(!good),
                         method = "linear", rule = 2, ties = "ordered")$y
    }
    v
  }

  safe_roc_stats <- function(response, predictor) {
    if (length(unique(response)) < 2) {
      return(list(auc = NA_real_, tpr = rep(NA_real_, length(fpr_grid))))
    }
    roc_obj <- try(
      pROC::roc(response = response, predictor = predictor, quiet = TRUE, direction = "auto"),
      silent = TRUE
    )
    if (inherits(roc_obj, "try-error")) {
      return(list(auc = NA_real_, tpr = rep(NA_real_, length(fpr_grid))))
    }
    auc_val <- try(as.numeric(pROC::auc(roc_obj)), silent = TRUE)
    if (inherits(auc_val, "try-error")) auc_val <- NA_real_

    sens_vec <- try(
      pROC::coords(roc_obj, x = spec_grid, input = "specificity", ret = "sensitivity", transpose = FALSE),
      silent = TRUE
    )
    if (inherits(sens_vec, "try-error")) {
      sens_vec <- rep(NA_real_, length(fpr_grid))
    } else {
      sens_vec <- fill_na_linear(as.numeric(unlist(sens_vec, use.names = FALSE)))
    }
    list(auc = as.numeric(auc_val), tpr = sens_vec)
  }

  safe_cindex <- function(df) {
    out <- try(
      survival::concordance(survival::Surv(time, event) ~ lp, data = df, reverse = TRUE)$concordance,
      silent = TRUE
    )
    if (inherits(out, "try-error")) return(NA_real_)
    as.numeric(out)
  }

  stratified_train_idx <- function(y, train_prop) {
    y <- as.integer(y)
    idx <- seq_along(y)
    cls <- sort(unique(y))
    train_idx <- integer(0)
    for (c in cls) {
      c_idx <- idx[y == c]
      n_c <- length(c_idx)
      n_tr <- floor(train_prop * n_c)
      if (n_tr < 1L || n_tr >= n_c) {
        return(sample.int(length(y), size = floor(train_prop * length(y))))
      }
      train_idx <- c(train_idx, sample(c_idx, n_tr))
    }
    train_idx
  }

  # Precompute k splits once, reuse across all model sizes
  set.seed(seed)
  split_list <- vector("list", k)
  for (i in seq_len(k)) {
    if (family == "binomial") {
      tr_idx <- stratified_train_idx(data$target, train_prop = 1 - test_ratio)
    } else if (family == "cox") {
      event_idx <- which(data$event == 1L)
      cens_idx  <- which(data$event == 0L)

      n_tr_event <- floor((1 - test_ratio) * length(event_idx))
      n_tr_cens  <- floor((1 - test_ratio) * length(cens_idx))

      if (n_tr_event < 1L || n_tr_cens < 1L) {
        n <- nrow(data)
        tr_idx <- sample.int(n, size = floor((1 - test_ratio) * n))
      } else {
        tr_idx <- c(sample(event_idx, size = n_tr_event),
                    sample(cens_idx,  size = n_tr_cens))
      }
    } else {
      n <- nrow(data)
      tr_idx <- sample.int(n, size = floor((1 - test_ratio) * n))
    }
    split_list[[i]] <- tr_idx
  }

  run_model <- function(m) {
    predictors_subset <- all_predictors[seq_len(m)]
    rhs <- paste(vapply(predictors_subset, quote_if_needed, character(1)), collapse = " + ")

    res_metric <- lapply(seq_along(split_list), function(i) {
      tr_idx <- split_list[[i]]
      train_data <- data[tr_idx, , drop = FALSE]
      test_data  <- data[-tr_idx, , drop = FALSE]

      if (family == "binomial") {
        form <- as.formula(paste("target ~", rhs))
        fit <- try(suppressWarnings(stats::glm(formula = form, data = train_data, family = stats::binomial())), silent = TRUE)
      } else if (family == "gaussian") {
        form <- as.formula(paste("target ~", rhs))
        fit <- try(stats::lm(formula = form, data = train_data), silent = TRUE)
      } else {
        form <- as.formula(paste("survival::Surv(time, event) ~", rhs))
        fit <- try(survival::coxph(formula = form, data = train_data, x = TRUE), silent = TRUE)
      }

      if (inherits(fit, "try-error")) {
        out <- list(beta_coefficients_names = predictors_subset,
                    beta_coefficients = rep(NA_real_, length(predictors_subset)),
                    metric_result = NA_real_,
                    metric_result_train = NA_real_)
        if (family == "binomial") {
          out$tpr <- rep(NA_real_, length(fpr_grid))
          out$tpr_train <- rep(NA_real_, length(fpr_grid))
        }
        return(out)
      }

      beta_coefficients <- stats::coef(fit)

      if (family == "binomial") {
        test_prob <- as.numeric(stats::predict(fit, newdata = test_data, type = "response"))
        roc_test <- safe_roc_stats(test_data$target, test_prob)

        train_prob <- as.numeric(stats::predict(fit, newdata = train_data, type = "response"))
        roc_train <- safe_roc_stats(train_data$target, train_prob)

        return(list(beta_coefficients_names = names(beta_coefficients),
                    beta_coefficients = as.vector(beta_coefficients),
                    metric_result = roc_test$auc,
                    metric_result_train = roc_train$auc,
                    tpr = roc_test$tpr,
                    tpr_train = roc_train$tpr))
      }

      if (family == "gaussian") {
        test_pred <- as.numeric(stats::predict(fit, newdata = test_data))
        metric_result_test <- suppressWarnings(stats::cor(test_data$target, test_pred, use = "complete.obs")^2)

        train_pred <- as.numeric(stats::predict(fit, newdata = train_data))
        metric_result_train <- suppressWarnings(stats::cor(train_data$target, train_pred, use = "complete.obs")^2)

        return(list(beta_coefficients_names = names(beta_coefficients),
                    beta_coefficients = as.vector(beta_coefficients),
                    metric_result = as.numeric(metric_result_test),
                    metric_result_train = as.numeric(metric_result_train)))
      }

      test_lp <- try(as.numeric(stats::predict(fit, newdata = test_data, type = "lp")), silent = TRUE)
      if (inherits(test_lp, "try-error")) test_lp <- rep(NA_real_, nrow(test_data))
      test_data$lp <- test_lp
      metric_result_test <- safe_cindex(test_data)

      train_lp <- try(as.numeric(stats::predict(fit, newdata = train_data, type = "lp")), silent = TRUE)
      if (inherits(train_lp, "try-error")) train_lp <- rep(NA_real_, nrow(train_data))
      train_data$lp <- train_lp
      metric_result_train <- safe_cindex(train_data)

      return(list(beta_coefficients_names = names(beta_coefficients),
                  beta_coefficients = as.vector(beta_coefficients),
                  metric_result = metric_result_test,
                  metric_result_train = metric_result_train))
    })

    coef_names <- unique(unlist(lapply(res_metric, function(z) z$beta_coefficients_names)))
    coef_names <- coef_names[!is.na(coef_names)]

    beta_mat <- do.call(rbind, lapply(res_metric, function(z) {
      b <- rep(NA_real_, length(coef_names))
      names(b) <- coef_names
      b[z$beta_coefficients_names] <- z$beta_coefficients
      unname(b)
    }))
    colnames(beta_mat) <- coef_names
    rownames(beta_mat) <- paste0("iter", seq_len(nrow(beta_mat)))

    metric_results <- as.numeric(sapply(res_metric, `[[`, "metric_result"))
    metric_results_train <- as.numeric(sapply(res_metric, `[[`, "metric_result_train"))

    final <- list(
      variable_names = coef_names,
      metric_mean = mean(metric_results, na.rm = TRUE),
      metric_sd   = stats::sd(metric_results, na.rm = TRUE),
      metric_lo   = stats::quantile(metric_results, probs = q_lo, na.rm = TRUE, names = FALSE),
      metric_hi   = stats::quantile(metric_results, probs = q_hi, na.rm = TRUE, names = FALSE),
      metric_mean_train = mean(metric_results_train, na.rm = TRUE),
      metric_sd_train   = stats::sd(metric_results_train, na.rm = TRUE),
      metric_lo_train   = stats::quantile(metric_results_train, probs = q_lo, na.rm = TRUE, names = FALSE),
      metric_hi_train   = stats::quantile(metric_results_train, probs = q_hi, na.rm = TRUE, names = FALSE),
      metric      = metric,
      metric_raw  = metric_results,
      metric_raw_train = metric_results_train,
      beta_coefficients_mean = colMeans(beta_mat, na.rm = TRUE),
      beta_coefficients_sd   = apply(beta_mat, 2, stats::sd, na.rm = TRUE),
      beta_coefficients      = beta_mat
    )

    if (family == "binomial") {
      tpr_mat <- do.call(rbind, lapply(res_metric, `[[`, "tpr"))
      if (!is.null(tpr_mat)) tpr_mat <- tpr_mat[rowSums(is.na(tpr_mat)) < ncol(tpr_mat), , drop = FALSE]
      if (!is.null(tpr_mat) && nrow(tpr_mat) > 0) {
        roc_ci_df <- data.frame(
          n_vars   = m,
          fpr      = fpr_grid,
          tpr_mean = colMeans(tpr_mat, na.rm = TRUE),
          tpr_lo   = apply(tpr_mat, 2, stats::quantile, probs = q_lo, na.rm = TRUE, names = FALSE),
          tpr_hi   = apply(tpr_mat, 2, stats::quantile, probs = q_hi, na.rm = TRUE, names = FALSE)
        )
      } else {
        roc_ci_df <- data.frame(n_vars = m, fpr = fpr_grid, tpr_mean = NA_real_, tpr_lo = NA_real_, tpr_hi = NA_real_)
      }
      final$roc_ci_df <- roc_ci_df

      tpr_mat_tr <- do.call(rbind, lapply(res_metric, `[[`, "tpr_train"))
      if (!is.null(tpr_mat_tr)) tpr_mat_tr <- tpr_mat_tr[rowSums(is.na(tpr_mat_tr)) < ncol(tpr_mat_tr), , drop = FALSE]
      if (!is.null(tpr_mat_tr) && nrow(tpr_mat_tr) > 0) {
        roc_ci_df_train <- data.frame(
          n_vars   = m,
          fpr      = fpr_grid,
          tpr_mean = colMeans(tpr_mat_tr, na.rm = TRUE),
          tpr_lo   = apply(tpr_mat_tr, 2, stats::quantile, probs = q_lo, na.rm = TRUE, names = FALSE),
          tpr_hi   = apply(tpr_mat_tr, 2, stats::quantile, probs = q_hi, na.rm = TRUE, names = FALSE)
        )
      } else {
        roc_ci_df_train <- data.frame(n_vars = m, fpr = fpr_grid, tpr_mean = NA_real_, tpr_lo = NA_real_, tpr_hi = NA_real_)
      }
      final$roc_ci_df_train <- roc_ci_df_train
    }

    final
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = ncores)
  details <- future.apply::future_lapply(seq_along(all_predictors), run_model, future.seed = TRUE)
  future::plan(future::sequential)

  added_vars <- all_predictors
  status_vec <- status_map[added_vars]

  metric_mean <- sapply(details, `[[`, "metric_mean")
  metric_sd   <- sapply(details, `[[`, "metric_sd")
  metric_lo   <- sapply(details, `[[`, "metric_lo")
  metric_hi   <- sapply(details, `[[`, "metric_hi")

  metric_mean_train <- sapply(details, `[[`, "metric_mean_train")
  metric_sd_train   <- sapply(details, `[[`, "metric_sd_train")
  metric_lo_train   <- sapply(details, `[[`, "metric_lo_train")
  metric_hi_train   <- sapply(details, `[[`, "metric_hi_train")

  variable_labels <- c(added_vars[1], paste0("+", added_vars[-1]))
  summary_output <- data.frame(
    n_vars   = seq_along(all_predictors),
    variable = factor(variable_labels, levels = variable_labels),
    metric_mean, metric_sd, metric_lo, metric_hi,
    metric_mean_train = metric_mean_train,
    metric_sd_train   = metric_sd_train,
    metric_lo_train   = metric_lo_train,
    metric_hi_train   = metric_hi_train,
    status = unname(status_vec)
  )

  roc_ci <- NULL
  roc_ci_train <- NULL
  if (family == "binomial") {
    roc_ci <- do.call(rbind, lapply(seq_along(details), function(i) {
      df <- details[[i]]$roc_ci_df
      if (is.null(df)) return(NULL)
      df$label <- as.character(summary_output$variable[i])
      df
    }))
    if (!is.null(roc_ci)) rownames(roc_ci) <- NULL

    roc_ci_train <- do.call(rbind, lapply(seq_along(details), function(i) {
      df <- details[[i]]$roc_ci_df_train
      if (is.null(df)) return(NULL)
      df$label <- as.character(summary_output$variable[i])
      df
    }))
    if (!is.null(roc_ci_train)) rownames(roc_ci_train) <- NULL
  }

  call <- list(
    test_ratio = test_ratio,
    k = k,
    metric = metric,
    ci_level = ci_level,
    roc_points = roc_points,
    family = family,
    n_predictors = n_predictors,
    seed = seed,
    store_data = store_data
  )
  res <- list(
    summary = summary_output,
    roc_ci = roc_ci,
    roc_ci_train = roc_ci_train,
    call = call,
    details = details,
    family = family,
    metric = metric,
    predictors = all_predictors,
    status = status_map,
    target_info = target_info
  )
  if (isTRUE(store_data)) {
    res$xdata <- xdata
    res$ydata <- ydata
  }
  class(res) <- c("sharpener_incremental", "list")

  if (verbose) {
    cat("\n--- Incremental Analysis Complete ---\n")
    cat(sprintf("Metric: %s\n", metric))
    cat(sprintf("Numbers of predictors tested: %d\n", nrow(summary_output)))
    if (family == "cox") cat("Cox evaluation uses survival::concordance with reverse=TRUE.\n")
    if (!is.null(roc_ci)) cat("Returned: roc_ci for empirical ROC bands on test data per combination.\n")
    if (!is.null(roc_ci_train)) cat("Returned: roc_ci_train for empirical ROC bands on training data per combination.\n")
  }
  res
}


# S3 helper and minimal-argument plotters ---------------------------------

print.sharpener_incremental <- function(x, ...) {
  cat("<sharpener_incremental>\n")
  cat("  family: ", x$family %||% "unknown", "\n", sep = "")
  if (!is.null(x$summary)) {
    cat("  increments: ", nrow(x$summary), "\n", sep = "")
  }
  invisible(x)
}

plot_incremental_roc <- function(incremental_res,
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
                                 ...) {
  if (!is.list(incremental_res) || is.null(incremental_res$summary)) {
    stop("incremental_res must be the output from IncrementalAnalysis().")
  }
  m_select <- match.arg(m_select)
  which <- match.arg(which)
  label_style <- match.arg(label_style)

  m_vec <- .resolve_m(incremental_res, m, m_select)
  label_vec <- .labels_for_m(
    m_vec,
    incremental_res,
    labels,
    label_predictors,
    label_style,
    label_collapse,
    label_include_m
  )

  ci_df_all <- if (which == "test") incremental_res$roc_ci else incremental_res$roc_ci_train
  if (is.null(ci_df_all)) stop("ROC data not available; binomial outcome required.")

  df_list <- lapply(seq_along(m_vec), function(i) {
    df <- ci_df_all[ci_df_all$n_vars == m_vec[i], , drop = FALSE]
    if (nrow(df) == 0L) stop("No ROC rows found for the requested m.")
    df$study <- label_vec[i]
    df
  })
  ci_df <- do.call(rbind, df_list)
  rownames(ci_df) <- NULL

  if (is.null(title)) {
    title <- if (which == "test") "Incremental ROC (test)" else "Incremental ROC (train)"
  }
  if (is.null(subtitle) && length(m_vec) == 1L) {
    subtitle <- label_vec[1]
  }

  plotROC_from_ci(ci_df, title = title, subtitle = subtitle, ...)
}

plot_incremental_survival <- function(incremental_res,
                                      m = NULL,
                                      m_select = c("best_test", "best_train"),
                                      ...) {
  plot_predicted_vs_actual_survival(
    incremental_res,
    xdata = NULL,
    ydata = NULL,
    m = m,
    m_select = m_select,
    ...
  )
}


# Predicted vs actual survival curve --------------------------------------

plot_predicted_vs_actual_survival <- function(incremental_res,
                                              xdata = NULL,
                                              ydata = NULL,
                                              m = NULL,
                                              m_select = c("best_test", "best_train"),
                                              labels = NULL,
                                              label_predictors = FALSE,
                                              label_style = c("full", "added"),
                                              label_collapse = ", ",
                                              label_include_m = TRUE,
                                              n_draws = 200,
                                              ci_level = NULL,
                                              times = NULL,
                                              curve_type = c("marginal", "mean_lp"),
                                              seed = 1,
                                              show_km_ci = FALSE,
                                              show_band = TRUE,
                                              band_alpha = 0.2,
                                              title = NULL,
                                              subtitle = NULL) {
  if (!requireNamespace("survival", quietly = TRUE)) stop("Package 'survival' is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")

  m_select <- match.arg(m_select)
  curve_type <- match.arg(curve_type)
  label_style <- match.arg(label_style)

  if (is.null(xdata) || is.null(ydata)) {
    if (!is.null(incremental_res$xdata) && !is.null(incremental_res$ydata)) {
      xdata <- incremental_res$xdata
      ydata <- incremental_res$ydata
    } else {
      stop("xdata/ydata not provided and not stored in incremental_res; rerun IncrementalAnalysis(store_data = TRUE) or pass data explicitly.")
    }
  }
  if (!inherits(ydata, "Surv")) stop("ydata must be a survival::Surv object for this plot.")
  if (!is.list(incremental_res) || is.null(incremental_res$details) || is.null(incremental_res$summary)) {
    stop("incremental_res must be the output from IncrementalAnalysis().")
  }

  xdf <- as.data.frame(xdata)
  df <- cbind.data.frame(
    time = as.numeric(ydata[, 1]),
    event = as.integer(ydata[, 2]),
    xdf
  )

  if (anyNA(df$time) || anyNA(df$event)) stop("Surv outcome contains NA values.")
  if (!all(df$event %in% c(0L, 1L))) stop("Surv event indicator must be coded 0/1.")

  if (is.null(ci_level)) {
    if (!is.null(incremental_res$call$ci_level)) ci_level <- incremental_res$call$ci_level
    else ci_level <- 0.95
  }
  alpha <- 1 - ci_level
  q_lo <- alpha / 2
  q_hi <- 1 - alpha / 2

  m_vec <- .resolve_m(incremental_res, m, m_select)
  predictor_order <- .predictor_order(incremental_res)
  if (is.null(predictor_order) || length(predictor_order) < max(m_vec)) {
    stop("Could not determine predictor order for requested m values.")
  }

  label_vec <- .labels_for_m(
    m_vec,
    incremental_res,
    labels,
    label_predictors,
    label_style,
    label_collapse,
    label_include_m
  )

  # Choose time grid
  km_fit <- survival::survfit(survival::Surv(time, event) ~ 1, data = df)
  if (is.null(times)) {
    times_grid <- km_fit$time
  } else {
    times_grid <- sort(unique(as.numeric(times)))
  }
  if (length(times_grid) < 2L) stop("times grid must contain at least 2 time points.")

  # Actual (Kaplan–Meier) curve
  km_sum <- summary(km_fit, times = times_grid, extend = TRUE)
  actual_surv <- as.numeric(km_sum$surv)
  actual_lo <- rep(NA_real_, length(times_grid))
  actual_hi <- rep(NA_real_, length(times_grid))
  if (show_km_ci) {
    actual_lo <- as.numeric(km_sum$lower)
    actual_hi <- as.numeric(km_sum$upper)
  }
  actual_df <- data.frame(
    time = times_grid,
    actual = actual_surv,
    actual_lo = actual_lo,
    actual_hi = actual_hi
  )

  safe_basehaz <- function(time_vec, event_vec, lp_vec) {
    fit_off <- try(
      survival::coxph(survival::Surv(time_vec, event_vec) ~ offset(lp_vec), ties = "breslow"),
      silent = TRUE
    )
    if (inherits(fit_off, "try-error")) return(NULL)
    bh <- try(survival::basehaz(fit_off, centered = FALSE), silent = TRUE)
    if (inherits(bh, "try-error")) return(NULL)
    bh
  }

  pred_list <- vector("list", length(m_vec))
  predictors_out <- vector("list", length(m_vec))
  draws_used <- integer(length(m_vec))

  for (i in seq_along(m_vec)) {
    m_i <- m_vec[i]
    predictor_set <- predictor_order[seq_len(m_i)]
    predictors_out[[i]] <- predictor_set

    det <- incremental_res$details[[m_i]]
    if (is.null(det$beta_coefficients) || !is.matrix(det$beta_coefficients)) {
      stop("incremental_res$details[[m]]$beta_coefficients must be a matrix.")
    }

    beta_mat <- det$beta_coefficients
    beta_cols <- colnames(beta_mat)
    if (is.null(beta_cols)) stop("beta_coefficients matrix must have column names.")

    beta_cols_use <- setdiff(beta_cols, "(Intercept)")
    beta_cols_use <- beta_cols_use[beta_cols_use %in% predictor_set]
    if (length(beta_cols_use) == 0L) stop("No overlapping coefficient columns found for predictor set.")

    beta_sub <- beta_mat[, beta_cols_use, drop = FALSE]

    ok_draw <- rowSums(is.na(beta_sub)) == 0L
    beta_sub <- beta_sub[ok_draw, , drop = FALSE]
    if (nrow(beta_sub) == 0L) stop("No valid coefficient draws available after removing NA draws.")

    set.seed(seed + m_i)
    if (!is.null(n_draws) && is.finite(n_draws)) {
      n_draws_i <- min(as.integer(n_draws), nrow(beta_sub))
      draw_idx <- sample.int(nrow(beta_sub), size = n_draws_i, replace = FALSE)
      beta_sub <- beta_sub[draw_idx, , drop = FALSE]
    }
    draws_used[i] <- nrow(beta_sub)

    X <- as.matrix(xdf[, beta_cols_use, drop = FALSE])
    mode(X) <- "numeric"

    predict_curve_one <- function(beta_vec) {
      lp <- as.numeric(X %*% beta_vec)
      bh <- safe_basehaz(df$time, df$event, lp)
      if (is.null(bh)) return(rep(NA_real_, length(times_grid)))

      H0 <- stats::approx(x = bh$time, y = bh$hazard, xout = times_grid, method = "linear", rule = 2)$y
      if (anyNA(H0)) return(rep(NA_real_, length(times_grid)))

      exp_lp <- exp(lp)

      if (curve_type == "mean_lp") {
        exp_lp_ref <- exp(mean(lp, na.rm = TRUE))
        out <- exp(-H0 * exp_lp_ref)
        return(as.numeric(out))
      }

      out <- vapply(H0, function(h) mean(exp(-h * exp_lp), na.rm = TRUE), numeric(1))
      as.numeric(out)
    }

    pred_mat <- matrix(NA_real_, nrow = nrow(beta_sub), ncol = length(times_grid))
    for (j in seq_len(nrow(beta_sub))) {
      pred_mat[j, ] <- predict_curve_one(as.numeric(beta_sub[j, ]))
    }

    pred_mat <- pred_mat[rowSums(is.na(pred_mat)) < ncol(pred_mat), , drop = FALSE]
    if (nrow(pred_mat) == 0L) stop("All predicted curves failed; check coefficient draws and data integrity.")

    pred_med <- apply(pred_mat, 2, stats::median, na.rm = TRUE)
    pred_lo  <- apply(pred_mat, 2, stats::quantile, probs = q_lo, na.rm = TRUE, names = FALSE)
    pred_hi  <- apply(pred_mat, 2, stats::quantile, probs = q_hi, na.rm = TRUE, names = FALSE)

    pred_list[[i]] <- data.frame(
      time = times_grid,
      pred = as.numeric(pred_med),
      pred_lo = as.numeric(pred_lo),
      pred_hi = as.numeric(pred_hi),
      m = m_i,
      label = label_vec[i]
    )
  }

  pred_df <- do.call(rbind, pred_list)
  pred_df$label <- factor(pred_df$label, levels = label_vec)

  if (is.null(title)) title <- "Predicted vs actual survival"
  if (is.null(subtitle) && length(m_vec) == 1L) {
    subtitle <- paste0(
      "Model size m = ", m_vec[1],
      " | Predicted CI from refit draws (n = ", draws_used[1], ")",
      " | Curve type: ", curve_type
    )
  }

  p <- ggplot2::ggplot(pred_df, ggplot2::aes(x = time, y = pred, colour = label)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      x = "Time",
      y = "Survival probability",
      title = title,
      subtitle = subtitle,
      colour = NULL,
      fill = NULL
    ) +
    ggplot2::theme_minimal()

  if (isTRUE(show_band)) {
    p <- p + ggplot2::geom_ribbon(
      data = pred_df,
      ggplot2::aes(ymin = pred_lo, ymax = pred_hi, fill = label),
      alpha = band_alpha,
      colour = NA
    )
  }

  p <- p + ggplot2::geom_step(
    data = actual_df,
    ggplot2::aes(x = time, y = actual),
    inherit.aes = FALSE,
    linewidth = 1,
    linetype = 2,
    colour = "black"
  )

  if (show_km_ci) {
    p <- p + ggplot2::geom_ribbon(
      data = actual_df,
      ggplot2::aes(x = time, ymin = actual_lo, ymax = actual_hi),
      inherit.aes = FALSE,
      alpha = 0.1
    )
  }

  list(
    plot = p,
    data = list(actual = actual_df, predicted = pred_df),
    m = m_vec,
    labels = label_vec,
    predictors = predictors_out,
    n_draws_used = draws_used,
    curve_type = curve_type
  )
}


# ROC plotting ------------------------------------------------------------
plotROC_from_ci <- function(ci_df,
                            title = "ROC",
                            subtitle = NULL,
                            ci_label = "95% CI",
                            digits = 3,
                            annotate_single = FALSE,
                            annotate_pos = c(0.65, 0.15)) {
  stopifnot(all(c("fpr", "tpr_mean") %in% names(ci_df)))

  trap_auc <- function(x, y) {
    o <- order(x)
    x <- x[o]; y <- y[o]
    dx <- diff(x)
    ym <- (y[-1] + y[-length(y)]) / 2
    sum(dx * ym, na.rm = TRUE)
  }

  df <- ci_df
  if (!"study" %in% names(df)) df$study <- "Pooled"

  has_lo <- "tpr_lo" %in% names(df)
  has_hi <- "tpr_hi" %in% names(df)

  stats <- df %>%
    dplyr::group_by(study) %>%
    dplyr::summarise(
      auc    = trap_auc(fpr, tpr_mean),
      auc_lo = if (has_lo) trap_auc(fpr, tpr_lo) else NA_real_,
      auc_hi = if (has_hi) trap_auc(fpr, tpr_hi) else NA_real_,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      lab = dplyr::if_else(
        !is.na(auc_lo) & !is.na(auc_hi),
        sprintf("%s  AUC %.3f [%s %.3f, %.3f]", study, round(auc, digits), ci_label, round(auc_lo, digits), round(auc_hi, digits)),
        sprintf("%s  AUC %.3f", study, round(auc, digits))
      )
    )

  df_plot <- df %>%
    dplyr::left_join(stats, by = "study") %>%
    dplyr::mutate(grp_label = factor(lab, levels = stats$lab)) %>%
    dplyr::mutate(tpr_mean = ifelse(fpr == 0.00, 0, tpr_mean))

  g <- ggplot2::ggplot(df_plot, ggplot2::aes(x = fpr)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = if (has_lo) tpr_lo else NA_real_,
                   ymax = if (has_hi) tpr_hi else NA_real_,
                   fill = grp_label),
      alpha = 0.12, na.rm = TRUE
    ) +
    ggplot2::geom_line(ggplot2::aes(y = tpr_mean, colour = grp_label)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_react_safe() +
    scale_colour_imperial_safe() +
    scale_fill_imperial_safe() +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "False positive rate",
      y = "True positive rate",
      colour = NULL, fill = NULL
    ) +
    ggplot2::theme(legend.position = c(0.65, 0.20))

  if (annotate_single && (length(levels(df_plot$grp_label)) == 1L)) {
    lab_txt <- as.character(levels(df_plot$grp_label))
    g <- g + ggplot2::annotate(
      "label", x = annotate_pos[1], y = annotate_pos[2],
      label = gsub("^Pooled\\s+", "", lab_txt),
      hjust = 0, vjust = 0, label.size = 0.2
    )
  }

  g
}
