# Univariable screening + volcano plotting utilities

.null_coalesce <- function(x, y) if (is.null(x)) y else x

.one_hot_xdata <- function(xdata, one_hot = TRUE) {
  if (!isTRUE(one_hot)) return(list(xdata = xdata, encoded = FALSE))

  if (is.data.frame(xdata)) {
    is_cat <- vapply(xdata, function(x) is.factor(x) || is.character(x) || is.logical(x), logical(1))
    if (any(is_cat)) {
      xdf <- xdata
      xdf[is_cat] <- lapply(xdf[is_cat], function(x) {
        if (is.character(x) || is.logical(x)) factor(x) else x
      })
      mm <- stats::model.matrix(~ . - 1, data = xdf)
      return(list(xdata = as.data.frame(mm), encoded = TRUE))
    }
    if (!all(vapply(xdata, is.numeric, logical(1)))) {
      stop("xdata must be numeric; one-hot encode categorical predictors before modeling.")
    }
    return(list(xdata = xdata, encoded = FALSE))
  }

  if (is.matrix(xdata)) {
    if (!is.numeric(xdata)) {
      stop("xdata must be a numeric matrix; one-hot encode categorical predictors before modeling.")
    }
    return(list(xdata = xdata, encoded = FALSE))
  }

  stop("xdata must be a data.frame or numeric matrix.")
}

.detect_family_y <- function(y, family = NULL) {
  if (!is.null(family)) {
    if (inherits(family, "family")) return(tolower(family$family))
    fam <- tolower(as.character(family))
    if (fam %in% c("coxph", "survival")) fam <- "cox"
    return(fam)
  }
  if (inherits(y, "Surv")) return("cox")
  if (is.factor(y)) return(if (nlevels(y) == 2L) "binomial" else "gaussian")
  if (is.logical(y)) return("binomial")
  if (is.numeric(y)) {
    u <- unique(stats::na.omit(y))
    if (length(u) == 2L && all(sort(u) %in% c(0, 1))) return("binomial")
  }
  "gaussian"
}

.resolve_selected_vars <- function(selected_raw, all_names) {
  if (is.null(selected_raw)) return(character(0))

  if (is.data.frame(selected_raw) || is.matrix(selected_raw)) {
    vec <- as.numeric(selected_raw)
    if (length(vec) == length(all_names) && all(vec %in% c(0, 1))) {
      return(all_names[vec > 0])
    }
    return(as.character(vec))
  }

  if (is.logical(selected_raw)) {
    if (length(selected_raw) == length(all_names)) {
      return(all_names[selected_raw])
    }
    return(as.character(which(selected_raw)))
  }

  if (is.numeric(selected_raw)) {
    if (!is.null(names(selected_raw)) && length(names(selected_raw)) == length(selected_raw)) {
      return(names(selected_raw)[selected_raw > 0])
    }
    if (length(selected_raw) == length(all_names) && all(selected_raw %in% c(0, 1))) {
      return(all_names[selected_raw > 0])
    }
    if (all(selected_raw %in% seq_len(length(all_names)))) {
      return(all_names[selected_raw])
    }
    return(as.character(selected_raw))
  }

  as.character(selected_raw)
}

.resolve_selection_threshold <- function(stability) {
  if (!inherits(stability, c("variable_selection", "structural_model"))) return(NA_real_)

  thr <- NA_real_
  argmax <- try(sharp::Argmax(stability), silent = TRUE)
  if (!inherits(argmax, "try-error")) {
    if (is.matrix(argmax) && "pi" %in% colnames(argmax)) {
      thr <- as.numeric(argmax[1, "pi"])
    } else if (is.numeric(argmax) && length(argmax) >= 2L) {
      thr <- as.numeric(argmax[2])
    }
  }

  if (!is.finite(thr)) {
    argmax_id <- try(sharp::ArgmaxId(stability), silent = TRUE)
    if (!inherits(argmax_id, "try-error") && !is.null(stability$params$pi_list)) {
      thr <- as.numeric(stability$params$pi_list[argmax_id[1, 2]])
    }
  }

  if (!is.finite(thr) && !is.null(stability$params$tau)) {
    thr <- as.numeric(stability$params$tau)
  }

  thr
}

.build_model_frame <- function(xdata, ydata, covariates = NULL, data = NULL) {
  xdf <- as.data.frame(xdata)
  if (nrow(xdf) == 0L) stop("xdata has zero rows.")
  if (is.null(colnames(xdf))) colnames(xdf) <- paste0("x", seq_len(ncol(xdf)))
  if (length(ydata) != nrow(xdf)) stop("ydata length must match nrow(xdata).")

  base_df <- if (is.null(data)) data.frame(.y = ydata) else data
  if (nrow(base_df) != nrow(xdf)) stop("data must have the same number of rows as xdata.")
  base_df$.y <- ydata

  cov_df <- NULL
  cov_names <- character(0)
  if (is.character(covariates)) {
    cov_names <- covariates
    missing_cov <- setdiff(cov_names, names(base_df))
    if (length(missing_cov)) stop("covariates not found in data: ", paste(missing_cov, collapse = ", "))
  } else if (is.data.frame(covariates) || is.matrix(covariates)) {
    cov_df <- as.data.frame(covariates)
    if (is.null(colnames(cov_df))) colnames(cov_df) <- paste0("cov", seq_len(ncol(cov_df)))
    cov_names <- colnames(cov_df)
  } else if (!is.null(covariates)) {
    stop("covariates must be NULL, a character vector, or a data.frame/matrix.")
  }

  all_names <- c(names(base_df), names(cov_df), names(xdf))
  if (anyDuplicated(all_names)) {
    dup <- unique(all_names[duplicated(all_names)])
    stop("Duplicate column names after combining data/covariates/xdata: ", paste(dup, collapse = ", "))
  }

  df <- if (is.null(cov_df)) {
    cbind(base_df, xdf)
  } else {
    cbind(base_df, cov_df, xdf)
  }

  list(df = df, predictors = colnames(xdf), covariates = cov_names)
}

.fit_univariable_single <- function(df,
                                    predictor,
                                    covariates = NULL,
                                    random = NULL,
                                    family = c("gaussian", "binomial", "cox"),
                                    verbose = TRUE) {
  family <- match.arg(family)

  cols <- unique(c(".y", predictor, covariates, random))
  cols <- cols[cols %in% names(df)]
  d <- df[stats::complete.cases(df[, cols, drop = FALSE]), , drop = FALSE]
  n <- nrow(d)
  if (n == 0L) {
    return(tibble::tibble(
      predictor = predictor,
      term = NA_character_,
      level = NA_character_,
      estimate = NA_real_,
      std_error = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      n = 0L,
      family = family,
      model = NA_character_,
      status = "fail",
      message = "all rows removed for missingness"
    ))
  }

  rhs_parts <- c(predictor, covariates)
  model <- NULL

  if (!is.null(random) && length(random)) {
    missing_re <- setdiff(random, names(d))
    if (length(missing_re)) stop("random effects not found in data: ", paste(missing_re, collapse = ", "))
    if (family == "cox" && verbose) {
      warning("Random effects are not supported for Cox models; ignoring random.")
    } else if (family != "cox") {
      rhs_parts <- c(rhs_parts, paste0("(1|", random, ")"))
    }
  }

  rhs <- if (length(rhs_parts)) paste(rhs_parts, collapse = " + ") else "1"
  form <- stats::as.formula(paste(".y ~", rhs))

  fit <- NULL
  fit_msg <- NA_character_

  if (family == "cox") {
    if (!requireNamespace("survival", quietly = TRUE)) stop("Package 'survival' is required for Cox models.")
    fit <- try(survival::coxph(form, data = d), silent = TRUE)
    model <- "coxph"
  } else if (!is.null(random) && length(random)) {
    if (!requireNamespace("lme4", quietly = TRUE)) stop("Package 'lme4' is required for random effects models.")
    if (family == "binomial") {
      fit <- try(lme4::glmer(form, data = d, family = stats::binomial()), silent = TRUE)
      model <- "glmer"
    } else {
      if (requireNamespace("lmerTest", quietly = TRUE)) {
        fit <- try(lmerTest::lmer(form, data = d), silent = TRUE)
        model <- "lmerTest"
      } else {
        fit <- try(lme4::lmer(form, data = d), silent = TRUE)
        model <- "lmer"
      }
    }
  } else if (family == "binomial") {
    fit <- try(stats::glm(form, data = d, family = stats::binomial()), silent = TRUE)
    model <- "glm"
  } else {
    fit <- try(stats::lm(form, data = d), silent = TRUE)
    model <- "lm"
  }

  if (inherits(fit, "try-error")) {
    fit_msg <- conditionMessage(attr(fit, "condition"))
    return(tibble::tibble(
      predictor = predictor,
      term = NA_character_,
      level = NA_character_,
      estimate = NA_real_,
      std_error = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      n = n,
      family = family,
      model = model,
      status = "fail",
      message = fit_msg
    ))
  }

  sm <- summary(fit)
  coef_mat <- sm$coefficients
  if (is.null(coef_mat)) {
    return(tibble::tibble(
      predictor = predictor,
      term = NA_character_,
      level = NA_character_,
      estimate = NA_real_,
      std_error = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      n = n,
      family = family,
      model = model,
      status = "fail",
      message = "no coefficients returned"
    ))
  }

  coef_df <- tibble::as_tibble(coef_mat, rownames = "term")
  stat_col <- dplyr::case_when(
    "t value" %in% names(coef_df) ~ "t value",
    "z value" %in% names(coef_df) ~ "z value",
    "z" %in% names(coef_df) ~ "z",
    TRUE ~ NA_character_
  )
  p_col <- grep("Pr\\(>.*\\)", names(coef_df), value = TRUE)[1]

  coef_df <- coef_df %>%
    dplyr::rename(estimate = Estimate, std_error = `Std. Error`) %>%
    dplyr::mutate(
      statistic = if (!is.na(stat_col)) .data[[stat_col]] else NA_real_,
      p_value = if (!is.na(p_col)) .data[[p_col]] else NA_real_
    )

  coef_df <- coef_df %>%
    dplyr::filter(grepl(paste0("^", predictor), .data$term) & !grepl(":", .data$term, fixed = TRUE))

  if (nrow(coef_df) == 0L) {
    return(tibble::tibble(
      predictor = predictor,
      term = NA_character_,
      level = NA_character_,
      estimate = NA_real_,
      std_error = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      n = n,
      family = family,
      model = model,
      status = "fail",
      message = "predictor term not found in model coefficients"
    ))
  }

  coef_df %>%
    dplyr::transmute(
      predictor = predictor,
      term,
      level = sub(paste0("^", predictor), "", term),
      estimate,
      std_error,
      statistic,
      p_value,
      n = n,
      family = family,
      model = model,
      status = "ok",
      message = NA_character_
    )
}

#' Run univariable screening across predictors
univariableScreen <- function(xdata,
                              ydata,
                              family = NULL,
                              covariates = NULL,
                              random = NULL,
                              data = NULL,
                              ncores = 1,
                              seed = 1,
                              verbose = TRUE,
                              one_hot = TRUE) {
  fam <- .detect_family_y(ydata, family)
  if (!fam %in% c("gaussian", "binomial", "cox")) {
    stop("family must be one of: gaussian, binomial, cox (or NULL for auto).")
  }
  if (fam == "cox" && !inherits(ydata, "Surv")) {
    stop("family='cox' requires ydata to be a survival::Surv object.")
  }

  prep <- .one_hot_xdata(xdata, one_hot = one_hot)
  build <- .build_model_frame(prep$xdata, ydata, covariates = covariates, data = data)
  df <- build$df
  predictors <- build$predictors
  covars <- build$covariates

  if (verbose) {
    message("univariableScreen: family=", fam, "; predictors=", length(predictors), "; ncores=", ncores)
  }

  fit_one <- function(p) {
    .fit_univariable_single(df, predictor = p, covariates = covars, random = random, family = fam, verbose = verbose)
  }

  if (ncores > 1L) {
    if (!requireNamespace("future", quietly = TRUE)) stop("Package 'future' is required for ncores > 1.")
    if (!requireNamespace("future.apply", quietly = TRUE)) stop("Package 'future.apply' is required for ncores > 1.")
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    future::plan(future::multisession, workers = ncores)
    res_list <- future.apply::future_lapply(predictors, fit_one, future.seed = seed)
    future::plan(future::sequential)
  } else {
    res_list <- lapply(predictors, fit_one)
  }

  res <- dplyr::bind_rows(res_list) %>%
    dplyr::mutate(
      p_adj = stats::p.adjust(p_value, method = "fdr")
    )

  summary <- res %>%
    dplyr::group_by(predictor) %>%
    dplyr::summarise(
      p_value_min = suppressWarnings(min(p_value, na.rm = TRUE)),
      p_adj_min = suppressWarnings(min(p_adj, na.rm = TRUE)),
      .groups = "drop"
    )

  out <- list(
    results = res,
    summary = summary,
    family = fam,
    call = list(ncores = ncores, seed = seed, one_hot = one_hot),
    encoded = prep$encoded
  )
  class(out) <- c("sharpener_univariable", "list")
  out
}

#' Fit a single univariable model (wrapper)
fitUnivariableModel <- function(xdata,
                                ydata,
                                predictor,
                                family = NULL,
                                covariates = NULL,
                                random = NULL,
                                data = NULL,
                                verbose = TRUE,
                                one_hot = TRUE) {
  fam <- .detect_family_y(ydata, family)
  prep <- .one_hot_xdata(xdata, one_hot = one_hot)
  build <- .build_model_frame(prep$xdata, ydata, covariates = covariates, data = data)
  if (!predictor %in% build$predictors) stop("predictor not found in xdata.")
  .fit_univariable_single(build$df, predictor = predictor, covariates = build$covariates, random = random, family = fam, verbose = verbose)
}

#' Run VariableSelection and univariableScreen on the same data
runVariableSelectionAndUnivariable <- function(xdata,
                                               ydata,
                                               family = NULL,
                                               variable_selection_args = list(),
                                               univariable_args = list(),
                                               one_hot = TRUE) {
  fam <- .detect_family_y(ydata, family)
  if (fam %in% c("coxph", "survival")) fam <- "cox"
  prep <- .one_hot_xdata(xdata, one_hot = one_hot)

  vs <- do.call(sharp::VariableSelection, c(list(xdata = prep$xdata, ydata = ydata, family = fam), variable_selection_args))
  uni <- do.call(univariableScreen, c(list(xdata = prep$xdata, ydata = ydata, family = fam, one_hot = FALSE), univariable_args))

  selprop_raw <- sharp::SelectionProportions(vs)
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
  if (is.null(names(selprop)) || all(names(selprop) == "")) {
    if (length(selprop) == ncol(prep$xdata)) {
      names(selprop) <- colnames(prep$xdata)
    }
  }

  all_names <- colnames(prep$xdata)
  if (is.null(all_names) || all(all_names == "")) {
    all_names <- names(selprop)
    if (is.null(all_names) || all(all_names == "")) {
      all_names <- paste0("V", seq_len(ncol(prep$xdata)))
    }
  }
  predictor_names <- names(selprop)
  if (is.null(predictor_names) || all(predictor_names == "")) {
    predictor_names <- all_names
    if (length(predictor_names) != length(selprop)) {
      predictor_names <- paste0("V", seq_along(selprop))
    }
  }

  selected <- .resolve_selected_vars(sharp::SelectedVariables(vs), all_names)

  comp <- tibble::tibble(
    predictor = predictor_names,
    selection_proportion = as.numeric(selprop),
    selected = predictor %in% selected
  ) %>%
    dplyr::left_join(uni$summary, by = "predictor")

  sel_threshold <- .resolve_selection_threshold(vs)
  plot_compare <- CompareVarSelUniv(
    comp,
    selection_threshold = sel_threshold,
    significance_alpha = 0.05,
    title = "Variable selection vs univariable"
  )

  list(
    variable_selection = vs,
    univariable = uni,
    comparison = comp,
    selection_threshold = sel_threshold,
    comparison_plot = plot_compare,
    encoded = prep$encoded
  )
}


#' Compare variable selection proportions vs univariable significance
CompareVarSelUniv <- function(data = NULL,
                              stability = NULL,
                              univariable = NULL,
                              pval_col = "p_value_min",
                              p_adj_col = "p_adj_min",
                              selprop_col = "selection_proportion",
                              selected_col = "selected",
                              label_col = "predictor",
                              significance_alpha = 0.05,
                              selection_threshold = NULL,
                              title = "Variable selection vs univariable",
                              subtitle = NULL,
                              text_size = 3) {
  if (!is.null(data) && inherits(data, "variable_selection")) {
    stability <- data
    data <- NULL
  }
  if (!is.null(data) && inherits(data, "sharpener_univariable")) {
    univariable <- data
    data <- NULL
  }

  if (is.list(data) && !is.null(data$comparison)) {
    if (is.null(selection_threshold) && !is.null(data$selection_threshold)) {
      selection_threshold <- data$selection_threshold
    }
    data <- data$comparison
  }

  if (is.null(data) && !is.null(stability) && !is.null(univariable)) {
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
    if (is.null(names(selprop)) || all(names(selprop) == "")) {
      if (!is.null(colnames(stability$selprop))) {
        names(selprop) <- colnames(stability$selprop)
      } else if (!is.null(dimnames(stability$Beta)) && length(dimnames(stability$Beta)) >= 2) {
        names(selprop) <- dimnames(stability$Beta)[[2]]
      }
    }

    all_names <- names(selprop)
    if (is.null(all_names) || all(all_names == "")) {
      all_names <- paste0("V", seq_along(selprop))
    }
    selected <- .resolve_selected_vars(sharp::SelectedVariables(stability), all_names)

    uni_summary <- if (is.list(univariable) && !is.null(univariable$summary)) {
      univariable$summary
    } else if (is.data.frame(univariable)) {
      univariable
    } else {
      stop("univariable must be the output from univariableScreen() or a summary data.frame.")
    }

    data <- tibble::tibble(
      predictor = all_names,
      selection_proportion = as.numeric(selprop),
      selected = predictor %in% selected
    ) %>%
      dplyr::left_join(uni_summary, by = "predictor")

    if (is.null(selection_threshold)) {
      selection_threshold <- .resolve_selection_threshold(stability)
    }
  }

  if (is.null(data)) {
    stop("Provide either a prepared data.frame, or both stability and univariable inputs.")
  }

  required_cols <- c(selprop_col, pval_col, label_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols)) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

  if (is.null(selection_threshold) || !is.finite(selection_threshold)) selection_threshold <- NA_real_

  if (p_adj_col %in% names(data)) {
    p_adj_vec <- data[[p_adj_col]]
  } else {
    p_adj_vec <- stats::p.adjust(data[[pval_col]], method = "fdr")
  }
  selected_exists <- selected_col %in% names(data)

  df <- data %>%
    dplyr::mutate(
      logp = -log10(.data[[pval_col]]),
      logp = dplyr::if_else(is.finite(logp), logp, NA_real_),
      p_adj_calc = p_adj_vec,
      sig_category = dplyr::case_when(
        !is.na(.data[[pval_col]]) & .data[[pval_col]] < (significance_alpha / n()) ~ "FWER<0.05",
        !is.na(p_adj_calc) & p_adj_calc < significance_alpha ~ "FDR<0.05",
        !is.na(.data[[pval_col]]) & .data[[pval_col]] < significance_alpha ~ "p<0.05",
        TRUE ~ "NS"
      ),
      sig_category = factor(sig_category, levels = c("FWER<0.05", "FDR<0.05", "p<0.05", "NS")),
      label_txt = .data[[label_col]],
      selprop_plot = .data[[selprop_col]],
      selected_flag = if (selected_exists) .data[[selected_col]] else FALSE
    )

  df$logp[is.na(df$logp)] <- 0

  df_label <- df %>%
    dplyr::filter(selected_flag | (!is.na(p_adj_calc) & p_adj_calc < significance_alpha))

  n_pred <- nrow(df)
  fwer_thr <- if (n_pred > 0L) significance_alpha / n_pred else NA_real_
  p_thr <- significance_alpha
  fdr_thr <- significance_alpha

  p <- ggplot2::ggplot(df, ggplot2::aes(x = selprop_plot, y = logp, color = sig_category)) +
    ggplot2::geom_point(alpha = 0.8, size = 1.6) +
    ggplot2::geom_point(
      data = df %>% dplyr::filter(selected_flag),
      ggplot2::aes(x = selprop_plot, y = logp),
      inherit.aes = FALSE,
      shape = 21,
      size = 2.4,
      stroke = 0.6,
      colour = "black",
      fill = NA_real_
    ) +
    ggplot2::geom_vline(
      xintercept = selection_threshold,
      linetype = "dotted",
      linewidth = 0.3,
      color = "firebrick3",
      na.rm = TRUE
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(fwer_thr),
      linetype = "dotted",
      linewidth = 0.3,
      color = "#E64B35",
      na.rm = TRUE
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(fdr_thr),
      linetype = "dotted",
      linewidth = 0.3,
      color = "orange2",
      na.rm = TRUE
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(p_thr),
      linetype = "dotted",
      linewidth = 0.3,
      color = "#3C5488",
      na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(values = c(
      "NS" = "grey50",
      "FWER<0.05" = "#E64B35",
      "FDR<0.05" = "orange2",
      "p<0.05" = "#3C5488"
    ), limits = c("FWER<0.05", "FDR<0.05", "p<0.05", "NS"), drop = FALSE) +
    ggplot2::labs(
      x = "Selection proportion (LASSO stability)",
      y = "-log10(p) (univariable)",
      color = NULL,
      title = title,
      subtitle = subtitle
    ) +
    theme_react_safe()

  if (nrow(df_label) > 0) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = df_label,
        ggplot2::aes(x = selprop_plot, y = logp, label = label_txt),
        size = text_size,
        show.legend = FALSE
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = df_label,
        ggplot2::aes(x = selprop_plot, y = logp, label = label_txt),
        size = text_size,
        vjust = -0.5,
        show.legend = FALSE
      )
    }
  }

  p
}


# Volcano plot ------------------------------------------------------------

plotVolcano <- function(data,
                        estimate_col = "estimate",
                        pval_col = "p_value",
                        p_adj_col = "p_adj",
                        label_col = "predictor",
                        significance_alpha = 0.05,
                        top_n_to_label = 10,
                        title = NULL,
                        text_size = 3) {
  if (is.list(data) && !is.null(data$results)) data <- data$results

  required_cols <- c(estimate_col, pval_col, label_col)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols)) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

  if (p_adj_col %in% names(data)) {
    p_adj_vec <- data[[p_adj_col]]
  } else {
    p_adj_vec <- stats::p.adjust(data[[pval_col]], method = "fdr")
  }

  df <- data %>%
    dplyr::mutate(
      logp = -log10(.data[[pval_col]]),
      logp = dplyr::if_else(is.finite(logp), logp, NA_real_),
      p_adj_calc = p_adj_vec,
      sig_category = dplyr::case_when(
        !is.na(.data[[pval_col]]) & .data[[pval_col]] < (significance_alpha / n()) ~ "FWER<0.05",
        !is.na(p_adj_calc) & p_adj_calc < significance_alpha ~ "FDR<0.05",
        !is.na(.data[[pval_col]]) & .data[[pval_col]] < significance_alpha ~ "p<0.05",
        TRUE ~ "NS"
      ),
      sig_category = factor(sig_category, levels = c("FWER<0.05", "FDR<0.05", "p<0.05", "NS")),
      label_txt = .data[[label_col]],
      estimate_plot = .data[[estimate_col]]
    )

  df$logp[is.na(df$logp)] <- 0

  df_to_label <- if (top_n_to_label > 0) {
    df %>% dplyr::arrange(dplyr::desc(logp)) %>% dplyr::slice(seq_len(min(top_n_to_label, nrow(df))))
  } else {
    df[FALSE, ]
  }

  legend_levels <- c("FWER<0.05", "FDR<0.05", "p<0.05", "NS")
  legend_df <- data.frame(
    estimate_plot = 0,
    logp = 0,
    sig_category = factor(legend_levels, levels = legend_levels)
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = estimate_plot, y = logp, color = sig_category)) +
    ggplot2::geom_point(alpha = 0.75, size = 1) +
    ggplot2::geom_point(
      data = legend_df,
      ggplot2::aes(x = estimate_plot, y = logp, color = sig_category),
      inherit.aes = FALSE,
      alpha = 0,
      size = 0,
      show.legend = TRUE
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3, color = "grey40") +
    ggplot2::geom_hline(yintercept = -log10(significance_alpha), linetype = "dotted", linewidth = 0.3, color = "grey40") +
    ggplot2::scale_color_manual(values = c(
      "NS" = "grey50",
      "FWER<0.05" = "#E64B35",
      "FDR<0.05" = "orange2",
      "p<0.05" = "#3C5488"
    ), limits = legend_levels, drop = FALSE) +
    ggplot2::labs(
      x = "Estimate",
      y = "-log10(p)",
      color = NULL,
      title = title
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1, size = 2))) +
    theme_react_safe()

  if (nrow(df_to_label) > 0) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = df_to_label,
        ggplot2::aes(label = label_txt),
        size = text_size,
        show.legend = FALSE
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = df_to_label,
        ggplot2::aes(label = label_txt),
        size = text_size,
        vjust = -0.5,
        show.legend = FALSE
      )
    }
  }

  p
}
