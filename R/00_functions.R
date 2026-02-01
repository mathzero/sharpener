

`%||%` <- function(a, b) if (!is.null(a)) a else b

.has_overreact <- requireNamespace("OverReact", quietly = TRUE)

mycols=c("#000080", "#00bfff", "#fa8072" , "#232333","#7b68ee"   ,"#006400",   "#c71585")

theme_react_safe <- function(...) {
  if (.has_overreact) OverReact::theme_react(...) else ggplot2::theme_minimal()
}

scale_colour_imperial_safe <- function() {
  # if (.has_overreact) OverReact::scale_color_imperial() else ggplot2::scale_colour_discrete()
  scale_colour_manual(values = mycols)
}

scale_fill_imperial_safe <- function() {
  # if (.has_overreact) OverReact::scale_fill_imperial() else ggplot2::scale_fill_discrete()
  scale_fill_manual(values = mycols)
  
}
# Utility functions -------------------------------------------------------

# Mean-impute missing values column-wise while preserving common types
mean_impute_df <- function(df) {
  stopifnot(is.data.frame(df))
  
  out <- df
  
  for (nm in names(out)) {
    x <- out[[nm]]
    
    # Handle Date columns by imputing with the mean calendar day
    if (inherits(x, "Date")) {
      x_num <- as.numeric(x)
      m <- mean(x_num, na.rm = TRUE)
      if (!is.na(m)) {
        m_date <- as.Date(round(m), origin = "1970-01-01")
        nas <- is.na(x)
        x[nas] <- m_date
        out[[nm]] <- x
      }
      next
    }
    
    # Handle POSIXct/POSIXlt by imputing with the mean timestamp
    if (inherits(x, c("POSIXct", "POSIXlt"))) {
      tz <- attr(x, "tzone")
      if (is.null(tz)) tz <- ""
      if (length(tz) > 0) tz <- tz[1]
      
      # Work in POSIXct for stability, then restore class if needed
      x_ct <- as.POSIXct(x, tz = tz)
      x_num <- as.numeric(x_ct)
      m <- mean(x_num, na.rm = TRUE)
      if (!is.na(m)) {
        m_ct <- as.POSIXct(m, origin = "1970-01-01", tz = tz)
        nas <- is.na(x_ct)
        x_ct[nas] <- m_ct
        # Restore original class
        if (inherits(x, "POSIXlt")) {
          out[[nm]] <- as.POSIXlt(x_ct, tz = tz)
        } else {
          out[[nm]] <- x_ct
        }
      }
      next
    }
    
    # Handle standard numeric columns (double or integer)
    if (is.numeric(x)) {
      m <- mean(x, na.rm = TRUE)
      if (!is.na(m)) {
        nas <- is.na(x)
        if (is.integer(x)) {
          # Preserve integer storage by rounding the mean
          m_int <- as.integer(round(m))
          x[nas] <- m_int
        } else {
          x[nas] <- m
        }
        out[[nm]] <- x
      }
      next
    }
    
    # Leave non-numeric/factor/character/logical columns unchanged
    # Potential future enhancement: add mode-imputation for factors/characters
  }
  
  out
}




plotSurvivalCurve <- function(
    mydata,
    rev_colours=F,
    testname   = "Lesion Severity Score",
    strat_col  = NULL,
    lesion_col = "lesion",
    id_col     = "subject_identifier_for_the_study",
    time_col   = "time",
    score_col  = "numeric_result_finding_in_standard_units",
    test_col   = "evaluation_test_name") {
  
  # ── Safety checks ------------------------------------------------------------
  if (!is.null(strat_col) && !is.character(strat_col)) {
    stop("`strat_col` must be a character string")
  }
  if (!is.null(strat_col) && !(strat_col %in% names(mydata))) {
    stop("Column '", strat_col, "' not found in `mydata`")
  }
  
  # ── 1) Keep rows for the target test ----------------------------------------
  event_df <- mydata %>%
    filter(.data[[test_col]] == testname) %>%
    mutate(
      lesion_uid = paste(.data[[id_col]], .data[[lesion_col]], sep = "_"),
      event_flag = .data[[score_col]] > 2 & .data[[time_col]] > 0
    )
  
  # ── 2) Lesion-level time & status (+ optional stratum) ----------------------
  surv_data <- event_df %>%
    arrange(lesion_uid, .data[[time_col]]) %>%
    group_by(lesion_uid) %>%
    summarise(
      time   = if (any(event_flag)) min(.data[[time_col]][event_flag]) else max(.data[[time_col]]),
      status = as.integer(any(event_flag)),
      across(all_of(strat_col), ~ dplyr::first(.x)),
      .groups = "drop"
    )
  
  # ── 3) Add an explicit grouping column for stratified fits ------------------
  if (is.null(strat_col)) {
    surv_data <- surv_data %>% mutate(.stratum = factor("All"))
    fit <- survfit(Surv(time, status) ~ 1, data = surv_data)
  } else {
    surv_data <- surv_data %>%
      mutate(.stratum = factor(.data[[strat_col]])) %>%
      droplevels()
    fit <- survfit(Surv(time, status) ~ .stratum, data = surv_data)
    if (!is.null(fit$strata)) {
      names(fit$strata) <- sub("^[^.]*\\.", "", names(fit$strata))
    }
  }
  
  # ── 4) Plot with survminer (pass the data explicitly) -----------------------
  # mycols=wesanderson::wes_palette(name = "Darjeeling1",n = 5)
  
  mycols= c("#3B3B3B", "#CD534C", "#7AA6DC", "#EFC000"  ,"#008280"  , "#631879")
  if(rev_colours){
    mycols=rev(mycols)
  }
  
  gp <- ggsurvplot(
    fit        = fit,
    data       = surv_data,
    palette = as.character(mycols),
    risk.table = TRUE,
    conf.int   = TRUE,
    ggtheme    = OverReact::theme_mw(),
    tables.theme = OverReact::theme_mw(),
    xlab       = "Follow-up time",
    ylab       = "Probability lesion severity ≤ 2",
    title      = "Time to Lesion-Severity > 2"
  )
  
  return((gp$plot & theme(legend.position = c(1, 1),
                          legend.justification = c(1, 1))) / 
           gp$table + patchwork::plot_layout(heights = c(8, 2)))
}


# km_plot() ── minimal-dependency Kaplan–Meier plot
# ------------------------------------------------------------------------------
# • Accepts either:
#     – a pre-computed survfit object,   or
#     – (formula, data) so the function can call survfit() for you.
# • Handles optional stratification (any number of factor levels).
# • Draws step curves with shaded 95 % CIs.
# • Returns a ggplot object you can print or further customise.
# ==============================================================================

km_plot <- function(fit       = NULL,
                    formula   = NULL,
                    data      = NULL,
                    conf_int  = TRUE,
                    palette   = NULL,
                    xlab      = "Time",
                    ylab      = "Survival probability",
                    title     = NULL) {
  
  stopifnot(
    xor(is.null(fit), is.null(formula)),   # must supply exactly one of them
    requireNamespace("survival", quietly = TRUE),
    requireNamespace("ggplot2",  quietly = TRUE)
  )
  
  library(survival)
  library(ggplot2)
  
  # ---------------------------------------------------------------------------
  # 1. Obtain a survfit object -------------------------------------------------
  # ---------------------------------------------------------------------------
  if (is.null(fit)) {
    fit <- survfit(formula, data = data)
  }
  
  # ---------------------------------------------------------------------------
  # 2. Tidy the survfit object -------------------------------------------------
  # ---------------------------------------------------------------------------
  # Convert to one data frame row per step per stratum
  tidy_survfit <- function(sf) {
    strata <- sf$strata
    if (is.null(strata)){
      strata <- ""
      strata_names <- ""
    }else{
      strata_names <- rep(names(strata), strata)
    }
    
    df <- data.frame(
      strata = factor(strata_names, levels = unique(strata_names)),
      time   = sf$time,
      surv   = sf$surv,
      upper  = sf$upper,
      lower  = sf$lower
    )
    
    # prepend time = 0  (survival = 1) for each stratum
    zeros <- df |> 
      dplyr::distinct(strata) |>
      dplyr::mutate(time = 0, surv = 1, upper = 1, lower = 1)
    
    dplyr::bind_rows(zeros, df) |>
      dplyr::arrange(strata, time)
  }
  
  df_plot <- tidy_survfit(fit)

  # ---------------------------------------------------------------------------
  # 3. Build ggplot ------------------------------------------------------------
  # ---------------------------------------------------------------------------
  n_strata <- length(levels(df_plot$strata))
  if (is.null(palette)) {
    palette <- scales::hue_pal()(max(1, n_strata))
  }
  
  p <- ggplot(df_plot, aes(time, surv, colour = strata)) +
    geom_step(size = 0.8) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    scale_colour_manual(values = palette, name = NULL) +
    labs(x = xlab, y = ylab, title = title) +
    OverReact::theme_mw() +
    ggsci::scale_color_jco()+
    ggsci::scale_fill_jco()+
    theme(legend.position = if (n_strata > 1) "right" else "none")
  
  if (conf_int) {
    p <- p +
      geom_ribbon(
        aes(ymin = lower, ymax = upper, fill = strata),
        alpha = 0.2,
        colour = NA
      ) +
      ggsci::scale_color_jco()+
      ggsci::scale_fill_jco() }
  
  return(p)
}




# Plot clinical descriptors -----------------------------------------------


plotBaselineCharacteristics <- function(data, cat="At what approximate age did atopic dermatitis begin?") {
  library(scales)
  data %>% 
    filter(
      subject_identifier_for_the_study%in%active_ids,
      evaluation_test_name==cat) %>%
    mutate(
      label = fct_reorder(character_result_finding_in_std_format,          # character to be reordered
                          numeric_result_finding_in_standard_units,           # numeric that carries the order
                          .fun  = min),   # **one number per group**
      label = factor(label, ordered = TRUE)
    ) %>% 
    ggplot(aes(y=label))+
    geom_histogram(stat="count") +
    OverReact::scale_color_imperial()+
    OverReact::theme_react() +
    labs(title=cat, y="")+
    scale_x_continuous(breaks = breaks_pretty())+
    theme(axis.text.x = element_text(size=15),
          axis.text.y = element_text(size=15),
          axis.title.x =  element_text(size=15))
}


plotScoresOverTime <- function(data, cat='Quantification of Staphylococcus aureus (log 10)') {
  library(scales)
   plotdat <- data %>% 
     filter(subject_identifier_for_the_study%in%active_ids,
            evaluation_test_name==cat) %>% 
     group_by(lesion,time,category_for_evaluation) %>%
     summarise(numeric_result_finding_in_standard_units=mean(numeric_result_finding_in_standard_units)) %>%
     left_join(df_relapse)

   plotdat |> 
     ggplot()+
    geom_point(aes(x=time,y=numeric_result_finding_in_standard_units,col=category_for_evaluation)) +
    geom_path(aes(x=time,y=numeric_result_finding_in_standard_units,col=category_for_evaluation)) +
    geom_point(data = plotdat |> group_by(lesion) |>arrange(time) |> slice_tail(n=1) |> ungroup(),
               aes(y=numeric_result_finding_in_standard_units,x=time_of_relapse),shape=1,size=3, col="red")+
    OverReact::scale_color_imperial()+
     scale_x_continuous(breaks = breaks_pretty())+
     
    facet_wrap(.~lesion)+
    OverReact::theme_react() +
    labs(title=cat, col="Evaluation \nmethod")
}



plotScoresOverTimeCompare <- function(data, 
                                      column_with_cats="evaluation_test_name",
                                      cat='Quantification of Staphylococcus aureus (log 10)') {
  library(scales)
  plotdat <- data %>% 
    filter(subject_identifier_for_the_study%in%active_ids,
           !!sym(column_with_cats)==cat) %>% 
    group_by(lesion,time,category_for_evaluation) %>%
    summarise(numeric_result_finding_in_standard_units=mean(numeric_result_finding_in_standard_units)) %>%
    left_join(df_relapse)
  
  plotdat |> 
    ggplot()+
    geom_point(aes(x=time,y=numeric_result_finding_in_standard_units,col=category_for_evaluation)) +
    geom_path(aes(x=time,y=numeric_result_finding_in_standard_units,col=category_for_evaluation)) +
    geom_point(data = plotdat |> group_by(lesion) |>arrange(time) |> slice_tail(n=1) |> ungroup(),
               aes(y=numeric_result_finding_in_standard_units,x=time_of_relapse),shape=1,size=3, col="red")+
    OverReact::scale_color_imperial()+
    scale_x_continuous(breaks = breaks_pretty())+
    
    facet_wrap(.~lesion)+
    OverReact::theme_react() +
    labs(title=cat, col="Evaluation \nmethod")
}




# Cox modelling functions -------------------------------------------------

# Optional aesthetics used in your KM code; these lines are safe if packages are absent
`%||%` <- function(x, y) if (is.null(x)) y else x

# ── 1. Helpers ────────────────────────────────────────────────────────────────

# First non-missing extractor for group summaries
first_non_missing <- function(x) {
  y <- x[!is.na(x)]
  if (length(y)) y[1] else NA
}

# Build lesion-level survival dataset using the SAME time/status definitions as KM
make_surv_lesion_data <- function(
    mydata,
    testname  = "Lesion Severity Score",
    lesion_col = "lesion",
    id_col     = "subject_identifier_for_the_study",
    time_col   = "time",
    score_col  = "numeric_result_finding_in_standard_units",
    test_col   = "evaluation_test_name",
    event_threshold = 2,
    event_time_min  = 0,
    extra_cols = NULL  # predictors/strata/cluster columns to carry to lesion level
) {
  # We construct a unique lesion id and compute event flags from the master test
  data_with_uid <- mydata %>%
    mutate(lesion_uid = paste(.data[[id_col]], .data[[lesion_col]], sep = "_"))
  
  event_df <- data_with_uid %>%
    filter(.data[[test_col]] == testname) %>%
    mutate(event_flag = .data[[score_col]] > event_threshold & .data[[time_col]] > event_time_min)
  
  # Lesion-level time-to-event summary (same logic as your KM function)
  surv_core <- event_df %>%
    arrange(lesion_uid, .data[[time_col]]) %>%
    group_by(lesion_uid) %>%
    summarise(
      time   = if (any(event_flag)) min(.data[[time_col]][event_flag]) else max(.data[[time_col]]),
      status = as.integer(any(event_flag)),
      # also keep the subject id at lesion level for clustering if needed
      !!id_col := first(.data[[id_col]]),
      .groups = "drop"
    )
  
  # Bring forward predictors/strata/cluster columns as first non-missing per lesion
  if (!is.null(extra_cols) && length(extra_cols)) {
    missing_cols <- setdiff(extra_cols, names(data_with_uid))
    if (length(missing_cols)) {
      stop("The following columns requested in `extra_cols` are not in `mydata`: ",
           paste(missing_cols, collapse = ", "))
    }
    
    extras <- data_with_uid %>%
      group_by(lesion_uid) %>%
      summarise(across(all_of(unique(c(extra_cols, id_col))), first_non_missing, .names = "{.col}"),
                .groups = "drop")
    
    surv_df <- surv_core %>% left_join(extras, by = "lesion_uid")
  } else {
    surv_df <- surv_core
  }
  
  surv_df
}


# ── Helpers shared by Cox/GLM forests ─────────────────────────────────────────

# Extract HR table from a coxph fit (Wald-style)
# ── Tables (unchanged) ────────────────────────────────────────────────────────
# Replace your hr_table() with this version
hr_table <- function(fit, conf_level = 0.95) {
  s  <- summary(fit)
  co <- s$coefficients
  ci <- s$conf.int
  
  # columns present in summary.coxph across standard/robust fits
  coef_col <- "coef"
  se_col   <- if ("robust se" %in% colnames(co)) "robust se" else "se(coef)"
  z_val    <- if ("z" %in% colnames(co)) co[, "z"] else co[, coef_col] / co[, se_col]
  pval     <- 2 * pnorm(abs(z_val), lower.tail = FALSE)
  
  lower_name <- grep("^lower", colnames(ci),  value = TRUE)[1]
  upper_name <- grep("^upper", colnames(ci),  value = TRUE)[1]
  hr_col     <- if ("exp(coef)" %in% colnames(ci)) "exp(coef)" else NULL
  hr         <- if (!is.null(hr_col)) ci[, hr_col] else exp(co[, coef_col])
  
  tibble::tibble(
    term     = gsub("`", "", rownames(co)),
    estimate = as.numeric(hr),
    lower    = as.numeric(ci[, lower_name]),
    upper    = as.numeric(ci[, upper_name]),
    pval     = as.numeric(pval),
    measure  = "HR"
  )
}

# (Optional) tighten up the GLM extractor too, to avoid regex surprises
or_table <- function(fit, conf_level = 0.95) {
  sm  <- summary(fit)
  co  <- sm$coefficients
  est <- co[, "Estimate"]
  se  <- co[, "Std. Error"]
  z   <- if ("z value" %in% colnames(co)) co[, "z value"] else est / se
  p   <- 2 * pnorm(abs(z), lower.tail = FALSE)
  
  zcrit <- qnorm(1 - (1 - conf_level)/2)
  tibble::tibble(
    term     = gsub("`", "", rownames(co)),
    estimate = exp(est),
    lower    = exp(est - zcrit * se),
    upper    = exp(est + zcrit * se),
    pval     = as.numeric(p),
    measure  = "OR"
  ) |>
    dplyr::filter(.data$term != "(Intercept)")
}


# ── Core ordering/plotting (FIXED) ────────────────────────────────────────────
.clean_label <- function(x) {
  x %>%
    gsub("(^|\\s)([A-Za-z0-9_.]+)\\((.+)\\)", "\\2: \\3", ., perl = TRUE) %>%
    gsub(":", " = ", .) %>%
    stringr::str_replace_all("_", " ")
}

.compute_and_order <- function(df,
                               order_by = c("none","pval","effect","abs_effect"),
                               var_order = NULL,          # RAW term names
                               label_wrap = 28) {
  order_by <- match.arg(order_by)
  
  # keep a clean label but DO NOT use it for matching
  df <- df %>%
    dplyr::mutate(
      term_clean = .clean_label(term),
      log_est    = log(estimate),
      effect_mag = abs(log_est)
    )
  
  # remove unusable rows (avoids NA/Inf in plot ranges)
  df <- df %>% dplyr::filter(is.finite(estimate), is.finite(lower), is.finite(upper))
  
  if (!is.null(var_order)) {
    # match by raw terms; append any extras not in var_order
    extras <- setdiff(df$term, var_order)
    lvls_raw <- c(var_order, extras)
    df <- df %>%
      dplyr::mutate(term = factor(term, levels = lvls_raw)) %>%
      dplyr::arrange(term) %>%
      dplyr::mutate(term_clean = factor(stringr::str_wrap(term_clean, width = label_wrap),
                                        levels = stringr::str_wrap(.clean_label(as.character(lvls_raw)),
                                                                   width = label_wrap)))
  } else {
    df <- switch(order_by,
                 pval       = dplyr::arrange(df, pval),
                 effect     = dplyr::arrange(df, log_est),
                 abs_effect = dplyr::arrange(df, dplyr::desc(effect_mag)),
                 none       = df
    )
    df <- df %>%
      dplyr::mutate(term_clean = factor(stringr::str_wrap(term_clean, width = label_wrap),
                                        levels = stringr::str_wrap(term_clean, width = label_wrap)))
  }
  
  df
}

.plot_forest_core <- function(df, ref_line = 1, point_size = 2.6, line_size = 0.8,
                              label_wrap = 28, title = "Effects", subtitle = NULL, caption = NULL,
                              x_lab = "Effect (log scale, 95% CI)") {
  ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = term_clean)) +
    ggplot2::geom_vline(xintercept = ref_line, linetype = 2) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lower, xmax = upper), height = 0, size = line_size) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    ggplot2::labs(x = x_lab, y = NULL, title = title, subtitle = subtitle, caption = caption) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor   = ggplot2::element_blank(),
                   axis.text.y        = ggplot2::element_text(hjust = 1)) +
    ggplot2::geom_text(
      ggplot2::aes(
        # x = max(upper, na.rm = TRUE) * 1.15,
                   x = upper* 1.15,
                   label = paste0("p = ", round(pval, 4))),
      hjust = 0, size = 3
    ) +
    ggplot2::coord_cartesian(xlim = c(min(df$lower, na.rm = TRUE) / 1.5,
                                      max(df$upper, na.rm = TRUE) * 1.5),
                             clip = "off") +
    ggplot2::theme(plot.margin = ggplot2::margin(5.5, 60, 5.5, 5.5, "pt")) +
    {if (requireNamespace("OverReact", quietly = TRUE)) OverReact::theme_mw() else ggplot2::theme()} +
    {if (requireNamespace("ggsci", quietly = TRUE)) ggsci::scale_color_jco() else ggplot2::scale_colour_discrete()} +
    {if (requireNamespace("ggsci", quietly = TRUE)) ggsci::scale_fill_jco() else ggplot2::scale_fill_discrete()}
}

# ── Public API (updated) ──────────────────────────────────────────────────────
plot_cox_forest <- function(fit_or_hr,
                            order_by = c("none","pval","effect","abs_effect"),
                            ref_line = 1, point_size = 2.6, line_size = 0.8,
                            label_wrap = 28,
                            title = "Hazard ratios", subtitle = NULL, caption = NULL,
                            var_order = NULL) {       # RAW term vector
  if (inherits(fit_or_hr, "coxph")) {
    df <- hr_table(fit_or_hr)
  } else {
    df <- fit_or_hr
    if (!("estimate" %in% names(df)) && "HR" %in% names(df)) df$estimate <- df$HR
    df$measure <- df$measure %||% "HR"
  }
  df <- .compute_and_order(df, order_by = match.arg(order_by),
                           var_order = var_order, label_wrap = label_wrap)
  .plot_forest_core(df, ref_line, point_size, line_size, label_wrap,
                    title, subtitle, caption,
                    x_lab = "Hazard ratio (log scale, 95% CI)")
}

plot_glm_forest <- function(fit_or_or,
                            order_by = c("none","pval","effect","abs_effect"),
                            ref_line = 1, point_size = 2.6, line_size = 0.8,
                            label_wrap = 28,
                            title = "Odds ratios", subtitle = NULL, caption = NULL,
                            var_order = NULL) {      # RAW term vector
  if (inherits(fit_or_or, "glm")) {
    if (!(family(fit_or_or)$family == "binomial" && family(fit_or_or)$link == "logit"))
      stop("plot_glm_forest expects glm(binomial, link='logit').")
    df <- or_table(fit_or_or)
  } else {
    df <- fit_or_or
    if (!("estimate" %in% names(df)) && "OR" %in% names(df)) df$estimate <- df$OR
    df$measure <- df$measure %||% "OR"
  }
  df <- .compute_and_order(df, order_by = match.arg(order_by),
                           var_order = var_order, label_wrap = label_wrap)
  .plot_forest_core(df, ref_line, point_size, line_size, label_wrap,
                    title, subtitle, caption,
                    x_lab = "Odds ratio (log scale, 95% CI)")
}

# Harmonised ordering helper: now returns RAW terms, not wrapped labels
get_forest_var_order <- function(x, order_by = c("none","pval","effect","abs_effect")) {
  if (inherits(x, "coxph")) {
    df <- hr_table(x)
  } else if (inherits(x, "glm") && family(x)$family == "binomial") {
    df <- or_table(x)
  } else if (is.data.frame(x)) {
    df <- x
    if (!("estimate" %in% names(df))) {
      if ("HR" %in% names(df)) df$estimate <- df$HR
      if ("OR" %in% names(df)) df$estimate <- df$OR
    }
  } else stop("Unsupported object for get_forest_var_order().")
  
  order_by <- match.arg(order_by)
  df <- df %>%
    dplyr::mutate(log_est = log(estimate),
                  effect_mag = abs(log_est)) %>%
    dplyr::filter(is.finite(estimate))
  
  df <- switch(order_by,
               pval       = dplyr::arrange(df, pval),
               effect     = dplyr::arrange(df, log_est),
               abs_effect = dplyr::arrange(df, dplyr::desc(effect_mag)),
               none       = df
  )
  df$term
}





fit_cox_model <- function(
    mydata,
    testname  = "Lesion Severity Score",
    predictors,                # character vector of RHS terms, e.g., c("age", "sex", "stage")
    lesion_col = "lesion",
    id_col     = "subject_identifier_for_the_study",
    time_col   = "time",
    score_col  = "numeric_result_finding_in_standard_units",
    test_col   = "evaluation_test_name",
    strata_vars = NULL,        # optional character vector of strata variables
    cluster_id_col = NULL,     # optional clustering variable for robust SE (e.g., subject id)
    ties = "efron",
    na_action = na.omit,       # we default to complete-case for model fitting
    event_threshold = 2,
    event_time_min  = 0
) {
  stopifnot(is.character(predictors), length(predictors) >= 1)
  
  # We need to carry any variables referenced by the formula to lesion level
  needed_cols <- unique(c(predictors, strata_vars %||% character(), cluster_id_col %||% character()))
  surv_df <- make_surv_lesion_data(
    mydata         = mydata,
    testname       = testname,
    lesion_col     = lesion_col,
    id_col         = id_col,
    time_col       = time_col,
    score_col      = score_col,
    test_col       = test_col,
    event_threshold = event_threshold,
    event_time_min  = event_time_min,
    extra_cols     = needed_cols
  )
  
  # Build formula: Surv(time, status) ~ predictors [+ strata(...)] [+ cluster(...)]
  rhs <- paste(predictors, collapse = " + ")
  
  if (!is.null(strata_vars) && length(strata_vars)) {
    # allow multiple strata variables; survival supports strata(a, b) or strata(a)+strata(b)
    strata_term <- paste0("strata(", paste(strata_vars, collapse = ", "), ")")
    rhs <- paste(rhs, strata_term, sep = " + ")
  }
  
  if (!is.null(cluster_id_col) && nzchar(cluster_id_col)) {
    rhs <- paste(rhs, sprintf("cluster(%s)", cluster_id_col), sep = " + ")
  }
  
  fml <- as.formula(paste0("Surv(time, status) ~ ", rhs))
  
  # Fit Cox model
  fit <- coxph(
    formula = fml,
    data    = surv_df,
    ties    = ties,
    na.action = na_action,
    model   = TRUE, x = TRUE, y = TRUE
  )
  
  # Schoenfeld PH test (useful to inspect)
  zph <- try(cox.zph(fit), silent = TRUE)
  
  # HR table
  hr <- hr_table(fit)
  
  list(
    fit     = fit,
    surv_df = surv_df,
    hr      = hr,
    zph     = zph,
    formula = fml
  )
}




# IncrementalAnalysis -----------------------------------------------------

IncrementalAnalysis <- function(xdata, ydata, n_predictors = NULL, stability = NULL,
                                test_ratio = 0.3, k = 100, seed = 123, ncores = 4, verbose = TRUE,
                                ci_level = 0.95, roc_points = 101, max_predictors = 40) {
  if (!requireNamespace("future", quietly = TRUE)) stop("Package 'future' is required.")
  if (!requireNamespace("future.apply", quietly = TRUE)) stop("Package 'future.apply' is required.")
  
  set.seed(seed)
  
  xdf <- as.data.frame(xdata)
  if (is.null(n_predictors)) n_predictors <- min(ncol(xdf), max_predictors)
  
  family <- NULL
  if (inherits(ydata, "Surv")) {
    family <- "cox"
  } else if (is.factor(ydata) || length(unique(na.omit(ydata))) == 2) {
    family <- "binomial"
  } else {
    family <- "gaussian"
  }
  
  metric <- if (family == "binomial") "AUC" else if (family == "cox") "C-index" else "R-squared"
  
  if (family == "binomial") {
    if (!requireNamespace("pROC", quietly = TRUE)) stop("Package 'pROC' is required for binomial outcomes.")
    if (!requireNamespace("caret", quietly = TRUE)) stop("Package 'caret' is required for binomial outcomes.")
    
    target_vec <- ydata
    if (is.factor(target_vec)) {
      target_vec <- as.integer(target_vec == levels(target_vec)[2])
    } else if (is.logical(target_vec)) {
      target_vec <- as.integer(target_vec)
    } else if (is.numeric(target_vec)) {
      u <- sort(unique(na.omit(as.numeric(target_vec))))
      if (identical(u, c(0, 1))) target_vec <- as.integer(target_vec)
      else if (identical(u, c(1, 2))) target_vec <- as.integer(target_vec - 1L)
      else stop("Binary outcome must be factor/logical/0-1/1-2.")
    } else stop("Unsupported ydata type for binomial.")
  }
  
  if (family == "gaussian") {
    target_vec <- as.numeric(ydata)
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
    if (is.na(v[1])) v[1] <- v[which(!is.na(v))[1]]
    if (is.na(v[length(v)])) v[length(v)] <- v[rev(which(!is.na(v)))[1]]
    na_idx <- which(is.na(v))
    if (length(na_idx) > 0) {
      good <- !is.na(v)
      v[!good] <- approx(x = which(good), y = v[good], xout = which(!good),
                         method = "linear", rule = 2, ties = "ordered")$y
    }
    v
  }
  
  safe_cindex <- function(df) {
    out <- try(
      survival::concordance(survival::Surv(time, event) ~ lp, data = df, reverse = TRUE)$concordance,
      silent = TRUE
    )
    if (inherits(out, "try-error")) return(NA_real_)
    as.numeric(out)
  }
  
  # Precompute k splits once, reuse across all model sizes
  set.seed(seed)
  split_list <- vector("list", k)
  for (i in seq_len(k)) {
    if (family == "binomial") {
      tr_idx <- caret::createDataPartition(factor(data$target), p = 1 - test_ratio, list = FALSE)
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
        roc_curve_test <- pROC::roc(response = test_data$target, predictor = test_prob, quiet = TRUE, direction = "auto")
        metric_result_test <- as.numeric(pROC::auc(roc_curve_test))
        sens_vec_test <- suppressWarnings(
          pROC::coords(roc_curve_test, x = spec_grid, input = "specificity", ret = "sensitivity", transpose = FALSE)
        )
        sens_vec_test <- fill_na_linear(as.numeric(unlist(sens_vec_test, use.names = FALSE)))
        
        train_prob <- as.numeric(stats::predict(fit, newdata = train_data, type = "response"))
        roc_curve_train <- pROC::roc(response = train_data$target, predictor = train_prob, quiet = TRUE, direction = "auto")
        metric_result_train <- as.numeric(pROC::auc(roc_curve_train))
        sens_vec_train <- suppressWarnings(
          pROC::coords(roc_curve_train, x = spec_grid, input = "specificity", ret = "sensitivity", transpose = FALSE)
        )
        sens_vec_train <- fill_na_linear(as.numeric(unlist(sens_vec_train, use.names = FALSE)))
        
        return(list(beta_coefficients_names = names(beta_coefficients),
                    beta_coefficients = as.vector(beta_coefficients),
                    metric_result = metric_result_test,
                    metric_result_train = metric_result_train,
                    tpr = sens_vec_test,
                    tpr_train = sens_vec_train))
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
  
  added_vars <- vapply(details, function(d) tail(d$variable_names, 1), character(1))
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
  
  call <- list(test_ratio = test_ratio, k = k, metric = metric, ci_level = ci_level, roc_points = roc_points)
  res <- list(summary = summary_output, roc_ci = roc_ci, roc_ci_train = roc_ci_train, call = call, details = details)
  
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






# Predicted vs actual survival curve --------------------------------------

plot_predicted_vs_actual_survival <- function(incremental_res,
                                              xdata,
                                              ydata,
                                              m = NULL,
                                              m_select = c("best_test", "best_train"),
                                              n_draws = 200,
                                              ci_level = NULL,
                                              times = NULL,
                                              curve_type = c("marginal", "mean_lp"),
                                              seed = 1,
                                              show_km_ci = FALSE) {
  if (!requireNamespace("survival", quietly = TRUE)) stop("Package 'survival' is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  
  m_select <- match.arg(m_select)
  curve_type <- match.arg(curve_type)
  
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
  
  summ <- incremental_res$summary
  if (is.null(m)) {
    if (m_select == "best_test") {
      m <- which.max(summ$metric_mean)
    } else {
      m <- which.max(summ$metric_mean_train)
    }
  }
  m <- as.integer(m)
  if (m < 1L || m > nrow(summ)) stop("m is out of range for incremental_res$summary.")
  
  # Reconstruct incremental predictor order from summary labels
  added_labels <- as.character(summ$variable)
  added_vars <- gsub("^\\+", "", added_labels)
  predictor_set <- added_vars[seq_len(m)]
  predictor_set <- predictor_set[predictor_set %in% colnames(xdf)]
  if (length(predictor_set) == 0L) stop("No predictors matched xdata columns for chosen m.")
  
  det <- incremental_res$details[[m]]
  if (is.null(det$beta_coefficients) || !is.matrix(det$beta_coefficients)) {
    stop("incremental_res$details[[m]]$beta_coefficients must be a matrix.")
  }
  
  beta_mat <- det$beta_coefficients
  beta_cols <- colnames(beta_mat)
  if (is.null(beta_cols)) stop("beta_coefficients matrix must have column names.")
  
  # Drop intercept column if present
  beta_cols_use <- setdiff(beta_cols, "(Intercept)")
  beta_cols_use <- beta_cols_use[beta_cols_use %in% predictor_set]
  if (length(beta_cols_use) == 0L) stop("No overlapping coefficient columns found for predictor set.")
  
  beta_sub <- beta_mat[, beta_cols_use, drop = FALSE]
  
  # Keep draws with complete coefficients
  ok_draw <- rowSums(is.na(beta_sub)) == 0L
  beta_sub <- beta_sub[ok_draw, , drop = FALSE]
  if (nrow(beta_sub) == 0L) stop("No valid coefficient draws available after removing NA draws.")
  
  set.seed(seed)
  if (!is.null(n_draws) && is.finite(n_draws)) {
    n_draws <- min(as.integer(n_draws), nrow(beta_sub))
    draw_idx <- sample.int(nrow(beta_sub), size = n_draws, replace = FALSE)
    beta_sub <- beta_sub[draw_idx, , drop = FALSE]
  }
  
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
  actual_lo <- actual_hi <- NULL
  if (show_km_ci) {
    actual_lo <- as.numeric(km_sum$lower)
    actual_hi <- as.numeric(km_sum$upper)
  }
  
  # Predicted curves from refit draws
  X <- as.matrix(xdf[, beta_cols_use, drop = FALSE])
  mode(X) <- "numeric"
  
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
  
  # Compute predicted marginal survival for a single draw
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
    
    # Marginal survival: average over subjects
    out <- vapply(H0, function(h) mean(exp(-h * exp_lp), na.rm = TRUE), numeric(1))
    as.numeric(out)
  }
  
  pred_mat <- matrix(NA_real_, nrow = nrow(beta_sub), ncol = length(times_grid))
  for (i in seq_len(nrow(beta_sub))) {
    pred_mat[i, ] <- predict_curve_one(as.numeric(beta_sub[i, ]))
  }
  
  pred_mat <- pred_mat[rowSums(is.na(pred_mat)) < ncol(pred_mat), , drop = FALSE]
  if (nrow(pred_mat) == 0L) stop("All predicted curves failed; check coefficient draws and data integrity.")
  
  pred_med <- apply(pred_mat, 2, stats::median, na.rm = TRUE)
  pred_lo  <- apply(pred_mat, 2, stats::quantile, probs = q_lo, na.rm = TRUE, names = FALSE)
  pred_hi  <- apply(pred_mat, 2, stats::quantile, probs = q_hi, na.rm = TRUE, names = FALSE)
  
  plot_df <- data.frame(
    time = times_grid,
    actual = actual_surv,
    pred = as.numeric(pred_med),
    pred_lo = as.numeric(pred_lo),
    pred_hi = as.numeric(pred_hi)
  )
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = time)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = pred_lo, ymax = pred_hi), alpha = 0.2) +
    ggplot2::geom_line(ggplot2::aes(y = pred), linewidth = 1) +
    ggplot2::geom_step(ggplot2::aes(y = actual), linewidth = 1, linetype = 2) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(
      x = "Time",
      y = "Survival probability",
      title = "Predicted vs actual survival",
      subtitle = paste0(
        "Model size m = ", m,
        " | Predicted CI from refit draws (n = ", nrow(pred_mat), ")",
        " | Curve type: ", curve_type
      )
    ) +
    ggplot2::theme_minimal()
  
  if (show_km_ci) {
    plot_df$actual_lo <- actual_lo
    plot_df$actual_hi <- actual_hi
    p <- p + ggplot2::geom_ribbon(
      data = plot_df,
      ggplot2::aes(ymin = actual_lo, ymax = actual_hi),
      alpha = 0.1
    )
  }
  
  list(
    plot = p,
    data = plot_df,
    m = m,
    predictors = predictor_set,
    n_draws_used = nrow(pred_mat),
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
    dplyr::mutate(grp_label = factor(lab, levels = stats$lab)) |> 
    mutate(tpr_mean=ifelse(fpr==0.00,0,tpr_mean))
  
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


