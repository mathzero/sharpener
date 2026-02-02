# Dataset-specific helper plots (kept here for reference)
# NOTE: These functions assume objects like `active_ids` and `df_relapse` exist.

plotBaselineCharacteristics <- function(data, cat = "At what approximate age did atopic dermatitis begin?") {
  library(scales)
  data %>%
    filter(
      subject_identifier_for_the_study %in% active_ids,
      evaluation_test_name == cat) %>%
    mutate(
      label = fct_reorder(character_result_finding_in_std_format,
                          numeric_result_finding_in_standard_units,
                          .fun  = min),
      label = factor(label, ordered = TRUE)
    ) %>%
    ggplot(aes(y = label)) +
    geom_histogram(stat = "count") +
    OverReact::scale_color_imperial() +
    OverReact::theme_react() +
    labs(title = cat, y = "") +
    scale_x_continuous(breaks = breaks_pretty()) +
    theme(axis.text.x = element_text(size = 15),
          axis.text.y = element_text(size = 15),
          axis.title.x = element_text(size = 15))
}

plotScoresOverTime <- function(data, cat = "Quantification of Staphylococcus aureus (log 10)") {
  library(scales)
  plotdat <- data %>%
    filter(subject_identifier_for_the_study %in% active_ids,
           evaluation_test_name == cat) %>%
    group_by(lesion, time, category_for_evaluation) %>%
    summarise(numeric_result_finding_in_standard_units = mean(numeric_result_finding_in_standard_units)) %>%
    left_join(df_relapse)

  plotdat |>
    ggplot() +
    geom_point(aes(x = time, y = numeric_result_finding_in_standard_units, col = category_for_evaluation)) +
    geom_path(aes(x = time, y = numeric_result_finding_in_standard_units, col = category_for_evaluation)) +
    geom_point(data = plotdat |> group_by(lesion) |> arrange(time) |> slice_tail(n = 1) |> ungroup(),
               aes(y = numeric_result_finding_in_standard_units, x = time_of_relapse), shape = 1, size = 3, col = "red") +
    OverReact::scale_color_imperial() +
    scale_x_continuous(breaks = breaks_pretty()) +
    facet_wrap(. ~ lesion) +
    OverReact::theme_react() +
    labs(title = cat, col = "Evaluation \nmethod")
}

plotScoresOverTimeCompare <- function(data,
                                      column_with_cats = "evaluation_test_name",
                                      cat = "Quantification of Staphylococcus aureus (log 10)") {
  library(scales)
  plotdat <- data %>%
    filter(subject_identifier_for_the_study %in% active_ids,
           !!sym(column_with_cats) == cat) %>%
    group_by(lesion, time, category_for_evaluation) %>%
    summarise(numeric_result_finding_in_standard_units = mean(numeric_result_finding_in_standard_units)) %>%
    left_join(df_relapse)

  plotdat |>
    ggplot() +
    geom_point(aes(x = time, y = numeric_result_finding_in_standard_units, col = category_for_evaluation)) +
    geom_path(aes(x = time, y = numeric_result_finding_in_standard_units, col = category_for_evaluation)) +
    geom_point(data = plotdat |> group_by(lesion) |> arrange(time) |> slice_tail(n = 1) |> ungroup(),
               aes(y = numeric_result_finding_in_standard_units, x = time_of_relapse), shape = 1, size = 3, col = "red") +
    OverReact::scale_color_imperial() +
    scale_x_continuous(breaks = breaks_pretty()) +
    facet_wrap(. ~ lesion) +
    OverReact::theme_react() +
    labs(title = cat, col = "Evaluation \nmethod")
}
