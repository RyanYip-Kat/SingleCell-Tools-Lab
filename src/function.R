library(compassR)
library(stringr)
library(magrittr)
library(argparse)
library(dplyr)

##################
cohens_d <- function(group_A_values, group_B_values) {
    n_A <- length(group_A_values)
    n_B <- length(group_B_values)
    mu_A <- mean(group_A_values)
    mu_B <- mean(group_B_values)
    var_A <- var(group_A_values)
    var_B <- var(group_B_values)
    pooled_sd <- sqrt(((n_A - 1) * var_A + (n_B - 1) * var_B) / (n_A + n_B - 2))
    cohen_d <- (mu_A - mu_B) / pooled_sd
    cohen_d
}

conduct_wilcoxon_test_1 = function(compass_data,group.by="status",for_metareactions = TRUE) {

  consistencies_matrix=as.data.frame(t(compass_data$reaction_consistencies))
  metareaction_ids <- rownames(compass_data$reaction_consistencies)
  cells=compass_data$cell_metadata$cells
  consistencies_matrix<-consistencies_matrix[cells,]
  groups=unique(compass_data$cell_metadata[[group.by]])
  #stopifnot(length(groups)==2)

  print(groups)
  cell_A=cells[compass_data$cell_metadata[[group.by]]%in%groups[1]]
  cell_B=cells[compass_data$cell_metadata[[group.by]]%in%groups[2]]

  if (0 < length(intersect(cell_A,cell_B))) {
    message("Groups A and B are not mutually exclusive. Continuing anyways ...")
  }

  group_A_values_per_metareaction=na.omit(consistencies_matrix[cell_A,])
  group_B_values_per_metareaction=na.omit(consistencies_matrix[cell_B,])
  wilcoxon_results <-
    purrr::map_dfr(
      metareaction_ids,
      function(metareaction_id) {
        group_A_values <- group_A_values_per_metareaction[,metareaction_id]
        group_B_values <- group_B_values_per_metareaction[,metareaction_id]
        wilcoxon_result_obj <- wilcox.test(group_A_values, group_B_values)
        wilcoxon_result_tbl <- data.frame(
          metareaction_id = metareaction_id,
          wilcoxon_statistic = wilcoxon_result_obj$statistic,
          cohens_d = cohens_d(group_A_values, group_B_values),
          p_value = wilcoxon_result_obj$p.value,
          stringsAsFactors = FALSE
        )
        wilcoxon_result_tbl
      }
    ) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(adjusted_p_value = p.adjust(dplyr::pull(., p_value), method = "BH"))
  if (!for_metareactions) {
    wilcoxon_results %<>% dplyr::rename(reaction_id = metareaction_id)
  }
  wilcoxon_results
}


conduct_wilcoxon_test_2 = function(compass_data, group_A_cell_ids, group_B_cell_ids, ..., for_metareactions = TRUE) {
  if (0 < length(intersect(group_A_cell_ids, group_B_cell_ids))) {
    message("Groups A and B are not mutually exclusive. Continuing anyways ...")
  }
  consistencies_matrix=compass_data$reaction_consistencies
  group_A_values_per_metareaction <-
    consistencies_matrix %>%
    t() %>%
    tibble::as_tibble(rownames = compass_data$settings$cell_id_col_name) %>%
    dplyr::right_join(
      tibble::tibble(!!compass_data$settings$cell_id_col_name := group_A_cell_ids),
      by = compass_data$settings$cell_id_col_name
    ) %>%
    as.data.frame()%>%na.omit()
  group_B_values_per_metareaction <-
    consistencies_matrix %>%
    t() %>%
    tibble::as_tibble(rownames = compass_data$settings$cell_id_col_name) %>%
    dplyr::right_join(
      tibble::tibble(!!compass_data$settings$cell_id_col_name := group_B_cell_ids),
      by = compass_data$settings$cell_id_col_name
    ) %>%
    as.data.frame()%>%na.omit()
  metareaction_ids <- rownames(consistencies_matrix)
  wilcoxon_results <-
    purrr::map_dfr(
      metareaction_ids,
      function(metareaction_id) {
        group_A_values <- group_A_values_per_metareaction[,metareaction_id]
        group_B_values <- group_B_values_per_metareaction[,metareaction_id]
        wilcoxon_result_obj <- wilcox.test(group_A_values, group_B_values)
        wilcoxon_result_tbl <- data.frame(
          metareaction_id = metareaction_id,
          wilcoxon_statistic = wilcoxon_result_obj$statistic,
          cohens_d = cohens_d(group_A_values, group_B_values),
          p_value = wilcoxon_result_obj$p.value,
          stringsAsFactors = FALSE
        )
        wilcoxon_result_tbl
      }
    ) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(adjusted_p_value = p.adjust(dplyr::pull(., p_value), method = "BH"))
  if (!for_metareactions) {
    wilcoxon_results %<>% dplyr::rename(reaction_id = metareaction_id)
  }
  wilcoxon_results
  }

