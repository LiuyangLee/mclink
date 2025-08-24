#' @title Process Plus-Separated KOs with Specified Scaling Method
#'
#' @description Handles plus-separated KOs by applying the specified scaling method (mean, min, or max).
#'              Processes multiple KOs separated by plus signs and aggregates them into a single row.
#'              Note: For mean calculation, uses the sum of existing KO abundances divided by total
#'              number of KOs (including those with zero abundance in all samples).
#'
#' @param module_abundance Data frame containing KO abundance data with required columns:
#'                         Orthology_Entry, Module_Entry, Definition
#' @param KOs Character string of plus-separated KO IDs (default: "K14126+K14127+K14128")
#' @param aggregrate_rowname Base name for row aggregation (default: 'step_1')
#' @param step_count Processing step counter (default: 1)
#' @param plus_scale_method Scaling method for plus-separated KOs ("mean", "min", or "max")
#'
#' @return List containing:
#'         - abundance_table: Processed data with aggregated values
#'         - step_count: Updated step counter
#'         - abundance_log: log
#' @export
process_step_plus <- function(module_abundance, KOs = c("K14126+K14127+K14128"), aggregrate_rowname,
                              step_count = 1, plus_scale_method) {
  # Process plus-separated KOs
  # For plus-separated entries, add new_step_name = paste0(aggregrate_rowname, '_', step_count)
  KOs_scale <- base::strsplit(KOs, "\\+")[[1]]
  #cat(paste0('\t\tRunning KOs plus: ', aggregrate_rowname, " = ", KOs_scale, '\n'))
  log_messages <- list(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']','    ','Running KOs plus: ', aggregrate_rowname, " = ", KOs_scale))
  # Prepare abundance table with selected KOs
  abundance_table = module_abundance %>%
    {
      rownames(.) = (.$Orthology_Entry)
      (.)
    } %>%
    {dplyr::select(., -c(Orthology_Entry, Module_Entry, Definition))} %>%
    {add_rows_if_not_exists(., add_rows = KOs_scale)} %>%
    {.[rownames(.) %in% KOs_scale, , drop = F]}

  # Apply specified scaling method
  if (plus_scale_method == "mean") {
    # Calculate mean: sum of values divided by total number of KOs (including zeros)
    abundance_table_scale = abundance_table %>% {t(colSums(.)/length(KOs_scale))}
  } else if (plus_scale_method == "min") {
    abundance_table_scale = abundance_table %>% {t(apply(., 2, min))}
  } else if (plus_scale_method == "max") {
    abundance_table_scale = abundance_table %>% {t(apply(., 2, max))}
  } else {
    stop(paste("Unknown plus scale method:", plus_scale_method))
  }

  # Format final result with metadata
  abundance_table = abundance_table_scale %>%
    {as.data.frame((.), row.names = aggregrate_rowname)} %>%
    {dplyr::mutate(.,
                   Orthology_Entry = rownames(.),
                   Module_Entry = unique(module_abundance$Module_Entry),
                   Definition = unique(module_abundance$Definition))}

  step_count = step_count + 1
  return(list(abundance_table = abundance_table, step_count = step_count, abundance_log = log_messages))
}
