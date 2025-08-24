#' @title Process Space-Separated KOs with Mean Calculation
#'
#' @description Handles space-separated KOs by calculating the mean abundance across all specified KOs.
#'              Processes multiple KOs separated by spaces and aggregates them into a single row.
#'              Note: For mean calculation, uses the sum of KO abundances divided by total number of KOs,
#'              including those with zero abundance in all samples.
#'
#' @param module_abundance Data frame containing KO abundance data with required columns:
#'                         Orthology_Entry, Module_Entry, Definition
#' @param KOs Character string of space-separated KO IDs (default: "K14126 K14127 K14128")
#' @param aggregrate_rowname Base name for row aggregation (default: 'step_1')
#' @param step_count Processing step counter (default: 1)
#'
#' @return List containing:
#'         - abundance_table: Processed data with mean values
#'         - step_count: Updated step counter
#'         - abundance_log: log
#' @export
process_step_space <- function(module_abundance, KOs = c("K14126 K14127 K14128"),
                               aggregrate_rowname, step_count = 1) {
  # Process space-separated KOs (does not enter loop, no new_step_name added)
  # Uses mean calculation similar to plus-separated KOs
  KOs_scale <- base::strsplit(KOs, " ")[[1]]
  #cat(paste0('\t\tRunning KOs space: ', aggregrate_rowname, " = ", KOs_scale, '\n'))
  log_messages <- list(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']','    ','Running KOs space: ', aggregrate_rowname, " = ", KOs_scale))
  # Prepare and process abundance data
  abundance_table = module_abundance %>%
    {
      rownames(.) = (.$Orthology_Entry)
      (.)
    } %>%
    {dplyr::select(., -c(Orthology_Entry, Module_Entry, Definition))} %>%
    {add_rows_if_not_exists(., add_rows = KOs_scale)} %>%
    {.[rownames(.) %in% KOs_scale, , drop = F]} %>%
    {t(colSums(.)/length(KOs_scale))} %>%  # Calculate mean across all KOs
    {as.data.frame((.), row.names = aggregrate_rowname)} %>%
    {dplyr::mutate(.,
                   Orthology_Entry = rownames(.),
                   Module_Entry = unique(module_abundance$Module_Entry),
                   Definition = unique(module_abundance$Definition))}

  step_count = step_count + 1
  return(list(abundance_table = abundance_table, step_count = step_count, abundance_log = log_messages))
}
