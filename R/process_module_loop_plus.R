#' @title Process Module Components with Plus Sign Handling
#'
#' @description Processes a vector of KOs, applying different handling methods depending on whether
#'              they contain plus signs or not. Handles pathway definitions with required components
#'              (plus-separated KOs representing complex subunits).
#'
#' @param KO_vector Character vector of KO identifiers to process
#' @param module_abundance Data frame containing KO abundance data
#' @param process_step_plus Function to handle plus-separated KOs (required components)
#' @param process_step_direct Function to handle individual KOs
#' @param aggregrate_rowname Base name for row aggregation (default: 'step_1')
#' @param step_count Counter for processing steps (default: 1)
#' @param plus_scale_method Scaling method for plus-separated KOs ("mean", "min", or "max")
#'
#' @return List containing:
#'         - abundance_table: Processed abundance values
#'         - step_count: Updated step counter
#'         - abundance_log: log
#' @export
process_module_loop_plus <- function(KO_vector,
                                     module_abundance,
                                     process_step_plus,
                                     process_step_direct,
                                     aggregrate_rowname = 'step_1',
                                     step_count = 1,
                                     plus_scale_method) {
  # Initialize empty data frames for results
  abundance_table = data.frame()
  abundance_table.tmp = data.frame()
  log_messages = list()

  # Process each KO in the input vector
  for (KOs in KO_vector) {
    # Check if KO string contains plus signs (required components)
    has_plus <- base::grepl("\\+", KOs)

    if (has_plus) {
      # Process plus-separated KOs using specified method
      result <- process_step_plus(
        module_abundance,
        KOs,
        aggregrate_rowname = paste0(aggregrate_rowname, '_', step_count),
        step_count,
        plus_scale_method
      )
      abundance_table.tmp = result[['abundance_table']]
      step_count = result[['step_count']]
      log_messages.tmp = result[['abundance_log']]
    }
    else {
      # Process single KO directly (no plus signs)
      abundance_list = process_step_direct(module_abundance, KOs)
      abundance_table.tmp = abundance_list[["abundance_table"]]
      log_messages.tmp = abundance_list[["abundance_log"]]

    }
    # Combine results while removing duplicates
    abundance_table <- rbind(abundance_table, abundance_table.tmp) %>% unique(.)
    log_messages = c(log_messages, log_messages.tmp)
  }
  return(list(abundance_table = abundance_table, step_count = step_count, abundance_log = log_messages))
}
