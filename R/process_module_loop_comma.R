#' @title Process Module Components with Comma Handling
#'
#' @description Processes a vector of KOs, applying different handling methods depending on whether
#'              they contain commas or not. Useful for processing complex pathway definitions with
#'              alternative KOs.
#'
#' @param KO_vector Character vector of KO identifiers to process
#' @param module_abundance Data frame containing KO abundance data
#' @param process_step_comma Function to handle comma-separated KOs (alternative forms)
#' @param process_step_direct Function to handle individual KOs
#' @param aggregrate_rowname Base name for row aggregation (default: 'step_1')
#' @param step_count Counter for processing steps (default: 1)
#' @param comma_scale_method Scaling method for comma-separated KOs ("sum" or "max")
#'
#' @return List containing:
#'         - abundance_table: Processed abundance Data frame
#'         - step_count: Updated step counter
#'         - abundance_log: log
#' @export
process_module_loop_comma <- function(KO_vector,
                                      module_abundance,
                                      process_step_comma,
                                      process_step_direct,
                                      aggregrate_rowname = 'step_1',
                                      step_count = 1,
                                      comma_scale_method) {
  # Initialize empty data frames for results
  abundance_table = data.frame()
  abundance_table.tmp = data.frame()
  log_messages = list()

  # Process each KO in the input vector
  for (KOs in KO_vector) {
    # Check if KO string contains commas (alternative KOs)
    has_comma <- base::grepl(",", KOs)

    if (has_comma) {
      # Process comma-separated KOs using specified method
      result = process_step_comma(
        module_abundance,
        KOs,
        aggregrate_rowname = paste0(aggregrate_rowname, '_', step_count),
        step_count,
        comma_scale_method
      )
      abundance_table.tmp = result[['abundance_table']]
      step_count = result[['step_count']]
      log_messages.tmp = result[['abundance_log']]
    }
    else {
      # Process single KO directly
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
