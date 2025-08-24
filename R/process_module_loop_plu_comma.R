#' @title Process Module Components with Plus and Comma Handling
#'
#' @description Processes a vector of KOs, applying different handling methods depending on whether
#'              they contain plus signs, commas, or both. Handles complex pathway definitions with
#'              both required components (plus-separated) and alternative forms (comma-separated).
#'
#' @param KO_vector Character vector of KO identifiers to process
#' @param module_abundance Data frame containing KO abundance data
#' @param process_step_plus Function to handle plus-separated KOs (required components)
#' @param process_step_comma Function to handle comma-separated KOs (alternative forms)
#' @param process_step_direct Function to handle individual KOs
#' @param aggregrate_rowname Base name for row aggregation (default: 'step_1')
#' @param step_count Counter for processing steps (default: 1)
#' @param plus_scale_method Scaling method for plus-separated KOs ("mean", "min", or "max")
#' @param comma_scale_method Scaling method for comma-separated KOs ("sum" or "max")
#'
#' @return List containing:
#'         - abundance_table: Processed abundance values
#'         - step_count: Updated step counter
#'         - abundance_log: log
#' @export
process_module_loop_plu_comma <- function(KO_vector,
                                          module_abundance,
                                          process_step_plus,
                                          process_step_comma,
                                          process_step_direct,
                                          aggregrate_rowname = 'step_1',
                                          step_count = 1,
                                          plus_scale_method,
                                          comma_scale_method) {
  # Initialize empty data frames for results
  abundance_table = data.frame()
  abundance_table.tmp = data.frame()
  abundance_log <- list()
  # Process each KO in the input vector
  for (KOs in KO_vector) {
    # Check for presence of operators
    has_comma <- base::grepl(",", KOs)
    has_plus <- base::grepl("\\+", KOs)

    # Case 1: Contains both '+' and ',' (e.g., "K14126+K14128,K22516+K00125,K00126")
    if (has_plus && has_comma) {
      # First process plus-separated components, then comma-separated alternatives
      KO_subvector <- base::strsplit(KOs, ",")[[1]]
      loop_plus_result <- process_module_loop_plus(
        KO_subvector,
        module_abundance,
        process_step_plus,
        process_step_direct,
        aggregrate_rowname,
        step_count,
        plus_scale_method
      )

      KO_scale <- paste(rownames(loop_plus_result[['abundance_table']]), collapse = ',')
      step_count <- loop_plus_result[['step_count']]
      abundance_log.tmp <- loop_plus_result[['abundance_log']]

      comma_result <- process_step_comma(
        loop_plus_result[['abundance_table']],
        KO_scale,
        aggregrate_rowname = paste0(aggregrate_rowname, '_', step_count),
        step_count,
        comma_scale_method
      )

      abundance_table.tmp <- comma_result[['abundance_table']]
      step_count <- comma_result[['step_count']]
      abundance_log.tmp <- c(abundance_log.tmp, comma_result[['abundance_log']])
    }
    # Case 2: Contains only commas (alternative forms, e.g., "K22516,K14126")
    else if (has_comma) {
      comma_result <- process_step_comma(
        module_abundance,
        KOs,
        aggregrate_rowname = paste0(aggregrate_rowname, '_', step_count),
        step_count,
        comma_scale_method
      )
      abundance_table.tmp <- comma_result[['abundance_table']]
      step_count <- comma_result[['step_count']]
      abundance_log.tmp <- comma_result[['abundance_log']]
    }
    # Case 3: Contains only plus signs (required components, e.g., "K22516+K00125")
    else if (has_plus) {
      plus_result <- process_step_plus(
        module_abundance,
        KOs,
        aggregrate_rowname = paste0(aggregrate_rowname, '_', step_count),
        step_count,
        plus_scale_method
      )
      abundance_table.tmp <- plus_result[['abundance_table']]
      step_count <- plus_result[['step_count']]
      abundance_log.tmp <- plus_result[['abundance_log']]
    }
    # Case 4: Single KO with no operators
    else {
      direct_result <- process_step_direct(module_abundance, KOs)
      abundance_table.tmp <- direct_result[['abundance_table']]
      abundance_log.tmp <- direct_result[['abundance_log']]
    }

    # Combine results while removing duplicates
    abundance_table <- rbind(abundance_table, abundance_table.tmp) %>% unique(.)
    abundance_log = c(abundance_log, abundance_log.tmp)
  }

  return(list(abundance_table = abundance_table, step_count = step_count, abundance_log = abundance_log))
}
