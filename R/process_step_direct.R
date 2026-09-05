#' @title Direct KO Processing Without Special Handling
#'
#' @description Processes KO abundances directly without any special scaling or aggregation.
#'              Simply extracts the specified KOs from the abundance table while maintaining
#'              the original module metadata.
#'
#' @param module_abundance Data frame containing KO abundance data with required columns:
#'                         Orthology_Entry, Module_Entry, Definition
#' @param KOs Character vector of KO IDs to extract (default: c("K14126","K14128","K14127"))
#' @return List containing:
#'         - abundance_table: A subset of the input data frame containing only the specified KOs,
#'         with original module metadata preserved
#'         - abundance_log: log
#' @export
process_step_direct <- function(module_abundance, KOs = c("K14126","K14128","K14127")) {
  # Process KOs directly without any special handling
  #cat(paste0('\t\tRunning direct KOs: ', KOs, '\n'))
  log_messages <- list(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']','    ','Running KOs direct: ', KOs))
  # Prepare abundance table with selected KOs (present KOs keep values, missing KOs become zero rows)
  sample_cols = setdiff(colnames(module_abundance), c("Orthology_Entry", "Module_Entry", "Definition"))
  hit = module_abundance$Orthology_Entry %in% KOs
  abundance_table = module_abundance[hit, sample_cols, drop = F]
  rownames(abundance_table) = module_abundance$Orthology_Entry[hit]
  abundance_table = add_rows_if_not_exists(abundance_table, add_rows = KOs)
  abundance_table$Orthology_Entry = rownames(abundance_table)
  abundance_table$Module_Entry = unique(module_abundance$Module_Entry)
  abundance_table$Definition = unique(module_abundance$Definition)
  return(list(abundance_table = abundance_table, abundance_log = log_messages))
}
