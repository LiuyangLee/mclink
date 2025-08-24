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
  abundance_table = module_abundance %>%
    {
      rownames(.) = (.$Orthology_Entry)
      (.)
    } %>%
    {dplyr::select(., -c(Orthology_Entry, Module_Entry, Definition))} %>%
    {add_rows_if_not_exists(., add_rows = KOs)} %>%
    {dplyr::mutate(.,
                   Orthology_Entry = rownames(.),
                   Module_Entry = unique(module_abundance$Module_Entry),
                   Definition = unique(module_abundance$Definition)
    )} %>%
    {.[rownames(.) %in% KOs, , drop = F]}
  return(list(abundance_table = abundance_table, abundance_log = log_messages))
}
