#' @title Merge Module Information with Module Table
#'
#' @description Merges pathway information with a module table to create a sample-by-module matrix
#'              with proper module names. Ensures all modules are represented in the output.
#' @param pathway_infor Data frame containing pathway information, see examples.
#' @param module_table Data frame containing module data with:
#'                     - Module_Entry: Matching module identifiers
#'                     - Orthology_Entry: KO identifiers
#'                     - Definition: Module definitions
#'                     - Sample columns with abundance values
#'
#' @return A data frame where:
#'         - Rows are module names (from Module_Name)
#'         - Columns are samples
#'         - All modules from pathway_infor are represented
#'         - Original row names are replaced with descriptive module names
#' @export
#'
#' @examples
#' \dontrun{
#' merge_module_name(pathway_infor, module_table)
#' }
merge_module_name <- function(pathway_infor, module_table) {
  module_infor <- pathway_infor %>%
    dplyr::select(Module_Entry, Module_Name) %>%
    {unique(.)}
  Module_Sample = module_table %>%
    base::merge(., module_infor, all.x = T, by = 'Module_Entry') %>%
    {dplyr::select(., -c(Module_Entry,Orthology_Entry,Definition))} %>%
    {
      rownames(.) = (.$Module_Name)
      (.)
    } %>%
    {dplyr::select(., -Module_Name)} %>%
    {add_rows_if_not_exists(., add_rows = unique(module_infor$Module_Name))}
  return(Module_Sample)
}
