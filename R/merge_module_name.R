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
merge_module_name <- function(pathway_infor, module_table) {
  module_infor <- unique(pathway_infor[, c("Module_Entry", "Module_Name")])

  # Same result as base::merge(module_table, module_infor, by = 'Module_Entry',
  # all.x = TRUE): rows sorted by Module_Entry, Module_Name attached by match
  mt <- module_table[order(module_table$Module_Entry), , drop = FALSE]
  sample_cols <- setdiff(colnames(mt), c("Module_Entry", "Orthology_Entry", "Definition"))
  Module_Sample <- mt[, sample_cols, drop = FALSE]
  rownames(Module_Sample) <- module_infor$Module_Name[match(mt$Module_Entry, module_infor$Module_Entry)]

  Module_Sample <- add_rows_if_not_exists(Module_Sample, add_rows = unique(module_infor$Module_Name))
  return(Module_Sample)
}
