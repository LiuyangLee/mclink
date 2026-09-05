#' @title Convert Abundance Values to Presence/Absence Indicators
#'
#' @description Transforms a numeric abundance matrix into a binary presence/absence matrix,
#'              where 1 indicates presence (abundance > 0) and 0 indicates absence.
#'              Preserves row names as Orthology_Entry column in the output.
#'
#' @param module_abundance A data frame containing KO abundance data, must include:
#'                         - Rows named by Orthology_Entry (KO identifiers)
#'                         - Numeric columns representing sample abundances
#'                         - An Orthology_Entry column
#'
#' @return A data frame with:
#'         - Binary values (1 = present, 0 = absent) for each sample
#'         - Original row names preserved in Orthology_Entry column
#'         - Same dimensions as input (excluding the Orthology_Entry column)
#' @export
convert_abundance_to_presence <- function(module_abundance) {
  value_cols = setdiff(colnames(module_abundance), "Orthology_Entry")
  # Short-circuit: values already binary (0/1) need no conversion
  already_binary = all(vapply(module_abundance[, value_cols, drop = FALSE],
                              function(x) all(x %in% c(0, 1)), logical(1)))
  if (already_binary) {
    rownames(module_abundance) = module_abundance$Orthology_Entry
    return(module_abundance)
  }
  ortho = module_abundance$Orthology_Entry
  m = as.matrix(module_abundance[, value_cols, drop = FALSE])
  out = as.data.frame((m > 0) * 1)
  rownames(out) = ortho
  out$Orthology_Entry = ortho
  out
}
