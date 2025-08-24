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
  module_abundance %>%
    {
      rownames(.) = (.$Orthology_Entry)
      (.)
    } %>%
    {dplyr::select(., -c(Orthology_Entry))} %>%
    # {apply(., 2, function(x) ifelse(x > 0, 1, 0))} %>%
    {ifelse(as.matrix(.) > 0, 1, 0)} %>%
    {as.data.frame(.)} %>%
    {dplyr::mutate(.,Orthology_Entry = rownames(.)
    )}
}
