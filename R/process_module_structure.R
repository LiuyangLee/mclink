#' @title Process Module Structure Data
#'
#' @description Filters and merges pathway information with sample KO data for a specific module.
#'              Returns a merged dataframe containing KO abundance data and pathway definitions.
#' @param pathway_infor Data frame containing pathway information, see examples.
#' @param Sample_KO Dataframe containing KO abundance data with KO IDs as row names
#' @param module Character string of the module ID to process (e.g. "M00563")
#'
#' @return A merged dataframe containing:
#'         - KO abundance data for the specified module
#'         - Corresponding pathway information
#'         - Empty dataframe if no matching KOs found
#' @export
process_module_structure <- function(pathway_infor, Sample_KO, module) {
  # Filter pathway information for the specified module
  each_pathway_infor = pathway_infor %>%
    {.[(.$Module_Entry) %in% module, c("Orthology_Entry","Module_Entry","Definition"),drop = F]} %>%
    {unique(.)}

  # Filter sample KO data for orthology entries in the module
  sub_Sample_KO = Sample_KO %>%
    {.[rownames(.) %in% unique(each_pathway_infor$Orthology_Entry),,drop = F]}

  # Check if any matching KOs were found
  if (length(sub_Sample_KO$Orthology_Entry)>0){
    # Merge KO data with pathway information
    sub_Sample_KO_pathway = sub_Sample_KO %>%
      base::merge(., each_pathway_infor, by = 'Orthology_Entry', all.x = T) %>%
      {
        rownames(.) = .$Orthology_Entry
        (.)
      }
  } else {
    sub_Sample_KO_pathway = data.frame()
  }
  return(sub_Sample_KO_pathway)
}
