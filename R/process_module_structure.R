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
  each_pathway_infor = unique(pathway_infor[pathway_infor$Module_Entry %in% module,
                                            c("Orthology_Entry","Module_Entry","Definition"),
                                            drop = FALSE])

  # Filter sample KO data for orthology entries in the module
  idx = match(unique(each_pathway_infor$Orthology_Entry), rownames(Sample_KO), nomatch = 0L)
  idx = idx[idx > 0L]

  # Check if any matching KOs were found
  if (length(idx) > 0){
    # Same result as base::merge(by = 'Orthology_Entry', all.x = TRUE): rows sorted by
    # Orthology_Entry, Orthology_Entry first column, pathway info attached by match
    sub_Sample_KO = Sample_KO[idx, , drop = F]
    sub_Sample_KO = sub_Sample_KO[order(sub_Sample_KO$Orthology_Entry), , drop = F]
    info_idx = match(sub_Sample_KO$Orthology_Entry, each_pathway_infor$Orthology_Entry)
    sub_Sample_KO_pathway = data.frame(
      Orthology_Entry = sub_Sample_KO$Orthology_Entry,
      sub_Sample_KO[, setdiff(colnames(sub_Sample_KO), "Orthology_Entry"), drop = F],
      each_pathway_infor[info_idx, c("Module_Entry", "Definition"), drop = F],
      stringsAsFactors = FALSE, check.names = FALSE
    )
    rownames(sub_Sample_KO_pathway) = sub_Sample_KO_pathway$Orthology_Entry
  } else {
    sub_Sample_KO_pathway = data.frame()
  }
  return(sub_Sample_KO_pathway)
}
