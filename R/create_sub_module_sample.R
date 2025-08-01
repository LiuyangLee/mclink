#' @title Create and Export Pathway-Specific Module Sample Files
#'
#' @description Processes pathway information and module sample data to create and export
#'              individual pathway-specific files containing scaled module data.
#'
#' @param pathway_infor Data frame containing pathway information, see examples.
#' @param Module_Sample_scale Data frame containing scaled module sample data with
#'                            module names as row names
#' @param out_DIR_Module_Sample_by_pathway Character string specifying output directory
#' @param plus_scale_method Scaling method for plus-separated KOs ("mean", "min", or "max")
#' @param comma_scale_method Scaling method for comma-separated KOs ("sum" or "max")
#'
#' @return None (writes files to disk)
#' @export
#'
#' @examples
#' \dontrun{
#' create_sub_module_sample(
#'     pathway_infor, Module_Sample_scale,
#'     out_DIR_Module_Sample_by_pathway, plus_scale_method, comma_scale_method)
#' }
create_sub_module_sample <- function(pathway_infor, Module_Sample_scale, out_DIR_Module_Sample_by_pathway, plus_scale_method, comma_scale_method) {
  cat(paste0('\n[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),'] Output Pathway: ',unique(pathway_infor$Level_2),'\n\n'))
  for (pathway in unique(pathway_infor$Level_2)) {
    pathway_scale <- gsub(' ', '_', pathway)
    out_file <- base::file.path(out_DIR_Module_Sample_by_pathway, paste0(pathway_scale,".tsv"))

    each_pathway_infor <- pathway_infor %>%
      {.[(.$Level_2) %in% pathway, c('Module_Name', 'Module_Entry'),drop = F]} %>%
      unique()
    sub_Module_Sample_out = Module_Sample_scale %>%
      {.[rownames(.) %in% unique(each_pathway_infor$Module_Name),]} %>%
      {tibble::rownames_to_column(., var = "Module_Name")}
    data.table::fwrite(sub_Module_Sample_out, file = out_file, sep = "\t", quote = FALSE)
  }
}
