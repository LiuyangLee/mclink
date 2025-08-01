#' @title Process All Pathways Analysis
#'
#' @description Processes module sample data across all pathways with specified scaling methods.
#'              Handles different comparison methods and outputs results by pathway.
#' @param pathway_infor Data frame containing pathway information, see examples.
#' @param Module_Sample Data frame of module sample data to process
#' @param out_DIR_Module_Sample_by_pathway Output directory for pathway-specific results
#' @param compare_method Comparison method to use: "log" (log10 transform),
#'                       "avg" (average calculation), or "round" (simple rounding)
#' @param plus_scale_method Scaling method for plus-separated KOs ("mean", "min", or "max")
#' @param comma_scale_method Scaling method for comma-separated KOs ("sum" or "max")
#'
#' @return Main outputs are written to:
#'         - Combined module file (All_modules.*.tsv)
#'         - Pathway-specific files (via create_sub_module_sample)
#' @export
#'
#' @examples
#' \dontrun{
#' process_all_pathways(
#'     pathway_infor, Module_Sample, out_DIR_Module_Sample_by_pathway,
#'     compare_method = c("log", "avg", "round"),
#'     plus_scale_method,comma_scale_method)
#' }
process_all_pathways <- function(pathway_infor, Module_Sample, out_DIR_Module_Sample_by_pathway, compare_method = c("log", "avg", "round"),
                                 plus_scale_method,comma_scale_method) {
  if (compare_method == "log") {
    Module_Sample_scale <- Module_Sample %>%
      {round(log10((.) + 1), 4)}
  } else if (compare_method == "avg") {
    Module_Sample_scale <- Module_Sample %>%
      {round(ata_cal(.), 4)}
  } else if (compare_method == "round") {
    Module_Sample_scale <- Module_Sample %>%
      {round(., 4)}
  } else {
    stop("Invalid compare_method. Choose anyone of the 'log', 'avg' or 'round'.")
  }

  out_file <- base::file.path(out_DIR_Module_Sample_by_pathway, paste0('All_modules',".tsv"))


  Module_Sample_out = Module_Sample_scale %>%
    {tibble::rownames_to_column(., var = "Module_Name")}

  data.table::fwrite(Module_Sample_out, file = out_file, sep = "\t", quote = FALSE)

  create_sub_module_sample(pathway_infor, Module_Sample_scale, out_DIR_Module_Sample_by_pathway, plus_scale_method,comma_scale_method)
}
