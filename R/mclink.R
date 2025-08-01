#' @importFrom dplyr %>%
NULL
#' @importFrom utils packageVersion data
#' @export KO_pathway_ref
#' @export KO_Sample_wide
NULL  # 仅用于生成NAMESPACE的占位符
#' Metabolic Pathway Coverage Analysis
#'
#' @description
#' Analyzes metabolic pathway completeness/abundance from (meta)genome KO presence/abundance data.
#' Can use either built-in KEGG datasets or user-provided data frames. Output includes pathway coverage
#' metrics and detected KOs in each pathway/module.
#'
#' @param ref Pathway information data frame. When `NULL` (default),
#'        uses the built-in \code{\link{KO_pathway_ref}} dataset. Must contain the
#'        same columns as the built-in dataset if providing custom data.
#' @param data Sample KO abundance data frame. When `NULL` (default),
#'        uses the built-in \code{\link{KO_Sample_wide}} dataset.
#' @param table_feature Analysis type, either:
#'        \itemize{
#'          \item "completeness" (binary presence/absence, default)
#'          \item "abundance" (weighted by KO abundance)
#'        }
#' @param plus_scale_method Scaling method for plus-separated KOs (subunits/complexes):
#'        \itemize{
#'          \item "mean" - Moderate approach (default)
#'          \item "min" - Rigorous/conservative estimate
#'          \item "max" - Liberal estimate
#'        }
#' @param comma_scale_method Scaling method for comma-separated KOs (isoforms/alternatives):
#'        \itemize{
#'          \item "max" - For completeness analysis (default)
#'          \item "sum" - For abundance analysis
#'        }
#' @param out_dir Output directory path. If `NULL` (default), results are only
#'        returned as R objects without writing files.
#' @param split_by_pathway Logical. If `TRUE`, splits results by pathway/module.
#'        Requires non-NULL `out_dir`. Default: `FALSE`.
#'
#' @return A list containing:
#' \itemize{
#'   \item coverage - Data frame with pathway coverage metrics
#'   \item detected_KOs - List of detected KOs per pathway/module
#' }
#'
#' If `out_dir` is specified, results are also written as TSV files.
#'
#' @export
#'
#' @examples
#' data(KO_pathway_ref)
#' data(KO_Sample_wide)
#' selected_modules <- c("M00176","M00165","M00173","M00374","M00375","M00376","M00377")
#' KO_pathway_ref_selected <- KO_pathway_ref[KO_pathway_ref$Module_Entry %in% selected_modules, ]
#' mc_list =
#'   mclink(ref = KO_pathway_ref_selected,
#'          data = KO_Sample_wide,
#'          table_feature = "completeness",
#'          plus_scale_method = "min",
#'          comma_scale_method = "max")
#' mc_coverage = mc_list[["coverage"]]
#' mc_detected_KOs = mc_list[["detected_KOs"]]
mclink <- function(ref = NULL,
                   data = NULL,
                   table_feature = "completeness",
                   plus_scale_method = "mean",
                   comma_scale_method = "max",
                   out_dir = NULL,
                   split_by_pathway = FALSE) {

  # Load default datasets if not provided
  if (is.null(ref)) {
    utils::data("KO_pathway_ref", envir = environment())
    ref <- KO_pathway_ref
  }

  if (is.null(data)) {
    utils::data("KO_Sample_wide", envir = environment())
    data <- KO_Sample_wide
  }

  # Input validation
  stopifnot(
    is.data.frame(ref),
    is.data.frame(data),
    table_feature %in% c("completeness", "abundance"),
    plus_scale_method %in% c("mean", "min", "max"),
    comma_scale_method %in% c("max", "sum")
  )

  message(paste0('\n[', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ']', ' mclink started!\n\n'))
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Input Sample-KO table type:',table_feature,'\n\n')
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Scale method for plus:',plus_scale_method,'\n\n')
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Scale method for comma:',comma_scale_method,'\n\n')

  ##########################      pathway infor input      ##########################
  # 'MASH_KEGG_Energy_metabolism.20240315.tsv'
  pathway_infor = read_and_process_pathway_infor(ref)
  module_level <- pathway_infor %>%
    dplyr::select(Module_Entry, Level_2, Level_3, Module_Name, Definition) %>%
    {unique(.)}
  #cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'dim(module_level)...,',dim(module_level),'\n\n')
  #print(head(module_level))
  ##########################      Sample-KO depth input      ##########################
  Sample_KO_abundance = read_and_process_KO_table(data, pathway_infor)
  ##############    convert presence     ##############
  if (table_feature == "completeness") {
    cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Converting abundance table to completeness table...\n\n')
    Sample_KO_abundance = convert_abundance_to_presence(Sample_KO_abundance)
  }
  ##############    Calculating coverage of respective modules     ##############
  module_table_coverage = process_all_modules(pathway_infor, Sample_KO = Sample_KO_abundance,
                                              plus_scale_method,
                                              comma_scale_method)
  Module_Sample_coverage = merge_module_name(pathway_infor, module_table = module_table_coverage)
  ##############    Output present KOs of respective modules    ##############
  cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),']'),'Summarizing present KOs of respective modules...\n\n')
  Module_Sample_KO_list = group_ko_by_module(pathway_infor, Sample_KO_abundance) %>%
    {.[match(rownames(Module_Sample_coverage), rownames(.)), match(colnames(Module_Sample_coverage), colnames(.))]}
  mc_detected_KOs = Module_Sample_KO_list %>%
    {tibble::rownames_to_column(., var = "Module_Name")} %>%
    {base::merge(module_level, ., by = "Module_Name", all.x = T)} %>%
    dplyr::arrange(match(Module_Name, module_level$Module_Name)) %>%
    dplyr::select(c(Module_Entry,Level_2,Level_3,Module_Name, Definition), dplyr::everything())
  ##############    Output coverage of respective modules for all pathways    ##############
  mc_coverage = Module_Sample_coverage %>%
    {tibble::rownames_to_column(., var = "Module_Name")} %>%
    {base::merge(module_level, ., by = "Module_Name", all.x = T)} %>%
    dplyr::arrange(match(Module_Name, module_level$Module_Name)) %>%
    dplyr::select(c(Module_Entry,Level_2,Level_3,Module_Name, Definition), dplyr::everything())

  mc_list <- list(coverage = mc_coverage, detected_KOs = mc_detected_KOs)

  if (is.null(out_dir)) {
    return(mc_list)
  } else {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }

    if (table_feature == "completeness") {
      # output detected KOs and coverage
      out_file <- base::file.path(out_dir, paste0('mc_detected_KOs', '.tsv'))
      data.table::fwrite(mc_detected_KOs, file = out_file, sep = "\t", quote = FALSE)

      out_file <- base::file.path(out_dir, paste0('mc_completeness', '.tsv'))
      data.table::fwrite(mc_coverage, file = out_file, sep = "\t", quote = FALSE)

      if (split_by_pathway) {
        out_dir_ko_sample_by_pathway_presence <- base::file.path(out_dir, 'mc_by_pathway')
        if (!dir.exists(out_dir_ko_sample_by_pathway_presence)) {
          dir.create(out_dir_ko_sample_by_pathway_presence, recursive = TRUE)
        }

        mc_coverage_arrange <- mc_coverage %>%
          dplyr::select(-c(Module_Entry, Level_2, Level_3, Definition)) %>%
          tibble::column_to_rownames(var = "Module_Name")

        Module_Sample_pathway_presence <- process_all_pathways(pathway_infor,
                                                               mc_coverage_arrange,
                                                               out_dir_ko_sample_by_pathway_presence,
                                                               compare_method = "round",
                                                               plus_scale_method,
                                                               comma_scale_method)
      }
    } else if (table_feature == "abundance") {
      # output detected KOs and coverage
      out_file <- base::file.path(out_dir, paste0('mc_detected_KOs','.tsv'))
      data.table::fwrite(mc_detected_KOs, file = out_file, sep = "\t", quote = FALSE)

      out_file <- base::file.path(out_dir, paste0('mc_abundance','.tsv'))
      data.table::fwrite(mc_coverage, file = out_file, sep = "\t", quote = FALSE)

      if (split_by_pathway) {
        out_dir_ko_sample_by_pathway_log <- base::file.path(out_dir, 'log_scale_by_pathway')
        if (!dir.exists(out_dir_ko_sample_by_pathway_log)) {
          dir.create(out_dir_ko_sample_by_pathway_log, recursive = TRUE)
        }

        out_dir_ko_sample_by_pathway_avg <- base::file.path(out_dir, 'avg_scale_by_pathway')
        if (!dir.exists(out_dir_ko_sample_by_pathway_avg)) {
          dir.create(out_dir_ko_sample_by_pathway_avg, recursive = TRUE)
        }

        mc_coverage_arrange <- mc_coverage %>%
          dplyr::select(-c(Module_Entry, Level_2, Level_3, Definition)) %>%
          tibble::column_to_rownames(var = "Module_Name")

        Module_Sample_pathway_log <- process_all_pathways(pathway_infor,
                                                          mc_coverage_arrange,
                                                          out_dir_ko_sample_by_pathway_log,
                                                          compare_method = "log",
                                                          plus_scale_method,
                                                          comma_scale_method)

        Module_Sample_pathway_avg <- process_all_pathways(pathway_infor,
                                                          mc_coverage_arrange,
                                                          out_dir_ko_sample_by_pathway_avg,
                                                          compare_method = "avg",
                                                          plus_scale_method,
                                                          comma_scale_method)
      }
    }
  }
  return(message(paste0('\n[', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ']', ' mclink finished!\n\n')))
}
