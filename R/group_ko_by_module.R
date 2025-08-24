#' @title Group KO Abundance Data by Module
#' @description Processes KO abundance data to group by metabolic modules, converting
#'              presence/absence data into module-level KO lists. Handles missing KOs
#'              and maintains sample-specific KO profiles.
#' @param pathway_infor Data frame containing pathway information, see examples.
#' @param Sample_KO_abundance Data frame of KO abundances with:
#'                            - Rows as KO identifiers
#'                            - Columns as samples
#'                            - Orthology_Entry column
#' @return A data frame where:
#'         - Rows are module names
#'         - Columns are samples
#'         - Cell values are space-separated lists of present KOs
#'         - Empty strings for modules with no detected KOs
#' @export
group_ko_by_module <- function(pathway_infor, Sample_KO_abundance) {
  replace_values <- function(df) {
    df[] <- lapply(df, function(x) {
      ifelse(x == 1, df$Orthology_Entry, ifelse(x == 0, "", x))
    })
    return(df)
  }
  Sample_ko_group_by_module = Sample_KO_abundance %>%
    {convert_abundance_to_presence(.)} %>%
    {dplyr::select(., -c(Orthology_Entry))} %>%
    {add_rows_if_not_exists(., add_rows = unique(pathway_infor$Orthology_Entry))} %>%
    {dplyr::mutate(., Orthology_Entry = rownames(.))} %>%
    {base::merge(., unique(pathway_infor[, c("Orthology_Entry","Module_Name")]), all.x = T, by = 'Orthology_Entry')} %>%
    {replace_values(.)} %>%
    {dplyr::select(., -Orthology_Entry)}  %>%
    dplyr::group_by(Module_Name) %>%
    dplyr::summarize(dplyr::across(dplyr::everything(), ~paste(.[. != ""], collapse = " "))) %>%
    dplyr::ungroup(.) %>%
    tibble::column_to_rownames(., 'Module_Name')
  return(Sample_ko_group_by_module)
}
