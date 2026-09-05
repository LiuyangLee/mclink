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
  pres = convert_abundance_to_presence(Sample_KO_abundance)
  sample_cols = setdiff(colnames(pres), "Orthology_Entry")
  m = as.matrix(pres[, sample_cols, drop = FALSE])
  rownames(m) = pres$Orthology_Entry

  # KO -> module map, sorted by Orthology_Entry to reproduce the row order of the
  # former base::merge() pipeline (KOs are pasted in that order within each module)
  map = unique(pathway_infor[, c("Orthology_Entry", "Module_Name")])
  map = map[order(map$Orthology_Entry), , drop = FALSE]

  row_idx = match(map$Orthology_Entry, rownames(m))
  found = !is.na(row_idx)

  modules = sort(unique(map$Module_Name))
  grp_idx = split(seq_len(nrow(map)), factor(map$Module_Name, levels = modules))

  out = matrix("", nrow = length(modules), ncol = length(sample_cols),
               dimnames = list(modules, sample_cols))
  ko_names = map$Orthology_Entry
  for (j in seq_along(sample_cols)) {
    v = numeric(nrow(map))
    v[found] = m[row_idx[found], j]
    is_one = v == 1
    s = character(nrow(map))
    s[is_one %in% TRUE] = ko_names[is_one %in% TRUE]
    s[is.na(is_one)] = "NA"  # NA abundances become the literal string "NA" in pasted lists
    out[, j] = vapply(grp_idx, function(ii) paste(s[ii][s[ii] != ""], collapse = " "),
                      character(1))
  }
  data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
}
