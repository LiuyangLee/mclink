#' KEGG Orthology (KO) Pathway Information Dataset
#'
#' A comprehensive dataset mapping KEGG Orthology (KO) entries to metabolic pathways,
#' including module hierarchy, definitions, and enzyme annotations.
#'
#' @format A data frame with 3846 rows (KO entries) and 10 variables:
#' \describe{
#'   \item{Orthology_Entry}{Character. KEGG Orthology ID (e.g., "K00844").}
#'   \item{Module_Type}{Character. Type of metabolic module (e.g., "Pathway modules").}
#'   \item{Level_2}{Character. Broad metabolic category (e.g., "Carbohydrate metabolism").}
#'   \item{Level_3}{Character. Specific metabolic subcategory (e.g., "Central carbohydrate metabolism").}
#'   \item{Module_Entry}{Character. KEGG Module ID (e.g., "M00001").}
#'   \item{Module_Name}{Character. Full name of the metabolic module (e.g., "Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate").}
#'   \item{Definition}{Character. KO composition of the module, truncated in display (e.g., "(K00844,K12407,...)").}
#'   \item{Orthology_Symbol}{Character. Short symbol for the KO (e.g., "HK" for hexokinase).}
#'   \item{Orthology_Name}{Character. Full enzyme name with EC number (e.g., "hexokinase \code{[EC:2.7.1.1]}").}
#'   \item{KO_Symbol}{Character. Combined KO ID and symbol (e.g., "K00844; HK").}
#' }
#' @export
"KO_pathway_ref"


#' KEGG Orthology (KO) Abundance/Presence Across Microbial Samples or Genomes
#'
#' A test dataset (wide-format) showing the relative abundance of KEGG Orthology (KO) entries
#' across multiple microbial samples. Values represent normalized abundance metrics
#' (e.g., TPM, RPKM, relative percentage, presence/absence).
#'
#' @format A data frame with 2495 rows (KO entries) and 5 variables:
#' \describe{
#'   \item{KO}{Character. KEGG Orthology ID (e.g., "K00001").}
#'   \item{Marinobacter salarius}{Numeric. Abundance in Genome "Marinobacter salarius".}
#'   \item{Pseudooceanicola nanhaiensis}{Numeric. Abundance in Genome "Pseudooceanicola nanhaiensis".}
#'   \item{Alteromonas australica}{Numeric. Abundance in Genome "Alteromonas australica".}
#'   \item{Henriciella pelagia}{Numeric. Abundance in Genome "Henriciella pelagia".}
#' }
#' @export
"KO_Sample_wide"

.reference_ko_data <- function() {
  list(
    ref1 = KO_pathway_ref,
    ref2 = KO_Sample_wide
  )
  invisible(NULL)
}
