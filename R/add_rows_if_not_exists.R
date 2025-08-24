#' @title Add Missing Rows to a Data Frame
#'
#' @description This function checks if specified rows exist in a data frame, and if not,
#'              adds them with all values set to 0. Useful for ensuring consistent KO
#'              representation across samples.
#'
#' @param module_abundance A data frame with row names representing KO identifiers
#'                         (e.g., K numbers) and numeric abundance values
#' @param add_rows A character vector of row names (KO identifiers) that should be
#'                 present in the output. Defaults to c("K14126","K14128","K14127")
#'
#' @return The original data frame with additional rows (if any were missing) where
#'         all values are set to 0. Row and column names are preserved.
#' @export
add_rows_if_not_exists <- function(module_abundance, add_rows = c("K14126","K14128","K14127")) {
  missing_rows <- base::setdiff(add_rows, rownames(module_abundance))
  absent_genome_KOs <- data.frame(
    matrix(0,
           nrow = length(missing_rows),
           ncol = ncol(module_abundance))
  )
  colnames(absent_genome_KOs) <- colnames(module_abundance)
  rownames(absent_genome_KOs) <- missing_rows
  result <- rbind(module_abundance, absent_genome_KOs)
  return(result)
}
