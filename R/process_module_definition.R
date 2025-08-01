#' @title Process Module Definition String
#'
#' @description Cleans and processes module definition strings by removing various patterns and formatting
#'              elements to extract core KO relationships. Handles special cases like negative KOs and
#'              parenthetical expressions.
#'
#' @param sub_Sample_KO_pathway Data frame containing module definitions in a 'Definition' column
#'
#' @return Character vector of cleaned module definition strings with:
#'         - Removed negative KO indicators
#'         - Simplified parentheses
#'         - Normalized space
#' @export
#'
#' @examples
#' \dontrun{
#' process_module_definition(sub_Sample_KO_pathway)
#' }
process_module_definition <- function(sub_Sample_KO_pathway) {
  # Use unique() to remove duplicate values in Definition column
  unique_definitions <- unique(sub_Sample_KO_pathway$Definition)
  # Use gsub() to remove specific patterns from strings:
  # Remove "-- " pattern
  step1 <- gsub("-- ", "", unique_definitions)
  # Remove " --" pattern
  step2 <- gsub(" --", "", step1)
  # Remove "-Knumber" pattern (negative KO indicators)
  step3 <- gsub("-K\\d+", "", step2)
  # Handle negative parenthetical case: M00011 -(K00242,K18859,K18860)
  step4 <- gsub("-\\(.*?\\)", "", step3)
  # Replace multiple consecutive spaces with single space
  step5 <- gsub("\\s+", " ", step4)
  # Replace "(Knumber)" pattern with "Knumber" (simplify parentheses)
  result <- gsub("\\((K\\d+)\\)", "\\1", step5)
  cat(paste0('\nProcessing module steps: ', unique_definitions, '\n\n'))
  cat(paste0('After omitting minus KOs: ', result, '\n\n'))
  return(result)
}
