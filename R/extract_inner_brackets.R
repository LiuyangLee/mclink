#' @title Extract Innermost Parentheses Content
#'
#' @description Extracts all text segments enclosed in the innermost level of parentheses from a string.
#'              This is useful for parsing hierarchical or nested parenthetical expressions.
#'
#' @param s A character string to process (can contain multiple parenthetical groups)
#'
#' @return A character vector containing all innermost parenthesized segments
#'         Returns empty character vector if no matches found
#' @export
extract_inner_brackets <- function(s) {
  # Regex explanation:
  # "\\(    " - Matches literal opening parenthesis
  # "[^()]* " - Matches any character except ( or ) zero or more times
  # "\\)    " - Matches literal closing parenthesis
  # This ensures we only match innermost parentheses (no nested cases)
  matches <- stringr::str_extract_all(s, "\\([^()]*\\)")[[1]]

  return(matches)
}
