#' @title Remove Outer Parentheses from String
#'
#' @description This function checks if a string is wrapped in outer parentheses
#'              and removes them if present.
#'
#' @param s A character string to be processed
#'
#' @return The input string with outer parentheses removed (if they existed),
#'         or the original string if no outer parentheses were found.
#' @export
#'
#' @examples
#' \dontrun{
#' remove_outer_brackets(s)
#' }
remove_outer_brackets <- function(s) {
  # Check if the string starts with '(' and ends with ')' using regex
  if (base::grepl("^\\(.*\\)$", s)) {
    # If outer parentheses exist, remove them using regex substitution
    # ^ matches start of string, $ matches end
    # \\( and \\) escape the parentheses
    # (.*) captures everything inside the parentheses
    s <- base::sub("^\\((.*)\\)$", "\\1", s)
  }
  return(s)
}
