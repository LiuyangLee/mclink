#' @title Escape Special Characters in a String
#'
#' @description Escapes all specified special characters in a string by adding a backslash before them.
#'              This is particularly useful for preparing strings for use in regular expressions or
#'              other contexts where special characters need to be treated as literals.
#'
#' @param s A character string to be processed
#'
#' @return A new string with all specified special characters escaped with backslashes
#' @export
escape_special_chars <- function(s) {
  # Define all special characters that need to be escaped
  special_chars <- c("(", ")", "+", "-")

  # Escape each special character by adding a backslash before it
  for (char in special_chars) {
    # Use fixed=TRUE for literal matching (faster than regex)
    # paste0 creates the escape sequence (e.g., "(" becomes "\\(")
    s <- gsub(char, paste0("\\", char), s, fixed = TRUE)
  }

  return(s)
}
