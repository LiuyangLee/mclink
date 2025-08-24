#' @title Process Module Steps with Complex KO String Handling
#'
#' @description Processes a KO string containing various combinations of KOs separated by
#'              different operators (commas, plus signs, or spaces). Handles complex pathway
#'              definitions with multiple types of relationships between KOs.
#'
#' @param module_abundance Data frame containing KO abundance data
#' @param KO_string String representation of KO relationships (default: "K03388,K03389+K03390+K14083,K14126+K14127,K14128")
#' @param aggregrate_rowname Base name for row aggregation (default: 'bracket_1')
#' @param step_count Counter for processing steps (default: 1)
#' @param plus_scale_method Scaling method for plus-separated KOs ("mean", "min", or "max")
#' @param comma_scale_method Scaling method for comma-separated KOs ("sum" or "max")
#'
#' @return List containing:
#'         - abundance_table: Processed abundance values
#'         - step_count: Updated step counter
#'         - abundance_log: log
#' @export
process_module_step <- function(module_abundance,
                                KO_string = "K03388,K03389+K03390+K14083,K14126+K14127,K14128",
                                aggregrate_rowname = 'bracket_1',
                                step_count = 1,
                                plus_scale_method,
                                comma_scale_method) {
  # Remove outer brackets if present
  KO_string = remove_outer_brackets(KO_string)

  # Check for different types of separators in the KO string
  has_comma <- base::grepl(",", KO_string)
  has_plus <- base::grepl("\\+", KO_string)
  has_space <- base::grepl(" ", KO_string)
  abundance_log.tmp = list()

  # Case 1: String contains plus, comma, and space (most complex case)
  # (K14126+K14128,K22516+K00125 K00126 K22516,K14126 K22516+K00125)
  if ((has_plus && has_comma && has_space)) {
    # Split by space first, then handle each component
    # (K14126+K14128,K22516+K00125 K00126) to (K14126+K14128,K22516+K00125) (K00126)
    KO_vector <- base::strsplit(KO_string, " ")[[1]]
    loop_plu_comma_list <- process_module_loop_plu_comma(
      KO_vector,
      module_abundance,
      process_step_plus,
      process_step_comma,
      process_step_direct,
      aggregrate_rowname,
      step_count,
      plus_scale_method,
      comma_scale_method
    )
    abundance_log.tmp = loop_plu_comma_list[['abundance_log']]
    KO_scale = paste(rownames(loop_plu_comma_list[['abundance_table']]), collapse = '+')
    result_list <- process_step_plus(
      loop_plu_comma_list[['abundance_table']],
      KO_scale,
      aggregrate_rowname,
      loop_plu_comma_list[['step_count']],
      plus_scale_method
    )
  }
  # Case 2: String contains both plus and comma (but no space)
  else if (has_plus && has_comma) {
    # Split by comma first, then handle each component
    KO_vector <- base::strsplit(KO_string, ",")[[1]]
    loop_plus_list <- process_module_loop_plus(
      KO_vector,
      module_abundance,
      process_step_plus,
      process_step_direct,
      aggregrate_rowname,
      step_count,
      plus_scale_method
    )
    abundance_log.tmp = loop_plus_list[['abundance_log']]
    KO_scale = paste(rownames(loop_plus_list[['abundance_table']]), collapse = ',')
    result_list <- process_step_comma(
      loop_plus_list[['abundance_table']],
      KO_scale,
      aggregrate_rowname,
      loop_plus_list[['step_count']],
      comma_scale_method
    )
  }
  # Case 3: String contains plus and space (but no comma)
  else if (has_plus && has_space) {
    # Split by space first, then handle each component
    KO_vector <- base::strsplit(KO_string, " ")[[1]]
    loop_plus_list <- process_module_loop_plus(
      KO_vector,
      module_abundance,
      process_step_plus,
      process_step_direct,
      aggregrate_rowname,
      step_count,
      plus_scale_method
    )
    abundance_log.tmp = loop_plus_list[['abundance_log']]
    KO_scale = paste(rownames(loop_plus_list[['abundance_table']]), collapse = ' ')
    result_list <- process_step_space(
      loop_plus_list[['abundance_table']],
      KO_scale,
      aggregrate_rowname,
      loop_plus_list[['step_count']]
    )
  }
  # Case 4: String contains comma and space (but no plus)
  else if (has_comma && has_space) {
    # Split by space first, then handle each component
    KO_vector <- base::strsplit(KO_string, " ")[[1]]
    loop_comma_list <- process_module_loop_comma(
      KO_vector,
      module_abundance,
      process_step_comma,
      process_step_direct,
      aggregrate_rowname,
      step_count,
      comma_scale_method
    )
    abundance_log.tmp = loop_comma_list[['abundance_log']]
    KO_scale = paste(rownames(loop_comma_list[['abundance_table']]), collapse = ' ')
    result_list <- process_step_space(
      loop_comma_list[['abundance_table']],
      KO_scale,
      aggregrate_rowname,
      loop_comma_list[['step_count']]
    )
  }
  # Case 5: String contains only commas (K14126,K22516)
  else if (has_comma) {
    result_list <- process_step_comma(
      module_abundance,
      KO_string,
      aggregrate_rowname,
      step_count,
      comma_scale_method
    )
  }
  # Case 6: String contains only plus signs (K14126+K14128)
  else if (has_plus) {
    result_list <- process_step_plus(
      module_abundance,
      KO_string,
      aggregrate_rowname,
      step_count,
      plus_scale_method
    )
  }
  # Case 7: String contains only spaces (K14126 K22516)
  else if (has_space) {
    result_list <- process_step_space(
      module_abundance,
      KO_string,
      aggregrate_rowname,
      step_count
    )
  }
  # Case 8: Single KO with no separators (K14126)
  else {
    result = process_step_direct(module_abundance, KO_string)
    abundance_table = result[['abundance_table']]
    abundance_log = result[['abundance_log']]
    result_list = list(abundance_table = abundance_table, step_count = step_count, abundance_log = abundance_log)
    #stop("Please check Module Definition!")
  }
  result_list = list(abundance_table = result_list[['abundance_table']],
                     step_count = result_list[['step_count']],
                     abundance_log = c(abundance_log.tmp ,result_list[['abundance_log']]))
  return(result_list)
}
