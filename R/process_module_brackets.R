#' @title Process Module Brackets Recursively
#'
#' @description Recursively processes nested brackets in module definitions to calculate
#'              pathway completeness or abundance scores. Handles complex pathway structures
#'              with multiple nesting levels.
#'
#' @param module_abundance Data frame containing KO abundance data for the module
#' @param module_steps_str String representation of module steps/structure
#' @param bracket_count Counter for tracking nested bracket levels
#' @param step_count Counter for tracking processing steps
#' @param module_name Name of the module being processed
#' @param raw_module_steps Original unprocessed module steps string
#' @param plus_scale_method Scaling method for plus-separated KOs ("mean", "min", or "max")
#' @param comma_scale_method Scaling method for comma-separated KOs ("sum" or "max")
#'
#' @return Data frame containing processed abundance values with:
#'         - Rows for each bracket level and final step
#'         - Consistent sample columns as input
#' @export
#'
#' @examples
#' \dontrun{
#' process_module_brackets(module_abundance = sub_Sample_KO_pathway,
#'                         module_steps_str = module_steps_str,
#'                         bracket_count = 1,
#'                         step_count = 1,
#'                         module_name = 'Module',
#'                         raw_module_steps = module_steps_str,
#'                         plus_scale_method,
#'                         comma_scale_method)
#' }
process_module_brackets <- function(
    module_abundance = sub_Sample_KO_pathway,
    module_steps_str = module_steps_str,
    bracket_count = 1,
    step_count = 1,
    module_name = 'Module',
    raw_module_steps = module_steps_str,
    plus_scale_method,
    comma_scale_method) {

  original_module_steps = raw_module_steps
  # Extract the innermost brackets
  brackets <- extract_inner_brackets(module_steps_str)
  # If no brackets are found, return the original string
  if (length(brackets) == 0) {
    cat('\n\tStart processing final step...\n\n')
    cat(paste0('\t\tAnalyzing ',module_name,': ',module_steps_str,'\n'))
    sample_module_steps = process_module_step(module_abundance,
                                              module_steps_str,
                                              aggregrate_rowname = module_name,
                                              step_count = bracket_count,
                                              plus_scale_method,
                                              comma_scale_method)[['abundance_table']]
    return((sample_module_steps))
  } else {
    cat(paste0('\n\tNested steps include: ',brackets,'\n\n'))
    cat(paste0('\tStart processing nested steps...\n\n'))
  }
  # Process each pair of brackets
  for (bracket in brackets) {
    bracket_name = paste0(module_name, '_', bracket_count)
    cat(paste0('\t\tAnalyzing ',bracket_name,': ',bracket,'\n'))
    cat(paste0('\t\tAnalyzing ',': ',bracket_count,'\n'))
    module_steps_str <- stringr::str_replace(module_steps_str, escape_special_chars(bracket), bracket_name)

    result = process_module_step(module_abundance, bracket, aggregrate_rowname = bracket_name, step_count = 1,
                                 plus_scale_method,comma_scale_method)
    sample_module_bracket = result[['abundance_table']]
    module_abundance = rbind(module_abundance, sample_module_bracket) %>% unique(.)

    bracket_count = bracket_count + 1
    cat('\n')
  }
  return(process_module_brackets(module_abundance = module_abundance,
                                 module_steps_str = module_steps_str,
                                 bracket_count = bracket_count,
                                 step_count = step_count,
                                 module_name = module_name,
                                 raw_module_steps = raw_module_steps,
                                 plus_scale_method = plus_scale_method,
                                 comma_scale_method = comma_scale_method))
}
