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
#' @param abundance_log A character vector of timestamped log messages
#' @return A list with two components:
#' \itemize{
#'   \item `data`: Data frame containing processed abundance values with:
#'         - Rows for each bracket level and final step
#'         - Consistent sample columns as input
#'   \item `log`: A character vector of timestamped log messages.
#' }
#' @export
process_module_brackets <- function(
    module_abundance = sub_Sample_KO_pathway,
    module_steps_str = module_steps_str,
    bracket_count = 1,
    step_count = 1,
    module_name = 'Module',
    raw_module_steps = module_steps_str,
    plus_scale_method,
    comma_scale_method,
    abundance_log = list()) {

  timestamp <- function() format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  log_entry <- function(msg) {
    paste0(timestamp(), " ", msg)
  }
  original_module_steps = raw_module_steps
  # Extract the innermost brackets
  brackets <- extract_inner_brackets(module_steps_str)
  # If no brackets are found, return the original string
  if (length(brackets) == 0) {
    #cat('\n\tStart processing final step...\n')
    #cat(paste0('\t\tAnalyzing ',module_name,': ',module_steps_str))
    abundance_log <- c(abundance_log, list(
      log_entry("Start processing final step..."),
      log_entry(paste0("Analyzing ", module_name, ": ", module_steps_str))
    ))
    module_step_list = process_module_step(module_abundance,
                                           module_steps_str,
                                           aggregrate_rowname = module_name,
                                           step_count = bracket_count,
                                           plus_scale_method,
                                           comma_scale_method
                                           )
    abundance_log <- c(abundance_log, module_step_list[['abundance_log']])
    return(list(
      data = module_step_list[['abundance_table']],
      log = abundance_log
    ))
  } else {
    #cat(paste0('\n\tNested steps include: ',brackets))
    #cat(paste0('\tStart processing nested steps...\n'))
    abundance_log <- c(abundance_log, list(
      log_entry(paste0("Nested steps include: ")),
      log_entry(paste0('    ',brackets)),
      log_entry("Start processing nested steps...")
    ))
  }
  # Process each pair of brackets
  for (bracket in brackets) {
    bracket_name = paste0(module_name, '_', bracket_count)
    #cat(paste0('\t\tAnalyzing ',bracket_name,': ',bracket))
    #cat(paste0('\t\tAnalyzing ',': ',bracket_count))
    abundance_log <- c(abundance_log, list(
        log_entry(paste0("Analyzing ", bracket_name, ": ", bracket)),
        log_entry(paste0("Bracket level: ", bracket_count))
      ))
    module_steps_str <- stringr::str_replace(module_steps_str, escape_special_chars(bracket), bracket_name)

    module_step_list = process_module_step(module_abundance, bracket, aggregrate_rowname = bracket_name, step_count = 1,
                                 plus_scale_method,comma_scale_method)
    sample_module_bracket = module_step_list[['abundance_table']]
    module_step_list.log = module_step_list[['abundance_log']]
    abundance_log <- c(abundance_log, module_step_list.log)
    module_abundance = rbind(module_abundance, sample_module_bracket) %>% unique(.)
    bracket_count = bracket_count + 1
  }
  return(process_module_brackets(module_abundance = module_abundance,
                                 module_steps_str = module_steps_str,
                                 bracket_count = bracket_count,
                                 step_count = step_count,
                                 module_name = module_name,
                                 raw_module_steps = raw_module_steps,
                                 plus_scale_method = plus_scale_method,
                                 comma_scale_method = comma_scale_method,
                                 abundance_log = abundance_log))
}
