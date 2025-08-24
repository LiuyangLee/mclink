#' @title Process All Modules in Pathway Information
#'
#' @description Processes all metabolic modules in pathway information, handling each module's structure,
#'              definition, and bracket components. Aggregates results across all modules.
#' @param pathway_infor Data frame containing pathway information, see examples.
#' @param Sample_KO Data frame containing KO (KEGG Orthology) sample data
#' @param plus_scale_method Scaling method for plus-separated KOs ("mean", "min", or "max")
#' @param comma_scale_method Scaling method for comma-separated KOs ("sum" or "max")
#' @param verbose Logical controlling console output:
#'        \itemize{
#'          \item \code{TRUE} (default): Print progress messages
#'          \item \code{FALSE}: Silent mode
#'        }
#' @return A list with two components:
#' \itemize{
#'   \item `data`: A data frame of processed results for all modules, with unique rows to avoid duplicates.
#'   \item `log`: A character vector of timestamped log messages.
#' }
#' @export
process_all_modules <- function(pathway_infor, Sample_KO,
                                plus_scale_method,
                                comma_scale_method,
                                verbose = TRUE) {

  # Initialize empty result dataframe with proper structure
  result <- data.frame()
  Module_log <- list()

  log_message <- function(msg) {
    timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
    entry <- paste(timestamp, msg)
    if (verbose) message(entry)
    return(entry)  # 返回标准化日志条目
  }

  # Process each module sequentially
  for (each_module in unique(pathway_infor[, 'Module_Entry'])) {

    Module_log <- c(Module_log, log_message(paste("Starting Module:", each_module)))

    # Process module structure to get KO-sample relationships
    sub_Sample_KO_pathway <- process_module_structure(
      pathway_infor = pathway_infor,
      Sample_KO = Sample_KO,
      module = each_module
    )
    # Check if module contains any KOs
    if (dim(sub_Sample_KO_pathway)[2] > 0) {
      # Process module definition string
      module_def_list <- process_module_definition(sub_Sample_KO_pathway)
      module_definition = module_def_list[["definition"]]
      module_definition_log = module_def_list[["log"]]
      Module_log <- c(Module_log, module_definition_log)

      # Process module brackets and get results
      brackets_list <- process_module_brackets(
        module_abundance = sub_Sample_KO_pathway,
        module_steps_str = module_definition,
        bracket_count = 1,
        step_count = 1,
        module_name = each_module,
        raw_module_steps = module_definition,
        plus_scale_method = plus_scale_method,
        comma_scale_method = comma_scale_method,
        abundance_log = list()
      )
      brackets_list.tmp = brackets_list[['data']]
      brackets_list.log = brackets_list[['log']]

      # Log successful processing
      #cat(paste0('\n[', format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      #          '] Module ', each_module, ' finished: ',
      #          unique(sub_Sample_KO_pathway$Definition), '\n'))
      Module_log <- c(Module_log, brackets_list.log, log_message(paste("Completed Module:", each_module)))

    } else {
      # Log warning for modules with no detected KOs
      #cat(paste0('\n[', format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      #          '] Warning: Orthology Entry of ', each_module,
      #          ' cannot be detected in this KO-Sample file!\n'))
      Module_log <- c(Module_log, log_message(paste("No KOs detected in module:", each_module)))
      brackets_list.tmp = data.frame()
    }
    # Combine results while removing duplicates
    result <- rbind(result, brackets_list.tmp) %>% unique(.)
  }

  return(list(
    data = result,
    log = Module_log
  ))
}
