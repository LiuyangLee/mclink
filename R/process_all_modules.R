#' @title Process All Modules in Pathway Information
#'
#' @description Processes all metabolic modules in pathway information, handling each module's structure,
#'              definition, and bracket components. Aggregates results across all modules.
#' @param pathway_infor Data frame containing pathway information, see examples.
#' @param Sample_KO Data frame containing KO (KEGG Orthology) sample data
#' @param plus_scale_method Scaling method for plus-separated KOs ("mean", "min", or "max")
#' @param comma_scale_method Scaling method for comma-separated KOs ("sum" or "max")
#'
#' @return A data frame containing processed results for all modules with:
#'         - Combined outputs from process_module_brackets for all modules
#'         - Unique rows to avoid duplicates
#' @export
#'
#' @examples
#' \dontrun{
#' process_all_modules(pathway_infor, Sample_KO,plus_scale_method,comma_scale_method)
#' }
process_all_modules <- function(pathway_infor, Sample_KO,
                                plus_scale_method,
                                comma_scale_method) {

  # Initialize empty result dataframe with proper structure
  result <- data.frame()

  # Process each module sequentially
  for (each_module in unique(pathway_infor[, 'Module_Entry'])) {

    # Process module structure to get KO-sample relationships
    sub_Sample_KO_pathway <- process_module_structure(
      pathway_infor = pathway_infor,
      Sample_KO = Sample_KO,
      module = each_module
    )

    # Check if module contains any KOs
    if (dim(sub_Sample_KO_pathway)[2] > 0) {
      # Process module definition string
      module_definition <- process_module_definition(sub_Sample_KO_pathway)

      # Process module brackets and get results
      result.tmp <- process_module_brackets(
        module_abundance = sub_Sample_KO_pathway,
        module_steps_str = module_definition,
        bracket_count = 1,
        step_count = 1,
        module_name = each_module,
        raw_module_steps = module_definition,
        plus_scale_method = plus_scale_method,
        comma_scale_method = comma_scale_method
      )

      # Log successful processing
      cat(paste0('\n[', format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                '] Module ', each_module, ' finished: ',
                unique(sub_Sample_KO_pathway$Definition), '\n\n'))

    } else {
      # Log warning for modules with no detected KOs
      cat(paste0('\n[', format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                '] Warning: Orthology Entry of ', each_module,
                ' cannot be detected in this KO-Sample file!\n\n'))
      result.tmp <- data.frame()
    }

    # Combine results while removing duplicates
    result <- rbind(result, result.tmp) %>% unique(.)
  }

  return(result)
}
