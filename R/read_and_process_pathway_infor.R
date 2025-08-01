#' @title Read and process pathway information file
#'
#' @description This function reads a tab-delimited file containing KEGG pathway information,
#' performs data validation and cleaning, and returns a processed data frame.
#'
#' @param in_KO_pathway_ref Character string specifying the path to the input file.
#'        The file should be a tab-delimited text file containing KEGG pathway information,
#'        with a header row and at least one column named "Module_Entry".
#'
#' @return A data frame containing the processed pathway information after removing
#'         empty/NA entries. The data frame will have the same columns as the input file.
#'         Returns NULL if the file cannot be read.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' read_and_process_pathway_infor(in_KO_pathway_ref)
#' }
read_and_process_pathway_infor <- function(in_KO_pathway_ref) {
  # Use tryCatch function to handle potential errors
  tryCatch({
    # Read the file
    #pathway_infor = read.csv(in_KO_pathway_ref, sep='\t',header = TRUE, check.names = F)
    pathway_infor = in_KO_pathway_ref

    # Check if the file was successfully imported
    if (!is.null(pathway_infor)) {
      cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),'] Pathway information file successfully imported.\n\n'))
    } else {
      stop(paste0('\n\n[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),'] *** Pathway information file is null, please check the input file path and format ***\n\n'))
    }
    # 删除列a中含有NA的行
    pathway_infor <- pathway_infor[!(pathway_infor$Module_Entry %in% 'NA'), , drop = F]
    pathway_infor <- pathway_infor[!(pathway_infor$Module_Entry %in% ''), , drop = F]
    pathway_infor <- pathway_infor[!is.na(pathway_infor$Module_Entry), , drop = F]

    Module_list <- unique(pathway_infor$Module_Entry)
    if (length(Module_list) > 0) {
      cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),'] There are ',length(Module_list),' Modules in the Pathway information file.\n\n'))
      cat(Module_list)
      cat('\n\n')
    } else {
      stop(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),'] *** No Module detected, please check the input file path and format ***\n\n'))
    }

    return(pathway_infor)
  }, error = function(e) {
    stop(paste0('\n\n[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),'] *** Pathway information file import failed, please check the input file path and format ***\n\n'))
  })
}
