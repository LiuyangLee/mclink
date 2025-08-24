#' @title Read and process Pathway information dataframe
#'
#' @description This function reads a tab-delimited file containing KEGG pathway information,
#' performs data validation and cleaning, and returns a processed data frame.
#'
#' @param in_KO_pathway_ref Character string specifying the path to the input file.
#'        The file should be a tab-delimited text file containing KEGG pathway information,
#'        with a header row and at least one column named "Module_Entry".
#'
#' @return A list with log and a data frame containing the processed pathway information after removing
#'         empty/NA entries. The data frame will have the same columns as the input file.
#'         Returns NULL if the file cannot be read.
#'
#' @export
read_and_process_pathway_infor <- function(in_KO_pathway_ref) {
  # Use tryCatch function to handle potential errors
  tryCatch({
    # Read the file
    #pathway_infor = read.csv(in_KO_pathway_ref, sep='\t',header = TRUE, check.names = F)
    pathway_infor = in_KO_pathway_ref

    timestamp <- function() format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
    log_entry <- function(msg) {
      paste0(timestamp(), " ", msg)
    }

    # Check if the file was successfully imported
    if (!is.null(pathway_infor)) {
      message(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),'] Pathway information dataframe successfully imported.'))
      pathway_log <- list(
          log_entry("Pathway information dataframe successfully imported.")
        )
    } else {
      pathway_log <- list(
          log_entry("*** Pathway information dataframe is null, please check the input file path and format ***")
        )
      stop(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),'] *** Pathway information dataframe is null, please check the input file path and format ***\n'))
    }
    # 删除列a中含有NA的行
    pathway_infor <- pathway_infor[!(pathway_infor$Module_Entry %in% 'NA'), , drop = F]
    pathway_infor <- pathway_infor[!(pathway_infor$Module_Entry %in% ''), , drop = F]
    pathway_infor <- pathway_infor[!is.na(pathway_infor$Module_Entry), , drop = F]

    Module_list <- unique(pathway_infor$Module_Entry)
    if (length(Module_list) > 0) {
      message(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),'] There are ',length(Module_list),' Modules in the Pathway information dataframe.'))
      pathway_log <- c(pathway_log, list(
          log_entry(paste0('There are ',length(Module_list),' Modules in the Pathway information dataframe: ', paste(Module_list, collapse = ' ')))
        ))
    } else {
      pathway_log <- c(pathway_log, list(
          log_entry('*** No Module detected, please check the input file path and format ***')
        ))
      stop(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),'] *** No Module detected, please check the input file path and format ***\n'))
    }
    return(list(data = pathway_infor, log = pathway_log))
  }, error = function(e) {
    stop(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),'] *** Pathway information dataframe import failed, please check the input file path and format ***\n'))
  })
}
