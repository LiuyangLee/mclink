#' @title Read and Process KO Sample Table with Pathway Information
#'
#' @description This function reads a KO sample table (wide format), validates its contents against
#' pathway information, and returns a filtered table containing only KOs present in both datasets.
#'
#' @param in_KO_Sample_wide Character string specifying the path to the input KO sample table file.
#'        Should be a tab-delimited file with KO identifiers as row names.
#' @param pathway_infor Data frame containing pathway information, see examples.
#'
#' @return A data frame containing the filtered sample data, with only rows that match KOs in the
#'         pathway information. Includes an added 'Orthology_Entry' column containing the row names.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' read_and_process_KO_table(in_KO_Sample_wide, pathway_infor)
#' }
read_and_process_KO_table <- function(in_KO_Sample_wide, pathway_infor) {
  # Use tryCatch function to handle potential errors
  tryCatch({
    # Read the file
    # Sample_KO = read.csv(in_KO_Sample_wide, sep='\t',header = TRUE, row.names = 1, check.names = F)
    # Sample_KO = data.table::fread(in_KO_Sample_wide, sep = "\t", header = TRUE, check.names = FALSE, data.table = FALSE)
    Sample_KO = in_KO_Sample_wide
    rownames(Sample_KO) = Sample_KO[,1]
    Sample_KO = Sample_KO[,-1]

    # Check if the file was successfully imported
    if (!is.null(Sample_KO)) {
      cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),'] File successfully imported.\n\n'))
    } else {
      stop(paste0('\n\n[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),'] *** File import failed, please check the input file path and format ***\n\n'))
    }

    KO_list <- unique(pathway_infor$Orthology_Entry)

    # Check the intersection of the first column and KO_list
    common_genes <- intersect(rownames(Sample_KO), KO_list)

    if (length(common_genes) > 0) {
      cat(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),'] There is an intersection between the first column of input file and KO list, with ',length(common_genes),' KOs:\n\n'))
      cat(common_genes)
      cat('\n\n')
    } else {
      stop(paste0('[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),'] *** There is no intersection between the first column of the input file and KO list, please check the input file path and format ***\n\n'))
    }

    # Use the mutate function to add a new column Orthology_Entry, the value of which is the row names of the data frame
    Sample_KO <- dplyr::mutate(Sample_KO, Orthology_Entry = rownames(Sample_KO))
    # Filter out the rows where Orthology_Entry is in pathway_infor$Orthology_Entry
    Sample_KO <- Sample_KO[rownames(Sample_KO) %in% unique(pathway_infor$Orthology_Entry), , drop = F]

    return(Sample_KO)
  }, error = function(e) {
    stop(paste0('\n\n[',format(Sys.time(), "%Y-%m-%d %H:%M:%S"),'] *** File import failed, please check the input file path and format ***\n\n'))
  })
}
