#' @title Calculate Log2 Ratio of Sample Values to Row Means
#'
#' @description This function calculates the log2 ratio of each value in a data frame
#' to its corresponding row mean. Zero values are preserved as zeros in the output.
#'
#' @param data A data frame or matrix containing numerical values to be processed.
#'             Rows represent features (e.g., KO terms) and columns represent samples.
#'             Default is KO_Sample_table.
#'
#' @return A data frame of the same dimensions as input, where each value is:
#'         - log2(sample_value/row_mean) when both sample_value and row_mean are non-zero
#'         - 0 when either sample_value or row_mean is zero
#'         The Mean_RA column used for calculations is removed from the output.
#'
#' @export
ata_cal = function(data = KO_Sample_table){
  data$Mean_RA = base::apply(data, 1, mean)
  scaled_df = data
  for(col in colnames(data)){
    scaled_df[[col]] <- ifelse(((data[[col]] != 0) & (data$Mean_RA != 0)), log2(data[[col]]/data$Mean_RA), 0)
  }
  scaled_df = dplyr::select(scaled_df, -Mean_RA)
  return(scaled_df)
}
