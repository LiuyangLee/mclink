% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_step_plus.R
\name{process_step_plus}
\alias{process_step_plus}
\title{Process Plus-Separated KOs with Specified Scaling Method}
\usage{
process_step_plus(
  module_abundance,
  KOs = c("K14126+K14127+K14128"),
  aggregrate_rowname,
  step_count = 1,
  plus_scale_method
)
}
\arguments{
\item{module_abundance}{Data frame containing KO abundance data with required columns:
Orthology_Entry, Module_Entry, Definition}

\item{KOs}{Character string of plus-separated KO IDs (default: "K14126+K14127+K14128")}

\item{aggregrate_rowname}{Base name for row aggregation (default: 'step_1')}

\item{step_count}{Processing step counter (default: 1)}

\item{plus_scale_method}{Scaling method for plus-separated KOs ("mean", "min", or "max")}
}
\value{
List containing:
- abundance_table: Processed data with aggregated values
- step_count: Updated step counter
}
\description{
Handles plus-separated KOs by applying the specified scaling method (mean, min, or max).
Processes multiple KOs separated by plus signs and aggregates them into a single row.
Note: For mean calculation, uses the sum of existing KO abundances divided by total
number of KOs (including those with zero abundance in all samples).
}
\examples{
\dontrun{
process_step_plus(module_abundance, KOs = c("K14126+K14127+K14128"), aggregrate_rowname,
                            step_count = 1, plus_scale_method)
}
}
