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
#' @param n_cores Number of CPU cores for module-level parallel processing.
#'        Default \code{1} (sequential). Values greater than 1 use
#'        \code{parallel::mclapply()}, which is only effective on Unix-like
#'        systems; on Windows processing falls back to sequential.
#' @return A list with two components:
#'   \itemize{
#'     \item `data`: A data frame of processed results for all modules, with unique rows to avoid duplicates.
#'     \item `log`: A character vector of timestamped log messages.
#'   }
#' @export
process_all_modules <- function(pathway_infor, Sample_KO,
                                plus_scale_method,
                                comma_scale_method,
                                verbose = TRUE,
                                n_cores = 1) {

  timestamp <- function() format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  log_message <- function(msg) {
    entry <- paste(timestamp(), msg)
    if (verbose) message(entry)
    return(entry)
  }

  # Numeric matrix view of the sample table (single conversion, shared by all modules)
  sample_cols <- setdiff(colnames(Sample_KO), "Orthology_Entry")
  Sample_mat <- as.matrix(Sample_KO[, sample_cols, drop = FALSE])
  rn <- rownames(Sample_KO)
  if (is.null(rn)) rn <- Sample_KO$Orthology_Entry
  rownames(Sample_mat) <- rn

  modules <- unique(pathway_infor[, 'Module_Entry'])
  n_modules <- length(modules)
  t_start <- Sys.time()

  run_one_module <- function(i) {
    each_module <- modules[i]
    module_log <- list(log_message(paste("Starting Module:", each_module)))

    # KOs annotated to this module, matched against the sample table
    each_pathway_infor <- unique(pathway_infor[pathway_infor$Module_Entry %in% each_module,
                                               c("Orthology_Entry", "Module_Entry", "Definition"),
                                               drop = FALSE])
    idx <- match(unique(each_pathway_infor$Orthology_Entry), rownames(Sample_mat), nomatch = 0L)
    idx <- idx[idx > 0L]

    if (length(idx) > 0) {
      sub_mat <- Sample_mat[idx, , drop = FALSE]
      # Definition metadata of the detected KOs only (mirrors the old merge-based path)
      def_source <- each_pathway_infor[each_pathway_infor$Orthology_Entry %in% rownames(sub_mat), , drop = FALSE]

      module_def_list <- process_module_definition(
        data.frame(Definition = def_source$Definition, stringsAsFactors = FALSE)
      )
      module_definition <- module_def_list[["definition"]]
      module_log <- c(module_log, module_def_list[["log"]])

      brackets_list <- .brackets_core(
        m = sub_mat,
        module_steps_str = module_definition,
        bracket_count = 1,
        step_count = 1,
        module_name = each_module,
        raw_module_steps = module_definition,
        plus_scale_method = plus_scale_method,
        comma_scale_method = comma_scale_method,
        abundance_log = list()
      )
      module_log <- c(module_log, brackets_list[['log']],
                      list(log_message(paste("Completed Module:", each_module))))

      module_df <- as.data.frame(brackets_list[['data']])
      module_df$Orthology_Entry <- rownames(module_df)
      module_df$Module_Entry <- unique(def_source$Module_Entry)
      module_df$Definition <- unique(def_source$Definition)
      list(data = module_df, log = module_log)
    } else {
      module_log <- c(module_log, list(log_message(paste("No KOs detected in module:", each_module))))
      list(data = NULL, log = module_log)
    }
  }

  use_parallel <- n_cores > 1 && n_modules > 1
  if (use_parallel && .Platform$OS.type == "windows") {
    message(paste(timestamp(), "Parallel processing is not supported on Windows; running sequentially."))
    use_parallel <- FALSE
  }
  if (use_parallel) {
    message(paste(timestamp(), paste0("Processing ", n_modules, " modules using ",
                                      min(n_cores, n_modules), " cores...")))
    per_module <- parallel::mclapply(seq_len(n_modules), run_one_module,
                                     mc.cores = min(n_cores, n_modules))
    message(paste(timestamp(), sprintf("All modules done (elapsed: %.1fs).",
                                       as.numeric(difftime(Sys.time(), t_start, units = "secs")))))
  } else {
    per_module <- vector("list", n_modules)
    progress_every <- max(1L, floor(n_modules / 10))
    for (i in seq_len(n_modules)) {
      per_module[[i]] <- run_one_module(i)
      if (i %% progress_every == 0 || i == n_modules) {
        message(paste(timestamp(),
                      sprintf("Module progress: %d/%d (%d%%, elapsed: %.1fs)",
                              i, n_modules, floor(100 * i / n_modules),
                              as.numeric(difftime(Sys.time(), t_start, units = "secs")))))
      }
    }
  }

  result_list <- lapply(per_module, function(x) x[['data']])
  result_list <- result_list[!vapply(result_list, is.null, logical(1))]
  Module_log <- do.call(c, lapply(per_module, function(x) x[['log']]))

  result <- if (length(result_list) > 0) unique(do.call(rbind, result_list)) else data.frame()

  return(list(
    data = result,
    log = Module_log
  ))
}
