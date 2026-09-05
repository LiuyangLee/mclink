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

  # Precompute KO -> row positions once (single hash of the table row names),
  # split by module; definitions restricted to detected KOs, mirroring the old
  # per-module merge of detected rows only.
  ko_module <- unique(pathway_infor[, c("Module_Entry", "Orthology_Entry")])
  ko_pos <- match(ko_module$Orthology_Entry, rn, nomatch = 0L)
  ko_module <- ko_module[ko_pos > 0L, , drop = FALSE]
  ko_pos <- ko_pos[ko_pos > 0L]
  mod_pos <- split(ko_pos, factor(ko_module$Module_Entry, levels = modules))

  ko_def <- unique(pathway_infor[, c("Module_Entry", "Orthology_Entry", "Definition")])
  ko_def <- ko_def[paste(ko_def$Module_Entry, ko_def$Orthology_Entry) %in%
                   paste(ko_module$Module_Entry, ko_module$Orthology_Entry), , drop = FALSE]
  mod_defs <- split(ko_def$Definition, factor(ko_def$Module_Entry, levels = modules))

  run_one_module <- function(i) {
    each_module <- modules[i]
    module_log <- list(log_message(paste("Starting Module:", each_module)))
    pos <- mod_pos[[each_module]]

    if (!is.null(pos) && length(pos) > 0) {
      sub_mat <- Sample_mat[pos, , drop = FALSE]

      module_def_list <- process_module_definition(
        data.frame(Definition = mod_defs[[each_module]], stringsAsFactors = FALSE)
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

      list(data = brackets_list[['data']],
           module_entry = each_module,
           definition = unique(mod_defs[[each_module]]),
           log = module_log)
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

  result_list <- per_module[!vapply(per_module, function(x) is.null(x[['data']]), logical(1))]
  Module_log <- do.call(c, lapply(per_module, function(x) x[['log']]))

  if (length(result_list) > 0) {
    n_rows <- vapply(result_list, function(x) nrow(x[['data']]), integer(1))
    res_mat <- do.call(rbind, lapply(result_list, function(x) x[['data']]))
    scalar_meta <- all(vapply(result_list, function(x) length(x[['definition']]), integer(1)) == 1)

    if (scalar_meta && !anyDuplicated(rownames(res_mat))) {
      # Fast path: row names are unique, so the old data.frame unique() could not
      # drop anything (its Orthology_Entry column equaled the row names); bind
      # matrices at C level and assemble the wide data frame exactly once.
      meta_entry <- rep(vapply(result_list, function(x) x[['module_entry']], character(1)), n_rows)
      meta_def <- rep(vapply(result_list, function(x) x[['definition']], character(1)), n_rows)
      result <- as.data.frame(res_mat)
      result$Orthology_Entry <- rownames(result)
      result$Module_Entry <- meta_entry
      result$Definition <- meta_def
    } else {
      # Fallback: original per-module data.frame assembly + full unique()
      result <- unique(do.call(rbind, lapply(result_list, function(x) {
        df <- as.data.frame(x[['data']])
        df$Orthology_Entry <- rownames(df)
        df$Module_Entry <- x[['module_entry']]
        df$Definition <- x[['definition']]
        df
      })))
    }
  } else {
    result <- data.frame()
  }

  return(list(
    data = result,
    log = Module_log
  ))
}
