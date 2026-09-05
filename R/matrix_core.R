# Internal matrix-based cores for the module recursion.
#
# The exported process_step_* / process_module_loop_* / process_module_step /
# process_module_brackets keep their data.frame interfaces unchanged. These cores
# do the same work on a plain numeric matrix (columns are samples, rows are named
# by KO or aggregate name), avoiding repeated data.frame dispatch, rbind factor
# checks and as.matrix() copies on wide tables. process_all_modules() runs the
# fast path through .brackets_core(); the exported wrappers remain for
# compatibility and must stay semantically in sync with these cores.

.step_log_ts <- function() paste0('[', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ']')

# Rows for KOs present in m (in m's row order) plus zero rows for missing KOs,
# mirroring subsetting + add_rows_if_not_exists().
.subset_ko_rows <- function(m, KOs) {
  hit <- rownames(m) %in% KOs
  sub <- m[hit, , drop = FALSE]
  missing_rows <- base::setdiff(KOs, rownames(sub))
  if (length(missing_rows) > 0) {
    zero <- matrix(0, nrow = length(missing_rows), ncol = ncol(m),
                   dimnames = list(missing_rows, colnames(m)))
    sub <- rbind(sub, zero)
  }
  sub
}

# Row-wise de-duplication equivalent to unique() on the old data.frame representation:
# the data frames carried an Orthology_Entry column equal to the row names, so rows
# could only be duplicates when both name and values matched.
.unique_rows <- function(m) {
  m[!duplicated(rownames(m)), , drop = FALSE]
}

.step_plus_core <- function(m, KOs, aggregrate_rowname, step_count = 1, plus_scale_method) {
  KOs_scale <- base::strsplit(KOs, "\\+")[[1]]
  log_messages <- list(paste0(.step_log_ts(), '    ', 'Running KOs plus: ', aggregrate_rowname, " = ", KOs_scale))
  sub <- .subset_ko_rows(m, KOs_scale)
  if (plus_scale_method == "mean") {
    scaled <- t(base::colSums(sub) / length(KOs_scale))
  } else if (plus_scale_method == "min") {
    scaled <- t(matrixStats::colMins(sub))
  } else if (plus_scale_method == "max") {
    scaled <- t(matrixStats::colMaxs(sub))
  } else {
    stop(paste("Unknown plus scale method:", plus_scale_method))
  }
  rownames(scaled) <- aggregrate_rowname
  list(abundance_matrix = scaled, step_count = step_count + 1, abundance_log = log_messages)
}

.step_comma_core <- function(m, KOs, aggregrate_rowname, step_count = 1, comma_scale_method) {
  KOs_scale <- base::strsplit(KOs, ",")[[1]]
  log_messages <- list(paste0(.step_log_ts(), '    ', 'Running KOs comma: ', aggregrate_rowname, " = ", KOs_scale))
  sub <- .subset_ko_rows(m, KOs_scale)
  if (comma_scale_method == "sum") {
    scaled <- t(base::colSums(sub))
  } else if (comma_scale_method == "max") {
    scaled <- t(matrixStats::colMaxs(sub))
  } else {
    stop(paste("Unknown comma scale method:", comma_scale_method))
  }
  rownames(scaled) <- aggregrate_rowname
  list(abundance_matrix = scaled, step_count = step_count + 1, abundance_log = log_messages)
}

.step_space_core <- function(m, KOs, aggregrate_rowname, step_count = 1) {
  KOs_scale <- base::strsplit(KOs, " ")[[1]]
  log_messages <- list(paste0(.step_log_ts(), '    ', 'Running KOs space: ', aggregrate_rowname, " = ", KOs_scale))
  sub <- .subset_ko_rows(m, KOs_scale)
  scaled <- t(base::colSums(sub) / length(KOs_scale))
  rownames(scaled) <- aggregrate_rowname
  list(abundance_matrix = scaled, step_count = step_count + 1, abundance_log = log_messages)
}

.step_direct_core <- function(m, KOs) {
  log_messages <- list(paste0(.step_log_ts(), '    ', 'Running KOs direct: ', KOs))
  sub <- .subset_ko_rows(m, KOs)
  list(abundance_matrix = sub, abundance_log = log_messages)
}

.loop_plus_core <- function(KO_vector, m, aggregrate_rowname = 'step_1', step_count = 1,
                            plus_scale_method) {
  mats <- list()
  log_messages <- list()
  for (KOs in KO_vector) {
    if (base::grepl("\\+", KOs)) {
      result <- .step_plus_core(m, KOs, paste0(aggregrate_rowname, '_', step_count),
                                step_count, plus_scale_method)
      step_count <- result[['step_count']]
    } else {
      result <- .step_direct_core(m, KOs)
    }
    mats[[length(mats) + 1]] <- result[['abundance_matrix']]
    log_messages <- c(log_messages, result[['abundance_log']])
  }
  abundance_matrix <- .unique_rows(do.call(rbind, mats))
  list(abundance_matrix = abundance_matrix, step_count = step_count, abundance_log = log_messages)
}

.loop_comma_core <- function(KO_vector, m, aggregrate_rowname = 'step_1', step_count = 1,
                             comma_scale_method) {
  mats <- list()
  log_messages <- list()
  for (KOs in KO_vector) {
    if (base::grepl(",", KOs)) {
      result <- .step_comma_core(m, KOs, paste0(aggregrate_rowname, '_', step_count),
                                 step_count, comma_scale_method)
      step_count <- result[['step_count']]
    } else {
      result <- .step_direct_core(m, KOs)
    }
    mats[[length(mats) + 1]] <- result[['abundance_matrix']]
    log_messages <- c(log_messages, result[['abundance_log']])
  }
  abundance_matrix <- .unique_rows(do.call(rbind, mats))
  list(abundance_matrix = abundance_matrix, step_count = step_count, abundance_log = log_messages)
}

.loop_plu_comma_core <- function(KO_vector, m, aggregrate_rowname = 'step_1', step_count = 1,
                                 plus_scale_method, comma_scale_method) {
  mats <- list()
  abundance_log <- list()
  for (KOs in KO_vector) {
    has_comma <- base::grepl(",", KOs)
    has_plus <- base::grepl("\\+", KOs)

    if (has_plus && has_comma) {
      KO_subvector <- base::strsplit(KOs, ",")[[1]]
      loop_plus_result <- .loop_plus_core(KO_subvector, m, aggregrate_rowname,
                                          step_count, plus_scale_method)
      KO_scale <- paste(rownames(loop_plus_result[['abundance_matrix']]), collapse = ',')
      step_count <- loop_plus_result[['step_count']]
      abundance_log.tmp <- loop_plus_result[['abundance_log']]

      comma_result <- .step_comma_core(loop_plus_result[['abundance_matrix']], KO_scale,
                                       paste0(aggregrate_rowname, '_', step_count),
                                       step_count, comma_scale_method)
      result_matrix <- comma_result[['abundance_matrix']]
      step_count <- comma_result[['step_count']]
      abundance_log.tmp <- c(abundance_log.tmp, comma_result[['abundance_log']])
    } else if (has_comma) {
      comma_result <- .step_comma_core(m, KOs, paste0(aggregrate_rowname, '_', step_count),
                                       step_count, comma_scale_method)
      result_matrix <- comma_result[['abundance_matrix']]
      step_count <- comma_result[['step_count']]
      abundance_log.tmp <- comma_result[['abundance_log']]
    } else if (has_plus) {
      plus_result <- .step_plus_core(m, KOs, paste0(aggregrate_rowname, '_', step_count),
                                     step_count, plus_scale_method)
      result_matrix <- plus_result[['abundance_matrix']]
      step_count <- plus_result[['step_count']]
      abundance_log.tmp <- plus_result[['abundance_log']]
    } else {
      direct_result <- .step_direct_core(m, KOs)
      result_matrix <- direct_result[['abundance_matrix']]
      abundance_log.tmp <- direct_result[['abundance_log']]
    }

    mats[[length(mats) + 1]] <- result_matrix
    abundance_log <- c(abundance_log, abundance_log.tmp)
  }
  abundance_matrix <- .unique_rows(do.call(rbind, mats))
  list(abundance_matrix = abundance_matrix, step_count = step_count, abundance_log = abundance_log)
}

.module_step_core <- function(m, KO_string, aggregrate_rowname = 'bracket_1', step_count = 1,
                              plus_scale_method, comma_scale_method) {
  KO_string <- remove_outer_brackets(KO_string)

  has_comma <- base::grepl(",", KO_string)
  has_plus <- base::grepl("\\+", KO_string)
  has_space <- base::grepl(" ", KO_string)
  abundance_log.tmp <- list()

  if ((has_plus && has_comma && has_space)) {
    KO_vector <- base::strsplit(KO_string, " ")[[1]]
    loop_list <- .loop_plu_comma_core(KO_vector, m, aggregrate_rowname, step_count,
                                      plus_scale_method, comma_scale_method)
    abundance_log.tmp <- loop_list[['abundance_log']]
    KO_scale <- paste(rownames(loop_list[['abundance_matrix']]), collapse = '+')
    result_list <- .step_plus_core(loop_list[['abundance_matrix']], KO_scale,
                                   aggregrate_rowname, loop_list[['step_count']],
                                   plus_scale_method)
  } else if (has_plus && has_comma) {
    KO_vector <- base::strsplit(KO_string, ",")[[1]]
    loop_list <- .loop_plus_core(KO_vector, m, aggregrate_rowname, step_count,
                                 plus_scale_method)
    abundance_log.tmp <- loop_list[['abundance_log']]
    KO_scale <- paste(rownames(loop_list[['abundance_matrix']]), collapse = ',')
    result_list <- .step_comma_core(loop_list[['abundance_matrix']], KO_scale,
                                    aggregrate_rowname, loop_list[['step_count']],
                                    comma_scale_method)
  } else if (has_plus && has_space) {
    KO_vector <- base::strsplit(KO_string, " ")[[1]]
    loop_list <- .loop_plus_core(KO_vector, m, aggregrate_rowname, step_count,
                                 plus_scale_method)
    abundance_log.tmp <- loop_list[['abundance_log']]
    KO_scale <- paste(rownames(loop_list[['abundance_matrix']]), collapse = ' ')
    result_list <- .step_space_core(loop_list[['abundance_matrix']], KO_scale,
                                    aggregrate_rowname, loop_list[['step_count']])
  } else if (has_comma && has_space) {
    KO_vector <- base::strsplit(KO_string, " ")[[1]]
    loop_list <- .loop_comma_core(KO_vector, m, aggregrate_rowname, step_count,
                                  comma_scale_method)
    abundance_log.tmp <- loop_list[['abundance_log']]
    KO_scale <- paste(rownames(loop_list[['abundance_matrix']]), collapse = ' ')
    result_list <- .step_space_core(loop_list[['abundance_matrix']], KO_scale,
                                    aggregrate_rowname, loop_list[['step_count']])
  } else if (has_comma) {
    result_list <- .step_comma_core(m, KO_string, aggregrate_rowname, step_count,
                                    comma_scale_method)
  } else if (has_plus) {
    result_list <- .step_plus_core(m, KO_string, aggregrate_rowname, step_count,
                                   plus_scale_method)
  } else if (has_space) {
    result_list <- .step_space_core(m, KO_string, aggregrate_rowname, step_count)
  } else {
    result <- .step_direct_core(m, KO_string)
    result_list <- list(abundance_matrix = result[['abundance_matrix']],
                        step_count = step_count,
                        abundance_log = result[['abundance_log']])
  }
  list(abundance_matrix = result_list[['abundance_matrix']],
       step_count = result_list[['step_count']],
       abundance_log = c(abundance_log.tmp, result_list[['abundance_log']]))
}

.brackets_core <- function(m, module_steps_str, bracket_count = 1, step_count = 1,
                           module_name = 'Module', raw_module_steps = module_steps_str,
                           plus_scale_method, comma_scale_method, abundance_log = list()) {
  timestamp <- function() format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  log_entry <- function(msg) {
    paste0(timestamp(), " ", msg)
  }
  brackets <- extract_inner_brackets(module_steps_str)
  if (length(brackets) == 0) {
    abundance_log <- c(abundance_log, list(
      log_entry("Start processing final step..."),
      log_entry(paste0("Analyzing ", module_name, ": ", module_steps_str))
    ))
    module_step_list <- .module_step_core(m, module_steps_str,
                                          aggregrate_rowname = module_name,
                                          step_count = bracket_count,
                                          plus_scale_method,
                                          comma_scale_method)
    abundance_log <- c(abundance_log, module_step_list[['abundance_log']])
    return(list(
      data = module_step_list[['abundance_matrix']],
      log = abundance_log
    ))
  } else {
    abundance_log <- c(abundance_log, list(
      log_entry(paste0("Nested steps include: ")),
      log_entry(paste0('    ', brackets)),
      log_entry("Start processing nested steps...")
    ))
  }
  for (bracket in brackets) {
    bracket_name <- paste0(module_name, '_', bracket_count)
    abundance_log <- c(abundance_log, list(
      log_entry(paste0("Analyzing ", bracket_name, ": ", bracket)),
      log_entry(paste0("Bracket level: ", bracket_count))
    ))
    module_steps_str <- stringr::str_replace(module_steps_str, escape_special_chars(bracket), bracket_name)
    module_step_list <- .module_step_core(m, bracket, aggregrate_rowname = bracket_name,
                                          step_count = 1,
                                          plus_scale_method, comma_scale_method)
    m <- rbind(m, module_step_list[['abundance_matrix']])
    abundance_log <- c(abundance_log, module_step_list[['abundance_log']])
    bracket_count <- bracket_count + 1
  }
  m <- .unique_rows(m)
  .brackets_core(m = m,
                 module_steps_str = module_steps_str,
                 bracket_count = bracket_count,
                 step_count = step_count,
                 module_name = module_name,
                 raw_module_steps = raw_module_steps,
                 plus_scale_method = plus_scale_method,
                 comma_scale_method = comma_scale_method,
                 abundance_log = abundance_log)
}
