# R/zzz.R

.onLoad <- function(libname, pkgname) {
  # 1. 初始化包选项（保持不变）
  op <- options()
  op.mclink <- list(
    mclink.path = tempdir(),
    mclink.verbose = FALSE
  )
  toset <- !(names(op.mclink) %in% names(op))
  if (any(toset)) options(op.mclink[toset])

  # 2. 安全加载数据（关键改进）
  ns <- asNamespace(pkgname)
  tryCatch({
    utils::data(list = c("KO_pathway_ref", "KO_Sample_wide"),
                package = pkgname,
                envir = ns)
  }, error = function(e) {
    warning("Data load failed: ", e$message, call. = FALSE)
  })

  invisible()
}

.onAttach <- function(libname, pkgname) {
  # 使用完全限定名调用（关键修复）
  version <- utils::packageVersion(pkgname)
  packageStartupMessage("Welcome to mclink v", version)
}
