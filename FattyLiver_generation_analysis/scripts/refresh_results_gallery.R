# Enable UTF-8 before sourcing support files containing Chinese path labels.
if (.Platform$OS.type == "windows" && !isTRUE(l10n_info()[["UTF-8"]]))
  if (!nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8")))) stop("UTF-8 locale required")
# Collect existing public outputs only; never executes analysis scripts.
.fl_gallery_root <- local({
  override <- Sys.getenv("FATTYLIVER_ROOT", unset = "")
  paths <- unlist(lapply(sys.frames(), function(f) {
    p <- if (!is.null(f$ofile)) f$ofile else f$file
    if (is.character(p) && length(p) == 1L) p else character()
  }), use.names = FALSE)
  cli <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  anchors <- c(dirname(c(paths, cli)), getwd())
  candidates <- if (nzchar(override)) override else unique(c(anchors, dirname(anchors)))
  valid <- vapply(candidates, function(p) file.exists(file.path(p, "R", "result_gallery.R")), logical(1))
  if (!any(valid)) stop("Cannot locate project. Set FATTYLIVER_ROOT.")
  normalizePath(candidates[which(valid)[1L]], winslash = "/", mustWork = TRUE)
})
source(file.path(.fl_gallery_root, "R", "result_gallery.R"), local = TRUE, encoding = "UTF-8")
fl_refresh_results_main <- function(args = character(), project_root = .fl_gallery_root) {
  old <- Sys.getlocale("LC_CTYPE"); on.exit(Sys.setlocale("LC_CTYPE", old))
  if (.Platform$OS.type == "windows") suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8"))
  if (length(args) && !identical(args, "--self-test")) stop("Only --self-test or no argument is supported")
  if (identical(args, "--self-test")) {
    source(file.path(project_root, "tests", "result_gallery_test.R"), local = TRUE, encoding = "UTF-8")
    return(invisible(fl_gallery_self_test(project_root)))
  }
  fl_refresh_result_gallery(project_root)
}
if (sys.nframe() == 0L) fl_refresh_results_main(commandArgs(TRUE)) else
  if (isTRUE(getOption("fatty_liver.autorun", TRUE))) fl_refresh_results_main()
