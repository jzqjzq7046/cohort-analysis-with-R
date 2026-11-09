# Viewing copies only: no clinical data or analysis functions are executed.
fl_gallery_url <- function(path) {
  vapply(strsplit(gsub("\\", "/", path, fixed = TRUE), "/", fixed = TRUE), function(parts)
    paste(vapply(parts, utils::URLencode, character(1), reserved = TRUE), collapse = "/"), character(1))
}

fl_gallery_escape <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE); x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub('"', "&quot;", x, fixed = TRUE)
}

fl_gallery_category <- function(path) {
  ext <- tolower(tools::file_ext(path)); name <- basename(path)
  if (ext %in% c("csv", "tsv", "xlsx", "xls", "xlsm", "ods")) return("Tables_表格")
  if (ext %in% c("html", "htm", "pdf") &&
      (grepl("^tables/", path) || grepl("^table[0-9_ .-]", name, ignore.case = TRUE))) return("Tables_表格")
  if (ext %in% c("png", "jpg", "jpeg", "tif", "tiff", "svg", "pdf", "eps", "webp", "gif", "bmp"))
    return("Figures_图片")
  if (ext %in% c("html", "htm", "md", "txt", "docx", "rtf", "pptx")) return("Reports_报告")
  NA_character_
}

# This guard is applied to the resolved absolute target before recursive cleanup.
fl_gallery_remove_owned <- function(path, project_root, prefix) {
  resolved <- normalizePath(path, winslash = "/", mustWork = TRUE)
  root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  if (!identical(tolower(dirname(resolved)), tolower(root)) ||
      !startsWith(basename(resolved), prefix)) stop("Unsafe gallery cleanup target: ", resolved)
  if (unlink(resolved, recursive = TRUE) != 0L) warning("Could not remove viewing-copy temporary folder: ", resolved)
  invisible(NULL)
}

fl_refresh_result_gallery_legacy <- function(project_root, gallery_name = "Results_表格与图片") {
  root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  if (length(gallery_name) != 1L || !nzchar(gallery_name) ||
      grepl("[/\\\\:]", gallery_name) || startsWith(gallery_name, ".")) stop("Use a simple gallery folder name")
  target <- file.path(root, gallery_name)
  marker <- ".managed_results_gallery"
  if (file.exists(target) && (!dir.exists(target) || !file.exists(file.path(target, marker))))
    stop("Refusing to replace an unmanaged folder: ", target)
  if (dir.exists(target) && !identical(tolower(dirname(normalizePath(target, winslash = "/", mustWork = TRUE))), tolower(root)))
    stop("Gallery resolves outside its project")
  output_root <- file.path(root, "outputs")
  modules <- list.dirs(output_root, recursive = FALSE, full.names = TRUE)
  modules <- sort(modules[grepl("^[0-9]+_", basename(modules))])
  files <- runs <- list(); k <- 0L
  categories <- c("Tables_表格", "Figures_图片", "Reports_报告")
  empty <- data.frame(module = character(), run_id = character(), category = character(),
    file = character(), viewing_path = character(), source_path = character(),
    source_md5 = character(), copied_md5 = character(), bytes = numeric(), stringsAsFactors = FALSE)
  for (module in modules) {
    adopted <- list.dirs(file.path(module, "current"), recursive = FALSE, full.names = TRUE)
    if (!length(adopted)) next
    if (length(adopted) != 1L) stop("More than one adopted run: ", module)
    run <- adopted[1L]; status <- file.path(run, "run_status.txt")
    if (!file.exists(status) || !identical(readLines(status, n = 1L, warn = FALSE), "COMPLETED"))
      stop("Adopted run is not complete: ", run)
    runs[[basename(module)]] <- data.frame(module = basename(module), run_id = basename(run))
    relative <- list.files(run, recursive = TRUE, full.names = FALSE)
    relative <- gsub("\\", "/", relative, fixed = TRUE)
    relative <- relative[!grepl("(^|/)(private|audit|history|_running)(/|$)", relative, ignore.case = TRUE) &
      !basename(relative) %in% c("run_status.txt", "STATUS", "CURRENT", "CURRENT.txt")]
    for (p in relative) {
      category <- fl_gallery_category(p)
      if (is.na(category)) next
      source <- normalizePath(file.path(run, p), winslash = "/", mustWork = TRUE)
      run_root <- paste0(normalizePath(run, winslash = "/", mustWork = TRUE), "/")
      if (!startsWith(tolower(source), tolower(run_root))) stop("Result resolves outside adopted run: ", source)
      # Strip only known result subdirectories, preserving nested names elsewhere.
      display <- sub("^(tables|figures|exploration|qc|preprocessing|clinical|similarity)/", "", p)
      k <- k + 1L
      files[[k]] <- data.frame(module = basename(module), run_id = basename(run), category = category,
        file = basename(p), viewing_path = paste(category, basename(module), display, sep = "/"),
        source_path = source, source_md5 = unname(tools::md5sum(source)), copied_md5 = "",
        bytes = unname(file.info(source)$size), stringsAsFactors = FALSE)
    }
  }
  manifest <- if (length(files)) do.call(rbind, files) else empty
  if (anyDuplicated(tolower(manifest$viewing_path))) stop("Colliding result names in viewing copies")
  stamp <- paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())
  stage_prefix <- ".result_gallery_build_"; backup_prefix <- ".result_gallery_previous_"
  stage <- file.path(root, paste0(stage_prefix, stamp)); backup <- file.path(root, paste0(backup_prefix, stamp))
  if (!dir.create(stage)) stop("Cannot create viewing-copy staging folder")
  committed <- FALSE
  on.exit(if (!committed && dir.exists(stage)) fl_gallery_remove_owned(stage, root, stage_prefix), add = TRUE)
  for (category in categories) dir.create(file.path(stage, category))
  for (i in seq_len(nrow(manifest))) {
    destination <- file.path(stage, manifest$viewing_path[i])
    dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(manifest$source_path[i], destination, copy.date = TRUE)) stop("Result copy failed")
    manifest$copied_md5[i] <- unname(tools::md5sum(destination))
  }
  if (!identical(manifest$source_md5, manifest$copied_md5)) stop("Viewing-copy checksum mismatch")
  # Keep user-added files, even though only generated copies appear in the index.
  if (dir.exists(target)) {
    old_manifest <- file.path(target, "file_manifest.csv")
    if (!file.exists(old_manifest)) stop("Managed gallery manifest missing")
    old <- utils::read.csv(old_manifest, fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
    owned <- c(old$viewing_path, "index.html", "README_zh.md", "README_en.md", "file_manifest.csv", marker)
    extra <- setdiff(gsub("\\", "/", list.files(target, recursive = TRUE, all.files = TRUE), fixed = TRUE), owned)
    # Office owner/lock files are transient process state, not user documents.
    # Leave them in the open target; never copy or remove the active locks.
    extra <- extra[!startsWith(basename(extra), "~$")]
    for (p in extra) {
      src <- normalizePath(file.path(target, p), winslash = "/", mustWork = TRUE)
      if (!startsWith(tolower(src), paste0(tolower(normalizePath(target, winslash = "/", mustWork = TRUE)), "/")))
        stop("User file resolves outside gallery")
      dest <- if (p %in% manifest$viewing_path) file.path(stage, "User_files_用户文件", stamp, p) else file.path(stage, p)
      dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
      if (!file.copy(src, dest, copy.date = TRUE)) stop("Cannot preserve user-added viewing-folder file: ", p)
    }
  }
  con <- file(file.path(stage, "file_manifest.csv"), "w", encoding = "UTF-8")
  writeLines("\ufeff", con, sep = ""); utils::write.table(manifest, con, sep = ",", row.names = FALSE, na = "")
  close(con)
  writeLines("managed viewing copies; authoritative analyses remain in outputs", file.path(stage, marker))
  write_doc <- function(name, text) {
    con <- file(file.path(stage, name), "w", encoding = "UTF-8"); writeLines(text, con); close(con)
  }
  write_doc("README_zh.md", c("# 表格、图片与报告统一查看目录", "",
    "双击index.html打开总目录，可搜索并点击表格、图片和报告。", "",
    "- Tables_表格：CSV/TSV/Excel及Table 1等可浏览表格，按分析模块归类。",
    "- Figures_图片：PNG/PDF/SVG等图片及可编辑版本，按分析模块归类。",
    "- Reports_报告：各模块HTML、Markdown及其他说明。",
    "- file_manifest.csv：来源运行、原文件路径和复制前后MD5，方便溯源。", "",
    "只汇集各模块current中的成功结果，未混入旧版或失败运行；历史结果仍在outputs各模块history中。",
    "不汇集private、audit、逐人数据包和原始数据；这些仍保留在原分析位置。",
    "这些是查看副本；分析成功后自动刷新。不要在自动生成的副本上保存需要长期保留的修改。用户新增文件会保留。",
    "报告中的相对链接按原始目录编写；如复制的报告内链接不能打开，请点击总目录中的“原始位置”。", "",
    "独立刷新：在RStudio Source scripts/refresh_results_gallery.R；只收集已有结果，不运行统计分析。",
    paste0("本次更新：", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")), "", "[English](README_en.md)"))
  write_doc("README_en.md", c("# Central result viewer", "", "Open index.html to browse/search adopted tables, figures and reports.",
    "Viewing copies are grouped by category and analysis module. Successful analysis publication refreshes this folder automatically.",
    "Only COMPLETED current runs are collected; historical, failed, private and audit outputs remain in their original locations.",
    "file_manifest.csv records source run IDs and matching source/copy MD5 hashes. Generated copies are refreshed; keep permanent edits elsewhere.",
    "User-added files are retained. For links inside copied reports, use the Original link in the index when a relative link no longer resolves.",
    "Source scripts/refresh_results_gallery.R to refresh copies without rerunning any statistical analysis.", "", "[中文](README_zh.md)"))
  html <- c("<!doctype html><html lang='zh-CN'><head><meta charset='utf-8'><title>表格与图片 · Results</title>",
    "<style>body{font:15px Arial,'Microsoft YaHei',sans-serif;max-width:1400px;margin:30px auto;padding:0 24px;color:#243746;line-height:1.6}a{color:#19618b}table{width:100%;border-collapse:collapse}th,td{border-bottom:1px solid #dfe5e8;padding:8px;text-align:left;word-break:break-word}th{background:#f0f4f6}input{padding:10px;width:95%;border:1px solid #aebdc6;border-radius:5px}small{color:#647784}nav{display:flex;gap:18px;flex-wrap:wrap}.tag{color:#516672;font-size:12px}details{margin:14px 0}img{max-width:95%;max-height:480px;object-fit:contain}.scroll{overflow:auto}h2{margin-top:30px}</style></head><body>",
    "<h1>表格与图片 · Results</h1>",
    paste0("<p>已汇集 ", length(runs), " 个分析模块，", nrow(manifest), " 个当前结果文件。<small>更新：", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "</small></p>"),
    "<p>查看副本按模块归类；每次成功运行自动更新。原始结果、历史记录和审计资料继续保留。</p>",
    "<nav><a href='统一配色_新版图/index.html'>统一配色新版图</a><a href='Results_重点表格/index.html'>重点表格 / Key tables</a><a href='#tables'>表格</a><a href='#figures'>图片</a><a href='#reports'>报告</a><a href='file_manifest.csv'>来源清单</a><a href='README_zh.md'>说明</a><a href='../outputs/RESULTS_zh.md'>研究结果总结</a></nav>",
    "<p><input id='q' type='search' placeholder='搜索模块、文件名，例如 Table1、Figure1、models、07…' oninput='filterRows()'></p>")
  for (j in seq_along(categories)) {
    category <- categories[j]; rows <- manifest[manifest$category == category, , drop = FALSE]
    html <- c(html, paste0("<section><h2 id='", c("tables", "figures", "reports")[j], "'>", category,
      " <small>(", nrow(rows), ")</small></h2><div class='scroll'><table><thead><tr><th>模块 / Module</th><th>文件 / File</th><th>来源 / Source</th></tr></thead><tbody>"))
    for (i in seq_len(nrow(rows))) {
      src_relative <- substring(rows$source_path[i], nchar(root) + 2L)
      html <- c(html, paste0("<tr class='result'><td>", fl_gallery_escape(rows$module[i]), "</td><td><a href='",
        fl_gallery_url(rows$viewing_path[i]), "'>", fl_gallery_escape(rows$file[i]), "</a></td><td><a href='",
        fl_gallery_url(paste0("../", src_relative)), "'>原始位置 / Original</a><br><small>",
        fl_gallery_escape(rows$run_id[i]), "</small></td></tr>"))
    }
    html <- c(html, "</tbody></table></div></section>")
  }
  html <- c(html, "<script>function filterRows(){const q=document.getElementById('q').value.toLowerCase();document.querySelectorAll('tr.result').forEach(r=>r.hidden=!r.textContent.toLowerCase().includes(q));}</script></body></html>")
  write_doc("index.html", html)
  renamed <- !dir.exists(target) || suppressWarnings(file.rename(target, backup))
  if (!renamed) {
    # Windows Explorer/current-directory handles may lock the folder name while
    # permitting file updates. For an additive refresh, retain the open folder,
    # copy only verified staged content, then publish the manifest/index last.
    # No old managed or user files are deleted in this fallback. If retirement
    # is needed, retain the old viewer and require the normal directory swap.
    obsolete <- setdiff(old$viewing_path, manifest$viewing_path)
    if (length(obsolete)) stop("Viewing folder is open and old managed files need retirement; close it and refresh")
    rel <- gsub("\\", "/", list.files(stage, recursive = TRUE, all.files = TRUE), fixed = TRUE)
    controls <- c("README_zh.md", "README_en.md", marker, "file_manifest.csv", "index.html")
    order <- c(setdiff(rel, controls), intersect(controls, rel))
    target_resolved <- paste0(tolower(normalizePath(target, winslash = "/", mustWork = TRUE)), "/")
    for (p in order) {
      src <- normalizePath(file.path(stage, p), winslash = "/", mustWork = TRUE)
      dest <- file.path(target, p)
      dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
      parent <- paste0(tolower(normalizePath(dirname(dest), winslash = "/", mustWork = TRUE)), "/")
      if (!startsWith(parent, target_resolved)) stop("Fallback destination resolves outside viewer")
      if (file.exists(dest)) {
        dest_resolved <- tolower(normalizePath(dest, winslash = "/", mustWork = TRUE))
        if (!startsWith(dest_resolved, target_resolved)) stop("Fallback file resolves outside viewer")
        if (identical(unname(tools::md5sum(src)), unname(tools::md5sum(dest)))) next
      }
      if (!file.copy(src, dest, overwrite = TRUE, copy.date = TRUE) ||
          !identical(unname(tools::md5sum(src)), unname(tools::md5sum(dest))))
        stop("Cannot update open viewing folder: ", p)
    }
    actual <- unname(tools::md5sum(file.path(target, manifest$viewing_path)))
    if (!identical(actual, manifest$source_md5)) stop("Open-folder copy checksum mismatch")
    fl_gallery_remove_owned(stage, root, stage_prefix)
    cat("Open Windows viewing folder refreshed in place; no old files removed.\n")
  } else if (!file.rename(stage, target)) {
    if (dir.exists(backup) && !file.rename(backup, target)) warning("Previous viewing copies remain recoverable at ", backup)
    stop("Cannot publish refreshed viewing copies")
  }
  committed <- TRUE
  if (dir.exists(backup)) fl_gallery_remove_owned(backup, root, backup_prefix)
  cat("Result viewer updated: ", file.path(target, "index.html"), "\n",
    nrow(manifest), " files from ", length(runs), " current modules; copy MD5 checks passed.\n", sep = "")
  invisible(manifest)
}

# The study-ordered viewer replaces simultaneous category/special-topic copies.
# Keep the original implementation for synthetic fixtures and projects without
# the explicit version-2 convention. Neither implementation runs statistics.
fl_refresh_result_gallery <- function(project_root, gallery_name = "Results_表格与图片") {
  settings <- file.path(project_root, "config", "results_viewer.json")
  if (!file.exists(settings)) return(fl_refresh_result_gallery_legacy(project_root, gallery_name))
  if (!identical(gallery_name, "Results_表格与图片")) stop("The organized viewer has one fixed result directory")
  python <- "C:/Users/jzqjzq/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe"
  if (!file.exists(python)) python <- Sys.which("python")
  if (!nzchar(python)) stop("Python runtime for viewing-copy organization is unavailable")
  status <- system2(python, c(shQuote(file.path(project_root, "tools", "organize_results.py")),
                             "--root", shQuote(normalizePath(project_root, winslash = "/", mustWork = TRUE))))
  if (status != 0L) stop("Result organization failed; previous viewer is preserved")
  manifest <- utils::read.csv(file.path(project_root, gallery_name, "file_manifest.csv"),
                             fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE)
  invisible(manifest)
}

fl_publish_organized_results <- function(cfg) {
  if (!file.exists(file.path(cfg$project_root, "config", "results_viewer.json"))) return(NULL)
  fl_refresh_result_gallery(cfg$project_root)
  invisible(file.path(cfg$project_root, "Results_表格与图片"))
}
