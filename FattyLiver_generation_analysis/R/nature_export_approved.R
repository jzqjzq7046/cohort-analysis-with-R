sys.source(file.path(Sys.getenv('FATTYLIVER_ROOT',unset='E:/YFS/workspace/FattyLiver_generation_analysis'),'R/model_display_labels.R'),envir=.GlobalEnv)
# R-only reusable export/layout helpers; structural adaptation of existing project figures.
nf_skill_root <- function() {
  p <- Sys.getenv("NATURE_FIGURE_SKILL", unset = file.path(Sys.getenv("USERPROFILE"), ".codex/skills/nature-figure"))
  if (!file.exists(file.path(p, "scripts/panel_alignment.R"))) stop("Nature figure QA helpers missing; set NATURE_FIGURE_SKILL")
  normalizePath(p, winslash = "/")
}

nf_theme <- function() ggplot2::theme_classic(base_size = 6.3, base_family = "Arial") +
  ggplot2::theme(
    axis.text = ggplot2::element_text(size = 5.7, colour = "black"),
    axis.title = ggplot2::element_text(size = 6.2),
    axis.line = ggplot2::element_line(linewidth = .25, colour = "black"),
    axis.ticks = ggplot2::element_line(linewidth = .25),
    axis.ticks.length = grid::unit(1, "mm"),
    panel.grid = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    plot.background = ggplot2::element_blank(),
    legend.position = "bottom", legend.title = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 5.7, colour = "black"),
    legend.key.size = grid::unit(3.0, "mm"),
    legend.margin = ggplot2::margin(0, 0, 0, 0, "pt"),
    legend.box.spacing = grid::unit(.6, "mm"),
    plot.title = ggplot2::element_text(size = 6.8, face = "plain", margin = ggplot2::margin(b = 3, l = 9, unit = "pt")),
    plot.tag = ggplot2::element_text(size = 8, face = "bold", colour = "black"),
    plot.tag.position = c(0, 1),
    plot.margin = ggplot2::margin(4, 4, 4, 4, "pt"))

nf_patch <- function(plots, design, heights = NULL, widths = NULL, collect_guides = FALSE) {
  # Apply themes before assembly: avoids the incompatible '&' operator in the
  # already-installed patchwork 1.3 / ggplot2 4 combination.
  plots <- lapply(seq_along(plots), function(i) plots[[i]] + nf_theme() + ggplot2::labs(tag = letters[i]))
  p <- patchwork::wrap_plots(plots, design = design, heights = heights, widths = widths,
    guides = if (collect_guides) "collect" else "keep")
  p
}


nf_alignment <- function(plot, stem, width_mm, height_mm, audit_dir, cfg) {
  plot <- fl_model_plot_labels(plot)
  source(file.path(nf_skill_root(), "scripts/panel_alignment.R"), local = TRUE)
  probe <- file.path(audit_dir, paste0(stem, ".layout-probe.pdf"))
  grDevices::cairo_pdf(probe, width = width_mm / 25.4, height = height_mm / 25.4, family = "Arial")
  open <- TRUE; on.exit(if (open) grDevices::dev.off())
  grob <- patchwork::patchworkGrob(plot)
  layout_restore <- getOption("fatty_liver.saved_layout_restore")
  if (is.function(layout_restore)) grob <- layout_restore(grob, stem)
  grob <- fl_model_label_grob(grob)
  grid::grid.newpage(); grid::grid.draw(grob); grid::grid.force()
  rows <- .nature_alignment_panel_rows(grob)
  # When a two-column row sits above a three-column row, patchwork allocates
  # some fixed axis space to the middle lower panel. Balance only flexible
  # layout columns, preserving all data, text, and fixed-size margins.
  row_key <- paste(rows$t, rows$b)
  row_sets <- split(seq_len(nrow(rows)), row_key)
  needs_balance <- any(vapply(row_sets, function(i) length(i) >= 3L && length(unique(rows$r[i] - rows$l[i])) == 1L, logical(1)))
  if (needs_balance) {
    nc <- length(grob$widths); nr <- length(grob$heights)
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nr, nc, widths = grob$widths, heights = grob$heights)))
    physical <- vapply(seq_len(nc), function(j) {
      grid::pushViewport(grid::viewport(layout.pos.col = j))
      w <- grid::convertWidth(grid::unit(1, "npc"), "inches", valueOnly = TRUE) * 72
      grid::popViewport(); w
    }, numeric(1))
    grid::popViewport()
    flexible <- which(grid::unitType(grob$widths) == "null")
    equations <- list(rep(1, nc))
    for (ii in row_sets) {
      if (length(ii) < 2L || length(unique(rows$r[ii] - rows$l[ii])) != 1L) next
      base <- rep(0, nc); base[rows$l[ii[1]]:rows$r[ii[1]]] <- 1
      for (j in ii[-1]) {
        eq <- base; eq[rows$l[j]:rows$r[j]] <- eq[rows$l[j]:rows$r[j]] - 1
        equations[[length(equations) + 1L]] <- eq
      }
    }
    A <- do.call(rbind, equations); target <- -drop(A %*% physical); target[1] <- 0
    M <- A[, flexible, drop = FALSE]; dec <- svd(M)
    keep <- dec$d > 1e-9
    delta <- drop(dec$v[, keep, drop = FALSE] %*% ((crossprod(dec$u[, keep, drop = FALSE], target)) / dec$d[keep]))
    if (max(abs(M %*% delta - target)) > 1e-6) stop("Cannot equalize the final comparable panel widths")
    physical[flexible] <- physical[flexible] + delta
    if (any(physical < -1e-8)) stop("Panel balancing would create a negative layout width")
    grob$widths <- grid::unit(pmax(physical, 0), "bigpts")
    grid::grid.newpage(); grid::grid.draw(grob); grid::grid.force()
  }
  # A gtable's flexible ('null') units do not resolve from its bare unit vectors.
  # Measure the drawn panel viewports, after the layout engine has resolved them.
  vp_names <- grid::grid.ls(viewports = TRUE, grobs = FALSE, print = FALSE)$name
  measure <- function(name, anchor = FALSE) {
    match <- vp_names[startsWith(vp_names, paste0(name, "."))]
    if (length(match) != 1L) stop("Cannot uniquely locate final viewport: ", name)
    grid::seekViewport(match)
    on.exit(grid::upViewport(0))
    low <- grid::deviceLoc(grid::unit(0, "npc"), grid::unit(0, "npc"), valueOnly = TRUE)
    high <- grid::deviceLoc(grid::unit(1, "npc"), grid::unit(1, "npc"), valueOnly = TRUE)
    if (anchor) return(72 * c(low$x, high$y))
    72 * c(low$x, low$y, high$x, high$y)
  }
  panels <- lapply(seq_len(nrow(rows)), function(i) {
    r <- rows[i, ]; number <- as.integer(sub("^panel-", "", r$name))
    panel <- list(id = letters[number], bbox_pt = measure(r$name),
      grid_id = "final-patchwork", row_start = r$t - 1L, row_stop = r$b, col_start = r$l - 1L, col_stop = r$r)
    tag <- grob$layout[grob$layout$name == paste0("tag-", number), , drop = FALSE]
    if (nrow(tag) == 1L) {
      panel$panel_label <- letters[number]
      panel$panel_label_anchor_pt <- measure(paste0("tag-", number), anchor = TRUE)
    }
    panel
  })
  layout <- list(schema_version = 1L, backend = "r-patchwork-cairo", measurement = "Final Cairo device; rendered gtable panel cells in PDF points",
    figure = list(width_pt = width_mm / 25.4 * 72, height_pt = height_mm / 25.4 * 72), panels = panels, exemptions = list())
  manifest <- file.path(audit_dir, paste0(stem, ".alignment-layout.json"))
  jsonlite::write_json(layout, manifest, auto_unbox = TRUE, pretty = TRUE, digits = NA)
  grDevices::dev.off(); open <- FALSE
  args <- c(shQuote(file.path(nf_skill_root(), "scripts/audit_panel_alignment.py")), shQuote(manifest),
    "--json-out", shQuote(file.path(audit_dir, paste0(stem, ".alignment.json"))), "--strict")
  status <- system2(cfg$python, args)
  if (status != 0L) stop("Rendered panel alignment gate blocked export: ", stem)
  invisible(list(layout = layout, grob = grob))
}

nf_qa_command <- function(cfg, script, args, log) {
  output <- suppressWarnings(system2(cfg$python, c(shQuote(file.path(nf_skill_root(), "scripts", script)), args), stdout = TRUE, stderr = TRUE))
  writeLines(output, log, useBytes = TRUE)
  status <- attr(output, "status"); if (is.null(status)) status <- 0L
  if (status != 0L) stop("Figure QA blocked delivery: ", script, "; see ", log)
  invisible(output)
}
nf_export <- function(spec, out, cfg) {
  audit <- file.path(out, "audit"); stem <- spec$name
  rendered <- nf_alignment(spec$plot, stem, spec$width, spec$height, audit, cfg)
  pdf <- file.path(out, "figures", paste0(stem, ".pdf"))
  grDevices::cairo_pdf(pdf, width = spec$width / 25.4, height = spec$height / 25.4, family = "Arial", bg = "white")
  tryCatch({ grid::grid.newpage(); grid::grid.draw(rendered$grob) }, finally = grDevices::dev.off())
  # Preserve unmodified lightweight reports. Their whole-trace/Tf assumptions
  # are checked against effective Cairo text geometry by the blocking adapter.
  raw_text <- suppressWarnings(system2(cfg$python, c(shQuote(file.path(nf_skill_root(), "scripts/audit_pdf_text.py")), shQuote(pdf), "--min-pt", "5", "--json"), stdout = TRUE, stderr = TRUE))
  writeLines(raw_text, file.path(audit, paste0(stem, ".raw-Tf-audit.json")))
  raw_collision <- suppressWarnings(system2(cfg$python, c(shQuote(file.path(nf_skill_root(), "scripts/audit_figure_collisions.py")), shQuote(pdf),
    "--json-out", shQuote(file.path(audit, paste0(stem, ".raw-trace-collision.json")))), stdout = TRUE, stderr = TRUE))
  writeLines(raw_collision, file.path(audit, paste0(stem, ".raw-collision-log.txt")))
  adapter <- file.path(cfg$project_root, "scripts/audit_nature_cairo_pdf.py")
  status <- system2(cfg$python, c(shQuote(adapter), shQuote(pdf), "--skill", shQuote(nf_skill_root()),
    "--prefix", shQuote(file.path(audit, stem))))
  if (status != 0L) stop("Effective rendered text/collision gate blocked export: ", stem)
  svglite::svglite(file.path(out, "figures", paste0(stem, ".svg")), width = spec$width / 25.4, height = spec$height / 25.4, bg = "white")
  tryCatch({ grid::grid.newpage(); grid::grid.draw(rendered$grob) }, finally = grDevices::dev.off())
  grDevices::png(file.path(out, "figures", paste0(stem, ".png")), width = spec$width, height = spec$height, units = "mm", res = 600, type = "cairo", bg = "white")
  tryCatch({ grid::grid.newpage(); grid::grid.draw(rendered$grob) }, finally = grDevices::dev.off())
  if (isTRUE(spec$export_tiff)) {
    grDevices::tiff(file.path(out, "figures", paste0(stem, ".tiff")), width = spec$width, height = spec$height,
      units = "mm", res = 600, type = "cairo", compression = "lzw", bg = "white")
    tryCatch({ grid::grid.newpage(); grid::grid.draw(rendered$grob) }, finally = grDevices::dev.off())
  }
  data.frame(figure = stem, nominal_width_mm = spec$width, nominal_height_mm = spec$height,
    body_text_pt = "5.5-6.8", panel_label_pt = 8, preview_dpi = 600,
    editable_main_file = paste0(stem, ".pdf"), formats = if (isTRUE(spec$export_tiff)) "PDF;SVG;PNG;TIFF" else "PDF;SVG;PNG", panel_count = spec$panels)
}


nf_baseline <- function(d) {
  nflow <- function(s) d$flow$n[match(s, d$flow$stage)]
  f <- d$composition$families[d$composition$scenario == "core_adjusted"][1]
  core <- d$samples[d$samples$scenario == "core_adjusted", ]
  adult <- d$samples[d$samples$scenario == "adult_BMI", ]
  boxes <- data.frame(x = rep(1, 5), y = 5:1,
    label = c(paste0("Clinical records: ", nflow("Clinical source records"), "\nMicrobiome samples: ", nflow("Microbiome participants")),
      paste0("Clinical-microbiome matched\nn = ", nflow("Matched participants")),
      paste0("Ultrasound assessed\nn = ", nflow("Final ultrasound cohort")),
      paste0("Three generations represented\n", f, " primary families; n = ", nflow("Three-generation family cohort")),
      paste0("Core: n = ", sum(core$n), "\nAdult common sample: n = ", sum(adult$n))))
  a <- ggplot2::ggplot(boxes, ggplot2::aes(x, y)) +
    ggplot2::geom_segment(data = data.frame(x = 1, y = 4.70:1.70),
      ggplot2::aes(x = 1, xend = 1, y = y, yend = y - .40), inherit.aes = FALSE,
      arrow = grid::arrow(length = grid::unit(1.5, "mm"), type = "closed"), linewidth = .3) +
    ggplot2::geom_rect(ggplot2::aes(xmin = .12, xmax = 1.88, ymin = y - .28, ymax = y + .28), fill = "white", colour = "#707070", linewidth = .3) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 2.15, lineheight = 1.08, family = "Arial") +
    ggplot2::coord_cartesian(xlim = c(0, 2), ylim = c(.5, 5.5), expand = FALSE) +
    ggplot2::labs(title = "Cohort") + ggplot2::theme_void(base_family = "Arial")
  age <- d$age; age$generation <- factor(age$generation, levels = c("G2", "G1", "G0"))
  b <- ggplot2::ggplot(age, ggplot2::aes(median, generation, fill = generation)) +
    ggplot2::geom_vline(xintercept = 18, linetype = 2, colour = "#777777", linewidth = .25) +
    ggplot2::geom_segment(ggplot2::aes(x = minimum, xend = maximum, yend = generation), linewidth = .18, colour = "#888888") +
    ggplot2::geom_segment(ggplot2::aes(x = q1, xend = q3, yend = generation, colour = generation), linewidth = .35) +
    ggplot2::geom_point(shape = 21, size = 1.8, stroke = .3) +
    ggplot2::scale_colour_manual(values = fl_viz_generation) + nf_palette_scale() +
    ggplot2::scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = ggplot2::expansion(mult = .02)) +
    ggplot2::labs(x = "Age (years)", y = NULL, title = "Age: median and IQR") + ggplot2::guides(fill = "none", colour = "none")
  p <- d$proportions; p <- p[p$group %in% c("G0", "G1", "G2"), ]; p$group <- factor(p$group, levels = c("G0", "G1", "G2"))
  c <- ggplot2::ggplot(p, ggplot2::aes(group, 100 * proportion, fill = group)) +
    ggplot2::geom_col(width = .55, colour = "#333333", linewidth = .25) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(cases, "/", n)), nudge_y = 3, size = 1.95, family = "Arial") +
    ggplot2::scale_fill_manual(values = fl_viz_generation) +
    ggplot2::scale_y_continuous(limits = c(0, 38), breaks = seq(0, 30, 10), expand = ggplot2::expansion(mult = c(0, .02))) +
    ggplot2::labs(x = "Generation", y = "Steatosis (%)", title = "Ultrasound steatosis") + ggplot2::guides(fill = "none")
  p2 <- d$proportions[match(c("G2 <18", "G2 >=18"), d$proportions$group), ]; p2$age <- factor(c("<18", ">=18"), levels = c("<18", ">=18"))
  e <- ggplot2::ggplot(p2, ggplot2::aes(age, 100 * proportion, fill = age)) +
    ggplot2::geom_col(width = .50, colour = "#333333", linewidth = .25) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(cases, "/", n)), nudge_y = .6, size = 1.95, family = "Arial") +
    ggplot2::scale_fill_manual(values = c("<18" = "#79CAFB", ">=18" = "#F7A24F")) +
    ggplot2::scale_y_continuous(limits = c(0, 8), breaks = c(0, 2, 4, 6, 8), expand = ggplot2::expansion(mult = c(0, .02))) +
    ggplot2::labs(x = "G2 age (years)", y = "Steatosis (%)", title = "G2 age groups") + ggplot2::guides(fill = "none")
  p <- nf_patch(list(a, b, c, e), "AB\nAC\nAD", widths = c(1, 1.15))
  # A is a flow diagram, so axes would misrepresent its non-quantitative coordinates.
  p[[1]] <- p[[1]] + ggplot2::theme(axis.line = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
    axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank())
  p
}

