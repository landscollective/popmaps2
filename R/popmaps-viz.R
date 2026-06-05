#' Plot modern POPMAPS raster outputs
#'
#' @description
#' `plot_popmaps()` draws cleaner, `terra`-based maps from the list returned by
#' [popmaps()]. It is intended as the modern plotting layer for popmaps2, while
#' [popmap_viz()] remains available for POPMAPS 1.03 compatibility.
#'
#' @param pop_raster_list A list returned by [popmaps()].
#' @param input_raster The raster used to create `pop_raster_list`. It may be a
#'   `terra::SpatRaster`, a legacy `raster::RasterLayer`, or a raster file path.
#' @param input_locs Optional POPMAPS location table. When supplied, empirical
#'   sampling locations are drawn as ancestry pies or dominant-ancestry points.
#' @param type Map type to draw. `"confidence"` plots the dominant ancestry
#'   confidence layer, `"boundary"` plots hard population assignments, and
#'   `"axis"` plots one weighted ancestry-axis layer. `"ancestry"` is accepted
#'   as a legacy alias for `"confidence"`.
#' @param style Plot style. `"modern"` uses the default popmaps2 map style.
#'   `"manuscript"` uses a Massatti and Winkler (2022)-inspired style with
#'   grayscale confidence, colored hard-boundary outlines, muted suitability
#'   background, and state outlines when available.
#' @param axis Ancestry axis used when `type = "axis"`. Supply either an integer
#'   axis number or a layer name such as `"axis_1"`.
#' @param sites How to draw empirical locations when `input_locs` is supplied.
#'   `"pies"` draws ancestry-proportion pies, `"points"` draws points colored by
#'   dominant ancestry, and `"none"` suppresses site symbols.
#' @param site_legend Logical. If `TRUE`, draw an ancestry-axis legend for site
#'   symbols.
#' @param boundaries Logical. If `TRUE`, hard-boundary outlines are overlaid on
#'   ancestry and axis maps.
#' @param palette Continuous color palette for ancestry and axis maps. Named
#'   palettes are matched to base R HCL palettes; common choices include
#'   `"Viridis"`, `"Mako"`, `"Inferno"`, `"YlGnBu"`, and `"Purple-Yellow"`.
#'   A custom color vector may also be supplied.
#' @param boundary_palette Categorical color palette for ancestry axes and hard
#'   boundaries. Named palettes are matched to base R HCL palettes, or a custom
#'   color vector may be supplied.
#' @param background_raster Optional suitability, habitat, or prediction raster
#'   drawn beneath the POPMAPS confidence layer. In `style = "manuscript"`, this
#'   defaults to `input_raster`.
#' @param background_threshold Optional numeric cutoff. When supplied, ancestry
#'   confidence or axis values are shown only where `background_raster` is
#'   greater than or equal to the cutoff; lower cells show only the background.
#' @param background_palette Color palette for `background_raster`.
#' @param state_lines Logical. If `TRUE`, overlay state boundaries from the
#'   `maps` package. This is intended for US longitude/latitude outputs.
#' @param state_col Color for state boundaries.
#' @param state_lwd Line width for state boundaries.
#' @param n Number of colors used for continuous raster maps.
#' @param legend Logical. If `TRUE`, draw raster and site legends.
#' @param axes Logical. If `TRUE`, draw map axes.
#' @param frame.plot Logical. If `TRUE`, draw a frame around the map panel.
#' @param main Optional plot title. If `NULL`, a title is chosen from `type`.
#' @param col_na Color used for `NA` raster cells.
#' @param boundary_alpha Fill transparency for hard-boundary maps.
#' @param boundary_lwd Line width for hard-boundary outlines.
#' @param pie_radius Radius for empirical-site pies in map units. When `NULL`,
#'   a conservative radius is chosen from the map extent.
#' @param point_cex Point size for `sites = "points"`.
#' @param add Logical. If `TRUE`, add the raster layer to the existing plot.
#' @param ... Additional arguments passed to `terra::plot()`.
#'
#' @return Invisibly returns the `terra::SpatRaster` created by [popmaps_rast()].
#'
#' @examples
#' ex_raster <- raster::aggregate(hija_raster, fact = 240)
#' aps <- popmaps(
#'   input_raster = ex_raster,
#'   input_locs = hija_struc,
#'   surface = "G",
#'   empirical_pt_dist = 0,
#'   num_sites = 5,
#'   num_tested = 2,
#'   popmod = -0.05,
#'   threshold = 0,
#'   ncore = 1
#' )
#' plot_popmaps(aps, ex_raster, hija_struc)
#'
#' @export
plot_popmaps <- function(pop_raster_list,
                         input_raster,
                         input_locs = NULL,
                         type = c("confidence", "boundary", "axis", "ancestry"),
                         style = c("modern", "manuscript"),
                         axis = 1,
                         sites = c("pies", "points", "none"),
                         site_legend = TRUE,
                         boundaries = TRUE,
                         palette = "Viridis",
                         boundary_palette = "Dark 3",
                         background_raster = NULL,
                         background_threshold = NULL,
                         background_palette = "manuscript_background",
                         state_lines = FALSE,
                         state_col = "grey15",
                         state_lwd = 0.8,
                         n = 100,
                         legend = TRUE,
                         axes = TRUE,
                         frame.plot = FALSE,
                         main = NULL,
                         col_na = "transparent",
                         boundary_alpha = 0.20,
                         boundary_lwd = 1.4,
                         pie_radius = NULL,
                         point_cex = 1,
                         add = FALSE,
                         ...) {
  type <- match.arg(type)
  style <- match.arg(style)
  sites <- match.arg(sites)
  if (style == "manuscript") {
    if (missing(palette)) {
      palette <- "manuscript_confidence"
    }
    if (missing(boundary_palette)) {
      boundary_palette <- "manuscript_axes"
    }
    if (is.null(background_raster)) {
      background_raster <- input_raster
    }
    if (missing(n)) {
      n <- 10
    }
    if (missing(boundary_lwd)) {
      boundary_lwd <- 2.2
    }
    if (missing(state_lines)) {
      state_lines <- TRUE
    }
    if (missing(site_legend)) {
      site_legend <- FALSE
    }
  }
  popmaps_check_logical_scalar(site_legend, "`site_legend`")
  popmaps_check_logical_scalar(boundaries, "`boundaries`")
  popmaps_check_logical_scalar(legend, "`legend`")
  popmaps_check_logical_scalar(state_lines, "`state_lines`")
  popmaps_check_logical_scalar(axes, "`axes`")
  popmaps_check_logical_scalar(frame.plot, "`frame.plot`")
  popmaps_check_logical_scalar(add, "`add`")
  popmaps_check_positive_plot_value(n, "`n`", whole = TRUE)
  popmaps_check_positive_plot_value(state_lwd, "`state_lwd`")
  popmaps_check_positive_plot_value(boundary_alpha, "`boundary_alpha`", allow_zero = TRUE)
  popmaps_check_positive_plot_value(boundary_lwd, "`boundary_lwd`")
  popmaps_check_positive_plot_value(point_cex, "`point_cex`")
  if (!is.null(background_threshold)) {
    popmaps_check_finite_scalar(background_threshold, "`background_threshold`")
  }
  if (!is.null(pie_radius)) {
    popmaps_check_positive_plot_value(pie_radius, "`pie_radius`")
  }

  raster_output <- popmaps_rast(pop_raster_list, input_raster)
  num_axes <- max(0L, terra::nlyr(raster_output) - 2L)
  if (num_axes < 1L) {
    stop("`pop_raster_list` must include at least one ancestry-axis layer.", call. = FALSE)
  }

  axis_colors <- popmaps_viz_palette(boundary_palette, num_axes)
  plot_type <- popmaps_viz_resolve_type(type)
  plot_layer <- popmaps_viz_select_layer(raster_output, plot_type, axis)
  plot_main <- popmaps_viz_title(plot_type, plot_layer, main)
  background_layer <- popmaps_prepare_viz_background(background_raster, raster_output)
  confidence_breaks <- if (style == "manuscript") {
    seq(0, 1, length.out = n + 1L)
  } else {
    NULL
  }
  locations <- NULL
  if (!is.null(input_locs) && sites != "none") {
    locations <- popmaps_prepare_locations(input_locs)
    popmaps_validate_site_axis_count(locations, num_axes)
  }

  if (type == "boundary") {
    boundary_layer <- popmaps_mask_viz_layer(
      raster_output[["hard_boundary"]],
      background_layer,
      background_threshold
    )
    popmaps_draw_background_map(
      background = background_layer,
      colors = popmaps_viz_palette(background_palette, max(n, 10L)),
      axes = axes,
      frame.plot = frame.plot,
      main = plot_main,
      col_na = col_na,
      add = add
    )
    popmaps_plot_boundary_map(
      boundary = boundary_layer,
      axis_colors = axis_colors,
      legend = legend,
      axes = if (is.null(background_layer)) axes else FALSE,
      frame.plot = if (is.null(background_layer)) frame.plot else FALSE,
      main = if (is.null(background_layer)) plot_main else "",
      col_na = col_na,
      boundary_alpha = boundary_alpha,
      boundary_lwd = boundary_lwd,
      add = add || !is.null(background_layer),
      ...
    )
  } else {
    display_layer <- popmaps_mask_viz_layer(
      raster_output[[plot_layer]],
      background_layer,
      background_threshold
    )
    popmaps_draw_background_map(
      background = background_layer,
      colors = popmaps_viz_palette(background_palette, max(n, 10L)),
      axes = axes,
      frame.plot = frame.plot,
      main = plot_main,
      col_na = col_na,
      add = add
    )
    popmaps_plot_continuous_map(
      layer = display_layer,
      colors = popmaps_viz_palette(palette, n),
      breaks = confidence_breaks,
      legend = legend,
      axes = if (is.null(background_layer)) axes else FALSE,
      frame.plot = if (is.null(background_layer)) frame.plot else FALSE,
      main = if (is.null(background_layer)) plot_main else "",
      col_na = col_na,
      add = add || !is.null(background_layer),
      ...
    )
    if (isTRUE(boundaries)) {
      popmaps_draw_boundary_outlines(
        boundary = raster_output[["hard_boundary"]],
        axis_colors = axis_colors,
        boundary_lwd = boundary_lwd
      )
    }
  }

  if (isTRUE(state_lines)) {
    popmaps_draw_state_lines(raster_output, state_col = state_col, state_lwd = state_lwd)
  }

  if (!is.null(locations) && sites != "none") {
    popmaps_draw_sites(
      locations = locations,
      axis_colors = axis_colors,
      site_style = sites,
      legend = site_legend,
      pie_radius = pie_radius,
      point_cex = point_cex,
      map_extent = terra::ext(raster_output)
    )
  }

  invisible(raster_output)
}

#' Write a modern POPMAPS map to PNG
#'
#' @description
#' `write_popmaps_plot()` is a small export wrapper around `plot_popmaps()`. It
#' writes a publication- or website-ready PNG while returning the written path.
#'
#' @inheritParams plot_popmaps
#' @param path Output PNG path.
#' @param width,height Plot dimensions in pixels.
#' @param res Output resolution in pixels per inch.
#' @param bg PNG background color.
#' @param overwrite Logical. If `FALSE`, the function stops before replacing an
#'   existing file.
#'
#' @return Normalized output path, invisibly.
#'
#' @examples
#' \dontrun{
#' write_popmaps_plot(aps, ex_raster, "hija-ancestry.png", input_locs = hija_struc)
#' }
#'
#' @export
write_popmaps_plot <- function(pop_raster_list,
                               input_raster,
	                               path,
	                               input_locs = NULL,
                               type = c("confidence", "boundary", "axis", "ancestry"),
	                               style = c("modern", "manuscript"),
	                               axis = 1,
	                               sites = c("pies", "points", "none"),
	                               site_legend = TRUE,
	                               boundaries = TRUE,
	                               palette = "Viridis",
	                               boundary_palette = "Dark 3",
	                               background_raster = NULL,
	                               background_threshold = NULL,
	                               background_palette = "manuscript_background",
	                               state_lines = FALSE,
	                               state_col = "grey15",
	                               state_lwd = 0.8,
	                               n = 100,
                               legend = TRUE,
                               axes = TRUE,
                               frame.plot = FALSE,
                               main = NULL,
                               col_na = "transparent",
                               boundary_alpha = 0.20,
                               boundary_lwd = 1.4,
                               pie_radius = NULL,
                               point_cex = 1,
                               width = 1800,
                               height = 1400,
                               res = 220,
                               bg = "white",
	                               overwrite = FALSE,
	                               ...) {
  style <- match.arg(style)
  if (style == "manuscript") {
    if (missing(palette)) {
      palette <- "manuscript_confidence"
    }
    if (missing(boundary_palette)) {
      boundary_palette <- "manuscript_axes"
    }
    if (missing(background_raster)) {
      background_raster <- NULL
    }
    if (missing(n)) {
      n <- 10
    }
    if (missing(boundary_lwd)) {
      boundary_lwd <- 2.2
    }
    if (missing(state_lines)) {
      state_lines <- TRUE
    }
    if (missing(site_legend)) {
      site_legend <- FALSE
    }
  }
  popmaps_check_character_scalar(path, "`path`")
  popmaps_check_png_dimension(width, "`width`")
  popmaps_check_png_dimension(height, "`height`")
  popmaps_check_png_dimension(res, "`res`")
  popmaps_check_logical_scalar(overwrite, "`overwrite`")
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop("`path` already exists. Use `overwrite = TRUE` to replace it.", call. = FALSE)
  }

  output_dir <- dirname(path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(output_dir)) {
    stop("Could not create the output directory for `path`.", call. = FALSE)
  }

  grDevices::png(path, width = width, height = height, res = res, bg = bg)
  on.exit(grDevices::dev.off(), add = TRUE)
  plot_popmaps(
    pop_raster_list = pop_raster_list,
    input_raster = input_raster,
	    input_locs = input_locs,
	    type = type,
	    style = style,
	    axis = axis,
	    sites = sites,
	    site_legend = site_legend,
	    boundaries = boundaries,
	    palette = palette,
	    boundary_palette = boundary_palette,
	    background_raster = background_raster,
	    background_threshold = background_threshold,
	    background_palette = background_palette,
	    state_lines = state_lines,
	    state_col = state_col,
	    state_lwd = state_lwd,
	    n = n,
    legend = legend,
    axes = axes,
    frame.plot = frame.plot,
    main = main,
    col_na = col_na,
    boundary_alpha = boundary_alpha,
    boundary_lwd = boundary_lwd,
    pie_radius = pie_radius,
    point_cex = point_cex,
    add = FALSE,
    ...
  )

  invisible(normalizePath(path, mustWork = FALSE))
}

popmaps_plot_continuous_map <- function(layer,
                                        colors,
                                        breaks = NULL,
                                        legend,
                                        axes,
                                        frame.plot,
                                        main,
                                        col_na,
                                        add,
                                        ...) {
  plot_args <- list(
    x = layer,
    col = colors,
    range = c(0, 1),
    colNA = col_na,
    axes = axes,
    legend = legend,
    frame = frame.plot,
    main = main,
    add = add
  )
  if (!is.null(breaks)) {
    plot_args$breaks <- breaks
  }
  extra <- list(...)
  plot_args[names(extra)] <- extra
  do.call(terra::plot, plot_args)
}

popmaps_draw_background_map <- function(background,
                                        colors,
                                        axes,
                                        frame.plot,
                                        main,
                                        col_na,
                                        add) {
  if (is.null(background)) {
    return(invisible(FALSE))
  }

  plot_args <- list(
    x = background,
    col = colors,
    range = c(0, 1),
    colNA = col_na,
    axes = axes,
    legend = FALSE,
    frame = frame.plot,
    main = main,
    add = add
  )
  do.call(terra::plot, plot_args)
  invisible(TRUE)
}

popmaps_plot_boundary_map <- function(boundary,
                                      axis_colors,
                                      legend,
                                      axes,
                                      frame.plot,
                                      main,
                                      col_na,
                                      boundary_alpha,
                                      boundary_lwd,
                                      add,
                                      ...) {
  num_axes <- length(axis_colors)
  fill_colors <- grDevices::adjustcolor(axis_colors, alpha.f = min(boundary_alpha, 1))
  plot_args <- list(
    x = boundary,
    col = fill_colors,
    breaks = seq(0.5, num_axes + 0.5, by = 1),
    colNA = col_na,
    axes = axes,
    legend = FALSE,
    frame = frame.plot,
    main = main,
    add = add
  )
  extra <- list(...)
  plot_args[names(extra)] <- extra
  do.call(terra::plot, plot_args)
  popmaps_draw_boundary_outlines(boundary, axis_colors, boundary_lwd = boundary_lwd)
  if (isTRUE(legend)) {
    popmaps_draw_axis_legend(axis_colors, terra::ext(boundary))
  }
}

popmaps_draw_boundary_outlines <- function(boundary, axis_colors, boundary_lwd) {
  polygons <- tryCatch(
    terra::as.polygons(boundary, dissolve = TRUE, values = TRUE, na.rm = TRUE),
    error = function(err) NULL
  )
  if (is.null(polygons) || nrow(polygons) < 1L) {
    return(invisible(NULL))
  }

  field <- names(boundary)[1]
  values <- polygons[[field]][, 1]
  for (axis_id in seq_along(axis_colors)) {
    keep <- !is.na(values) & as.integer(values) == axis_id
    if (any(keep)) {
      terra::plot(
        polygons[keep, ],
        add = TRUE,
        col = NA,
        border = axis_colors[[axis_id]],
        lwd = boundary_lwd
      )
    }
  }

  invisible(NULL)
}

popmaps_draw_sites <- function(locations,
                               axis_colors,
                               site_style,
                               legend,
                               pie_radius,
                               point_cex,
                               map_extent) {
  ancestry <- as.matrix(locations[, seq.int(4, ncol(locations)), drop = FALSE])

  if (site_style == "pies") {
    radius <- pie_radius %||% popmaps_default_pie_radius(map_extent)
    for (idx in seq_len(nrow(locations))) {
      slices <- ancestry[idx, ]
      if (sum(slices) <= 0) {
        slices <- rep(1, length(slices))
      }
      plotrix::floating.pie(
        xpos = locations$V2[[idx]],
        ypos = locations$V3[[idx]],
        x = slices,
        radius = radius,
        col = axis_colors,
        border = "grey15"
      )
    }
  } else if (site_style == "points") {
    dominant <- max.col(ancestry, ties.method = "first")
    purity <- apply(ancestry, 1, max, na.rm = TRUE)
    purity[!is.finite(purity)] <- 0
    graphics::points(
      locations$V2,
      locations$V3,
      pch = 21,
      bg = axis_colors[dominant],
      col = "grey10",
      cex = point_cex * (0.85 + 0.45 * pmin(purity, 1)),
      lwd = 0.8
    )
  }

  if (isTRUE(legend)) {
    popmaps_draw_axis_legend(axis_colors, map_extent)
  }
  invisible(NULL)
}

popmaps_draw_axis_legend <- function(axis_colors, map_extent) {
  x_span <- terra::xmax(map_extent) - terra::xmin(map_extent)
  y_span <- terra::ymax(map_extent) - terra::ymin(map_extent)
  legend_x <- terra::xmin(map_extent) + 0.025 * x_span
  legend_y <- terra::ymax(map_extent) - 0.025 * y_span
  graphics::legend(
    x = legend_x,
    y = legend_y,
    legend = paste("Axis", seq_along(axis_colors)),
    pt.bg = axis_colors,
    col = "grey10",
    pch = 21,
    pt.cex = 1.1,
    bg = grDevices::adjustcolor("white", alpha.f = 0.82),
    box.col = "transparent",
    xjust = 0,
    yjust = 1,
    cex = 0.78,
    x.intersp = 0.7
  )
}

popmaps_viz_select_layer <- function(raster_output, type, axis) {
  if (type == "confidence") {
    return("dominant_ancestry_confidence")
  }
  if (type == "boundary") {
    return("hard_boundary")
  }

  if (is.numeric(axis) && length(axis) == 1L && is.finite(axis)) {
    axis <- as.integer(axis)
    layer <- paste0("axis_", axis)
  } else if (is.character(axis) && length(axis) == 1L && !is.na(axis) && nzchar(axis)) {
    layer <- axis
  } else {
    stop("`axis` must be one ancestry-axis number or layer name.", call. = FALSE)
  }

  axis_layers <- names(raster_output)[startsWith(names(raster_output), "axis_")]
  if (!layer %in% axis_layers) {
    stop("`axis` does not identify an available ancestry-axis layer.", call. = FALSE)
  }
  layer
}

popmaps_viz_title <- function(type, layer, main) {
  if (!is.null(main)) {
    return(main)
  }
  switch(
    type,
    confidence = "Dominant ancestry confidence",
    boundary = "Hard ancestry boundary",
    axis = paste("Ancestry", gsub("_", " ", layer))
  )
}

popmaps_viz_resolve_type <- function(type) {
  if (identical(type, "ancestry")) {
    return("confidence")
  }

  type
}

popmaps_viz_palette <- function(palette, n) {
  if (is.character(palette) && length(palette) == 1L && !is.na(palette) && nzchar(palette)) {
    special <- popmaps_special_viz_palette(palette, n)
    if (!is.null(special)) {
      return(special)
    }
    palette_name <- popmaps_match_hcl_palette(palette)
    return(grDevices::hcl.colors(n, palette_name))
  }
  if (!is.character(palette) || length(palette) < 1L || anyNA(palette)) {
    stop("Palette values must be color names, hex colors, or one HCL palette name.", call. = FALSE)
  }
  grDevices::colorRampPalette(palette)(n)
}

popmaps_match_hcl_palette <- function(palette) {
  hcl_palettes <- grDevices::hcl.pals()
  exact <- match(tolower(palette), tolower(hcl_palettes))
  if (!is.na(exact)) {
    return(hcl_palettes[[exact]])
  }

  aliases <- c(
    viridis = "Viridis",
    magma = "Inferno",
    inferno = "Inferno",
    plasma = "Plasma",
    mako = "Mako",
    rocket = "Rocket",
    batlow = "Batlow",
    blueyellow = "Blue-Yellow",
    purpleyellow = "Purple-Yellow",
    bluegreen = "BluGrn",
    dark = "Dark 3"
  )
  key <- gsub("[^a-z0-9]", "", tolower(palette))
  if (!is.na(aliases[[key]])) {
    return(aliases[[key]])
  }

  stop("Unknown HCL palette `", palette, "`.", call. = FALSE)
}

popmaps_special_viz_palette <- function(palette, n) {
  key <- gsub("[^a-z0-9]", "", tolower(palette))
  if (key %in% c(
    "manuscriptconfidence",
    "confidencegray",
    "confidencegrey",
    "manuscriptprobability",
    "probabilitygray",
    "probabilitygrey"
  )) {
    return(grDevices::gray.colors(n, start = 1, end = 0.08))
  }
  if (key %in% c("manuscriptbackground", "habitatbackground")) {
    return(grDevices::colorRampPalette(c("#d8be8d", "#6e665a"))(n))
  }
  if (key %in% c("manuscriptaxes", "popmapsmanuscript")) {
    manuscript_colors <- c("#f5df2e", "#4b155f", "#159b9b")
    return(rep(manuscript_colors, length.out = n))
  }
  NULL
}

popmaps_prepare_viz_background <- function(background_raster, raster_output) {
  if (is.null(background_raster)) {
    return(NULL)
  }

  background <- popmaps_prepare_raster(background_raster)$rast
  template <- raster_output[[1]]
  if (!terra::compareGeom(background, template, stopOnError = FALSE)) {
    background <- terra::resample(background, template, method = "bilinear")
  }
  names(background) <- "background"
  background
}

popmaps_mask_viz_layer <- function(layer, background, threshold) {
  if (is.null(threshold)) {
    return(layer)
  }
  if (is.null(background)) {
    stop("`background_threshold` requires `background_raster`.", call. = FALSE)
  }
  terra::ifel(background >= threshold, layer, NA)
}

popmaps_draw_state_lines <- function(raster_output, state_col, state_lwd) {
  extent <- terra::ext(raster_output)
  maps::map(
    "state",
    xlim = c(terra::xmin(extent), terra::xmax(extent)),
    ylim = c(terra::ymin(extent), terra::ymax(extent)),
    add = TRUE,
    col = state_col,
    lwd = state_lwd,
    fill = FALSE
  )
  invisible(NULL)
}

popmaps_validate_site_axis_count <- function(locations, num_axes) {
  if (ncol(locations) - 3L != num_axes) {
    stop(
      "`input_locs` ancestry columns must match the number of ancestry-axis layers in `pop_raster_list`.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

popmaps_default_pie_radius <- function(map_extent) {
  x_span <- terra::xmax(map_extent) - terra::xmin(map_extent)
  y_span <- terra::ymax(map_extent) - terra::ymin(map_extent)
  0.012 * max(x_span, y_span)
}

popmaps_check_positive_plot_value <- function(x, label, allow_zero = FALSE, whole = FALSE) {
  popmaps_check_finite_scalar(x, label)
  if (whole && abs(x - round(x)) > sqrt(.Machine$double.eps)) {
    stop(label, " must be a whole number.", call. = FALSE)
  }
  if (allow_zero) {
    if (x < 0) {
      stop(label, " must be non-negative.", call. = FALSE)
    }
  } else if (x <= 0) {
    stop(label, " must be positive.", call. = FALSE)
  }
  invisible(TRUE)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
