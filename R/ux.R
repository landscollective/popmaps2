popmaps_check_quiet <- function(quiet) {
  if (!is.logical(quiet) || length(quiet) != 1 || is.na(quiet)) {
    stop("`quiet` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  invisible(TRUE)
}

popmaps_inform <- function(..., quiet = TRUE) {
  if (!isTRUE(quiet)) {
    message(..., appendLF = TRUE)
  }

  invisible(TRUE)
}

popmaps_progress_message <- function(index, total, label, quiet = TRUE) {
  if (isTRUE(quiet) || total < 1) {
    return(invisible(TRUE))
  }

  interval <- max(1L, floor(total / 10))
  if (index == 1L || index == total || index %% interval == 0L) {
    message(label, " ", index, " of ", total, ".", appendLF = TRUE)
  }

  invisible(TRUE)
}

popmaps_site_label <- function(locations, idx) {
  site <- as.character(locations$V1[idx])
  paste0(site, " (row ", idx, ")")
}

popmaps_collapse_examples <- function(values, max_values = 5) {
  values <- as.character(values)
  shown <- utils::head(values, max_values)
  suffix <- if (length(values) > max_values) {
    paste0(", and ", length(values) - max_values, " more")
  } else {
    ""
  }

  paste0(paste(shown, collapse = ", "), suffix)
}
