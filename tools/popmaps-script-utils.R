truthy_env <- function(x) {
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

popmaps_numeric_env <- function(env, default = NA_real_) {
  value <- suppressWarnings(as.numeric(Sys.getenv(env, unset = as.character(default))))
  if (length(value) != 1 || is.na(value) || !is.finite(value)) {
    return(default)
  }

  value
}

popmaps_integer_env <- function(env, default = NA_integer_) {
  value <- popmaps_numeric_env(env, default)
  if (length(value) != 1 || is.na(value) || value != floor(value)) {
    return(default)
  }

  as.integer(value)
}

popmaps_first_positive_integer <- function(values) {
  values <- suppressWarnings(as.integer(values))
  values <- values[is.finite(values) & values > 0]
  if (length(values) < 1) {
    return(NA_integer_)
  }

  values[[1]]
}

popmaps_command_output <- function(command, args = character()) {
  output <- tryCatch(
    suppressWarnings(system2(command, args = args, stdout = TRUE, stderr = FALSE)),
    error = function(err) character()
  )
  output <- trimws(output)
  output[nzchar(output)]
}

popmaps_available_threads <- function() {
  env_candidates <- c(
    "POPMAPS_AVAILABLE_THREADS",
    "SLURM_CPUS_PER_TASK",
    "PBS_NP",
    "NSLOTS",
    "NUMBER_OF_PROCESSORS"
  )
  env_values <- vapply(env_candidates, Sys.getenv, character(1), unset = "")
  env_threads <- popmaps_first_positive_integer(env_values[nzchar(env_values)])
  if (!is.na(env_threads)) {
    return(env_threads)
  }

  sysname <- Sys.info()[["sysname"]]
  sysname <- if (is.na(sysname)) "" else tolower(sysname)
  command_values <- character()

  if (identical(sysname, "darwin")) {
    command_values <- c(command_values, popmaps_command_output("sysctl", c("-n", "hw.logicalcpu")))
  }
  if (!identical(sysname, "windows")) {
    command_values <- c(
      command_values,
      popmaps_command_output("getconf", "_NPROCESSORS_ONLN"),
      popmaps_command_output("nproc")
    )
  }

  command_threads <- popmaps_first_positive_integer(command_values)
  if (!is.na(command_threads)) {
    return(command_threads)
  }

  1L
}

popmaps_read_thread_request <- function(env_names, available_threads) {
  for (env in env_names) {
    raw_value <- Sys.getenv(env, unset = "")
    if (!nzchar(raw_value)) {
      next
    }

    value <- tolower(trimws(raw_value))
    if (value %in% c("auto", "default")) {
      return(NULL)
    }
    if (value %in% c("all", "max")) {
      return(available_threads)
    }

    requested <- popmaps_first_positive_integer(value)
    if (!is.na(requested)) {
      return(min(requested, available_threads))
    }

    stop(env, " must be a positive whole number, 'auto', or 'all'.", call. = FALSE)
  }

  NULL
}

popmaps_default_threads <- function(available_threads) {
  fraction <- popmaps_numeric_env("POPMAPS_THREAD_FRACTION", 0.75)
  if (!is.finite(fraction) || fraction <= 0 || fraction > 1) {
    stop("POPMAPS_THREAD_FRACTION must be > 0 and <= 1.", call. = FALSE)
  }

  if (available_threads <= 2) {
    return(1L)
  }

  max(1L, min(available_threads - 1L, floor(available_threads * fraction)))
}

popmaps_set_env_if_unset <- function(env, value) {
  if (!nzchar(Sys.getenv(env, unset = ""))) {
    do.call(Sys.setenv, as.list(stats::setNames(as.character(value), env)))
  }
}

popmaps_configure_script_resources <- function(script_prefix = NULL, quiet = FALSE) {
  available_threads <- popmaps_available_threads()
  prefix_thread_env <- if (is.null(script_prefix)) character() else paste0(script_prefix, "_THREADS")
  requested_threads <- popmaps_read_thread_request(
    env_names = c(prefix_thread_env, "POPMAPS_THREADS"),
    available_threads = available_threads
  )
  threads <- if (is.null(requested_threads)) {
    popmaps_default_threads(available_threads)
  } else {
    requested_threads
  }

  terra_memfrac <- popmaps_numeric_env("POPMAPS_TERRA_MEMFRAC", 0.70)
  if (!is.finite(terra_memfrac) || terra_memfrac <= 0 || terra_memfrac > 0.95) {
    stop("POPMAPS_TERRA_MEMFRAC must be > 0 and <= 0.95.", call. = FALSE)
  }

  temp_dir <- Sys.getenv("POPMAPS_TMPDIR", unset = tempdir())
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  temp_dir <- normalizePath(temp_dir, mustWork = TRUE)

  thread_env <- c(
    "OMP_NUM_THREADS",
    "OMP_THREAD_LIMIT",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "RCPP_PARALLEL_NUM_THREADS",
    "GDAL_NUM_THREADS"
  )
  for (env in thread_env) {
    popmaps_set_env_if_unset(env, threads)
  }
  for (env in c("TMPDIR", "TMP", "TEMP")) {
    popmaps_set_env_if_unset(env, temp_dir)
  }

  options(Ncpus = threads)

  if (requireNamespace("terra", quietly = TRUE)) {
    terra::terraOptions(memfrac = terra_memfrac, tempdir = temp_dir)
  }

  config <- list(
    available_threads = available_threads,
    threads = threads,
    thread_fraction = popmaps_numeric_env("POPMAPS_THREAD_FRACTION", 0.75),
    terra_memfrac = terra_memfrac,
    temp_dir = temp_dir
  )
  class(config) <- "popmaps_script_resources"

  if (!isTRUE(quiet)) {
    message(
      "Resource configuration: using ",
      threads,
      " of ",
      available_threads,
      " detected logical processors; terra memfrac = ",
      terra_memfrac,
      "; tempdir = ",
      temp_dir
    )
  }

  config
}

popmaps_resource_row <- function(resource_config) {
  data.frame(
    available_threads = resource_config$available_threads,
    threads = resource_config$threads,
    thread_fraction = resource_config$thread_fraction,
    terra_memfrac = resource_config$terra_memfrac,
    temp_dir = resource_config$temp_dir,
    stringsAsFactors = FALSE
  )
}
