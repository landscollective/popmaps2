dir.create("inst/extdata/validation", recursive = TRUE, showWarnings = FALSE)

required_packages <- c("doParallel", "foreach", "raster", "sp")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop(
    "Missing packages needed to create validation reference: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

for (pkg in required_packages) {
  library(pkg, character.only = TRUE)
}

load("data/hija_raster.rda")
load("data/hija_struc.rda")

baseline_source <- system2(
  "git",
  c("show", "7bf9812:R/popmaps.R"),
  stdout = TRUE,
  stderr = TRUE
)
if (!length(baseline_source)) {
  stop("Could not read the baseline popmaps() source from git.", call. = FALSE)
}

baseline_source <- sub(
  "^popmaps <- function",
  ".popmaps103 <- function",
  baseline_source
)
cluster_line <- grep("parallel::makeCluster", baseline_source, fixed = TRUE)
if (length(cluster_line) != 1) {
  stop("Could not find the baseline cluster setup block.", call. = FALSE)
}
baseline_source <- c(
  baseline_source[seq_len(cluster_line - 1)],
  "  foreach::registerDoSEQ()",
  baseline_source[(cluster_line + 6):length(baseline_source)]
)

baseline_file <- tempfile(fileext = ".R")
writeLines(baseline_source, baseline_file)
source(baseline_file, local = TRUE)

params <- list(
  aggregate_fact = 240,
  surface = "G",
  empirical_pt_dist = 0,
  num_sites = 5,
  num_tested = 2,
  popmod = -0.05,
  threshold = 0
)

ex_raster <- raster::aggregate(hija_raster, fact = params$aggregate_fact)
reference_result <- .popmaps103(
  input_raster = ex_raster,
  input_locs = hija_struc,
  surface = params$surface,
  empirical_pt_dist = params$empirical_pt_dist,
  num_sites = params$num_sites,
  num_tested = params$num_tested,
  popmod = params$popmod,
  threshold = params$threshold,
  ncore = 1
)

reference <- list(
  metadata = list(
    source = "POPMAPS 1.03 popmaps() implementation from commit 7bf9812",
    execution = "Serial foreach backend used to avoid platform-specific PSOCK cluster setup.",
    created_with = R.version.string,
    created_at = as.character(Sys.time()),
    dataset = "Bundled Hilaria jamesii example data"
  ),
  params = params,
  result = reference_result
)

saveRDS(
  reference,
  file = "inst/extdata/validation/popmaps-1.03-hija-small-reference.rds",
  version = 2
)
