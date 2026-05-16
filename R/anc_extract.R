#' @title Extract estimated ancestry coefficients
#'
#' @description This function extracts the estimated ancestry coefficients of a specified 
#'     geographic location from an object resulting from the function popmaps(). 
#' @param input_raster An R RasterLayer object defining the geographic extent of the
#'     interpolation recorded in pop_raster_list. 
#' @param input_locs An R object (rows = total # empirical sites, columns = total # genetic axes + 3) 
#'     with column 1: site name; column 2: decimal longitude; column 3: decimal latitude; 
#'     column 4…column x: ancestry coefficients for genetic axis 1…genetic axis x. 
#'     Function depends on this precise format – see example data hija_struc.
#' @param pop_raster_list An R object resulting from executing the function popmaps().
#' @param dec_lat A float specifying the decimal latitude from which to extract
#'     estimated ancestry coefficients.
#' @param dec_long A float specifying the decimal longitude from which to extract
#'     estimated ancestry coefficients.
#' @references Massatti R & Winkler DE. (2022) Spatially explicit management of genetic diversity using 
#'     ancestry probability surfaces. Methods in Ecology and Evolution. http://dx.doi.org/10.1111/2041-210X.13902
#' @author Rob Massatti
#' @examples
#'     ex_raster <- raster::aggregate(hija_raster,fact=16)  #Cells in embedded raster are aggregated to reduce computation time
#'     pp <- popmaps(input_raster=ex_raster,input_locs=hija_struc,empirical_pt_dist=5,num_sites=15,num_tested=4,popmod=-0.05,threshold=0,surface='G')
#'     anc_extract(pop_raster_list=pp, input_raster=ex_raster, input_locs=hija_struc, dec_lat=39.46522,dec_long=-110.9525)
#' @export

anc_extract <- function(pop_raster_list='',input_raster='', input_locs='', dec_lat='', dec_long='') {

  if (!is.list(pop_raster_list) || length(pop_raster_list) < 3) {
    stop("`pop_raster_list` must be a list returned by popmaps().", call. = FALSE)
  }

  raster_input <- popmaps_prepare_raster(input_raster)
  species_data <- popmaps_prepare_locations(input_locs)
  num_axes <- length(species_data[1,]) - 3

  if (length(pop_raster_list) < (2 + num_axes)) {
    stop("`pop_raster_list` does not contain enough ancestry coefficient surfaces for `input_locs`.", call. = FALSE)
  }

  point <- popmaps_prepare_point(dec_long = dec_long, dec_lat = dec_lat)
  cell <- terra::cellFromXY(raster_input$rast, point)
  if (is.na(cell)) {
    stop("The requested coordinate falls outside `input_raster`.", call. = FALSE)
  }

  row_col <- terra::rowColFromCell(raster_input$rast, cell)
  i <- row_col[, 1]
  j <- row_col[, 2]

  anc_vec <- vapply(
    seq.int(3, 2 + num_axes),
    function(k) pop_raster_list[[k]][i, j],
    numeric(1)
  )
  anc_sum <- sum(anc_vec)

  if (!is.finite(anc_sum) || anc_sum == 0) {
    return(rep(NA_real_, num_axes))
  }

  anc_vec / anc_sum
}
  
