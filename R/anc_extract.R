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
  
  if(is.character(input_raster) == F) {
    raster_surface <- input_raster
  } else {
    raster_surface <- raster::raster(input_raster)
  }
  cell_size <- raster::res(raster_surface)[1]
  species_data <- input_locs
  num_axes <- length(species_data[1,])-3
  
  nrows <- raster_surface@nrows
  ncols <- raster_surface@ncols
  ymax <- raster_surface@extent@ymax
  xmin <- raster_surface@extent@xmin
  
  i <- (ymax-dec_lat)/cell_size 
  if (i<1) {
    i <- 1
  } else {
    i <- trunc(i)
  }
  j <- (dec_long-xmin)/cell_size
  if (j<1) {
    j <-1
  } else {
    j <- trunc(j)
  }
 
  anc_sum <- 0
  for(m in 3:length(pop_raster_list)) {
    anc_sum <- anc_sum + pop_raster_list[[m]][i,j]
  }
  anc_vec <- NULL
  for(k in 3:length(pop_raster_list)) {
   anc_vec <- c(anc_vec,pop_raster_list[[k]][i,j]/anc_sum)
  }
  return(anc_vec)
}
  
