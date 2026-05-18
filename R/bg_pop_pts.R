#' @title Create background points for genetically defined populations
#'
#' @description This function generates random points across a raster surface for the purpose
#'     of exploring the environmental variation across a defined geographic extent. Returns
#'     a list with as many elements as genetically defined populations. Each element is a SpatialPoints
#'     object containing the random points that fall within one population as defined by popmaps().
#' @param input_raster An R RasterLayer object defining the geographic extent for the 
#'     spatial interpolation. 
#' @param input_locs An R object (rows = total # empirical sites, columns = total # genetic axes + 3) 
#'     with column 1: site name; column 2: decimal longitude; column 3: decimal latitude; 
#'     column 4…column x: ancestry coefficients for genetic axis 1…genetic axis x. 
#'     Function depends on this precise format – see example data hija_struc.
#' @param pop_raster_list An R object resulting from executing the function popmaps().
#' @param bg_pts The number of background points to generate across the entire surface.
#' @param crs A string defining a mapping projection. The default defines the Albers Equal Area 
#'     Conic projection suitable for the contiguous United States.
#' @references Massatti R & Winkler DE. (2022) Spatially explicit management of genetic diversity using
#'     ancestry probability surfaces. Methods in Ecology and Evolution. http://dx.doi.org/10.1111/2041-210X.13902
#' @author Rob Massatti
#' @examples
#' \dontrun{
#'     ex_raster <- raster::aggregate(hija_raster,fact=16)  #Cells in embedded raster are aggregated to reduce computation time
#'     pp <- popmaps(input_raster=ex_raster,input_locs=hija_struc,empirical_pt_dist=5,num_sites=15,num_tested=4,popmod=-0.05,threshold=0,surface='G')
#'     bg_pts <- bg_pop_pts(pop_raster_list = pp, input_locs = hija_struc, input_raster = ex_raster, bg_pts=1000)
#'     plot(ex_raster)
#'     points(bg_pts$bg1, col='red', pch=19)
#'     points(bg_pts$bg2, col='blue', pch=19)
#'     points(bg_pts$bg3, col='yellow', pch=19)
#' }
#' @export

bg_pop_pts <- function(pop_raster_list='', input_locs='', input_raster='', bg_pts=1000, crs="+init=epsg:5070") {
  raster_surface <- popmaps_prepare_raster(input_raster)$raster
  input_locs <- popmaps_prepare_locations(input_locs)
  popmaps_check_whole_number(bg_pts, "`bg_pts`")

  cell_size <- raster::res(raster_surface)[1]
  nrows <- raster_surface@nrows
  ncols <- raster_surface@ncols
  ymax <- raster_surface@extent@ymax
  xmin <- raster_surface@extent@xmin
  
  h_boundary <- raster::raster(pop_raster_list[[1]],xmn=xmin,xmx=xmin+(cell_size*ncols),ymn=ymax-(cell_size*nrows),ymx=ymax,crs= sp::CRS(crs))
#  hard_boundary <- raster::rasterToPolygons(h_boundary,dissolve=T)
#  hard_boundary <- sp::spTransform(hard_boundary,CRSobj=CRS(crs))
 
  num_axes <- length(input_locs[1,])-3
#  colors_axes <- viridis::viridis(num_axes,begin=0,end=1,direction=-1)
  
  bg <- sp::spsample(as(raster_surface@extent, "SpatialPolygons"),bg_pts,type='random')
 # raster::extract(h_boundary,bg)
  
  bg_list <- list()
  for(w in 1:num_axes) {
    yy<-paste("bg",w,sep='')
    bg_list[[yy]]<- 0
  }
  
  for(j in 1:num_axes) {
    temp <- which(raster::extract(h_boundary, bg) == j)
    bg_pts_temp <- bg[temp,]
    bg_list[[j]] <- bg_pts_temp
  }
  
  return(bg_list)
}
  
