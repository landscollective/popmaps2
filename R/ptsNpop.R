#' @title Partition points by genetically defined populations
#'
#' @description This function assigns points contained in a data frame to a genetically defined population as 
#'     estimated by the function popmaps(). It returns a list with as many elements as genetically defined populations;
#'     each element is a data frame containing the points that are assigned to one of the genetically defined populations.
#' @param input_raster An R RasterLayer object defining the geographic extent for the 
#'     spatial interpolation. 
#' @param input_locs An R object (rows = total # empirical sites, columns = total # genetic axes + 3) 
#'     with column 1: site name; column 2: decimal longitude; column 3: decimal latitude; 
#'     column 4…column x: ancestry coefficients for genetic axis 1…genetic axis x. 
#'     Function depends on this precise format – see example data hija_struc.
#' @param pop_raster_list An R object resulting from executing the function popmaps().
#' @param sampling_pts A data frame of sampling points that will be partitioned by population. For example,
#'     a list of herbarium voucher specimen locations - see example data hija_herb.
#' @param crs A string defining a mapping projection. The default defines the Albers Equal Area 
#'     Conic projection suitable for the contiguous United States.
#' @references Massatti R & Winkler DE. (2022) Spatially explicit management of genetic diversity using 
#'     ancestry probability surfaces. Methods in Ecology and Evolution. http://dx.doi.org/10.1111/2041-210X.13902
#' @author Rob Massatti
#' @examples
#'    ex_raster <- raster::aggregate(hija_raster,fact=16)   #Cells in embedded raster are aggregated to reduce computation time
#'    pp <- popmaps(input_raster=ex_raster,input_locs=hija_struc,empirical_pt_dist=5,num_sites=15,num_tested=4,popmod=-0.05,threshold=0,surface='G')
#'    samp_per_pop <- ptsNpop(pop_raster_list=pp, input_locs=hija_struc, input_raster=ex_raster,sampling_pts=hija_herb)
#'    plot(ex_raster)
#'    points(samp_per_pop$pop1, pch=24, col= 'black', bg='red')
#'    points(samp_per_pop$pop2, pch=24, col= 'black', bg='blue')
#'    points(samp_per_pop$pop3, pch=24, col= 'black', bg='yellow')
#' @export

ptsNpop <- function(pop_raster_list='', input_locs='', input_raster='', sampling_pts='', crs="+init=epsg:5070") {
  if(is.character(input_raster) == F) {
    raster_surface <- input_raster
  } else {
    raster_surface <- raster::raster(input_raster)
  }

  cell_size <- raster::res(raster_surface)[1]
  nrows <- raster_surface@nrows
  ncols <- raster_surface@ncols
  ymax <- raster_surface@extent@ymax
  xmin <- raster_surface@extent@xmin
  
  h_boundary <- raster::raster(pop_raster_list[[1]],xmn=xmin,xmx=xmin+(cell_size*ncols),ymn=ymax-(cell_size*nrows),ymx=ymax,crs= sp::CRS(crs))
  #  hard_boundary <- raster::rasterToPolygons(h_boundary,dissolve=T)
  #  hard_boundary <- sp::spTransform(hard_boundary,CRSobj=CRS(crs))
  
  num_axes <- length(input_locs[1,])-3

    sampling_list <- list()
  for(w in 1:num_axes) {
    yy<-paste("pop",w,sep='')
    sampling_list[[yy]]<- 0
  }
  
  for(j in 1:num_axes) {
    
    temp <- which(raster::extract(h_boundary, sampling_pts) == j)
    sampling_pts_temp <- sampling_pts[temp,]
    sampling_list[[j]] <- sampling_pts_temp
  }
  
  return(sampling_list)
}
