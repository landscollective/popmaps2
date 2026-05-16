#' @title Visualize an ancestry probability surface
#'
#' @description This function visualizes the results from the function popmaps(). Maps
#'     can be drawn with hard boundaries only or with hard boundaries plus estimations
#'     of ancestry probabilities (see maptype). Pie charts representing genetic patterns
#'     at empirical sampling sites are drawn in all maps. 
#' @param input_raster An R RasterLayer object defining the geographic extent for the 
#'     spatial interpolation. 
#' @param input_locs An R object (rows = total # empirical sites, columns = total # genetic axes + 3) 
#'     with column 1: site name; column 2: decimal longitude; column 3: decimal latitude; 
#'     column 4…column x: ancestry coefficients for genetic axis 1…genetic axis x. 
#'     Function depends on this precise format – see example data hija_struc.
#' @param pop_raster_list An R object resulting from executing the function popmaps().
#' @param maptype A string (either ‘bound’ or ‘ancestry’) that defines the type of map to be drawn. 
#'     Specifying ‘bound’ will draw the hard boundaries only, while specifying ‘ancestry’ will 
#'     draw the hard boundaries and ancestry probability surface.
#' @param pie_radius A float modifying the size of the pie charts depicting empirical ancestry 
#'     coefficients drawn on top of the probability surface (see Fig. 3 in Massatti & Winkler 2022).
#' @param boundary_width A float modifying the hard boundaries drawn on top of a probability surface. 
#'     The scale of the analysis may require this variable to be modified so that adjacent hard 
#'     boundaries touch one another (see Fig. 3 in Massatti & Winkler 2022).
#' @param crs A string defining a mapping projection. The default defines the Albers Equal Area 
#'     Conic projection suitable for the contiguous United States.
#' @references Massatti R & Winkler DE. (2022) Spatially explicit management of genetic diversity using 
#'     ancestry probability surfaces. Methods in Ecology and Evolution. http://dx.doi.org/10.1111/2041-210X.13902
#' @author Rob Massatti
#' @examples
#'     ex_raster <- raster::aggregate(hija_raster,fact=16)  #Cells in embedded raster are aggregated to reduce computation time
#'     pp <- popmaps(input_raster=ex_raster,input_locs=hija_struc,empirical_pt_dist=5,num_sites=15,num_tested=4,popmod=-0.05,threshold=0,surface='G')
#'     popmap_viz(pop_raster_list=pp,input_locs=hija_struc,input_raster=ex_raster,maptype='ancestry',crs="+init=epsg:5070",boundary_width=-0.01)
#' @export

popmap_viz <- function(pop_raster_list='',input_locs='',input_raster='',maptype=c('bound','ancestry'), pie_radius= 0.15, boundary_width= -0.015,crs="+init=epsg:5070") {

  maptype <- match.arg(maptype)
  if (!requireNamespace("rgeos", quietly = TRUE)) {
    stop("The 'rgeos' package is required for legacy boundary plotting. This will be replaced in a future sf/terra plotting path.", call. = FALSE)
  }
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("The 'viridis' package is required for popmap_viz().", call. = FALSE)
  }

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
	hard_boundary <- raster::rasterToPolygons(h_boundary,dissolve=T)
	hard_boundary <- sp::spTransform(hard_boundary,CRSobj=sp::CRS(crs))

	ancest_surface <- raster::raster(pop_raster_list[[2]],xmn=xmin,xmx=xmin+(cell_size*ncols),ymn=ymax-(cell_size*nrows),ymx=ymax,crs= sp::CRS(crs))
	
	num_axes <- length(input_locs[1,])-3
	colors_axes <- viridis::viridis(num_axes,begin=0,end=1,direction=-1)

	if(maptype == 'bound') {
		temp_raster <- h_boundary
		colors <- grey.colors(num_axes,start=1,end=0.3)
		raster::plot(temp_raster,box=F,col=colors,legend=F)
		maps::map("state", xlim=c(temp_raster@extent[1],temp_raster@extent[2]), ylim=c(temp_raster@extent[3],temp_raster@extent[4]), add=T, col='tan',lwd = 2, fill=T)
		raster::plot(temp_raster,box=F,col=colors,legend=F, alpha=0.7,add=T)
		for(n in 1:length(hard_boundary@data[,1])) {
			raster::plot(rgeos::gBuffer(hard_boundary[n,],byid=T,width=boundary_width),add=T,border=colors_axes[n],cex=1,lwd=3)
		}
		for(i in 1:length(input_locs$V2)) { 
			plotrix::floating.pie(input_locs[i,2],input_locs[i,3],as.numeric(input_locs[i,4:(4+(num_axes-1))]),radius=pie_radius,col=colors_axes)
		}	
	} else if (maptype == 'ancestry'){
		temp_raster <- ancest_surface
		#temp_raster[temp_raster[]==-9999] <-0
		breakpoints <- c(0,.05,.10,.15,.20,.25,.30,.35,.40,.45,.50,.55,.60,.65,.70,.75,.80,.85,.90,.95,1.00)
		colors <- grey.colors(20,start=1,end=0.3)
		raster::plot(temp_raster,box=F,breaks=breakpoints,col=colors)
		#plot(temp_raster,box=F,col=colors)
		maps::map("state", xlim=c(temp_raster@extent[1],temp_raster@extent[2]), ylim=c(temp_raster@extent[3],temp_raster@extent[4]), add=T, col='tan',lwd = 2, fill=T)
		raster::plot(temp_raster,box=F,breaks=breakpoints,col=colors,alpha=0.85,add=T)
		#plot(temp_raster,box=F,col=colors,alpha=0.85,add=T)
		for(n in 1:length(hard_boundary@data[,1])) {
			raster::plot(rgeos::gBuffer(hard_boundary[n,],byid=T,width=boundary_width),add=T,border=colors_axes[n],cex=1,lwd=3)
		}
		for(i in 1:length(input_locs$V2)) { 
			plotrix::floating.pie(input_locs[i,2],input_locs[i,3],as.numeric(input_locs[i,4:(4+(num_axes-1))]),radius=pie_radius,col=colors_axes)
		}
	}
}
