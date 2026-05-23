#' @title Legacy POPMAPS 1.03 ancestry-surface visualization
#'
#' @description This function preserves the original POPMAPS plotting interface
#'     for compatibility with POPMAPS 1.03-era scripts. New analyses should use
#'     [plot_popmaps()] or [write_popmaps_plot()], which provide terra-based map
#'     output, manuscript-style defaults, optional background rasters, and PNG
#'     export. `popmap_viz()` can draw hard boundaries only or hard boundaries
#'     with ancestry probabilities. Pie charts representing empirical ancestry
#'     patterns are drawn at sampling sites.
#' @param input_raster An R RasterLayer object defining the geographic extent for
#'     the spatial interpolation.
#' @param input_locs An R object with rows as empirical sites and columns as:
#'     column 1, site name; column 2, decimal longitude or x-coordinate; column
#'     3, decimal latitude or y-coordinate; columns 4 through n, ancestry
#'     coefficients for each genetic axis. The function depends on this precise
#'     column order; see the example data `hija_struc`.
#' @param pop_raster_list An R object resulting from executing [popmaps()].
#' @param maptype A string, either `"bound"` or `"ancestry"`, that defines the
#'     map type. `"bound"` draws hard boundaries only. `"ancestry"` draws hard
#'     boundaries and the maximum ancestry-probability surface.
#' @param pie_radius A numeric value modifying the size of the pie charts
#'     depicting empirical ancestry coefficients drawn on top of the probability
#'     surface.
#' @param boundary_width Retained for compatibility with POPMAPS 1.03. Boundary buffering
#'     previously depended on retired spatial packages and is no longer applied.
#' @param crs A string defining a mapping projection. The historical default is
#'     an Albers Equal Area Conic projection used by older POPMAPS examples.
#' @references Massatti R & Winkler DE. (2022) Spatially explicit management of genetic diversity using 
#'     ancestry probability surfaces. Methods in Ecology and Evolution. http://dx.doi.org/10.1111/2041-210X.13902
#' @author Rob Massatti
#' @examples
#' \dontrun{
#'     ex_raster <- raster::aggregate(hija_raster,fact=16)  #Cells in embedded raster are aggregated to reduce computation time
#'     pp <- popmaps(input_raster=ex_raster,input_locs=hija_struc,empirical_pt_dist=5,num_sites=15,num_tested=4,popmod=-0.05,threshold=0,surface='G')
#'     popmap_viz(pop_raster_list=pp,input_locs=hija_struc,input_raster=ex_raster,maptype='ancestry',crs="+init=epsg:5070",boundary_width=-0.01)
#' }
#' @export

popmap_viz <- function(pop_raster_list='',input_locs='',input_raster='',maptype=c('bound','ancestry'), pie_radius= 0.15, boundary_width= -0.015,crs="+init=epsg:5070") {

  maptype <- match.arg(maptype)
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("The 'viridis' package is required for popmap_viz().", call. = FALSE)
  }

  raster_surface <- popmaps_prepare_raster(input_raster)$raster
  input_locs <- popmaps_prepare_locations(input_locs)

	cell_size <- raster::res(raster_surface)[1]
	nrows <- raster_surface@nrows
	ncols <- raster_surface@ncols
	ymax <- raster_surface@extent@ymax
	xmin <- raster_surface@extent@xmin

	h_boundary <- raster::raster(pop_raster_list[[1]],xmn=xmin,xmx=xmin+(cell_size*ncols),ymn=ymax-(cell_size*nrows),ymx=ymax,crs= sp::CRS(crs))

	ancest_surface <- raster::raster(pop_raster_list[[2]],xmn=xmin,xmx=xmin+(cell_size*ncols),ymn=ymax-(cell_size*nrows),ymx=ymax,crs= sp::CRS(crs))
	
	num_axes <- length(input_locs[1,])-3
	colors_axes <- viridis::viridis(num_axes,begin=0,end=1,direction=-1)

	plot_boundaries <- function() {
		hard_boundary <- raster::rasterToPolygons(h_boundary,dissolve=F)
		pop_ids <- as.integer(hard_boundary@data[[1]])
		for(n in seq_along(pop_ids)) {
			pop_id <- pop_ids[n]
			if(is.na(pop_id) || pop_id < 1 || pop_id > num_axes) {
				next
			}
			raster::plot(hard_boundary[n,],add=T,border=colors_axes[pop_id],col=NA,lwd=3)
		}
	}

	if(maptype == 'bound') {
		temp_raster <- h_boundary
		colors <- grey.colors(num_axes,start=1,end=0.3)
		raster::plot(temp_raster,box=F,col=colors,legend=F)
		maps::map("state", xlim=c(temp_raster@extent[1],temp_raster@extent[2]), ylim=c(temp_raster@extent[3],temp_raster@extent[4]), add=T, col='tan',lwd = 2, fill=T)
		raster::plot(temp_raster,box=F,col=colors,legend=F, alpha=0.7,add=T)
		plot_boundaries()
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
		plot_boundaries()
		for(i in 1:length(input_locs$V2)) { 
			plotrix::floating.pie(input_locs[i,2],input_locs[i,3],as.numeric(input_locs[i,4:(4+(num_axes-1))]),radius=pie_radius,col=colors_axes)
		}
	}
}
