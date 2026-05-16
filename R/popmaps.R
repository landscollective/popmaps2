#' @title Estimate an ancestry probability surface
#'
#' @description This function estimates an ancestry probability surface when supplied
#'     a geospatial layer and empirical genetic data describing patterns of ancestry coefficients
#'     across sampling locations. 
#' @param input_raster An R RasterLayer object defining the geographic extent for the spatial 
#'     interpolation. Values in the cells will be used to calculate distance used in the 
#'     dist_prob_func if surface = ‘C’. See example data hija_raster.
#' @param input_locs An R object (rows = total # empirical sites, columns = total # genetic axes + 3]) 
#'     with column 1: site name; column 2: decimal longitude; column 3: decimal latitude; 
#'     column 4…column x: ancestry coefficients for genetic axis 1…genetic axis x. 
#'     Functions depend on this precise format – see example data hija_struc.
#' @param surface A string (either ‘G’ or ‘C’) that determines how input_raster is used to calculate 
#'     ancestry coefficients. If 'G', the spatial attributes of input_raster will be used to calculate 
#'     geographic (i.e., Euclidean) distances between empirical sites and inference sites. If ‘C’, 
#'     least-cost distances are calculated using the values contained in the raster cells.
#' @param empirical_pt_dist An integer representing the minimum distance that all sites in num_tested 
#'     must be separated by. If num_tested = 3, the closest two empirical sites to the cell are 
#'     determined first. If the two sites are closer than this value, the closest will be kept and the 
#'     second discarded. The next closest site will then be selected and compared to the first; this 
#'     process repeats until the two sites are farther than empirical_pt_dist. Subsequently, the third 
#'     empirical site will be selected using the same process. This rarefaction procedure may reduce 
#'     high spatial autocorrelation expected in genetic patterns among proximate empirical sites that 
#'     may mislead estimations of ancestry coefficients. The parameter should be informed by the 
#'     distribution of empirical sites and the biology of the focal species.
#' @param num_tested An integer representing the number of empirical sampling sites that will contribute 
#'     to the estimation of ancestry coefficients.
#' @param popmod A float used in dist_prob_func that modifies the relationship between the contribution 
#'     of an empirical site to the estimation of ancestry coefficients at a cell and distance (see Fig. 2 
#'     in Massatti & Winkler 2022). Using the default dist_prob_func, values ranging from -0.00001 and -0.9 
#'     have the most influence on the shape of the curve. Modifying dist_prob_func will require users to 
#'     determine how to set popmod accordingly.
#' @param threshold A float (0.0-1.0) determined during species distribution modeling inference. If > 0, 
#'     ancestry coefficients will not be estimated in raster cells with values < threshold, which saves 
#'     computation time by avoiding estimating ancestry coefficients for cells in which the focal species 
#'     would never be expected to occur.
#' @param ncore An integer allowing the user to specify how many cores on the local machine should be 
#'     used to facilitate spatial interpolation.
#' @param dist_prob_func A function defining the relationship between distance and the contribution 
#'     of an empirical site’s ancestry coefficients to the estimation of ancestry coefficients at 
#'     an inference cell. The default equation defines the relationship in Fig. 2 of Massatti & Winkler (2022).
#' @param num_sites An integer representing the pool of empirical sites considered when selecting 
#'     num_tested sites to estimate ancestry coefficients. If empirical sites are highly clustered 
#'     and rarefaction due to empirical_pt_dist causes many to be discarded, this variable will likely 
#'     need to be increased at the cost of computing time.
#' @references Massatti R & Winkler DE. (2022) Spatially explicit management of genetic diversity using 
#'     ancestry probability surfaces. Methods in Ecology and Evolution. http://dx.doi.org/10.1111/2041-210X.13902
#' @author Rob Massatti
#' @examples
#'     ex_raster <- raster::aggregate(hija_raster,fact=16)  #Cells in embedded raster are aggregated to reduce computation time
#'     pp <- popmaps(input_raster=ex_raster,input_locs=hija_struc,empirical_pt_dist=5,num_sites=15,num_tested=4,popmod=-0.05,threshold=0,surface='G')
#' @export


popmaps <- function(input_raster='', input_locs='', surface='G', empirical_pt_dist=5, num_sites=10, num_tested=3, popmod=-0.001, ncore=4, threshold=0, dist_prob_func=function(popmod_temp,distance) {exp(popmod_temp*distance)}) {

  surface <- match.arg(surface, c("G", "C"))
  if (surface == "C" && !requireNamespace("gdistance", quietly = TRUE)) {
    stop("The 'gdistance' package is required when surface = 'C'.", call. = FALSE)
  }
  foreach_packages <- if (surface == "C") {
    c("raster", "gdistance")
  } else {
    "raster"
  }

  cl<-parallel::makeCluster(ncore)
	doParallel::registerDoParallel(cl)
	on.exit({
	  parallel::stopCluster(cl)
	  doParallel::stopImplicitCluster()
	}, add = TRUE)
	
	if(is.character(input_raster) == F) {
		raster_surface <- input_raster
	} else {
		raster_surface <- raster::raster(input_raster)
	}
	species_data <- input_locs
	species_pts <- sp::SpatialPointsDataFrame(species_data[,2:3], species_data)
	sampling_loc_coords <- cbind(species_pts@data$V2,species_pts@data$V3)
	num_emp_sites <- length(sampling_loc_coords[,1])
	num_axes <- length(species_data[1,])-3
	cell_size <- raster::res(raster_surface)[1]
	
	if(surface=='G') {
		species_surface <- raster_surface
	}
	if(surface=='C') {
		species_surface <- gdistance::transition(raster_surface, transitionFunction=mean,directions=8)
		species_surface <- gdistance::geoCorrection(species_surface,type="c",scl=T)
	}

	nrows <- species_surface@nrows
	ncols <- species_surface@ncols
	ymax <- species_surface@extent@ymax
	xmin <- species_surface@extent@xmin
#	pty <- pretty(1:nrows,n=10)
#	pty_iter <- 2
		
	x <- foreach::foreach(i = 1:nrows, .combine='rbind') %:%	#nrows
	  foreach::foreach(j = 1:ncols, .combine='c',.packages=foreach_packages,.inorder=T) %dopar% {	
			
			pop_function <- function(i='',j='') {
				earth.dist<-function(lat1,long1,lat2,long2){
					rad <- pi/180
					a1 <- lat1 * rad
					a2 <- long1 * rad
					b1 <- lat2 * rad
					b2 <- long2 * rad
					dlat <- b1-a1
					dlon<- b2-a2
					a <- (sin(dlat/2))^2 +cos(a1)*cos(b1)*(sin(dlon/2))^2
					c <- 2*atan2(sqrt(a),sqrt(1-a))
					R <- 6378.145
					dist <- R *c
					return(dist)
				}

				y<- ymax+(cell_size/2)-(i*cell_size)
				if(i==nrows) {
					y<- ymax+(cell_size/2)-((nrows-1)*cell_size)
				}

				x<- xmin+(cell_size/2)+(j*cell_size)
				if(j==ncols){
					x <- xmin+(cell_size/2)+((ncols-1)*cell_size)
				}
				point <- cbind(x,y)
				distance_vec<-NULL

				for(a in 1:num_emp_sites) {
					ed <- earth.dist(point[1,2],point[1,1],sampling_loc_coords[a,2],sampling_loc_coords[a,1])
					distance_vec <- c(distance_vec,ed)
				}

				temp_site_data <- species_pts@data
				temp_site_data <- cbind(temp_site_data,distance_vec)
				dist_vec_position <- length(temp_site_data[1,])
				temp_site_data <- temp_site_data[order(temp_site_data[,dist_vec_position]),]
				temp_site_data <- temp_site_data[1:num_sites,]

				######################################################################################
				###If clause that is executed when lanscape-based distances are chosen for inference
				#######################################################################################
				if(surface=="C") {
					if(is.na(raster::extract(raster_surface, cbind(point[,1],point[,2]))) == T) {
					  h_bound <- NA
					  ancest_surface <- NA
					  cell_prob_list <- list()
					  for(w in 1:num_axes) {
					    zz<-paste("cell_prob",w,sep='')
					    cell_prob_list[[zz]]<- NA
					  }
					} else if(raster::extract(raster_surface, cbind(point[,1],point[,2])) < threshold) {

						##########################
						######This section ensures that the hard boundries will not be affected by the use of the threshold value
						site_locs_temp <- temp_site_data[,2:3]
						temp_data <- rbind(as.matrix(point),as.matrix(site_locs_temp))
						cosDist <- gdistance::costDistance(species_surface, temp_data)
						cost_vec <- cosDist[1:num_sites]
						temp_site_data <- cbind(temp_site_data,cost_vec)
						temp_site_data <- temp_site_data[order(temp_site_data[,(dist_vec_position+1)]),]
						temp_site_data <- cbind(temp_site_data,seq(1:num_sites))
						k_list<-list()
						for(w in 1:num_axes) {
							xx<-paste("k",w,sep='')
							k_list[[xx]]<- 0
						}
						match <- as.numeric(temp_site_data[1,2:3])
						k_iter <- 4
						for(u in 1:length(k_list)) {
							k_list[[u]] <- temp_site_data[,k_iter][as.numeric(match[1])==temp_site_data$V2]
							k_iter <- k_iter + 1
						}
						max_axis <- max(unlist(k_list))
						max_axis <- which(k_list==max_axis)
						##########################
						##########################
						
						h_bound <- max_axis
						ancest_surface <- NA
						cell_prob_list <- list()
						for(w in 1:num_axes) {
						  zz<-paste("cell_prob",w,sep='')
						  cell_prob_list[[zz]]<- NA
						}

					} else {
						site_locs_temp <- temp_site_data[,2:3]
						temp_data <- rbind(as.matrix(point),as.matrix(site_locs_temp))
						cosDist <- gdistance::costDistance(species_surface, temp_data)
						cost_vec <- cosDist[1:num_sites]

						temp_site_data <- cbind(temp_site_data,cost_vec)
						temp_site_data <- temp_site_data[order(temp_site_data[,(dist_vec_position+1)]),]
						temp_site_data <- cbind(temp_site_data,seq(1:num_sites))

						match <- temp_site_data[temp_site_data[,(dist_vec_position+2)]==1,2:3]
						match2 <- temp_site_data[temp_site_data[,(dist_vec_position+2)]==2,2:3]

						cosDist_loc_pair <- as.matrix(rbind(match,match2))
						pt_dist <- gdistance::costDistance(species_surface,cosDist_loc_pair)

						iterator <- 3
						while(pt_dist < empirical_pt_dist) {
							match2 <- temp_site_data[temp_site_data[,(dist_vec_position+2)]==iterator,2:3]
							cosDist_loc_pair <- as.matrix(rbind(match,match2))
							pt_dist <- gdistance::costDistance(species_surface,cosDist_loc_pair)
							iterator <- iterator + 1
						}
						empirical_site_list <- rbind(match,match2)

						while(length(empirical_site_list[,1]) < num_tested) {
							match_temp <- temp_site_data[temp_site_data[,(dist_vec_position+2)]==iterator,2:3]

							keep <- T
							for(g in 1:length(empirical_site_list[,1])) {
								match <- empirical_site_list[g,]
								cosDist_loc_pair <- as.matrix(rbind(match,match_temp))
								pt_dist <- gdistance::costDistance(species_surface,cosDist_loc_pair)

								if(pt_dist < empirical_pt_dist) {
									keep <- F
								}
							}
							if(keep == T) {
								empirical_site_list <- rbind(empirical_site_list,match_temp)
							}
							iterator <- iterator + 1
						}

						k_list<-list()
						axis_list <- list()
						cell_prob_list <- list()
						temp_list <- list()
						for(w in 1:num_axes) {
							xx<-paste("k",w,sep='')
							k_list[[xx]]<- 0

							yy<-paste("axis",w,sep='')
							axis_list[[yy]]<- 0

							zz<-paste("cell_prob",w,sep='')
							cell_prob_list[[zz]]<- 0

							aa<-paste("temp_list",w,sep='')
							temp_list[[aa]]<- 0
						}

						match <- as.numeric(temp_site_data[1,2:3])

						k_iter <- 4
						for(u in 1:length(k_list)) {
							k_list[[u]] <- temp_site_data[,k_iter][as.numeric(match[1])==temp_site_data$V2]
							k_iter <- k_iter + 1
						}

						max_axis <- max(unlist(k_list))
						max_axis <- which(k_list==max_axis)
						h_bound <- max_axis###

						site_weight <- NULL
						for(h in 1:num_tested) {
							match <- empirical_site_list[h,]
							temp_list_iter <- 4
							for(o in 1:num_axes) {
								temp_list[[o]] <- temp_site_data[,temp_list_iter][as.numeric(match[1])==temp_site_data$V2]
								temp_list_iter <- temp_list_iter + 1
							}
							for(s in 1:num_axes) {
								axis_list[[s]] <- c(axis_list[[s]], temp_list[[s]])
							}

							distance <- temp_site_data[,(dist_vec_position+1)][as.numeric(match[1])==temp_site_data$V2]
							dist_prob <- dist_prob_func(popmod,distance)
							site_weight <- c(site_weight, dist_prob)
							
						}
						for(t in 1:length(axis_list)) {
							axis_list[[t]] <- axis_list[[t]][2:(num_tested+1)]
						}
						for(aa in 1:length(cell_prob_list)) {
							cell_prob_list[[aa]] <- sum(axis_list[[aa]]*site_weight/num_tested)
						}

						max_avg <- max(unlist(cell_prob_list))/sum(unlist(cell_prob_list))
						if(is.na(max_avg)==T) {
							max_avg<-0
						}

						if(max_avg < .5) {
							max_avg <- 0
						} else {
							max_avg <- ((max_avg) - .5) *2
						}
						ancest_surface <- max_avg
					}	#end if else loop	
				} #end surface == C loop

				#######################################################################################
				###If clause that is executed when Euclidean distances are chosen for inference
				#######################################################################################	
				if(surface=="G") {
					if(is.na(raster::extract(raster_surface, cbind(point[,1],point[,2]))) == T) {
						h_bound <- NA
						ancest_surface <- NA
						cell_prob_list <- list()
						for(w in 1:num_axes) {
						  zz<-paste("cell_prob",w,sep='')
						  cell_prob_list[[zz]]<- NA
						}
					} else if(raster::extract(raster_surface, cbind(point[,1],point[,2])) < threshold) {

						##########################
						######This section ensures that the hard boundries will not be affected by the use of the threshold value
						k_list<-list()
						for(w in 1:num_axes) {
							xx<-paste("k",w,sep='')
							k_list[[xx]]<- 0
						}
						match <- as.numeric(temp_site_data[1,2:3])
						k_iter <- 4
						for(u in 1:length(k_list)) {
							k_list[[u]] <- temp_site_data[,k_iter][as.numeric(match[1])==temp_site_data$V2]
							k_iter <- k_iter + 1
						}
						max_axis <- max(unlist(k_list))
						max_axis <- which(k_list==max_axis)
						##########################
						##########################
						
						h_bound <- max_axis
						ancest_surface <- NA
						cell_prob_list <- list()
						for(w in 1:num_axes) {
						  zz<-paste("cell_prob",w,sep='')
						  cell_prob_list[[zz]]<- NA
						}
					} else {
						site_locs_temp <- temp_site_data[,2:3]
						match <- as.numeric(temp_site_data[1,2:3])
						match2 <- as.numeric(temp_site_data[2,2:3])
						pt_dist <- earth.dist(match[2],match[1],match2[2],match2[1])

						iterator <- 3
						while(pt_dist < empirical_pt_dist) {
							match2 <- as.numeric(temp_site_data[iterator,2:3])
							pt_dist <- earth.dist(match[2],match[1],match2[2],match2[1])
							iterator <- iterator + 1
						}
						empirical_site_list <- rbind(match,match2)

						while(length(empirical_site_list[,1]) < num_tested) {
							match_temp <- as.numeric(temp_site_data[iterator,2:3])

							keep <- T
							for(g in 1:length(empirical_site_list[,1])) {
								match <- empirical_site_list[g,]
								pt_dist <- earth.dist(match[2],match[1],match_temp[2],match_temp[1])

								if(pt_dist < empirical_pt_dist) {
									keep <- F
								}
							}
							if(keep == T) {
								empirical_site_list <- rbind(empirical_site_list,match_temp)
							}
							iterator <- iterator + 1
						}

						k_list<-list()
						axis_list <- list()
						cell_prob_list <- list()
						temp_list <- list()
						for(w in 1:num_axes) {
							xx<-paste("k",w,sep='')
							k_list[[xx]]<- 0

							yy<-paste("axis",w,sep='')
							axis_list[[yy]]<- 0

							zz<-paste("cell_prob",w,sep='')
							cell_prob_list[[zz]]<- 0

							aa<-paste("temp_list",w,sep='')
							temp_list[[aa]]<- 0
						}

						match <- as.numeric(temp_site_data[1,2:3])

						k_iter <- 4
						for(u in 1:length(k_list)) {
							k_list[[u]] <- temp_site_data[,k_iter][as.numeric(match[1])==temp_site_data$V2]
							k_iter <- k_iter + 1
						}

						max_axis <- max(unlist(k_list))
						max_axis <- which(k_list==max_axis)
						h_bound <- max_axis###

						site_weight <- NULL
						for(h in 1:num_tested) {
							match <- empirical_site_list[h,]
							temp_list_iter <- 4
							for(o in 1:num_axes) {
								temp_list[[o]] <- temp_site_data[,temp_list_iter][as.numeric(match[1])==temp_site_data$V2]
								temp_list_iter <- temp_list_iter + 1
							}
							for(s in 1:num_axes) {
								axis_list[[s]] <- c(axis_list[[s]], temp_list[[s]])
							}

							distance <- temp_site_data[,(dist_vec_position)][as.numeric(match[1])==temp_site_data$V2]
							dist_prob <- dist_prob_func(popmod,distance)
							site_weight <- c(site_weight, dist_prob)
							
						}
						for(t in 1:length(axis_list)) {
							axis_list[[t]] <- axis_list[[t]][2:(num_tested+1)]
						}
						for(aa in 1:length(cell_prob_list)) {
							cell_prob_list[[aa]] <- sum(axis_list[[aa]]*site_weight/num_tested)
						}

						max_avg <- max(unlist(cell_prob_list))/sum(unlist(cell_prob_list))
						if(is.na(max_avg)==T) {
							max_avg<-0
						}

						if(max_avg < .5) {
							max_avg <- 0
						} else {
							max_avg <- ((max_avg) - .5) *2
						}
						ancest_surface <- max_avg						
					}	#end if else loop	
				} #end surface == G loop
				
			out_data <- c(h_bound,ancest_surface)
			for(rm in 1:length(cell_prob_list)) {
			  out_data <- c(out_data,unlist(cell_prob_list[rm]))
			}
			return(out_data)
			}#end pop_function

			pop_function(i,j)

#			}#else bracket
	}#dopar end bracket

#	cat("Almost done")
#	cat('\n')
	h_bound <- x[,seq(from=1,to=(ncols*(2+num_axes)),by=2+num_axes)]
	ancest_surface <- x[,seq(from=2,to=(ncols*(2+num_axes)),by=2+num_axes)]
	name_list <- list()
	for(fl in 1:num_axes) {
	  xx<-paste("axis",fl,sep='')
	  name_list[[xx]]<- 0
	}
	for(am in 1:length(name_list)) {
	  name_list[[am]] <- x[,seq(from=(2+am),to=(ncols*(2+num_axes)),by=2+num_axes)]
	}
	y <- list(h_bound,ancest_surface)
	y <- append(y, name_list)
	
	return(y)	
} #function end bracket
