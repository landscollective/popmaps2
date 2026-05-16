#' @title Estimate ancestry coefficents for empirical sites
#'
#' @description This function tests how parameter combinations (num_tested & popmod)
#'     influence the ability of the popmaps algorithm to estimate known values of
#'     ancestry coefficients using a leave-one-out approach.
#' @param input_raster An R RasterLayer object defining the geographic extent for the 
#'     spatial interpolation. Values in the cells will be used to calculate distance used 
#'     in the dist_prob_func if surface = ‘C’. See example data hija_raster.
#' @param input_locs An R object (rows = total # empirical sites, columns = total # genetic axes + 3) 
#'     with column 1: site name; column 2: decimal longitude; column 3: decimal latitude; 
#'     column 4…column x: ancestry coefficients for genetic axis 1…genetic axis x. 
#'     Function depends on this precise format – see example data hija_struc.
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
#' @param num_sites An integer representing the pool of empirical sites considered when selecting 
#'     num_tested sites to estimate ancestry coefficients. If empirical sites are highly clustered 
#'     and rarefaction due to empirical_pt_dist causes many to be discarded, this variable will likely 
#'     need to be increased at the cost of computing time.
#' @param num_tested_vec A vector of integers representing the values of num_tested that should be 
#'     tested during jackknifing. All pairwise combinations with popmod_vec will be tested.
#' @param popmod_vec A vector of floats representing the values of popmod that should be tested 
#'     during jackknifing. All pairwise combinations with num_tested_vec will be tested.
#' @param dist_prob_func A function defining the relationship between distance and the contribution 
#'     of an empirical site’s ancestry coefficients to the estimation of ancestry coefficients at 
#'     an inference cell. The default equation defines the relationship in Fig. 2 of Massatti & Winkler (2022).
#' @references Massatti R & Winkler DE. (2022) Spatially explicit management of genetic diversity using 
#'     ancestry probability surfaces. Methods in Ecology and Evolution. http://dx.doi.org/10.1111/2041-210X.13902
#' @author Rob Massatti
#' @examples 
#'     ex_raster <- raster::aggregate(hija_raster,fact=16)  #Cells in embedded raster are aggregated to reduce computation time
#'     jack_data <- jackknife(input_raster=ex_raster,input_locs=hija_struc,surface="G")
#'
#' @export

jackknife <- function(input_raster="", input_locs="", surface='G', empirical_pt_dist=5, num_sites=10, num_tested_vec=c(2,3,4,5,6,7,8), popmod_vec=c(-0.001,-0.01,-0.05,-0.1,-0.15),dist_prob_func=function(popmod_temp,distance) {exp(popmod_temp*distance)}) {

  if (!is.numeric(num_tested_vec) || length(num_tested_vec) < 1 ||
      any(!is.finite(num_tested_vec)) || any(num_tested_vec < 1) ||
      any(num_tested_vec != floor(num_tested_vec))) {
    stop("`num_tested_vec` must contain positive whole numbers.", call. = FALSE)
  }
  if (!is.numeric(popmod_vec) || length(popmod_vec) < 1 || any(!is.finite(popmod_vec))) {
    stop("`popmod_vec` must contain finite numeric values.", call. = FALSE)
  }
  if (!is.function(dist_prob_func)) {
    stop("`dist_prob_func` must be a function.", call. = FALSE)
  }

  prepared <- popmaps_prepare_inputs(
    input_raster = input_raster,
    input_locs = input_locs,
    surface = surface,
    num_sites = num_sites,
    num_tested = max(num_tested_vec),
    empirical_pt_dist = empirical_pt_dist,
    jackknife = TRUE
  )
  surface <- prepared$surface

  raster_surface <- prepared$raster
	cell_size <- raster::res(raster_surface)[1]
	species_data <- prepared$locations
	species_pts <- sp::SpatialPointsDataFrame(species_data[,2:3], species_data)
	sampling_loc_coords <- cbind(species_pts@data$V2,species_pts@data$V3)
	num_emp_sites <- length(sampling_loc_coords[,1])
	num_axes <- length(species_data[1,])-3

	###############################################
	#####Modifies species_surface according to surface setting
	###############################################
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

	###############################################
	#####For each of the empirical sampling sites, determines the raster cell (in terms of cell location in the raster)
	###############################################
	ij_mat <- matrix(nrow=num_emp_sites,ncol=2)
	for(d in 1:length(sampling_loc_coords[,1])){
		i <- (ymax-sampling_loc_coords[d,2])/cell_size
		j <- (sampling_loc_coords[d,1]-xmin)/cell_size
		ij_mat[d,1] <- round(i)
		ij_mat[d,2] <- round(j)
	}

	###############################################
	#####Function to calculate geographic distance between locations
	###############################################
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

	###############################################
	#####Main body of function
	###############################################

	jack_data <- NULL
	cat("Completed: ")
	#num_tested_vec iterator
	for(m in 1:length(num_tested_vec)) {
	  num_tested_temp <- num_tested_vec[m]

		#popmod_vec iterator
		for(n in 1:length(popmod_vec)) {
			popmod_temp <- popmod_vec[n]

			#emp. site iterator
			for(p in 1:length(ij_mat[,1])) {
				i<-ij_mat[p,1]
				j<-ij_mat[p,2]
				y<- ymax+(cell_size/2)-(i*cell_size)
				x<- xmin+(cell_size/2)+(j*cell_size)
				point <- cbind(x,y)
				distance_vec<-NULL

				for(a in 1:num_emp_sites) {
					ed <- earth.dist(point[1,2],point[1,1],sampling_loc_coords[a,2],sampling_loc_coords[a,1])
					distance_vec <- c(distance_vec,ed)
				}

				#drop current point from distance_vec
				distance_vec <- distance_vec[-p]

				temp_site_data <- species_pts@data
				temp_site_data <- temp_site_data[-p,]
				temp_site_data <- cbind(temp_site_data,distance_vec)
				dist_vec_position <- length(temp_site_data[1,])
				temp_site_data <- temp_site_data[order(temp_site_data[,dist_vec_position]),]
				temp_site_data <- temp_site_data[1:num_sites,]
				#temp_site_data <- temp_site_data[temp_site_data[,(4+num_axes)]>empirical_pt_dist,]

				if(surface=="G") {
				    temp_site_data <- cbind(temp_site_data,seq(1:length(temp_site_data[,1])))
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

						while(length(empirical_site_list[,1]) < num_tested_temp) {
							match_temp <- as.numeric(temp_site_data[iterator,2:3])

							keep <- T
							for(k in 1:length(empirical_site_list[,1])) {
								match <- empirical_site_list[k,]
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
						axis_list <- list()
						cell_prob_list <- list()
						temp_list <- list()
						for(w in 1:num_axes) {
							yy<-paste("axis",w,sep='')
							axis_list[[yy]]<- 0

							zz<-paste("cell_prob",w,sep='')
							cell_prob_list[[zz]]<- 0

							aa<-paste("temp_list",w,sep='')
							temp_list[[aa]]<- 0
						}
						site_weight <- NULL
						for(b in 1:num_tested_temp) {
							match <- empirical_site_list[b,]

							temp_list_iter <- 4
							for(c in 1:num_axes) {
								temp_list[[c]] <- temp_site_data[,temp_list_iter][as.numeric(match[1])==temp_site_data$V2]
								temp_list_iter <- temp_list_iter + 1
							}
							for(f in 1:length(temp_list)) {
								axis_list[[f]] <- c(axis_list[[f]], temp_list[[f]])
							}
							distance <- temp_site_data[,dist_vec_position][as.numeric(match[1])==temp_site_data$V2]
							dist_prob <- dist_prob_func(popmod_temp, distance)
							if(dist_prob > 1) {
							  dist_prob <- 1
							} else if (dist_prob < 0) {
							  dist_prob <- 0 
							}
							site_weight <- c(site_weight, dist_prob)
						}
						for(q in 1:length(axis_list)) {
							axis_list[[q]] <- axis_list[[q]][2:(num_tested_temp+1)]
						}
						for(r in 1:length(cell_prob_list)) {
							cell_prob_list[[r]] <- sum(axis_list[[r]]*site_weight/num_tested_temp)
						}

						data_vec <- NULL
						data_vec <- c(as.numeric(num_tested_temp),as.numeric(popmod_temp))

						for(u in 1:num_axes) {
							data_vec <- c(data_vec,as.numeric(cell_prob_list[[u]]/sum(unlist(cell_prob_list))))
						}
						for(v in 1:num_axes) {
							data_vec <- c(data_vec,unlist(species_pts@data[p,][v+3]))
						}
				}

				if(surface=="C") {
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

						while(length(empirical_site_list[,1]) < num_tested_temp) {
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
						axis_list <- list()
						cell_prob_list <- list()
						temp_list <- list()
						for(w in 1:num_axes) {
							yy<-paste("axis",w,sep='')
							axis_list[[yy]]<- 0

							zz<-paste("cell_prob",w,sep='')
							cell_prob_list[[zz]]<- 0

							aa<-paste("temp_list",w,sep='')
							temp_list[[aa]]<- 0
						}
						site_weight <- NULL
						for(h in 1:num_tested_temp) {
							match <- empirical_site_list[h,]
							temp_list_iter <- 4
							for(o in 1:num_axes) {
								temp_list[[o]] <- temp_site_data[,temp_list_iter][as.numeric(match[1])==temp_site_data$V2]
								temp_list_iter <- temp_list_iter + 1
							}
							for(s in 1:length(temp_list)) {
								axis_list[[s]] <- c(axis_list[[s]], temp_list[[s]])
							}
							distance <- temp_site_data[,(dist_vec_position+1)][as.numeric(match[1])==temp_site_data$V2]
							dist_prob <- dist_prob_func(popmod_temp,distance)
							if(dist_prob > 1) {
							  dist_prob <- 1
							} else if (dist_prob < 0) {
							  dist_prob <- 0 
							}
							site_weight <- c(site_weight, dist_prob)
						}
						for(t in 1:length(axis_list)) {
							axis_list[[t]] <- axis_list[[t]][2:(num_tested_temp+1)]
						}
						for(dd in 1:length(cell_prob_list)) {
							cell_prob_list[[dd]] <- sum(axis_list[[dd]]*site_weight/num_tested_temp)
						}

						data_vec <- NULL
						data_vec <- c(as.numeric(num_tested_temp),as.numeric(popmod_temp))

						for(bb in 1:num_axes) {
							data_vec <- c(data_vec,as.numeric(cell_prob_list[[bb]]/sum(unlist(cell_prob_list))))
						}
						for(cc in 1:num_axes) {
							data_vec <- c(data_vec,unlist(species_pts@data[p,][cc+3]))
						}
				}
			jack_data <- rbind(jack_data,data_vec)
			}
		}
	  cat(floor(m/length(num_tested_vec)*100))
	  cat('%')
	  cat(' ')
	}
return(jack_data)
}
