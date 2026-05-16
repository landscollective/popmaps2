#' @title Visualize jackknife results
#'
#' @description This function takes results from jackknife() and uses linear models
#'    to calculate R-squared values for each genetic axis for each combination of the parameters num_tested and popmod.
#'    In other words, lm(known ancestry coefficients across sampling locations) ~ 
#'    (predicted ancestry coefficients across sampling locations for a parameter combination).
#'    Results of the linear models can be visualized per axis or averaged across axes (default).
#' @param input_locs An R object (rows = total # empirical sites, columns = total # genetic axes + 3) 
#'    with column 1: site name; column 2: decimal longitude; column 3: decimal latitude; 
#'    column 4…column x: ancestry coefficients for genetic axis 1…genetic axis x. 
#'    Function depends on this precise format – see example data hija_struc.
#' @param jackknife_data An R object resulting from executing function jackknife(). The matrix has the columns:
#'    num_tested parameter; popmod parameter; estimated ancestry coefficient axis 1...estimated ancestry
#'    coefficient axis x; empirical ancestry coefficient axis 1...empirical ancestry coefficient axis x. 
#' @param axis An integer indicating which genetic axis to visualize R2 values as a heatmap for 
#'    all parameter combinations. The default value averages R2 values across all axes for visualization.
#'@references Massatti R & Winkler DE. (2022) Spatially explicit management of genetic diversity using 
#'    ancestry probability surfaces. Methods in Ecology and Evolution. http://dx.doi.org/10.1111/2041-210X.13902
#' @author Rob Massatti
#' @examples
#'    ex_raster <- raster::aggregate(hija_raster,fact=16)   #Cells in embedded raster are aggregated to reduce computation time
#'    jack_data <- jackknife(input_raster=ex_raster,input_locs=hija_struc,surface="G")
#'    jackknife_viz(input_locs=hija_struc, jackknife_data=jack_data,axis=0)
#'
#' @export

jackknife_viz <- function(input_locs="", jackknife_data="",axis=0) {
	if (!requireNamespace("gplots", quietly = TRUE)) {
		stop("The 'gplots' package is required for jackknife_viz().", call. = FALSE)
	}
	################################################
	# Define outfiles & variables
	################################################
	data <- jackknife_data
	n_tested <- length(unique(data[,1]))
	exp <- length(unique(data[,2]))
	species_loc <- input_locs
	num_axes <- length(species_loc[1,])-3
	n_samp <- length(species_loc[,1])
	n_comb <- length(data[,1])/n_samp
	
	rsquare <- NULL 
	for(i in 1:n_comb) {
		start <- (i*n_samp-(n_samp-1))
		fin <- i*n_samp
		data_temp <- data[start:fin,]
		k_obs_list <- list()
		k_pred_list <- list()
		for(w in 1:num_axes) {
			xx<-paste("kobs_",w,sep='')
			k_obs_list[[xx]]<- 0	
			ww<-paste("kpred_",w,sep='')
			k_pred_list[[ww]]<- 0	
		}      
		for(q in 1:num_axes) {
		  for(rtm in 1:n_samp) {
		    if(data_temp[,(2+num_axes+q)][rtm] == 0) {
		      data_temp[,(2+num_axes+q)][rtm] <- 0.00001
		    } else if (data_temp[,(2+num_axes+q)][rtm] == 1) {
		      data_temp[,(2+num_axes+q)][rtm] <- .999999
		    }
		    if(data_temp[,(2+q)][rtm] == 0) {
		      data_temp[,(2+q)][rtm] <- 0.00001
		    } else if (data_temp[,(2+q)][rtm] == 1) {
		      data_temp[,(2+q)][rtm] <- .999999
		    }
		  }
		  k_obs_list[[q]] <- data_temp[,(2+num_axes+q)]
		  k_pred_list[[q]] <- data_temp[,(2+q)]
		}  
		
		for(p in 1:num_axes) {
		  data_vec <- NULL
		  temp_model <- lm(gtools::logit(unlist(k_pred_list[p]))~gtools::logit(unlist(k_obs_list[p])))
		  temp_r <- summary(temp_model)$r.squared		
		  data_vec <- c(data_temp[1,1],data_temp[1,2],temp_r)
		  rsquare <- rbind(rsquare,data_vec)
		}	
	}

	
	out_t_list <- list()
	for(e in 1:num_axes) {
		xx<-paste("out",e,"_t",sep='')
		out_t_list[[xx]]<- rsquare[seq(from=e,to=length(rsquare[,1]),by=num_axes),]
	}
	out_trun_list <- list()
	for(d in 1:num_axes) {
		xx<-paste("out",d,"_trun",sep='')
		out_trun_list[[xx]]<- 0		
	}
	for(i in 1:n_tested) {
		start <- (i*exp-(exp-1))
		fin <- i*exp
		out_t_temp <- list()
		for(d in 1:num_axes) {
			xx<-paste("out",d,"_t_temp",sep='')
			out_t_temp[[xx]]<- out_t_list[[d]][start:fin,]
		}
		for(e in 1:num_axes) {
			xx<-paste("out",e,"_trun",sep='')
			out_trun_list[[xx]]<- cbind(out_trun_list[[xx]],out_t_temp[[e]][,3])
		}	
	}
	for(q in 1:num_axes) {
		out_trun_list[[q]] <- out_trun_list[[q]][,2:(n_tested+1)]
	}        
	col_names <- (unique(data[,1]))
	row_names <- out_t_temp[[1]][,2]
	if(axis > 0) {
		r_temp <- out_trun_list[[axis]]
		suppressWarnings(gplots::heatmap.2(as.matrix(r_temp),Rowv=F,Colv=F,scale='none',trace='none',cellnote=as.matrix(round(r_temp,digits=3)),notecol='black',col=gplots::bluered(20),labCol=col_names,labRow=row_names))
	} else if(axis==0) {

	temp_data <- NULL
	for(k in 1:exp) {
		temp_vec <- NULL
		for(l in 1:n_tested) {	
			temp <- 0
			for(b in 1:num_axes) {
				temp <- temp + out_trun_list[[b]][k,l]
			}
			temp <- temp / num_axes
			temp_vec <- c(temp_vec, temp)
		}
		temp_data <- rbind(temp_data, temp_vec)
	}
	suppressWarnings(gplots::heatmap.2(as.matrix(temp_data), Rowv=F,Colv=F,scale='none',trace='none',cellnote=as.matrix(round(temp_data,digits=3)),notecol='black',col=gplots::bluered(20),labCol=col_names,labRow=row_names))
	}
}
