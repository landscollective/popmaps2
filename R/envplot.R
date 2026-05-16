#' @title Visualizing the environmental space of genetically defined populations
#'
#' @description This function visualizes the environmental variability across a defined geographic area and within the estimated 
#'     genetically defined populations of a species. After estimating an ancestry probability surface using popmaps(), a user would 
#'     generate random background points for each genetically defined population using bg_pop_pts(). In addition, points representing the focal
#'     species (i.e., the 'focal points,' such as may be downloaded from a herbarium or biodiversity database) are assigned to the estimated 
#'     populations using ptsNpop(). Finally, a principal components analysis (PCA) performed using popmap_pca() summarizes environmental 
#'     variation across user-supplied environmental layers. After extracting PC data for the background and focal points, data can be 
#'     visualized in PC space or as box plots. Data may be visualized for any pair of PC axes resulting from the PCA.
#' @param input_raster An R RasterLayer object defining the geographic extent for the 
#'     spatial interpolation. 
#' @param bg_env An R object resulting from extracting PC data for spatial points defined by bg_pop_pts(). After PC data are extracted by
#'     genetically defined population (see example below), data per axis are combined into a list, with each element representing a matrix of PC data 
#'     for the background points falling within each genetically defined population. 
#' @param pop_env An R object resulting from extracting PC data for spatial points defined by ptsNpop(). After PC data are extracted by
#'     genetically defined population (see example below), data per axis are combined into a list, with each element representing a matrix of PC data 
#'     for the points falling within each genetically defined population. 
#' @param pt_env This parameter is not currently implemented.
#' @param axis1 This will be the x-axis, or environmental variable 1, in the output from envplot(). Any number from 1 to the number of PC axes (i.e.,
#'     the number of environmental layers used in the PCA) may be used, though the first few PC axes generally explain the majority of environmental
#'     variation across the study area.
#' @param axis2 This will be the y-axis, or environmental variable 2, in the output from envplot().Any number from 1 to the number of PC axes (i.e.,
#'     the number of environmental layers used in the PCA) may be used, though the first few PC axes generally explain the majority of environmental
#'     variation across the study area.
#' @param plot_type Determines the graphical output that results from executing envplot(). If 'env', the output will be a graph displaying PC space
#'     for axis1 x axis2. Points (both background and focal) are colored uniquely by their genetically defined population. Background points are
#'     smaller and partially transparent. Focal points are larger and have contour lines that represent the density of points in environmental space. 
#'     If 'box', the output will contain multiple boxplots showing the range of environmental variation for the focal points according to the 
#'     genetically defined populations, as well as axis1 and axis2. 
#' @references Massatti R & Winkler DE. (2022) Spatially explicit management of genetic diversity using 
#'     ancestry probability surfaces. Methods in Ecology and Evolution. http://dx.doi.org/10.1111/2041-210X.13902
#' @author Rob Massatti
#' @examples
#'     ex_raster <- raster::aggregate(hija_raster,fact=16)  #Cells in embedded raster are aggregated to reduce computation time
#'     pp <- popmaps(input_raster=ex_raster,input_locs=hija_struc,empirical_pt_dist=5,num_sites=15,num_tested=4,popmod=-0.05,threshold=0,surface='G')
#'     bg_pts <- bg_pop_pts(pop_raster_list = pp, input_locs = hija_struc, input_raster = hija_raster, bg_pts=2000)
#'     samp_per_pop <- ptsNpop(pop_raster_list=pp, input_locs=hija_struc, input_raster=hija_raster,sampling_pts=hija_herb)
#'     
#'     #This function requires a pathway to a folder containing environmental data layers that are not included in the POPMAPS package.
#'     pca <- popmap_pca(input_raster=hija_raster, bio_dir='./wc2.1_30s_bio/')
#'     
#'     xx <- raster::extract(pca,bg_pts[[1]])
#'     xx <- xx[!is.na(xx[,1]),]
#'     yy <- raster::extract(pca,bg_pts[[2]])
#'     yy <- yy[!is.na(yy[,1]),]
#'     zz <- raster::extract(pca,bg_pts[[3]])
#'     zz <- zz[!is.na(zz[,1]),]
#'     bg_env <- list(xx,yy,zz)
#'     
#'     pop1_pca <- raster::extract(pca,samp_per_pop[[1]])
#'     pop2_pca <- raster::extract(pca, samp_per_pop[[2]])
#'     pop3_pca <- raster::extract(pca, samp_per_pop[[3]])
#'     pop_env <- list(pop1_pca,pop2_pca,pop3_pca)
#'     
#'     envplot(bg_env=bg_env,pop_env=pop_env,input_raster=hija_raster,plot_type='env')
#'     
#' @export

envplot <- function(bg_env='', pop_env='', pt_env=NULL, input_raster='', axis1=1, axis2=2, plot_type=c('env','box')) {

  plot_type <- match.arg(plot_type)
  if (!requireNamespace("viridis", quietly = TRUE)) {
    stop("The 'viridis' package is required for envplot().", call. = FALSE)
  }

  if(is.character(input_raster) == F) {
    raster_surface <- input_raster
  } else {
    raster_surface <- raster::raster(input_raster)
  }
  
  t_col <- function(color, percent = 50, name = NULL) {
    rgb.val <- col2rgb(color)
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                 maxColorValue = 255,
                 alpha = (100 - percent) * 255 / 100,
                 names = name)
    invisible(t.col)
  }

  num_axes <- length(bg_env)
  colors_axes <- viridis::viridis(num_axes,begin=0,end=1,direction=-1)
  col_axis <- colors_axes

  xmin <- NA
  xmax <- NA
  ymin <- NA
  ymax <- NA
  for(k in 1:num_axes) {
    xmin_temp <- min(bg_env[[k]][,axis1])
    xmax_temp <- max(bg_env[[k]][,axis1])
    ymin_temp <- min(bg_env[[k]][,axis2])
    ymax_temp <- max(bg_env[[k]][,axis2])
    
    if(is.na(xmin) == TRUE) {
      xmin <- xmin_temp
    } else if(xmin_temp < xmin) {
      xmin <- xmin_temp
    }
    
    if(is.na(xmax) == TRUE) {
      xmax <- xmax_temp
    } else if(xmax_temp > xmax) {
      xmax <- xmax_temp
    }
    
    if(is.na(ymin) == TRUE) {
      ymin <- ymin_temp
    } else if(ymin_temp < ymin) {
      ymin <- ymin_temp
    }
    
    if(is.na(ymax) == TRUE) {
      ymax <- ymax_temp
    } else if(ymax_temp > ymax) {
      ymax <- ymax_temp
    }
  }

  par(mfrow=c(1,1))
  par(oma=c(3,3,1,1))
  par(mar=c(1,2,1,2))
  
  if(plot_type == 'env') {
    
    plot(bg_env[[1]][,axis1],bg_env[[1]][,axis2],col=t_col(color=col_axis[1],percent=30), pch=16,cex=0.5,xlim=c(round(xmin)-1, round(xmax)+1), ylim=c(round(ymin)-1, round(ymax)+1))
    for(i in 2:num_axes) {
      points(bg_env[[i]][,axis1],bg_env[[i]][,axis2],col=t_col(color=col_axis[i],percent=70), pch=16,cex=0.5)
    }
    for(j in 1:num_axes) {
    
      zzz <- MASS::kde2d(pop_env[[j]][,axis1],pop_env[[j]][,axis2],n=50,lims=c(round(xmin)-2, round(xmax)+2, round(ymin)-2, round(ymax)+2))
      graphics::contour(zzz,lwd=1, drawlabels=F, add=T, col=colors_axes[j])
      points(pop_env[[j]][,axis1],pop_env[[j]][,axis2],bg=colors_axes[j],pch=21)
    }

    mtext(side=1, text='Env. axis 1',line=2.5, cex=1)
    mtext(side=2, text='Env. axis 2', line=2.8, cex=1)
  } else if(plot_type == 'box') {
  
    color_vec <- c(colors_axes, colors_axes)
    at_vec <- c(seq(1:num_axes),seq(from=(num_axes+2),to=(num_axes*2+1)))
    name_vec <- NULL
    for(t in 1:num_axes) {
      name_vec <- c(name_vec, paste('pop',t,sep=''))
    }
    name_vec <- c(name_vec,name_vec)
  
    f_list <- list()
    for(f in 1:num_axes) {
      f_list[[f]] <- pop_env[[f]][,1]
    }
    for(g in 1:num_axes) {
      f_list[[g+num_axes]] <- pop_env[[g]][,2]
    }
    
    boxplot(f_list,
          las=1, 
          names=name_vec, 
          notch=F, 
          boxwex=0.9, 
          at = at_vec, 
          col=color_vec, 
          border=T)
    mtext(side=1, text='Env. axis 1',line=2.5, cex=1, at=num_axes/2+0.5)
    mtext(side=1, text='Env. axis 2',line=2.5, cex=1, at=(num_axes/2+(num_axes+2)-0.5))
  }
}
