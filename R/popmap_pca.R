#' @title Calculate an environmental PCA
#'
#' @description This function performs a principal components analysis on environmental variables across the geographic 
#'     extent of a popmaps() inference. The function crops environmental layers to the extent of the raster layer used in
#'     other popmaps analyses (i.e., input_raster).
#' @param input_raster An R RasterLayer object defining the geographic extent for the spatial 
#'     interpolation.  
#' @param bio_dir A string defining the pathway to a folder containing the environmental variables 
#'     that will be used in a PCA. Environmental layers must be able to be loaded by raster::raster(). For example,
#'     the 19 biolclimatically informative variable available form WorldClim (https://www.worldclim.org/) work well. Layers must
#'     be the same size or larger than input_raster - if larger, popmap_pca() will crop layers to the dimensions of input_raster.
#' @references Massatti R & Winkler DE. (2022) Spatially explicit management of genetic diversity using 
#'     ancestry probability surfaces. Methods in Ecology and Evolution. http://dx.doi.org/10.1111/2041-210X.13902
#' @author Rob Massatti
#' @examples
#' \dontrun{
#'     #This function requires a pathway to a folder containing environmental data layers that are not included in the POPMAPS package.
#'     pca <- popmap_pca(input_raster=hija_raster, bio_dir='./wc2.1_30s_bio/')
#'     plot(pca$Comp.1, ext = hija_raster@extent) 
#'     points(hija_struc[,2:3],pch=19, col='black')
#' }
#' @export

popmap_pca <- function(input_raster='', bio_dir='') {
  raster_surface <- popmaps_prepare_raster(input_raster)$raster
  
  file_list <- list.files(path=bio_dir,recursive=F)
  num_files <- length(list.files(path=bio_dir,recursive=F))
  predictors <- NULL
  for(i in 1:num_files) {
    temp <- raster::crop(raster::raster(paste(bio_dir,file_list[i],sep='')),raster_surface@extent)
    predictors <- c(predictors, temp)
  }
  predictors <- raster::stack(predictors)
  pca <- stats::princomp(stats::na.omit(predictors[]),cor=T)
#  summary(pca)
#  pca$loadings
  pcaRast <- raster::predict(predictors,pca,index=1:num_files)
  
  return(pcaRast)
}
