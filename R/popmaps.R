#' @title Estimate an ancestry probability surface
#'
#' @description This function estimates an ancestry probability surface when supplied
#'     a geospatial layer and empirical genetic data describing patterns of ancestry coefficients
#'     across sampling locations. 
#' @param input_raster An R RasterLayer object defining the geographic extent for the spatial
#'     interpolation. Values in the cells will be used to calculate distance used in the
#'     dist_prob_func if surface = 'C'. See example data hija_raster.
#' @param input_locs An R object (rows = total # empirical sites, columns = total # genetic axes + 3]) 
#'     with column 1: site name; column 2: decimal longitude; column 3: decimal latitude; 
#'     column 4…column x: ancestry coefficients for genetic axis 1…genetic axis x. 
#'     Functions depend on this precise format – see example data hija_struc.
#' @param surface A string (either 'G' or 'C') that determines how input_raster is used to calculate
#'     ancestry coefficients. If 'G', the spatial attributes of input_raster will be used to calculate
#'     geographic (i.e., Euclidean) distances between empirical sites and inference sites. If 'C',
#'     least-cost distances are calculated using the values contained in the raster cells.
#' @param surface_values Character. Meaning of raster values when `surface = "C"`.
#'     `"suitability"` and `"conductance"` use values directly; `"resistance"`
#'     converts values to conductance using an inverse transform.
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
#' @param ncore An integer retained for compatibility with POPMAPS 1.03 calls.
#'     The optimized geographic- and least-cost-distance paths currently
#'     precompute distances and run sequentially.
#' @param dist_prob_func A function defining the relationship between distance and the contribution 
#'     of an empirical site’s ancestry coefficients to the estimation of ancestry coefficients at 
#'     an inference cell. The default equation defines the relationship in Fig. 2 of Massatti & Winkler (2022).
#' @param num_sites An integer representing the pool of empirical sites considered when selecting 
#'     num_tested sites to estimate ancestry coefficients. If empirical sites are highly clustered 
#'     and rarefaction due to empirical_pt_dist causes many to be discarded, this variable will likely 
#'     need to be increased at the cost of computing time.
#' @param rescale_conductance Logical. If `TRUE`, rescale non-missing
#'     conductance values to the range 0-1 after any resistance conversion.
#' @param resistance_epsilon Positive numeric scalar added to resistance values
#'     before inversion when `surface_values = "resistance"`.
#' @references Massatti R & Winkler DE. (2022) Spatially explicit management of genetic diversity using 
#'     ancestry probability surfaces. Methods in Ecology and Evolution. http://dx.doi.org/10.1111/2041-210X.13902
#' @author Rob Massatti
#' @examples
#' \dontrun{
#'     ex_raster <- raster::aggregate(hija_raster,fact=16)  #Cells in embedded raster are aggregated to reduce computation time
#'     pp <- popmaps(input_raster=ex_raster,input_locs=hija_struc,empirical_pt_dist=5,num_sites=15,num_tested=4,popmod=-0.05,threshold=0,surface='G')
#' }
#' @export


popmaps <- function(input_raster='',
                    input_locs='',
                    surface='G',
                    empirical_pt_dist=5,
                    num_sites=10,
                    num_tested=3,
                    popmod=-0.001,
                    ncore=4,
                    threshold=0,
                    dist_prob_func=function(popmod_temp,distance) {exp(popmod_temp*distance)},
                    surface_values = c("suitability", "conductance", "resistance"),
                    rescale_conductance = FALSE,
                    resistance_epsilon = sqrt(.Machine$double.eps)) {

  surface_values <- match.arg(surface_values)
  if (!is.logical(rescale_conductance) || length(rescale_conductance) != 1 || is.na(rescale_conductance)) {
    stop("`rescale_conductance` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  prepared <- popmaps_prepare_inputs(
    input_raster = input_raster,
    input_locs = input_locs,
    surface = surface,
    num_sites = num_sites,
    num_tested = num_tested,
    threshold = threshold,
    ncore = ncore,
    empirical_pt_dist = empirical_pt_dist,
    popmod = popmod,
    require_legacy_c = FALSE
  )
  surface <- prepared$surface
  if (!is.function(dist_prob_func)) {
    stop("`dist_prob_func` must be a function.", call. = FALSE)
  }

	raster_surface <- prepared$raster
	species_data <- prepared$locations

  if (surface == "G") {
    return(popmaps_geographic_surface(
      raster_surface = raster_surface,
      species_data = species_data,
      empirical_pt_dist = empirical_pt_dist,
      num_sites = num_sites,
      num_tested = num_tested,
      popmod = popmod,
      threshold = threshold,
      dist_prob_func = dist_prob_func
    ))
  }

  popmaps_cost_surface(
    raster_surface = raster_surface,
    species_data = species_data,
    empirical_pt_dist = empirical_pt_dist,
    num_sites = num_sites,
    num_tested = num_tested,
    popmod = popmod,
    threshold = threshold,
    dist_prob_func = dist_prob_func,
    surface_values = surface_values,
    rescale_conductance = rescale_conductance,
    resistance_epsilon = resistance_epsilon
  )
} #function end bracket
