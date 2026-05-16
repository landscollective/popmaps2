#'@title Pairwise Fst values between empirical sampling locations
#'
#'@description A data set consisting of Fst values between all pairs of Hilaria
#'    jamesii empirical sampling locations contained in hija_struc. 
#'
#'@format An object of class "dist" with 16 rows and 16 columns. Floats below the diagonal
#'    indicate pairwise Fst values. Row and column names indicate sampling site - 
#'    see hija_struc.
#'    
#'@references Massatti R & Knowles LL. (2020) The historical context of contemporary
#'    climatic adaptation: a case study in the climatically dynamic and environmentally complex
#'    southwestern United States. Ecography 43(5) 735-746.
#'
#'@seealso hija_struc, hija_raster
#'
#'@examples 
#'    location <- as.matrix(hija_struc[,2:3])
#'    geoDist <- raster::pointDistance(location,longlat=T)
#'    geoDist <- as.dist(geoDist)
#'    transition_layer <- gdistance::transition(hija_raster, transitionFunction=mean, directions=8)
#'    corrected_transC <- gdistance::geoCorrection(transition_layer,type="c",scl=T)
#'    costDist <- gdistance::costDistance(corrected_transC,location)
#'    cor(hija_fst,geoDist)
#'    cor(hija_fst,costDist) 
"hija_fst"
#'
#'@title Hilaria jamesii empirical genetic data
#'
#'@description A data set containing location information and ancestry coefficients
#'    for Hilaria jamesii sampling locations.
#'
#'@format A data frame with 16 rows and 6 columns. Columns include: 1 - sampling location
#'    name; 2 - decimal longitude; 3 - decimal latitude; 4 to the number of genetic 
#'    axes - ancestry coefficients per genetic axis.
#'    
#'@references Massatti R & Knowles LL. (2020) The historical context of contemporary
#'    climatic adaptation: a case study in the climatically dynamic and environmentally complex
#'    southwestern United States. Ecography 43(5) 735-746.
#'
#'@seealso hija_fst, hija_raster
#'
#'@examples 
#'    ex_raster <- raster::aggregate(hija_raster,fact=16)  #Cells in embedded raster are aggregated to reduce computation time
#'    pp <- popmaps(input_raster=ex_raster,input_locs=hija_struc,empirical_pt_dist=5,num_sites=15,num_tested=4,popmod=-0.05,threshold=0,surface='G')
#'    popmap_viz(pop_raster_list=pp,input_locs=hija_struc,input_raster=ex_raster,maptype='ancestry',crs="+init=epsg:5070",boundary_width=-0.01)
"hija_struc"
#'
#'@title Geospatial layer for ancestry probability surface estimation
#'
#'@description A raster layer created from a species distribution model (asci file) for 
#'    Hilaria jamesii (see Massatti & Winkler 2022). 
#'
#'@format A RasterLayer object with 480 rows and 540 columns. Values in the cells
#'    are logistic probabilities of occurrence calculated by Maxent. 
#'    
#'@references Massatti R & Knowles LL. (2020) The historical context of contemporary
#'    climatic adaptation: a case study in the climatically dynamic and environmentally complex
#'    southwestern United States. Ecography 43(5) 735-746.
#'    
#'    Massatti R & Winkler DE. (2022) Spatially explicit management of genetic diversity using 
#'    ancestry probability surfaces. Methods in Ecology and Evolution. http://dx.doi.org/10.1111/2041-210X.13902
#'
#'@seealso hija_struc
#'
#'@examples 
#'    plot(hija_raster)
#'    points(as.matrix(hija_struc[,2:3]))
"hija_raster"
#'
#'@title Locations of herbarium specimens for Hilaria jamesii.
#'
#'@description A data set consisting of 214 verified and georeferenced herbarium voucher specimens for Hilaria jamesii occurring 
#'    within the extent of hija_raster. These herbarium specimens were downloaded from the SEINet data portal (https://swbiodiversity.org/seinet/).
#'
#'@format A data frame consisting of 2 columns and 214 rows. Columns contain longitude and latitude in decimal degrees and 
#'    rows contain unique herbarium specimens.
#'    
#'@references Massatti R & Knowles LL. (2020) The historical context of contemporary
#'    climatic adaptation: a case study in the climatically dynamic and environmentally complex
#'    southwestern United States. Ecography 43(5) 735-746.
#'
#'@seealso hija_struc
#'
#'@examples 
#'    plot(hija_raster)
#'    points(hija_herb, pch=24, col= 'black', bg='red')
"hija_herb"

