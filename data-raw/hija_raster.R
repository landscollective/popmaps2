# Hilaria jamesii ascii

library(usethis)
library(raster)

hija_raster <- raster::raster('./hija_raster.asc')
hija_raster <- readAll(hija_raster)

usethis::use_data(hija_raster, overwrite=T)
