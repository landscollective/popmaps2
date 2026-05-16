# Hilaria jamesii structure and location data

library(usethis)

hija_struc <- read.csv(file='hija_struc.txt',header=F, sep='\t',stringsAsFactors=F)

usethis::use_data(hija_struc, overwrite=TRUE)