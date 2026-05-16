# Hilaria jamesii herbarium points

library(usethis)

hija_herb <- read.csv(file='./hija_herb.txt',header=T, sep='\t',stringsAsFactors=F)
hija_herb <- hija_herb[,2:3]

usethis::use_data(hija_herb, overwrite=TRUE)
