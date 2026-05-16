# Hilaria jamesii Fst table

library(usethis)

hija_fst <- read.table(file='hija_fst.txt',header=T,sep='\t')

row.names(hija_fst) <- hija_fst[,1]
hija_fst <- hija_fst[,-1]
hija_fst <- as.matrix(hija_fst)
hija_fst <- as.dist(m=hija_fst)

usethis::use_data(hija_fst, overwrite=TRUE)