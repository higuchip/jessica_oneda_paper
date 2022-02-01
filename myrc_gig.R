library(rgbif)

#GBIF
key.mg <- name_backbone(name='Myrcianthes gigantea')$speciesKey
key.mg
dat.mg <- occ_search(taxonKey=key.mg, return='data', limit=100)
gbifmap(dat.mg)
coord.gbif.mg<-cbind(dat.mg$decimalLongitude, dat.mg$decimalLatitude)
head(coord.gbif.mg)

plot(wrld_simpl, xlim=c(-55,-51), ylim=c(-35,-10), axes=TRUE, col="gray")
points(coord.gbif.mg, pch=19, cex=.51, col="black")
