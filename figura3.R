

tiff(filename = "plot3_new_new.tiff",
     width = 1750, height = 3800, units = "px",
     bg = "white", 
     res = 300)
par(mfrow=c(7,3),mar=c(4.5,4.5,1,1), oma=c(4.5,3,1,1)) 


#Group 1
# gym_klo

plot(myBiomodRangeSize.drimys.2070.85.gym.klo$Diff.By.Pixel, col= c("darkred","seagreen", "gray", "white"), 
     legend=FALSE, xlim=c(-60,-35))

plot(wrld_simpl,  axes="TRUE", add=T, xlim=c(-60,-35))
text(-40,-25, "-45.9 %", cex=.8)
text(-46,-35.5,"G. klotzschiana", font=3, cex=.8)
#polygon(poligono_mata_atlantica)



# all.edu

plot(myBiomodRangeSize.drimys.2070.85.all.edu$Diff.By.Pixel, col= c("darkred","seagreen", "gray", "white"), 
     legend=FALSE)

plot(wrld_simpl,  axes="TRUE", add=T)
text(-40,-25, "-49.6 %", cex=.8)
text(-46,-35.5,"A. edulis", font=3, cex=.8)
#polygon(poligono_mata_atlantica)


# vit.meg


plot(myBiomodRangeSize.drimys.2070.85.vit.meg$Diff.By.Pixel, col= c("darkred","seagreen", "gray", "white"), 
     legend=FALSE)

plot(wrld_simpl,  axes="TRUE", add=T)
text(-40,-25, "-57.5 %", cex=.8)
text(-46,-35.5,"V. megapotamica", font=3, cex=.8)
#polygon(poligono_mata_atlantica)




# cam.xan


plot(myBiomodRangeSize.drimys.2070.85.cam.xan.new$Diff.By.Pixel, col= c("darkred","seagreen", "gray", "white"), 
     legend=FALSE)

plot(wrld_simpl,  axes="TRUE", add=T)
text(-40,-25, "-52.4 %", cex=.8)
text(-46,-35.5,"C. xanthocarpa", font=3, cex=.8)
# polygon(poligono_mata_atlantica)


# eug.uni

plot(myBiomodRangeSize.drimys.2070.85.eug.uni$Diff.By.Pixel, col= c("darkred","seagreen", "gray", "white"), 
     legend=FALSE)

plot(wrld_simpl,  axes="TRUE", add=T)
text(-40,-25, "-34.8 %", cex=.8)
text(-46,-35.5,"E. uniflora", font=3,  cex=.8)
#polygon(poligono_mata_atlantica)


#Group 2

# and.fra

plot(myBiomodRangeSize.drimys.2070.85.and.fra$Diff.By.Pixel, col= c("darkred","seagreen", "gray", "white"), 
     legend=FALSE)

plot(wrld_simpl,  axes="TRUE", add=T)
text(-40,-25, "85.2 %", cex=.8)
text(-46,-35.5,"A. fraxinifolia", font=3,  cex=.8)
#polygon(poligono_mata_atlantica)

# per.gla

plot(myBiomodRangeSize.2070.85.per.gla$Diff.By.Pixel, col= c("darkred","seagreen", "gray", "white"), 
     legend=FALSE)

plot(wrld_simpl,  axes="TRUE", add=T)
text(-40,-25, "-76.3 %", cex=.8)
text(-46,-35.5,"P. glabrata", font=3,  cex=.8)
# polygon(poligono_mata_atlantica)



# ani.fir

plot(myBiomodRangeSize.2070.85.ani.fir$Diff.By.Pixel,col= c("darkred","seagreen", "gray", "white"), 
     legend=FALSE)

plot(wrld_simpl,  axes="TRUE", add=T)
text(-40,-25, "-60.3  %", cex=.8)
text(-46,-35.5,"A. firmula", font=3, cex=.8)
# polygon(poligono_mata_atlantica)

# cec.gla

plot(myBiomodRangeSize.2070.85.cec.gla$Diff.By.Pixel, col= c("darkred","seagreen", "gray", "white"), 
     legend=FALSE)

plot(wrld_simpl,  axes="TRUE", add=T)
text(-40,-25, "-86.6  %", cex=.8)
text(-46,-35.5,"C. glaziovii", font=3,  cex=.8)
# polygon(poligono_mata_atlantica)


# hye.alc

plot(myBiomodRangeSize.2070.85.hye.alc$Diff.By.Pixel, col= c("darkred","seagreen", "gray", "white"), 
     legend=FALSE)

plot(wrld_simpl,  axes="TRUE", add=T)
text(-40,-25, "-33.1  %", cex=.8)
text(-46,-35.5,"H. alchorneoides", font=3,  cex=.8)
# polygon(poligono_mata_atlantica)


# nec.opp

plot(myBiomodRangeSize.2070.85.nec.opp$Diff.By.Pixel, col= c("darkred","seagreen", "gray", "white"), 
     legend=FALSE)

plot(wrld_simpl,  axes="TRUE", add=T)
text(-40,-25, "-58.4  %", cex=.8)
text(-46,-35.5,"N. oppositifolia", font=3,  cex=.8)
# polygon(poligono_mata_atlantica)


# gua.aus

plot(myBiomodRangeSize.2070.85.gua.aus$Diff.By.Pixel, col= c("darkred","seagreen", "gray", "white"), 
     legend=FALSE)

plot(wrld_simpl,  axes="TRUE", add=T)
text(-40,-25, "-83.39   %", cex=.8)
text(-46,-35.5,"G. australis", font=3,  cex=.8)
# polygon(poligono_mata_atlantica)



# lec.pis
plot(myBiomodRangeSize.2070.85.lec.pis$Diff.By.Pixel, col= c("darkred","seagreen", "gray", "white"), 
     legend=FALSE)

plot(wrld_simpl,  axes="TRUE", add=T)
text(-40,-25, "-15.3  %", cex=.8)
text(-46,-35.5,"L. pisonis", font=3, cex=.8)
# polygon(poligono_mata_atlantica)

# myr.rac
plot(myBiomodRangeSize.2070.85.myr.rac$Diff.By.Pixel, col= c("darkred","seagreen", "gray", "white"), 
     legend=FALSE)

plot(wrld_simpl,  axes="TRUE", add=T)
text(-40,-25, "-2.1 %", cex=.8)
text(-46,-35.5,"M. racemosa", font=3, cex=.8)
# polygon(poligono_mata_atlantica)

# pse.con
plot(myBiomodRangeSize.2070.85.pse.con$Diff.By.Pixel, col= c("darkred","seagreen", "gray", "white"), 
     legend=FALSE)

plot(wrld_simpl,  axes="TRUE", add=T)
text(-40,-25, "-78.2  %", cex=.8)
text(-46,-35.5,"P. contorta", font=3,  cex=.8)
# polygon(poligono_mata_atlantica)

# ann.dol

plot(myBiomodRangeSize.2070.85.ann.dol$Diff.By.Pixel, col= c("darkred","seagreen", "gray", "white"), 
     legend=FALSE)

plot(wrld_simpl,  axes="TRUE", add=T)
text(-40,-25, " -70.0 %", cex=.8)
text(-46,-35.5, "A. dolabripetala", font=3,  cex=.8)
# polygon(poligono_mata_atlantica)


# cou.mic
plot(myBiomodRangeSize.2070.85.cou.mic$Diff.By.Pixel, col= c("darkred","seagreen", "gray", "white"), 
     legend=FALSE)

plot(wrld_simpl,  axes="TRUE", add=T)
text(-40,-25, " -4.1 %", cex=.8)
text(-46,-35.5, "C. microcarpa", font=3,  cex=.8)
# polygon(poligono_mata_atlantica)

# eut.edu

plot(myBiomodRangeSize.2070.85.eut.edu$Diff.By.Pixel, col= c("darkred","seagreen", "gray", "white"), 
     legend=FALSE)

plot(wrld_simpl,  axes="TRUE", add=T)
text(-40,-25, "-77.4 %", cex=.8)
text(-46,-35.5, "E. edulis", font=3, cex=.8)
# polygon(poligono_mata_atlantica)

# gar.gar

plot(myBiomodRangeSize.2070.85.gar.gar$Diff.By.Pixel, col= c("darkred","seagreen", "gray", "white"), 
     legend=FALSE)

plot(wrld_simpl,  axes="TRUE", add=T)
text(-40,-25, "-35.1 %", cex=.8)
text(-46,-35.5, "G. gardneriana", font=3,  cex=.8)
# polygon(poligono_mata_atlantica)


# pse.gra


plot(myBiomodRangeSize.2070.85.pse.gra$Diff.By.Pixel, col= c("darkred","seagreen", "gray", "white"), 
     legend=FALSE)

plot(wrld_simpl,  axes="TRUE", add=T)
text(-40,-25, "-67.1 %", cex=.8)
text(-46,-35.5, "P. grandiflorum", font=3,  cex=.8)
# polygon(poligono_mata_atlantica)


# alc.tri

plot(myBiomodRangeSize.drimys.2070.85.alc.tri$Diff.By.Pixel,col= c("darkred","seagreen", "gray", "white"), 
     legend=FALSE)

plot(wrld_simpl,  axes="TRUE", add=T)
text(-40,-25, "-59.9 %", cex=.8)
text(-46,-35.5, "A. triplinervia", font=3, cex=.8)
# polygon(poligono_mata_atlantica)

dev.off()

