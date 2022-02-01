####ADEQUACAO
install.packages("ff")
library(biomod2)
library(splancs)
library(SDMTools)
library(dplyr)
library(tidyr)
library(gridExtra) 
library(raster)
library(spThin)
library(BIEN)

BioC<- getData('worldclim', var='bio', res=10) # Bioclima
variaveis <- stack(BioC)

plot(wrld_simpl, xlim=c(-130,-35), ylim=c(-45,35), axes=TRUE, col="gray")
plot(wrld_simpl, xlim=c(-70,-45), ylim=c(-35,-5), axes=TRUE, col="gray")
points(coord.ing.thi, pch=19, cex=.51, col="black")


##Poligono MA

require(rgdal)
data.shape<-readOGR(dsn="/media/dendrologia/423EF1863EF172F1/User/Documents/orientados/Guilherme/IDS_09_Desm_Bioma_Mata_Atlantica_1_2", 
                    layer="IDS_09_Desm_Bioma_Mata_Atlantica_1_2")

poligono_mata_atlantica<-data.shape@polygons[[1]]@Polygons[[1]]@coords

polygon(poligono_mata_atlantica)
###Poligono especie


poligono_ing.thi<-getpoly()

colnames(poligono_ing.thi) = c("Long", "Lat")
polygon(poligono_ing.thi) 



########REMOCAO DUVIDAS

out.ing.thi= pnt.in.poly(coord.ing.thi,poligono_ing.thi)
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
plot(data.shape, add=T, col="lightgreen")
points(out.ing.thi[which(out.ing.thi$pip==1),1:2],pch=20, cex=.5, col="black")
coord.ing.thi.limpo<-out.ing.thi[which(out.ing.thi$pip==1),1:2]
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
points(coord.ing.thi.limpo,pch=20, cex=.5, col="black")


###FILTRAGEM ESPACIAl

head(coord.ing.thi.limpo)
coord.ing.thi.limpo<-na.omit(coord.ing.thi.limpo)

dim(coord.ing.thi.limpo)
ing.thi_names<-c(rep("ing.thi",  dim(coord.ing.thi.limpo)[1]))
coord.ing.thi.limpo.names<-cbind(ing.thi_names, coord.ing.thi.limpo)

head(coord.ing.thi.limpo.names)
is.data.frame(coord.ing.thi.limpo.names)
colnames(coord.ing.thi.limpo.names)<-c("SPEC", "LONG", "LAT")
# 
# ?write.table
# write.table(coord.ing.thi.limpo.names, "tap_gui.csv", sep=";", dec=",")

coord.ing.thi.names_thinned_dataset_full <-
  thin( loc.data =coord.ing.thi.limpo.names, 
        lat.col = "LAT", 
        long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 26.27, reps =1, #10 min resolution
        locs.thinned.list.return = TRUE, 
        write.files = T, 
        max.files = 5, 
        out.dir = "coord.ing.thi.limpo.names_full/", out.base = "coord.ing.thi.limpo.names_thinned", 
        write.log.file = TRUE,
        log.file = "coord.ing.thi.limpo.names_thinned_full_log_file.txt" )



dim(coord.ing.thi.limpo.names)
dim(coord.ing.thi.names_thinned_dataset_full[[1]])



###ENVIRONMENTAL BACKGROUND

reg.variaveis.ing.thi= crop(variaveis, poligono_ing.thi)
plot(reg.variaveis.ing.thi, 1)
points(coord.ing.thi.limpo, pch=20, cex=.5)

#remoção variaveis nao colineares
vari.n.colineares.ing.thi<-vifstep(reg.variaveis.ing.thi)
vari.n.colineares.ing.thi


reg.predictors.ing.thi<- stack(reg.variaveis.ing.thi[[2]],
                               reg.variaveis.ing.thi[[3]],
                               reg.variaveis.ing.thi[[4]],
                               reg.variaveis.ing.thi[[8]],
                               reg.variaveis.ing.thi[[13]],
                               reg.variaveis.ing.thi[[14]],
                               reg.variaveis.ing.thi[[15]],
                               reg.variaveis.ing.thi[[18]],
                               reg.variaveis.ing.thi[[19]])



names(reg.predictors.ing.thi)


###MATA ATLANTICA projection presente

ma.presente = crop(variaveis, poligono_mata_atlantica)
plot(ma.presente, 1)
polygon(poligono_mata_atlantica)
ma.presente.poly <- mask(ma.presente, data.shape)
plot(ma.presente.poly, 1 )


ma.presente.poly.ing.thi<-stack(ma.presente.poly[[2]],
                                ma.presente.poly[[3]],
                                ma.presente.poly[[4]],
                                ma.presente.poly[[8]],
                                ma.presente.poly[[13]],
                                ma.presente.poly[[14]],
                                ma.presente.poly[[15]],
                                ma.presente.poly[[18]],
                                ma.presente.poly[[19]])




names(ma.presente.poly.ing.thi)


myRespName.ing.thi <-  'ing.thi'
myResp.ing.thi<-rep(1, nrow(coord.ing.thi.names_thinned_dataset_full[[1]]))
myResp.ing.thi
length(myResp.ing.thi)
myRespXY.ing.thi<-coord.ing.thi.names_thinned_dataset_full[[1]]

# function to define the intersect of rasters
intersect_mask <- function(x){
  values_x <- getValues(x)
  inter_x <- values_x %*% rep(1,nlayers(x))
  mask <- setValues(subset(x,1),values = (inter_x>0))
  return(mask)
}
# keep only all cells that are defined for all layers

myExpl.ing.thi<-stack(reg.predictors.ing.thi)

myExpl.ing.thi <- stack(mask(myExpl.ing.thi, intersect_mask(myExpl.ing.thi)))


##Number of pseudo-absence: Lobo & Tognelli sugerem 100 pontos de pseudo-
#ausencia para cada ponto de presenca.



myBiomodData.ing.thi<- BIOMOD_FormatingData(resp.var = myResp.ing.thi,
                                            expl.var = myExpl.ing.thi,
                                            resp.xy = myRespXY.ing.thi,
                                            resp.name = myRespName.ing.thi,
                                            PA.nb.rep = 1,
                                            PA.nb.absences = 10000,
                                            PA.strategy = 'disk',
                                            PA.dist.max=500000) # com base em http://www.sciencedirect.com/science/article/pii/S0304380008005486
myBiomodData.ing.thi
save.image()
###MODELING

# 1. Defining Models Options using default options.
myBiomodOption.ing.thi <- BIOMOD_ModelingOptions()

# 3. Computing the models

myBiomodModelOut.ing.thi<- BIOMOD_Modeling(myBiomodData.ing.thi,
                                           models = c("MAXENT.Phillips"),
                                           models.options = myBiomodOption.ing.thi,
                                           NbRunEval=5, #conforme http://onlinelibrary.wiley.com/doi/10.1111/j.1523-1739.2006.00354.x/full
                                           DataSplit=70,
                                           Prevalence=0.5,
                                           VarImport=1,
                                           models.eval.meth = c('TSS'),
                                           ing.thiveObj = TRUE,
                                           rescal.all.models = TRUE,
                                           do.full.models = FALSE,
                                           modeling.id = paste(myRespName.ing.thi,"FirstModeling",sep=""))
(myBiomodModelOut.ing.thi)
save.image()
# get all models evaluation

#Excelente TSS>0.75
#Bom 0.40>=TSS<0.75
#Ruim TSS<0.40


myBiomodModelEval.ing.thi <- get_evaluations(myBiomodModelOut.ing.thi)

myBiomodModelEval.ing.thi


#####STOP NON TSS > 0.4




citations.ing.thi<-BIEN_metadata_citation(dataframe=ing.thi.bien)#If you are referencing occurrence data
citations.ing.thi
