####ADEQUACAO

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
points(coord.gar.gar, pch=19, cex=.51, col="black")


##Poligono MA

require(rgdal)
data.shape<-readOGR(dsn="/media/dendrologia/423EF1863EF172F1/User/Documents/orientados/Guilherme/IDS_09_Desm_Bioma_Mata_Atlantica_1_2", 
                    layer="IDS_09_Desm_Bioma_Mata_Atlantica_1_2")

poligono_mata_atlantica<-data.shape@polygons[[1]]@Polygons[[1]]@coords

polygon(poligono_mata_atlantica)
###Poligono especie


poligono_gar.gar<-getpoly()

colnames(poligono_gar.gar) = c("Long", "Lat")
polygon(poligono_gar.gar) 



########REMOCAO DUVIDAS

out.gar.gar= pnt.in.poly(coord.gar.gar,poligono_gar.gar)
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
plot(data.shape, add=T, col="lightgreen")
points(out.gar.gar[which(out.gar.gar$pip==1),1:2],pch=20, cex=.5, col="black")
coord.gar.gar.limpo<-out.gar.gar[which(out.gar.gar$pip==1),1:2]
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
points(coord.gar.gar.limpo,pch=20, cex=.5, col="black")


###FILTRAGEM ESPACIAl

head(coord.gar.gar.limpo)
coord.gar.gar.limpo<-na.omit(coord.gar.gar.limpo)

dim(coord.gar.gar.limpo)
gar.gar_names<-c(rep("gar.gar",  dim(coord.gar.gar.limpo)[1]))
coord.gar.gar.limpo.names<-cbind(gar.gar_names, coord.gar.gar.limpo)

head(coord.gar.gar.limpo.names)
is.data.frame(coord.gar.gar.limpo.names)
colnames(coord.gar.gar.limpo.names)<-c("SPEC", "LONG", "LAT")
# 
# ?write.table
# write.table(coord.gar.gar.limpo.names, "tap_gui.csv", sep=";", dec=",")




coord.gar.gar.names_thinned_dataset_full <-
  thin( loc.data =coord.gar.gar.limpo.names, 
        lat.col = "LAT", 
        long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 26.27, reps =1, #10 min resolution
        locs.thinned.list.return = TRUE, 
        write.files = T, 
        max.files = 5, 
        out.dir = "coord.gar.gar.limpo.names_full/", out.base = "coord.gar.gar.limpo.names_thinned", 
        write.log.file = TRUE,
        log.file = "coord.gar.gar.limpo.names_thinned_full_log_file.txt" )



dim(coord.gar.gar.limpo.names)
dim(coord.gar.gar.names_thinned_dataset_full[[1]])



###ENVIRONMENTAL BACKGROUND

reg.variaveis.gar.gar= crop(variaveis, poligono_gar.gar)
plot(reg.variaveis.gar.gar, 1)
points(coord.gar.gar.limpo, pch=20, cex=.5)

#remoção variaveis nao colineares
vari.n.colineares.gar.gar<-vifstep(reg.variaveis.gar.gar)
vari.n.colineares.gar.gar


reg.predictors.gar.gar<- stack(reg.variaveis.gar.gar[[2]],
                               reg.variaveis.gar.gar[[3]],
                               reg.variaveis.gar.gar[[8]],
                               reg.variaveis.gar.gar[[9]],
                               reg.variaveis.gar.gar[[13]],
                               reg.variaveis.gar.gar[[14]],
                               reg.variaveis.gar.gar[[15]],
                               reg.variaveis.gar.gar[[18]],
                               reg.variaveis.gar.gar[[19]])



names(reg.predictors.gar.gar)


###MATA ATLANTICA projection presente

ma.presente = crop(variaveis, poligono_mata_atlantica)
plot(ma.presente, 1)
polygon(poligono_mata_atlantica)
ma.presente.poly <- mask(ma.presente, data.shape)
plot(ma.presente.poly, 1 )


ma.presente.poly.gar.gar<-stack(ma.presente.new[[2]],
                                ma.presente.new[[3]],
                                ma.presente.new[[8]],
                                ma.presente.new[[9]],
                                ma.presente.new[[13]],
                                ma.presente.new[[14]],
                                ma.presente.new[[15]],
                                ma.presente.new[[18]],
                                ma.presente.new[[19]])




names(ma.presente.poly.gar.gar)


myRespName.gar.gar <-  'gar.gar'
myResp.gar.gar<-rep(1, nrow(coord.gar.gar.names_thinned_dataset_full[[1]]))
myRespXY.gar.gar<-coord.gar.gar.names_thinned_dataset_full[[1]]

# function to define the intersect of rasters
intersect_mask <- function(x){
  values_x <- getValues(x)
  inter_x <- values_x %*% rep(1,nlayers(x))
  mask <- setValues(subset(x,1),values = (inter_x>0))
  return(mask)
}
# keep only all cells that are defined for all layers

myExpl.gar.gar<-stack(reg.predictors.gar.gar)

myExpl.gar.gar <- stack(mask(myExpl.gar.gar, intersect_mask(myExpl.gar.gar)))


##Number of pseudo-absence: Lobo & Tognelli sugerem 100 pontos de pseudo-
#ausencia para cada ponto de presenca.



myBiomodData.gar.gar<- BIOMOD_FormatingData(resp.var = myResp.gar.gar,
                                            expl.var = myExpl.gar.gar,
                                            resp.xy = myRespXY.gar.gar,
                                            resp.name = myRespName.gar.gar,
                                            PA.nb.rep = 1,
                                            PA.nb.absences = 10000,
                                            PA.strategy = 'disk',
                                            PA.dist.max=500000) # com base em http://www.sciencedirect.com/science/article/pii/S0304380008005486
myBiomodData.gar.gar
save.image()
###MODELING

# 1. Defining Models Options using default options.
myBiomodOption.gar.gar <- BIOMOD_ModelingOptions()

# 3. Computing the models

myBiomodModelOut.gar.gar<- BIOMOD_Modeling(myBiomodData.gar.gar,
                                           models = c("MAXENT.Phillips"),
                                           models.options = myBiomodOption.gar.gar,
                                           NbRunEval=5, #conforme http://onlinelibrary.wiley.com/doi/10.1111/j.1523-1739.2006.00354.x/full
                                           DataSplit=70,
                                           Prevalence=0.5,
                                           VarImport=1,
                                           models.eval.meth = c('TSS'),
                                           gar.garveObj = TRUE,
                                           rescal.all.models = TRUE,
                                           do.full.models = FALSE,
                                           modeling.id = paste(myRespName.gar.gar,"FirstModeling",sep=""))
(myBiomodModelOut.gar.gar)
save.image()
# get all models evaluation

#Excelente TSS>0.75
#Bom 0.40>=TSS<0.75
#Ruim TSS<0.40


myBiomodModelEval.gar.gar <- get_evaluations(myBiomodModelOut.gar.gar)

myBiomodModelEval.gar.gar




get_variables_importance(myBiomodModelOut.gar.gar)
get_evaluations_matrix_gar.gar<-cbind(get_variables_importance(myBiomodModelOut.gar.gar)[1:9], #RUN 1
                                      get_variables_importance(myBiomodModelOut.gar.gar)[10:18], #RUN 2
                                      get_variables_importance(myBiomodModelOut.gar.gar)[19:27], #RUN 3
                                      get_variables_importance(myBiomodModelOut.gar.gar)[28:36], #RUN 4
                                      get_variables_importance(myBiomodModelOut.gar.gar)[37:45]) #RUN 5
colnames(get_evaluations_matrix_gar.gar)<-c("RUN1", "RUN2", "RUN3", "RUN4", "RUN5")
names(reg.predictors.gar.gar)
rownames(get_evaluations_matrix_gar.gar)<-names(reg.predictors.gar.gar)
get_evaluations_matrix_gar.gar_mean<-apply(get_evaluations_matrix_gar.gar, 1, mean)
which.max(get_evaluations_matrix_gar.gar_mean)



mymaxent.gar.gar <- BIOMOD_LoadModels(myBiomodModelOut.gar.gar, models='MAXENT.Phillips')
myRespPlot2D.gar.gar <- response.plot2(models  = mymaxent.gar.gar,
                                       Data = get_formal_data(myBiomodModelOut.gar.gar,'expl.var'), 
                                       show.variables= get_formal_data(myBiomodModelOut.gar.gar,'expl.var.names'),
                                       do.bivariate = FALSE,
                                       fixed.var.metric = 'mean',
                                       col = c("blue", "red"),
                                       legend = F,
                                       data_species = get_formal_data(myBiomodModelOut.gar.gar,'resp.var'),
                                       plot=F)

save.image()

bio_3_response_gar.gar<-gather(myRespPlot2D.gar.gar$bio3,model,Probability,
                                gar.gar_PA1_RUN1_MAXENT.Phillips:gar.gar_PA1_RUN5_MAXENT.Phillips)



p_gar.gar_bio3<-ggplot(data =bio_3_response_gar.gar) + 
  geom_smooth(mapping = aes(x = bio3, y = Probability))+
  labs(x = bquote('Bio 3'))
p_gar.gar_bio3



###Ensemble modeling
myBiomodEM.gar.gar<-BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut.gar.gar,
                                             chosen.models = 'all',
                                             em.by = 'all',
                                             eval.metric = "all",
                                             eval.metric.quality.threshold = c(0.40),
                                             models.eval.meth = c('TSS'),
                                             prob.mean = TRUE,
                                             prob.cv = TRUE,
                                             prob.ci = FALSE,
                                             prob.ci.alpha = 0.05,
                                             prob.median = FALSE,
                                             committee.averaging = FALSE,
                                             prob.mean.weight = FALSE,
                                             prob.mean.weight.decay = 'proportional',
                                             VarImport = 3)

# Individual models projections on current environmental conditions
plot(ma.presente.poly.gar.gar,1)

myBiomodProjection.gar.gar <- BIOMOD_Projection(modeling.output = myBiomodModelOut.gar.gar,
                                                new.env = stack(ma.presente.poly.gar.gar),
                                                proj.name = 'current_MA',
                                                selected.models = "all",
                                                binary.meth = 'TSS',
                                                compress = FALSE,
                                                build.clamping.mask = FALSE)

#Ensemble Forcasting

myBiomodEF.gar.gar.presente <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM.gar.gar,
  projection.output = myBiomodProjection.gar.gar,binary.meth = "TSS")



plot(myBiomodEF.gar.gar.presente)

currentPred.ensemble.gar.gar <- stack("gar.gar/proj_current_MA/proj_current_MA_gar.gar_ensemble.grd")
plot(currentPred.ensemble.gar.gar[[1]])

plot(currentPred.ensemble.gar.gar[[1]],col = gray.colors(12, start = 0.9, end = 0.0))
plot(wrld_simpl, xlim=c(-60,-40), ylim=c(-34,-15), axes=TRUE, add=T)
polygon(poligono_mata_atlantica)
points(coord.gar.gar.names_thinned_dataset_full[[1]], cex=.6,col="black")

plot(currentPred.ensemble.gar.gar[[2]],col = gray.colors(12, start = 0.9, end = 0.0))
points(coord.gar.gar.names_thinned_dataset_full[[1]], cex=.6,col="black")
polygon(poligono_mata_atlantica)



citations.gar.gar<-BIEN_metadata_citation(dataframe=gar.gar.bien)#If you are referencing occurrence data
citations.gar.gar

save.image()


######################### 2070 ##########################
################################################################

ma.2070.poly.gar.gar<-stack(ma.2070.new[[2]],
                                ma.2070.new[[3]],
                                ma.2070.new[[8]],
                                ma.2070.new[[9]],
                                ma.2070.new[[13]],
                                ma.2070.new[[14]],
                                ma.2070.new[[15]],
                                ma.2070.new[[18]],
                                ma.2070.new[[19]])

names(ma.2070.poly.gar.gar)<-names(ma.presente.poly.gar.gar)


myExplFuture.2070.85.ma<-stack(ma.2070.poly.gar.gar) #modelo HE


projection.2070.85.ma.gar.gar <-BIOMOD_Projection(modeling.output=myBiomodModelOut.gar.gar,
                                                  new.env= myExplFuture.2070.85.ma,
                                                  proj.name='gar.gar.2070.85',
                                                  selected.models =  "all",
                                                  Bin.trans=TRUE,
                                                  binary.meth ='TSS',
                                                  compress = 'gzip',
                                                  build.clamping.mask	= FALSE,
                                                  SaveObj=TRUE,
                                                  keep.in.memory=F,
                                                  do.stack=T,
                                                  silent=F)


gar.gar.2070.85<-get_predictions(projection.2070.85.ma.gar.gar)
gar.gar.2070.85
summary(gar.gar.2070.85)





################## ENSAMBLE MODELLING 2070  ################
#################################################################



EnsambleForecast.2070.85.ma.gar.gar <-BIOMOD_EnsembleForecasting (
  projection.output = projection.2070.85.ma.gar.gar,
  EM.output=myBiomodEM.gar.gar,
  selected.models = "all",
  binary.meth ='TSS') 




EnsambleForecast.2070.85.ma.gar.gar
plot(EnsambleForecast.2070.85.ma.gar.gar)



ensemble.2070.HE.85.ma.gar.gar <- stack("gar.gar/proj_gar.gar.2070.85/proj_gar.gar.2070.85_gar.gar_ensemble.grd")
ensemble.2070.HE.85.ma.gar.gar



plot(ensemble.2070.HE.85.ma.gar.gar[[1]],col = gray.colors(12, start = 0.9, end = 0))
points(coord.gar.gar.names_thinned_dataset_full[[1]],pch=20, col="black", cex=.7)

plot(ensemble.2070.HE.85.ma.gar.gar[[2]],col = gray.colors(12, start = 0.9, end = 0))



futurePred.ensemble.2070.85.gar.gar <- stack("gar.gar/proj_gar.gar.2070.85/proj_gar.gar.2070.85_gar.gar_ensemble.grd")



plot(currentPred.ensemble.gar.gar[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.2070.85.gar.gar[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS


currentPred.ensemble.bin.gar.gar<- stack("gar.gar/proj_current_MA/proj_current_MA_gar.gar_ensemble_TSSbin.grd")
futurePred.ensemble.bin.2070.85.gar.gar <- stack("gar.gar/proj_gar.gar.2070.85/proj_gar.gar.2070.85_gar.gar_ensemble_TSSbin.gri")


plot(currentPred.ensemble.bin.gar.gar[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.bin.2070.85.gar.gar[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS



myBiomodRangeSize.2070.85.gar.gar <- BIOMOD_RangeSize(
  CurrentPred=currentPred.ensemble.bin.gar.gar[[1]],
  FutureProj=futurePred.ensemble.bin.2070.85.gar.gar[[1]])
myBiomodRangeSize.2070.85.gar.gar$Compt.By.Models





plot(myBiomodRangeSize.2070.85.gar.gar$Diff.By.Pixel, col= gray.colors(4, start = 0, end = 1, gamma = 1, alpha = NULL))
plot(wrld_simpl,  axes="TRUE", add=T)

polygon(poligono_mata_atlantica)








save.image()










