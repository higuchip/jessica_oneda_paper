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
library(ff)

BioC<- getData('worldclim', var='bio', res=10) # Bioclima
variaveis <- stack(BioC)

plot(wrld_simpl, xlim=c(-130,-35), ylim=c(-45,35), axes=TRUE, col="gray")
plot(wrld_simpl, xlim=c(-70,-45), ylim=c(-35,-5), axes=TRUE, col="gray")
points(coord.lec.pis, pch=19, cex=.51, col="black")


##Poligono MA

require(rgdal)
data.shape<-readOGR(dsn="/media/dendrologia/423EF1863EF172F1/User/Documents/orientados/Guilherme/IDS_09_Desm_Bioma_Mata_Atlantica_1_2", 
                    layer="IDS_09_Desm_Bioma_Mata_Atlantica_1_2")

poligono_mata_atlantica<-data.shape@polygons[[1]]@Polygons[[1]]@coords

polygon(poligono_mata_atlantica)
###Poligono especie


poligono_lec.pis<-getpoly()

colnames(poligono_lec.pis) = c("Long", "Lat")
polygon(poligono_lec.pis) 



########REMOCAO DUVIDAS

out.lec.pis= pnt.in.poly(coord.lec.pis,poligono_lec.pis)
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
plot(data.shape, add=T, col="lightgreen")
points(out.lec.pis[which(out.lec.pis$pip==1),1:2],pch=20, cex=.5, col="black")
coord.lec.pis.limpo<-out.lec.pis[which(out.lec.pis$pip==1),1:2]
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
points(coord.lec.pis.limpo,pch=20, cex=.5, col="black")


###FILTRAGEM ESPACIAl

head(coord.lec.pis.limpo)
coord.lec.pis.limpo<-na.omit(coord.lec.pis.limpo)

dim(coord.lec.pis.limpo)
lec.pis_names<-c(rep("lec.pis",  dim(coord.lec.pis.limpo)[1]))
coord.lec.pis.limpo.names<-cbind(lec.pis_names, coord.lec.pis.limpo)

head(coord.lec.pis.limpo.names)
is.data.frame(coord.lec.pis.limpo.names)
colnames(coord.lec.pis.limpo.names)<-c("SPEC", "LONG", "LAT")
# 
# ?write.table
# write.table(coord.lec.pis.limpo.names, "tap_gui.csv", sep=";", dec=",")

coord.lec.pis.names_thinned_dataset_full <-
  thin( loc.data =coord.lec.pis.limpo.names, 
        lat.col = "LAT", 
        long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 26.27, reps =1, #10 min resolution
        locs.thinned.list.return = TRUE, 
        write.files = T, 
        max.files = 5, 
        out.dir = "coord.lec.pis.limpo.names_full/", out.base = "coord.lec.pis.limpo.names_thinned", 
        write.log.file = TRUE,
        log.file = "coord.lec.pis.limpo.names_thinned_full_log_file.txt" )



dim(coord.lec.pis.limpo.names)
dim(coord.lec.pis.names_thinned_dataset_full[[1]])



###ENVIRONMENTAL BACKGROUND

reg.variaveis.lec.pis= crop(variaveis, poligono_lec.pis)
plot(reg.variaveis.lec.pis, 1)
points(coord.lec.pis.limpo, pch=20, cex=.5)

#remoção variaveis nao colineares
vari.n.colineares.lec.pis<-vifstep(reg.variaveis.lec.pis)
vari.n.colineares.lec.pis


reg.predictors.lec.pis<- stack(reg.variaveis.lec.pis[[2]],
                               reg.variaveis.lec.pis[[3]],
                               reg.variaveis.lec.pis[[4]],
                               reg.variaveis.lec.pis[[8]],
                               reg.variaveis.lec.pis[[13]],
                               reg.variaveis.lec.pis[[14]],
                               reg.variaveis.lec.pis[[15]],
                               reg.variaveis.lec.pis[[18]],
                               reg.variaveis.lec.pis[[19]])



names(reg.predictors.lec.pis)


###MATA ATLANTICA projection presente

ma.presente = crop(variaveis, poligono_mata_atlantica)
plot(ma.presente, 1)
polygon(poligono_mata_atlantica)
ma.presente.poly <- mask(ma.presente, data.shape)
plot(ma.presente.poly, 1 )


ma.presente.poly.lec.pis<-stack(ma.presente.new[[2]],
                                ma.presente.new[[3]],
                                ma.presente.new[[4]],
                                ma.presente.new[[8]],
                                ma.presente.new[[13]],
                                ma.presente.new[[14]],
                                ma.presente.new[[15]],
                                ma.presente.new[[18]],
                                ma.presente.new[[19]])




names(ma.presente.poly.lec.pis)


myRespName.lec.pis <-  'lec.pis'
myResp.lec.pis<-rep(1, nrow(coord.lec.pis.names_thinned_dataset_full[[1]]))
myResp.lec.pis
length(myResp.lec.pis)
myRespXY.lec.pis<-coord.lec.pis.names_thinned_dataset_full[[1]]

# function to define the intersect of rasters
intersect_mask <- function(x){
  values_x <- getValues(x)
  inter_x <- values_x %*% rep(1,nlayers(x))
  mask <- setValues(subset(x,1),values = (inter_x>0))
  return(mask)
}
# keep only all cells that are defined for all layers

myExpl.lec.pis<-stack(reg.predictors.lec.pis)

myExpl.lec.pis <- stack(mask(myExpl.lec.pis, intersect_mask(myExpl.lec.pis)))


##Number of pseudo-absence: Lobo & Tognelli sugerem 100 pontos de pseudo-
#ausencia para cada ponto de presenca.



myBiomodData.lec.pis<- BIOMOD_FormatingData(resp.var = myResp.lec.pis,
                                            expl.var = myExpl.lec.pis,
                                            resp.xy = myRespXY.lec.pis,
                                            resp.name = myRespName.lec.pis,
                                            PA.nb.rep = 1,
                                            PA.nb.absences = 10000,
                                            PA.strategy = 'disk',
                                            PA.dist.max=500000) # com base em http://www.sciencedirect.com/science/article/pii/S0304380008005486
myBiomodData.lec.pis
save.image()
###MODELING

# 1. Defining Models Options using default options.
myBiomodOption.lec.pis <- BIOMOD_ModelingOptions()

# 3. Computing the models

myBiomodModelOut.lec.pis<- BIOMOD_Modeling(myBiomodData.lec.pis,
                                           models = c("MAXENT.Phillips"),
                                           models.options = myBiomodOption.lec.pis,
                                           NbRunEval=5, #conforme http://onlinelibrary.wiley.com/doi/10.1111/j.1523-1739.2006.00354.x/full
                                           DataSplit=70,
                                           Prevalence=0.5,
                                           VarImport=1,
                                           models.eval.meth = c('TSS'),
                                           lec.pisveObj = TRUE,
                                           rescal.all.models = TRUE,
                                           do.full.models = FALSE,
                                           modeling.id = paste(myRespName.lec.pis,"FirstModeling",sep=""))
(myBiomodModelOut.lec.pis)
save.image()
# get all models evaluation

#Excelente TSS>0.75
#Bom 0.40>=TSS<0.75
#Ruim TSS<0.40


myBiomodModelEval.lec.pis <- get_evaluations(myBiomodModelOut.lec.pis)

myBiomodModelEval.lec.pis






get_variables_importance(myBiomodModelOut.lec.pis)
get_evaluations_matrix_lec.pis<-cbind(get_variables_importance(myBiomodModelOut.lec.pis)[1:9], #RUN 1
                                      get_variables_importance(myBiomodModelOut.lec.pis)[10:18], #RUN 2
                                      get_variables_importance(myBiomodModelOut.lec.pis)[19:27], #RUN 3
                                      get_variables_importance(myBiomodModelOut.lec.pis)[28:36], #RUN 4
                                      get_variables_importance(myBiomodModelOut.lec.pis)[37:45]) #RUN 5
colnames(get_evaluations_matrix_lec.pis)<-c("RUN1", "RUN2", "RUN3", "RUN4", "RUN5")
names(reg.predictors.lec.pis)
rownames(get_evaluations_matrix_lec.pis)<-names(reg.predictors.lec.pis)
get_evaluations_matrix_lec.pis_mean<-apply(get_evaluations_matrix_lec.pis, 1, mean)
which.max(get_evaluations_matrix_lec.pis_mean)



mymaxent.lec.pis <- BIOMOD_LoadModels(myBiomodModelOut.lec.pis, models='MAXENT.Phillips')
myRespPlot2D.lec.pis <- response.plot2(models  = mymaxent.lec.pis,
                                       Data = get_formal_data(myBiomodModelOut.lec.pis,'expl.var'), 
                                       show.variables= get_formal_data(myBiomodModelOut.lec.pis,'expl.var.names'),
                                       do.bivariate = FALSE,
                                       fixed.var.metric = 'mean',
                                       col = c("blue", "red"),
                                       legend = F,
                                       data_species = get_formal_data(myBiomodModelOut.lec.pis,'resp.var'),
                                       plot=F)

save.image()

bio_14_response_lec.pis<-gather(myRespPlot2D.lec.pis$bio14,model,Probability,
                                lec.pis_PA1_RUN1_MAXENT.Phillips:lec.pis_PA1_RUN5_MAXENT.Phillips)



p_lec.pis_bio14<-ggplot(data =bio_14_response_lec.pis) + 
  geom_smooth(mapping = aes(x = bio14, y = Probability))+
  labs(x = bquote('Bio 14'))
p_lec.pis_bio14



###Ensemble modeling
myBiomodEM.lec.pis<-BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut.lec.pis,
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

myBiomodProjection.lec.pis <- BIOMOD_Projection(modeling.output = myBiomodModelOut.lec.pis,
                                                new.env = stack(ma.presente.poly.lec.pis),
                                                proj.name = 'current_MA',
                                                selected.models = "all",
                                                binary.meth = 'TSS',
                                                compress = FALSE,
                                                build.clamping.mask = FALSE)

#Ensemble Forcasting

myBiomodEF.lec.pis.presente <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM.lec.pis,
  projection.output = myBiomodProjection.lec.pis,binary.meth = "TSS")

save.image()

plot(myBiomodEF.lec.pis.presente)

currentPred.ensemble.lec.pis <- stack("lec.pis/proj_current_MA/proj_current_MA_lec.pis_ensemble.grd")
plot(currentPred.ensemble.lec.pis[[1]])

plot(currentPred.ensemble.lec.pis[[1]],col = gray.colors(12, start = 0.9, end = 0.0))
plot(wrld_simpl, xlim=c(-60,-40), ylim=c(-34,-15), axes=TRUE, add=T)
polygon(poligono_mata_atlantica)
points(coord.lec.pis.names_thinned_dataset_full[[1]], cex=.6,col="black")

plot(currentPred.ensemble.lec.pis[[2]],col = gray.colors(12, start = 0.9, end = 0.0))
points(coord.lec.pis.names_thinned_dataset_full[[1]], cex=.6,col="black")
polygon(poligono_mata_atlantica)


?BIEN_metadata_citation

citations.lec.pis<-BIEN_metadata_citation(dataframe=lec.pis.bien)#If you are referencing occurrence data
citations.lec.pis


########################## 2070 ##########################
################################################################



ma.2070.poly.lec.pis<-stack(ma.2070.new[[2]],
                                ma.2070.new[[3]],
                                ma.2070.new[[4]],
                                ma.2070.new[[8]],
                                ma.2070.new[[13]],
                                ma.2070.new[[14]],
                                ma.2070.new[[15]],
                                ma.2070.new[[18]],
                                ma.2070.new[[19]])




names(ma.2070.poly.lec.pis)<-names(ma.presente.poly.lec.pis)


myExplFuture.2070.85.ma<-stack(ma.2070.poly.lec.pis) #modelo HE


projection.2070.85.ma.lec.pis <-BIOMOD_Projection(modeling.output=myBiomodModelOut.lec.pis,
                                                  new.env= myExplFuture.2070.85.ma,
                                                  proj.name='lec.pis.2070.85',
                                                  selected.models =  "all",
                                                  Bin.trans=TRUE,
                                                  binary.meth ='TSS',
                                                  compress = 'gzip',
                                                  build.clamping.mask	= FALSE,
                                                  SaveObj=TRUE,
                                                  keep.in.memory=F,
                                                  do.stack=T,
                                                  silent=F)


lec.pis.2070.85<-get_predictions(projection.2070.85.ma.lec.pis)
lec.pis.2070.85
summary(lec.pis.2070.85)





################## ENSAMBLE MODELLING 2070  ################
#################################################################



EnsambleForecast.2070.85.ma.lec.pis <-BIOMOD_EnsembleForecasting (
  projection.output = projection.2070.85.ma.lec.pis,
  EM.output=myBiomodEM.lec.pis,
  selected.models = "all",
  binary.meth ='TSS') 




EnsambleForecast.2070.85.ma.lec.pis
plot(EnsambleForecast.2070.85.ma.lec.pis)



ensemble.2070.HE.85.ma.lec.pis <- stack("lec.pis/proj_lec.pis.2070.85/proj_lec.pis.2070.85_lec.pis_ensemble.grd")
ensemble.2070.HE.85.ma.lec.pis



plot(ensemble.2070.HE.85.ma.lec.pis[[1]],col = gray.colors(12, start = 0.9, end = 0))
points(coord.gym_klo.names_thinned_dataset_full[[1]],pch=20, col="black", cex=.7)

plot(ensemble.2070.HE.85.ma.lec.pis[[2]],col = gray.colors(12, start = 0.9, end = 0))



futurePred.ensemble.2070.85.lec.pis <- stack("lec.pis/proj_lec.pis.2070.85/proj_lec.pis.2070.85_lec.pis_ensemble.grd")



plot(currentPred.ensemble.lec.pis[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.2070.85.lec.pis[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS


currentPred.ensemble.bin.lec.pis<- stack("lec.pis/proj_current_MA/proj_current_MA_lec.pis_ensemble_TSSbin.grd")
futurePred.ensemble.bin.2070.85.lec.pis <- stack("lec.pis/proj_lec.pis.2070.85/proj_lec.pis.2070.85_lec.pis_ensemble_TSSbin.gri")


plot(currentPred.ensemble.bin.lec.pis[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.bin.2070.85.lec.pis[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS



myBiomodRangeSize.2070.85.lec.pis <- BIOMOD_RangeSize(
  CurrentPred=currentPred.ensemble.bin.lec.pis[[1]],
  FutureProj=futurePred.ensemble.bin.2070.85.lec.pis[[1]])
myBiomodRangeSize.2070.85.lec.pis$Compt.By.Models





plot(myBiomodRangeSize.2070.85.lec.pis$Diff.By.Pixel, col= gray.colors(4, start = 0, end = 1, gamma = 1, alpha = NULL))
plot(wrld_simpl,  axes="TRUE", add=T)

polygon(poligono_mata_atlantica)








save.image()




