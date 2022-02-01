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
points(coord.cou.mic, pch=19, cex=.51, col="black")


##Poligono MA

require(rgdal)
data.shape<-readOGR(dsn="/media/dendrologia/423EF1863EF172F1/User/Documents/orientados/Guilherme/IDS_09_Desm_Bioma_Mata_Atlantica_1_2", 
                    layer="IDS_09_Desm_Bioma_Mata_Atlantica_1_2")

poligono_mata_atlantica<-data.shape@polygons[[1]]@Polygons[[1]]@coords

polygon(poligono_mata_atlantica)
###Poligono especie


poligono_cou.mic<-getpoly()

colnames(poligono_cou.mic) = c("Long", "Lat")
polygon(poligono_cou.mic) 



########REMOCAO DUVIDAS

out.cou.mic= pnt.in.poly(coord.cou.mic,poligono_cou.mic)
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
plot(data.shape, add=T, col="lightgreen")
points(out.cou.mic[which(out.cou.mic$pip==1),1:2],pch=20, cex=.5, col="black")
coord.cou.mic.limpo<-out.cou.mic[which(out.cou.mic$pip==1),1:2]
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
points(coord.cou.mic.limpo,pch=20, cex=.5, col="black")


###FILTRAGEM ESPACIAl

head(coord.cou.mic.limpo)
coord.cou.mic.limpo<-na.omit(coord.cou.mic.limpo)

dim(coord.cou.mic.limpo)
cou.mic_names<-c(rep("cou.mic",  dim(coord.cou.mic.limpo)[1]))
coord.cou.mic.limpo.names<-cbind(cou.mic_names, coord.cou.mic.limpo)

head(coord.cou.mic.limpo.names)
is.data.frame(coord.cou.mic.limpo.names)
colnames(coord.cou.mic.limpo.names)<-c("SPEC", "LONG", "LAT")
# 
# ?write.table
# write.table(coord.cou.mic.limpo.names, "tap_gui.csv", sep=";", dec=",")




coord.cou.mic.names_thinned_dataset_full <-
  thin( loc.data =coord.cou.mic.limpo.names, 
        lat.col = "LAT", 
        long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 26.27, reps =1, #10 min resolution
        locs.thinned.list.return = TRUE, 
        write.files = T, 
        max.files = 5, 
        out.dir = "coord.cou.mic.limpo.names_full/", out.base = "coord.cou.mic.limpo.names_thinned", 
        write.log.file = TRUE,
        log.file = "coord.cou.mic.limpo.names_thinned_full_log_file.txt" )



dim(coord.cou.mic.limpo.names)
dim(coord.cou.mic.names_thinned_dataset_full[[1]])



###ENVIRONMENTAL BACKGROUND

reg.variaveis.cou.mic= crop(variaveis, poligono_cou.mic)
plot(reg.variaveis.cou.mic, 1)
points(coord.cou.mic.limpo, pch=20, cex=.5)

#remoção variaveis nao colineares
vari.n.colineares.cou.mic<-vifstep(reg.variaveis.cou.mic)
vari.n.colineares.cou.mic


reg.predictors.cou.mic<- stack(reg.variaveis.cou.mic[[2]],
                               reg.variaveis.cou.mic[[3]],
                               reg.variaveis.cou.mic[[8]],
                               reg.variaveis.cou.mic[[9]],
                               reg.variaveis.cou.mic[[13]],
                               reg.variaveis.cou.mic[[14]],
                               reg.variaveis.cou.mic[[15]],
                               reg.variaveis.cou.mic[[18]],
                               reg.variaveis.cou.mic[[19]])



names(reg.predictors.cou.mic)


###MATA ATLANTICA projection presente

ma.presente = crop(variaveis, poligono_mata_atlantica)
plot(ma.presente, 1)
polygon(poligono_mata_atlantica)
ma.presente.poly <- mask(ma.presente, data.shape)
plot(ma.presente.poly, 1 )


ma.presente.poly.cou.mic<-stack(ma.presente.new[[2]],
                                ma.presente.new[[3]],
                                ma.presente.new[[8]],
                                ma.presente.new[[9]],
                                ma.presente.new[[13]],
                                ma.presente.new[[14]],
                                ma.presente.new[[15]],
                                ma.presente.new[[18]],
                                ma.presente.new[[19]])




names(ma.presente.poly.cou.mic)


myRespName.cou.mic <-  'cou.mic'
myResp.cou.mic<-rep(1, nrow(coord.cou.mic.names_thinned_dataset_full[[1]]))
myResp.cou.mic
length(myResp.cou.mic)
myRespXY.cou.mic<-coord.cou.mic.names_thinned_dataset_full[[1]]

# function to define the intersect of rasters
intersect_mask <- function(x){
  values_x <- getValues(x)
  inter_x <- values_x %*% rep(1,nlayers(x))
  mask <- setValues(subset(x,1),values = (inter_x>0))
  return(mask)
}
# keep only all cells that are defined for all layers

myExpl.cou.mic<-stack(reg.predictors.cou.mic)

myExpl.cou.mic <- stack(mask(myExpl.cou.mic, intersect_mask(myExpl.cou.mic)))


##Number of pseudo-absence: Lobo & Tognelli sugerem 100 pontos de pseudo-
#ausencia para cada ponto de presenca.



myBiomodData.cou.mic<- BIOMOD_FormatingData(resp.var = myResp.cou.mic,
                                            expl.var = myExpl.cou.mic,
                                            resp.xy = myRespXY.cou.mic,
                                            resp.name = myRespName.cou.mic,
                                            PA.nb.rep = 1,
                                            PA.nb.absences = 10000,
                                            PA.strategy = 'disk',
                                            PA.dist.max=500000) # com base em http://www.sciencedirect.com/science/article/pii/S0304380008005486
myBiomodData.cou.mic
save.image()
###MODELING

# 1. Defining Models Options using default options.
myBiomodOption.cou.mic <- BIOMOD_ModelingOptions()

# 3. Computing the models

myBiomodModelOut.cou.mic<- BIOMOD_Modeling(myBiomodData.cou.mic,
                                           models = c("MAXENT.Phillips"),
                                           models.options = myBiomodOption.cou.mic,
                                           NbRunEval=5, #conforme http://onlinelibrary.wiley.com/doi/10.1111/j.1523-1739.2006.00354.x/full
                                           DataSplit=70,
                                           Prevalence=0.5,
                                           VarImport=1,
                                           models.eval.meth = c('TSS'),
                                           cou.micveObj = TRUE,
                                           rescal.all.models = TRUE,
                                           do.full.models = FALSE,
                                           modeling.id = paste(myRespName.cou.mic,"FirstModeling",sep=""))
(myBiomodModelOut.cou.mic)
save.image()
# get all models evaluation

#Excelente TSS>0.75
#Bom 0.40>=TSS<0.75
#Ruim TSS<0.40


myBiomodModelEval.cou.mic <- get_evaluations(myBiomodModelOut.cou.mic)

myBiomodModelEval.cou.mic




get_variables_importance(myBiomodModelOut.cou.mic)
get_evaluations_matrix_cou.mic<-cbind(get_variables_importance(myBiomodModelOut.cou.mic)[1:9], #RUN 1
                                      get_variables_importance(myBiomodModelOut.cou.mic)[10:18], #RUN 2
                                      get_variables_importance(myBiomodModelOut.cou.mic)[19:27], #RUN 3
                                      get_variables_importance(myBiomodModelOut.cou.mic)[28:36], #RUN 4
                                      get_variables_importance(myBiomodModelOut.cou.mic)[37:45]) #RUN 5
colnames(get_evaluations_matrix_cou.mic)<-c("RUN1", "RUN2", "RUN3", "RUN4", "RUN5")
names(reg.predictors.cou.mic)
rownames(get_evaluations_matrix_cou.mic)<-names(reg.predictors.cou.mic)
get_evaluations_matrix_cou.mic_mean<-apply(get_evaluations_matrix_cou.mic, 1, mean)
which.max(get_evaluations_matrix_cou.mic_mean)



mymaxent.cou.mic <- BIOMOD_LoadModels(myBiomodModelOut.cou.mic, models='MAXENT.Phillips')
myRespPlot2D.cou.mic <- response.plot2(models  = mymaxent.cou.mic,
                                       Data = get_formal_data(myBiomodModelOut.cou.mic,'expl.var'), 
                                       show.variables= get_formal_data(myBiomodModelOut.cou.mic,'expl.var.names'),
                                       do.bivariate = FALSE,
                                       fixed.var.metric = 'mean',
                                       col = c("blue", "red"),
                                       legend = F,
                                       data_species = get_formal_data(myBiomodModelOut.cou.mic,'resp.var'),
                                       plot=F)

save.image()

bio_3_response_cou.mic<-gather(myRespPlot2D.cou.mic$bio3,model,Probability,
                                cou.mic_PA1_RUN1_MAXENT.Phillips:cou.mic_PA1_RUN5_MAXENT.Phillips)



p_cou.mic_bio3<-ggplot(data =bio_3_response_cou.mic) + 
  geom_smooth(mapping = aes(x = bio3, y = Probability))+
  labs(x = bquote('Bio 3'))
p_cou.mic_bio3



###Ensemble modeling
myBiomodEM.cou.mic<-BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut.cou.mic,
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

myBiomodProjection.cou.mic <- BIOMOD_Projection(modeling.output = myBiomodModelOut.cou.mic,
                                                new.env = stack(ma.presente.poly.cou.mic),
                                                proj.name = 'current_MA',
                                                selected.models = "all",
                                                binary.meth = 'TSS',
                                                compress = FALSE,
                                                build.clamping.mask = FALSE)

#Ensemble Forcasting

myBiomodEF.cou.mic.presente <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM.cou.mic,
  projection.output = myBiomodProjection.cou.mic,binary.meth = "TSS")



plot(myBiomodEF.cou.mic.presente)

currentPred.ensemble.cou.mic <- stack("cou.mic/proj_current_MA/proj_current_MA_cou.mic_ensemble.grd")
plot(currentPred.ensemble.cou.mic[[1]])

plot(currentPred.ensemble.cou.mic[[1]],col = gray.colors(12, start = 0.9, end = 0.0))
plot(wrld_simpl, xlim=c(-60,-40), ylim=c(-34,-15), axes=TRUE, add=T)
polygon(poligono_mata_atlantica)
points(coord.cou.mic.names_thinned_dataset_full[[1]], cex=.6,col="black")

plot(currentPred.ensemble.cou.mic[[2]],col = gray.colors(12, start = 0.9, end = 0.0))
points(coord.cou.mic.names_thinned_dataset_full[[1]], cex=.6,col="black")
polygon(poligono_mata_atlantica)



citations.cou.mic<-BIEN_metadata_citation(dataframe=cou.mic.bien)#If you are referencing occurrence data
citations.cou.mic

save.image()



########################## 2070 ##########################
################################################################

ma.2070.poly.cou.mic<-stack(ma.2070.new[[2]],
                                ma.2070.new[[3]],
                                ma.2070.new[[8]],
                                ma.2070.new[[9]],
                                ma.2070.new[[13]],
                                ma.2070.new[[14]],
                                ma.2070.new[[15]],
                                ma.2070.new[[18]],
                                ma.2070.new[[19]])




names(ma.2070.poly.cou.mic)<-names(ma.presente.poly.cou.mic)


myExplFuture.2070.85.ma<-stack(ma.2070.poly.cou.mic) #modelo HE


projection.2070.85.ma.cou.mic <-BIOMOD_Projection(modeling.output=myBiomodModelOut.cou.mic,
                                                  new.env= myExplFuture.2070.85.ma,
                                                  proj.name='cou.mic.2070.85',
                                                  selected.models =  "all",
                                                  Bin.trans=TRUE,
                                                  binary.meth ='TSS',
                                                  compress = 'gzip',
                                                  build.clamping.mask	= FALSE,
                                                  SaveObj=TRUE,
                                                  keep.in.memory=F,
                                                  do.stack=T,
                                                  silent=F)


cou.mic.2070.85<-get_predictions(projection.2070.85.ma.cou.mic)
cou.mic.2070.85
summary(cou.mic.2070.85)





################## ENSAMBLE MODELLING 2070  ################
#################################################################



EnsambleForecast.2070.85.ma.cou.mic <-BIOMOD_EnsembleForecasting (
  projection.output = projection.2070.85.ma.cou.mic,
  EM.output=myBiomodEM.cou.mic,
  selected.models = "all",
  binary.meth ='TSS') 




EnsambleForecast.2070.85.ma.cou.mic
plot(EnsambleForecast.2070.85.ma.cou.mic)



ensemble.2070.HE.85.ma.cou.mic <- stack("cou.mic/proj_cou.mic.2070.85/proj_cou.mic.2070.85_cou.mic_ensemble.grd")
ensemble.2070.HE.85.ma.cou.mic



plot(ensemble.2070.HE.85.ma.cou.mic[[1]],col = gray.colors(12, start = 0.9, end = 0))
points(coord.cou.mic.names_thinned_dataset_full[[1]],pch=20, col="black", cex=.7)

plot(ensemble.2070.HE.85.ma.cou.mic[[2]],col = gray.colors(12, start = 0.9, end = 0))



futurePred.ensemble.2070.85.cou.mic <- stack("cou.mic/proj_cou.mic.2070.85/proj_cou.mic.2070.85_cou.mic_ensemble.grd")



plot(currentPred.ensemble.cou.mic[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.2070.85.cou.mic[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS


currentPred.ensemble.bin.cou.mic<- stack("cou.mic/proj_current_MA/proj_current_MA_cou.mic_ensemble_TSSbin.grd")
futurePred.ensemble.bin.2070.85.cou.mic <- stack("cou.mic/proj_cou.mic.2070.85/proj_cou.mic.2070.85_cou.mic_ensemble_TSSbin.gri")


plot(currentPred.ensemble.bin.cou.mic[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.bin.2070.85.cou.mic[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS



myBiomodRangeSize.2070.85.cou.mic <- BIOMOD_RangeSize(
  CurrentPred=currentPred.ensemble.bin.cou.mic[[1]],
  FutureProj=futurePred.ensemble.bin.2070.85.cou.mic[[1]])
myBiomodRangeSize.2070.85.cou.mic$Compt.By.Models





plot(myBiomodRangeSize.2070.85.cou.mic$Diff.By.Pixel, col= gray.colors(4, start = 0, end = 1, gamma = 1, alpha = NULL))
plot(wrld_simpl,  axes="TRUE", add=T)

polygon(poligono_mata_atlantica)








save.image()








