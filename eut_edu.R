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
points(coord.eut.edu, pch=19, cex=.51, col="black")


##Poligono MA

require(rgdal)
data.shape<-readOGR(dsn="/media/dendrologia/423EF1863EF172F1/User/Documents/orientados/Guilherme/IDS_09_Desm_Bioma_Mata_Atlantica_1_2", 
                    layer="IDS_09_Desm_Bioma_Mata_Atlantica_1_2")

poligono_mata_atlantica<-data.shape@polygons[[1]]@Polygons[[1]]@coords

polygon(poligono_mata_atlantica)
###Poligono especie


poligono_eut.edu<-getpoly()

colnames(poligono_eut.edu) = c("Long", "Lat")
polygon(poligono_eut.edu) 



########REMOCAO DUVIDAS

out.eut.edu= pnt.in.poly(coord.eut.edu,poligono_eut.edu)
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
plot(data.shape, add=T, col="lightgreen")
points(out.eut.edu[which(out.eut.edu$pip==1),1:2],pch=20, cex=.5, col="black")
coord.eut.edu.limpo<-out.eut.edu[which(out.eut.edu$pip==1),1:2]
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
points(coord.eut.edu.limpo,pch=20, cex=.5, col="black")


###FILTRAGEM ESPACIAl

head(coord.eut.edu.limpo)
coord.eut.edu.limpo<-na.omit(coord.eut.edu.limpo)

dim(coord.eut.edu.limpo)
eut.edu_names<-c(rep("eut.edu",  dim(coord.eut.edu.limpo)[1]))
coord.eut.edu.limpo.names<-cbind(eut.edu_names, coord.eut.edu.limpo)

head(coord.eut.edu.limpo.names)
is.data.frame(coord.eut.edu.limpo.names)
colnames(coord.eut.edu.limpo.names)<-c("SPEC", "LONG", "LAT")
# 
# ?write.table
# write.table(coord.eut.edu.limpo.names, "tap_gui.csv", sep=";", dec=",")




coord.eut.edu.names_thinned_dataset_full <-
  thin( loc.data =coord.eut.edu.limpo.names, 
        lat.col = "LAT", 
        long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 26.27, reps =1, #10 min resolution
        locs.thinned.list.return = TRUE, 
        write.files = T, 
        max.files = 5, 
        out.dir = "coord.eut.edu.limpo.names_full/", out.base = "coord.eut.edu.limpo.names_thinned", 
        write.log.file = TRUE,
        log.file = "coord.eut.edu.limpo.names_thinned_full_log_file.txt" )



dim(coord.eut.edu.limpo.names)
dim(coord.eut.edu.names_thinned_dataset_full[[1]])



###ENVIRONMENTAL BACKGROUND

reg.variaveis.eut.edu= crop(variaveis, poligono_eut.edu)
plot(reg.variaveis.eut.edu, 1)
points(coord.eut.edu.limpo, pch=20, cex=.5)

#remoção variaveis nao colineares
vari.n.colineares.eut.edu<-vifstep(reg.variaveis.eut.edu)
vari.n.colineares.eut.edu


reg.predictors.eut.edu<- stack(reg.variaveis.eut.edu[[2]],
                               reg.variaveis.eut.edu[[3]],
                               reg.variaveis.eut.edu[[8]],
                               reg.variaveis.eut.edu[[9]],
                               reg.variaveis.eut.edu[[10]],
                               reg.variaveis.eut.edu[[13]],
                               reg.variaveis.eut.edu[[14]],
                               reg.variaveis.eut.edu[[15]],
                               reg.variaveis.eut.edu[[18]],
                               reg.variaveis.eut.edu[[19]])



names(reg.predictors.eut.edu)


###MATA ATLANTICA projection presente

ma.presente = crop(variaveis, poligono_mata_atlantica)
plot(ma.presente, 1)
polygon(poligono_mata_atlantica)
ma.presente.poly <- mask(ma.presente, data.shape)
plot(ma.presente.poly, 1 )


ma.presente.poly.eut.edu<-stack(ma.presente.new[[2]],
                                ma.presente.new[[3]],
                                ma.presente.new[[8]],
                                ma.presente.new[[9]],
                                ma.presente.new[[10]],
                                ma.presente.new[[13]],
                                ma.presente.new[[14]],
                                ma.presente.new[[15]],
                                ma.presente.new[[18]],
                                ma.presente.new[[19]])




names(ma.presente.poly.eut.edu)


myRespName.eut.edu <-  'eut.edu'
myResp.eut.edu<-rep(1, nrow(coord.eut.edu.names_thinned_dataset_full[[1]]))
myResp.eut.edu
length(myResp.eut.edu)
myRespXY.eut.edu<-coord.eut.edu.names_thinned_dataset_full[[1]]

# function to define the intersect of rasters
intersect_mask <- function(x){
  values_x <- getValues(x)
  inter_x <- values_x %*% rep(1,nlayers(x))
  mask <- setValues(subset(x,1),values = (inter_x>0))
  return(mask)
}
# keep only all cells that are defined for all layers

myExpl.eut.edu<-stack(reg.predictors.eut.edu)

myExpl.eut.edu <- stack(mask(myExpl.eut.edu, intersect_mask(myExpl.eut.edu)))


##Number of pseudo-absence: Lobo & Tognelli sugerem 100 pontos de pseudo-
#ausencia para cada ponto de presenca.



myBiomodData.eut.edu<- BIOMOD_FormatingData(resp.var = myResp.eut.edu,
                                            expl.var = myExpl.eut.edu,
                                            resp.xy = myRespXY.eut.edu,
                                            resp.name = myRespName.eut.edu,
                                            PA.nb.rep = 1,
                                            PA.nb.absences = 10000,
                                            PA.strategy = 'disk',
                                            PA.dist.max=500000) # com base em http://www.sciencedirect.com/science/article/pii/S0304380008005486
myBiomodData.eut.edu
save.image()
###MODELING

# 1. Defining Models Options using default options.
myBiomodOption.eut.edu <- BIOMOD_ModelingOptions()

# 3. Computing the models

myBiomodModelOut.eut.edu<- BIOMOD_Modeling(myBiomodData.eut.edu,
                                           models = c("MAXENT.Phillips"),
                                           models.options = myBiomodOption.eut.edu,
                                           NbRunEval=5, #conforme http://onlinelibrary.wiley.com/doi/10.1111/j.1523-1739.2006.00354.x/full
                                           DataSplit=70,
                                           Prevalence=0.5,
                                           VarImport=1,
                                           models.eval.meth = c('TSS'),
                                           eut.eduveObj = TRUE,
                                           rescal.all.models = TRUE,
                                           do.full.models = FALSE,
                                           modeling.id = paste(myRespName.eut.edu,"FirstModeling",sep=""))
(myBiomodModelOut.eut.edu)
save.image()
# get all models evaluation

#Excelente TSS>0.75
#Bom 0.40>=TSS<0.75
#Ruim TSS<0.40


myBiomodModelEval.eut.edu <- get_evaluations(myBiomodModelOut.eut.edu)

myBiomodModelEval.eut.edu




get_variables_importance(myBiomodModelOut.eut.edu)
get_evaluations_matrix_eut.edu<-cbind(get_variables_importance(myBiomodModelOut.eut.edu)[1:10], #RUN 1
                                      get_variables_importance(myBiomodModelOut.eut.edu)[11:20], #RUN 2
                                      get_variables_importance(myBiomodModelOut.eut.edu)[21:30], #RUN 3
                                      get_variables_importance(myBiomodModelOut.eut.edu)[31:40], #RUN 4
                                      get_variables_importance(myBiomodModelOut.eut.edu)[41:50]) #RUN 5
colnames(get_evaluations_matrix_eut.edu)<-c("RUN1", "RUN2", "RUN3", "RUN4", "RUN5")
names(reg.predictors.eut.edu)
rownames(get_evaluations_matrix_eut.edu)<-names(reg.predictors.eut.edu)
get_evaluations_matrix_eut.edu_mean<-apply(get_evaluations_matrix_eut.edu, 1, mean)
which.max(get_evaluations_matrix_eut.edu_mean)



mymaxent.eut.edu <- BIOMOD_LoadModels(myBiomodModelOut.eut.edu, models='MAXENT.Phillips')
myRespPlot2D.eut.edu <- response.plot2(models  = mymaxent.eut.edu,
                                       Data = get_formal_data(myBiomodModelOut.eut.edu,'expl.var'), 
                                       show.variables= get_formal_data(myBiomodModelOut.eut.edu,'expl.var.names'),
                                       do.bivariate = FALSE,
                                       fixed.var.metric = 'mean',
                                       col = c("blue", "red"),
                                       legend = F,
                                       data_species = get_formal_data(myBiomodModelOut.eut.edu,'resp.var'),
                                       plot=F)

save.image()

bio_14_response_eut.edu<-gather(myRespPlot2D.eut.edu$bio14,model,Probability,
                               eut.edu_PA1_RUN1_MAXENT.Phillips:eut.edu_PA1_RUN5_MAXENT.Phillips)



p_eut.edu_bio14<-ggplot(data =bio_14_response_eut.edu) + 
  geom_smooth(mapping = aes(x = bio14, y = Probability))+
  labs(x = bquote('Bio 14'))
p_eut.edu_bio14



###Ensemble modeling
myBiomodEM.eut.edu<-BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut.eut.edu,
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

myBiomodProjection.eut.edu <- BIOMOD_Projection(modeling.output = myBiomodModelOut.eut.edu,
                                                new.env = stack(ma.presente.poly.eut.edu),
                                                proj.name = 'current_MA',
                                                selected.models = "all",
                                                binary.meth = 'TSS',
                                                compress = FALSE,
                                                build.clamping.mask = FALSE)

#Ensemble Forcasting

myBiomodEF.eut.edu.presente <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM.eut.edu,
  projection.output = myBiomodProjection.eut.edu,binary.meth = "TSS")

save.image()

plot(myBiomodEF.eut.edu.presente)

currentPred.ensemble.eut.edu <- stack("eut.edu/proj_current_MA/proj_current_MA_eut.edu_ensemble.grd")
plot(currentPred.ensemble.eut.edu[[1]])

plot(currentPred.ensemble.eut.edu[[1]],col = gray.colors(12, start = 0.9, end = 0.0))
plot(wrld_simpl, xlim=c(-60,-40), ylim=c(-34,-15), axes=TRUE, add=T)
polygon(poligono_mata_atlantica)
points(coord.eut.edu.names_thinned_dataset_full[[1]], cex=.6,col="black")

plot(currentPred.ensemble.eut.edu[[2]],col = gray.colors(12, start = 0.9, end = 0.0))
points(coord.eut.edu.names_thinned_dataset_full[[1]], cex=.6,col="black")
polygon(poligono_mata_atlantica)



citations.eut.edu<-BIEN_metadata_citation(dataframe=eut.edu.bien)#If you are referencing occurrence data
citations.eut.edu


######################### 2070 ##########################
################################################################

ma.2070.poly.eut.edu<-stack(ma.2070.new[[2]],
                                ma.2070.new[[3]],
                                ma.2070.new[[8]],
                                ma.2070.new[[9]],
                                ma.2070.new[[10]],
                                ma.2070.new[[13]],
                                ma.2070.new[[14]],
                                ma.2070.new[[15]],
                                ma.2070.new[[18]],
                                ma.2070.new[[19]])





names(ma.2070.poly.eut.edu)<-names(ma.presente.poly.eut.edu)


myExplFuture.2070.85.ma<-stack(ma.2070.poly.eut.edu) #modelo HE


projection.2070.85.ma.eut.edu <-BIOMOD_Projection(modeling.output=myBiomodModelOut.eut.edu,
                                                  new.env= myExplFuture.2070.85.ma,
                                                  proj.name='eut.edu.2070.85',
                                                  selected.models =  "all",
                                                  Bin.trans=TRUE,
                                                  binary.meth ='TSS',
                                                  compress = 'gzip',
                                                  build.clamping.mask	= FALSE,
                                                  SaveObj=TRUE,
                                                  keep.in.memory=F,
                                                  do.stack=T,
                                                  silent=F)


eut.edu.2070.85<-get_predictions(projection.2070.85.ma.eut.edu)
eut.edu.2070.85
summary(eut.edu.2070.85)





################## ENSAMBLE MODELLING 2070  ################
#################################################################



EnsambleForecast.2070.85.ma.eut.edu <-BIOMOD_EnsembleForecasting (
  projection.output = projection.2070.85.ma.eut.edu,
  EM.output=myBiomodEM.eut.edu,
  selected.models = "all",
  binary.meth ='TSS') 




EnsambleForecast.2070.85.ma.eut.edu
plot(EnsambleForecast.2070.85.ma.eut.edu)



ensemble.2070.HE.85.ma.eut.edu <- stack("eut.edu/proj_eut.edu.2070.85/proj_eut.edu.2070.85_eut.edu_ensemble.grd")
ensemble.2070.HE.85.ma.eut.edu



plot(ensemble.2070.HE.85.ma.eut.edu[[1]],col = gray.colors(12, start = 0.9, end = 0))
points(coord.eut.edu.names_thinned_dataset_full[[1]],pch=20, col="black", cex=.7)

plot(ensemble.2070.HE.85.ma.eut.edu[[2]],col = gray.colors(12, start = 0.9, end = 0))



futurePred.ensemble.2070.85.eut.edu <- stack("eut.edu/proj_eut.edu.2070.85/proj_eut.edu.2070.85_eut.edu_ensemble.grd")



plot(currentPred.ensemble.eut.edu[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.2070.85.eut.edu[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS


currentPred.ensemble.bin.eut.edu<- stack("eut.edu/proj_current_MA/proj_current_MA_eut.edu_ensemble_TSSbin.grd")
futurePred.ensemble.bin.2070.85.eut.edu <- stack("eut.edu/proj_eut.edu.2070.85/proj_eut.edu.2070.85_eut.edu_ensemble_TSSbin.gri")


plot(currentPred.ensemble.bin.eut.edu[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.bin.2070.85.eut.edu[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS



myBiomodRangeSize.2070.85.eut.edu <- BIOMOD_RangeSize(
  CurrentPred=currentPred.ensemble.bin.eut.edu[[1]],
  FutureProj=futurePred.ensemble.bin.2070.85.eut.edu[[1]])
myBiomodRangeSize.2070.85.eut.edu$Compt.By.Models





plot(myBiomodRangeSize.2070.85.eut.edu$Diff.By.Pixel, col= gray.colors(4, start = 0, end = 1, gamma = 1, alpha = NULL))
plot(wrld_simpl,  axes="TRUE", add=T)

polygon(poligono_mata_atlantica)








save.image()









