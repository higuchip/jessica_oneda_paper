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
points(coord.gua.aus, pch=19, cex=.51, col="black")


##Poligono MA

require(rgdal)
data.shape<-readOGR(dsn="/media/dendrologia/423EF1863EF172F1/User/Documents/orientados/Guilherme/IDS_09_Desm_Bioma_Mata_Atlantica_1_2", 
                    layer="IDS_09_Desm_Bioma_Mata_Atlantica_1_2")

poligono_mata_atlantica<-data.shape@polygons[[1]]@Polygons[[1]]@coords

polygon(poligono_mata_atlantica)
###Poligono especie


poligono_gua.aus<-getpoly()

colnames(poligono_gua.aus) = c("Long", "Lat")
polygon(poligono_gua.aus) 



########REMOCAO DUVIDAS

out.gua.aus= pnt.in.poly(coord.gua.aus,poligono_gua.aus)
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
plot(data.shape, add=T, col="lightgreen")
points(out.gua.aus[which(out.gua.aus$pip==1),1:2],pch=20, cex=.5, col="black")
coord.gua.aus.limpo<-out.gua.aus[which(out.gua.aus$pip==1),1:2]
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
points(coord.gua.aus.limpo,pch=20, cex=.5, col="black")


###FILTRAGEM ESPACIAl

head(coord.gua.aus.limpo)
coord.gua.aus.limpo<-na.omit(coord.gua.aus.limpo)

dim(coord.gua.aus.limpo)
gua.aus_names<-c(rep("gua.aus",  dim(coord.gua.aus.limpo)[1]))
coord.gua.aus.limpo.names<-cbind(gua.aus_names, coord.gua.aus.limpo)

head(coord.gua.aus.limpo.names)
is.data.frame(coord.gua.aus.limpo.names)
colnames(coord.gua.aus.limpo.names)<-c("SPEC", "LONG", "LAT")
# 
# ?write.table
# write.table(coord.gua.aus.limpo.names, "tap_gui.csv", sep=";", dec=",")

coord.gua.aus.names_thinned_dataset_full <-
  thin( loc.data =coord.gua.aus.limpo.names, 
        lat.col = "LAT", 
        long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 26.27, reps =1, #10 min resolution
        locs.thinned.list.return = TRUE, 
        write.files = T, 
        max.files = 5, 
        out.dir = "coord.gua.aus.limpo.names_full/", out.base = "coord.gua.aus.limpo.names_thinned", 
        write.log.file = TRUE,
        log.file = "coord.gua.aus.limpo.names_thinned_full_log_file.txt" )



dim(coord.gua.aus.limpo.names)
dim(coord.gua.aus.names_thinned_dataset_full[[1]])



###ENVIRONMENTAL BACKGROUND

reg.variaveis.gua.aus= crop(variaveis, poligono_gua.aus)
plot(reg.variaveis.gua.aus, 1)
points(coord.gua.aus.limpo, pch=20, cex=.5)

#remoção variaveis nao colineares
vari.n.colineares.gua.aus<-vifstep(reg.variaveis.gua.aus)
vari.n.colineares.gua.aus


reg.predictors.gua.aus<- stack(reg.variaveis.gua.aus[[2]],
                               reg.variaveis.gua.aus[[3]],
                               reg.variaveis.gua.aus[[8]],
                               reg.variaveis.gua.aus[[9]],
                               reg.variaveis.gua.aus[[13]],
                               reg.variaveis.gua.aus[[14]],
                               reg.variaveis.gua.aus[[15]],
                               reg.variaveis.gua.aus[[18]],
                               reg.variaveis.gua.aus[[19]])



names(reg.predictors.gua.aus)


###MATA ATLANTICA projection presente

ma.presente = crop(variaveis, poligono_mata_atlantica)
plot(ma.presente, 1)
polygon(poligono_mata_atlantica)
ma.presente.poly <- mask(ma.presente, data.shape)
plot(ma.presente.poly, 1 )


ma.presente.new.gua.aus<-stack(ma.presente.new[[2]],
                                ma.presente.new[[3]],
                                ma.presente.new[[8]],
                                ma.presente.new[[9]],
                                ma.presente.new[[13]],
                                ma.presente.new[[14]],
                                ma.presente.new[[15]],
                                ma.presente.new[[18]],
                                ma.presente.new[[19]])


plot(ma.presente.new.gua.aus,2)

names(ma.presente.poly.gua.aus)


myRespName.gua.aus <-  'gua.aus'
myResp.gua.aus<-rep(1, nrow(coord.gua.aus.names_thinned_dataset_full[[1]]))
myResp.gua.aus
length(myResp.gua.aus)
myRespXY.gua.aus<-coord.gua.aus.names_thinned_dataset_full[[1]]

# function to define the intersect of rasters
intersect_mask <- function(x){
  values_x <- getValues(x)
  inter_x <- values_x %*% rep(1,nlayers(x))
  mask <- setValues(subset(x,1),values = (inter_x>0))
  return(mask)
}
# keep only all cells that are defined for all layers

myExpl.gua.aus<-stack(reg.predictors.gua.aus)

myExpl.gua.aus <- stack(mask(myExpl.gua.aus, intersect_mask(myExpl.gua.aus)))


##Number of pseudo-absence: Lobo & Tognelli sugerem 100 pontos de pseudo-
#ausencia para cada ponto de presenca.



myBiomodData.gua.aus<- BIOMOD_FormatingData(resp.var = myResp.gua.aus,
                                            expl.var = myExpl.gua.aus,
                                            resp.xy = myRespXY.gua.aus,
                                            resp.name = myRespName.gua.aus,
                                            PA.nb.rep = 1,
                                            PA.nb.absences = 10000,
                                            PA.strategy = 'disk',
                                            PA.dist.max=500000) # com base em http://www.sciencedirect.com/science/article/pii/S0304380008005486
myBiomodData.gua.aus
save.image()
###MODELING

# 1. Defining Models Options using default options.
myBiomodOption.gua.aus <- BIOMOD_ModelingOptions()

# 3. Computing the models

myBiomodModelOut.gua.aus<- BIOMOD_Modeling(myBiomodData.gua.aus,
                                           models = c("MAXENT.Phillips"),
                                           models.options = myBiomodOption.gua.aus,
                                           NbRunEval=5, #conforme http://onlinelibrary.wiley.com/doi/10.1111/j.1523-1739.2006.00354.x/full
                                           DataSplit=70,
                                           Prevalence=0.5,
                                           VarImport=1,
                                           models.eval.meth = c('TSS'),
                                           gua.ausveObj = TRUE,
                                           rescal.all.models = TRUE,
                                           do.full.models = FALSE,
                                           modeling.id = paste(myRespName.gua.aus,"FirstModeling",sep=""))
(myBiomodModelOut.gua.aus)
save.image()
# get all models evaluation

#Excelente TSS>0.75
#Bom 0.40>=TSS<0.75
#Ruim TSS<0.40


myBiomodModelEval.gua.aus <- get_evaluations(myBiomodModelOut.gua.aus)

myBiomodModelEval.gua.aus





get_variables_importance(myBiomodModelOut.gua.aus)
get_evaluations_matrix_gua.aus<-cbind(get_variables_importance(myBiomodModelOut.gua.aus)[1:9], #RUN 1
                                      get_variables_importance(myBiomodModelOut.gua.aus)[10:18], #RUN 2
                                      get_variables_importance(myBiomodModelOut.gua.aus)[19:27], #RUN 3
                                      get_variables_importance(myBiomodModelOut.gua.aus)[28:36], #RUN 4
                                      get_variables_importance(myBiomodModelOut.gua.aus)[37:45]) #RUN 5
colnames(get_evaluations_matrix_gua.aus)<-c("RUN1", "RUN2", "RUN3", "RUN4", "RUN5")
names(reg.predictors.gua.aus)
rownames(get_evaluations_matrix_gua.aus)<-names(reg.predictors.gua.aus)
get_evaluations_matrix_gua.aus_mean<-apply(get_evaluations_matrix_gua.aus, 1, mean)
which.max(get_evaluations_matrix_gua.aus_mean)



mymaxent.gua.aus <- BIOMOD_LoadModels(myBiomodModelOut.gua.aus, models='MAXENT.Phillips')
myRespPlot2D.gua.aus <- response.plot2(models  = mymaxent.gua.aus,
                                       Data = get_formal_data(myBiomodModelOut.gua.aus,'expl.var'), 
                                       show.variables= get_formal_data(myBiomodModelOut.gua.aus,'expl.var.names'),
                                       do.bivariate = FALSE,
                                       fixed.var.metric = 'mean',
                                       col = c("blue", "red"),
                                       legend = F,
                                       data_species = get_formal_data(myBiomodModelOut.gua.aus,'resp.var'),
                                       plot=F)

save.image()

bio_14_response_gua.aus<-gather(myRespPlot2D.gua.aus$bio14,model,Probability,
                               gua.aus_PA1_RUN1_MAXENT.Phillips:gua.aus_PA1_RUN5_MAXENT.Phillips)



p_gua.aus_bio14<-ggplot(data =bio_14_response_gua.aus) + 
  geom_smooth(mapping = aes(x = bio14, y = Probability))+
  labs(x = bquote('Bio 14'))
p_gua.aus_bio14



###Ensemble modeling
myBiomodEM.gua.aus<-BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut.gua.aus,
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

myBiomodProjection.gua.aus <- BIOMOD_Projection(modeling.output = myBiomodModelOut.gua.aus,
                                                new.env = stack(ma.presente.new.gua.aus),
                                                proj.name = 'current_MA',
                                                selected.models = "all",
                                                binary.meth = 'TSS',
                                                compress = FALSE,
                                                build.clamping.mask = FALSE)

plot(myBiomodProjection.gua.aus )

#Ensemble Forcasting

myBiomodEF.gua.aus.presente <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM.gua.aus,
  projection.output = myBiomodProjection.gua.aus,binary.meth = "TSS")

save.image()

plot(myBiomodEF.gua.aus.presente)

currentPred.ensemble.gua.aus <- stack("gua.aus/proj_current_MA/proj_current_MA_gua.aus_ensemble.grd")
plot(currentPred.ensemble.gua.aus[[1]])

plot(currentPred.ensemble.gua.aus[[1]],col = gray.colors(12, start = 0.9, end = 0.0))
plot(wrld_simpl, xlim=c(-60,-40), ylim=c(-34,-15), axes=TRUE, add=T)
polygon(poligono_mata_atlantica)
points(coord.gua.aus.names_thinned_dataset_full[[1]], cex=.6,col="black")

plot(currentPred.ensemble.gua.aus[[2]],col = gray.colors(12, start = 0.9, end = 0.0))
points(coord.gua.aus.names_thinned_dataset_full[[1]], cex=.6,col="black")
polygon(poligono_mata_atlantica)


?BIEN_metadata_citation

citations.gua.aus<-BIEN_metadata_citation(dataframe=gua.aus.bien)#If you are referencing occurrence data
citations.gua.aus


########################## 2070 ##########################
################################################################



ma.2070.poly.gua.aus<-stack(ma.2070.new[[2]],
                                ma.2070.new[[3]],
                                ma.2070.new[[8]],
                                ma.2070.new[[9]],
                                ma.2070.new[[13]],
                                ma.2070.new[[14]],
                                ma.2070.new[[15]],
                                ma.2070.new[[18]],
                                ma.2070.new[[19]])



names(ma.2070.poly.gua.aus)<-names(ma.presente.poly.gua.aus)


myExplFuture.2070.85.ma<-stack(ma.2070.poly.gua.aus) #modelo HE


projection.2070.85.ma.gua.aus <-BIOMOD_Projection(modeling.output=myBiomodModelOut.gua.aus,
                                                  new.env= myExplFuture.2070.85.ma,
                                                  proj.name='gua.aus.2070.85',
                                                  selected.models =  "all",
                                                  Bin.trans=TRUE,
                                                  binary.meth ='TSS',
                                                  compress = 'gzip',
                                                  build.clamping.mask	= FALSE,
                                                  SaveObj=TRUE,
                                                  keep.in.memory=F,
                                                  do.stack=T,
                                                  silent=F)


gua.aus.2070.85<-get_predictions(projection.2070.85.ma.gua.aus)
gua.aus.2070.85
summary(gua.aus.2070.85)





################## ENSAMBLE MODELLING 2070  ################
#################################################################



EnsambleForecast.2070.85.ma.gua.aus <-BIOMOD_EnsembleForecasting (
  projection.output = projection.2070.85.ma.gua.aus,
  EM.output=myBiomodEM.gua.aus,
  selected.models = "all",
  binary.meth ='TSS') 




EnsambleForecast.2070.85.ma.gua.aus
plot(EnsambleForecast.2070.85.ma.gua.aus)



ensemble.2070.HE.85.ma.gua.aus <- stack("gua.aus/proj_gua.aus.2070.85/proj_gua.aus.2070.85_gua.aus_ensemble.grd")
ensemble.2070.HE.85.ma.gua.aus



plot(ensemble.2070.HE.85.ma.gua.aus[[1]],col = gray.colors(12, start = 0.9, end = 0))
points(coord.gym_klo.names_thinned_dataset_full[[1]],pch=20, col="black", cex=.7)

plot(ensemble.2070.HE.85.ma.gua.aus[[2]],col = gray.colors(12, start = 0.9, end = 0))



futurePred.ensemble.2070.85.gua.aus <- stack("gua.aus/proj_gua.aus.2070.85/proj_gua.aus.2070.85_gua.aus_ensemble.grd")



plot(currentPred.ensemble.gua.aus[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.2070.85.gua.aus[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS


currentPred.ensemble.bin.gua.aus<- stack("gua.aus/proj_current_MA/proj_current_MA_gua.aus_ensemble_TSSbin.grd")
futurePred.ensemble.bin.2070.85.gua.aus <- stack("gua.aus/proj_gua.aus.2070.85/proj_gua.aus.2070.85_gua.aus_ensemble_TSSbin.gri")


plot(currentPred.ensemble.bin.gua.aus[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.bin.2070.85.gua.aus[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS



myBiomodRangeSize.2070.85.gua.aus <- BIOMOD_RangeSize(
  CurrentPred=currentPred.ensemble.bin.gua.aus[[1]],
  FutureProj=futurePred.ensemble.bin.2070.85.gua.aus[[1]])
myBiomodRangeSize.2070.85.gua.aus$Compt.By.Models





plot(myBiomodRangeSize.2070.85.gua.aus$Diff.By.Pixel, col= gray.colors(4, start = 0, end = 1, gamma = 1, alpha = NULL))
plot(wrld_simpl,  axes="TRUE", add=T)

polygon(poligono_mata_atlantica)








save.image()

