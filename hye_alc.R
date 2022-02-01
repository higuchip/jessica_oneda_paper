####ADEQUACAO

library(biomod2)
library(splancs)
library(SDMTools)
library(dplyr)
library(tidyr)
library(gridExtra) 
library(raster)
library(spThin)

BioC<- getData('worldclim', var='bio', res=10) # Bioclima
variaveis <- stack(BioC)

plot(wrld_simpl, xlim=c(-130,-35), ylim=c(-45,35), axes=TRUE, col="gray")
plot(wrld_simpl, xlim=c(-70,-45), ylim=c(-35,-5), axes=TRUE, col="gray")
points(coord.hye.alc, pch=19, cex=.51, col="black")


##Poligono MA

require(rgdal)
data.shape<-readOGR(dsn="/media/dendrologia/423EF1863EF172F1/User/Documents/orientados/Guilherme/IDS_09_Desm_Bioma_Mata_Atlantica_1_2", 
                    layer="IDS_09_Desm_Bioma_Mata_Atlantica_1_2")

poligono_mata_atlantica<-data.shape@polygons[[1]]@Polygons[[1]]@coords

polygon(poligono_mata_atlantica)
###Poligono especie


poligono_hye.alc<-getpoly()

colnames(poligono_hye.alc) = c("Long", "Lat")
polygon(poligono_hye.alc) 



########REMOCAO DUVIDAS

out.hye.alc= pnt.in.poly(coord.hye.alc,poligono_hye.alc)
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
plot(data.shape, add=T, col="lightgreen")
points(out.hye.alc[which(out.hye.alc$pip==1),1:2],pch=20, cex=.5, col="black")
coord.hye.alc.limpo<-out.hye.alc[which(out.hye.alc$pip==1),1:2]
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
points(coord.hye.alc.limpo,pch=20, cex=.5, col="black")


###FILTRAGEM ESPACIAl

head(coord.hye.alc.limpo)
coord.hye.alc.limpo<-na.omit(coord.hye.alc.limpo)

dim(coord.hye.alc.limpo)
hye.alc_names<-c(rep("hye.alc",  dim(coord.hye.alc.limpo)[1]))
coord.hye.alc.limpo.names<-cbind(hye.alc_names, coord.hye.alc.limpo)

head(coord.hye.alc.limpo.names)
is.data.frame(coord.hye.alc.limpo.names)
colnames(coord.hye.alc.limpo.names)<-c("SPEC", "LONG", "LAT")

coord.hye.alc.names_thinned_dataset_full <-
  thin( loc.data =coord.hye.alc.limpo.names, 
        lat.col = "LAT", 
        long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 26.27, reps =1, #10 min resolution
        locs.thinned.list.return = TRUE, 
        write.files = T, 
        max.files = 5, 
        out.dir = "coord.hye.alc.limpo.names_full/", out.base = "coord.hye.alc.limpo.names_thinned", 
        write.log.file = TRUE,
        log.file = "coord.hye.alc.limpo.names_thinned_full_log_file.txt" )




dim(coord.hye.alc.names_thinned_dataset_full[[1]])



###ENVIRONMENTAL BACKGROUND

reg.variaveis.hye.alc= crop(variaveis, poligono_hye.alc)
plot(reg.variaveis.hye.alc, 1)
points(coord.hye.alc.limpo, pch=20, cex=.5)

#remoção variaveis nao colineares
vari.n.colineares.hye.alc<-vifstep(reg.variaveis.hye.alc)
vari.n.colineares.hye.alc


reg.predictors.hye.alc<- stack(reg.variaveis.hye.alc[[2]],
                               reg.variaveis.hye.alc[[3]],
                               reg.variaveis.hye.alc[[8]],
                               reg.variaveis.hye.alc[[9]],
                               reg.variaveis.hye.alc[[13]],
                               reg.variaveis.hye.alc[[14]],
                               reg.variaveis.hye.alc[[15]],
                               reg.variaveis.hye.alc[[18]],
                               reg.variaveis.hye.alc[[19]])



names(reg.predictors.hye.alc)


###MATA ATLANTICA projection presente

ma.presente = crop(variaveis, poligono_mata_atlantica)
plot(ma.presente, 1)
polygon(poligono_mata_atlantica)
ma.presente.poly <- mask(ma.presente, data.shape)
plot(ma.presente.poly, 1 )


ma.presente.poly.hye.alc<-stack(ma.presente.new[[2]],
                                ma.presente.new[[3]],
                                ma.presente.new[[8]],
                                ma.presente.new[[9]],
                                ma.presente.new[[13]],
                                ma.presente.new[[14]],
                                ma.presente.new[[15]],
                                ma.presente.new[[18]],
                                ma.presente.new[[19]])




names(ma.presente.poly.hye.alc)


myRespName.hye.alc <-  'hye.alc'
myResp.hye.alc<-rep(1, nrow(coord.hye.alc.names_thinned_dataset_full[[1]]))
myResp.hye.alc
length(myResp.hye.alc)
myRespXY.hye.alc<-coord.hye.alc.names_thinned_dataset_full[[1]]

# function to define the intersect of rasters
intersect_mask <- function(x){
  values_x <- getValues(x)
  inter_x <- values_x %*% rep(1,nlayers(x))
  mask <- setValues(subset(x,1),values = (inter_x>0))
  return(mask)
}
# keep only all cells that are defined for all layers

myExpl.hye.alc<-stack(reg.predictors.hye.alc)

myExpl.hye.alc <- stack(mask(myExpl.hye.alc, intersect_mask(myExpl.hye.alc)))


##Number of pseudo-absence: Lobo & Tognelli sugerem 100 pontos de pseudo-
#ausencia para cada ponto de presenca.



myBiomodData.hye.alc<- BIOMOD_FormatingData(resp.var = myResp.hye.alc,
                                            expl.var = myExpl.hye.alc,
                                            resp.xy = myRespXY.hye.alc,
                                            resp.name = myRespName.hye.alc,
                                            PA.nb.rep = 1,
                                            PA.nb.absences = 10000,
                                            PA.strategy = 'disk',
                                            PA.dist.max=500000) # com base em http://www.sciencedirect.com/science/article/pii/S0304380008005486
myBiomodData.hye.alc
save.image()
###MODELING

# 1. Defining Models Options using default options.
myBiomodOption.hye.alc <- BIOMOD_ModelingOptions()

# 3. Computing the models

myBiomodModelOut.hye.alc<- BIOMOD_Modeling(myBiomodData.hye.alc,
                                           models = c("MAXENT.Phillips"),
                                           models.options = myBiomodOption.hye.alc,
                                           NbRunEval=5, #conforme http://onlinelibrary.wiley.com/doi/10.1111/j.1523-1739.2006.00354.x/full
                                           DataSplit=70,
                                           Prevalence=0.5,
                                           VarImport=1,
                                           models.eval.meth = c('TSS'),
                                           hye.alcveObj = TRUE,
                                           rescal.all.models = TRUE,
                                           do.full.models = FALSE,
                                           modeling.id = paste(myRespName.hye.alc,"FirstModeling",sep=""))
(myBiomodModelOut.hye.alc)
save.image()
# get all models evaluation

#Excelente TSS>0.75
#Bom 0.40>=TSS<0.75
#Ruim TSS<0.40


myBiomodModelEval.hye.alc <- get_evaluations(myBiomodModelOut.hye.alc)

myBiomodModelEval.hye.alc




get_variables_importance(myBiomodModelOut.hye.alc)
get_evaluations_matrix_hye.alc<-cbind(get_variables_importance(myBiomodModelOut.hye.alc)[1:9], #RUN 1
                                      get_variables_importance(myBiomodModelOut.hye.alc)[10:18], #RUN 2
                                      get_variables_importance(myBiomodModelOut.hye.alc)[19:27], #RUN 3
                                      get_variables_importance(myBiomodModelOut.hye.alc)[28:36], #RUN 4
                                      get_variables_importance(myBiomodModelOut.hye.alc)[37:45]) #RUN 5
colnames(get_evaluations_matrix_hye.alc)<-c("RUN1", "RUN2", "RUN3", "RUN4", "RUN5")
names(reg.predictors.hye.alc)
rownames(get_evaluations_matrix_hye.alc)<-names(reg.predictors.hye.alc)
get_evaluations_matrix_hye.alc_mean<-apply(get_evaluations_matrix_hye.alc, 1, mean)
which.max(get_evaluations_matrix_hye.alc_mean)



mymaxent.hye.alc <- BIOMOD_LoadModels(myBiomodModelOut.hye.alc, models='MAXENT.Phillips')
myRespPlot2D.hye.alc <- response.plot2(models  = mymaxent.hye.alc,
                                       Data = get_formal_data(myBiomodModelOut.hye.alc,'expl.var'), 
                                       show.variables= get_formal_data(myBiomodModelOut.hye.alc,'expl.var.names'),
                                       do.bivariate = FALSE,
                                       fixed.var.metric = 'mean',
                                       col = c("blue", "red"),
                                       legend = F,
                                       data_species = get_formal_data(myBiomodModelOut.hye.alc,'resp.var'),
                                       plot=F)

save.image()

bio_9_response_hye.alc<-gather(myRespPlot2D.hye.alc$bio9,model,Probability,
                                hye.alc_PA1_RUN1_MAXENT.Phillips:hye.alc_PA1_RUN5_MAXENT.Phillips)



p_hye.alc_bio9<-ggplot(data =bio_9_response_hye.alc) + 
  geom_smooth(mapping = aes(x = bio9/10, y = Probability))+
  labs(x = bquote('Bio 9'))
p_hye.alc_bio9



###Ensemble modeling
myBiomodEM.hye.alc<-BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut.hye.alc,
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

myBiomodProjection.hye.alc <- BIOMOD_Projection(modeling.output = myBiomodModelOut.hye.alc,
                                                new.env = stack(ma.presente.poly.hye.alc),
                                                proj.name = 'current_MA',
                                                selected.models = "all",
                                                binary.meth = 'TSS',
                                                compress = FALSE,
                                                build.clamping.mask = FALSE)

#Ensemble Forcasting

myBiomodEF.hye.alc.presente <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM.hye.alc,
  projection.output = myBiomodProjection.hye.alc,binary.meth = "TSS")

save.image()

plot(myBiomodEF.hye.alc.presente)

currentPred.ensemble.hye.alc <- stack("hye.alc/proj_current_MA/proj_current_MA_hye.alc_ensemble.grd")
plot(currentPred.ensemble.hye.alc[[1]])

plot(currentPred.ensemble.hye.alc[[1]],col = gray.colors(12, start = 0.9, end = 0.0))
plot(wrld_simpl, xlim=c(-60,-40), ylim=c(-34,-15), axes=TRUE, add=T)
polygon(poligono_mata_atlantica)


plot(currentPred.ensemble.hye.alc[[2]],col = gray.colors(12, start = 0.9, end = 0.0))
points(coord.hye.alc.names_thinned_dataset_full[[1]], cex=1,col="black", pch=20)
polygon(poligono_mata_atlantica)


########################## 2070 ##########################
################################################################



ma.2070.poly.hye.alc<-stack(ma.2070.new[[2]],
                                ma.2070.new[[3]],
                                ma.2070.new[[8]],
                                ma.2070.new[[9]],
                                ma.2070.new[[13]],
                                ma.2070.new[[14]],
                                ma.2070.new[[15]],
                                ma.2070.new[[18]],
                                ma.2070.new[[19]])



names(ma.2070.poly.hye.alc)<-names(ma.presente.poly.hye.alc)


myExplFuture.2070.85.ma<-stack(ma.2070.poly.hye.alc) #modelo HE


projection.2070.85.ma.hye.alc <-BIOMOD_Projection(modeling.output=myBiomodModelOut.hye.alc,
                                                  new.env= myExplFuture.2070.85.ma,
                                                  proj.name='hye.alc.2070.85',
                                                  selected.models =  "all",
                                                  Bin.trans=TRUE,
                                                  binary.meth ='TSS',
                                                  compress = 'gzip',
                                                  build.clamping.mask	= FALSE,
                                                  SaveObj=TRUE,
                                                  keep.in.memory=F,
                                                  do.stack=T,
                                                  silent=F)


hye.alc.2070.85<-get_predictions(projection.2070.85.ma.hye.alc)
hye.alc.2070.85
summary(hye.alc.2070.85)





################## ENSAMBLE MODELLING 2070  ################
#################################################################



EnsambleForecast.2070.85.ma.hye.alc <-BIOMOD_EnsembleForecasting (
  projection.output = projection.2070.85.ma.hye.alc,
  EM.output=myBiomodEM.hye.alc,
  selected.models = "all",
  binary.meth ='TSS') 




EnsambleForecast.2070.85.ma.hye.alc
plot(EnsambleForecast.2070.85.ma.hye.alc)



ensemble.2070.HE.85.ma.hye.alc <- stack("hye.alc/proj_hye.alc.2070.85/proj_hye.alc.2070.85_hye.alc_ensemble.grd")
ensemble.2070.HE.85.ma.hye.alc



plot(ensemble.2070.HE.85.ma.hye.alc[[1]],col = gray.colors(12, start = 0.9, end = 0))
points(coord.gym_klo.names_thinned_dataset_full[[1]],pch=20, col="black", cex=.7)

plot(ensemble.2070.HE.85.ma.hye.alc[[2]],col = gray.colors(12, start = 0.9, end = 0))



futurePred.ensemble.2070.85.hye.alc <- stack("hye.alc/proj_hye.alc.2070.85/proj_hye.alc.2070.85_hye.alc_ensemble.grd")



plot(currentPred.ensemble.hye.alc[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.2070.85.hye.alc[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS



currentPred.ensemble.bin.hye.alc<- stack("hye.alc/proj_current_MA/proj_current_MA_hye.alc_ensemble_TSSbin.grd")
futurePred.ensemble.bin.2070.85.hye.alc <- stack("hye.alc/proj_hye.alc.2070.85/proj_hye.alc.2070.85_hye.alc_ensemble_TSSbin.gri")


plot(currentPred.ensemble.bin.hye.alc[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.bin.2070.85.hye.alc[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS



myBiomodRangeSize.2070.85.hye.alc <- BIOMOD_RangeSize(
  CurrentPred=currentPred.ensemble.bin.hye.alc[[1]],
  FutureProj=futurePred.ensemble.bin.2070.85.hye.alc[[1]])
myBiomodRangeSize.2070.85.hye.alc$Compt.By.Models





plot(myBiomodRangeSize.2070.85.hye.alc$Diff.By.Pixel, col= gray.colors(4, start = 0, end = 1, gamma = 1, alpha = NULL))
plot(wrld_simpl,  axes="TRUE", add=T)

polygon(poligono_mata_atlantica)








