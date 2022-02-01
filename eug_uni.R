####ADEQUACAO

library(biomod2)
library(splancs)
library(SDMTools)
library(dplyr)
library(tidyr)
library(gridExtra) 
library(raster)
library(spThin)


plot(wrld_simpl, xlim=c(-130,-35), ylim=c(-45,35), axes=TRUE, col="gray")
plot(wrld_simpl, xlim=c(-70,-45), ylim=c(-35,-5), axes=TRUE, col="gray")
points(coord.eug.uni, pch=19, cex=.51, col="black")


##Poligono MA

require(rgdal)
data.shape<-readOGR(dsn="/media/dendrologia/423EF1863EF172F1/User/Documents/orientados/Guilherme/IDS_09_Desm_Bioma_Mata_Atlantica_1_2", 
                    layer="IDS_09_Desm_Bioma_Mata_Atlantica_1_2")

poligono_mata_atlantica<-data.shape@polygons[[1]]@Polygons[[1]]@coords

polygon(poligono_mata_atlantica)
###Poligono especie


poligono_eug.uni<-getpoly()

colnames(poligono_eug.uni) = c("Long", "Lat")
polygon(poligono_eug.uni) 



########REMOCAO DUVIDAS

out.eug.uni= pnt.in.poly(coord.eug.uni,poligono_eug.uni)
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
plot(data.shape, add=T, col="lightgreen")
points(out.eug.uni[which(out.eug.uni$pip==1),1:2],pch=20, cex=.5, col="black")
coord.eug.uni.limpo<-out.eug.uni[which(out.eug.uni$pip==1),1:2]
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
points(coord.eug.uni.limpo,pch=20, cex=.5, col="black")


###FILTRAGEM ESPACIAl

head(coord.eug.uni.limpo)
coord.eug.uni.limpo<-na.omit(coord.eug.uni.limpo)

dim(coord.eug.uni.limpo)
eug.uni_names<-c(rep("eug.uni",  dim(coord.eug.uni.limpo)[1]))
coord.eug.uni.limpo.names<-cbind(eug.uni_names, coord.eug.uni.limpo)

head(coord.eug.uni.limpo.names)
is.data.frame(coord.eug.uni.limpo.names)
colnames(coord.eug.uni.limpo.names)<-c("SPEC", "LONG", "LAT")

coord.eug.uni.names_thinned_dataset_full <-
  thin( loc.data =coord.eug.uni.limpo.names, 
        lat.col = "LAT", 
        long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 26.27, reps =1, #10 min resolution
        locs.thinned.list.return = TRUE, 
        write.files = T, 
        max.files = 5, 
        out.dir = "coord.eug.uni.limpo.names_full/", out.base = "coord.eug.uni.limpo.names_thinned", 
        write.log.file = TRUE,
        log.file = "coord.eug.uni.limpo.names_thinned_full_log_file.txt" )




dim(coord.eug.uni.names_thinned_dataset_full[[1]])



###ENVIRONMENTAL BACKGROUND

reg.variaveis.eug.uni= crop(variaveis, poligono_eug.uni)
plot(reg.variaveis.eug.uni, 1)
points(coord.eug.uni.limpo, pch=20, cex=.5)

#remoção variaveis nao colineares
vari.n.colineares.eug.uni<-vifstep(reg.variaveis.eug.uni)
vari.n.colineares.eug.uni


reg.predictors.eug.uni<- stack(reg.variaveis.eug.uni[[2]],
                               reg.variaveis.eug.uni[[3]],
                               reg.variaveis.eug.uni[[8]],
                               reg.variaveis.eug.uni[[9]],
                               reg.variaveis.eug.uni[[13]],
                               reg.variaveis.eug.uni[[14]],
                               reg.variaveis.eug.uni[[15]],
                               reg.variaveis.eug.uni[[18]],
                               reg.variaveis.eug.uni[[19]])



names(reg.predictors.eug.uni)


###MATA ATLANTICA projection presente

ma.presente = crop(variaveis, poligono_mata_atlantica)
plot(ma.presente, 1)
polygon(poligono_mata_atlantica)
ma.presente.poly <- mask(ma.presente, data.shape)
plot(ma.presente.poly, 1 )


ma.presente.poly.eug.uni<-stack(ma.presente.new[[2]],
                                ma.presente.new[[3]],
                                ma.presente.new[[8]],
                                ma.presente.new[[9]],
                                ma.presente.new[[13]],
                                ma.presente.new[[14]],
                                ma.presente.new[[15]],
                                ma.presente.new[[18]],
                                ma.presente.new[[19]])




names(ma.presente.poly.eug.uni)


myRespName.eug.uni <-  'eug.uni'
myResp.eug.uni<-rep(1, nrow(coord.eug.uni.names_thinned_dataset_full[[1]]))
myResp.eug.uni
length(myResp.eug.uni)
myRespXY.eug.uni<-coord.eug.uni.names_thinned_dataset_full[[1]]

# function to define the intersect of rasters
intersect_mask <- function(x){
  values_x <- getValues(x)
  inter_x <- values_x %*% rep(1,nlayers(x))
  mask <- setValues(subset(x,1),values = (inter_x>0))
  return(mask)
}
# keep only all cells that are defined for all layers

myExpl.eug.uni<-stack(reg.predictors.eug.uni)

myExpl.eug.uni <- stack(mask(myExpl.eug.uni, intersect_mask(myExpl.eug.uni)))


##Number of pseudo-absence: Lobo & Tognelli sugerem 100 pontos de pseudo-
#ausencia para cada ponto de presenca.



myBiomodData.eug.uni<- BIOMOD_FormatingData(resp.var = myResp.eug.uni,
                                            expl.var = myExpl.eug.uni,
                                            resp.xy = myRespXY.eug.uni,
                                            resp.name = myRespName.eug.uni,
                                            PA.nb.rep = 1,
                                            PA.nb.absences = 10000,
                                            PA.strategy = 'disk',
                                            PA.dist.max=500000) # com base em http://www.sciencedirect.com/science/article/pii/S0304380008005486
myBiomodData.eug.uni
save.image()
###MODELING

# 1. Defining Models Options using default options.
myBiomodOption.eug.uni <- BIOMOD_ModelingOptions()

# 3. Computing the models

myBiomodModelOut.eug.uni<- BIOMOD_Modeling(myBiomodData.eug.uni,
                                           models = c("MAXENT.Phillips"),
                                           models.options = myBiomodOption.eug.uni,
                                           NbRunEval=5, #conforme http://onlinelibrary.wiley.com/doi/10.1111/j.1523-1739.2006.00354.x/full
                                           DataSplit=70,
                                           Prevalence=0.5,
                                           VarImport=1,
                                           models.eval.meth = c('TSS'),
                                           eug.univeObj = TRUE,
                                           rescal.all.models = TRUE,
                                           do.full.models = FALSE,
                                           modeling.id = paste(myRespName.eug.uni,"FirstModeling",sep=""))
(myBiomodModelOut.eug.uni)
save.image()
# get all models evaluation

#Excelente TSS>0.75
#Bom 0.40>=TSS<0.75
#Ruim TSS<0.40


myBiomodModelEval.eug.uni <- get_evaluations(myBiomodModelOut.eug.uni)

myBiomodModelEval.eug.uni




get_variables_importance(myBiomodModelOut.eug.uni)
get_evaluations_matrix_eug.uni<-cbind(get_variables_importance(myBiomodModelOut.eug.uni)[1:9], #RUN 1
                                      get_variables_importance(myBiomodModelOut.eug.uni)[10:18], #RUN 2
                                      get_variables_importance(myBiomodModelOut.eug.uni)[19:27], #RUN 3
                                      get_variables_importance(myBiomodModelOut.eug.uni)[28:36], #RUN 4
                                      get_variables_importance(myBiomodModelOut.eug.uni)[37:45]) #RUN 5
colnames(get_evaluations_matrix_eug.uni)<-c("RUN1", "RUN2", "RUN3", "RUN4", "RUN5")
names(reg.predictors.eug.uni)
rownames(get_evaluations_matrix_eug.uni)<-names(reg.predictors.eug.uni)
get_evaluations_matrix_eug.uni_mean<-apply(get_evaluations_matrix_eug.uni, 1, mean)
which.max(get_evaluations_matrix_eug.uni_mean)



mymaxent.eug.uni <- BIOMOD_LoadModels(myBiomodModelOut.eug.uni, models='MAXENT.Phillips')
myRespPlot2D.eug.uni <- response.plot2(models  = mymaxent.eug.uni,
                                       Data = get_formal_data(myBiomodModelOut.eug.uni,'expl.var'), 
                                       show.variables= get_formal_data(myBiomodModelOut.eug.uni,'expl.var.names'),
                                       do.bivariate = FALSE,
                                       fixed.var.metric = 'mean',
                                       col = c("blue", "red"),
                                       legend = F,
                                       data_species = get_formal_data(myBiomodModelOut.eug.uni,'resp.var'),
                                       plot=F)

save.image()

bio_3_response_eug.uni<-gather(myRespPlot2D.eug.uni$bio3,model,Probability,
                               eug.uni_PA1_RUN1_MAXENT.Phillips:eug.uni_PA1_RUN5_MAXENT.Phillips)



p_eug.uni_bio3<-ggplot(data =bio_3_response_eug.uni) + 
  geom_smooth(mapping = aes(x = bio3, y = Probability))+
  labs(x = bquote('Bio 3'))
p_eug.uni_bio3



###Ensemble modeling
myBiomodEM.eug.uni<-BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut.eug.uni,
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

myBiomodProjection.eug.uni <- BIOMOD_Projection(modeling.output = myBiomodModelOut.eug.uni,
                                                new.env = stack(ma.presente.poly.eug.uni),
                                                proj.name = 'current_MA',
                                                selected.models = "all",
                                                binary.meth = 'TSS',
                                                compress = FALSE,
                                                build.clamping.mask = FALSE)

#Ensemble Forcasting

myBiomodEF.eug.uni.presente <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM.eug.uni,
  projection.output = myBiomodProjection.eug.uni,binary.meth = "TSS")

save.image()

plot(myBiomodEF.eug.uni.presente)

currentPred.ensemble.eug.uni <- stack("eug.uni/proj_current_MA/proj_current_MA_eug.uni_ensemble.grd")
plot(currentPred.ensemble.eug.uni[[1]])

plot(currentPred.ensemble.eug.uni[[1]],col = gray.colors(12, start = 0.9, end = 0.0))
plot(wrld_simpl, xlim=c(-60,-40), ylim=c(-34,-15), axes=TRUE, add=T)
polygon(poligono_mata_atlantica)


plot(currentPred.ensemble.eug.uni[[2]],col = gray.colors(12, start = 0.9, end = 0.0))
points(coord.eug.uni.names_thinned_dataset_full[[1]], cex=.6,col="black")
polygon(poligono_mata_atlantica)



########################## 2070 ##########################
################################################################


ma.2070.poly.eug.uni<-stack(ma.2070.new[[2]],
                                ma.2070.new[[3]],
                                ma.2070.new[[8]],
                                ma.2070.new[[9]],
                                ma.2070.new[[13]],
                                ma.2070.new[[14]],
                                ma.2070.new[[15]],
                                ma.2070.new[[18]],
                                ma.2070.new[[19]])



names(ma.2070.poly.eug.uni)<-names(ma.presente.poly.eug.uni)


myExplFuture.2070.85.ma<-stack(ma.2070.poly.eug.uni) #modelo HE


projection.2070.85.ma.eug.uni <-BIOMOD_Projection(modeling.output=myBiomodModelOut.eug.uni,
                                                  new.env= myExplFuture.2070.85.ma,
                                                  proj.name='eug.uni.2070.85',
                                                  selected.models =  "all",
                                                  Bin.trans=TRUE,
                                                  binary.meth ='TSS',
                                                  compress = 'gzip',
                                                  build.clamping.mask	= FALSE,
                                                  SaveObj=TRUE,
                                                  keep.in.memory=F,
                                                  do.stack=T,
                                                  silent=F)


eug.uni.2070.85<-get_predictions(projection.2070.85.ma.eug.uni)
eug.uni.2070.85
summary(eug.uni.2070.85)





################## ENSAMBLE MODELLING 2070  ################
#################################################################



EnsambleForecast.2070.85.ma.eug.uni <-BIOMOD_EnsembleForecasting (
  projection.output = projection.2070.85.ma.eug.uni,
  EM.output=myBiomodEM.eug.uni,
  selected.models = "all",
  binary.meth ='TSS') 




EnsambleForecast.2070.85.ma.eug.uni
plot(EnsambleForecast.2070.85.ma.eug.uni)



ensemble.2070.HE.85.ma.eug.uni <- stack("eug.uni/proj_eug.uni.2070.85/proj_eug.uni.2070.85_eug.uni_ensemble.grd")
ensemble.2070.HE.85.ma.eug.uni



plot(ensemble.2070.HE.85.ma.eug.uni[[1]],col = gray.colors(12, start = 0.9, end = 0))
points(coord.gym_klo.names_thinned_dataset_full[[1]],pch=20, col="black", cex=.7)

plot(ensemble.2070.HE.85.ma.eug.uni[[2]],col = gray.colors(12, start = 0.9, end = 0))



futurePred.ensemble.2070.85.eug.uni <- stack("eug.uni/proj_eug.uni.2070.85/proj_eug.uni.2070.85_eug.uni_ensemble.grd")



plot(currentPred.ensemble.eug.uni[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.2070.85.eug.uni[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS


currentPred.ensemble.bin.eug.uni<- stack("eug.uni/proj_current_MA/proj_current_MA_eug.uni_ensemble_TSSbin.grd")
futurePred.ensemble.bin.2070.85.eug.uni <- stack("eug.uni/proj_eug.uni.2070.85/proj_eug.uni.2070.85_eug.uni_ensemble_TSSbin.gri")


plot(currentPred.ensemble.bin.eug.uni[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.bin.2070.85.eug.uni[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS



myBiomodRangeSize.drimys.2070.85.eug.uni <- BIOMOD_RangeSize(
  CurrentPred=currentPred.ensemble.bin.eug.uni[[1]],
  FutureProj=futurePred.ensemble.bin.2070.85.eug.uni[[1]])
myBiomodRangeSize.drimys.2070.85.eug.uni$Compt.By.Models





plot(myBiomodRangeSize.drimys.2070.85.eug.uni$Diff.By.Pixel, col= gray.colors(4, start = 0, end = 1, gamma = 1, alpha = NULL))
plot(wrld_simpl,  axes="TRUE", add=T)

polygon(poligono_mata_atlantica)

dev.off()
save.image()
