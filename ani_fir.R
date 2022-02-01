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
points(coord.ani.fir, pch=19, cex=.51, col="black")


##Poligono MA

require(rgdal)
data.shape<-readOGR(dsn="/media/dendrologia/423EF1863EF172F1/User/Documents/orientados/Guilherme/IDS_09_Desm_Bioma_Mata_Atlantica_1_2", 
                    layer="IDS_09_Desm_Bioma_Mata_Atlantica_1_2")

poligono_mata_atlantica<-data.shape@polygons[[1]]@Polygons[[1]]@coords

polygon(poligono_mata_atlantica)
###Poligono especie


poligono_ani.fir<-getpoly()

colnames(poligono_ani.fir) = c("Long", "Lat")
polygon(poligono_ani.fir) 



########REMOCAO DUVIDAS

out.ani.fir= pnt.in.poly(coord.ani.fir,poligono_ani.fir)
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
plot(data.shape, add=T, col="lightgreen")
points(out.ani.fir[which(out.ani.fir$pip==1),1:2],pch=20, cex=.5, col="black")
coord.ani.fir.limpo<-out.ani.fir[which(out.ani.fir$pip==1),1:2]
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
points(coord.ani.fir.limpo,pch=20, cex=.5, col="black")


###FILTRAGEM ESPACIAl

head(coord.ani.fir.limpo)
coord.ani.fir.limpo<-na.omit(coord.ani.fir.limpo)

dim(coord.ani.fir.limpo)
ani.fir_names<-c(rep("ani.fir",  dim(coord.ani.fir.limpo)[1]))
coord.ani.fir.limpo.names<-cbind(ani.fir_names, coord.ani.fir.limpo)

head(coord.ani.fir.limpo.names)
is.data.frame(coord.ani.fir.limpo.names)
colnames(coord.ani.fir.limpo.names)<-c("SPEC", "LONG", "LAT")

coord.ani.fir.names_thinned_dataset_full <-
  thin( loc.data =coord.ani.fir.limpo.names, 
        lat.col = "LAT", 
        long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 26.27, reps =1, #10 min resolution
        locs.thinned.list.return = TRUE, 
        write.files = T, 
        max.files = 5, 
        out.dir = "coord.ani.fir.limpo.names_full/", out.base = "coord.ani.fir.limpo.names_thinned", 
        write.log.file = TRUE,
        log.file = "coord.ani.fir.limpo.names_thinned_full_log_file.txt" )




dim(coord.ani.fir.names_thinned_dataset_full[[1]])



###ENVIRONMENTAL BACKGROUND

reg.variaveis.ani.fir= crop(variaveis, poligono_ani.fir)
plot(reg.variaveis.ani.fir, 1)
points(coord.ani.fir.limpo, pch=20, cex=.5)

#remoção variaveis nao colineares
vari.n.colineares.ani.fir<-vifstep(reg.variaveis.ani.fir)
vari.n.colineares.ani.fir


reg.predictors.ani.fir<- stack(reg.variaveis.ani.fir[[2]],
                               reg.variaveis.ani.fir[[3]],
                               reg.variaveis.ani.fir[[8]],
                               reg.variaveis.ani.fir[[9]],
                               reg.variaveis.ani.fir[[13]],
                               reg.variaveis.ani.fir[[14]],
                               reg.variaveis.ani.fir[[15]],
                               reg.variaveis.ani.fir[[18]],
                               reg.variaveis.ani.fir[[19]])



names(reg.predictors.ani.fir)


###MATA ATLANTICA projection presente

ma.presente = crop(variaveis, poligono_mata_atlantica)
plot(ma.presente, 1)
polygon(poligono_mata_atlantica)
ma.presente.poly <- mask(ma.presente, data.shape)
plot(ma.presente.poly, 1 )


ma.presente.poly.ani.fir<-stack(ma.presente.new[[2]],
                                ma.presente.new[[3]],
                                ma.presente.new[[8]],
                                ma.presente.new[[9]],
                                ma.presente.new[[13]],
                                ma.presente.new[[14]],
                                ma.presente.new[[15]],
                                ma.presente.new[[18]],
                                ma.presente.new[[19]])




names(ma.presente.poly.ani.fir)


myRespName.ani.fir <-  'ani.fir'
myResp.ani.fir<-rep(1, nrow(coord.ani.fir.names_thinned_dataset_full[[1]]))
myResp.ani.fir
length(myResp.ani.fir)
myRespXY.ani.fir<-coord.ani.fir.names_thinned_dataset_full[[1]]

# function to define the intersect of rasters
intersect_mask <- function(x){
  values_x <- getValues(x)
  inter_x <- values_x %*% rep(1,nlayers(x))
  mask <- setValues(subset(x,1),values = (inter_x>0))
  return(mask)
}
# keep only all cells that are defined for all layers

myExpl.ani.fir<-stack(reg.predictors.ani.fir)

myExpl.ani.fir <- stack(mask(myExpl.ani.fir, intersect_mask(myExpl.ani.fir)))


##Number of pseudo-absence: Lobo & Tognelli sugerem 100 pontos de pseudo-
#ausencia para cada ponto de presenca.



myBiomodData.ani.fir<- BIOMOD_FormatingData(resp.var = myResp.ani.fir,
                                            expl.var = myExpl.ani.fir,
                                            resp.xy = myRespXY.ani.fir,
                                            resp.name = myRespName.ani.fir,
                                            PA.nb.rep = 1,
                                            PA.nb.absences = 10000,
                                            PA.strategy = 'disk',
                                            PA.dist.max=500000) # com base em http://www.sciencedirect.com/science/article/pii/S0304380008005486
myBiomodData.ani.fir
save.image()
###MODELING

# 1. Defining Models Options using default options.
myBiomodOption.ani.fir <- BIOMOD_ModelingOptions()

# 3. Computing the models

myBiomodModelOut.ani.fir<- BIOMOD_Modeling(myBiomodData.ani.fir,
                                           models = c("MAXENT.Phillips"),
                                           models.options = myBiomodOption.ani.fir,
                                           NbRunEval=5, #conforme http://onlinelibrary.wiley.com/doi/10.1111/j.1523-1739.2006.00354.x/full
                                           DataSplit=70,
                                           Prevalence=0.5,
                                           VarImport=1,
                                           models.eval.meth = c('TSS'),
                                           ani.firveObj = TRUE,
                                           rescal.all.models = TRUE,
                                           do.full.models = FALSE,
                                           modeling.id = paste(myRespName.ani.fir,"FirstModeling",sep=""))
(myBiomodModelOut.ani.fir)
save.image()
# get all models evaluation

#Excelente TSS>0.75
#Bom 0.40>=TSS<0.75
#Ruim TSS<0.40


myBiomodModelEval.ani.fir <- get_evaluations(myBiomodModelOut.ani.fir)

myBiomodModelEval.ani.fir




get_variables_importance(myBiomodModelOut.ani.fir)
get_evaluations_matrix_ani.fir<-cbind(get_variables_importance(myBiomodModelOut.ani.fir)[1:9], #RUN 1
                                      get_variables_importance(myBiomodModelOut.ani.fir)[10:18], #RUN 2
                                      get_variables_importance(myBiomodModelOut.ani.fir)[19:27], #RUN 3
                                      get_variables_importance(myBiomodModelOut.ani.fir)[28:36], #RUN 4
                                      get_variables_importance(myBiomodModelOut.ani.fir)[37:45]) #RUN 5
colnames(get_evaluations_matrix_ani.fir)<-c("RUN1", "RUN2", "RUN3", "RUN4", "RUN5")
names(reg.predictors.ani.fir)
rownames(get_evaluations_matrix_ani.fir)<-names(reg.predictors.ani.fir)
get_evaluations_matrix_ani.fir_mean<-apply(get_evaluations_matrix_ani.fir, 1, mean)
which.max(get_evaluations_matrix_ani.fir_mean)



mymaxent.ani.fir <- BIOMOD_LoadModels(myBiomodModelOut.ani.fir, models='MAXENT.Phillips')
myRespPlot2D.ani.fir <- response.plot2(models  = mymaxent.ani.fir,
                                       Data = get_formal_data(myBiomodModelOut.ani.fir,'expl.var'), 
                                       show.variables= get_formal_data(myBiomodModelOut.ani.fir,'expl.var.names'),
                                       do.bivariate = FALSE,
                                       fixed.var.metric = 'mean',
                                       col = c("blue", "red"),
                                       legend = F,
                                       data_species = get_formal_data(myBiomodModelOut.ani.fir,'resp.var'),
                                       plot=F)

save.image()

bio_14_response_ani.fir<-gather(myRespPlot2D.ani.fir$bio14,model,Probability,
                               ani.fir_PA1_RUN1_MAXENT.Phillips:ani.fir_PA1_RUN5_MAXENT.Phillips)



p_ani.fir_bio14<-ggplot(data =bio_14_response_ani.fir) + 
  geom_smooth(mapping = aes(x = bio14, y = Probability))+
  labs(x = bquote('Bio 14'))
p_ani.fir_bio14



###Ensemble modeling
myBiomodEM.ani.fir<-BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut.ani.fir,
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

myBiomodProjection.ani.fir <- BIOMOD_Projection(modeling.output = myBiomodModelOut.ani.fir,
                                                new.env = stack(ma.presente.poly.ani.fir),
                                                proj.name = 'current_MA',
                                                selected.models = "all",
                                                binary.meth = 'TSS',
                                                compress = FALSE,
                                                build.clamping.mask = FALSE)

#Ensemble Forcasting

myBiomodEF.ani.fir.presente <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM.ani.fir,
  projection.output = myBiomodProjection.ani.fir,binary.meth = "TSS")

save.image()

plot(myBiomodEF.ani.fir.presente)

currentPred.ensemble.ani.fir <- stack("ani.fir/proj_current_MA/proj_current_MA_ani.fir_ensemble.grd")
plot(currentPred.ensemble.ani.fir[[1]])

plot(currentPred.ensemble.ani.fir[[1]],col = gray.colors(12, start = 0.9, end = 0.0))
plot(wrld_simpl, xlim=c(-60,-40), ylim=c(-34,-15), axes=TRUE, add=T)
polygon(poligono_mata_atlantica)


plot(currentPred.ensemble.ani.fir[[2]],col = gray.colors(12, start = 0.9, end = 0.0))
points(coord.ani.fir.names_thinned_dataset_full[[1]], cex=.6,col="black")
polygon(poligono_mata_atlantica)

########################## 2070 ##########################
################################################################



ma.2070.poly.ani.fir<-stack(ma.2070.new[[2]],
                                ma.2070.new[[3]],
                                ma.2070.new[[8]],
                                ma.2070.new[[9]],
                                ma.2070.new[[13]],
                                ma.2070.new[[14]],
                                ma.2070.new[[15]],
                                ma.2070.new[[18]],
                                ma.2070.new[[19]])


names(ma.2070.poly.ani.fir)<-names(ma.presente.poly.ani.fir)


myExplFuture.2070.85.ma<-stack(ma.2070.poly.ani.fir) #modelo HE


projection.2070.85.ma.ani.fir <-BIOMOD_Projection(modeling.output=myBiomodModelOut.ani.fir,
                                                  new.env= myExplFuture.2070.85.ma,
                                                  proj.name='ani.fir.2070.85',
                                                  selected.models =  "all",
                                                  Bin.trans=TRUE,
                                                  binary.meth ='TSS',
                                                  compress = 'gzip',
                                                  build.clamping.mask	= FALSE,
                                                  SaveObj=TRUE,
                                                  keep.in.memory=F,
                                                  do.stack=T,
                                                  silent=F)


ani.fir.2070.85<-get_predictions(projection.2070.85.ma.ani.fir)
ani.fir.2070.85
summary(ani.fir.2070.85)





################## ENSAMBLE MODELLING 2070  ################
#################################################################



EnsambleForecast.2070.85.ma.ani.fir <-BIOMOD_EnsembleForecasting (
  projection.output = projection.2070.85.ma.ani.fir,
  EM.output=myBiomodEM.ani.fir,
  selected.models = "all",
  binary.meth ='TSS') 




EnsambleForecast.2070.85.ma.ani.fir
plot(EnsambleForecast.2070.85.ma.ani.fir)



ensemble.2070.HE.85.ma.ani.fir <- stack("ani.fir/proj_ani.fir.2070.85/proj_ani.fir.2070.85_ani.fir_ensemble.grd")
ensemble.2070.HE.85.ma.ani.fir



plot(ensemble.2070.HE.85.ma.ani.fir[[1]],col = gray.colors(12, start = 0.9, end = 0))
points(coord.gym_klo.names_thinned_dataset_full[[1]],pch=20, col="black", cex=.7)

plot(ensemble.2070.HE.85.ma.ani.fir[[2]],col = gray.colors(12, start = 0.9, end = 0))



futurePred.ensemble.2070.85.ani.fir <- stack("ani.fir/proj_ani.fir.2070.85/proj_ani.fir.2070.85_ani.fir_ensemble.grd")



plot(currentPred.ensemble.ani.fir[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.2070.85.ani.fir[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS


currentPred.ensemble.bin.ani.fir<- stack("ani.fir/proj_current_MA/proj_current_MA_ani.fir_ensemble_TSSbin.grd")
futurePred.ensemble.bin.2070.85.ani.fir <- stack("ani.fir/proj_ani.fir.2070.85/proj_ani.fir.2070.85_ani.fir_ensemble_TSSbin.gri")


plot(currentPred.ensemble.bin.ani.fir[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.bin.2070.85.ani.fir[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS



myBiomodRangeSize.2070.85.ani.fir <- BIOMOD_RangeSize(
  CurrentPred=currentPred.ensemble.bin.ani.fir[[1]],
  FutureProj=futurePred.ensemble.bin.2070.85.ani.fir[[1]])
myBiomodRangeSize.2070.85.ani.fir$Compt.By.Models





plot(myBiomodRangeSize.2070.85.ani.fir$Diff.By.Pixel, col= gray.colors(4, start = 0, end = 1, gamma = 1, alpha = NULL))
plot(wrld_simpl,  axes="TRUE", add=T)

polygon(poligono_mata_atlantica)

save.image()






