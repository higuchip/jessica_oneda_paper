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
points(coord.per.gla, pch=19, cex=.51, col="black")


##Poligono MA

require(rgdal)
data.shape<-readOGR(dsn="/media/dendrologia/423EF1863EF172F1/User/Documents/orientados/Guilherme/IDS_09_Desm_Bioma_Mata_Atlantica_1_2", 
                    layer="IDS_09_Desm_Bioma_Mata_Atlantica_1_2")

poligono_mata_atlantica<-data.shape@polygons[[1]]@Polygons[[1]]@coords

polygon(poligono_mata_atlantica)
###Poligono especie


poligono_per.gla<-getpoly()

colnames(poligono_per.gla) = c("Long", "Lat")
polygon(poligono_per.gla) 



########REMOCAO DUVIDAS

out.per.gla= pnt.in.poly(coord.per.gla,poligono_per.gla)
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
plot(data.shape, add=T, col="lightgreen")
points(out.per.gla[which(out.per.gla$pip==1),1:2],pch=20, cex=.5, col="black")
coord.per.gla.limpo<-out.per.gla[which(out.per.gla$pip==1),1:2]
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
points(coord.per.gla.limpo,pch=20, cex=.5, col="black")


###FILTRAGEM ESPACIAl

head(coord.per.gla.limpo)
coord.per.gla.limpo<-na.omit(coord.per.gla.limpo)

dim(coord.per.gla.limpo)
per.gla_names<-c(rep("per.gla",  dim(coord.per.gla.limpo)[1]))
coord.per.gla.limpo.names<-cbind(per.gla_names, coord.per.gla.limpo)

head(coord.per.gla.limpo.names)
is.data.frame(coord.per.gla.limpo.names)
colnames(coord.per.gla.limpo.names)<-c("SPEC", "LONG", "LAT")

coord.per.gla.names_thinned_dataset_full <-
  thin( loc.data =coord.per.gla.limpo.names, 
        lat.col = "LAT", 
        long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 26.27, reps =1, #10 min resolution
        locs.thinned.list.return = TRUE, 
        write.files = T, 
        max.files = 5, 
        out.dir = "coord.per.gla.limpo.names_full/", out.base = "coord.per.gla.limpo.names_thinned", 
        write.log.file = TRUE,
        log.file = "coord.per.gla.limpo.names_thinned_full_log_file.txt" )




dim(coord.per.gla.names_thinned_dataset_full[[1]])



###ENVIRONMENTAL BACKGROUND

reg.variaveis.per.gla= crop(variaveis, poligono_per.gla)
plot(reg.variaveis.per.gla, 1)
points(coord.per.gla.limpo, pch=20, cex=.5)

#remoção variaveis nao colineares
vari.n.colineares.per.gla<-vifstep(reg.variaveis.per.gla)
vari.n.colineares.per.gla


reg.predictors.per.gla<- stack(reg.variaveis.per.gla[[2]],
                               reg.variaveis.per.gla[[3]],
                               reg.variaveis.per.gla[[8]],
                               reg.variaveis.per.gla[[9]],
                               reg.variaveis.per.gla[[13]],
                               reg.variaveis.per.gla[[14]],
                               reg.variaveis.per.gla[[15]],
                               reg.variaveis.per.gla[[18]],
                               reg.variaveis.per.gla[[19]])



names(reg.predictors.per.gla)


###MATA ATLANTICA projection presente

ma.presente = crop(variaveis, poligono_mata_atlantica)
plot(ma.presente, 1)
polygon(poligono_mata_atlantica)
ma.presente.poly <- mask(ma.presente, data.shape)
plot(ma.presente.poly, 1 )


ma.presente.poly.per.gla<-stack(ma.presente.new[[2]],
                                ma.presente.new[[3]],
                                ma.presente.new[[8]],
                                ma.presente.new[[9]],
                                ma.presente.new[[13]],
                                ma.presente.new[[14]],
                                ma.presente.new[[15]],
                                ma.presente.new[[18]],
                                ma.presente.new[[19]])




names(ma.presente.poly.per.gla)


myRespName.per.gla <-  'per.gla'
myResp.per.gla<-rep(1, nrow(coord.per.gla.names_thinned_dataset_full[[1]]))
myResp.per.gla
length(myResp.per.gla)
myRespXY.per.gla<-coord.per.gla.names_thinned_dataset_full[[1]]

# function to define the intersect of rasters
intersect_mask <- function(x){
  values_x <- getValues(x)
  inter_x <- values_x %*% rep(1,nlayers(x))
  mask <- setValues(subset(x,1),values = (inter_x>0))
  return(mask)
}
# keep only all cells that are defined for all layers

myExpl.per.gla<-stack(reg.predictors.per.gla)

myExpl.per.gla <- stack(mask(myExpl.per.gla, intersect_mask(myExpl.per.gla)))


##Number of pseudo-absence: Lobo & Tognelli sugerem 100 pontos de pseudo-
#ausencia para cada ponto de presenca.



myBiomodData.per.gla<- BIOMOD_FormatingData(resp.var = myResp.per.gla,
                                            expl.var = myExpl.per.gla,
                                            resp.xy = myRespXY.per.gla,
                                            resp.name = myRespName.per.gla,
                                            PA.nb.rep = 1,
                                            PA.nb.absences = 10000,
                                            PA.strategy = 'disk',
                                            PA.dist.max=500000) # com base em http://www.sciencedirect.com/science/article/pii/S0304380008005486
myBiomodData.per.gla
save.image()
###MODELING

# 1. Defining Models Options using default options.
myBiomodOption.per.gla <- BIOMOD_ModelingOptions()

# 3. Computing the models

myBiomodModelOut.per.gla<- BIOMOD_Modeling(myBiomodData.per.gla,
                                           models = c("MAXENT.Phillips"),
                                           models.options = myBiomodOption.per.gla,
                                           NbRunEval=5, #conforme http://onlinelibrary.wiley.com/doi/10.1111/j.1523-1739.2006.00354.x/full
                                           DataSplit=70,
                                           Prevalence=0.5,
                                           VarImport=1,
                                           models.eval.meth = c('TSS'),
                                           per.glaveObj = TRUE,
                                           rescal.all.models = TRUE,
                                           do.full.models = FALSE,
                                           modeling.id = paste(myRespName.per.gla,"FirstModeling",sep=""))
(myBiomodModelOut.per.gla)
save.image()
# get all models evaluation

#Excelente TSS>0.75
#Bom 0.40>=TSS<0.75
#Ruim TSS<0.40


myBiomodModelEval.per.gla <- get_evaluations(myBiomodModelOut.per.gla)

myBiomodModelEval.per.gla




get_variables_importance(myBiomodModelOut.per.gla)
get_evaluations_matrix_per.gla<-cbind(get_variables_importance(myBiomodModelOut.per.gla)[1:9], #RUN 1
                                      get_variables_importance(myBiomodModelOut.per.gla)[10:18], #RUN 2
                                      get_variables_importance(myBiomodModelOut.per.gla)[19:27], #RUN 3
                                      get_variables_importance(myBiomodModelOut.per.gla)[28:36], #RUN 4
                                      get_variables_importance(myBiomodModelOut.per.gla)[37:45]) #RUN 5
colnames(get_evaluations_matrix_per.gla)<-c("RUN1", "RUN2", "RUN3", "RUN4", "RUN5")
names(reg.predictors.per.gla)
rownames(get_evaluations_matrix_per.gla)<-names(reg.predictors.per.gla)
get_evaluations_matrix_per.gla_mean<-apply(get_evaluations_matrix_per.gla, 1, mean)
which.max(get_evaluations_matrix_per.gla_mean)



mymaxent.per.gla <- BIOMOD_LoadModels(myBiomodModelOut.per.gla, models='MAXENT.Phillips')
myRespPlot2D.per.gla <- response.plot2(models  = mymaxent.per.gla,
                                       Data = get_formal_data(myBiomodModelOut.per.gla,'expl.var'), 
                                       show.variables= get_formal_data(myBiomodModelOut.per.gla,'expl.var.names'),
                                       do.bivariate = FALSE,
                                       fixed.var.metric = 'mean',
                                       col = c("blue", "red"),
                                       legend = F,
                                       data_species = get_formal_data(myBiomodModelOut.per.gla,'resp.var'),
                                       plot=F)

save.image()

bio_8_response_per.gla<-gather(myRespPlot2D.per.gla$bio8,model,Probability,
                                per.gla_PA1_RUN1_MAXENT.Phillips:per.gla_PA1_RUN5_MAXENT.Phillips)



p_per.gla_bio8<-ggplot(data =bio_8_response_per.gla) + 
  geom_smooth(mapping = aes(x = bio8/10, y = Probability))+
  labs(x = bquote('Bio 8'))
p_per.gla_bio8



###Ensemble modeling
myBiomodEM.per.gla<-BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut.per.gla,
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

myBiomodProjection.per.gla <- BIOMOD_Projection(modeling.output = myBiomodModelOut.per.gla,
                                                new.env = stack(ma.presente.poly.per.gla),
                                                proj.name = 'current_MA',
                                                selected.models = "all",
                                                binary.meth = 'TSS',
                                                compress = FALSE,
                                                build.clamping.mask = FALSE)

#Ensemble Forcasting

myBiomodEF.per.gla.presente <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM.per.gla,
  projection.output = myBiomodProjection.per.gla,binary.meth = "TSS")

save.image()

plot(myBiomodEF.per.gla.presente)

currentPred.ensemble.per.gla <- stack("per.gla/proj_current_MA/proj_current_MA_per.gla_ensemble.grd")
plot(currentPred.ensemble.per.gla[[1]])

plot(currentPred.ensemble.per.gla[[1]],col = gray.colors(12, start = 0.9, end = 0.0))
plot(wrld_simpl, xlim=c(-60,-40), ylim=c(-34,-15), axes=TRUE, add=T)
polygon(poligono_mata_atlantica)


plot(currentPred.ensemble.per.gla[[2]],col = gray.colors(12, start = 0.9, end = 0.0))
points(coord.per.gla.names_thinned_dataset_full[[1]], cex=.6,col="black")
polygon(poligono_mata_atlantica)




########################## 2070 ##########################
################################################################



ma.2070.poly.per.gla<-stack(ma.2070.new[[2]],
                                ma.2070.new[[3]],
                                ma.2070.new[[8]],
                                ma.2070.new[[9]],
                                ma.2070.new[[13]],
                                ma.2070.new[[14]],
                                ma.2070.new[[15]],
                                ma.2070.new[[18]],
                                ma.2070.new[[19]])





names(ma.2070.poly.per.gla)<-names(ma.presente.poly.per.gla)


myExplFuture.2070.85.ma<-stack(ma.2070.poly.per.gla) #modelo HE


projection.2070.85.ma.per.gla <-BIOMOD_Projection(modeling.output=myBiomodModelOut.per.gla,
                                                  new.env= myExplFuture.2070.85.ma,
                                                  proj.name='per.gla.2070.85',
                                                  selected.models =  "all",
                                                  Bin.trans=TRUE,
                                                  binary.meth ='TSS',
                                                  compress = 'gzip',
                                                  build.clamping.mask	= FALSE,
                                                  SaveObj=TRUE,
                                                  keep.in.memory=F,
                                                  do.stack=T,
                                                  silent=F)


per.gla.2070.85<-get_predictions(projection.2070.85.ma.per.gla)
per.gla.2070.85
summary(per.gla.2070.85)





################## ENSAMBLE MODELLING 2070  ################
#################################################################



EnsambleForecast.2070.85.ma.per.gla <-BIOMOD_EnsembleForecasting (
  projection.output = projection.2070.85.ma.per.gla,
  EM.output=myBiomodEM.per.gla,
  selected.models = "all",
  binary.meth ='TSS') 




EnsambleForecast.2070.85.ma.per.gla
plot(EnsambleForecast.2070.85.ma.per.gla)



ensemble.2070.HE.85.ma.per.gla <- stack("per.gla/proj_per.gla.2070.85/proj_per.gla.2070.85_per.gla_ensemble.grd")
ensemble.2070.HE.85.ma.per.gla



plot(ensemble.2070.HE.85.ma.per.gla[[1]],col = gray.colors(12, start = 0.9, end = 0))
points(coord.gym_klo.names_thinned_dataset_full[[1]],pch=20, col="black", cex=.7)

plot(ensemble.2070.HE.85.ma.per.gla[[2]],col = gray.colors(12, start = 0.9, end = 0))



futurePred.ensemble.2070.85.per.gla <- stack("per.gla/proj_per.gla.2070.85/proj_per.gla.2070.85_per.gla_ensemble.grd")



plot(currentPred.ensemble.per.gla[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.2070.85.per.gla[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS


currentPred.ensemble.bin.per.gla<- stack("per.gla/proj_current_MA/proj_current_MA_per.gla_ensemble_TSSbin.grd")
futurePred.ensemble.bin.2070.85.per.gla <- stack("per.gla/proj_per.gla.2070.85/proj_per.gla.2070.85_per.gla_ensemble_TSSbin.gri")


plot(currentPred.ensemble.bin.per.gla[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.bin.2070.85.per.gla[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS



myBiomodRangeSize.2070.85.per.gla <- BIOMOD_RangeSize(
  CurrentPred=currentPred.ensemble.bin.per.gla[[1]],
  FutureProj=futurePred.ensemble.bin.2070.85.per.gla[[1]])
myBiomodRangeSize.2070.85.per.gla$Compt.By.Models





plot(myBiomodRangeSize.2070.85.per.gla$Diff.By.Pixel, col= gray.colors(4, start = 0, end = 1, gamma = 1, alpha = NULL))
plot(wrld_simpl,  axes="TRUE", add=T)

polygon(poligono_mata_atlantica)



