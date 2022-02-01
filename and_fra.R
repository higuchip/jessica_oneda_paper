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
points(coord.and.fra, pch=19, cex=.51, col="black")


##Poligono MA

require(rgdal)
data.shape<-readOGR(dsn="/media/dendrologia/423EF1863EF172F1/User/Documents/orientados/Guilherme/IDS_09_Desm_Bioma_Mata_Atlantica_1_2", 
                    layer="IDS_09_Desm_Bioma_Mata_Atlantica_1_2")

poligono_mata_atlantica<-data.shape@polygons[[1]]@Polygons[[1]]@coords

polygon(poligono_mata_atlantica)
###Poligono especie


poligono_and.fra<-getpoly()

colnames(poligono_and.fra) = c("Long", "Lat")
polygon(poligono_and.fra) 



########REMOCAO DUVIDAS

out.and.fra= pnt.in.poly(coord.and.fra,poligono_and.fra)
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
plot(data.shape, add=T, col="lightgreen")
points(out.and.fra[which(out.and.fra$pip==1),1:2],pch=20, cex=.5, col="black")
coord.and.fra.limpo<-out.and.fra[which(out.and.fra$pip==1),1:2]
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
points(coord.and.fra.limpo,pch=20, cex=.5, col="black")


###FILTRAGEM ESPACIAl

head(coord.and.fra.limpo)
coord.and.fra.limpo<-na.omit(coord.and.fra.limpo)

dim(coord.and.fra.limpo)
and.fra_names<-c(rep("and.fra",  dim(coord.and.fra.limpo)[1]))
coord.and.fra.limpo.names<-cbind(and.fra_names, coord.and.fra.limpo)

head(coord.and.fra.limpo.names)
is.data.frame(coord.and.fra.limpo.names)
colnames(coord.and.fra.limpo.names)<-c("SPEC", "LONG", "LAT")

coord.and.fra.names_thinned_dataset_full <-
  thin( loc.data =coord.and.fra.limpo.names, 
        lat.col = "LAT", 
        long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 26.27, reps =1, #10 min resolution
        locs.thinned.list.return = TRUE, 
        write.files = T, 
        max.files = 5, 
        out.dir = "coord.and.fra.limpo.names_full/", out.base = "coord.and.fra.limpo.names_thinned", 
        write.log.file = TRUE,
        log.file = "coord.and.fra.limpo.names_thinned_full_log_file.txt" )




dim(coord.and.fra.names_thinned_dataset_full[[1]])



###ENVIRONMENTAL BACKGROUND

reg.variaveis.and.fra= crop(variaveis, poligono_and.fra)
plot(reg.variaveis.and.fra, 1)
points(coord.and.fra.limpo, pch=20, cex=.5)

#remoção variaveis nao colineares
library(usdm)
vari.n.colineares.and.fra<-vifstep(reg.variaveis.and.fra)
vari.n.colineares.and.fra


reg.predictors.and.fra<- stack(reg.variaveis.and.fra[[2]],
                               reg.variaveis.and.fra[[3]],
                               reg.variaveis.and.fra[[8]],
                               reg.variaveis.and.fra[[9]],
                               reg.variaveis.and.fra[[10]],
                               reg.variaveis.and.fra[[13]],
                               reg.variaveis.and.fra[[14]],
                               reg.variaveis.and.fra[[15]],
                               reg.variaveis.and.fra[[18]],
                               reg.variaveis.and.fra[[19]])



plot(reg.predictors.and.fra,1)


###MATA ATLANTICA projection presente

ma.presente = crop(variaveis, poligono_mata_atlantica)
plot(ma.presente, 1)
polygon(poligono_mata_atlantica)
ma.presente.poly <- mask(ma.presente, data.shape)
plot(ma.presente.poly, 1 )


ma.presente.poly.and.fra<-stack(ma.presente.new[[2]],
                                ma.presente.new[[3]],
                                ma.presente.new[[8]],
                                ma.presente.new[[9]],
                                ma.presente.new[[10]],
                                ma.presente.new[[13]],
                                ma.presente.new[[14]],
                                ma.presente.new[[15]],
                                ma.presente.new[[18]],
                                ma.presente.new[[19]])


plot(ma.presente.poly.and.fra,2)
plot(ma.2070.poly.and.fra,2)

names(ma.presente.poly.and.fra)


myRespName.and.fra <-  'and.fra'
myResp.and.fra<-rep(1, nrow(coord.and.fra.names_thinned_dataset_full[[1]]))
myResp.and.fra
length(myResp.and.fra)
myRespXY.and.fra<-coord.and.fra.names_thinned_dataset_full[[1]]

# function to define the intersect of rasters
intersect_mask <- function(x){
  values_x <- getValues(x)
  inter_x <- values_x %*% rep(1,nlayers(x))
  mask <- setValues(subset(x,1),values = (inter_x>0))
  return(mask)
}
# keep only all cells that are defined for all layers

myExpl.and.fra<-stack(reg.predictors.and.fra)
plot(reg.predictors.and.fra,1)
myExpl.and.fra <- stack(mask(myExpl.and.fra, intersect_mask(myExpl.and.fra)))


myExpl.and.fra

##Number of pseudo-absence: Lobo & Tognelli sugerem 100 pontos de pseudo-
#ausencia para cada ponto de presenca.



myBiomodData.and.fra<- BIOMOD_FormatingData(resp.var = myResp.and.fra,
                                            expl.var = myExpl.and.fra,
                                            resp.xy = myRespXY.and.fra,
                                            resp.name = myRespName.and.fra,
                                            PA.nb.rep = 1,
                                            PA.nb.absences = 10000,
                                            PA.strategy = 'disk',
                                            PA.dist.max=500000) # com base em http://www.sciencedirect.com/science/article/pii/S0304380008005486
myBiomodData.and.fra
save.image()
###MODELING

# 1. Defining Models Options using default options.
myBiomodOption.and.fra <- BIOMOD_ModelingOptions()

# 3. Computing the models

myBiomodModelOut.and.fra<- BIOMOD_Modeling(myBiomodData.and.fra,
                                           models = c("MAXENT.Phillips"),
                                           models.options = myBiomodOption.and.fra,
                                           NbRunEval=5, #conforme http://onlinelibrary.wiley.com/doi/10.1111/j.1523-1739.2006.00354.x/full
                                           DataSplit=70,
                                           Prevalence=0.5,
                                           VarImport=1,
                                           models.eval.meth = c('TSS'),
                                           and.fraveObj = TRUE,
                                           rescal.all.models = TRUE,
                                           do.full.models = FALSE,
                                           modeling.id = paste(myRespName.and.fra,"FirstModeling",sep=""))
(myBiomodModelOut.and.fra)
save.image()
# get all models evaluation

#Excelente TSS>0.75
#Bom 0.40>=TSS<0.75
#Ruim TSS<0.40


myBiomodModelEval.and.fra <- get_evaluations(myBiomodModelOut.and.fra)

myBiomodModelEval.and.fra




get_variables_importance(myBiomodModelOut.and.fra)
get_evaluations_matrix_and.fra<-cbind(get_variables_importance(myBiomodModelOut.and.fra)[1:10], #RUN 1
                                      get_variables_importance(myBiomodModelOut.and.fra)[11:20], #RUN 2
                                      get_variables_importance(myBiomodModelOut.and.fra)[21:30], #RUN 3
                                      get_variables_importance(myBiomodModelOut.and.fra)[31:40], #RUN 4
                                      get_variables_importance(myBiomodModelOut.and.fra)[41:50]) #RUN 5
colnames(get_evaluations_matrix_and.fra)<-c("RUN1", "RUN2", "RUN3", "RUN4", "RUN5")
names(reg.predictors.and.fra)
rownames(get_evaluations_matrix_and.fra)<-names(reg.predictors.and.fra)
get_evaluations_matrix_and.fra_mean<-apply(get_evaluations_matrix_and.fra, 1, mean)
which.max(get_evaluations_matrix_and.fra_mean)



mymaxent.and.fra <- BIOMOD_LoadModels(myBiomodModelOut.and.fra, models='MAXENT.Phillips')
myRespPlot2D.and.fra <- response.plot2(models  = mymaxent.and.fra,
                                       Data = get_formal_data(myBiomodModelOut.and.fra,'expl.var'), 
                                       show.variables= get_formal_data(myBiomodModelOut.and.fra,'expl.var.names'),
                                       do.bivariate = FALSE,
                                       fixed.var.metric = 'mean',
                                       col = c("blue", "red"),
                                       legend = F,
                                       data_species = get_formal_data(myBiomodModelOut.and.fra,'resp.var'),
                                       plot=F)

save.image()

bio_15_response_and.fra<-gather(myRespPlot2D.and.fra$bio15,model,Probability,
                               and.fra_PA1_RUN1_MAXENT.Phillips:and.fra_PA1_RUN5_MAXENT.Phillips)



p_and.fra_bio15<-ggplot(data =bio_15_response_and.fra) + 
  geom_smooth(mapping = aes(x = bio15, y = Probability))+
  labs(x = bquote('Bio 15'))
p_and.fra_bio15



###Ensemble modeling
myBiomodEM.and.fra<-BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut.and.fra,
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

myBiomodProjection.and.fra <- BIOMOD_Projection(modeling.output = myBiomodModelOut.and.fra,
                                                new.env = stack(ma.presente.poly.and.fra),
                                                proj.name = 'current_MA',
                                                selected.models = "all",
                                                binary.meth = 'TSS',
                                                compress = FALSE,
                                                build.clamping.mask = FALSE)

#Ensemble Forcasting

myBiomodEF.and.fra.presente <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM.and.fra,
  projection.output = myBiomodProjection.and.fra,binary.meth = "TSS")

save.image()

plot(myBiomodEF.and.fra.presente)

currentPred.ensemble.and.fra <- stack("and.fra/proj_current_MA/proj_current_MA_and.fra_ensemble.grd")
plot(currentPred.ensemble.and.fra[[1]])

plot(currentPred.ensemble.and.fra[[1]],col = gray.colors(12, start = 0.9, end = 0.0))
plot(wrld_simpl, xlim=c(-60,-40), ylim=c(-34,-15), axes=TRUE, add=T)
polygon(poligono_mata_atlantica)


plot(currentPred.ensemble.and.fra[[2]],col = gray.colors(12, start = 0.9, end = 0.0))
points(coord.and.fra.names_thinned_dataset_full[[1]], cex=1,col="black", pch=20)
polygon(poligono_mata_atlantica)


########################## 2070 ##########################
################################################################



ma.2070.poly.and.fra<-stack(ma.2070.new[[2]],
                                ma.2070.new[[3]],
                                ma.2070.new[[8]],
                                ma.2070.new[[9]],
                                ma.2070.new[[10]],
                                ma.2070.new[[13]],
                                ma.2070.new[[14]],
                                ma.2070.new[[15]],
                                ma.2070.new[[18]],
                                ma.2070.new[[19]])

ma.presente.poly.and.fra
plot(ma.2070.poly.and.fra,1)
names(ma.2070.poly.and.fra)<-names(ma.presente.poly.and.fra)


myExplFuture.2070.85.ma<-stack(ma.2070.poly.and.fra) #modelo HE


projection.2070.85.ma.and.fra <-BIOMOD_Projection(modeling.output=myBiomodModelOut.and.fra,
                                                  new.env= myExplFuture.2070.85.ma,
                                                  proj.name='and.fra.2070.85',
                                                  selected.models =  "all",
                                                  Bin.trans=TRUE,
                                                  binary.meth ='TSS',
                                                  compress = 'gzip',
                                                  build.clamping.mask	= FALSE,
                                                  SaveObj=TRUE,
                                                  keep.in.memory=F,
                                                  do.stack=T,
                                                  silent=F)

and.fra.2070.85<-get_predictions(projection.2070.85.ma.and.fra)
and.fra.2070.85
summary(and.fra.2070.85)






################## ENSAMBLE MODELLING 2070  ################
#################################################################




EnsambleForecast.2070.85.ma.and.fra <-BIOMOD_EnsembleForecasting (
  projection.output = projection.2070.85.ma.and.fra,
  EM.output=myBiomodEM.and.fra,
  selected.models = "all",
  binary.meth ='TSS') 




EnsambleForecast.2070.85.ma.and.fra
plot(EnsambleForecast.2070.85.ma.and.fra)



ensemble.2070.HE.85.ma.and.fra <- stack("and.fra/proj_and.fra.2070.85/proj_and.fra.2070.85_and.fra_ensemble.grd")
ensemble.2070.HE.85.ma.and.fra



plot(ensemble.2070.HE.85.ma.and.fra[[1]],col = gray.colors(12, start = 0.9, end = 0))
points(coord.gym_klo.names_thinned_dataset_full[[1]],pch=20, col="black", cex=.7)

plot(ensemble.2070.HE.85.ma.and.fra[[2]],col = gray.colors(12, start = 0.9, end = 0))



futurePred.ensemble.2070.85.and.fra <- stack("and.fra/proj_and.fra.2070.85/proj_and.fra.2070.85_and.fra_ensemble.grd")



plot(currentPred.ensemble.and.fra[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.2070.85.and.fra[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS


currentPred.ensemble.bin.and.fra<- stack("and.fra/proj_current_MA/proj_current_MA_and.fra_ensemble_TSSbin.grd")
futurePred.ensemble.bin.2070.85.and.fra <- stack("and.fra/proj_and.fra.2070.85/proj_and.fra.2070.85_and.fra_ensemble_TSSbin.gri")


plot(currentPred.ensemble.bin.and.fra[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.bin.2070.85.and.fra[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS



myBiomodRangeSize.drimys.2070.85.and.fra <- BIOMOD_RangeSize(
  CurrentPred=currentPred.ensemble.bin.and.fra[[1]],
  FutureProj=futurePred.ensemble.bin.2070.85.and.fra[[1]])
myBiomodRangeSize.drimys.2070.85.and.fra$Compt.By.Models





plot(myBiomodRangeSize.drimys.2070.85.and.fra$Diff.By.Pixel, col= gray.colors(4, start = 0, end = 1, gamma = 1, alpha = NULL))
plot(wrld_simpl,  axes="TRUE", add=T)

polygon(poligono_mata_atlantica)



