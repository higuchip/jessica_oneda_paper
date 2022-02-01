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
library(ff)

BioC<- getData('worldclim', var='bio', res=10) # Bioclima
variaveis <- stack(BioC)

plot(wrld_simpl, xlim=c(-130,-35), ylim=c(-45,35), axes=TRUE, col="gray")
plot(wrld_simpl, xlim=c(-70,-45), ylim=c(-35,-5), axes=TRUE, col="gray")
points(coord.myr.rac, pch=19, cex=.51, col="black")


##Poligono MA

require(rgdal)
data.shape<-readOGR(dsn="/media/dendrologia/423EF1863EF172F1/User/Documents/orientados/Guilherme/IDS_09_Desm_Bioma_Mata_Atlantica_1_2", 
                    layer="IDS_09_Desm_Bioma_Mata_Atlantica_1_2")

poligono_mata_atlantica<-data.shape@polygons[[1]]@Polygons[[1]]@coords

polygon(poligono_mata_atlantica)
###Poligono especie


poligono_myr.rac<-getpoly()

colnames(poligono_myr.rac) = c("Long", "Lat")
polygon(poligono_myr.rac) 



########REMOCAO DUVIDAS

out.myr.rac= pnt.in.poly(coord.myr.rac,poligono_myr.rac)
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
plot(data.shape, add=T, col="lightgreen")
points(out.myr.rac[which(out.myr.rac$pip==1),1:2],pch=20, cex=.5, col="black")
coord.myr.rac.limpo<-out.myr.rac[which(out.myr.rac$pip==1),1:2]
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
points(coord.myr.rac.limpo,pch=20, cex=.5, col="black")


###FILTRAGEM ESPACIAl

head(coord.myr.rac.limpo)
coord.myr.rac.limpo<-na.omit(coord.myr.rac.limpo)

dim(coord.myr.rac.limpo)
myr.rac_names<-c(rep("myr.rac",  dim(coord.myr.rac.limpo)[1]))
coord.myr.rac.limpo.names<-cbind(myr.rac_names, coord.myr.rac.limpo)

head(coord.myr.rac.limpo.names)
is.data.frame(coord.myr.rac.limpo.names)
colnames(coord.myr.rac.limpo.names)<-c("SPEC", "LONG", "LAT")
# 
# ?write.table
# write.table(coord.myr.rac.limpo.names, "tap_gui.csv", sep=";", dec=",")

coord.myr.rac.names_thinned_dataset_full <-
  thin( loc.data =coord.myr.rac.limpo.names, 
        lat.col = "LAT", 
        long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 26.27, reps =1, #10 min resolution
        locs.thinned.list.return = TRUE, 
        write.files = T, 
        max.files = 5, 
        out.dir = "coord.myr.rac.limpo.names_full/", out.base = "coord.myr.rac.limpo.names_thinned", 
        write.log.file = TRUE,
        log.file = "coord.myr.rac.limpo.names_thinned_full_log_file.txt" )



dim(coord.myr.rac.limpo.names)
dim(coord.myr.rac.names_thinned_dataset_full[[1]])



###ENVIRONMENTAL BACKGROUND

reg.variaveis.myr.rac= crop(variaveis, poligono_myr.rac)
plot(reg.variaveis.myr.rac, 1)
points(coord.myr.rac.limpo, pch=20, cex=.5)

#remoção variaveis nao colineares
vari.n.colineares.myr.rac<-vifstep(reg.variaveis.myr.rac)
vari.n.colineares.myr.rac


reg.predictors.myr.rac<- stack(reg.variaveis.myr.rac[[2]],
                               reg.variaveis.myr.rac[[3]],
                               reg.variaveis.myr.rac[[8]],
                               reg.variaveis.myr.rac[[9]],
                               reg.variaveis.myr.rac[[10]],
                               reg.variaveis.myr.rac[[13]],
                               reg.variaveis.myr.rac[[14]],
                               reg.variaveis.myr.rac[[15]],
                               reg.variaveis.myr.rac[[18]],
                               reg.variaveis.myr.rac[[19]])



names(reg.predictors.myr.rac)


###MATA ATLANTICA projection presente

ma.presente = crop(variaveis, poligono_mata_atlantica)
plot(ma.presente, 1)
polygon(poligono_mata_atlantica)
ma.presente.poly <- mask(ma.presente, data.shape)
plot(ma.presente.poly, 1 )


ma.presente.poly.myr.rac<-stack(ma.presente.new[[2]],
                                ma.presente.new[[3]],
                                ma.presente.new[[8]],
                                ma.presente.new[[9]],
                                ma.presente.new[[10]],
                                ma.presente.new[[13]],
                                ma.presente.new[[14]],
                                ma.presente.new[[15]],
                                ma.presente.new[[18]],
                                ma.presente.new[[19]])




names(ma.presente.poly.myr.rac)


myRespName.myr.rac <-  'myr.rac'
myResp.myr.rac<-rep(1, nrow(coord.myr.rac.names_thinned_dataset_full[[1]]))
myResp.myr.rac
length(myResp.myr.rac)
myRespXY.myr.rac<-coord.myr.rac.names_thinned_dataset_full[[1]]

# function to define the intersect of rasters
intersect_mask <- function(x){
  values_x <- getValues(x)
  inter_x <- values_x %*% rep(1,nlayers(x))
  mask <- setValues(subset(x,1),values = (inter_x>0))
  return(mask)
}
# keep only all cells that are defined for all layers

myExpl.myr.rac<-stack(reg.predictors.myr.rac)

myExpl.myr.rac <- stack(mask(myExpl.myr.rac, intersect_mask(myExpl.myr.rac)))


##Number of pseudo-absence: Lobo & Tognelli sugerem 100 pontos de pseudo-
#ausencia para cada ponto de presenca.



myBiomodData.myr.rac<- BIOMOD_FormatingData(resp.var = myResp.myr.rac,
                                            expl.var = myExpl.myr.rac,
                                            resp.xy = myRespXY.myr.rac,
                                            resp.name = myRespName.myr.rac,
                                            PA.nb.rep = 1,
                                            PA.nb.absences = 10000,
                                            PA.strategy = 'disk',
                                            PA.dist.max=500000) # com base em http://www.sciencedirect.com/science/article/pii/S0304380008005486
myBiomodData.myr.rac
save.image()
###MODELING

# 1. Defining Models Options using default options.
myBiomodOption.myr.rac <- BIOMOD_ModelingOptions()

# 3. Computing the models

myBiomodModelOut.myr.rac<- BIOMOD_Modeling(myBiomodData.myr.rac,
                                           models = c("MAXENT.Phillips"),
                                           models.options = myBiomodOption.myr.rac,
                                           NbRunEval=5, #conforme http://onlinelibrary.wiley.com/doi/10.1111/j.1523-1739.2006.00354.x/full
                                           DataSplit=70,
                                           Prevalence=0.5,
                                           VarImport=1,
                                           models.eval.meth = c('TSS'),
                                           myr.racveObj = TRUE,
                                           rescal.all.models = TRUE,
                                           do.full.models = FALSE,
                                           modeling.id = paste(myRespName.myr.rac,"FirstModeling",sep=""))
(myBiomodModelOut.myr.rac)
save.image()
# get all models evaluation

#Excelente TSS>0.75
#Bom 0.40>=TSS<0.75
#Ruim TSS<0.40


myBiomodModelEval.myr.rac <- get_evaluations(myBiomodModelOut.myr.rac)

myBiomodModelEval.myr.rac






get_variables_importance(myBiomodModelOut.myr.rac)
get_evaluations_matrix_myr.rac<-cbind(get_variables_importance(myBiomodModelOut.myr.rac)[1:10], #RUN 1
                                      get_variables_importance(myBiomodModelOut.myr.rac)[11:20], #RUN 2
                                      get_variables_importance(myBiomodModelOut.myr.rac)[21:30], #RUN 3
                                      get_variables_importance(myBiomodModelOut.myr.rac)[31:40], #RUN 4
                                      get_variables_importance(myBiomodModelOut.myr.rac)[41:50]) #RUN 5
colnames(get_evaluations_matrix_myr.rac)<-c("RUN1", "RUN2", "RUN3", "RUN4", "RUN5")
names(reg.predictors.myr.rac)
rownames(get_evaluations_matrix_myr.rac)<-names(reg.predictors.myr.rac)
get_evaluations_matrix_myr.rac_mean<-apply(get_evaluations_matrix_myr.rac, 1, mean)
which.max(get_evaluations_matrix_myr.rac_mean)



mymaxent.myr.rac <- BIOMOD_LoadModels(myBiomodModelOut.myr.rac, models='MAXENT.Phillips')
myRespPlot2D.myr.rac <- response.plot2(models  = mymaxent.myr.rac,
                                       Data = get_formal_data(myBiomodModelOut.myr.rac,'expl.var'), 
                                       show.variables= get_formal_data(myBiomodModelOut.myr.rac,'expl.var.names'),
                                       do.bivariate = FALSE,
                                       fixed.var.metric = 'mean',
                                       col = c("blue", "red"),
                                       legend = F,
                                       data_species = get_formal_data(myBiomodModelOut.myr.rac,'resp.var'),
                                       plot=F)

save.image()

bio_2_response_myr.rac<-gather(myRespPlot2D.myr.rac$bio2,model,Probability,
                               myr.rac_PA1_RUN1_MAXENT.Phillips:myr.rac_PA1_RUN5_MAXENT.Phillips)



p_myr.rac_bio2<-ggplot(data =bio_2_response_myr.rac) + 
  geom_smooth(mapping = aes(x = bio2/10, y = Probability))+
  labs(x = bquote('Bio 2'))
p_myr.rac_bio2



###Ensemble modeling
myBiomodEM.myr.rac<-BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut.myr.rac,
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

myBiomodProjection.myr.rac <- BIOMOD_Projection(modeling.output = myBiomodModelOut.myr.rac,
                                                new.env = stack(ma.presente.poly.myr.rac),
                                                proj.name = 'current_MA',
                                                selected.models = "all",
                                                binary.meth = 'TSS',
                                                compress = FALSE,
                                                build.clamping.mask = FALSE)

#Ensemble Forcasting

myBiomodEF.myr.rac.presente <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM.myr.rac,
  projection.output = myBiomodProjection.myr.rac,binary.meth = "TSS")

save.image()

plot(myBiomodEF.myr.rac.presente)

currentPred.ensemble.myr.rac <- stack("myr.rac/proj_current_MA/proj_current_MA_myr.rac_ensemble.grd")
plot(currentPred.ensemble.myr.rac[[1]])

plot(currentPred.ensemble.myr.rac[[1]],col = gray.colors(12, start = 0.9, end = 0.0))
plot(wrld_simpl, xlim=c(-60,-40), ylim=c(-34,-15), axes=TRUE, add=T)
polygon(poligono_mata_atlantica)
points(coord.myr.rac.names_thinned_dataset_full[[1]], cex=.6,col="black")

plot(currentPred.ensemble.myr.rac[[2]],col = gray.colors(12, start = 0.9, end = 0.0))
points(coord.myr.rac.names_thinned_dataset_full[[1]], cex=.6,col="black")
polygon(poligono_mata_atlantica)


?BIEN_metadata_citation

citations.myr.rac<-BIEN_metadata_citation(dataframe=myr.rac.bien)#If you are referencing occurrence data
citations.myr.rac



########################## 2070 ##########################
################################################################



ma.2070.poly.myr.rac<-stack(ma.2070.new[[2]],
                                ma.2070.new[[3]],
                                ma.2070.new[[8]],
                                ma.2070.new[[9]],
                                ma.2070.new[[10]],
                                ma.2070.new[[13]],
                                ma.2070.new[[14]],
                                ma.2070.new[[15]],
                                ma.2070.new[[18]],
                                ma.2070.new[[19]])


names(ma.2070.poly.myr.rac)<-names(ma.presente.poly.myr.rac)


myExplFuture.2070.85.ma<-stack(ma.2070.poly.myr.rac) #modelo HE


projection.2070.85.ma.myr.rac <-BIOMOD_Projection(modeling.output=myBiomodModelOut.myr.rac,
                                                  new.env= myExplFuture.2070.85.ma,
                                                  proj.name='myr.rac.2070.85',
                                                  selected.models =  "all",
                                                  Bin.trans=TRUE,
                                                  binary.meth ='TSS',
                                                  compress = 'gzip',
                                                  build.clamping.mask	= FALSE,
                                                  SaveObj=TRUE,
                                                  keep.in.memory=F,
                                                  do.stack=T,
                                                  silent=F)


myr.rac.2070.85<-get_predictions(projection.2070.85.ma.myr.rac)
myr.rac.2070.85
summary(myr.rac.2070.85)





################## ENSAMBLE MODELLING 2070  ################
#################################################################



EnsambleForecast.2070.85.ma.myr.rac <-BIOMOD_EnsembleForecasting (
  projection.output = projection.2070.85.ma.myr.rac,
  EM.output=myBiomodEM.myr.rac,
  selected.models = "all",
  binary.meth ='TSS') 




EnsambleForecast.2070.85.ma.myr.rac
plot(EnsambleForecast.2070.85.ma.myr.rac)



ensemble.2070.HE.85.ma.myr.rac <- stack("myr.rac/proj_myr.rac.2070.85/proj_myr.rac.2070.85_myr.rac_ensemble.grd")
ensemble.2070.HE.85.ma.myr.rac



plot(ensemble.2070.HE.85.ma.myr.rac[[1]],col = gray.colors(12, start = 0.9, end = 0))
points(coord.gym_klo.names_thinned_dataset_full[[1]],pch=20, col="black", cex=.7)

plot(ensemble.2070.HE.85.ma.myr.rac[[2]],col = gray.colors(12, start = 0.9, end = 0))



futurePred.ensemble.2070.85.myr.rac <- stack("myr.rac/proj_myr.rac.2070.85/proj_myr.rac.2070.85_myr.rac_ensemble.grd")



plot(currentPred.ensemble.myr.rac[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.2070.85.myr.rac[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS


currentPred.ensemble.bin.myr.rac<- stack("myr.rac/proj_current_MA/proj_current_MA_myr.rac_ensemble_TSSbin.grd")
futurePred.ensemble.bin.2070.85.myr.rac <- stack("myr.rac/proj_myr.rac.2070.85/proj_myr.rac.2070.85_myr.rac_ensemble_TSSbin.gri")


plot(currentPred.ensemble.bin.myr.rac[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.bin.2070.85.myr.rac[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS



myBiomodRangeSize.2070.85.myr.rac <- BIOMOD_RangeSize(
  CurrentPred=currentPred.ensemble.bin.myr.rac[[1]],
  FutureProj=futurePred.ensemble.bin.2070.85.myr.rac[[1]])
myBiomodRangeSize.2070.85.myr.rac$Compt.By.Models





plot(myBiomodRangeSize.2070.85.myr.rac$Diff.By.Pixel, col= gray.colors(4, start = 0, end = 1, gamma = 1, alpha = NULL))
plot(wrld_simpl,  axes="TRUE", add=T)

polygon(poligono_mata_atlantica)








save.image()




