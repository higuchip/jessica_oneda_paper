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
points(coord.cec.gla, pch=19, cex=.51, col="black")


##Poligono MA

require(rgdal)
data.shape<-readOGR(dsn="/media/dendrologia/423EF1863EF172F1/User/Documents/orientados/Guilherme/IDS_09_Desm_Bioma_Mata_Atlantica_1_2", 
                    layer="IDS_09_Desm_Bioma_Mata_Atlantica_1_2")

poligono_mata_atlantica<-data.shape@polygons[[1]]@Polygons[[1]]@coords

polygon(poligono_mata_atlantica)
###Poligono especie


poligono_cec.gla<-getpoly()

colnames(poligono_cec.gla) = c("Long", "Lat")
polygon(poligono_cec.gla) 



########REMOCAO DUVIDAS

out.cec.gla= pnt.in.poly(coord.cec.gla,poligono_cec.gla)
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
plot(data.shape, add=T, col="lightgreen")
points(out.cec.gla[which(out.cec.gla$pip==1),1:2],pch=20, cex=.5, col="black")
coord.cec.gla.limpo<-out.cec.gla[which(out.cec.gla$pip==1),1:2]
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
points(coord.cec.gla.limpo,pch=20, cex=.5, col="black")


###FILTRAGEM ESPACIAl

head(coord.cec.gla.limpo)
coord.cec.gla.limpo<-na.omit(coord.cec.gla.limpo)

dim(coord.cec.gla.limpo)
cec.gla_names<-c(rep("cec.gla",  dim(coord.cec.gla.limpo)[1]))
coord.cec.gla.limpo.names<-cbind(cec.gla_names, coord.cec.gla.limpo)

head(coord.cec.gla.limpo.names)
is.data.frame(coord.cec.gla.limpo.names)
colnames(coord.cec.gla.limpo.names)<-c("SPEC", "LONG", "LAT")

coord.cec.gla.names_thinned_dataset_full <-
  thin( loc.data =coord.cec.gla.limpo.names, 
        lat.col = "LAT", 
        long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 26.27, reps =1, #10 min resolution
        locs.thinned.list.return = TRUE, 
        write.files = T, 
        max.files = 5, 
        out.dir = "coord.cec.gla.limpo.names_full/", out.base = "coord.cec.gla.limpo.names_thinned", 
        write.log.file = TRUE,
        log.file = "coord.cec.gla.limpo.names_thinned_full_log_file.txt" )




dim(coord.cec.gla.names_thinned_dataset_full[[1]])



###ENVIRONMENTAL BACKGROUND

reg.variaveis.cec.gla= crop(variaveis, poligono_cec.gla)
plot(reg.variaveis.cec.gla, 1)
points(coord.cec.gla.limpo, pch=20, cex=.5)

#remoção variaveis nao colineares
vari.n.colineares.cec.gla<-vifstep(reg.variaveis.cec.gla)
vari.n.colineares.cec.gla


reg.predictors.cec.gla<- stack(reg.variaveis.cec.gla[[2]],
                               reg.variaveis.cec.gla[[3]],
                               reg.variaveis.cec.gla[[8]],
                               reg.variaveis.cec.gla[[9]],
                               reg.variaveis.cec.gla[[10]],
                               reg.variaveis.cec.gla[[12]],
                               reg.variaveis.cec.gla[[14]],
                               reg.variaveis.cec.gla[[18]],
                               reg.variaveis.cec.gla[[19]])



names(reg.predictors.cec.gla)


###MATA ATLANTICA projection presente

ma.presente = crop(variaveis, poligono_mata_atlantica)
plot(ma.presente, 1)
polygon(poligono_mata_atlantica)
ma.presente.poly <- mask(ma.presente, data.shape)
plot(ma.presente.poly, 1 )


ma.presente.poly.cec.gla<-stack(ma.presente.new[[2]],
                                ma.presente.new[[3]],
                                ma.presente.new[[8]],
                                ma.presente.new[[9]],
                                ma.presente.new[[10]],
                                ma.presente.new[[12]],
                                ma.presente.new[[14]],
                                ma.presente.new[[18]],
                                ma.presente.new[[19]])




names(ma.presente.poly.cec.gla)


myRespName.cec.gla <-  'cec.gla'
myResp.cec.gla<-rep(1, nrow(coord.cec.gla.names_thinned_dataset_full[[1]]))
myResp.cec.gla
length(myResp.cec.gla)
myRespXY.cec.gla<-coord.cec.gla.names_thinned_dataset_full[[1]]

# function to define the intersect of rasters
intersect_mask <- function(x){
  values_x <- getValues(x)
  inter_x <- values_x %*% rep(1,nlayers(x))
  mask <- setValues(subset(x,1),values = (inter_x>0))
  return(mask)
}
# keep only all cells that are defined for all layers

myExpl.cec.gla<-stack(reg.predictors.cec.gla)

myExpl.cec.gla <- stack(mask(myExpl.cec.gla, intersect_mask(myExpl.cec.gla)))


##Number of pseudo-absence: Lobo & Tognelli sugerem 100 pontos de pseudo-
#ausencia para cada ponto de presenca.



myBiomodData.cec.gla<- BIOMOD_FormatingData(resp.var = myResp.cec.gla,
                                            expl.var = myExpl.cec.gla,
                                            resp.xy = myRespXY.cec.gla,
                                            resp.name = myRespName.cec.gla,
                                            PA.nb.rep = 1,
                                            PA.nb.absences = 10000,
                                            PA.strategy = 'disk',
                                            PA.dist.max=500000) # com base em http://www.sciencedirect.com/science/article/pii/S0304380008005486
myBiomodData.cec.gla
save.image()
###MODELING

# 1. Defining Models Options using default options.
myBiomodOption.cec.gla <- BIOMOD_ModelingOptions()

# 3. Computing the models

myBiomodModelOut.cec.gla<- BIOMOD_Modeling(myBiomodData.cec.gla,
                                           models = c("MAXENT.Phillips"),
                                           models.options = myBiomodOption.cec.gla,
                                           NbRunEval=5, #conforme http://onlinelibrary.wiley.com/doi/10.1111/j.1523-1739.2006.00354.x/full
                                           DataSplit=70,
                                           Prevalence=0.5,
                                           VarImport=1,
                                           models.eval.meth = c('TSS'),
                                           cec.glaveObj = TRUE,
                                           rescal.all.models = TRUE,
                                           do.full.models = FALSE,
                                           modeling.id = paste(myRespName.cec.gla,"FirstModeling",sep=""))
(myBiomodModelOut.cec.gla)
save.image()
# get all models evaluation

#Excelente TSS>0.75
#Bom 0.40>=TSS<0.75
#Ruim TSS<0.40


myBiomodModelEval.cec.gla <- get_evaluations(myBiomodModelOut.cec.gla)

myBiomodModelEval.cec.gla




get_variables_importance(myBiomodModelOut.cec.gla)
get_evaluations_matrix_cec.gla<-cbind(get_variables_importance(myBiomodModelOut.cec.gla)[1:9], #RUN 1
                                      get_variables_importance(myBiomodModelOut.cec.gla)[10:18], #RUN 2
                                      get_variables_importance(myBiomodModelOut.cec.gla)[19:27], #RUN 3
                                      get_variables_importance(myBiomodModelOut.cec.gla)[28:36], #RUN 4
                                      get_variables_importance(myBiomodModelOut.cec.gla)[37:45]) #RUN 5
colnames(get_evaluations_matrix_cec.gla)<-c("RUN1", "RUN2", "RUN3", "RUN4", "RUN5")
names(reg.predictors.cec.gla)
rownames(get_evaluations_matrix_cec.gla)<-names(reg.predictors.cec.gla)
get_evaluations_matrix_cec.gla_mean<-apply(get_evaluations_matrix_cec.gla, 1, mean)
which.max(get_evaluations_matrix_cec.gla_mean)



mymaxent.cec.gla <- BIOMOD_LoadModels(myBiomodModelOut.cec.gla, models='MAXENT.Phillips')
myRespPlot2D.cec.gla <- response.plot2(models  = mymaxent.cec.gla,
                                       Data = get_formal_data(myBiomodModelOut.cec.gla,'expl.var'), 
                                       show.variables= get_formal_data(myBiomodModelOut.cec.gla,'expl.var.names'),
                                       do.bivariate = FALSE,
                                       fixed.var.metric = 'mean',
                                       col = c("blue", "red"),
                                       legend = F,
                                       data_species = get_formal_data(myBiomodModelOut.cec.gla,'resp.var'),
                                       plot=F)

save.image()

bio_3_response_cec.gla<-gather(myRespPlot2D.cec.gla$bio3,model,Probability,
                               cec.gla_PA1_RUN1_MAXENT.Phillips:cec.gla_PA1_RUN5_MAXENT.Phillips)



p_cec.gla_bio3<-ggplot(data =bio_3_response_cec.gla) + 
  geom_smooth(mapping = aes(x = bio3, y = Probability))+
  labs(x = bquote('Bio 3'))
p_cec.gla_bio3



###Ensemble modeling
myBiomodEM.cec.gla<-BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut.cec.gla,
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

myBiomodProjection.cec.gla <- BIOMOD_Projection(modeling.output = myBiomodModelOut.cec.gla,
                                                new.env = stack(ma.presente.poly.cec.gla),
                                                proj.name = 'current_MA',
                                                selected.models = "all",
                                                binary.meth = 'TSS',
                                                compress = FALSE,
                                                build.clamping.mask = FALSE)

#Ensemble Forcasting

myBiomodEF.cec.gla.presente <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM.cec.gla,
  projection.output = myBiomodProjection.cec.gla,binary.meth = "TSS")

save.image()

plot(myBiomodEF.cec.gla.presente)

currentPred.ensemble.cec.gla <- stack("cec.gla/proj_current_MA/proj_current_MA_cec.gla_ensemble.grd")
plot(currentPred.ensemble.cec.gla[[1]])

plot(currentPred.ensemble.cec.gla[[1]],col = gray.colors(12, start = 0.9, end = 0.0))
plot(wrld_simpl, xlim=c(-60,-40), ylim=c(-34,-15), axes=TRUE, add=T)
polygon(poligono_mata_atlantica)


plot(currentPred.ensemble.cec.gla[[2]],col = gray.colors(12, start = 0.9, end = 0.0))
points(coord.cec.gla.names_thinned_dataset_full[[1]], cex=1,col="black", pch=20)
polygon(poligono_mata_atlantica)


?BIEN_metadata_citation
library(BIEN)
BIEN_metadata_citation()
cec.gla.bien
citations.cec.gla<-BIEN_metadata_citation(dataframe=cec.gla.bien, acknowledgement_file = "citations_cec_gla.txt")#If you are referencing occurrence data
citations.cec.gla


citations.alc.tri
Xanthium_data<-BIEN_occurrence_species("Xanthium strumarium")
citations.teste<-BIEN_metadata_citation(dataframe=Xanthium_data)#If you are referencing occurrence data
citations.teste

########################## 2070 ##########################
################################################################



ma.2070.poly.cec.gla<-stack(ma.2070.new[[2]],
                                ma.2070.new[[3]],
                                ma.2070.new[[8]],
                                ma.2070.new[[9]],
                                ma.2070.new[[10]],
                                ma.2070.new[[12]],
                                ma.2070.new[[14]],
                                ma.2070.new[[18]],
                                ma.2070.new[[19]])


names(ma.2070.poly.cec.gla)<-names(ma.presente.poly.cec.gla)


myExplFuture.2070.85.ma<-stack(ma.2070.poly.cec.gla) #modelo HE


projection.2070.85.ma.cec.gla <-BIOMOD_Projection(modeling.output=myBiomodModelOut.cec.gla,
                                                  new.env= myExplFuture.2070.85.ma,
                                                  proj.name='cec.gla.2070.85',
                                                  selected.models =  "all",
                                                  Bin.trans=TRUE,
                                                  binary.meth ='TSS',
                                                  compress = 'gzip',
                                                  build.clamping.mask	= FALSE,
                                                  SaveObj=TRUE,
                                                  keep.in.memory=F,
                                                  do.stack=T,
                                                  silent=F)


cec.gla.2070.85<-get_predictions(projection.2070.85.ma.cec.gla)
cec.gla.2070.85
summary(cec.gla.2070.85)





################## ENSAMBLE MODELLING 2070  ################
#################################################################



EnsambleForecast.2070.85.ma.cec.gla <-BIOMOD_EnsembleForecasting (
  projection.output = projection.2070.85.ma.cec.gla,
  EM.output=myBiomodEM.cec.gla,
  selected.models = "all",
  binary.meth ='TSS') 




EnsambleForecast.2070.85.ma.cec.gla
plot(EnsambleForecast.2070.85.ma.cec.gla)



ensemble.2070.HE.85.ma.cec.gla <- stack("cec.gla/proj_cec.gla.2070.85/proj_cec.gla.2070.85_cec.gla_ensemble.grd")
ensemble.2070.HE.85.ma.cec.gla



plot(ensemble.2070.HE.85.ma.cec.gla[[1]],col = gray.colors(12, start = 0.9, end = 0))
points(coord.gym_klo.names_thinned_dataset_full[[1]],pch=20, col="black", cex=.7)

plot(ensemble.2070.HE.85.ma.cec.gla[[2]],col = gray.colors(12, start = 0.9, end = 0))



futurePred.ensemble.2070.85.cec.gla <- stack("cec.gla/proj_cec.gla.2070.85/proj_cec.gla.2070.85_cec.gla_ensemble.grd")



plot(currentPred.ensemble.cec.gla[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.2070.85.cec.gla[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS


currentPred.ensemble.bin.cec.gla<- stack("cec.gla/proj_current_MA/proj_current_MA_cec.gla_ensemble_TSSbin.grd")
futurePred.ensemble.bin.2070.85.cec.gla <- stack("cec.gla/proj_cec.gla.2070.85/proj_cec.gla.2070.85_cec.gla_ensemble_TSSbin.gri")


plot(currentPred.ensemble.bin.cec.gla[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.bin.2070.85.cec.gla[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS



myBiomodRangeSize.2070.85.cec.gla <- BIOMOD_RangeSize(
  CurrentPred=currentPred.ensemble.bin.cec.gla[[1]],
  FutureProj=futurePred.ensemble.bin.2070.85.cec.gla[[1]])
myBiomodRangeSize.2070.85.cec.gla$Compt.By.Models





plot(myBiomodRangeSize.2070.85.cec.gla$Diff.By.Pixel, col= gray.colors(4, start = 0, end = 1, gamma = 1, alpha = NULL))
plot(wrld_simpl,  axes="TRUE", add=T)

polygon(poligono_mata_atlantica)






