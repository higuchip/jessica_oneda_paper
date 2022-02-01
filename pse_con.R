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
points(coord.pse.con, pch=19, cex=.51, col="black")


##Poligono MA

require(rgdal)
data.shape<-readOGR(dsn="/media/dendrologia/423EF1863EF172F1/User/Documents/orientados/Guilherme/IDS_09_Desm_Bioma_Mata_Atlantica_1_2", 
                    layer="IDS_09_Desm_Bioma_Mata_Atlantica_1_2")

poligono_mata_atlantica<-data.shape@polygons[[1]]@Polygons[[1]]@coords

polygon(poligono_mata_atlantica)
###Poligono especie


poligono_pse.con<-getpoly()

colnames(poligono_pse.con) = c("Long", "Lat")
polygon(poligono_pse.con) 



########REMOCAO DUVIDAS

out.pse.con= pnt.in.poly(coord.pse.con,poligono_pse.con)
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
plot(data.shape, add=T, col="lightgreen")
points(out.pse.con[which(out.pse.con$pip==1),1:2],pch=20, cex=.5, col="black")
coord.pse.con.limpo<-out.pse.con[which(out.pse.con$pip==1),1:2]
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
points(coord.pse.con.limpo,pch=20, cex=.5, col="black")


###FILTRAGEM ESPACIAl

head(coord.pse.con.limpo)
coord.pse.con.limpo<-na.omit(coord.pse.con.limpo)

dim(coord.pse.con.limpo)
pse.con_names<-c(rep("pse.con",  dim(coord.pse.con.limpo)[1]))
coord.pse.con.limpo.names<-cbind(pse.con_names, coord.pse.con.limpo)

head(coord.pse.con.limpo.names)
is.data.frame(coord.pse.con.limpo.names)
colnames(coord.pse.con.limpo.names)<-c("SPEC", "LONG", "LAT")
# 
# ?write.table
# write.table(coord.pse.con.limpo.names, "tap_gui.csv", sep=";", dec=",")




coord.pse.con.names_thinned_dataset_full <-
  thin( loc.data =coord.pse.con.limpo.names, 
        lat.col = "LAT", 
        long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 26.27, reps =1, #10 min resolution
        locs.thinned.list.return = TRUE, 
        write.files = T, 
        max.files = 5, 
        out.dir = "coord.pse.con.limpo.names_full/", out.base = "coord.pse.con.limpo.names_thinned", 
        write.log.file = TRUE,
        log.file = "coord.pse.con.limpo.names_thinned_full_log_file.txt" )



dim(coord.pse.con.limpo.names)
dim(coord.pse.con.names_thinned_dataset_full[[1]])



###ENVIRONMENTAL BACKGROUND

reg.variaveis.pse.con= crop(variaveis, poligono_pse.con)
plot(reg.variaveis.pse.con, 1)
points(coord.pse.con.limpo, pch=20, cex=.5)

#remoção variaveis nao colineares
vari.n.colineares.pse.con<-vifstep(reg.variaveis.pse.con)
vari.n.colineares.pse.con


reg.predictors.pse.con<- stack(reg.variaveis.pse.con[[2]],
                               reg.variaveis.pse.con[[3]],
                               reg.variaveis.pse.con[[4]],
                               reg.variaveis.pse.con[[8]],
                               reg.variaveis.pse.con[[10]],
                               reg.variaveis.pse.con[[13]],
                               reg.variaveis.pse.con[[14]],
                               reg.variaveis.pse.con[[15]],
                               reg.variaveis.pse.con[[18]],
                               reg.variaveis.pse.con[[19]])



names(reg.predictors.pse.con)


###MATA ATLANTICA projection presente

ma.presente = crop(variaveis, poligono_mata_atlantica)
plot(ma.presente, 1)
polygon(poligono_mata_atlantica)
ma.presente.poly <- mask(ma.presente, data.shape)
plot(ma.presente.poly, 1 )


ma.presente.poly.pse.con<-stack(ma.presente.new[[2]],
                                ma.presente.new[[3]],
                                ma.presente.new[[4]],
                                ma.presente.new[[8]],
                                ma.presente.new[[10]],
                                ma.presente.new[[13]],
                                ma.presente.new[[14]],
                                ma.presente.new[[15]],
                                ma.presente.new[[18]],
                                ma.presente.new[[19]])




names(ma.presente.poly.pse.con)


myRespName.pse.con <-  'pse.con'
myResp.pse.con<-rep(1, nrow(coord.pse.con.names_thinned_dataset_full[[1]]))
myResp.pse.con
length(myResp.pse.con)
myRespXY.pse.con<-coord.pse.con.names_thinned_dataset_full[[1]]

# function to define the intersect of rasters
intersect_mask <- function(x){
  values_x <- getValues(x)
  inter_x <- values_x %*% rep(1,nlayers(x))
  mask <- setValues(subset(x,1),values = (inter_x>0))
  return(mask)
}
# keep only all cells that are defined for all layers

myExpl.pse.con<-stack(reg.predictors.pse.con)

myExpl.pse.con <- stack(mask(myExpl.pse.con, intersect_mask(myExpl.pse.con)))


##Number of pseudo-absence: Lobo & Tognelli sugerem 100 pontos de pseudo-
#ausencia para cada ponto de presenca.



myBiomodData.pse.con<- BIOMOD_FormatingData(resp.var = myResp.pse.con,
                                            expl.var = myExpl.pse.con,
                                            resp.xy = myRespXY.pse.con,
                                            resp.name = myRespName.pse.con,
                                            PA.nb.rep = 1,
                                            PA.nb.absences = 10000,
                                            PA.strategy = 'disk',
                                            PA.dist.max=500000) # com base em http://www.sciencedirect.com/science/article/pii/S0304380008005486
myBiomodData.pse.con
save.image()
###MODELING

# 1. Defining Models Options using default options.
myBiomodOption.pse.con <- BIOMOD_ModelingOptions()

# 3. Computing the models

myBiomodModelOut.pse.con<- BIOMOD_Modeling(myBiomodData.pse.con,
                                           models = c("MAXENT.Phillips"),
                                           models.options = myBiomodOption.pse.con,
                                           NbRunEval=5, #conforme http://onlinelibrary.wiley.com/doi/10.1111/j.1523-1739.2006.00354.x/full
                                           DataSplit=70,
                                           Prevalence=0.5,
                                           VarImport=1,
                                           models.eval.meth = c('TSS'),
                                           pse.conveObj = TRUE,
                                           rescal.all.models = TRUE,
                                           do.full.models = FALSE,
                                           modeling.id = paste(myRespName.pse.con,"FirstModeling",sep=""))
(myBiomodModelOut.pse.con)
save.image()
# get all models evaluation

#Excelente TSS>0.75
#Bom 0.40>=TSS<0.75
#Ruim TSS<0.40


myBiomodModelEval.pse.con <- get_evaluations(myBiomodModelOut.pse.con)

myBiomodModelEval.pse.con






get_variables_importance(myBiomodModelOut.pse.con)
get_evaluations_matrix_pse.con<-cbind(get_variables_importance(myBiomodModelOut.pse.con)[1:10], #RUN 1
                                      get_variables_importance(myBiomodModelOut.pse.con)[11:20], #RUN 2
                                      get_variables_importance(myBiomodModelOut.pse.con)[21:30], #RUN 3
                                      get_variables_importance(myBiomodModelOut.pse.con)[31:40], #RUN 4
                                      get_variables_importance(myBiomodModelOut.pse.con)[41:50]) #RUN 5
colnames(get_evaluations_matrix_pse.con)<-c("RUN1", "RUN2", "RUN3", "RUN4", "RUN5")
names(reg.predictors.pse.con)
rownames(get_evaluations_matrix_pse.con)<-names(reg.predictors.pse.con)
get_evaluations_matrix_pse.con_mean<-apply(get_evaluations_matrix_pse.con, 1, mean)
which.max(get_evaluations_matrix_pse.con_mean)



mymaxent.pse.con <- BIOMOD_LoadModels(myBiomodModelOut.pse.con, models='MAXENT.Phillips')
myRespPlot2D.pse.con <- response.plot2(models  = mymaxent.pse.con,
                                       Data = get_formal_data(myBiomodModelOut.pse.con,'expl.var'), 
                                       show.variables= get_formal_data(myBiomodModelOut.pse.con,'expl.var.names'),
                                       do.bivariate = FALSE,
                                       fixed.var.metric = 'mean',
                                       col = c("blue", "red"),
                                       legend = F,
                                       data_species = get_formal_data(myBiomodModelOut.pse.con,'resp.var'),
                                       plot=F)

save.image()

bio_14_response_pse.con<-gather(myRespPlot2D.pse.con$bio14,model,Probability,
                               pse.con_PA1_RUN1_MAXENT.Phillips:pse.con_PA1_RUN5_MAXENT.Phillips)



p_pse.con_bio14<-ggplot(data =bio_14_response_pse.con) + 
  geom_smooth(mapping = aes(x = bio14, y = Probability))+
  labs(x = bquote('Bio 14'))
p_pse.con_bio14



###Ensemble modeling
myBiomodEM.pse.con<-BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut.pse.con,
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

myBiomodProjection.pse.con <- BIOMOD_Projection(modeling.output = myBiomodModelOut.pse.con,
                                                new.env = stack(ma.presente.poly.pse.con),
                                                proj.name = 'current_MA',
                                                selected.models = "all",
                                                binary.meth = 'TSS',
                                                compress = FALSE,
                                                build.clamping.mask = FALSE)

#Ensemble Forcasting

myBiomodEF.pse.con.presente <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM.pse.con,
  projection.output = myBiomodProjection.pse.con,binary.meth = "TSS")

save.image()

plot(myBiomodEF.pse.con.presente)

currentPred.ensemble.pse.con <- stack("pse.con/proj_current_MA/proj_current_MA_pse.con_ensemble.grd")
plot(currentPred.ensemble.pse.con[[1]])

plot(currentPred.ensemble.pse.con[[1]],col = gray.colors(12, start = 0.9, end = 0.0))
plot(wrld_simpl, xlim=c(-60,-40), ylim=c(-34,-15), axes=TRUE, add=T)
polygon(poligono_mata_atlantica)
points(coord.pse.con.names_thinned_dataset_full[[1]], cex=.6,col="black")

plot(currentPred.ensemble.pse.con[[2]],col = gray.colors(12, start = 0.9, end = 0.0))
points(coord.pse.con.names_thinned_dataset_full[[1]], cex=.6,col="black")
polygon(poligono_mata_atlantica)


?BIEN_metadata_citation

citations.pse.con<-BIEN_metadata_citation(dataframe=pse.con.bien)#If you are referencing occurrence data
citations.pse.con



########################## 2070 ##########################
################################################################



ma.2070.poly.pse.con<-stack(ma.2070.new[[2]],
                                ma.2070.new[[3]],
                                ma.2070.new[[4]],
                                ma.2070.new[[8]],
                                ma.2070.new[[10]],
                                ma.2070.new[[13]],
                                ma.2070.new[[14]],
                                ma.2070.new[[15]],
                                ma.2070.new[[18]],
                                ma.2070.new[[19]])





names(ma.2070.poly.pse.con)<-names(ma.presente.poly.pse.con)


myExplFuture.2070.85.ma<-stack(ma.2070.poly.pse.con) #modelo HE


projection.2070.85.ma.pse.con <-BIOMOD_Projection(modeling.output=myBiomodModelOut.pse.con,
                                                  new.env= myExplFuture.2070.85.ma,
                                                  proj.name='pse.con.2070.85',
                                                  selected.models =  "all",
                                                  Bin.trans=TRUE,
                                                  binary.meth ='TSS',
                                                  compress = 'gzip',
                                                  build.clamping.mask	= FALSE,
                                                  SaveObj=TRUE,
                                                  keep.in.memory=F,
                                                  do.stack=T,
                                                  silent=F)


pse.con.2070.85<-get_predictions(projection.2070.85.ma.pse.con)
pse.con.2070.85
summary(pse.con.2070.85)





################## ENSAMBLE MODELLING 2070  ################
#################################################################



EnsambleForecast.2070.85.ma.pse.con <-BIOMOD_EnsembleForecasting (
  projection.output = projection.2070.85.ma.pse.con,
  EM.output=myBiomodEM.pse.con,
  selected.models = "all",
  binary.meth ='TSS') 




EnsambleForecast.2070.85.ma.pse.con
plot(EnsambleForecast.2070.85.ma.pse.con)



ensemble.2070.HE.85.ma.pse.con <- stack("pse.con/proj_pse.con.2070.85/proj_pse.con.2070.85_pse.con_ensemble.grd")
ensemble.2070.HE.85.ma.pse.con



plot(ensemble.2070.HE.85.ma.pse.con[[1]],col = gray.colors(12, start = 0.9, end = 0))
points(coord.gym_klo.names_thinned_dataset_full[[1]],pch=20, col="black", cex=.7)

plot(ensemble.2070.HE.85.ma.pse.con[[2]],col = gray.colors(12, start = 0.9, end = 0))



futurePred.ensemble.2070.85.pse.con <- stack("pse.con/proj_pse.con.2070.85/proj_pse.con.2070.85_pse.con_ensemble.grd")



plot(currentPred.ensemble.pse.con[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.2070.85.pse.con[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS


currentPred.ensemble.bin.pse.con<- stack("pse.con/proj_current_MA/proj_current_MA_pse.con_ensemble_TSSbin.grd")
futurePred.ensemble.bin.2070.85.pse.con <- stack("pse.con/proj_pse.con.2070.85/proj_pse.con.2070.85_pse.con_ensemble_TSSbin.gri")


plot(currentPred.ensemble.bin.pse.con[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.bin.2070.85.pse.con[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS



myBiomodRangeSize.2070.85.pse.con <- BIOMOD_RangeSize(
  CurrentPred=currentPred.ensemble.bin.pse.con[[1]],
  FutureProj=futurePred.ensemble.bin.2070.85.pse.con[[1]])
myBiomodRangeSize.2070.85.pse.con$Compt.By.Models





plot(myBiomodRangeSize.2070.85.pse.con$Diff.By.Pixel, col= gray.colors(4, start = 0, end = 1, gamma = 1, alpha = NULL))
plot(wrld_simpl,  axes="TRUE", add=T)

polygon(poligono_mata_atlantica)








save.image()




