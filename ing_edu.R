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
points(coord.ing.edu, pch=19, cex=.51, col="black")


##Poligono MA

require(rgdal)
data.shape<-readOGR(dsn="/media/dendrologia/423EF1863EF172F1/User/Documents/orientados/Guilherme/IDS_09_Desm_Bioma_Mata_Atlantica_1_2", 
                    layer="IDS_09_Desm_Bioma_Mata_Atlantica_1_2")

poligono_mata_atlantica<-data.shape@polygons[[1]]@Polygons[[1]]@coords

polygon(poligono_mata_atlantica)
###Poligono especie


poligono_ing.edu<-getpoly()

colnames(poligono_ing.edu) = c("Long", "Lat")
polygon(poligono_ing.edu) 



########REMOCAO DUVIDAS

out.ing.edu= pnt.in.poly(coord.ing.edu,poligono_ing.edu)
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
plot(data.shape, add=T, col="lightgreen")
points(out.ing.edu[which(out.ing.edu$pip==1),1:2],pch=20, cex=.5, col="black")
coord.ing.edu.limpo<-out.ing.edu[which(out.ing.edu$pip==1),1:2]
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
points(coord.ing.edu.limpo,pch=20, cex=.5, col="black")


###FILTRAGEM ESPACIAl

head(coord.ing.edu.limpo)
coord.ing.edu.limpo<-na.omit(coord.ing.edu.limpo)

dim(coord.ing.edu.limpo)
ing.edu_names<-c(rep("ing.edu",  dim(coord.ing.edu.limpo)[1]))
coord.ing.edu.limpo.names<-cbind(ing.edu_names, coord.ing.edu.limpo)

head(coord.ing.edu.limpo.names)
is.data.frame(coord.ing.edu.limpo.names)
colnames(coord.ing.edu.limpo.names)<-c("SPEC", "LONG", "LAT")

coord.ing.edu.names_thinned_dataset_full <-
  thin( loc.data =coord.ing.edu.limpo.names, 
        lat.col = "LAT", 
        long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 26.27, reps =1, #10 min resolution
        locs.thinned.list.return = TRUE, 
        write.files = T, 
        max.files = 5, 
        out.dir = "coord.ing.edu.limpo.names_full/", out.base = "coord.ing.edu.limpo.names_thinned", 
        write.log.file = TRUE,
        log.file = "coord.ing.edu.limpo.names_thinned_full_log_file.txt" )




dim(coord.ing.edu.names_thinned_dataset_full[[1]])



###ENVIRONMENTAL BACKGROUND

reg.variaveis.ing.edu= crop(variaveis, poligono_ing.edu)
plot(reg.variaveis.ing.edu, 1)
points(coord.ing.edu.limpo, pch=20, cex=.5)

#remoção variaveis nao colineares
vari.n.colineares.ing.edu<-vifstep(reg.variaveis.ing.edu)
vari.n.colineares.ing.edu


reg.predictors.ing.edu<- stack(reg.variaveis.ing.edu[[2]],
                               reg.variaveis.ing.edu[[3]],
                               reg.variaveis.ing.edu[[8]],
                               reg.variaveis.ing.edu[[9]],
                               reg.variaveis.ing.edu[[13]],
                               reg.variaveis.ing.edu[[14]],
                               reg.variaveis.ing.edu[[15]],
                               reg.variaveis.ing.edu[[18]],
                               reg.variaveis.ing.edu[[19]])



names(reg.predictors.ing.edu)


###MATA ATLANTICA projection presente

ma.presente = crop(variaveis, poligono_mata_atlantica)
plot(ma.presente, 1)
polygon(poligono_mata_atlantica)
ma.presente.poly <- mask(ma.presente, data.shape)
plot(ma.presente.poly, 1 )


ma.presente.poly.ing.edu<-stack(ma.presente.poly[[2]],
                                ma.presente.poly[[3]],
                                ma.presente.poly[[8]],
                                ma.presente.poly[[9]],
                                ma.presente.poly[[13]],
                                ma.presente.poly[[14]],
                                ma.presente.poly[[15]],
                                ma.presente.poly[[18]],
                                ma.presente.poly[[19]])




names(ma.presente.poly.ing.edu)


myRespName.ing.edu <-  'ing.edu'
myResp.ing.edu<-rep(1, nrow(coord.ing.edu.names_thinned_dataset_full[[1]]))
myResp.ing.edu
length(myResp.ing.edu)
myRespXY.ing.edu<-coord.ing.edu.names_thinned_dataset_full[[1]]

# function to define the intersect of rasters
intersect_mask <- function(x){
  values_x <- getValues(x)
  inter_x <- values_x %*% rep(1,nlayers(x))
  mask <- setValues(subset(x,1),values = (inter_x>0))
  return(mask)
}
# keep only all cells that are defined for all layers

myExpl.ing.edu<-stack(reg.predictors.ing.edu)

myExpl.ing.edu <- stack(mask(myExpl.ing.edu, intersect_mask(myExpl.ing.edu)))


##Number of pseudo-absence: Lobo & Tognelli sugerem 100 pontos de pseudo-
#ausencia para cada ponto de presenca.



myBiomodData.ing.edu<- BIOMOD_FormatingData(resp.var = myResp.ing.edu,
                                            expl.var = myExpl.ing.edu,
                                            resp.xy = myRespXY.ing.edu,
                                            resp.name = myRespName.ing.edu,
                                            PA.nb.rep = 1,
                                            PA.nb.absences = 10000,
                                            PA.strategy = 'disk',
                                            PA.dist.max=500000) # com base em http://www.sciencedirect.com/science/article/pii/S0304380008005486
myBiomodData.ing.edu
save.image()
###MODELING

# 1. Defining Models Options using default options.
myBiomodOption.ing.edu <- BIOMOD_ModelingOptions()

# 3. Computing the models

myBiomodModelOut.ing.edu<- BIOMOD_Modeling(myBiomodData.ing.edu,
                                           models = c("MAXENT.Phillips"),
                                           models.options = myBiomodOption.ing.edu,
                                           NbRunEval=5, #conforme http://onlinelibrary.wiley.com/doi/10.1111/j.1523-1739.2006.00354.x/full
                                           DataSplit=70,
                                           Prevalence=0.5,
                                           VarImport=1,
                                           models.eval.meth = c('TSS'),
                                           ing.eduveObj = TRUE,
                                           rescal.all.models = TRUE,
                                           do.full.models = FALSE,
                                           modeling.id = paste(myRespName.ing.edu,"FirstModeling",sep=""))
(myBiomodModelOut.ing.edu)
save.image()
# get all models evaluation

#Excelente TSS>0.75
#Bom 0.40>=TSS<0.75
#Ruim TSS<0.40


myBiomodModelEval.ing.edu <- get_evaluations(myBiomodModelOut.ing.edu)

myBiomodModelEval.ing.edu

####STOP - NO TSS > 0.4


get_variables_importance(myBiomodModelOut.ing.edu)
get_evaluations_matrix_ing.edu<-cbind(get_variables_importance(myBiomodModelOut.ing.edu)[1:9], #RUN 1
                                      get_variables_importance(myBiomodModelOut.ing.edu)[10:18], #RUN 2
                                      get_variables_importance(myBiomodModelOut.ing.edu)[19:27], #RUN 3
                                      get_variables_importance(myBiomodModelOut.ing.edu)[28:36], #RUN 4
                                      get_variables_importance(myBiomodModelOut.ing.edu)[37:45]) #RUN 5
colnames(get_evaluations_matrix_ing.edu)<-c("RUN1", "RUN2", "RUN3", "RUN4", "RUN5")
names(reg.predictors.ing.edu)
rownames(get_evaluations_matrix_ing.edu)<-names(reg.predictors.ing.edu)
get_evaluations_matrix_ing.edu_mean<-apply(get_evaluations_matrix_ing.edu, 1, mean)
which.max(get_evaluations_matrix_ing.edu_mean)



mymaxent.ing.edu <- BIOMOD_LoadModels(myBiomodModelOut.ing.edu, models='MAXENT.Phillips')
myRespPlot2D.ing.edu <- response.plot2(models  = mymaxent.ing.edu,
                                       Data = get_formal_data(myBiomodModelOut.ing.edu,'expl.var'), 
                                       show.variables= get_formal_data(myBiomodModelOut.ing.edu,'expl.var.names'),
                                       do.bivariate = FALSE,
                                       fixed.var.metric = 'mean',
                                       col = c("blue", "red"),
                                       legend = F,
                                       data_species = get_formal_data(myBiomodModelOut.ing.edu,'resp.var'),
                                       plot=F)

save.image()

bio_14_response_ing.edu<-gather(myRespPlot2D.ing.edu$bio14,model,Probability,
                                ing.edu_PA1_RUN1_MAXENT.Phillips:ing.edu_PA1_RUN5_MAXENT.Phillips)



p_ing.edu_bio14<-ggplot(data =bio_14_response_ing.edu) + 
  geom_smooth(mapping = aes(x = bio14, y = Probability))+
  labs(x = bquote('Precipitation of Driest Month (mm)'))
p_ing.edu_bio14



###Ensemble modeling
myBiomodEM.ing.edu<-BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut.ing.edu,
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

myBiomodProjection.ing.edu <- BIOMOD_Projection(modeling.output = myBiomodModelOut.ing.edu,
                                                new.env = stack(ma.presente.poly.ing.edu),
                                                proj.name = 'current_MA',
                                                selected.models = "all",
                                                binary.meth = 'TSS',
                                                compress = FALSE,
                                                build.clamping.mask = FALSE)

#Ensemble Forcasting

myBiomodEF.ing.edu.presente <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM.ing.edu,
  projection.output = myBiomodProjection.ing.edu,binary.meth = "TSS")

save.image()

plot(myBiomodEF.ing.edu.presente)

currentPred.ensemble.ing.edu <- stack("ing.edu/proj_current_MA/proj_current_MA_ing.edu_ensemble.grd")
plot(currentPred.ensemble.ing.edu[[1]])

plot(currentPred.ensemble.ing.edu[[1]],col = gray.colors(12, start = 0.9, end = 0.0))
plot(wrld_simpl, xlim=c(-60,-40), ylim=c(-34,-15), axes=TRUE, add=T)
polygon(poligono_mata_atlantica)


plot(currentPred.ensemble.ing.edu[[2]],col = gray.colors(12, start = 0.9, end = 0.0))
points(coord.ing.edu.names_thinned_dataset_full[[1]], cex=.6,col="black")
polygon(poligono_mata_atlantica)


?BIEN_metadata_citation

citations<-BIEN_metadata_citation(dataframe=ing.edu.bien)#If you are referencing occurrence data
citations


########################## 2070 ##########################
################################################################


###NO TSS > 0.4


citations.ing.edu<-BIEN_metadata_citation(dataframe=ing.edu.bien)#If you are referencing occurrence data
citations.ing.edu
