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
points(coord.cal.bra, pch=19, cex=.51, col="black")


##Poligono MA

require(rgdal)
data.shape<-readOGR(dsn="/media/dendrologia/423EF1863EF172F1/User/Documents/orientados/Guilherme/IDS_09_Desm_Bioma_Mata_Atlantica_1_2", 
                    layer="IDS_09_Desm_Bioma_Mata_Atlantica_1_2")

poligono_mata_atlantica<-data.shape@polygons[[1]]@Polygons[[1]]@coords

polygon(poligono_mata_atlantica)
###Poligono especie


poligono_cal.bra<-getpoly()

colnames(poligono_cal.bra) = c("Long", "Lat")
polygon(poligono_cal.bra) 



########REMOCAO DUVIDAS

out.cal.bra= pnt.in.poly(coord.cal.bra,poligono_cal.bra)
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
plot(data.shape, add=T, col="lightgreen")
points(out.cal.bra[which(out.cal.bra$pip==1),1:2],pch=20, cex=.5, col="black")
coord.cal.bra.limpo<-out.cal.bra[which(out.cal.bra$pip==1),1:2]
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
points(coord.cal.bra.limpo,pch=20, cex=.5, col="black")


###FILTRAGEM ESPACIAl

head(coord.cal.bra.limpo)
coord.cal.bra.limpo<-na.omit(coord.cal.bra.limpo)

dim(coord.cal.bra.limpo)
cal.bra_names<-c(rep("cal.bra",  dim(coord.cal.bra.limpo)[1]))
coord.cal.bra.limpo.names<-cbind(cal.bra_names, coord.cal.bra.limpo)

head(coord.cal.bra.limpo.names)
is.data.frame(coord.cal.bra.limpo.names)
colnames(coord.cal.bra.limpo.names)<-c("SPEC", "LONG", "LAT")
# 
# ?write.table
# write.table(coord.cal.bra.limpo.names, "tap_gui.csv", sep=";", dec=",")

coord.cal.bra.names_thinned_dataset_full <-
  thin( loc.data =coord.cal.bra.limpo.names, 
        lat.col = "LAT", 
        long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 26.27, reps =1, #10 min resolution
        locs.thinned.list.return = TRUE, 
        write.files = T, 
        max.files = 5, 
        out.dir = "coord.cal.bra.limpo.names_full/", out.base = "coord.cal.bra.limpo.names_thinned", 
        write.log.file = TRUE,
        log.file = "coord.cal.bra.limpo.names_thinned_full_log_file.txt" )



dim(coord.cal.bra.limpo.names)
dim(coord.cal.bra.names_thinned_dataset_full[[1]])



###ENVIRONMENTAL BACKGROUND

reg.variaveis.cal.bra= crop(variaveis, poligono_cal.bra)
plot(reg.variaveis.cal.bra, 1)
points(coord.cal.bra.limpo, pch=20, cex=.5)

#remoção variaveis nao colineares
vari.n.colineares.cal.bra<-vifstep(reg.variaveis.cal.bra)
vari.n.colineares.cal.bra


reg.predictors.cal.bra<- stack(reg.variaveis.cal.bra[[2]],
                               reg.variaveis.cal.bra[[3]],
                               reg.variaveis.cal.bra[[4]],
                               reg.variaveis.cal.bra[[8]],
                               reg.variaveis.cal.bra[[13]],
                               reg.variaveis.cal.bra[[14]],
                               reg.variaveis.cal.bra[[15]],
                               reg.variaveis.cal.bra[[18]],
                               reg.variaveis.cal.bra[[19]])



names(reg.predictors.cal.bra)


###MATA ATLANTICA projection presente

ma.presente = crop(variaveis, poligono_mata_atlantica)
plot(ma.presente, 1)
polygon(poligono_mata_atlantica)
ma.presente.poly <- mask(ma.presente, data.shape)
plot(ma.presente.poly, 1 )


ma.presente.poly.cal.bra<-stack(ma.presente.poly[[2]],
                                ma.presente.poly[[3]],
                                ma.presente.poly[[4]],
                                ma.presente.poly[[8]],
                                ma.presente.poly[[13]],
                                ma.presente.poly[[14]],
                                ma.presente.poly[[15]],
                                ma.presente.poly[[18]],
                                ma.presente.poly[[19]])




names(ma.presente.poly.cal.bra)


myRespName.cal.bra <-  'cal.bra'
myResp.cal.bra<-rep(1, nrow(coord.cal.bra.names_thinned_dataset_full[[1]]))
myResp.cal.bra
length(myResp.cal.bra)
myRespXY.cal.bra<-coord.cal.bra.names_thinned_dataset_full[[1]]

# function to define the intersect of rasters
intersect_mask <- function(x){
  values_x <- getValues(x)
  inter_x <- values_x %*% rep(1,nlayers(x))
  mask <- setValues(subset(x,1),values = (inter_x>0))
  return(mask)
}
# keep only all cells that are defined for all layers

myExpl.cal.bra<-stack(reg.predictors.cal.bra)

myExpl.cal.bra <- stack(mask(myExpl.cal.bra, intersect_mask(myExpl.cal.bra)))


##Number of pseudo-absence: Lobo & Tognelli sugerem 100 pontos de pseudo-
#ausencia para cada ponto de presenca.



myBiomodData.cal.bra<- BIOMOD_FormatingData(resp.var = myResp.cal.bra,
                                            expl.var = myExpl.cal.bra,
                                            resp.xy = myRespXY.cal.bra,
                                            resp.name = myRespName.cal.bra,
                                            PA.nb.rep = 1,
                                            PA.nb.absences = 10000,
                                            PA.strategy = 'disk',
                                            PA.dist.max=500000) # com base em http://www.sciencedirect.com/science/article/pii/S0304380008005486
myBiomodData.cal.bra
save.image()
###MODELING

# 1. Defining Models Options using default options.
myBiomodOption.cal.bra <- BIOMOD_ModelingOptions()

# 3. Computing the models

myBiomodModelOut.cal.bra<- BIOMOD_Modeling(myBiomodData.cal.bra,
                                           models = c("MAXENT.Phillips"),
                                           models.options = myBiomodOption.cal.bra,
                                           NbRunEval=5, #conforme http://onlinelibrary.wiley.com/doi/10.1111/j.1523-1739.2006.00354.x/full
                                           DataSplit=70,
                                           Prevalence=0.5,
                                           VarImport=1,
                                           models.eval.meth = c('TSS'),
                                           cal.braveObj = TRUE,
                                           rescal.all.models = TRUE,
                                           do.full.models = FALSE,
                                           modeling.id = paste(myRespName.cal.bra,"FirstModeling",sep=""))
(myBiomodModelOut.cal.bra)
save.image()
# get all models evaluation

#Excelente TSS>0.75
#Bom 0.40>=TSS<0.75
#Ruim TSS<0.40


myBiomodModelEval.cal.bra <- get_evaluations(myBiomodModelOut.cal.bra)

myBiomodModelEval.cal.bra


#####STOP NON TSS > 0.4




get_variables_importance(myBiomodModelOut.cal.bra)
get_evaluations_matrix_cal.bra<-cbind(get_variables_importance(myBiomodModelOut.cal.bra)[1:9], #RUN 1
                                      get_variables_importance(myBiomodModelOut.cal.bra)[10:18], #RUN 2
                                      get_variables_importance(myBiomodModelOut.cal.bra)[19:27], #RUN 3
                                      get_variables_importance(myBiomodModelOut.cal.bra)[28:36], #RUN 4
                                      get_variables_importance(myBiomodModelOut.cal.bra)[37:45]) #RUN 5
colnames(get_evaluations_matrix_cal.bra)<-c("RUN1", "RUN2", "RUN3", "RUN4", "RUN5")
names(reg.predictors.cal.bra)
rownames(get_evaluations_matrix_cal.bra)<-names(reg.predictors.cal.bra)
get_evaluations_matrix_cal.bra_mean<-apply(get_evaluations_matrix_cal.bra, 1, mean)
which.max(get_evaluations_matrix_cal.bra_mean)



mymaxent.cal.bra <- BIOMOD_LoadModels(myBiomodModelOut.cal.bra, models='MAXENT.Phillips')
myRespPlot2D.cal.bra <- response.plot2(models  = mymaxent.cal.bra,
                                       Data = get_formal_data(myBiomodModelOut.cal.bra,'expl.var'), 
                                       show.variables= get_formal_data(myBiomodModelOut.cal.bra,'expl.var.names'),
                                       do.bivariate = FALSE,
                                       fixed.var.metric = 'mean',
                                       col = c("blue", "red"),
                                       legend = F,
                                       data_species = get_formal_data(myBiomodModelOut.cal.bra,'resp.var'),
                                       plot=F)

save.image()

bio_3_response_cal.bra<-gather(myRespPlot2D.cal.bra$bio3,model,Probability,
                               cal.bra_PA1_RUN1_MAXENT.Phillips:cal.bra_PA1_RUN5_MAXENT.Phillips)



p_cal.bra_bio3<-ggplot(data =bio_3_response_cal.bra) + 
  geom_smooth(mapping = aes(x = bio3, y = Probability))+
  labs(x = bquote('Isothermality (%)'))
p_cal.bra_bio3



###Ensemble modeling
myBiomodEM.cal.bra<-BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut.cal.bra,
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

myBiomodProjection.cal.bra <- BIOMOD_Projection(modeling.output = myBiomodModelOut.cal.bra,
                                                new.env = stack(ma.presente.poly.cal.bra),
                                                proj.name = 'current_MA',
                                                selected.models = "all",
                                                binary.meth = 'TSS',
                                                compress = FALSE,
                                                build.clamping.mask = FALSE)

#Ensemble Forcasting

myBiomodEF.cal.bra.presente <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM.cal.bra,
  projection.output = myBiomodProjection.cal.bra,binary.meth = "TSS")

save.image()

plot(myBiomodEF.cal.bra.presente)

currentPred.ensemble.cal.bra <- stack("cal.bra/proj_current_MA/proj_current_MA_cal.bra_ensemble.grd")
plot(currentPred.ensemble.cal.bra[[1]])

plot(currentPred.ensemble.cal.bra[[1]],col = gray.colors(12, start = 0.9, end = 0.0))
plot(wrld_simpl, xlim=c(-60,-40), ylim=c(-34,-15), axes=TRUE, add=T)
polygon(poligono_mata_atlantica)
points(coord.cal.bra.names_thinned_dataset_full[[1]], cex=.6,col="black")

plot(currentPred.ensemble.cal.bra[[2]],col = gray.colors(12, start = 0.9, end = 0.0))
points(coord.cal.bra.names_thinned_dataset_full[[1]], cex=.6,col="black")
polygon(poligono_mata_atlantica)


?BIEN_metadata_citation

citations.cal.bra<-BIEN_metadata_citation(dataframe=cal.bra.bien)#If you are referencing occurrence data
citations.cal.bra


########################## 2070 ##########################
################################################################


#NO TSS > 0.70
