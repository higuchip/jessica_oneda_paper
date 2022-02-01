
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
points(coord.all.edu, pch=19, cex=.51, col="black")


##Poligono MA

require(rgdal)
data.shape<-readOGR(dsn="/media/dendrologia/423EF1863EF172F1/User/Documents/orientados/Guilherme/IDS_09_Desm_Bioma_Mata_Atlantica_1_2", 
                    layer="IDS_09_Desm_Bioma_Mata_Atlantica_1_2")

poligono_mata_atlantica<-data.shape@polygons[[1]]@Polygons[[1]]@coords

polygon(poligono_mata_atlantica)


###Poligono especie


poligono_all.edu<-getpoly()

colnames(poligono_all.edu) = c("Long", "Lat")
polygon(poligono_all.edu) 



########REMOCAO DUVIDAS

out.all.edu= pnt.in.poly(coord.all.edu,poligono_all.edu)
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
plot(data.shape, add=T, col="lightgreen")
points(out.all.edu[which(out.all.edu$pip==1),1:2],pch=20, cex=.5, col="black")
coord.all.edu.limpo<-out.all.edu[which(out.all.edu$pip==1),1:2]
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
points(coord.all.edu.limpo,pch=20, cex=.5, col="black")


###FILTRAGEM ESPACIAl

head(coord.all.edu.limpo)
dim(coord.all.edu.limpo)
all.edu_names<-c(rep("all.edu",  dim(coord.all.edu.limpo)[1]))
coord.all.edu.limpo.names<-cbind(all.edu_names, coord.all.edu.limpo)

head(coord.all.edu.limpo.names)
is.data.frame(coord.all.edu.limpo.names)
colnames(coord.all.edu.limpo.names)<-c("SPEC", "LONG", "LAT")

coord.all.edu.names_thinned_dataset_full <-
  thin( loc.data =coord.all.edu.limpo.names, 
        lat.col = "LAT", 
        long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 26.27, reps =1, #10 min resolution
        locs.thinned.list.return = TRUE, 
        write.files = T, 
        max.files = 5, 
        out.dir = "coord.all.edu.limpo.names_full/", out.base = "coord.all.edu.limpo.names_thinned", 
        write.log.file = TRUE,
        log.file = "coord.all.edu.limpo.names_thinned_full_log_file.txt" )




dim(coord.all.edu.names_thinned_dataset_full[[1]])



###ENVIRONMENTAL BACKGROUND

reg.variaveis.all.edu= crop(variaveis, poligono_all.edu)
plot(reg.variaveis.all.edu, 1)
points(coord.all.edu.limpo, pch=20, cex=.5)

#remoção variaveis nao colineares
vari.n.colineares.all.edu<-vifstep(reg.variaveis.all.edu)
vari.n.colineares.all.edu


reg.predictors.all.edu<- stack(reg.variaveis.all.edu[[2]],
                               reg.variaveis.all.edu[[3]],
                               reg.variaveis.all.edu[[8]],
                               reg.variaveis.all.edu[[9]],
                               reg.variaveis.all.edu[[13]],
                               reg.variaveis.all.edu[[14]],
                               reg.variaveis.all.edu[[15]],
                               reg.variaveis.all.edu[[18]],
                               reg.variaveis.all.edu[[19]])



names(reg.predictors.all.edu)


###MATA ATLANTICA projection presente

ma.presente = crop(variaveis, poligono_mata_atlantica)
plot(ma.presente, 1)
polygon(poligono_mata_atlantica, add=T)
ma.presente.poly <- mask(MA.presente, data.shape)
plot(ma.presente.poly, 1 )

ma.presente.poly.all.edu <- stack(ma.presente.new[[2]],
                                  ma.presente.new[[3]],
                                  ma.presente.new[[8]],
                                  ma.presente.new[[9]],
                                  ma.presente.new[[13]],
                                  ma.presente.new[[14]],
                                  ma.presente.new[[15]],
                                  ma.presente.new[[18]],
                                  ma.presente.new[[19]])



names(ma.presente.poly.all.edu)


myRespName.all.edu <-  'all.edu'
myResp.all.edu<-rep(1, nrow(coord.all.edu.names_thinned_dataset_full[[1]]))
myResp.all.edu
length(myResp.all.edu)
myRespXY.all.edu<-coord.all.edu.names_thinned_dataset_full[[1]]

# function to define the intersect of rasters
intersect_mask <- function(x){
  values_x <- getValues(x)
  inter_x <- values_x %*% rep(1,nlayers(x))
  mask <- setValues(subset(x,1),values = (inter_x>0))
  return(mask)
}
# keep only all cells that are defined for all layers

myExpl.all.edu<-stack(reg.predictors.all.edu)

myExpl.all.edu <- stack(mask(myExpl.all.edu, intersect_mask(myExpl.all.edu)))


##Number of pseudo-absence: Lobo & Tognelli sugerem 100 pontos de pseudo-
#ausencia para cada ponto de presenca.



myBiomodData.all.edu<- BIOMOD_FormatingData(resp.var = myResp.all.edu,
                                            expl.var = myExpl.all.edu,
                                            resp.xy = myRespXY.all.edu,
                                            resp.name = myRespName.all.edu,
                                            PA.nb.rep = 1,
                                            PA.nb.absences = 10000,
                                            PA.strategy = 'disk',
                                            PA.dist.max=500000) # com base em http://www.sciencedirect.com/science/article/pii/S0304380008005486
myBiomodData.all.edu
save.image()
###MODELING

# 1. Defining Models Options using default options.
myBiomodOption.all.edu <- BIOMOD_ModelingOptions()

# 3. Computing the models

myBiomodModelOut.all.edu<- BIOMOD_Modeling(myBiomodData.all.edu,
                                           models = c("MAXENT.Phillips"),
                                           models.options = myBiomodOption.all.edu,
                                           NbRunEval=5, #conforme http://onlinelibrary.wiley.com/doi/10.1111/j.1523-1739.2006.00354.x/full
                                           DataSplit=70,
                                           Prevalence=0.5,
                                           VarImport=1,
                                           models.eval.meth = c('TSS'),
                                           all.eduveObj = TRUE,
                                           rescal.all.models = TRUE,
                                           do.full.models = FALSE,
                                           modeling.id = paste(myRespName.all.edu,"FirstModeling",sep=""))
(myBiomodModelOut.all.edu)
save.image()
# get all models evaluation

#Excelente TSS>0.75
#Bom 0.40>=TSS<0.75
#Ruim TSS<0.40


myBiomodModelEval.all.edu <- get_evaluations(myBiomodModelOut.all.edu)

myBiomodModelEval.all.edu




get_variables_importance(myBiomodModelOut.all.edu)
get_evaluations_matrix_all.edu<-cbind(get_variables_importance(myBiomodModelOut.all.edu)[1:9], #RUN 1
                                      get_variables_importance(myBiomodModelOut.all.edu)[10:18], #RUN 2
                                      get_variables_importance(myBiomodModelOut.all.edu)[19:27], #RUN 3
                                      get_variables_importance(myBiomodModelOut.all.edu)[28:36], #RUN 4
                                      get_variables_importance(myBiomodModelOut.all.edu)[37:45]) #RUN 5
colnames(get_evaluations_matrix_all.edu)<-c("RUN1", "RUN2", "RUN3", "RUN4", "RUN5")
names(reg.predictors.all.edu)
rownames(get_evaluations_matrix_all.edu)<-c("bio2", "bio3", "bio8", "bio9", "bio13", "bio14", "bio15", "bio18", "bio19")
get_evaluations_matrix_all.edu_mean<-apply(get_evaluations_matrix_all.edu, 1, mean)
which.max(get_evaluations_matrix_all.edu_mean)



mymaxent.all.edu <- BIOMOD_LoadModels(myBiomodModelOut.all.edu, models='MAXENT.Phillips')
myRespPlot2D.all.edu <- response.plot2(models  = mymaxent.all.edu,
                                       Data = get_formal_data(myBiomodModelOut.all.edu,'expl.var'), 
                                       show.variables= get_formal_data(myBiomodModelOut.all.edu,'expl.var.names'),
                                       do.bivariate = FALSE,
                                       fixed.var.metric = 'mean',
                                       col = c("blue", "red"),
                                       legend = F,
                                       data_species = get_formal_data(myBiomodModelOut.all.edu,'resp.var'),
                                       plot=F)

save.image()

bio_3_response_all.edu<-gather(myRespPlot2D.all.edu$bio3,model,Probability,
                               all.edu_PA1_RUN1_MAXENT.Phillips:all.edu_PA1_RUN5_MAXENT.Phillips)



p_all.edu_bio3<-ggplot(data =bio_3_response_all.edu) + 
  geom_smooth(mapping = aes(x = bio3, y = Probability))+
  labs(x = bquote('Bio 3'))
p_all.edu_bio3



###Ensemble modeling
myBiomodEM.all.edu<-BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut.all.edu,
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

myBiomodProjection.all.edu <- BIOMOD_Projection(modeling.output = myBiomodModelOut.all.edu,
                                                new.env = stack(ma.presente.poly.all.edu),
                                                proj.name = 'current_MA',
                                                selected.models = "all",
                                                binary.meth = 'TSS',
                                                compress = FALSE,
                                                build.clamping.mask = FALSE)

#Ensemble Forcasting

myBiomodEF.all.edu.presente <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM.all.edu,
  projection.output = myBiomodProjection.all.edu,binary.meth = "TSS")

save.image()

plot(myBiomodEF.all.edu.presente)

currentPred.ensemble.all.edu <- stack("all.edu/proj_current_MA/proj_current_MA_all.edu_ensemble.grd")
plot(currentPred.ensemble.all.edu[[1]])

plot(currentPred.ensemble.all.edu[[1]],col = gray.colors(12, start = 0.9, end = 0.0))
plot(wrld_simpl, xlim=c(-60,-40), ylim=c(-34,-15), axes=TRUE, add=T)
polygon(poligono_mata_atlantica)


plot(currentPred.ensemble.all.edu[[2]],col = gray.colors(12, start = 0.9, end = 0.0))
points(coord.all.edu.names_thinned_dataset_full[[1]], cex=.6,col="black")
polygon(poligono_mata_atlantica)



########################## 2070 ##########################
################################################################


#install.packages("rgdal")
#library(rgdal)
#Worst cenario In RCP 4.5, emissions continue to rise throughout the 21st century
#GCM HadGEM2-SE
BioC.fut.HE.85<-getData('CMIP5', var='bio', res=10, rcp=85, model='HE', year=70,download=T) 


#plot(BioC.fut.HE.85, 1)
variaveis.fut.HE.2070.85 <- stack(BioC.fut.HE.85)

###MATA ATLANTICA projection 2070 RCP 85

ma.2070 = crop(variaveis.fut.HE.2070.85, poligono_mata_atlantica)
plot(ma.2070, 1,  axes=TRUE)
polygon(poligono_mata_atlantica, add=T)
ma.2070.poly <- mask(ma.2070, data.shape)
ma.2070.poly <- mask(ma.2070, )

plot(ma.2070.poly, 1 )
plot(ma.presente.poly,1)



ma.2070.poly.all.edu <- stack(ma.2070.new[[2]],
                                  ma.2070.new[[3]],
                                  ma.2070.new[[8]],
                                  ma.2070.new[[9]],
                                  ma.2070.new[[13]],
                                  ma.2070.new[[14]],
                                  ma.2070.new[[15]],
                                  ma.2070.new[[18]],
                                  ma.2070.new[[19]])



names(ma.2070.poly.all.edu)<-names(ma.presente.poly.all.edu)


myExplFuture.2070.85.ma<-stack(ma.2070.poly.all.edu) #modelo HE


projection.2070.85.ma.all.edu <-BIOMOD_Projection(modeling.output=myBiomodModelOut.all.edu,
                                                   new.env= myExplFuture.2070.85.ma,
                                                   proj.name='all.edu.2070.85',
                                                   selected.models =  "all",
                                                   Bin.trans=TRUE,
                                                   binary.meth ='TSS',
                                                   compress = 'gzip',
                                                   build.clamping.mask	= FALSE,
                                                   SaveObj=TRUE,
                                                   keep.in.memory=F,
                                                   do.stack=T,
                                                   silent=F)


all.edu.2070.85<-get_predictions(projection.2070.85.ma.all.edu)
all.edu.2070.85
summary(all.edu.2070.85)





################## ENSAMBLE MODELLING 2070  ################
#################################################################



EnsambleForecast.2070.85.ma.all.edu <-BIOMOD_EnsembleForecasting (
  projection.output = projection.2070.85.ma.all.edu,
  EM.output=myBiomodEM.all.edu,
  selected.models = "all",
  binary.meth ='TSS') 



save.image()


EnsambleForecast.2070.85.ma.all.edu
plot(EnsambleForecast.2070.85.ma.all.edu)



ensemble.2070.HE.85.ma.all.edu <- stack("all.edu/proj_all.edu.2070.85/proj_all.edu.2070.85_all.edu_ensemble.grd")
ensemble.2070.HE.85.ma.all.edu



plot(ensemble.2070.HE.85.ma.all.edu[[1]],col = gray.colors(12, start = 0.9, end = 0))
points(coord.gym_klo.names_thinned_dataset_full[[1]],pch=20, col="black", cex=.7)

plot(ensemble.2070.HE.85.ma.all.edu[[2]],col = gray.colors(12, start = 0.9, end = 0))



futurePred.ensemble.2070.85.all.edu <- stack("all.edu/proj_all.edu.2070.85/proj_all.edu.2070.85_all.edu_ensemble.grd")



plot(currentPred.ensemble.all.edu[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.2070.85.all.edu[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS


currentPred.ensemble.bin.all.edu<- stack("all.edu/proj_current_MA/proj_current_MA_all.edu_ensemble_TSSbin.grd")
futurePred.ensemble.bin.2070.85.all.edu <- stack("all.edu/proj_all.edu.2070.85/proj_all.edu.2070.85_all.edu_ensemble_TSSbin.gri")


plot(currentPred.ensemble.bin.all.edu[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.bin.2070.85.all.edu[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS


# plot(futurePred.ensemble.bin.85.ma,col = gray.colors(12, start = 0.9, end = 0), ext=e_ma, ylim=c(-40,-3))#TSS
# plot(wrld_simpl,  axes="TRUE", add=T)


#  projection(futurePred.ensemble.bin.2070.85.all.edu[[1]]) <- projection(currentPred.ensemble.bin.all.edu[[1]])
# plot(futurePred.ensemble.bin.2070.85.all.edu[[1]])
# plot(currentPred.ensemble.bin.all.edu[[1]])
# 
#  
# distpobResamp.all.edu.2070.85 <- resample(futurePred.ensemble.bin.2070.85.all.edu[[1]], 
# currentPred.ensemble.bin.all.edu[[1]], resample='bilinear')
# plot(distpobResamp.all.edu.2070.85)
#  s.all.edu.2070.85<-stack(distpobResamp.all.edu.2070.85, currentPred.ensemble.bin.all.edu[[1]]) #creates stack, 's' 
# plot( s.all.edu.2070.85[[2]])



myBiomodRangeSize.drimys.2070.85.all.edu <- BIOMOD_RangeSize(
  CurrentPred=currentPred.ensemble.bin.all.edu[[1]],
  FutureProj=futurePred.ensemble.bin.2070.85.all.edu[[1]])
myBiomodRangeSize.drimys.2070.85.all.edu$Compt.By.Models





plot(myBiomodRangeSize.drimys.2070.85.all.edu$Diff.By.Pixel, col= gray.colors(4, start = 0, end = 1, gamma = 1, alpha = NULL))
plot(wrld_simpl,  axes="TRUE", add=T)

polygon(poligono_mata_atlantica)









