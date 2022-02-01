library(biomod2)
library(splancs)
library(SDMTools)
library(dplyr)
library(tidyr)
library(gridExtra) 
library(raster)
library(spThin)

data(world)

plot(wrld_simpl, xlim=c(-130,-35), ylim=c(-45,35), axes=TRUE, col="gray")
points(coord.gym.klo, pch=19, cex=.51, col="black")


##Poligono MA

require(rgdal)
data.shape<-readOGR(dsn="/media/dendrologia/423EF1863EF172F1/User/Documents/orientados/Guilherme/IDS_09_Desm_Bioma_Mata_Atlantica_1_2", 
                    layer="IDS_09_Desm_Bioma_Mata_Atlantica_1_2")

poligono_mata_atlantica<-data.shape@polygons[[1]]@Polygons[[1]]@coords

polygon(poligono_mata_atlantica)


###Poligono especie


poligono_gym.klo<-getpoly()

colnames(poligono_gym.klo) = c("Long", "Lat")
polygon(poligono_gym.klo) 



########REMOCAO DUVIDAS

out.gym.klo= pnt.in.poly(coord.gym.klo,poligono_gym.klo)
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
plot(data.shape, add=T, col="lightgreen")
points(out.gym.klo[which(out.gym.klo$pip==1),1:2],pch=20, cex=.5, col="black")
coord.gym.klo.limpo<-out.gym.klo[which(out.gym.klo$pip==1),1:2]
plot(wrld_simpl, xlim=c(-120,-51), ylim=c(-45,35), axes=TRUE, col="white")
points(coord.gym.klo.limpo,pch=20, cex=.5, col="black")


###FILTRAGEM ESPACIAl

head(coord.gym.klo.limpo)
dim(coord.gym.klo.limpo)
gym_klo_names<-c(rep("gym_klo",  dim(coord.gym.klo.limpo)[1]))
coord.gym_klo.limpo.names<-cbind(gym_klo_names, coord.gym.klo.limpo)

head(coord.gym_klo.limpo.names)
is.data.frame(coord.gym_klo.limpo.names)
colnames(coord.gym_klo.limpo.names)<-c("SPEC", "LONG", "LAT")

coord.gym_klo.names_thinned_dataset_full <-
  thin( loc.data =coord.gym_klo.limpo.names, 
        lat.col = "LAT", 
        long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 26.27, reps =1, #10 min resolution
        locs.thinned.list.return = TRUE, 
        write.files = T, 
        max.files = 5, 
        out.dir = "coord.gym_klo.limpo.names_full/", out.base = "coord.gym_klo.limpo.names_thinned", 
        write.log.file = TRUE,
        log.file = "coord.gym_klo.limpo.names_thinned_full_log_file.txt" )




dim(coord.gym_klo.names_thinned_dataset_full[[1]])



###ENVIRONMENTAL BACKGROUND

reg.variaveis.gym.klo= crop(variaveis, poligono_gym.klo)
plot(reg.variaveis.gym.klo, 1)
points(coord.gym.klo.limpo, pch=20, cex=.5)

#remoção variaveis nao colineares
vari.n.colineares.gym.klo<-vifstep(reg.variaveis.gym.klo)
vari.n.colineares.gym.klo


reg.predictors.gym.klo<- stack(reg.variaveis.gym.klo[[2]],
                          reg.variaveis.gym.klo[[3]],
                          reg.variaveis.gym.klo[[8]],
                          reg.variaveis.gym.klo[[9]],
                          reg.variaveis.gym.klo[[13]],
                          reg.variaveis.gym.klo[[14]],
                          reg.variaveis.gym.klo[[15]],
                          reg.variaveis.gym.klo[[18]],
                          reg.variaveis.gym.klo[[19]])



names(reg.predictors.gym.klo)


###MATA ATLANTICA projection presente

ma.presente = crop(variaveis, poligono_mata_atlantica)
plot(ma.presente, 1)
polygon(poligono_mata_atlantica, add=T)
ma.presente.poly <- mask(MA.presente, data.shape)
plot(ma.presente.poly, 1 )

ma.presente.poly.gym.klo <- stack(ma.presente.poly[[2]],
                         ma.presente.poly[[3]],
                         ma.presente.poly[[8]],
                         ma.presente.poly[[9]],
                         ma.presente.poly[[13]],
                         ma.presente.poly[[14]],
                         ma.presente.poly[[15]],
                         ma.presente.poly[[18]],
                         ma.presente.poly[[19]])



names(ma.presente.poly.gym.klo)




ma.presente.poly.gym.klo.new <- stack(ma.presente.new[[2]],
                                  ma.presente.new[[3]],
                                  ma.presente.new[[8]],
                                  ma.presente.new[[9]],
                                  ma.presente.new[[13]],
                                  ma.presente.new[[14]],
                                  ma.presente.new[[15]],
                                  ma.presente.new[[18]],
                                  ma.presente.new[[19]])



names(ma.presente.poly.gym.klo.new)


###DADOS AMBIENTAIS

# 
# coord.gk.climate<-extract(reg.predictors,coord.gk.limpo)
# coord.gk.climate<-cbind(coord.gk.climate, coord.gk.limpo)
# #summary(coord.ilex.climate)
# coord.gk.climate[is.na(coord.gk.climate)]
# coord.gk.climate<-na.omit(coord.gk.climate)
# colnames(coord.gk.climate)
# 
# clim.pca.gk<-rda(coord.gk.climate[,c(1:9)], scale=T) # The scale argument
# 
# clim.pca.gk
# scores(clim.pca.gk, choices = 1:3, display = "species", scaling = 0)
# 
# eixosPCA.clim.gk<-scores(clim.pca.gk, choices = 1:2, display = "sites")
# clim.PCA1.gk<-eixosPCA.clim.gk[,1]
# is.vector(clim.PCA1.gk)
# clim.PCA2.gk<-eixosPCA.clim.gk[,2]
# 
# screeplot(clim.pca.gk,  bstick = TRUE, type = "lines")
# biplot(clim.pca.gk, scaling = 3, display="species", type="points", col="black", 
#        font=6, cex.lab=1.5, cex.axis=2, ylim=c(-3,3), xlim=c(-3,3), font.lab=6)
# text(clim.pca.gk, display = "species", 
#      scaling = 3, cex = 1.3, col = "darkgray",
#      font=6)
# points(clim.pca.gk, "sites", cex=1.2, col="black", bg="black")
# 
# 
# #Filtragem ambiental, para remover bias de amostragem, ugym.klondo metodologia de varela
# colnames(coord.gk.limpo)
# gk.env.filtered<-envgym.klomple (coord.gk.climate[,10:11], filters=list (clim.PCA1.gk, clim.PCA2.gk), 
#                             res=list (.1, .1), do.plot=TRUE)
# plot(wrld_simpl, xlim=c(-55,-51), ylim=c(-35,30), axes=TRUE, col="white")
# points(gk.env.filtered,pch=20, col="black")
# 
##PREPARANDO DADOS PARA BIOMOD2

myRespName.gym.klo <-  'gym.klo'
myResp.gym.klo<-rep(1, nrow(coord.gym_klo.names_thinned_dataset_full[[1]]))
myResp.gym.klo
length(myResp.gym.klo)
myRespXY.gym.klo<-coord.gym_klo.names_thinned_dataset_full[[1]]
myExpl.gym.klo<-stack(reg.predictors.gym.klo)


##Number of pseudo-absence: Lobo & Tognelli sugerem 100 pontos de pseudo-
#ausencia para cada ponto de presenca.



myBiomodData.gym.klo<- BIOMOD_FormatingData(resp.var = myResp.gym.klo,
                                         expl.var = myExpl.gym.klo,
                                         resp.xy = myRespXY.gym.klo,
                                         resp.name = myRespName.gym.klo,
                                         PA.nb.rep = 1,
                                         PA.nb.absences = 10000,
                                         PA.strategy = 'disk',
                                         PA.dist.max=500000) # com base em http://www.sciencedirect.com/science/article/pii/S0304380008005486
myBiomodData.gym.klo
save.image()
###MODELING

# 1. Defining Models Options using default options.
myBiomodOption.gym.klo <- BIOMOD_ModelingOptions()

# 3. Computing the models

myBiomodModelOut.gym.klo<- BIOMOD_Modeling(myBiomodData.gym.klo,
                                         models = c("MAXENT.Phillips"),
                                         models.options = myBiomodOption.gym.klo,
                                         NbRunEval=5, #conforme http://onlinelibrary.wiley.com/doi/10.1111/j.1523-1739.2006.00354.x/full
                                         DataSplit=70,
                                         Prevalence=0.5,
                                         VarImport=1,
                                         models.eval.meth = c('TSS'),
                                         gym.kloveObj = TRUE,
                                         rescal.all.models = TRUE,
                                         do.full.models = FALSE,
                                         modeling.id = paste(myRespName.gym.klo,"FirstModeling",sep=""))
(myBiomodModelOut.gym.klo)
save.image()
# get all models evaluation

#Excelente TSS>0.75
#Bom 0.40>=TSS<0.75
#Ruim TSS<0.40


myBiomodModelEval.gym.klo <- get_evaluations(myBiomodModelOut.gym.klo)

myBiomodModelEval.gym.klo



 
get_variables_importance(myBiomodModelOut.gym.klo)
get_evaluations_matrix_gym.klo<-cbind(get_variables_importance(myBiomodModelOut.gym.klo)[1:9], #RUN 1
get_variables_importance(myBiomodModelOut.gym.klo)[10:18], #RUN 2
get_variables_importance(myBiomodModelOut.gym.klo)[19:27], #RUN 3
get_variables_importance(myBiomodModelOut.gym.klo)[28:36], #RUN 4
get_variables_importance(myBiomodModelOut.gym.klo)[37:45]) #RUN 5
colnames(get_evaluations_matrix_gym.klo)<-c("RUN1", "RUN2", "RUN3", "RUN4", "RUN5")
names(reg.predictors.gym.klo)
rownames(get_evaluations_matrix_gym.klo)<-c("bio2", "bio3", "bio8", "bio9", "bio13", "bio14", "bio15", "bio18", "bio19")
get_evaluations_matrix_gym.klo_mean<-apply(get_evaluations_matrix_gym.klo, 1, mean)
which.max(get_evaluations_matrix_gym.klo_mean)
  


mymaxent.gym.klo <- BIOMOD_LoadModels(myBiomodModelOut.gym.klo, models='MAXENT.Phillips')
myRespPlot2D.gym.klo <- response.plot2(models  = mymaxent.gym.klo,
                                        Data = get_formal_data(myBiomodModelOut.gym.klo,'expl.var'), 
                                        show.variables= get_formal_data(myBiomodModelOut.gym.klo,'expl.var.names'),
                                        do.bivariate = FALSE,
                                        fixed.var.metric = 'mean',
                                        col = c("blue", "red"),
                                        legend = F,
                                        data_species = get_formal_data(myBiomodModelOut.gym.klo,'resp.var'),
                                        plot=F)

save.image()

bio_8_response_gym.klo<-gather(myRespPlot2D.gym.klo$bio8,model,Probability,
                               gym.klo_PA1_RUN1_MAXENT.Phillips:gym.klo_PA1_RUN5_MAXENT.Phillips)



p_gym.klo_bio8<-ggplot(data =bio_8_response_gym.klo) + 
  geom_smooth(mapping = aes(x = bio8/10, y = Probability))+
  labs(x = bquote('Bio 8'))
p_gym.klo_bio8



###Ensemble modeling
myBiomodEM.gym.klo<-BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut.gym.klo,
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

myBiomodProjection.gym.klo <- BIOMOD_Projection(modeling.output = myBiomodModelOut.gym.klo,
                                                 new.env = stack(ma.presente.poly.gym.klo.new),
                                                 proj.name = 'current_MA',
                                                 selected.models = "all",
                                                 binary.meth = 'TSS',
                                                 compress = FALSE,
                                                 build.clamping.mask = FALSE)

#Ensemble Forcasting

myBiomodEF.gym.klo.presente <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM.gym.klo,
  projection.output = myBiomodProjection.gym.klo,binary.meth = "TSS")

save.image()

plot(myBiomodEF.gym.klo.presente)

currentPred.ensemble.gym.klo <- stack("gym.klo/proj_current_MA/proj_current_MA_gym.klo_ensemble.grd")
plot(currentPred.ensemble.gym.klo[[1]])

plot(currentPred.ensemble.gym.klo[[1]],col = gray.colors(12, start = 0.9, end = 0.0))
plot(wrld_simpl, xlim=c(-60,-40), ylim=c(-34,-15), axes=TRUE, add=T)
polygon(poligono_mata_atlantica)


plot(currentPred.ensemble.gym.klo[[2]],col = gray.colors(12, start = 0.9, end = 0.0))
points(coord.gym_klo.names_thinned_dataset_full[[1]], cex=.6,col="black")
polygon(poligono_mata_atlantica)


########################### 2070 ##########################
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
plot(ma.2070.poly, 1 )
plot(ma.presente.poly,1)


ma.2070.poly.gym.klo <- stack(ma.2070.poly[[2]],
                                  ma.2070.poly[[3]],
                                  ma.2070.poly[[8]],
                                  ma.2070.poly[[9]],
                                  ma.2070.poly[[13]],
                                  ma.2070.poly[[14]],
                                  ma.2070.poly[[15]],
                                  ma.2070.poly[[18]],
                                  ma.2070.poly[[19]])

ma.2070.poly.gym.klo <- stack(ma.2070.new[[2]],
                              ma.2070.new[[3]],
                              ma.2070.new[[8]],
                              ma.2070.new[[9]],
                              ma.2070.new[[13]],
                              ma.2070.new[[14]],
                              ma.2070.new[[15]],
                              ma.2070.new[[18]],
                              ma.2070.new[[19]])


names(ma.2070.poly.gym.klo)<-names(ma.presente.poly.gym.klo)


myExplFuture.2070.85.ma<-stack(ma.2070.poly.gym.klo) #modelo HE


projection.2070.85.MAS.gym.klo <-BIOMOD_Projection(modeling.output=myBiomodModelOut.gym.klo,
                                              new.env= myExplFuture.2070.85.ma,
                                              proj.name='gym.klo.2070.85',
                                              selected.models =  "all",
                                              Bin.trans=TRUE,
                                              binary.meth ='TSS',
                                              compress = 'gzip',
                                              build.clamping.mask	= FALSE,
                                              SaveObj=TRUE,
                                              keep.in.memory=F,
                                              do.stack=T,
                                              silent=F)

projection.2070.85.ma.gym.klo <- projection.2070.85.MAS.gym.klo 
plot(projection.2070.85.ma.gym.klo)
save.image()


gym.klo.2070.85<-get_predictions(projection.2070.85.ma.gym.klo)
gym.klo.2070.85
summary(gym.klo.2070.85)



gk.2070.85.pred.ma<-stack("gym.klo/proj_gym.klo.2070.85/proj_gym.klo.2070.85_gym.klo_TSSbin.grd")
gk.2070.85.pred.ma



################## ENSAMBLE MODELLING 2070  ################
#################################################################



EnsambleForecast.2070.85.ma.gym.klo <-BIOMOD_EnsembleForecasting (
  projection.output = projection.2070.85.ma.gym.klo,
  EM.output=myBiomodEM.gym.klo,
  selected.models = "all",
  binary.meth ='TSS') 



save.image()


EnsambleForecast.2070.85.ma.gym.klo
plot(EnsambleForecast.2070.85.ma.gym.klo)



ensemble.2070.HE.85.ma.gym.klo <- stack("gym.klo/proj_gym.klo.2070.85/proj_gym.klo.2070.85_gym.klo_ensemble.grd")
ensemble.2070.HE.85.ma.gym.klo



plot(ensemble.2070.HE.85.ma.gym.klo[[1]],col = gray.colors(12, start = 0.9, end = 0))
points(coord.gym_klo.names_thinned_dataset_full[[1]],pch=20, col="black", cex=.7)

plot(ensemble.2070.HE.85.ma.gym.klo[[2]],col = gray.colors(12, start = 0.9, end = 0))



futurePred.ensemble.2070.85.gym.klo <- stack("gym.klo/proj_gym.klo.2070.85/proj_gym.klo.2070.85_gym.klo_ensemble.grd")



plot(currentPred.ensemble.gym.klo[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.2070.85.gym.klo[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS


currentPred.ensemble.bin.gym.klo<- stack("gym.klo/proj_current_MA/proj_current_MA_gym.klo_ensemble_TSSbin.grd")
futurePred.ensemble.bin.2070.85.gym.klo <- stack("gym.klo/proj_gym.klo.2070.85/proj_gym.klo.2070.85_gym.klo_ensemble_TSSbin.gri")


plot(currentPred.ensemble.bin.gym.klo[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.bin.2070.85.gym.klo[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS


# plot(futurePred.ensemble.bin.85.ma,col = gray.colors(12, start = 0.9, end = 0), ext=e_ma, ylim=c(-40,-3))#TSS
# plot(wrld_simpl,  axes="TRUE", add=T)


#  projection(futurePred.ensemble.bin.2070.85.gym.klo[[1]]) <- projection(currentPred.ensemble.bin.gym.klo[[1]])
# plot(futurePred.ensemble.bin.2070.85.gym.klo[[1]])
# plot(currentPred.ensemble.bin.gym.klo[[1]])
# 
#  
# distpobResamp.gym.klo.2070.85 <- resample(futurePred.ensemble.bin.2070.85.gym.klo[[1]], 
# currentPred.ensemble.bin.gym.klo[[1]], resample='bilinear')
# plot(distpobResamp.gym.klo.2070.85)
#  s.gym.klo.2070.85<-stack(distpobResamp.gym.klo.2070.85, currentPred.ensemble.bin.gym.klo[[1]]) #creates stack, 's' 
# plot( s.gym.klo.2070.85[[2]])



myBiomodRangeSize.drimys.2070.85.gym.klo <- BIOMOD_RangeSize(
  CurrentPred=currentPred.ensemble.bin.gym.klo[[1]],
  FutureProj=futurePred.ensemble.bin.2070.85.gym.klo[[1]])
myBiomodRangeSize.drimys.2070.85.gym.klo$Compt.By.Models





plot(myBiomodRangeSize.drimys.2070.85.gym.klo$Diff.By.Pixel, col= gray.colors(4, start = 0, end = 1, gamma = 1, alpha = NULL))
plot(wrld_simpl,  axes="TRUE", add=T)

polygon(poligono_mata_atlantica)





citations.gym.klo<-BIEN_metadata_citation(dataframe=gym.klo.bien)#If you are referencing occurrence data
citations.alc.tri


