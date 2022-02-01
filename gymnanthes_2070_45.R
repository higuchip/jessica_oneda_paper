########################### 2070 ##########################
################################################################
#install.packages("rgdal")
#library(rgdal)
#Worst cenario In RCP 4.5, emissions continue to rise throughout the 21st century
#GCM HadGEM2-SE
BioC.fut.HE.45<-getData('CMIP5', var='bio', res=2.5, rcp=45, model='HE', year=70,download=T) 
#plot(BioC.fut.HE.45, 1)
variaveis.fut.HE.2070.45 <- stack(BioC.fut.HE.45)
reg.variaveis.2070.45.MAS<-crop(variaveis.fut.HE.2070.45, YbrevRange.MAS)
plot(reg.variaveis.2070.45.MAS, 1,  axes=TRUE)



#points(coord.limpo, pch=19, cex=.5, col="black")


reg.predictors.MAS.2070.45<- stack(reg.variaveis.2070.45.MAS[[2]],
                                   reg.variaveis.2070.45.MAS[[3]],
                                   reg.variaveis.2070.45.MAS[[8]],
                                   reg.variaveis.2070.45.MAS[[9]],
                                   reg.variaveis.2070.45.MAS[[13]], 
                                   reg.variaveis.2070.45.MAS[[14]],
                                   reg.variaveis.2070.45.MAS[[15]],
                                   reg.variaveis.2070.45.MAS[[18]],
                                   reg.variaveis.2070.45.MAS[[19]])



names(reg.predictors.MAS.2070.45)<-names(reg.predictors)


plot(reg.predictors.MAS.2070.45, 3,  axes=TRUE)



myExplFuture.2070.45.MAS<-stack(reg.predictors.MAS.2070.45) #modelo HE



# projection.2070.45.MAS.gk <-BIOMOD_Projection(modeling.output=myBiomodModelOut.gk,
#                                           new.env= myExplFuture.2070.45.MAS,
#                                           proj.name='2070.45',
#                                           selected.models =  c('Gymnanthes_PA1_RUN1_MAXENT.Phillips',
#                                                                'Gymnanthes_PA1_RUN2_MAXENT.Phillips',
#                                                                "Gymnanthes_PA2_RUN1_MAXENT.Phillips",
#                                                                "Gymnanthes_PA2_RUN2_MAXENT.Phillips",
#                                                                "Gymnanthes_PA3_RUN1_MAXENT.Phillips",
#                                                                "Gymnanthes_PA3_RUN2_MAXENT.Phillips"),
#                                           Bin.trans=TRUE,
#                                           binary.meth ='TSS',
#                                           compress = 'gzip',
#                                           build.clamping.mask	= FALSE,
#                                           SaveObj=TRUE,
#                                           keep.in.memory=F,
#                                           do.stack=T,
#                                           silent=F)

projection.2070.45.MAS.gk <-BIOMOD_Projection(modeling.output=myBiomodModelOut.gk,
                                              new.env= myExplFuture.2070.45.MAS,
                                              proj.name='2070.45',
                                              selected.models =  "all",
                                              Bin.trans=TRUE,
                                              binary.meth ='TSS',
                                              compress = 'gzip',
                                              build.clamping.mask	= FALSE,
                                              SaveObj=TRUE,
                                              keep.in.memory=F,
                                              do.stack=T,
                                              silent=F)


save.image(file="jessica.RData")


gk.2070.45.sa<-get_predictions(projection.2070.45.MAS.gk)
gk.2070.45.sa
summary(gk.2070.45.sa)



gk.2070.45.pred.MAS<-stack("Gymnanthes/proj_2070.45/proj_2070.45_Gymnanthes_TSSbin.grd")
gk.2070.45.pred.MAS



################## ENSAMBLE MODELLING 2070  ################
#################################################################



EnsambleForecast.2070.45.MAS.gk <-BIOMOD_EnsembleForecasting (
  projection.output = projection.2070.45.MAS.gk,
  EM.output=myBiomodEM.gk,
  selected.models = "all",
  binary.meth ='TSS') 



save.image(file="jessica.RData")


EnsambleForecast.2070.45.MAS.gk
plot(EnsambleForecast.2070.45.MAS.gk)



ensemble.2070.HE.45.MAS.gk <- stack("Gymnanthes/proj_2070.45/proj_2070.45_Gymnanthes_ensemble.grd")
ensemble.2070.HE.45.MAS.gk



plot(ensemble.2070.HE.45.MAS.gk[[1]],col = gray.colors(12, start = 0.9, end = 0))
points(coord.gk.limpo,pch=20, col="black", cex=.7)

plot(ensemble.2070.HE.45.MAS.gk[[2]],col = gray.colors(12, start = 0.9, end = 0))
points(coord.gk.limpo,pch=20, col="black", cex=.7)



currentPred.ensemble.gk <- stack("Gymnanthes/proj_current_MAS/proj_current_MAS_Gymnanthes_ensemble.grd")
futurePred.ensemble.2070.45.gk <- stack("Gymnanthes/proj_2070.45/proj_2070.45_Gymnanthes_ensemble.grd")

plot(currentPred.ensemble.gk[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.2070.45.gk[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS

currentPred.ensemble.bin.gk<- stack("Gymnanthes/proj_current_MAS/proj_current_MAS_Gymnanthes_ensemble_TSSbin.grd")
futurePred.ensemble.bin.2070.45.gk <- stack("Gymnanthes/proj_2070.45/proj_2070.45_Gymnanthes_ensemble_TSSbin.grd")

plot(currentPred.ensemble.bin.gk[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.bin.2070.45.gk[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS


# plot(futurePred.ensemble.bin.45.ma,col = gray.colors(12, start = 0.9, end = 0), ext=e_ma, ylim=c(-40,-3))#TSS
# plot(wrld_simpl,  axes="TRUE", add=T)




myBiomodRangeSize.drimys.2070.45.gk <- BIOMOD_RangeSize(
  CurrentPred=currentPred.ensemble.bin.gk[[1]],
  FutureProj=futurePred.ensemble.bin.2070.45.gk[[1]])
myBiomodRangeSize.drimys.2070.45.gk





plot(myBiomodRangeSize.drimys.2070.45.gk$Diff.By.Pixel, col= gray.colors(4, start = 0, end = 1, gamma = 1, alpha = NULL))
plot(wrld_simpl,  axes="TRUE", add=T)




