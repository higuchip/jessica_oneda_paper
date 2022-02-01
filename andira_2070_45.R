projection.2070.45.MAS.af <-BIOMOD_Projection(modeling.output=myBiomodModelOut.af,
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


af.2070.45.sa<-get_predictions(projection.2070.45.MAS.af)
af.2070.45.sa
summary(af.2070.45.sa)



af.2070.45.pred.MAS<-stack("Andira/proj_2070.45/proj_2070.45_Andira_TSSbin.grd")
af.2070.45.pred.MAS



################## ENSAMBLE MODELLING 2070  ################
#################################################################



EnsambleForecast.2070.45.MAS.af <-BIOMOD_EnsembleForecasting (
  projection.output = projection.2070.45.MAS.af,
  EM.output=myBiomodEM.af,
  selected.models = "all",
  binary.meth ='TSS') 



save.image(file="jessica.RData")


EnsambleForecast.2070.45.MAS.af
plot(EnsambleForecast.2070.45.MAS.af)



ensemble.2070.HE.45.MAS.af <- stack("Andira/proj_2070.45/proj_2070.45_Andira_ensemble.grd")
ensemble.2070.HE.45.MAS.af



plot(ensemble.2070.HE.45.MAS.af[[1]],col = gray.colors(12, start = 0.9, end = 0))
points(coord.af.limpo,pch=20, col="black", cex=.7)

plot(ensemble.2070.HE.45.MAS.af[[2]],col = gray.colors(12, start = 0.9, end = 0))
points(coord.af.limpo,pch=20, col="black", cex=.7)



currentPred.ensemble.af <- stack("Andira/proj_current_MAS/proj_current_MAS_Andira_ensemble.grd")
futurePred.ensemble.2070.45.af <- stack("Andira/proj_2070.45/proj_2070.45_Andira_ensemble.grd")

plot(currentPred.ensemble.af[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.2070.45.af[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS

currentPred.ensemble.bin.af<- stack("Andira/proj_current_MAS/proj_current_MAS_Andira_ensemble_TSSbin.grd")
futurePred.ensemble.bin.2070.45.af <- stack("Andira/proj_2070.45/proj_2070.45_Andira_ensemble_TSSbin.grd")

plot(currentPred.ensemble.bin.af[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.bin.2070.45.af[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS


# plot(futurePred.ensemble.bin.45.ma,col = gray.colors(12, start = 0.9, end = 0), ext=e_ma, ylim=c(-40,-3))#TSS
# plot(wrld_simpl,  axes="TRUE", add=T)




myBiomodRangeSize.drimys.2070.45.af <- BIOMOD_RangeSize(
  CurrentPred=currentPred.ensemble.bin.af[[1]],
  FutureProj=futurePred.ensemble.bin.2070.45.af[[1]])
myBiomodRangeSize.drimys.2070.45.af





plot(myBiomodRangeSize.drimys.2070.45.af$Diff.By.Pixel, col= gray.colors(4, start = 0, end = 1, gamma = 1, alpha = NULL))
plot(wrld_simpl,  axes="TRUE", add=T)




