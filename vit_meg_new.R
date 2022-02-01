library(biomod2)

dim(coord.vit.meg.names_thinned_dataset_full[[1]])


myBiomodData.vit.meg<- BIOMOD_FormatingData(resp.var = myResp.vit.meg,
                                            expl.var = myExpl.vit.meg,
                                            resp.xy = myRespXY.vit.meg,
                                            resp.name = myRespName.vit.meg,
                                            PA.nb.rep = 1,
                                            PA.nb.absences = 10000,
                                            PA.strategy = 'disk',
                                            PA.dist.max=500000) # com base em http://www.sciencedirect.com/science/article/pii/S0304380008005486
myBiomodData.vit.meg
save.image()


# 1. Defining Models Options using default options.
myBiomodOption.vit.meg <- BIOMOD_ModelingOptions()

# 3. Computing the models

myBiomodModelOut.vit.meg<- BIOMOD_Modeling(myBiomodData.vit.meg,
                                           models = c("MAXENT.Phillips"),
                                           models.options = myBiomodOption.vit.meg,
                                           NbRunEval=5, #conforme http://onlinelibrary.wiley.com/doi/10.1111/j.1523-1739.2006.00354.x/full
                                           DataSplit=70,
                                           Prevalence=0.5,
                                           VarImport=1,
                                           models.eval.meth = c('TSS'),
                                           vit.megveObj = TRUE,
                                           rescal.all.models = TRUE,
                                           do.full.models = FALSE,
                                           modeling.id = paste(myRespName.vit.meg,"FirstModeling",sep=""))
(myBiomodModelOut.vit.meg)

# get all models evaluation

#Excelente TSS>0.75
#Bom 0.40>=TSS<0.75
#Ruim TSS<0.40


myBiomodModelEval.vit.meg <- get_evaluations(myBiomodModelOut.vit.meg)

myBiomodModelEval.vit.meg




get_variables_importance(myBiomodModelOut.vit.meg)
get_evaluations_matrix_vit.meg<-cbind(get_variables_importance(myBiomodModelOut.vit.meg)[1:10], #RUN 1
                                      get_variables_importance(myBiomodModelOut.vit.meg)[11:20], #RUN 2
                                      get_variables_importance(myBiomodModelOut.vit.meg)[21:30], #RUN 3
                                      get_variables_importance(myBiomodModelOut.vit.meg)[31:40], #RUN 4
                                      get_variables_importance(myBiomodModelOut.vit.meg)[41:50]) #RUN 5
colnames(get_evaluations_matrix_vit.meg)<-c("RUN1", "RUN2", "RUN3", "RUN4", "RUN5")
names(reg.predictors.vit.meg)
rownames(get_evaluations_matrix_vit.meg)<-c("bio2", "bio3", "bio8", "bio9", "bio10", "bio13", "bio14", "bio15", "bio18", "bio19")
get_evaluations_matrix_vit.meg_mean<-apply(get_evaluations_matrix_vit.meg, 1, mean)
which.max(get_evaluations_matrix_vit.meg_mean)

save.image()



mymaxent.vit.meg <- BIOMOD_LoadModels(myBiomodModelOut.vit.meg, models='MAXENT.Phillips')
myRespPlot2D.vit.meg <- response.plot2(models  = mymaxent.vit.meg,
                                       Data = get_formal_data(myBiomodModelOut.vit.meg,'expl.var'), 
                                       show.variables= get_formal_data(myBiomodModelOut.vit.meg,'expl.var.names'),
                                       do.bivariate = FALSE,
                                       fixed.var.metric = 'mean',
                                       col = c("blue", "red"),
                                       legend = F,
                                       data_species = get_formal_data(myBiomodModelOut.vit.meg,'resp.var'),
                                       plot=F)

save.image()

bio_10_response_vit.meg<-gather(myRespPlot2D.vit.meg$bio10,model,Probability,
                               vit.meg_PA1_RUN1_MAXENT.Phillips:vit.meg_PA1_RUN5_MAXENT.Phillips)



p_vit.meg_bio10<-ggplot(data =bio_10_response_vit.meg) + 
  geom_smooth(mapping = aes(x = bio10, y = Probability))+
  labs(x = bquote('Bio 10'))
p_vit.meg_bio10



###Ensemble modeling
myBiomodEM.vit.meg<-BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut.vit.meg,
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

plot(ma.presente.poly.vit.meg, 1)
ma.presente.poly.vit.meg
ma.presente.poly.vit.meg <- stack(ma.presente.new[[2]],
                                  ma.presente.new[[3]],
                                  ma.presente.new[[8]],
                                  ma.presente.new[[9]],
                                  ma.presente.new[[10]],
                                  ma.presente.new[[13]],
                                  ma.presente.new[[14]],
                                  ma.presente.new[[15]],
                                  ma.presente.new[[18]],
                                  ma.presente.new[[19]])


ma.2070.poly.vit.meg <- stack(ma.2070.new[[2]],
                              ma.2070.new[[3]],
                              ma.2070.new[[8]],
                              ma.2070.new[[9]],
                              ma.2070.new[[10]],
                              ma.2070.new[[13]],
                              ma.2070.new[[14]],
                              ma.2070.new[[15]],
                              ma.2070.new[[18]],
                              ma.2070.new[[19]])




myBiomodProjection.vit.meg <- BIOMOD_Projection(modeling.output = myBiomodModelOut.vit.meg,
                                                new.env = stack(ma.presente.poly.vit.meg),
                                                proj.name = 'current_MA',
                                                selected.models = "all",
                                                binary.meth = 'TSS',
                                                compress = FALSE,
                                                build.clamping.mask = FALSE)

#Ensemble Forcasting

myBiomodEF.vit.meg.presente <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM.vit.meg,
  projection.output = myBiomodProjection.vit.meg,binary.meth = "TSS")

save.image()

plot(myBiomodEF.vit.meg.presente)

currentPred.ensemble.vit.meg <- stack("vit.meg/proj_current_MA/proj_current_MA_vit.meg_ensemble.grd")
plot(currentPred.ensemble.vit.meg[[1]])

plot(currentPred.ensemble.vit.meg[[1]],col = gray.colors(12, start = 0.9, end = 0.0))
plot(wrld_simpl, xlim=c(-60,-40), ylim=c(-34,-15), axes=TRUE, add=T)
polygon(poligono_mata_atlantica)


plot(currentPred.ensemble.vit.meg[[2]],col = gray.colors(12, start = 0.9, end = 0.0))
points(coord.vit.meg.names_thinned_dataset_full[[1]], cex=.6,col="black")
polygon(poligono_mata_atlantica)



########################## 2070 ##########################
################################################################


#install.packages("rgdal")
#library(rgdal)
#Worst cenario In RCP 4.5, emissions continue to rise throughout the 21st century
#GCM HadGEM2-SE
# BioC.fut.HE.85<-getData('CMIP5', var='bio', res=10, rcp=85, model='HE', year=70,download=T) 
# 
# 
# #plot(BioC.fut.HE.85, 1)
# variaveis.fut.HE.2070.85 <- stack(BioC.fut.HE.85)
# 
# ###MATA ATLANTICA projection 2070 RCP 85
# 
# ma.2070 = crop(variaveis.fut.HE.2070.85, poligono_mata_atlantica)
# plot(ma.2070, 1,  axes=TRUE)
# polygon(poligono_mata_atlantica, add=T)
# ma.2070.poly <- mask(ma.2070, data.shape)
# ma.2070.poly <- mask(ma.2070, )
# 
# plot(ma.2070.poly, 1 )
# plot(ma.presente.poly,1)
# 
# 
# 
# ma.2070.poly.vit.meg <- stack(ma.2070.new[[2]],
#                               ma.2070.new[[3]],
#                               ma.2070.new[[8]],
#                               ma.2070.new[[9]],
#                               ma.2070.new[[13]],
#                               ma.2070.new[[14]],
#                               ma.2070.new[[15]],
#                               ma.2070.new[[18]],
#                               ma.2070.new[[19]])
# 


names(ma.2070.poly.vit.meg)<-names(ma.presente.poly.vit.meg)


myExplFuture.2070.85.ma<-stack(ma.2070.poly.vit.meg) #modelo HE


projection.2070.85.ma.vit.meg <-BIOMOD_Projection(modeling.output=myBiomodModelOut.vit.meg,
                                                  new.env= myExplFuture.2070.85.ma,
                                                  proj.name='vit.meg.2070.85',
                                                  selected.models =  "all",
                                                  Bin.trans=TRUE,
                                                  binary.meth ='TSS',
                                                  compress = 'gzip',
                                                  build.clamping.mask	= FALSE,
                                                  SaveObj=TRUE,
                                                  keep.in.memory=F,
                                                  do.stack=T,
                                                  silent=F)


vit.meg.2070.85<-get_predictions(projection.2070.85.ma.vit.meg)
vit.meg.2070.85
summary(vit.meg.2070.85)





################## ENSAMBLE MODELLING 2070  ################
#################################################################



EnsambleForecast.2070.85.ma.vit.meg <-BIOMOD_EnsembleForecasting (
  projection.output = projection.2070.85.ma.vit.meg,
  EM.output=myBiomodEM.vit.meg,
  selected.models = "all",
  binary.meth ='TSS') 



save.image()


EnsambleForecast.2070.85.ma.vit.meg
plot(EnsambleForecast.2070.85.ma.vit.meg)



ensemble.2070.HE.85.ma.vit.meg <- stack("vit.meg/proj_vit.meg.2070.85/proj_vit.meg.2070.85_vit.meg_ensemble.grd")
ensemble.2070.HE.85.ma.vit.meg



plot(ensemble.2070.HE.85.ma.vit.meg[[1]],col = gray.colors(12, start = 0.9, end = 0))
points(coord.gym_klo.names_thinned_dataset_full[[1]],pch=20, col="black", cex=.7)

plot(ensemble.2070.HE.85.ma.vit.meg[[2]],col = gray.colors(12, start = 0.9, end = 0))



futurePred.ensemble.2070.85.vit.meg <- stack("vit.meg/proj_vit.meg.2070.85/proj_vit.meg.2070.85_vit.meg_ensemble.grd")



plot(currentPred.ensemble.vit.meg[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.2070.85.vit.meg[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS


currentPred.ensemble.bin.vit.meg<- stack("vit.meg/proj_current_MA/proj_current_MA_vit.meg_ensemble_TSSbin.grd")
futurePred.ensemble.bin.2070.85.vit.meg <- stack("vit.meg/proj_vit.meg.2070.85/proj_vit.meg.2070.85_vit.meg_ensemble_TSSbin.gri")


plot(currentPred.ensemble.bin.vit.meg[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS
plot(futurePred.ensemble.bin.2070.85.vit.meg[[1]],col = gray.colors(12, start = 0.9, end = 0))#TSS


# plot(futurePred.ensemble.bin.85.ma,col = gray.colors(12, start = 0.9, end = 0), ext=e_ma, ylim=c(-40,-3))#TSS
# plot(wrld_simpl,  axes="TRUE", add=T)


#  projection(futurePred.ensemble.bin.2070.85.vit.meg[[1]]) <- projection(currentPred.ensemble.bin.vit.meg[[1]])
# plot(futurePred.ensemble.bin.2070.85.vit.meg[[1]])
# plot(currentPred.ensemble.bin.vit.meg[[1]])
# 
#  
# distpobResamp.vit.meg.2070.85 <- resample(futurePred.ensemble.bin.2070.85.vit.meg[[1]], 
# currentPred.ensemble.bin.vit.meg[[1]], resample='bilinear')
# plot(distpobResamp.vit.meg.2070.85)
#  s.vit.meg.2070.85<-stack(distpobResamp.vit.meg.2070.85, currentPred.ensemble.bin.vit.meg[[1]]) #creates stack, 's' 
# plot( s.vit.meg.2070.85[[2]])



myBiomodRangeSize.drimys.2070.85.vit.meg <- BIOMOD_RangeSize(
  CurrentPred=currentPred.ensemble.bin.vit.meg[[1]],
  FutureProj=futurePred.ensemble.bin.2070.85.vit.meg[[1]])
myBiomodRangeSize.drimys.2070.85.vit.meg$Compt.By.Models





plot(myBiomodRangeSize.drimys.2070.85.vit.meg$Diff.By.Pixel, col= gray.colors(4, start = 0, end = 1, gamma = 1, alpha = NULL))
plot(wrld_simpl,  axes="TRUE", add=T)

polygon(poligono_mata_atlantica)

save.image()










