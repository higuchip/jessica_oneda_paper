library(rgbif)
library(biomod2)

#GBIF
key.gk <- name_backbone(name='Gymnanthes klotzschiana')$speciesKey
key.gk
dat.gk <- occ_search(taxonKey=key.gk, return='data', limit=1182)
gbifmap(dat.gk)
coord.gbif.gk<-cbind(dat.gk$decimalLongitude, dat.gk$decimalLatitude)
head(coord.gbif.gk)

plot(wrld_simpl, xlim=c(-55,-51), ylim=c(-35,10), axes=TRUE, col="gray")
points(coord.gbif.gk, pch=19, cex=.51, col="black")

head(coord.gbif.gk)

coord


# row.names(coord)
# coord.allo.bd<-coord[c(1:6,8,11:15,21:24,30,32,34,36:38,43,44,47,49,50),]
# colnames(coord.allo.bd)
# colnames(coord.gbif.ae)<-colnames(coord.allo.bd)
# coord.allo.total<-rbind(coord.gbif.ae, coord.allo.bd)
# plot(wrld_simpl, xlim=c(-55,-51), ylim=c(-35,30), axes=TRUE, col="gray")
# points(coord.allo.total, pch=19, cex=.51, col="black")
# save.image(file="jessica.RData")
 poligono<-getpoly()
 polygon(poligono)

colnames(poligono) = c("Long", "Lat")


#create check which points fall within the polygon
out.gk= pnt.in.poly(coord.gbif.gk,poligono)
plot(wrld_simpl, xlim=c(-55,-51), ylim=c(-35,30), axes=TRUE, col="gray")
points(out.gk[which(out.gk$pip==1),1:2],pch=20, col="black")


coord.gk.limpo<-out.gk[which(out.gk$pip==1),1:2]

###DADOS AMBIENTAIS



BioC<- getData('worldclim', var='bio', res=2.5) # Bioclima
Alt<- getData('worldclim', var='alt', res=2.5) # Alt #Sem altitude por indicacao do revisor

variaveis <- stack(BioC)



YbrevRange.MAS = extent(-62,-35,-34,-10) # define the extent conforme quali mariele
reg.variaveis.MAS= crop(variaveis, YbrevRange.MAS)
plot(reg.variaveis.MAS, 1, xlim=c(-62,-30))
summary(reg.variaveis.MAS)


YbrevRange=extent(-100,-34,-35,20) #define the extent
reg.variaveis= crop(variaveis, YbrevRange)
plot(reg.variaveis, 1)


#remoção variaveis nao colineares
vari.n.colineares.gk<-vifstep(reg.variaveis)
(vari.n.colineares.gk)

names(reg.variaveis)

reg.predictors.MAS<- stack(reg.variaveis.MAS[[2]],
                           reg.variaveis.MAS[[3]],
                           reg.variaveis.MAS[[8]],
                           reg.variaveis.MAS[[9]],
                           reg.variaveis.MAS[[13]], 
                           reg.variaveis.MAS[[14]],
                           reg.variaveis.MAS[[15]],
                           reg.variaveis.MAS[[18]],
                           reg.variaveis.MAS[[19]])

reg.predictors<- stack(reg.variaveis[[2]],
                           reg.variaveis[[3]],
                           reg.variaveis[[8]],
                           reg.variaveis[[9]],
                           reg.variaveis[[13]], 
                           reg.variaveis[[14]],
                           reg.variaveis[[15]],
                           reg.variaveis[[18]],
                           reg.variaveis[[19]])

names(reg.predictors)
dim(coord.gk.limpo)
coord.gk.climate<-extract(reg.predictors,coord.gk.limpo)
coord.gk.climate<-cbind(coord.gk.climate, coord.gk.limpo)
#summary(coord.ilex.climate)
coord.gk.climate[is.na(coord.gk.climate)]
coord.gk.climate<-na.omit(coord.gk.climate)
colnames(coord.gk.climate)

clim.pca.gk<-rda(coord.gk.climate[,c(1:9)], scale=T) # The scale argument

clim.pca.gk
scores(clim.pca.gk, choices = 1:3, display = "species", scaling = 0)

eixosPCA.clim.gk<-scores(clim.pca.gk, choices = 1:2, display = "sites")
clim.PCA1.gk<-eixosPCA.clim.gk[,1]
is.vector(clim.PCA1.gk)
clim.PCA2.gk<-eixosPCA.clim.gk[,2]

screeplot(clim.pca.gk,  bstick = TRUE, type = "lines")
biplot(clim.pca.gk, scaling = 3, display="species", type="points", col="black", 
       font=6, cex.lab=1.5, cex.axis=2, ylim=c(-3,3), xlim=c(-3,3), font.lab=6)
text(clim.pca.gk, display = "species", 
     scaling = 3, cex = 1.3, col = "darkgray",
     font=6)
points(clim.pca.gk, "sites", cex=1.2, col="black", bg="black")


#Filtragem ambiental, para remover bias de amostragem, usando metodologia de varela
colnames(coord.gk.limpo)
gk.env.filtered<-envSample (coord.gk.climate[,10:11], filters=list (clim.PCA1.gk, clim.PCA2.gk), 
                            res=list (.1, .1), do.plot=TRUE)
plot(wrld_simpl, xlim=c(-55,-51), ylim=c(-35,30), axes=TRUE, col="white")
points(gk.env.filtered,pch=20, col="black")

##PREPARANDO DADOS PARA BIOMOD2

myRespName.gk <-  'Gymnanthes'
myResp.gk<-rep(1, nrow(gk.env.filtered))
myResp.gk
length(myResp.gk)
myRespXY.gk<-gk.env.filtered
myExpl.gk<-stack(reg.predictors)
myExpl.MAS.gk<-stack(reg.predictors.MAS)


##Number of pseudo-absence: Lobo & Tognelli sugerem 100 pontos de pseudo-
#ausencia para cada ponto de presenca.
Sys.getenv('R_USER')
tempdir()

#Brazil
myBiomodData.gk<- BIOMOD_FormatingData(resp.var = myResp.gk,
                                         expl.var = myExpl.gk,
                                         resp.xy = myRespXY.gk,
                                         resp.name = myRespName.gk,
                                         PA.nb.rep = 3,
                                         PA.nb.absences = 9400,
                                         PA.strategy = 'disk',
                                         PA.dist.max=500000) # com base em http://www.sciencedirect.com/science/article/pii/S0304380008005486
myBiomodData.gk
save.image(file="jessica.RData")
###MODELING

# 1. Defining Models Options using default options.
myBiomodOption.gk <- BIOMOD_ModelingOptions()

# 3. Computing the models

myBiomodModelOut.gk <- BIOMOD_Modeling(myBiomodData.gk,
                                         models = c("MAXENT.Phillips"),
                                         models.options = myBiomodOption.allo,
                                         NbRunEval=5,
                                         DataSplit=70,
                                         Prevalence=0.5,
                                         VarImport=3,
                                         models.eval.meth = c('TSS'),
                                         SaveObj = TRUE,
                                         rescal.all.models = TRUE,
                                         do.full.models = FALSE,
                                         modeling.id = paste(myRespName.allo,"FirstModeling",sep=""))
(myBiomodModelOut.gk)
save.image(file="jessica.RData")
# get all models evaluation
myBiomodModelEval.gk <- get_evaluations(myBiomodModelOut.gk)
myBiomodModelEval.gk

# , , MAXENT.Phillips, RUN1, PA3
# 
# Testing.data Cutoff Sensitivity Specificity
# TSS        0.727    454      96.429      76.241

get_variables_importance(myBiomodModelOut.gk)
  
##Modelo melhor ajustado
mymaxent.gk <- BIOMOD_LoadModels(myBiomodModelOut.gk, models='MAXENT.Phillips')
myRespPlot2D.gk <- response.plot2(models  = mymaxent.gk,
                                        Data = get_formal_data(myBiomodModelOut.gk,'expl.var'), 
                                        show.variables= get_formal_data(myBiomodModelOut.gk,'expl.var.names'),
                                        do.bivariate = FALSE,
                                        fixed.var.metric = 'mean',
                                        col = c("blue", "red"),
                                        legend = TRUE,
                                        data_species = get_formal_data(myBiomodModelOut.gk,'resp.var'),
                                        plot=F)

jpeg(filename = "response_gk.jpg",
     width = 2000, height = 4000, units = "px",
     quality = 100,
     bg = "white", 
     res = 300,
     family="times")
par(mfcol=c(2,1),mar=c(3.5,4,1.5,0), oma=c(1.5,1,1.5,1), 
    mgp = c(2, 1, 0), family= "serif")


plot(myRespPlot2D.gk$bio8$bio8,
     myRespPlot2D.gk$bio8$Gymnanthes_PA3_RUN1_MAXENT.Phillips,
     type = "l", xlab = 'Temperatura no trimestre mais úmido',
     ylab = 'Estimativa de Adequabilidade Climática', lty=1, lwd=2, cex.lab=1.5, cex.axis=1.5)

plot(myRespPlot2D.gk$bio14$bio14,
     myRespPlot2D.gk$bio14$Gymnanthes_PA3_RUN1_MAXENT.Phillips,
     type = "l", xlab = 'Precipitação no mês mais seco',
     ylab = 'Estimativa de Adequabilidade Climática', lty=1, lwd=2, cex.lab=1.5, cex.axis=1.5)



dev.off()

###Ensemble modeling
myBiomodEM.gk<-BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut.gk,
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

myBiomodProjection.gk <- BIOMOD_Projection(modeling.output = myBiomodModelOut.gk,
                                                 new.env = stack(reg.predictors.MAS),
                                                 proj.name = 'current_MAS',
                                                 selected.models = "all",
                                                 binary.meth = 'TSS',
                                                 compress = FALSE,
                                                 build.clamping.mask = FALSE)

#Ensemble Forcasting

myBiomodEF.gk <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM.gk,
  projection.output = myBiomodProjection.gk,binary.meth = "TSS")

save.image(file="jessica.RData")

plot(myBiomodEF.gk)

currentPred.ensemble.gk <- stack("Gymnanthes/proj_current_MAS/proj_current_MAS_Gymnanthes_ensemble.grd")
plot(currentPred.ensemble.gk[[1]])
points(coord.gbif.gk,pch=20, col="black", cex=0.5)

plot(currentPred.ensemble.gk[[1]],col = gray.colors(12, start = 0.9, end = 0.0))
points(coord.gbif.gk,pch=19, cex=.6,col="black")
plot(wrld_simpl, xlim=c(-60,-40), ylim=c(-34,-15), axes=TRUE, add=T)

plot(currentPred.ensemble.gk[[2]],col = gray.colors(12, start = 0.9, end = 0.0))
points(coord.gbif.gk,pch=19, cex=.6,col="black")
plot(wrld_simpl, xlim=c(-60,-40), ylim=c(-34,-15), axes=TRUE, add=T)

