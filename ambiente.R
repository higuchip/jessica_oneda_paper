###DADOS AMBIENTAIS
library(raster)
library(usdm)

BioC<- getData('worldclim', var='bio', res=10) # Bioclima
1.85 * 1.42 * 10


#Alt<- getData('worldclim', var='alt', res=2.5) # Alt #Sem altitude por indicacao do revisor

variaveis <- stack(BioC)
plot(variaveis,1)




YbrevRange.MAS = extent(-62,-35,-34,-10) # define the extent conforme quali mariele
YbrevRange.MA = extent(poligono_mata_atlantica)
YbrevRange.SA = extent(-130,-35,-45,35) # define the extent conforme quali mariele



reg.variaveis.MA= crop(variaveis, YbrevRange.MA)
reg.variaveis.SA= crop(variaveis, YbrevRange.SA)


plot(reg.variaveis.SA, 1)
polygon(poligono_mata_atlantica)

summary(reg.variaveis.MA)




# YbrevRange=extent(-100,-34,-35,20) #define the extent
# reg.variaveis= crop(variaveis, YbrevRange)
# plot(reg.variaveis, 1)

reg.variaveis.MA
#remoção variaveis nao colineares
vari.n.colineares<-vifstep(reg.variaveis.SA)
vari.n.colineares

plot(reg.predictors.SA$bio2)

reg.predictors.SA<- stack(reg.variaveis.SA[[2]],
                           reg.variaveis.SA[[3]],
                           reg.variaveis.SA[[8]],
                           reg.variaveis.SA[[9]],
                           reg.variaveis.SA[[13]],
                           reg.variaveis.SA[[14]],
                           reg.variaveis.SA[[15]],
                           reg.variaveis.SA[[18]],
                          reg.variaveis.SA[[19]])


names(reg.predictors.SA)
