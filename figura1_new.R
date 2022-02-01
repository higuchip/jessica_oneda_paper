library(grid)
library(gridExtra)
library(ggplot2)
library(ggthemes)
library(maps)
library(maptools)


gg1 <- ggplot() + 
  geom_polygon(data = wrld_simpl, aes(x=long, y = lat, group = group), fill = "white", color = "black") + 
  coord_fixed(1.3)  +
  coord_cartesian(xlim = c(-125, -30), ylim = c(-40, 25)) + theme_map(base_size = 5)

gg1

# Relacao de spp com > 100 pontos de ocorrencias filtrados e TSS > 0.4, disponivel no BIEN

#Group 1
# gym_klo

gym.klo.spatial.plot<-gg1 + geom_point(data = coord.gym_klo.limpo.names, 
                                       aes(x = LONG, y = LAT), 
                                       shape = 20, 
                                       color = "blue", 
                                       fill = "blue", 
                                       size = .5)



gym.klo.spatial.plot1 <- arrangeGrob(gym.klo.spatial.plot, 
                                     top = textGrob("A", x = unit(1, "npc"), 
                                       y   = unit(-4, "npc"), just=c("right"),,
                                      gp=gpar(col="black", 
                                      fontsize=12, 
                                      fontfamily="Times Roman")))



plot(gym.klo.spatial.plot1 )



# all.edu

all.edu.spatial.plot<-gg1 + geom_point(data = coord.all.edu.limpo.names, 
                                       aes(x = LONG, y = LAT), 
                                       shape = 20, 
                                       color = "blue", 
                                       fill = "blue", 
                                       size = .5)



all.edu.spatial.plot1 <- arrangeGrob(all.edu.spatial.plot, 
                                     top = textGrob("B", x = unit(1, "npc"), 
                                                     y   = unit(-4, "npc"), just=c("right"),,
                                                    gp=gpar(col="black", 
                                                            fontsize=12, 
                                                            fontfamily="Times Roman")))



plot(all.edu.spatial.plot1 )


# vit.meg

vit.meg.spatial.plot<-gg1 + geom_point(data = coord.vit.meg.limpo.names, 
                                       aes(x = LONG, y = LAT), 
                                       shape = 20, 
                                       color = "blue", 
                                       fill = "blue", 
                                       size = .5)



vit.meg.spatial.plot1 <- arrangeGrob(vit.meg.spatial.plot, 
                                     top = textGrob("C", x = unit(1, "npc"), 
                                                     y   = unit(-4, "npc"), just=c("right"),,
                                                    gp=gpar(col="black", 
                                                            fontsize=12, 
                                                            fontfamily="Times Roman")))



plot(vit.meg.spatial.plot1 )



# cam.xan

cam.xan.spatial.plot<-gg1 + geom_point(data = coord.cam.xan.limpo.names, 
                                       aes(x = LONG, y = LAT), 
                                       shape = 20, 
                                       color = "blue", 
                                       fill = "blue", 
                                       size = .5)



cam.xan.spatial.plot1 <- arrangeGrob(cam.xan.spatial.plot, 
                                     top = textGrob("D", x = unit(1, "npc"), 
                                                     y   = unit(-4, "npc"), just=c("right"),,
                                                    gp=gpar(col="black", 
                                                            fontsize=12, 
                                                            fontfamily="Times Roman")))



plot(cam.xan.spatial.plot1 )


# eug.uni

eug.uni.spatial.plot<-gg1 + geom_point(data = coord.eug.uni.limpo.names, 
                                       aes(x = LONG, y = LAT), 
                                       shape = 20, 
                                       color = "blue", 
                                       fill = "blue", 
                                       size = .5)

eug.uni.spatial.plot1 <- arrangeGrob(eug.uni.spatial.plot, 
                                     top = textGrob("E", x = unit(1, "npc"), 
                                                     y   = unit(-4, "npc"), just=c("right"),,
                                                    gp=gpar(col="black", 
                                                            fontsize=12, 
                                                            fontfamily="Times Roman")))
plot(eug.uni.spatial.plot1 )

#Group 2

# and.fra

and.fra.spatial.plot<-gg1 + geom_point(data = coord.and.fra.limpo.names, 
                                       aes(x = LONG, y = LAT), 
                                       shape = 20, 
                                       color = "green", 
                                       fill = "green", 
                                       size = .5)

and.fra.spatial.plot1 <- arrangeGrob(and.fra.spatial.plot, 
                                     top = textGrob("F", x = unit(1, "npc"), 
                                                     y   = unit(-4, "npc"), just=c("right"),,
                                                    gp=gpar(col="black", 
                                                            fontsize=12, 
                                                            fontfamily="Times Roman")))
plot(and.fra.spatial.plot1 )


# per.gla

per.gla.spatial.plot<-gg1 + geom_point(data = coord.per.gla.limpo.names, 
                                       aes(x = LONG, y = LAT), 
                                       shape = 20, 
                                       color = "green", 
                                       fill = "green", 
                                       size = .5)

per.gla.spatial.plot1 <- arrangeGrob(per.gla.spatial.plot, 
                                     top = textGrob("G", x = unit(1, "npc"), 
                                                     y   = unit(-4, "npc"), just=c("right"),,
                                                    gp=gpar(col="black", 
                                                            fontsize=12, 
                                                            fontfamily="Times Roman")))
plot(per.gla.spatial.plot1 )
# ani.fir

ani.fir.spatial.plot<-gg1 + geom_point(data = coord.ani.fir.limpo.names, 
                                       aes(x = LONG, y = LAT), 
                                       shape = 20, 
                                       color = "green", 
                                       fill = "green", 
                                       size = .5)

ani.fir.spatial.plot1 <- arrangeGrob(ani.fir.spatial.plot, 
                                     top = textGrob("H", x = unit(1, "npc"), 
                                                     y   = unit(-4, "npc"), just=c("right"),,
                                                    gp=gpar(col="black", 
                                                            fontsize=12, 
                                                            fontfamily="Times Roman")))
plot(ani.fir.spatial.plot1 )

# cec.gla

cec.gla.spatial.plot<-gg1 + geom_point(data = coord.cec.gla.limpo.names, 
                                       aes(x = LONG, y = LAT), 
                                       shape = 20, 
                                       color = "green", 
                                       fill = "green", 
                                       size = .5)

cec.gla.spatial.plot1 <- arrangeGrob(cec.gla.spatial.plot, 
                                     top = textGrob("I", x = unit(1, "npc"), 
                                                     y   = unit(-4, "npc"), just=c("right"),,
                                                    gp=gpar(col="black", 
                                                            fontsize=12, 
                                                            fontfamily="Times Roman")))
plot(cec.gla.spatial.plot1 )

# hye.alc

hye.alc.spatial.plot<-gg1 + geom_point(data = coord.hye.alc.limpo.names, 
                                       aes(x = LONG, y = LAT), 
                                       shape = 20, 
                                       color = "green", 
                                       fill = "green", 
                                       size = .5)

hye.alc.spatial.plot1 <- arrangeGrob(hye.alc.spatial.plot, 
                                     top = textGrob("J", x = unit(1, "npc"), 
                                                     y   = unit(-4, "npc"), just=c("right"),,
                                                    gp=gpar(col="black", 
                                                            fontsize=12, 
                                                            fontfamily="Times Roman")))
plot(hye.alc.spatial.plot1 )

# nec.opp

nec.opp.spatial.plot<-gg1 + geom_point(data = coord.nec.opp.limpo.names, 
                                       aes(x = LONG, y = LAT), 
                                       shape = 20, 
                                       color = "green", 
                                       fill = "green", 
                                       size = .5)

nec.opp.spatial.plot1 <- arrangeGrob(nec.opp.spatial.plot, 
                                     top = textGrob("K", x = unit(1, "npc"), 
                                                     y   = unit(-4, "npc"), just=c("right"),,
                                                    gp=gpar(col="black", 
                                                            fontsize=12, 
                                                            fontfamily="Times Roman")))
plot(nec.opp.spatial.plot1 )


# gua.aus

gua.aus.spatial.plot<-gg1 + geom_point(data = coord.gua.aus.limpo.names, 
                                       aes(x = LONG, y = LAT), 
                                       shape = 20, 
                                       color = "green", 
                                       fill = "green", 
                                       size = .5)

gua.aus.spatial.plot1 <- arrangeGrob(gua.aus.spatial.plot, 
                                     top = textGrob("L", x = unit(1, "npc"), 
                                                     y   = unit(-4, "npc"), just=c("right"),,
                                                    gp=gpar(col="black", 
                                                            fontsize=12, 
                                                            fontfamily="Times Roman")))
plot(gua.aus.spatial.plot1 )
# lec.pis

lec.pis.spatial.plot<-gg1 + geom_point(data = coord.lec.pis.limpo.names, 
                                       aes(x = LONG, y = LAT), 
                                       shape = 20, 
                                       color = "green", 
                                       fill = "green", 
                                       size = .5)

lec.pis.spatial.plot1 <- arrangeGrob(lec.pis.spatial.plot, 
                                     top = textGrob("M", x = unit(1, "npc"), 
                                                     y   = unit(-4, "npc"), just=c("right"),,
                                                    gp=gpar(col="black", 
                                                            fontsize=12, 
                                                            fontfamily="Times Roman")))
plot(lec.pis.spatial.plot1 )

# myr.rac

myr.rac.spatial.plot<-gg1 + geom_point(data = coord.myr.rac.limpo.names, 
                                       aes(x = LONG, y = LAT), 
                                       shape = 20, 
                                       color = "green", 
                                       fill = "green", 
                                       size = .5)

myr.rac.spatial.plot1 <- arrangeGrob(myr.rac.spatial.plot, 
                                     top = textGrob("N", x = unit(1, "npc"), 
                                                     y   = unit(-4, "npc"), just=c("right"),,
                                                    gp=gpar(col="black", 
                                                            fontsize=12, 
                                                            fontfamily="Times Roman")))
plot(myr.rac.spatial.plot1 )


# pse.con

pse.con.spatial.plot<-gg1 + geom_point(data = coord.pse.con.limpo.names, 
                                       aes(x = LONG, y = LAT), 
                                       shape = 20, 
                                       color = "green", 
                                       fill = "green", 
                                       size = .5)

pse.con.spatial.plot1 <- arrangeGrob(pse.con.spatial.plot, 
                                     top = textGrob("O", x = unit(1, "npc"), 
                                                     y   = unit(-4, "npc"), just=c("right"),,
                                                    gp=gpar(col="black", 
                                                            fontsize=12, 
                                                            fontfamily="Times Roman")))
plot(pse.con.spatial.plot1 )

# ann.dol

ann.dol.spatial.plot<-gg1 + geom_point(data = coord.ann.dol.limpo.names, 
                                       aes(x = LONG, y = LAT), 
                                       shape = 20, 
                                       color = "green", 
                                       fill = "green", 
                                       size = .5)

ann.dol.spatial.plot1 <- arrangeGrob(ann.dol.spatial.plot, 
                                     top = textGrob("P", x = unit(1, "npc"), 
                                                     y   = unit(-4, "npc"), just=c("right"),,
                                                    gp=gpar(col="black", 
                                                            fontsize=12, 
                                                            fontfamily="Times Roman")))
plot(ann.dol.spatial.plot1 )

# cou.mic

cou.mic.spatial.plot<-gg1 + geom_point(data = coord.cou.mic.limpo.names, 
                                       aes(x = LONG, y = LAT), 
                                       shape = 20, 
                                       color = "green", 
                                       fill = "green", 
                                       size = .5)

cou.mic.spatial.plot1 <- arrangeGrob(cou.mic.spatial.plot, 
                                     top = textGrob("Q", x = unit(1, "npc"), 
                                                     y   = unit(-4, "npc"), just=c("right"),,
                                                    gp=gpar(col="black", 
                                                            fontsize=12, 
                                                            fontfamily="Times Roman")))
plot(cou.mic.spatial.plot1 )


# eut.edu

eut.edu.spatial.plot<-gg1 + geom_point(data = coord.eut.edu.limpo.names, 
                                       aes(x = LONG, y = LAT), 
                                       shape = 20, 
                                       color = "green", 
                                       fill = "green", 
                                       size = .5)

eut.edu.spatial.plot1 <- arrangeGrob(eut.edu.spatial.plot, 
                                     top = textGrob("R", x = unit(1, "npc"), 
                                                     y   = unit(-4, "npc"), just=c("right"),,
                                                    gp=gpar(col="black", 
                                                            fontsize=12, 
                                                            fontfamily="Times Roman")))
plot(eut.edu.spatial.plot1 )


# gar.gar
gar.gar.spatial.plot<-gg1 + geom_point(data = coord.gar.gar.limpo.names, 
                                       aes(x = LONG, y = LAT), 
                                       shape = 20, 
                                       color = "green", 
                                       fill = "green", 
                                       size = .5)

gar.gar.spatial.plot1 <- arrangeGrob(gar.gar.spatial.plot, 
                                     top = textGrob("S", x = unit(1, "npc"), 
                                                     y   = unit(-4, "npc"), just=c("right"),,
                                                    gp=gpar(col="black", 
                                                            fontsize=12, 
                                                            fontfamily="Times Roman")))
plot(gar.gar.spatial.plot1 )
# pse.gra
pse.gra.spatial.plot<-gg1 + geom_point(data = coord.pse.gra.limpo.names, 
                                       aes(x = LONG, y = LAT), 
                                       shape = 20, 
                                       color = "green", 
                                       fill = "green", 
                                       size = .5)

pse.gra.spatial.plot1 <- arrangeGrob(pse.gra.spatial.plot, 
                                     top = textGrob("T", x = unit(1, "npc"), 
                                                     y   = unit(-4, "npc"), just=c("right"),,
                                                    gp=gpar(col="black", 
                                                            fontsize=12, 
                                                            fontfamily="Times Roman")))
plot(pse.gra.spatial.plot1 )
# alc.tri
alc.tri.spatial.plot<-gg1 + geom_point(data = coord.alc.tri.limpo.names, 
                                       aes(x = LONG, y = LAT), 
                                       shape = 20, 
                                       color = "green", 
                                       fill = "green", 
                                       size = .5)

alc.tri.spatial.plot1 <- arrangeGrob(alc.tri.spatial.plot, 
                                     top = textGrob("U", x = unit(1, "npc"), 
                                                     y   = unit(-4, "npc"), just=c("right"),,
                                                    gp=gpar(col="black", 
                                                            fontsize=12, 
                                                            fontfamily="Times Roman")))
plot(alc.tri.spatial.plot1 )

plot(and.fra.spatial.plot1)

tiff(filename = "plot1_new.tiff",
     width = 16.5, height = 22, units = "cm",
     bg = "white", 
     res = 300)

grid.arrange(gym.klo.spatial.plot1 , all.edu.spatial.plot1 , vit.meg.spatial.plot1, cam.xan.spatial.plot1, 
             eug.uni.spatial.plot1,and.fra.spatial.plot1,per.gla.spatial.plot1,ani.fir.spatial.plot1,
             cec.gla.spatial.plot1, hye.alc.spatial.plot1, nec.opp.spatial.plot1, gua.aus.spatial.plot1,
             lec.pis.spatial.plot1,myr.rac.spatial.plot1,pse.con.spatial.plot1 ,ann.dol.spatial.plot1,
             cou.mic.spatial.plot1, eut.edu.spatial.plot1,
             gar.gar.spatial.plot1, pse.gra.spatial.plot1,alc.tri.spatial.plot1 ,nrow = 6, ncol=4)

dev.off()

tiff(filename = "plot1_new.tiff",
     width = 16.5, height = 22, units = "cm",
     bg = "white", 
     res = 300)

grid.arrange(gym.klo.spatial.plot1 , all.edu.spatial.plot1 , vit.meg.spatial.plot1 ,
             cam.xan.spatial.plot1 ,eug.uni.spatial.plot1,and.fra.spatial.plot1,
             per.gla.spatial.plot1,ani.fir.spatial.plot1, cec.gla.spatial.plot1 ,
             hye.alc.spatial.plot1, nec.opp.spatial.plot1, gua.aus.spatial.plot1,
             lec.pis.spatial.plot1,myr.rac.spatial.plot1,pse.con.spatial.plot1 ,
             ann.dol.spatial.plot1 ,cou.mic.spatial.plot1, eut.edu.spatial.plot1,
             gar.gar.spatial.plot1, pse.gra.spatial.plot1,alc.tri.spatial.plot1 ,  nrow = 6, ncol=4)


dev.off()

