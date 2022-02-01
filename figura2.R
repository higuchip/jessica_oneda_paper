library(grid)
library(gridExtra)



# Relacao de spp com > 100 pontos de ocorrencias filtrados e TSS > 0.4, disponivel no BIEN

#Group 1
# gym_klo




gym.klo.bio <- arrangeGrob(p_gym.klo_bio8, 
                                     top = textGrob("A", x = unit(1, "npc"), 
                                                    y   = unit(-8, "npc"), just=c("right"),
                                                    gp=gpar(col="black", 
                                                            fontsize=18, 
                                                            fontfamily="Times Roman")))



plot(gym.klo.bio )



# all.edu



all.edu.bio <- arrangeGrob(p_all.edu_bio3, 
                           top = textGrob("B", x = unit(1, "npc"), 
                                          y   = unit(-8, "npc"), just=c("right"),
                                          gp=gpar(col="black", 
                                                  fontsize=18, 
                                                  fontfamily="Times Roman")))



plot(all.edu.bio )


# vit.meg

vit.meg.bio <- arrangeGrob(p_vit.meg_bio10, 
                           top = textGrob("C", x = unit(1, "npc"), 
                                          y   = unit(-8, "npc"), just=c("right"),
                                          gp=gpar(col="black", 
                                                  fontsize=18, 
                                                  fontfamily="Times Roman")))



plot(vit.meg.bio )

# cam.xan

cam.xan.bio <- arrangeGrob(p_cam.xan_bio8, 
                           top = textGrob("D", x = unit(1, "npc"), 
                                          y   = unit(-8, "npc"), just=c("right"),
                                          gp=gpar(col="black", 
                                                  fontsize=18, 
                                                  fontfamily="Times Roman")))



plot(cam.xan.bio )

# eug.uni

eug.uni.bio <- arrangeGrob(p_eug.uni_bio3, 
                           top = textGrob("E", x = unit(1, "npc"), 
                                          y   = unit(-8, "npc"), just=c("right"),
                                          gp=gpar(col="black", 
                                                  fontsize=18, 
                                                  fontfamily="Times Roman")))



plot(eug.uni.bio )

#Group 2

# and.fra

and.fra.bio <- arrangeGrob(p_and.fra_bio15, 
                           top = textGrob("F", x = unit(1, "npc"), 
                                          y   = unit(-8, "npc"), just=c("right"),
                                          gp=gpar(col="black", 
                                                  fontsize=18, 
                                                  fontfamily="Times Roman")))



plot(and.fra.bio )


# per.gla

per.gla.bio <- arrangeGrob(p_per.gla_bio8, 
                           top = textGrob("G", x = unit(1, "npc"), 
                                          y   = unit(-8, "npc"), just=c("right"),
                                          gp=gpar(col="black", 
                                                  fontsize=18, 
                                                  fontfamily="Times Roman")))



plot(per.gla.bio )


# ani.fir

ani.fir.bio <- arrangeGrob(p_ani.fir_bio14, 
                           top = textGrob("H", x = unit(1, "npc"), 
                                          y   = unit(-8, "npc"), just=c("right"),
                                          gp=gpar(col="black", 
                                                  fontsize=18, 
                                                  fontfamily="Times Roman")))



plot(ani.fir.bio )

# cec.gla

cec.gla.bio <- arrangeGrob(p_cec.gla_bio3, 
                           top = textGrob("I", x = unit(1, "npc"), 
                                          y   = unit(-8, "npc"), just=c("right"),
                                          gp=gpar(col="black", 
                                                  fontsize=18, 
                                                  fontfamily="Times Roman")))



plot(cec.gla.bio )

# hye.alc

hye.alc.bio <- arrangeGrob(p_hye.alc_bio9, 
                           top = textGrob("J", x = unit(1, "npc"), 
                                          y   = unit(-8, "npc"), just=c("right"),
                                          gp=gpar(col="black", 
                                                  fontsize=18, 
                                                  fontfamily="Times Roman")))



plot(hye.alc.bio )

# nec.opp

nec.opp.bio <- arrangeGrob(p_nec.opp_bio3, 
                           top = textGrob("K", x = unit(1, "npc"), 
                                          y   = unit(-8, "npc"), just=c("right"),
                                          gp=gpar(col="black", 
                                                  fontsize=18, 
                                                  fontfamily="Times Roman")))



plot(nec.opp.bio )


# gua.aus

gua.aus.bio <- arrangeGrob(p_gua.aus_bio14, 
                           top = textGrob("L", x = unit(1, "npc"), 
                                          y   = unit(-8, "npc"), just=c("right"),
                                          gp=gpar(col="black", 
                                                  fontsize=18, 
                                                  fontfamily="Times Roman")))



plot(gua.aus.bio )
# lec.pis

lec.pis.bio <- arrangeGrob(p_lec.pis_bio14, 
                           top = textGrob("M", x = unit(1, "npc"), 
                                          y   = unit(-8, "npc"), just=c("right"),
                                          gp=gpar(col="black", 
                                                  fontsize=18, 
                                                  fontfamily="Times Roman")))



plot(lec.pis.bio )
# myr.rac

myr.rac.bio <- arrangeGrob(p_myr.rac_bio2, 
                           top = textGrob("N", x = unit(1, "npc"), 
                                          y   = unit(-8, "npc"), just=c("right"),
                                          gp=gpar(col="black", 
                                                  fontsize=18, 
                                                  fontfamily="Times Roman")))



plot(myr.rac.bio )

# pse.con

pse.con.bio <- arrangeGrob(p_pse.con_bio14, 
                           top = textGrob("O", x = unit(1, "npc"), 
                                          y   = unit(-8, "npc"), just=c("right"),
                                          gp=gpar(col="black", 
                                                  fontsize=18, 
                                                  fontfamily="Times Roman")))



plot(pse.con.bio )

# ann.dol

ann.dol.bio <- arrangeGrob(p_ann.dol_bio14, 
                           top = textGrob("P", x = unit(1, "npc"), 
                                          y   = unit(-8, "npc"), just=c("right"),
                                          gp=gpar(col="black", 
                                                  fontsize=18, 
                                                  fontfamily="Times Roman")))



plot(ann.dol.bio )

# cou.mic

cou.mic.bio <- arrangeGrob(p_cou.mic_bio3, 
                           top = textGrob("Q", x = unit(1, "npc"), 
                                          y   = unit(-8, "npc"), just=c("right"),
                                          gp=gpar(col="black", 
                                                  fontsize=18, 
                                                  fontfamily="Times Roman")))



plot(cou.mic.bio )

# eut.edu

eut.edu.bio <- arrangeGrob(p_eut.edu_bio14, 
                           top = textGrob("R", x = unit(1, "npc"), 
                                          y   = unit(-8, "npc"), just=c("right"),
                                          gp=gpar(col="black", 
                                                  fontsize=18, 
                                                  fontfamily="Times Roman")))



plot(eut.edu.bio )


# gar.gar
gar.gar.bio <- arrangeGrob(p_gar.gar_bio3, 
                           top = textGrob("S", x = unit(1, "npc"), 
                                          y   = unit(-8, "npc"), just=c("right"),
                                          gp=gpar(col="black", 
                                                  fontsize=18, 
                                                  fontfamily="Times Roman")))



plot(gar.gar.bio )

# pse.gra
pse.gra.bio <- arrangeGrob(p_pse.gra_bio14, 
                           top = textGrob("T", x = unit(1, "npc"), 
                                          y   = unit(-8, "npc"), just=c("right"),
                                          gp=gpar(col="black", 
                                                  fontsize=18, 
                                                  fontfamily="Times Roman")))



plot(pse.gra.bio )
# alc.tri
alc.tri.bio <- arrangeGrob(p_alc.tri_bio8, 
                           top = textGrob("U", x = unit(1, "npc"), 
                                          y   = unit(-8, "npc"), just=c("right"),
                                          gp=gpar(col="black", 
                                                  fontsize=18, 
                                                  fontfamily="Times Roman")))



plot(alc.tri.bio )



tiff(filename = "plot2_new.tiff",
     width = 3300, height = 4800, units = "px",
     bg = "white", 
     res = 300)

grid.arrange(gym.klo.bio , all.edu.bio , vit.meg.bio ,
             cam.xan.bio ,eug.uni.bio,and.fra.bio,
             per.gla.bio,ani.fir.bio, cec.gla.bio ,
             hye.alc.bio, nec.opp.bio, gua.aus.bio,
             lec.pis.bio,myr.rac.bio,pse.con.bio ,
             ann.dol.bio ,cou.mic.bio, eut.edu.bio,
             gar.gar.bio, pse.gra.bio,alc.tri.bio ,  nrow = 6, ncol=4)


dev.off()

