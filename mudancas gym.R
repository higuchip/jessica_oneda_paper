jpeg(filename = "mudancas_gk.jpg",
     width = 2000, height = 3500, units = "px",
     quality = 100,
     bg = "white", 
     res = 300,
     family="times")
par(mfcol=c(2,1),mar=c(3.5,4,1.5,0), oma=c(1.5,1,1.5,1), 
    mgp = c(2, 1, 0), family= "serif")
plot(myBiomodRangeSize.drimys.2070.45.gk$Diff.By.Pixel, col= gray.colors(4, start = 0, end = 1, gamma = 1, alpha = NULL))
plot(wrld_simpl,  axes="TRUE", add=T)
text(-50,-20, "BRASIL", cex=1)
text(-45,-30, "OCEANO ATLÂNTICO", cex=1)
mtext("a", side=1, line=3, adj=1, font=6, cex=1.5)

plot(myBiomodRangeSize.drimys.2070.85.gk$Diff.By.Pixel, col= gray.colors(4, start = 0, end = 1, gamma = 1, alpha = NULL))
plot(wrld_simpl,  axes="TRUE", add=T)
text(-50,-20, "BRASIL", cex=1)
text(-45,-30, "OCEANO ATLÂNTICO", cex=1)
mtext("b", side=1, line=3, adj=1, font=6, cex=1.5)

dev.off()


jpeg(filename = "mudancas_af.jpg",
     width = 2000, height = 3500, units = "px",
     quality = 100,
     bg = "white", 
     res = 300,
     family="times")
par(mfcol=c(2,1),mar=c(3.5,4,1.5,0), oma=c(1.5,1,1.5,1), 
    mgp = c(2, 1, 0), family= "serif")
plot(myBiomodRangeSize.drimys.2070.45.af$Diff.By.Pixel, col= gray.colors(4, start = 0, end = 1, gamma = 1, alpha = NULL))
plot(wrld_simpl,  axes="TRUE", add=T)
text(-50,-20, "BRASIL", cex=1)
text(-45,-30, "OCEANO ATLÂNTICO", cex=1)
mtext("a", side=1, line=3, adj=1, font=6, cex=1.5)

plot(myBiomodRangeSize.drimys.2070.85.af$Diff.By.Pixel, col= gray.colors(4, start = 0, end = 1, gamma = 1, alpha = NULL))
plot(wrld_simpl,  axes="TRUE", add=T)
text(-50,-20, "BRASIL", cex=1)
text(-45,-30, "OCEANO ATLÂNTICO", cex=1)
mtext("b", side=1, line=3, adj=1, font=6, cex=1.5)

dev.off()