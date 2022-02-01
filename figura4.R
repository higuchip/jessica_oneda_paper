library(ggplot2)
quantiles_95 <- function(x) {
  r <- quantile(x, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

efeito_clima_new<-c(-45.904, -49.638, -57.539, -52.394 , -34.792,
                    85.244, -76.316, -60.317,  -86.639, -33.066,
                    -58.418,  -83.39, -15.285, -2.048, -78.149, 
                    -69.951, -4.059, -77.375, -35.047, -67.055,
                    -59.858 )
length(efeito_clima_new)


df$efeito_clima_new<-c(-45.904, -49.638, -57.539, -52.394 , -34.792,
                       85.244, -76.316, -60.317,  -86.639, -33.066,
                       -58.418,  -83.39, -15.285, -2.048, -78.149, 
                       -69.951, -4.059, -77.375, -35.047, -67.055,
                       -59.858 )
dim(df)
df

tiff(filename = "plot4_new1.tiff",
     width = 10, height = 10, units = "cm",
     bg = "white", 
     res = 300)


p1<-ggplot(df, aes(x=bacia, y=efeito_clima_new, fill=bacia)) +
  guides(fill=F) +
  stat_summary(fun.data = quantiles_95, geom="boxplot") + xlab("Basin") + ylab("Decrease in area of climate adequacy (%)")

p1 + geom_jitter(shape=16, position=position_jitter(0.2))

dev.off()

shapiro.test(df$efeito_clima_new)
hist(df$efeito_clima_new)
wilcox.test(df$efeito_clima_new~df$bacia)
wilcox.test(df$efeito_clima~df$bacia)
