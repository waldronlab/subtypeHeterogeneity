############################################################
# 
# author: Ludwig Geistlinger
# date: 2018-10-17 17:01:07
# 
# descr: 
# 
############################################################

library(ggplot2)

x <- rep(c("DIF", "IMR", "MES", "PRO"), each=100000)
y <- c(rnorm(100000, mean=0), rnorm(100000, mean=3), 
            rnorm(100000, mean=4), rnorm(100000, mean=7))
df <- data.frame(x=x, y=y)
ggplot(df, aes(y, fill=x, colour=x)) + 
    geom_density(alpha = 0.3) + 
    theme_classic() + 
    ylim(0,4) +
    scale_fill_manual(values=c(cb.orange, cb.pink, cb.green, cb.lightblue)) +
    scale_colour_manual(values=c(cb.orange, cb.pink, cb.green, cb.lightblue))  
