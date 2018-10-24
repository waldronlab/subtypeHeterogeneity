############################################################
# 
# author: Ludwig Geistlinger
# date: 2018-09-06 16:35:38
# 
# descr: plotly helper figures
# 
############################################################

library(plotly)

## TUMOR COMPOSITION

subcl <- rep(0, 40)
subcl2 <- cumsum(1:12)
subcl2 <- seq(subcl2[1], subcl2[12], length=60) 
subcl <- c(subcl, subcl2)

subcl[45:95] <- subcl[45:95] + runif(51, min=-15, max=15)
subcl[subcl < 0] <- 1 
subcl[subcl > 80] <- 80 
#ind <- seq(60, 100, by=3)


plot_ly(x = 1:100, y = 100, type = 'scatter', name = "Founder clone",
                    mode = 'none', fill = 'tozeroy', fillcolor = '#F5FF8D') %>%
  add_trace(y = subcl , name = 'Subclone', fillcolor = '#50CB86') %>%
  layout(title = "",
         xaxis = list(visible = FALSE),
         yaxis = list(visible = FALSE))


## Correlation subtype association / subclonality

plot_ly(x = 1:100, y = 100:1, type="scatter", mode="markers") %>%
    layout(
         xaxis = list(title = "Subtype association",
                        showticklabels = FALSE,
                        showgrid = FALSE),
         yaxis = list(title = "Subclonality",
                        showticklabels = FALSE,
                        showgrid = FALSE))

barplot(c(25,25,25,25) + runif(4, min=-5, max=5), 
    names=names(stcols)[c(3,4,2,1)], 
    col=stcols[c(3,4,2,1)], 
    ylab="%tumors", ylim=c(0,100))   

barplot(c(15,15,15,55) + runif(4, min=-5, max=5), 
    names=names(stcols)[c(3,4,2,1)], 
    col=stcols[c(3,4,2,1)], 
    ylab="%tumors", ylim=c(0,100))   

## pie charts
cb.pink <- "#CC79A7"
cb.red <- "#D55E00"
cb.blue <- "#0072B2"
cb.yellow <- "#F0E442"
cb.green <- "#009E73"
cb.lightblue <- "#56B4E9"
cb.orange <- "#E69F00"

stcols <- c(cb.lightblue, cb.green, cb.orange, cb.pink) 
names(stcols) <- c("PRO", "MES", "DIF", "IMR")

# 99genes, all cells, epithelial: c(1,3,4), c(0.027, 0.054, 0.919) 
# 99genes, all cells, stromal: c(2,4), c(0.207, 0.793) 
# 99genes, top50 cells, epithelial: 4, 1
# 99genes, top50 cells, stromal: c(2,4), c(0.25, 0.75)

# 92genes, all cells, epithelial: c(1,3,4), c(0.054, 0.324, 0.622) 
# 92genes, all cells, stromal: 1:4, c(0.069, 0.31, 0.207, 0.414) 
# 92genes, top50 cells, epithelial: 3:4, c(0.263, 0.737)
# 92genes, top50 cells, stromal: 1:4, c(0.067, 0.333, 0.133, 0.467)

p <- plot_ly(labels = names(stcols)[3:4], values = c(0.263, 0.737), type = 'pie',
       textposition = 'inside',
       textinfo = 'label+percent',
       insidetextfont = list(color = '#FFFFFF', size=30),
        pull = 0.01,
        hole = 0.01, 
       marker = list(colors = stcols[3:4],
                     line = list(color = '#FFFFFF', width = 1)),
       showlegend = FALSE) %>%
 layout(
        xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
        yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

Sys.setenv('MAPBOX_TOKEN' = "pk.eyJ1IjoibHVkd2lnZyIsImEiOiJjamx3cXFwMm8xOGJ3M2tvZGV5amozNG5rIn0.-EpukkiU2QxcbUvFtq1QOw")
orca(p, "pie-plot.pdf")


