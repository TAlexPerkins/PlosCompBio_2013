if(!require(RColorBrewer)){install.packages('RColorBrewer'); library(RColorBrewer)}
if(!require(spatstat)){install.packages('spatstat'); library(spatstat)}
load('../data/parametersCommon.RData')
source('functionsPhylogram.R')



pdf(file = '../output/heterogeneity.pdf', width = 5.5, height = 7)
par()
layout(mat = cbind(c(1, 4 : 7), c(2, 8 : 11), c(3, 12 : 15)), heights = c(rep(c(1, 3, 3, 3, 3), 3)))
par(oma = c(4, 5, 3, 3), mar = c(0, 0, 2, 0))



plot(-1, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, frame = FALSE, xlab = '', ylab = '')
rect(seq(0, .9, .1), rep(0, 10), seq(.1, 1, .1), rep(1, 10), border = c(rev(brewer.pal(8, 'Blues'))[1 : 5], brewer.pal(8, 'Reds')[4 : 8]), col = c(rev(brewer.pal(8, 'Blues'))[1 : 5], brewer.pal(8, 'Reds')[4 : 8]))
axis(3, at = c(0, .2, .4, .6, .8, 1), labels = c('0', '0.2', '0.4', '0.6', '0.8', '1'))
plot(-1, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, frame = FALSE, xlab = '', ylab = '')
rect(seq(0, .9, .1), rep(0, 10), seq(.1, 1, .1), rep(1, 10), border = c(rev(brewer.pal(8, 'Blues'))[1 : 5], brewer.pal(8, 'Reds')[4 : 8]), col = c(rev(brewer.pal(8, 'Blues'))[1 : 5], brewer.pal(8, 'Reds')[4 : 8]))
axis(3, at = c(0, .2, .4, .6, .8, 1), labels = c('0', '0.2', '0.4', '0.6', '0.8', '1'))
mtext('Heterogeneity', side = 3, line = 2.3, cex = 1)
plot(-1, xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, frame = FALSE, xlab = '', ylab = '')
rect(seq(0, .9, .1), rep(0, 10), seq(.1, 1, .1), rep(1, 10), border = c(rev(brewer.pal(8, 'Blues'))[1 : 5], brewer.pal(8, 'Reds')[4 : 8]), col = c(rev(brewer.pal(8, 'Blues'))[1 : 5], brewer.pal(8, 'Reds')[4 : 8]))
axis(3, at = c(0, .2, .4, .6, .8, 1), labels = c('0', '0.2', '0.4', '0.6', '0.8', '1'))



par(mar = rep(0, 4))

load('../../fig6/data/overdispersed.RData')
source('functionsPhylogram.R')
source('calculateMetrics.R')

p = as.phylo(hclust(dist(cbind(XFcoord, YFcoord))))
plotHetBySize(p, S)
mtext('Patch size', side = 2, line = 3, cex = 1)
axis(side = 2, line = 0, las = 1)

p = as.phylo(hclust(dist(cbind(XFcoord, YFcoord))))
plotHetByDist(p, S)
mtext('Spatial extent', side = 2, line = 3, cex = 1)
axis(side = 2, line = 0, las = 1)

p = as.phylo(hclust(dist(cbind(XFcoord, YFcoord))))
plotPhyloHet(p, S)

par(mar = c(0, 0, 1, 0))
plot(disc(radius = sqrt(1 / pi)), main = '')
points(XFcoord, YFcoord, pch = 19)
mtext('Overdispersed', side = 1, line = .5, cex = 1)
axis(side = 2, line = 0, las = 2)
mtext(side = 2, line = 3, 'Distance', cex = 1)



par(mar = rep(0, 4))

load('../../fig6/data/poisson.RData')
source('functionsPhylogram.R')
source('calculateMetrics.R')

p = as.phylo(hclust(dist(cbind(XFcoord, YFcoord))))
plotHetBySize(p, S)

p = as.phylo(hclust(dist(cbind(XFcoord, YFcoord))))
plotHetByDist(p, S)

p = as.phylo(hclust(dist(cbind(XFcoord, YFcoord))))
plotPhyloHet(p, S)

par(mar = c(0, 0, 1, 0))
plot(disc(radius = sqrt(1 / pi)), main = '')
points(XFcoord, YFcoord, pch = 19)
mtext('Poisson', side = 1, line = .5, cex = 1)
mtext('Spatial arrangement of blood-feeding habitats', side = 1, line = 2.5, cex = 1)



par(mar = rep(0, 4))

load('../../fig6/data/clustered.RData')
source('functionsPhylogram.R')
source('calculateMetrics.R')

p = as.phylo(hclust(dist(cbind(XFcoord, YFcoord))))
plotHetBySize(p, S)

p = as.phylo(hclust(dist(cbind(XFcoord, YFcoord))))
plotHetByDist(p, S)

p = as.phylo(hclust(dist(cbind(XFcoord, YFcoord))))
plotPhyloHet(p, S)

par(mar = c(0, 0, 1, 0))
plot(disc(radius = sqrt(1 / pi)), main = '')
points(XFcoord, YFcoord, pch = 19)
mtext('Clustered', side = 1, line = .5, cex = 1)

dev.off()
