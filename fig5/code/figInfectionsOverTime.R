if(!require(MASS)){install.packages('MASS'); library(MASS)}
if(!require(RColorBrewer)){install.packages('RColorBrewer'); library(RColorBrewer)}

fileList = sapply(1 : length(list.files('../data/')), function(ii)
  substr(list.files('../data/')[ii], start = 1, stop = nchar(list.files('../data/')[ii]) - nchar('.RData')))


nFiles = length(fileList) - 1

for(ff in 1 : nFiles)
{
  load(paste('../data/', fileList[ff], '.RData', sep = ''))
} # end for


cumHumanInfections_000 = matrix(0, Nreps, T + 1)
humanInfections_000 = matrix(0, Nreps, T + 1)
cumHumanInfections_110 = matrix(0, Nreps, T + 1)
humanInfections_110 = matrix(0, Nreps, T + 1)
for(ii in 1 : Nreps)
{
  cumHumanInfections_000[ii, ] = data_000[1, ii][[1]]
  humanInfections_000[ii, ] = data_000[2, ii][[1]]
  cumHumanInfections_110[ii, ] = data_110[1, ii][[1]]
  humanInfections_110[ii, ] = data_110[2, ii][[1]]
} # end for



pdf(file = '../output/fig5.pdf', width = 6, height = 3.1)

# figure set up
par(oma = c(2.5, 2.5, .25, 0), mar = c(1, 1, 1.25, 1))
layout(rbind(c(1, 2), c(3, 2)), widths = c(.45, .55))
branch.colors = c(rev(brewer.pal(8, 'Blues'))[1 : 5], brewer.pal(8, 'Reds')[4 : 8])
lineColor_000 = rgb(col2rgb(branch.colors[10])[1], col2rgb(branch.colors[10])[2], col2rgb(branch.colors[10])[3], alpha = 60, maxColorValue = 255)
lineColor_110 = rgb(col2rgb(branch.colors[1])[1], col2rgb(branch.colors[1])[2], col2rgb(branch.colors[1])[3], alpha = 60, maxColorValue = 255)

# plot panel a
plot(secondaryInfections_000[1, ], type = 'l', lwd = 5, col = lineColor_000, ylim = c(0, max(secondaryInfections_110)), axes = FALSE)
axis(2, at = 0 : max(secondaryInfections_110), las = 1)
axis(1, at = seq(0, dim(secondaryInfections_000)[2], 5), labels = FALSE)
for(ii in (2 : dim(secondaryInfections_000)[1])[-1])
{
  lines(secondaryInfections_000[ii, ], type = 'l', lwd = 5, col = lineColor_000)
}
lines(secondaryInfections_000[1, ], type = 'l', lwd = 5, col = lineColor_000)
mtext('(a)', 3, line = .1, adj = 0)

# plot panel c main
par(mar = par()$mar + c(0, 2.5, 0, 0))
plot(humanInfections_000[1, 1 : 250], type = 'l', col = lineColor_000, ylim = c(0, .8 * Nh), xlab = 'Time', ylab = '', las = 1)
lines(humanInfections_110[1, 1 : 250], type = 'l', col = lineColor_110)
for(ii in 2 : Nreps)
{
  lines(humanInfections_000[ii, 1 : 250], type = 'l', col = lineColor_000)
  lines(humanInfections_110[ii, 1 : 250], type = 'l', col = lineColor_110)
} # end for
mtext('Infected hosts', 2, line = 2.75)
mtext('Time', 1, line = 2.5)
mtext('(c)', 3, line = .1, adj = 0)

# plot panel c inset
par(new = T, mar = c(9, 9, 1.25, 1))
plot(cumHumanInfections_000[1, 1 : 250], type = 'l', col = lineColor_000, ylim = c(0, Nh), frame.plot = TRUE, axes = FALSE, xlab = '', ylab = '')
lines(cumHumanInfections_110[1, 1 : 250], type = 'l', col = lineColor_110)
for(ii in 2 : Nreps)
{
  lines(cumHumanInfections_000[ii, 1 : 250], type = 'l', col = lineColor_000)
  lines(cumHumanInfections_110[ii, 1 : 250], type = 'l', col = lineColor_110)
} # end for
mtext('Time', 1, line = 0, cex = .8)
mtext('Cumulative\nnew infections', 2, line = 0, cex = .8)

# plot panel b
par(new = FALSE, mar = c(1, 1, 1.25, 1))
plot(secondaryInfections_110[humanMeanR, ], type = 'l', lwd = 5, col = lineColor_110, ylim = c(0, max(secondaryInfections_110)), axes = FALSE)
axis(2, at = 0 : max(secondaryInfections_110), las = 1)
axis(1, at = seq(0, dim(secondaryInfections_110)[2], 5))
for (ii in (2 : dim(secondaryInfections_110)[1])[-humanMeanR])
{
  lines(secondaryInfections_110[ii, ], type = 'l', lwd = 5, col = lineColor_110)
}
lines(secondaryInfections_110[humanMeanR, ], type = 'l', lwd = 5, col = lineColor_110)
mtext('Cumulative infectious bites per host', 2, line = 1, at = .477, outer = TRUE)
mtext('Time', 1, line = 2.5)
mtext('(b)', 3, line = .1, adj = 0)

dev.off()
