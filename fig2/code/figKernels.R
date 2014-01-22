source('parametersCommon.R')
if(!require(plotrix)){install.packages('plotrix'); library(plotrix)}


pF = points.poisson(Nf)
XFcoord = pF$x # blood-feeding habitat coordinate X
YFcoord = pF$y # blood-feeding habitat coordinate Y
XLcoord = rnorm(Nl, XFcoord[larvalLocations], .01) # larval habitat coordinate X
YLcoord = rnorm(Nl, YFcoord[larvalLocations], .01) # larval habitat coordinate Y

L = matrix(0, Nf, Nl)
F = matrix(0, Nl, Nf)
L = mosquitoMovement(XFcoord, YFcoord, XLcoord, YLcoord, surv = survL, mixing = 'poor')
F = mosquitoMovement(XLcoord, YLcoord, XFcoord, YFcoord, surv = survF, mixing = 'poor')
H = makeHpoorHostMixing(Nf, Nh, Hh)
U = makeUrossMacdonald(H)
Nh = dim(H)[1]
source('variables.R')
source('calculateMetrics.R')
source('calculateKernels.R')



pdf(file = '../output/fig2.pdf', width = 6, height = 6)

par(oma = c(3, 5, 1, 0), mar = c(1.5, 1.5, .5, 1.5), las = 1)
layout(matrix(1 : 4, 2, 2))

plot(-1,
		 xlim = range(bins), ylim = c(0, .03),
		 xlab = '', ylab = '',
		 axes = FALSE)
maxQ = max(Q_pdf.emp) - max(Q_pdf.emp) %% 0.01
meanDQ = sum(Q_pdf.emp * bins)
Q_pdf.emp[which.max(Q_pdf.emp)] = .03
segments(bins, rep(0, length(bins)), bins, Q_pdf.emp, lwd = 3, col = 'gray')
lines(Q_pdf.sth$x, Q_pdf.sth$y, lwd = 3)
segments(meanDQ, 0, meanDQ, 0.03, col = 1, lty = 2, lwd = 3)
axis(1, seq(0, 1, .2), labels = FALSE)
axis(2, seq(0, .03, .01), labels = c('0.0', '0.01', '0.02', as.character(maxQ)))
axis.break(2, .028, style = 'slash', brw = .03)
mtext('(a)', side = 3, adj = 0, line = -.3, cex = 1)
text(x = meanDQ + .01, y = .0285, pos = 4,
  labels = paste('E{    } =', as.character(meanDQ - meanDQ %% 0.01)))
text(x = meanDQ + .01 + .08, y = .0282, pos = 4,
  labels = expression(delta[Q]))



plot(-1,
		 xlim = range(bins), ylim = c(0, .03),
		 xlab = '', ylab = '',
		 axes = FALSE)
maxZ = max(Z_pdf.emp) - max(Z_pdf.emp) %% 0.01
meanDZ = sum(Z_pdf.emp * bins)
Z_pdf.emp[which.max(Z_pdf.emp)] = .03
segments(bins, rep(0, length(bins)), bins, Z_pdf.emp, lwd = 3, col = 'gray')
lines(Z_pdf.sth$x, Z_pdf.sth$y, lwd = 3)
segments(meanDZ, 0, meanDZ, 0.03, col = 1, lty = 2, lwd = 3)
axis(1, seq(0, 1, .2))
axis(2, seq(0, .03, .01), labels = c('0.0', '0.01', '0.02', as.character(maxZ)))
axis.break(2, .028, style = 'slash', brw = .03)
mtext('(c)', side = 3, adj = 0, line = -.3, cex = 1)
mtext('Distance', side = 1, adj = 0, line = 2.75, at = 1.01)
par(las = 0)
mtext('Kernel density', side = 2, adj = 0, line = 2.75, at = .42, outer = TRUE)
par(las = 1)
text(x = meanDZ + .01, y = .0285, pos = 4,
  labels = paste('E{    } =', as.character(meanDZ - meanDZ %% 0.01)))
text(x = meanDZ + .01 + .08, y = .0282, pos = 4,
  labels = expression(delta[Z]))



plot(-1,
		 xlim = range(bins), ylim = c(0, .03),
		 xlab = '', ylab = '',
		 axes = FALSE)
maxP = max(P_pdf.emp) - max(P_pdf.emp) %% 0.01
meanDP = sum(P_pdf.emp * bins)
P_pdf.emp[which.max(P_pdf.emp)] = .03
segments(bins, rep(0, length(bins)), bins, P_pdf.emp, lwd = 3, col = 'gray')
lines(P_pdf.sth$x, P_pdf.sth$y, lwd = 3)
segments(meanDP, 0, meanDP, 0.03, col = 1, lty = 2, lwd = 3)
axis(1, seq(0, 1, .2), labels = FALSE)
axis(2, seq(0, .03, .01), labels = c('0.0', '0.01', '0.02', as.character(maxP)))
axis.break(2, .028, style = 'slash', brw = .03)
mtext('(b)', side = 3, adj = 0, line = -.3, cex = 1)
text(x = meanDP + .01, y = .0285, pos = 4,
  labels = paste('E{    } =', as.character(meanDP - meanDP %% 0.01)))
text(x = meanDP + .01 + .08, y = .0282, pos = 4,
  labels = expression(delta[P]))



plot(-1,
		 xlim = range(bins), ylim = c(0, .03),
		 xlab = '', ylab = '',
		 axes = FALSE)
maxS = max(S_pdf.emp) - max(S_pdf.emp) %% 0.01
meanDS = sum(S_pdf.emp * bins)
S_pdf.emp[which.max(S_pdf.emp)] = .03
segments(bins, rep(0, length(bins)), bins, S_pdf.emp, lwd = 3, col = 'gray')
lines(S_pdf.sth$x, S_pdf.sth$y, lwd = 3)
segments(meanDS, 0, meanDS, 0.03, col = 1, lty = 2, lwd = 3)
axis(1, seq(0, 1, .2))
axis(2, seq(0, .03, .01), labels = c('0.0', '0.01', '0.02', as.character(maxS)))
axis.break(2, .028, style = 'slash', brw = .03)
mtext('(d)', side = 3, adj = 0, line = -.3, cex = 1)
text(x = meanDS + .01, y = .0285, pos = 4,
  labels = paste('E{    } =', as.character(meanDS - meanDS %% 0.01)))
text(x = meanDS + .01 + .08, y = .0282, pos = 4,
  labels = expression(delta[S]))

dev.off()
