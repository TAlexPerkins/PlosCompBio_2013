load('../data/parametersCommon.RData')



pF = points.clustered(Nf)
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

save(list = ls(), file = '../data/clustered.RData')



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

save(list = ls(), file = '../data/poisson.RData')



pF = points.overdispersed(Nf)
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

save(list = ls(), file = '../data/overdispersed.RData')

