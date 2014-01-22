library('stats4')

source('functions.R')

T = 400 # number of time steps to run
Nf = 100 # number of feeding stations
Nl = 100 # number of larval stations
v = 10 # number of mosquito eggs per capita
xi = 4 # larval stages for maturation

survL = 0.93
survF = 0.93

sigma = 3 # mosquito extrinsic incubation period
tau = 3 # host intrinsic incubation period
rho_max = 40 # rho is host recovery and these lines define its distribution
rho_pmf = dnbinom(0 : (rho_max - 1), 3, .2)
rho_mean = sum((1 : rho_max) * rho_pmf)
rho_fail = rho_pmf / (1 - c(0, pnbinom(0 : (rho_max - 2), 3, .2)))
c = 0.15 # host-to-mosquito transmission efficiency
b = 0.55 # mosquito-to-host transmission efficiency

pF = points.poisson(Nf)
larvalLocations = sample(1 : Nf, Nl, replace = TRUE)
XFcoord = pF$x # feeding station coordinate X
YFcoord = pF$y # feeding station coordinate Y
XLcoord = rnorm(Nl, XFcoord[larvalLocations], .01) # larval station coordinate X
YLcoord = rnorm(Nl, YFcoord[larvalLocations], .01) # larval station coordinate Y

betaParam1 = 14; betaParam2 = 6
betaMean = betaParam1 / (betaParam1 + betaParam2)
betaVar = (betaParam1 * betaParam2) / ((betaParam1 + betaParam2) ^ 2 *(betaParam1 + betaParam2 + 1))
alpha = matrix(betaMean, nrow = Nl, ncol = 1)

bites = read.table('bitingDeBenedictis2003.txt', sep = ' ')[, 5]
bites = bites / sum(bites)
bitingRate = coef(mle(function(p1 = 1 / mean(bites)) -sum(log(dexp(bites, rate = p1)))))

hSize = 3.5
Hh = 2 + rpois(Nf, hSize - 2) # Poisson-like distribution of people across feeding stations, with two-person minimum
Nh = sum(Hh) # total number of hosts

L = matrix(0, Nf, Nl)
F = matrix(0, Nl, Nf)
L = mosquitoMovement(XFcoord, YFcoord, XLcoord, YLcoord, surv = survL, mixing = 'poor')
F = mosquitoMovement(XLcoord, YLcoord, XFcoord, YFcoord, surv = survF, mixing = 'poor')

H = makeHpoorHostMixing(Nf, Nh, Hh)
U = makeUrossMacdonald(H)

source('variablesNew.R')
source('calculateMetrics.R')
