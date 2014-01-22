source('functionsModel.R')
if(!require(stats4)){install.packages('stats4'); library(stats4)}


T = 400 # number of time steps to run
Nf = 50 # number of blood-feeding habitats
Nl = 50 # number of larval habitats
v = 10 # number of mosquito eggs per capita
xi = 4 # larval stages for maturation

survL = 0.9
survF = 0.9

sigma = 3 # mosquito incubation period
tau = 3 # host incubation period
rho_max = 40 # rho is host recovery and these lines define its distribution
rho_pmf = normalize(dnbinom(0 : (rho_max - 1), 3, .2))
rho_mean = sum((1 : rho_max) * rho_pmf)
rho_fail = rho_pmf / (1 - c(0, pnbinom(0 : (rho_max - 2), 3, .2)))
c = 0.15 # host-to-mosquito transmission efficiency
b = 0.55 # mosquito-to-host transmission efficiency

fac = matrix(seq(0.9, 1, length.out = xi), byrow = T, nrow = Nl, ncol = xi)
betaParam1 = 14; betaParam2 = 6
betaMean = betaParam1 / (betaParam1 + betaParam2)
betaVar = (betaParam1 * betaParam2) / ((betaParam1 + betaParam2) ^ 2 *(betaParam1 + betaParam2 + 1))
alphaEven = matrix(betaMean, nrow = Nl, ncol = xi) * fac
alphaUneven = matrix(rbeta(Nl, betaParam1, betaParam2), nrow = Nl, ncol = xi) * fac
alpha = alphaEven
larvalLocations = sample(1 : Nf, Nl, replace = TRUE)

bites = read.table('../data/bitingDeBenedictis2003.txt', sep = ' ')[, 5]
bites = bites / sum(bites)
bitingRate = coef(mle(function(p1 = 1 / mean(bites)) -sum(log(dexp(bites, rate = p1)))))

hSize = 5.5
Hh = 2 + rpois(Nf, hSize - 2) # Poisson-like distribution of hosts across feeding stations, with two-host minimum
Nh = sum(Hh) # total number of hosts



save(list = ls(), file = '../data/parametersCommon.RData')