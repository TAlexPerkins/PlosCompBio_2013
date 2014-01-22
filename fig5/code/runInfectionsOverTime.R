load('../data/parametersCommon.RData')
if(!require(snowfall)){install.packages('snowfall'); library(snowfall)}
source('functions.R')
source('simulateEpidemic.R')
source('simulateSecondaryInfections.R')
Nreps = 100



alpha = alphaEven
L = matrix(0, Nf, Nl)
F = matrix(0, Nl, Nf)
L = mosquitoMovement(XFcoord, YFcoord, XLcoord, YLcoord, surv = survL, mixing = 'well')
F = mosquitoMovement(XLcoord, YLcoord, XFcoord, YFcoord, surv = survF, mixing = 'well')
H = HrossMacdonald
U = UrossMacdonaldRossMacdonald
Nh = dim(H)[1]
M_larvae = matrix(0, nrow = Nl, ncol = xi)
M_larvaeNew = matrix(0, nrow = Nl, ncol = xi)
M_larvaeNew = round(estimateLarvalEquilib())
SmlNew = round(estimateMosquitoEquilibLarval())
SmfNew = round(estimateMosquitoEquilibFeeding())
source('variables.R')
source('calculateMetrics.R')

secondaryInfections_000 = matrix(0, Nh, round(rho_mean))
while(colSums(secondaryInfections_000)[dim(secondaryInfections_000)[2]] != round(mean(rowSums(R_nonlinear))))
{
  secondaryInfections_000 = simulateSecondaryInfections(1)
} # end while

data_000 = vector('list', Nreps)

sfInit(parallel = TRUE, cpus = 16)
stopifnot(sfCpus() >= 4)
stopifnot(sfParallel() == TRUE)

sfExportAll()
data_000 = sfSapply(1 : Nreps, function(ii) simulateEpidemic(1))
save(list = ls(), file = '../data/000.RData')



alpha = alphaEven
L = matrix(0, Nf, Nl)
F = matrix(0, Nl, Nf)
L = mosquitoMovement(XFcoord, YFcoord, XLcoord, YLcoord, surv = survL, mixing = 'poor')
F = mosquitoMovement(XLcoord, YLcoord, XFcoord, YFcoord, surv = survF, mixing = 'poor')
H = HpoorHostMixing
U = UpoorHostMixingRossMacdonald
Nh = dim(H)[1]
source('variables.R')
source('calculateMetrics.R')

humanMaxR = which.max(rowSums(R_nonlinear))
humanMeanR = which.min(abs(rowSums(R_nonlinear) - mean(rowSums(R_nonlinear))))

secondaryInfections_110 = matrix(0, Nh, round(rho_mean))
while((
	colSums(secondaryInfections_110)[dim(secondaryInfections_110)[2]] !=
	round(mean(rowSums(R_nonlinear)))) | (max(secondaryInfections_110[, dim(secondaryInfections_110)[2]]) < 3))
{
  secondaryInfections_110 = simulateSecondaryInfections(humanMeanR)
} # end while

data_110 = vector('list', Nreps)

sfInit(parallel = TRUE, cpus = 16)
stopifnot(sfCpus() >= 4)
stopifnot(sfParallel() == TRUE)

sfExportAll()
data_110 = sfSapply(1 : Nreps, function(ii) simulateEpidemic(humanMeanR))
save(list = ls(), file = '../data/110.RData')
