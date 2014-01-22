# load function for simulating an epidemic
source('simulateEpidemic.R')



# define function for estimating force of infection given simulated data
estimateFOI = function()
{
  # allocate vectors to track new and cumulative infections
  newInfections = numeric()
  cumInfections = numeric()

  # simulate Nreps epidemics initiated by the individual with the most average individual R0
  hh = which.min(abs(rowSums(R_nonlinear) - mean(rowSums(R_nonlinear))))
  Nreps = 100
  for(rr in 1 : Nreps)
  {
    incidence = simulateEpidemic(hh, T = 50)
    newInfections = c(newInfections, incidence / Nh)
    cumInfections = c(cumInfections, cumsum(incidence) / Nh)
  } # end for
    
  # find the slope of a linear regression of new infections against cumulative infections
  slope = sapply(sort(unique(cumInfections))[min(which(sort(unique(cumInfections)) > .01)) : length(unique(cumInfections))],
  							 function(ii) lm(newInfections[which(cumInfections < ii)] ~ cumInfections[which(cumInfections < ii)])$coefficients[2])
  r.squared = sapply(sort(unique(cumInfections))[min(which(sort(unique(cumInfections)) > .01)) : length(unique(cumInfections))],
  									 function(ii) summary(lm(newInfections[which(cumInfections < ii)] ~ cumInfections[which(cumInfections < ii)]))$r.squared)
  
  return(slope[which.max(r.squared)])
} # end estimateFOI function definition



# estimate force of infection and calculate R0 under different formulations for doing so
# low transmission setting with high mosquito mortality
source('parameters_low.R')
FOI = estimateFOI()
R0_FOI1.low = 1 + FOI / rho_mean # Marques et al. 1994
R0_FOI2.low = (1 + FOI / rho_mean) * (1 + FOI / -log(survL * survF)) # Massad et al. 2001
R0_FOI3.low = (1 + FOI / rho_mean) * (1 + FOI / -log(survL * survF)) * exp(FOI * (sigma + tau)) # Favier et al. 2006
save(list = ls(), file = '../data/low.RData')

# high transmission setting with low mosquito mortality
source('parameters_high.R')
FOI = estimateFOI()
R0_FOI1.high = 1 + FOI / rho_mean # Marques et al. 1994
R0_FOI2.high = (1 + FOI / rho_mean) * (1 + FOI / -log(survL * survF)) # Massad et al. 2001
R0_FOI3.high = (1 + FOI / rho_mean) * (1 + FOI / -log(survL * survF)) * exp(FOI * (sigma + tau)) # Favier et al. 2006
save(list = ls(), file = '../data/high.RData')
rm(list = ls())



# compile estimates of R0 from above and save in ../output/table.RData
tab = matrix(0, 6, 6)

load('../data/low.RData')
tab[1, 2] = R0_rossMacdonald
tab[2, 2] = R0_rossMacdonald * (1 + (sd(rowSums(R_linear)) / mean(rowSums(R_linear))) ^ 2)
tab[3, 2] = R0_rossMacdonald * (1 + (sd(rowSums(R_nonlinear)) / mean(rowSums(R_nonlinear))) ^ 2)
tab[4, 2] = R0_linear
tab[5, 2] = R0_nonlinear
tab[6, 2] = R0_FOI3.low

tab[, 1] = tab[, 2] / tab[5, 2]

tab[, 3] = 1 - 1 / tab[, 2]


load('../data/high.RData')
tab[1, 5] = R0_rossMacdonald
tab[2, 5] = R0_rossMacdonald * (1 + (sd(rowSums(R_linear)) / mean(rowSums(R_linear))) ^ 2)
tab[3, 5] = R0_rossMacdonald * (1 + (sd(rowSums(R_nonlinear)) / mean(rowSums(R_nonlinear))) ^ 2)
tab[4, 5] = R0_linear
tab[5, 5] = R0_nonlinear
tab[6, 5] = R0_FOI3.high

tab[, 4] = tab[, 5] / tab[5, 5]

tab[, 6] = 1 - 1 / tab[, 5]


save(tab, file = '../output/table.RData')