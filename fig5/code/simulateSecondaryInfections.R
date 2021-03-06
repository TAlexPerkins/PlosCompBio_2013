simulateSecondaryInfections = function(hh)
{
	source('variables.R')
	
	Sh[hh] = 0
	Ih = matrix(0, Nh, round(rho_mean))
	Ih[hh, 1] = 1
	IIh = as.numeric(rowSums(Ih) > 0)
	
	# loop over time
	tt = 1
	while((sum(Ih) + sum(Emf) + sum(Imf)) > 0){
		# mosquito movement from blood feeding habitats to larval habitats
		Sml = M_larvae[,xi] + rmultinomMatrix(Smf, L)
		for(ii in 1 : sigma)
			Eml[, ii] = matrix(rmultinomMatrix(Emf[, ii], L), nrow = Nl, ncol = 1)
		Iml = rmultinomMatrix(Imf, L)
		
		# egg laying and advancement through larval stages
		M_larvae = cbind(rpois(Nl, v * (Sml + rowSums(Eml) + Iml)), M_larvae[, c(1 : xi)[-xi]])
		M_larvae = matrix(rbinom(Nl * xi, M_larvae, (M_larvae + 1) ^ (matrix(rep(alpha, xi), Nl, xi) - 1)), Nl, xi)
		
		# mosquito movement from larval habitats to blood feeding habitats
		Smf = rmultinomMatrix(Sml, F)
		for(ii in 1 : (sigma - 1))
			Emf[, ii + 1] = matrix(rmultinomMatrix(Eml[, ii], F), nrow = Nf, ncol = 1)
		Imf = rmultinomMatrix(Iml, F) + rmultinomMatrix(Eml[, sigma], F)
		
		# transmission from infectious hosts to susceptible mosquitoes
		if(sum(Ih) > 1)
			Emf[, 1] = rbinom(Nf, sapply(1 : Nf, function(ii) sum(rmultinom(1, Smf[ii], U[, ii])[which(rowSums(Ih) > 0)])), c)
		if(sum(Ih) == 1)
			Emf[, 1] = rbinom(Nf, sapply(1 : Nf, function(ii) sum(rmultinom(1, Smf[ii], U[, ii])[which(rowSums(Ih) > 0)])), c)
		if(sum(Ih) == 0)
			Emf[, 1] = 0
		Smf = Smf - Emf[, 1]
		
		# infectious host recovery
		IhOld = Ih
		for(ii in seq(round(rho_mean) - 1, 1, -1))
			Ih[, ii + 1] = rbinom(Nh, Ih[, ii], 1 - rho_fail[ii])
		Rh = Rh + rowSums(IhOld) - rowSums(Ih[, 2 : round(rho_mean)])
		rm(IhOld)
		
		# host progression through the pathogen incubation period
		Ih[, 1] = Eh[, tau]
		for(ii in seq(tau - 1, 1, -1))
			Eh[, ii + 1] = Eh[, ii]
		
		# transmission from infectious mosquitoes to susceptible hosts
		secBites_h = rowSums(sapply(1 : Nf, function(ii) rmultinom(1, Imf[ii], U[, ii])))
		
		if(tt > 1){IIh = cbind(IIh, IIh[, tt] + rbinom(Nh, Sh, 1 - (1 - b) ^ secBites_h))}
		else{IIh = cbind(IIh, IIh + rbinom(Nh, Sh, 1 - (1 - b) ^ secBites_h))}
		tt = tt + 1
	}
	
	return(IIh)
} # end function simulatedSecondaryInfections definition
