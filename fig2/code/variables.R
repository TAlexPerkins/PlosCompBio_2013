# initialize variables

M_larvae = matrix(0, nrow = Nl, ncol = xi)
MM = array(0, dim = c(Nl, xi, T))
M_larvae = round(estimateLarvalEquilib())

Sml = round(estimateMosquitoEquilibLarval())
Smf = round(estimateMosquitoEquilibFeeding())
Eml = matrix(0, nrow = Nl, ncol = sigma + 1)
Emf = matrix(0, nrow = Nf, ncol = sigma + 1)
Iml = rep(0, Nl)
Imf = rep(0, Nf)

Sh = rep(0, Nh)
Eh = matrix(0, nrow = Nf, ncol = tau + 1)
Ih = rep(0, Nh)
Rh = rep(0, Nh)

SSml = matrix(0, nrow = Nl, ncol = T)
SSmf = matrix(0, nrow = Nf, ncol = T)
EEml = array(0, dim = c(Nl, sigma + 1, T))
EEmf = array(0, dim = c(Nf, sigma + 1, T))
IIml = matrix(0, nrow = Nl, ncol = T)
IImf = matrix(0, nrow = Nf, ncol = T)

SSh = matrix(0, nrow = Nh, ncol = T)
EEh = array(0, dim = c(Nf, tau + 1, T))
IIh = matrix(0, nrow = Nh, ncol = T)
RRh = matrix(0, nrow = Nh, ncol = T)