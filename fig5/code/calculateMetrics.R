if(!require(Biodem)){install.packages('Biodem'); library(Biodem)}
if(!require(popbio)){install.packages('popbio'); library(popbio)}


lambda_equilib = estimateLarvalEquilib()[, xi]

nSum = 100
Q = matrix(data = 0, nrow = Nf, ncol = Nf)
for (qq in sigma : nSum)
{
  Q = Q + mtx.exp((L %*% F), qq)
} # summing over time since infection in mosquitoes

nSum = 100
M = rep(0, Nf)
for (bb in 0 : (sigma + nSum))
{
  M = M + F %*% mtx.exp((L %*% F), bb)
} # summing over adult mosquito ages
M = lambda_equilib %*% M

B = diag(as.vector(M), Nf, Nf) %*% t(U)

W = t(sapply(1 : Nh, function(kk) normalize(B[, kk])))

P = t(U) %*% W

Z = (1 - exp(-b * Q %*% t(U))) %*% (c * t(B) * rho_mean)

V = t(B) %*% Q %*% t(U)

R_linear = b * c * t(B) %*% Q %*% t(U) * rho_mean
diag(R_linear) = 0

R_nonlinear = 1 - exp(-b * c * t(B) %*% Q %*% t(U) * rho_mean)
diag(R_nonlinear) = 0

R0_rossMacdonald = b * c * sum(B) / Nh * rho_mean * sum(Q) / Nf

R0_linear = lambda(R_linear)

R0_nonlinear = lambda(R_nonlinear)