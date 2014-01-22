if(!require(KernSmooth)){install.packages('KernSmooth'); library(KernSmooth)}


# function to numerically differentiate a vector y on x
numDiff = function(x, y){
	sapply(3 : (length(x) - 2), function(ii)
		(-y[ii + 2] + 8 * y[ii + 1] - 8 * y[ii - 1] + y[ii - 2]) / (12 * (x[2] - x[1])))
}
	
# calculate distance matrix for blood-feeding habitats
dist.F = matrix(0, Nf, Nf)
for(ii in 1 : Nf){
	for(jj in 1 : Nf){
		dist.F[ii, jj] = sqrt((XFcoord[ii] - XFcoord[jj]) ^ 2 + (YFcoord[ii] - YFcoord[jj]) ^ 2)
	}
}

# vectorize distances
distVec.F = as.vector(dist.F)



# compute pdf and cdf of kernel for Q matrix

# vectorize the matrix
Q_vec = as.vector(Q)

# make bins according to unique distances
bins = unique(as.vector(dist.F))

# empirical pdf and cdf
Q_pdf.emp = normalize(sapply(bins, function(dd) sum(Q_vec[which(distVec.F == dd)])))
Q_cdf.emp = sapply(bins, function(dd) sum(Q_vec[which(distVec.F <= dd)]))
Q_cdf.emp = Q_cdf.emp / max(Q_cdf.emp)

# smoothed cdf
Q_cdf.sth = ksmooth(bins[2 : length(bins)], Q_cdf.emp[2 : length(bins)], kernel = 'normal', bandwidth = 5 * dpill(bins, Q_cdf.emp))

# differentiate the smoothed cdf to obtain a smoothed pdf
Q_pdf.sth = list(x = Q_cdf.sth$x[3 : (length(Q_cdf.sth$x) - 2)],
								 y = normalize(numDiff(Q_cdf.sth$x, Q_cdf.sth$y)))

# redefine the smoothed cdf to be exactly consistent with the smoothed pdf
Q_cdf.sth$x = Q_pdf.sth$x
Q_cdf.sth$y = cumsum(Q_pdf.sth$y)



# compute pdf and cdf of kernel for P matrix

# vectorize the matrix
P_vec = as.vector(P)

# empirical pdf and cdf
P_pdf.emp = normalize(sapply(bins, function(dd) sum(P_vec[which(distVec.F == dd)])))
P_cdf.emp = sapply(bins, function(dd) sum(P_vec[which(distVec.F <= dd)]))
P_cdf.emp = P_cdf.emp / max(P_cdf.emp)

# smoothed cdf
P_cdf.sth = ksmooth(bins[2 : length(bins)], P_cdf.emp[2 : length(bins)], kernel = 'normal', bandwidth = 5 * dpill(bins, P_cdf.emp))

# differentiate the smoothed cdf to obtain a smoothed pdf
P_pdf.sth = list(x = P_cdf.sth$x[3 : (length(P_cdf.sth$x) - 2)],
								 y = normalize(numDiff(P_cdf.sth$x, P_cdf.sth$y)))

# redefine the smoothed cdf to be exactly consistent with the smoothed pdf
P_cdf.sth$x = P_pdf.sth$x
P_cdf.sth$y = cumsum(P_pdf.sth$y)



# compute pdf and cdf of kernel for Z matrix

# vectorize the matrix
Z_vec = as.vector(Z)

# empirical pdf and cdf
Z_pdf.emp = normalize(sapply(bins, function(dd) sum(Z_vec[which(distVec.F == dd)])))
Z_cdf.emp = sapply(bins, function(dd) sum(Z_vec[which(distVec.F <= dd)]))
Z_cdf.emp = Z_cdf.emp / max(Z_cdf.emp)

# smoothed cdf
Z_cdf.sth = ksmooth(bins[2 : length(bins)], Z_cdf.emp[2 : length(bins)], kernel = 'normal', bandwidth = 5 * dpill(bins, Z_cdf.emp))

# differentiate the smoothed cdf to obtain a smoothed pdf
Z_pdf.sth = list(x = Z_cdf.sth$x[3 : (length(Z_cdf.sth$x) - 2)],
								 y = normalize(numDiff(Z_cdf.sth$x, Z_cdf.sth$y)))

# redefine the smoothed cdf to be exactly consistent with the smoothed pdf
Z_cdf.sth$x = Z_pdf.sth$x
Z_cdf.sth$y = cumsum(Z_pdf.sth$y)



# compute pdf and cdf of kernel for S matrix

# vectorize the matrix
S_vec = as.vector(S)

# empirical pdf and cdf
S_pdf.emp = normalize(sapply(bins, function(dd) sum(S_vec[which(distVec.F == dd)])))
S_cdf.emp = sapply(bins, function(dd) sum(S_vec[which(distVec.F <= dd)]))
S_cdf.emp = S_cdf.emp / max(S_cdf.emp)

# smoothed cdf
S_cdf.sth = ksmooth(bins[2 : length(bins)], S_cdf.emp[2 : length(bins)], kernel = 'normal', bandwidth = 5 * dpill(bins, S_cdf.emp))

# differentiate the smoothed cdf to obtain a smoothed pdf
S_pdf.sth = list(x = S_cdf.sth$x[3 : (length(S_cdf.sth$x) - 2)],
								 y = normalize(numDiff(S_cdf.sth$x, S_cdf.sth$y)))

# redefine the smoothed cdf to be exactly consistent with the smoothed pdf
S_cdf.sth$x = S_pdf.sth$x
S_cdf.sth$y = cumsum(S_pdf.sth$y)
