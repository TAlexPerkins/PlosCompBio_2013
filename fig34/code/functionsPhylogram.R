if(!require(ape)){install.packages('ape'); library(ape)}
if(!require(cluster)){install.packages('cluster'); library(cluster)}
if(!require(RColorBrewer)){install.packages('RColorBrewer'); library(RColorBrewer)}
if(!require(colorspace)){install.packages('colorspace'); library(colorspace)}


evenness = function(x, alpha = 1)
{
  x = x / sum(x)
  
  if(alpha == 1)
  {
    return(exp(-sum(x[which(x > 0)] * log(x[which(x > 0)]))) / length(x))
  } # end if
  
  if(alpha != 1)
  {
    return(sum(x ^ alpha) ^ (1 / (1 - alpha)) / length(x))
  } # end if
} # end function HillEvenness definition



getDescendants<-function(tree,node,curr=NULL){
  # by Liam Revell
  # http://phytools.blogspot.com/2012/01/function-to-get-descendant-node-numbers.html
  if(is.null(curr)) curr<-vector()
  daughters<-tree$edge[which(tree$edge[,1]==node),2]
  curr<-c(curr,daughters)
  w<-which(daughters>=length(tree$tip))
  if(length(w)>0) for(i in 1:length(w)) 
    curr<-getDescendants(tree,daughters[w[i]],curr)
  return(curr)
} # end function getDescendants definition



makePhylo = function(m)
{
  # Perform agglomerative nesting algorithm across rows of matrix and convert to phylo class.
  return(as.phylo(as.hclust(agnes(m, metric = 'euclidean', method = 'average'))))
} # end function makePhylo definition



phylo.het = function(p, m, matrixIndex = 1 : length(p$tip.label))
{
  hetList = list()
	
	for(ii in 1 : dim(p$edge)[1])
	{
		desc = getDescendants(p, ii)[which(getDescendants(p, ii) < length(p$tip))]
		if(length(desc)){
      if(length(desc) == 1){hetList[[ii]] = .001}
      else{hetList[[ii]] = sd(rowSums(m[desc, desc])) / mean(rowSums(m[desc, desc]))}
		}
		rm(desc)
	} # end for
	
	return(hetList)
} # end function phylo.het definition



phylo.mix = function(p, m, matrixIndex = 1 : length(p$tip.label))
{
	mixList = list()
	
	for(ii in 1 : dim(p$edge)[1])
	{
		desc = getDescendants(p, ii)[which(getDescendants(p, ii) < length(p$tip))]
		if(length(desc)){
			mixList[[ii]] = evenness(as.vector(m[desc, desc]))
		}
		rm(desc)
	} # end for
	
	return(mixList)
} # end function phylo.mix definition



phylo.magnitude = function(p, m)
{
  magnitudeList = list()
  
  for(ii in 1 : dim(p$edge)[1])
  {
    desc = getDescendants(p, ii)[which(getDescendants(p, ii) < length(p$tip))]
    magnitudeList[[ii]] = sapply(desc, function(ii) sum(sapply(1 : length(desc), function(jj) sum(m[ii, jj]))))
    rm(desc)
  } # end for
  
  return(magnitudeList)
} # end function phylo.magnitude definiton



plotPhyloHet = function(p, m, matrixIndex = 1 : length(p$tip.label), direction = 'downwards')
{
  # Colors
  branch.colors = c(rev(brewer.pal(8, 'Blues'))[1 : 5], brewer.pal(8, 'Reds')[4 : 8])
  
  # Calculate heterogeneity of matrix subset of the descendants of each non-tip node.
  ph = phylo.het(p, m, matrixIndex)
  ph = unlist(sapply(ph, function(pp) pp / max(range(ph))))
  
  # Painted phylograms
  plot.phylo(p, show.tip.label = FALSE,
    edge.color = (c(rep('grey', length(ph) + 2), branch.colors[ceiling(ph * 10)]))[rank(p$edge[, 2])],
    edge.width = 2, type = 'phylogram', direction = direction, no.margin = TRUE,
    x.lim = c(-2, 52), y.lim = c(0, .6))
  par(mar = rep(2, 4))
} # end function plotPhyloHet definition



plotPhyloMix = function(p, m, matrixIndex = 1 : length(p$tip.label), direction = 'downwards')
{
  # Colors
  branch.colors = c(rev(brewer.pal(8, 'Blues'))[1 : 5], brewer.pal(8, 'Reds')[4 : 8])
  
  # Calculate evenness of matrix subset of the descendants of each non-tip node.
  pm = phylo.mix(p, m, matrixIndex)
  pm = unlist(sapply(pm, function(pp) pp / max(range(pm))))
  
  # Painted phylograms
  plot.phylo(p, show.tip.label = FALSE,
    edge.color = (c(rep('grey', length(pm) + 2), branch.colors[ceiling(pm * 10)]))[rank(p$edge[, 2])],
    edge.width = 2, type = 'phylogram', direction = direction, no.margin = TRUE,
    x.lim = c(-2, 52), y.lim = c(0, .6))
  par(mar = rep(2, 4))
} # end function plotPhyloMix definition



plotHetByDist = function(p, m)
{
  point.colors = c(rev(brewer.pal(8, 'Blues'))[1 : 5], brewer.pal(8, 'Reds')[4 : 8])
  
  het = unlist(phylo.het(p, m))
  het = het / max(het)
  
  btp = branching.times(p)[2 : p$Nnode]
  plot(het, btp,
    col = point.colors[ceiling(het * 10)],
    xlim = c(0, max(het)),
    ylim = c(0, max(branching.times(p))),
    pch = 15, xlab = '', ylab = '', xaxt = 'n', yaxt = 'n') 
} # end function plotHetByDist definition



plotMixByDist = function(p, m)
{
  point.colors = c(rev(brewer.pal(8, 'Blues'))[1 : 5], brewer.pal(8, 'Reds')[4 : 8])
  
  mix = unlist(phylo.mix(p, m))
  mix = mix / max(mix)
  
  btp = branching.times(p)[2 : p$Nnode]
  plot(mix, btp,
    col = point.colors[ceiling(mix * 10)],
    xlim = c(0, max(mix)),
    ylim = c(0, max(branching.times(p))),
    pch = 15, xlab = '', ylab = '', xaxt = 'n', yaxt = 'n') 
} # end function plotMixByDist definition



plotHetBySize = function(p, m)
{
  point.colors = c(rev(brewer.pal(8, 'Blues'))[1 : 5], brewer.pal(8, 'Reds')[4 : 8])
  
  het = unlist(phylo.het(p, m))
  het = het / max(het)
  
  numberOfDescendants = sapply((p$Nnode + 2) : (2 * p$Nnode), function(ii) length(getDescendants(p, ii)))
  plot(het, numberOfDescendants,
    col = point.colors[ceiling(het * 10)],
    xlim = c(0, max(het)),
    ylim = c(0, max(numberOfDescendants)),
    pch = 15, xlab = '', ylab = '', xaxt = 'n', yaxt = 'n') 
} # end function plotHetBySize definition



plotMixBySize = function(p, m)
{
  point.colors = c(rev(brewer.pal(8, 'Blues'))[1 : 5], brewer.pal(8, 'Reds')[4 : 8])
  
  mix = unlist(phylo.mix(p, m))
  mix = mix / max(mix)
  
  numberOfDescendants = sapply((p$Nnode + 2) : (2 * p$Nnode), function(ii) length(getDescendants(p, ii)))
  plot(mix, numberOfDescendants,
    col = point.colors[ceiling(mix * 10)],
    xlim = c(0, max(mix)),
    ylim = c(0, max(numberOfDescendants)),
    pch = 15, xlab = '', ylab = '', xaxt = 'n', yaxt = 'n') 
} # end function plotMixBySize definition
