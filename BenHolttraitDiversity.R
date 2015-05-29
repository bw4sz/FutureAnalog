MNND <- function(A,B,sp.list,dists){
# Based on Holt et al. 2013. An Update of Wallace’s Zoogeographic Regions of the World. Science 339: 74-78.
# Supplemental Material : http://www.jeanphilippelessard.com/Jean-Philippe_Lessard/Publications_files/Holt.SM.pdf
# Inputs a row iterator (A), and a row index (B), the species list (which could be traits or phylo) and a distance matrix
# Outputs mean nearest neighbor taxon distance

Asp     <- sp.list[[A]]
Bsp     <- sp.list[[B]]
compmat <- dists[Asp,Bsp]
Ann     <- apply(as.matrix(compmat),1,min)
Bnn     <- apply(as.matrix(compmat),2,min)
Dnn     <- mean(c(Ann, Bnn))
#turn    <- min(c(mean(Ann),mean(Bnn)))
#nest    <- Dnn - turn
#res <- c(Dnn,turn,nest)
res<-c(Dnn)
#names(res) <- c("MNND","MNNDturn","MNNDnest")
names(res) <- c("MNND")

res
}

func.dist.mat <- function(sppXsite, sp.list, dists){
  # get all combinations of site pairs
  site.pairs <- data.frame(t(combn(rownames(sppXsite), 2)), stringsAsFactors = FALSE)
  
  # calculate teh mean nearest neighbour taxon distance
  site.pairs$MNND <- apply(site.pairs, 1, function(x) {
    MNND(x[1], x[2], sp.list = sp.list, dists = dists)
  })
  
  # the first and last sites will not have a row and column respectively unless
  # their measurements to themselves is given
  first <- rownames(sppXsite)[1]
  last <- rownames(sppXsite)[nrow(sppXsite)]
  first.last <- data.frame(X1 = c(first, last), X2 = c(first, last), MNND = c(0, 0))
  site.pairs <- rbind(site.pairs, first.last)
  
  # site numbers need to be numeric for the ordering to work
  site.pairs$X1 <- as.numeric(site.pairs$X1)
  site.pairs$X2 <- as.numeric(site.pairs$X2)
  
  # create matrix, set diag to zero and copy the bottom to the top triangle
  func.dist.mat <- spread(site.pairs, "X1", "MNND")
  rownames(func.dist.mat) <- func.dist.mat[,1]
  func.dist.mat <- func.dist.mat[-1]
  diag(func.dist.mat) <- 0
  func.dist.mat[upper.tri(func.dist.mat)] <- 
    t(func.dist.mat)[upper.tri(func.dist.mat)]
  return(func.dist.mat)
}
