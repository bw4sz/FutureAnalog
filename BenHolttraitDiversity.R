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
