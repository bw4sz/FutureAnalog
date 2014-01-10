MNND <- function(A,B,sp.list,dists)
{

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
