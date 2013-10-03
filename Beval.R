# hierachical clustering eval

Beval.meth <- function(m.list, dismat, s.all,maxclust = 100)
{
  
  cl <- makeCluster(4)
  registerDoSNOW(cl)
  
  Beval <- function (clusts, dismat, s.all) {
    bBeta <- foreach(i = unique(clusts), .combine = "c")%do%{dismat[names(clusts[clusts %in% unique(clusts)[-c(1:which(unique(clusts) == i))]]),names(clusts[clusts == i])]}
    eval <- 1-(sum(bBeta)/s.all)
    return(eval)
  }
  
  all.eval <- foreach(m = m.list, .packages = "foreach", .combine = "rbind")%dopar%
{
  foreach(i = 2:maxclust, .combine = "c")%do%{gc();Beval(cutree(m,i),dismat,s.all)}
  
}
  stopCluster(cl)
 return(all.eval) 
}


# non-hierachical clustering eval

Beval.pam <- function( dis, dismat, s.all,maxclust = 100)
{
  
  cl <- makeCluster(4)
  registerDoSNOW(cl)
  
  Beval <- function (clusts, dismat, s.all) {
    bBeta <- foreach(i = unique(clusts)[-length(unique(clusts))], .combine = "c")%do%{dismat[names(clusts[clusts %in% unique(clusts)[-c(1:which(unique(clusts) == i))]]),names(clusts[clusts == i])]}
    eval <- 1-(sum(bBeta)/s.all)
    return(eval)
  }
  
  all.eval <- foreach(i = 2:(maxclust-1), .combine = "c", .packages = c("cluster","foreach"))%do%{gc();Beval(pam(dis,i,diss = T)$clustering,dismat,s.all)}
  

  stopCluster(cl)
  return(all.eval) 
}

