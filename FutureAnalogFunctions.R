fnBetaDiv <- function(mod){
  # get niche for the GCM 
  niche <- niche.crops[grep(mod,niche.crops,value=FALSE)]
  
  #Create siteXspp table from input rasters, function is from
  #AlphaMappingFunctions.R, sourced at the top.
  siteXspps <- tableFromRaster(niche, threshold = 0.05)
  
  #Remove NAs from siteXspps 
  fails <- na.test(siteXspps)
  siteXspps <- siteXspps[,!colnames(siteXspps) %in% fails]
  
  # PART II: Calculate beta diversity (between time) -----------------------------
  # compare current with each future scenario
  
  # Step 1) TAXONOMIC BETA DIVERSITY ---------------------------------------------
  beta.time.taxa <- analogue::distance(current, siteXspps, "bray")
  
  # Step 2) PHYLO BETA DIVERSITY -------------------------------------------------
  # For phylobeta, there needs to be more than 2 species for a rooted tree
  phylo.dat <- siteXspps[,colnames(siteXspps) %in% trx$tip.label]
  phylo.dat <- phylo.dat[!rowSums(phylo.dat)<=2,]   
  
  strt <- Sys.time()
  beta.time.phylo <- matpsim.pairwise(phyl = trx, 
                                      com.x = current.phylo, 
                                      com.y = phylo.dat)
  Sys.time()-strt
  
  # Step 3) FUNC BETA DIVERSITY ---------------------------------------------------
  func.dat <- siteXspps[,colnames(siteXspps) %in% colnames(fco)]
  func.dat <- func.dat[!apply(func.dat,1,sum)<=2,]  
  
  sp.list_current <- lapply(rownames(current.func), function(k){
    g <- current.func[k,]
    names(g[which(g==1)])
  })
  
  names(sp.list_current) <- rownames(current.func)
  
  sp.list_future <- lapply(func.dat, function(k){
    g <- func.dat[k,]
    names(g[which(g==1)])
  })
  
  names(sp.list_future) <- rownames(func.dat)
  
  #Get distances from the cophenetic matrix?
  dists <- as.matrix(fco)
  rownames(dists) <- rownames(fco)
  colnames(dists) <- rownames(fco)
  
  cl <- makeCluster(clust) # make another cluster
  registerDoSNOW(cl)
  
  beta.time.func <- foreach(fu=rownames(func.dat)) %dopar%{
    sapply(rownames(current.func), function(cur){
      MNND_fc(fu, cur, sp.list_current, sp.list_future, dists)
    })}
  
  
  stopCluster(cl)

 res <- list(beta.time.taxa, beta.time.phylo, beta.time.func)
 names(res) <- c("beta.time.taxa", "beta.time.phylo", "beta.time.func")
 return(res)
 save(res, file = paste0(out_path, "/beta_diversity_", mod, ".rda"))
  return(beta.time.taxa)
}

fnCurrent2Future <- function(betadiv, arb.thresh) {
  n.analogs <- sapply(rownames(betadiv), function(x){
    sum(betadiv[rownames(betadiv) %in% x,] <= arb.thresh) 
    #counts the number of assemblages with beta div values less than arb.thresh
  })
  
  current_to_future.analog <- data.frame(rownames(betadiv), n.analogs)
  colnames(current_to_future.analog) <- c("cell.number", "numberofanalogs")  
  
  c_f <- cellVis(cells=current_to_future.analog$cell.number, 
                 value=current_to_future.analog$numberofanalogs)
  hist(current_to_future.analog$numberofanalogs)
  return(c_f)
}

fnFuture2Current <- function(betadiv, arb.thresh){
  n.analogs <- sapply(colnames(betadiv), function(x){
    sum(betadiv[,colnames(betadiv) %in% x] <= arb.thresh)
  })
  
  future_to_current.analog <- data.frame(colnames(betadiv), n.analogs)
  colnames(future_to_current.analog) <- c("cell.number", "numberofanalogs")  
  
  f_c <- cellVis(cell=future_to_current.analog$cell.number, 
                 value=future_to_current.analog$numberofanalogs)
  hist(future_to_current.analog$numberofanalogs)
  return(f_c)
}

na.test <-  function (x) {
  w <- apply(x, 2, function(x)all(is.na(x)))
  if (any(w)) {
    fails <- names(which(w))
    print(paste("All NA in columns", paste(names(which(w)), collapse=", ")))
    return(fails)
  }
}
