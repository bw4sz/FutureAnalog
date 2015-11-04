
#Write a function that gets the siteXspp matrix from an input list of niche_model outputs
tableFromRaster<-function(fil_list,threshold, clust = 7){
  
  cl <- makeCluster(clust) # create parellel clusters
  registerDoSNOW(cl)
  
  trial<-foreach(mod=1:length(fil_list),.combine=cbind, .packages = c("raster", "dplyr")) %dopar% {
    #read in the niche models from file
    niche.m<-lapply(fil_list[[mod]],function(x) raster(x))
    
    niche_ens<-niche.m[[1]][[1]]
    
    #get species name
    #This portion is a bit sensitive to file structure, if error calls, check here first!
    # Fixed so not so sensitive to file structure
    us.pos <- gregexpr(pattern = "_", niche_ens@file@name) # pos of underscores
    sp.start <- us.pos[[1]][[length(us.pos[[1]]) - 1]] # pos of penultimate us
    sp.end <- us.pos[[1]][[length(us.pos[[1]])]] # pos of last us
    sp.n <- substr(niche_ens@file@name, sp.start + 1, sp.end - 1)
    
    #Lets go get the presence data on hummingbird distributions
    Niche_locals <- read.csv("InputData/MASTER_POINTLOCALITYarcmap_review.csv")

    #Just take the columns you want. 
    PAdat <- select(Niche_locals, RECORD_ID, SPECIES, COUNTRY, LOCALITY, LATDECDEG,
                  LONGDECDEG, Decision, SpatialCheck, MapDecision)
    
    gooddata <- c("ok", "Ok", "OK", "Y") #note that blanks ("") and REJECT data are excluded
    loc_clean <- filter(PAdat, SpatialCheck=="Y", MapDecision %in% gooddata)
    
    sp.loc <- filter(loc_clean, SPECIES %in%  gsub("\\."," ",sp.n))
    sp.loc<-SpatialPointsDataFrame(sp.loc[,c("LONGDECDEG","LATDECDEG")],sp.loc)
    
    #draw suitability for occurence points
    site_suit<-raster::extract(niche_ens,sp.loc)
    suit_cut<-stats::quantile(na.omit(site_suit),threshold)
    
    #plot to file
    #print(qplot(site_suit) + geom_vline(aes(xintercept=suit_cut),linetype="dashed",col="Red"))
    
    paste("suitability threshold:",suit_cut,"")
    #Predicted Presence absence Column
    A_list<-values(niche_ens > suit_cut)*1
  }
  
  stopCluster(cl)
  
  #use a regular expression to extract names
  species <- str_match(fil_list,pattern=paste(cell,"(\\w+.\\w+)/proj_",sep="/"))[,2]
    
  #Name the rows and columns. 
  colnames(trial)<-species
  rownames(trial)<- as.factor(1:nrow(trial))
  
  #How many species in a each row?
  sum.row<-apply(trial,1,sum,na.rm=TRUE)
  
  #Get rid of cells with no species
  siteXspp<-trial[!is.na(sum.row),]
  siteXspp.raster<-siteXspp[apply(siteXspp,1,sum,na.rm=TRUE)>0,]
  #return the raster for each siteXspp
  return(siteXspp.raster)
}

#write a function that chunks vectors into pieces
chunk<-function(dat,max){
  j <- seq_along(dat)
  ind<-split(dat, ceiling(j/max))
  return(ind)
}

AlphaPhylo<-function(siteXspp.raster){
  #remove any species not in the phylogeny
  siteXspp.raster.Phylo<-siteXspp.raster[ , colnames(siteXspp.raster) %in% colnames(co)]
  
  #Get sum rows
  rowS.Phylo<-apply(siteXspp.raster.Phylo,1,sum,na.rm=TRUE)
  
  #Trim siteXspp
  siteXspp.Phylo<-siteXspp.raster.Phylo[rowS.Phylo > 1,]
  
  #remove the row sum
  rm(rowS.Phylo)
  
  #remove untrimmed siteXspp data
  rm(siteXspp.raster.Phylo)
  
  #remove any species with na
  toremove<-apply(siteXspp.Phylo,2,sum)
  
  siteXspp.Phylo<-siteXspp.Phylo[,colnames(siteXspp.Phylo) %in% names(which(!is.na(toremove)))]
  
  #mean phylogenetic branch length at each community
  a<-mpd(siteXspp.Phylo,cophenetic(trx))
  
  #Add in cell names and write to file
  MPDcell<-data.frame(rownames(siteXspp.Phylo),a)
  colnames(MPDcell)<-c("Cell","MPD")
  return(MPDcell)
}

# Needs to be changed with new understanding of traits (no dendrograms!)
AlphaFunc<-function(siteXspp.raster){
  
  #Only keep communities with species with functional information and richness > 1
  siteXspp.raster.Func<-siteXspp.raster[,colnames(siteXspp.raster) %in% colnames(fco)]
  
  rowS.Func<-apply(siteXspp.raster.Func,1,sum,na.rm=TRUE)
  siteXspp.Func<-siteXspp.raster.Func[rowS.Func > 1,]
  
  #remove any species with na
  toremove<-apply(siteXspp.Func,2,sum)
  
  siteXspp.Func<-siteXspp.Func[ ,colnames(siteXspp.Func) %in% names(which(!is.na(toremove)))]
  
  #find mean func neighbor at each assemblage
  a<-mpd(siteXspp.Func,fco)
    
  #Put the cell numbers so it can be identified as each community
  MFDCell<-data.frame(rownames(siteXspp.Func),a)
  colnames(MFDCell)<-c("Cell","MFD")
  
  return(MFDCell)
}

#Start with the three trait data, means we need more than 2 species 

#Create FD compute function
# Delete
AlphaFunc.FD<-function(siteXspp.raster,traits){
    
  #Need more species than traits
  siteXspp.FD<-siteXspp.raster[,colnames(siteXspp.raster) %in% gsub(" ","\\.",traits$Species)]
  rowS.Func<-apply(siteXspp.FD,1,sum)
  siteXspp.FD<-siteXspp.raster[rowS.Func > ncol(traits),]
  
  #Run FD metric on each row in the raster
  cl<-makeCluster(8,"SOCK")
  registerDoSNOW(cl)
  
  trait_FD<-foreach(h=1:nrow(siteXspp.FD),.packages="FD",.errorhandling="pass",.combine=rbind) %dopar% {
    g<-siteXspp.FD[h,] 
    
    #Which species are present in the community
    f<-names(g[which(g==1)])
    
    #For each community apply the functional metrics
    trait.site<-traits[traits$Species %in% gsub("\\."," ",f),]
    
    #make species rownames
    rownames(trait.site)<-trait.site[,1]
    trait.site<-trait.site[,-1]
    
    trait.output<-dbFD(trait.site,calc.CWM=FALSE,messages=TRUE)
    trait.out<-as.data.frame(trait.output)
    
    return(trait.out)}
  
  stopCluster(cl)
  
  rownames(trait_FD)<-rownames(siteXspp.FD)
  return(trait_FD)}

cellVis<-function(cells,value){
  #create a empty raster the exact dimensions of the model, this needs to match extent and grain first
  alpha_structure<-raster(blank)
  values(alpha_structure)<-NA
  
  alpha_structure[as.numeric(as.character(cells))]<-value
  names(alpha_structure)<-NULL
  plot(alpha_structure)
  
  return(alpha_structure)
}

#wrapper function to be called for all types of climate layers and time periods
cellVisuals<-function(inp.name){
  
  #Taxonomic richness
  richnessvalues<-apply(siteXspps[names(siteXspps) %in% inp.name][[1]],1,sum,na.rm=TRUE)
  cellP<-rownames(siteXspps[names(siteXspps) %in% inp.name][[1]])
  richness<-cellVis(cells=cellP,value=richnessvalues)
  
  #Phylogenetic richness
  prichness<-MPDs[names(MPDs) %in% inp.name][[1]]$MPD
  pcellP<-MPDs[names(MPDs) %in% inp.name][[1]]$Cell
  MPD.vis<-cellVis(cells=pcellP,value=prichness)
  
  #Func Tree
  fcellP<-MFDs[names(MFDs) %in% inp.name][[1]]$Cell
  frichness<-MFDs[names(MFDs) %in% inp.name][[1]]$MFD
  MFD.vis<-cellVis(fcellP,frichness)
  
  out<-list(richness,MPD.vis,MFD.vis)
  names(out)<-c("Richness","Phylogenetic","Trait")
  return(out)}

###Parallel Phylosor

#create parallel phylobeta functions
PD.par<-function (samp, tree, include.root = TRUE) 
{
  if (is.null(tree$edge.length)) {
    stop("Tree has no branch lengths, cannot compute pd")
  }
  species <- colnames(samp)
  SR <- rowSums(ifelse(samp > 0, 1, 0))
  nlocations = dim(samp)[1]
  nspecies = dim(samp)[2]
  PDs = NULL
  PDs<-foreach (i=1:nlocations,.combine=c,.packages="picante") %dopar% {
    present <- species[samp[i, ] > 0]
    treeabsent <- tree$tip.label[which(!(tree$tip.label %in% 
                                           present))]
    if (length(present) == 0) {
      PD <- 0
    }
    else if (length(present) == 1) {
      if (!is.rooted(tree) || !include.root) {
        warning("Rooted tree and include.root=TRUE argument required to calculate PD of single-species sampunities. Single species sampunity assigned PD value of NA.")
        PD <- NA
      }
      else {
        PD <- node.age(tree)$ages[which(tree$edge[, 
                                                  2] == which(tree$tip.label == present))]
      }
    }
    else if (length(treeabsent) == 0) {
      PD <- sum(tree$edge.length)
    }
    else {
      sub.tree <- drop.tip(tree, treeabsent)
      if (include.root) {
        if (!is.rooted(tree)) {
          stop("Rooted tree required to calculate PD with include.root=TRUE argument")
        }
        sub.tree.depth <- max(node.age(sub.tree)$ages)
        orig.tree.depth <- max(node.age(tree)$ages[which(tree$edge[, 
                                                                   2] %in% which(tree$tip.label %in% present))])
        PD <- sum(sub.tree$edge.length) + (orig.tree.depth - 
                                             sub.tree.depth)
      }
      else {
        PD <- sum(sub.tree$edge.length)
      }
    }
    return(PD)
  }
  PDout <- data.frame(PD = PDs, SR = SR)
  rownames(PDout) <- rownames(samp)
  return(PDout)
}

par.phylosor<-function (samp, tree) 
{
  if (is.null(tree$edge.length)) {
    stop("Tree has no branch lengths, cannot compute pd")
  }
  if (!is.rooted(tree)) {
    stop("Rooted phylogeny required for phylosor calculation")
  }
  samp <- as.matrix(samp)
  s <- nrow(samp)
  phylodist <- matrix(NA, s, s)
  rownames(phylodist) <- rownames(samp)
  colnames(phylodist) <- rownames(samp)
  samp_comb <- matrix(NA, s * (s - 1)/2, ncol(samp))
  colnames(samp_comb) <- colnames(samp)
  i <- 1
  for (l in 1:(s - 1)) {
    for (k in (l + 1):s) {
      samp_comb[i, ] <- samp[l, ] + samp[k, ]
      i <- i + 1
    }
  }
  pdsamp <- PD.par(samp, tree)
  pdsamp_comb <- PD.par(samp_comb, tree)
  i <- 1
  for (l in 1:(s - 1)) {
    pdl <- pdsamp[l, "PD"]
    for (k in (l + 1):s) {
      pdk <- pdsamp[k, "PD"]
      pdcomb <- pdsamp_comb[i, "PD"]
      pdsharedlk <- pdl + pdk - pdcomb
      phylodist[k, l] = 2 * pdsharedlk/(pdl + pdk)
      i <- i + 1
    }
  }
  return(as.dist(phylodist))
}

# phyl = phylo object, com = site by species matrix

# phyl = phylo object, com = site by species matrix

matpsim <- function(phyl, com, clust = 7) # make sure nodes are labelled and that com and phyl species match
{
  
  require(phylobase)  # detail information for all phylo branches
  new <- phylo4(phyl)
  dat <- as.data.frame(print(new))
  allbr <- dat$edge.length
  names(allbr) <- getEdge(new)
  spp <- colnames(com)
  
  require(foreach)
  require(doSNOW)
  
  cl <- makeCluster(clust) # create parellel clusters
  registerDoSNOW(cl)
  
  # create a list of phy branches for each species
  
  brs <-  foreach(i = spp, .packages = "phylobase") %dopar% #this loop makes a list of branches for each species
{  
  print(which(spp == i)/length(spp))
  print(date())
  brsp <- vector()
  br   <- as.numeric(rownames(dat[which(dat$label==i),]))
  repeat{
    brsn <- getEdge(new,br)
    brsl <- dat[names(brsn),"edge.length"]
    names(brsl) <- brsn
    brsp <- c(brsp, brsl)
    br   <- dat[br,3]
    if(br == 0) {break}
  }
  
  brsp
}
  names(brs) <- spp
  stopCluster(cl)
  
    
  print("brs")
  
  # create a species by phy branch matrix
  
  spp_br <- matrix(0,nrow = length(spp), ncol = length(allbr))
  rownames(spp_br) <- spp
  colnames(spp_br) <- names(allbr)
  
  for(i in spp)
  {
    spp_br[i,names(brs[[i]])] <- brs[[i]]
  }
  
  spp_br <- spp_br[,!colSums(spp_br) %in% c(0,nrow(spp_br))]  # here take out all common branches instead
  spp_br <- spp_br[,-(ncol(com)+1)] # removes root 
  
  
  print("spp_br")
  
  spp_br <<- spp_br
  
  # function to give the pres/abs of phy branches within cell i
  
  cellbr <- function(i,spp_br, com)
  {
    i_spp <- rownames(spp_br)[com[i,]>0]
    if(length(i_spp) > 1)
    {i_br  <- as.numeric(apply(spp_br[i_spp,],2,max))}
    else  {i_br  <- as.numeric(spp_br[i_spp,]) }
    names(i_br) <- colnames(spp_br)
    return(i_br)
  }
  
  require(foreach)
  require(doSNOW)
  
  cl <- makeCluster(clust) # create parellel clusters
  registerDoSNOW(cl)
  
  tcellbr <- foreach(j = rownames(com), .combine = "rbind") %dopar% {cellbr(j,spp_br,com)}
  
  stopCluster(cl)
  
  print("cell_br")
  rownames(tcellbr) <- rownames(com)
  tcellbr <<- tcellbr
  
  # function to calculate phylobsim between cell_a and cell_b
  
  nmatsim <- function(cell_b) # samp = grid cell of interest
  {
    a_br  <- tcellbr[cell_a,]
    b_br <- tcellbr[cell_b,]
    s_br <- rbind(a_br,b_br)
    s_br <- as.matrix(s_br[,colSums(s_br) >0])
    pa_br <- s_br > 0
    both <- s_br[1,colSums(pa_br > 0)==2]
    ubr <- s_br[,colSums(pa_br > 0)==1]
    a_ubr <- as.matrix(ubr)[1,]
    b_ubr <- as.matrix(ubr)[2,]
    psim <- 1 - (sum(both,na.rm=T)/(min(sum(a_ubr,na.rm=T),sum(b_ubr,na.rm=T))+sum(both,na.rm=T)))
    return(psim)
  }
  
  # calculate full cell by cell phylobsim matrix
  
  cl <- makeCluster(clust) # make another cluster
  registerDoSNOW(cl)
  
  psim <- foreach(j = rownames(tcellbr)) %dopar% {
    print(j);cell_a <- j; 
    unlist(
      lapply(
        rownames(tcellbr)[1:(which(rownames(tcellbr) == j))], 
        nmatsim
      )
    )
  }
  
  stopCluster(cl)
  print("psim")
  psim <- do.call("rbind", psim)
  rownames(psim) <- rownames(tcellbr)
  colnames(psim) <- rownames(tcellbr)    
  psim[upper.tri(psim)] <- psim[lower.tri(psim)]
  psim <- as.dist(psim)
  print("its over")
  return(psim)
  }

######################
#Non parallel version
#####################
matpsim.nopar <- function(phyl, com, clust = 7) # make sure nodes are labelled and that com and phyl species match
{
  
  require(phylobase)  # detail information for all phylo branches
  new <- phylo4(phyl)
  dat <- as.data.frame(print(new))
  allbr <- dat$edge.length
  names(allbr) <- getEdge(new)
  spp <- colnames(com)
  
  require(foreach)
  require(doSNOW)
  
  #cl <- makeCluster(clust) # create parellel clusters
  #registerDoSNOW(cl)
  
  # create a list of phy branches for each species
  
  brs <-  foreach(i = spp, .packages = "phylobase") %do% #this loop makes a list of branches for each species
{  
  print(which(spp == i)/length(spp))
  print(date())
  brsp <- vector()
  br   <- as.numeric(rownames(dat[which(dat$label==i),]))
  repeat{
    brsn <- getEdge(new,br)
    brsl <- dat[names(brsn),"edge.length"]
    names(brsl) <- brsn
    brsp <- c(brsp, brsl)
    br   <- dat[br,3]
    if(br == 0) {break}
  }
  
  brsp
}
  names(brs) <- spp
  #stopCluster(cl)
  
  
  print("brs")
  
  # create a species by phy branch matrix
  
  spp_br <- matrix(0,nrow = length(spp), ncol = length(allbr))
  rownames(spp_br) <- spp
  colnames(spp_br) <- names(allbr)
  
  for(i in spp)
  {
    spp_br[i,names(brs[[i]])] <- brs[[i]]
  }
  
  spp_br <- spp_br[,!colSums(spp_br) %in% c(0,nrow(spp_br))]  # here take out all common branches instead
  spp_br <- spp_br[,-(ncol(com)+1)] # removes root
  
  
  print("spp_br")
  
  spp_br <<- spp_br
  
  # function to give the pres/abs of phy branches withi cell i
  
  cellbr <- function(i,spp_br, com)
  {
    i_spp <- rownames(spp_br)[com[i,]>0]
    if(length(i_spp) > 1)
    {i_br  <- as.numeric(apply(spp_br[i_spp,],2,max))}
    else  {i_br  <- as.numeric(spp_br[i_spp,]) }
    names(i_br) <- colnames(spp_br)
    return(i_br)
  }
  
  require(foreach)
  require(doSNOW)
  
  #cl <- makeCluster(clust) # create parellel clusters
  #registerDoSNOW(cl)
  
  tcellbr <- foreach(j = rownames(com), .combine = "rbind") %do% {cellbr(j,spp_br,com)}
  
  #stopCluster(cl)
  
  print("cell_br")
  rownames(tcellbr) <- rownames(com)
  tcellbr <<- tcellbr
  
  # function to calculate phylobsim between cell_a and cell_b
  
  nmatsim <- function(cell_b) # samp = grid cell of interest
  {
    a_br  <- tcellbr[cell_a,]
    b_br <- tcellbr[cell_b,]
    s_br <- rbind(a_br,b_br)
    s_br <- as.matrix(s_br[,colSums(s_br) >0])
    pa_br <- s_br > 0
    both <- s_br[1,colSums(pa_br > 0)==2]
    ubr <- s_br[,colSums(pa_br > 0)==1]
    a_ubr <- as.matrix(ubr)[1,]
    b_ubr <- as.matrix(ubr)[2,]
    psim <- 1 - (sum(both,na.rm=T)/(min(sum(a_ubr,na.rm=T),sum(b_ubr,na.rm=T))+sum(both,na.rm=T)))
    return(psim)
  }
  
  # calculate full cell by cell phylobsim matrix
  
  #cl <- makeCluster(clust) # make another cluster
  #registerDoSNOW(cl)
  
  psim <- foreach(j = rownames(tcellbr)) %do% {
    print(j);cell_a <- j; 
    unlist(
      lapply(
        rownames(tcellbr)[1:(which(rownames(tcellbr) == j))], 
        nmatsim
      )
    )
  }
  
  #stopCluster(cl)
  print("psim")
  psim <- do.call("rbind", psim)
  rownames(psim) <- rownames(tcellbr)
  colnames(psim) <- rownames(tcellbr)    
  psim[upper.tri(psim)] <- psim[lower.tri(psim)]
  psim <- as.dist(psim)
  print("its over")
  return(psim)
}


##################
#Pairwise Version
##################

matpsim.pairwise <- function(phyl, com.x, com.y, clust = 7) # make sure nodes are labelled and that com and phyl species match
{
  
  require(phylobase)  # detail information for all phylo branches
  new <- phylo4(phyl)
  dat <- as.data.frame(print(new))
  allbr <- dat$edge.length
  names(allbr) <- getEdge(new)
  #Assumes same species in the comm matrices
  spp <- colnames(com.x)
  
  require(foreach)
  require(doSNOW)
  
  cl <- makeCluster(clust) # create parellel clusters
  registerDoSNOW(cl)
  
  # create a list of phy branches for each species
  
  brs <-  foreach(i = spp, .packages = "phylobase") %dopar% #this loop makes a list of branches for each species
{  
  print(which(spp == i)/length(spp))
  print(date())
  brsp <- vector()
  br   <- as.numeric(rownames(dat[which(dat$label==i),]))
  repeat{
    brsn <- getEdge(new,br)
    brsl <- dat[names(brsn),"edge.length"]
    names(brsl) <- brsn
    brsp <- c(brsp, brsl)
    br   <- dat[br,3]
    if(br == 0) {break}
  }
  
  brsp
}
  names(brs) <- spp
  stopCluster(cl)

  # create a species by phy branch matrix
  
  spp_br <- matrix(0,nrow = length(spp), ncol = length(allbr))
  rownames(spp_br) <- spp
  colnames(spp_br) <- names(allbr)
  
  for(i in spp)
  {
    spp_br[i,names(brs[[i]])] <- brs[[i]]
  }
  
  spp_br <- spp_br[,!colSums(spp_br) %in% c(0,nrow(spp_br))]  # here take out all common branches instead
  spp_br <- spp_br[,-(ncol(com.x)+1)] # removes root
  print("spp_br")
  
  spp_br <<- spp_br
  
  # function to give the pres/abs of phy branches withi cell i
  
  cellbr <- function(i,spp_br, com)
  {
    i_spp <- rownames(spp_br)[com[i,]>0]
    if(length(i_spp) > 1)
    {i_br  <- as.numeric(apply(spp_br[i_spp,],2,max))}
    else  {i_br  <- as.numeric(spp_br[i_spp,]) }
    names(i_br) <- colnames(spp_br)
    return(i_br)
  }
  
  require(foreach)
  require(doSNOW)
  
  cl <- makeCluster(clust) # create parellel clusters
  registerDoSNOW(cl)
  
  #Looks like we need to duplicate this function! Once for the current matrix, and one for the future matrix
  
  tcellbr.x <- foreach(j = rownames(com.x), .combine = "rbind") %dopar% {cellbr(j,spp_br,com.x)}
  tcellbr.y <- foreach(j = rownames(com.y), .combine = "rbind") %dopar% {cellbr(j,spp_br,com.y)}
  
  stopCluster(cl)
  
  rownames(tcellbr.x) <- rownames(com.x)
  rownames(tcellbr.y) <- rownames(com.y)
  #tcellbr <<- tcellbr
  
  # function to calculate phylobsim between cell_a and cell_b
  
  nmatsim <- function(cell_a,cell_b) # samp = grid cell of interest
  {
    a_br  <- tcellbr.x[cell_a,]
    b_br <- tcellbr.y[cell_b,]
    s_br <- rbind(a_br,b_br)
    s_br <- as.matrix(s_br[,colSums(s_br) >0])
    pa_br <- s_br > 0
    both <- s_br[1,colSums(pa_br > 0)==2]
    ubr <- s_br[,colSums(pa_br > 0)==1]
    a_ubr <- as.matrix(ubr)[1,]
    b_ubr <- as.matrix(ubr)[2,]
    res <- (sum(a_ubr,na.rm=T) + sum(b_ubr,na.rm=T))/(2*(sum(both, na.rm=T)) + sum(a_ubr,na.rm=T) + sum(b_ubr,na.rm=T)
    return(res)
  }
  
  # calculate full cell by cell phylobsim matrix
  
  cl <- makeCluster(clust) # make another cluster
  registerDoSNOW(cl)
  
  psor<-foreach(j=rownames(tcellbr.x)) %dopar%{
    sapply(rownames(tcellbr.y),function(k){
      nmatsim(j,k)
    })}
  
  stopCluster(cl)
  
  psor <- do.call("rbind", psor)
  rownames(psor) <- rownames(tcellbr.x)
  colnames(psor) <- rownames(tcellbr.y)    
  
  return(psor)
}

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


###Between time func

MNND_fc <- function(fu,cur,sp.list_current,sp.list_future,dists)
{
  Asp     <- sp.list_current[[fu]]
  Bsp     <- sp.list_future[[cur]]
  
  if(is.null(Asp) | is.null(Bsp)){
    res <- 1
    #names(res) <- "MNND"
  } else {
    compmat <- dists[Asp,Bsp]
    Ann     <- apply(as.matrix(compmat),1,min)
    Bnn     <- apply(as.matrix(compmat),2,min)
    Dnn     <- mean(c(Ann, Bnn))
    #turn    <- min(c(mean(Ann),mean(Bnn)))
    #nest    <- Dnn - turn
    #res <- c(Dnn,turn,nest)
    res<-c(Dnn)
    #names(res) <- c("MNND","MNNDturn","MNNDnest")
    #names(res) <- c("MNND")  
  }
  
  return(res)
}


trait.betadiv <- function(cur, fu, current.func, func.dat, beta.measure,
                          sp.list_current, sp.list_future, traits) {
  species.dat <- unique(unlist(c(sp.list_current[cur], sp.list_future[fu])))
  comm.dat <- data.frame(cbind(current.func[cur,], func.dat[fu,]))
  comm.dat <- comm.dat[rownames(comm.dat) %in% species.dat,]
  comm.dat <- t(comm.dat)
  trait.dat <- as.matrix(traits[rownames(traits) %in% species.dat,])
  betadiv <- functional.beta.pair(comm.dat, trait.dat)
  res <- switch(beta.measure, "sor" = betadiv$funct.beta.sor, "sim" = betadiv$funct.beta.sim)
  return(as.vector(res))
}

