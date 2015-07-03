#Define Function

require(picante)
require(foreach)
require(doSNOW)
require(phylobase)
#Set dropbox path

droppath<-"/home1/02443/bw4sz/GlobalMammal/"


#Read in species matrix
siteXspp <- read.csv(paste(droppath,"ninexdat.csv",sep=""))

#Just get the species data, starts on column 33 for this example
siteXspp<-siteXspp[,33:ncol(siteXspp)]

#Remove lines with less than 2 species
richness<-apply(siteXspp,1,sum)
keep<-which(richness > 2)
siteXspp<-siteXspp[keep,]

#Get entire species list
splist<-colnames(siteXspp)

#Read in phylogeny
tree<-read.tree(paste(droppath,"Sep19_InterpolatedMammals_ResolvedPolytomies.nwk",sep=""))

#remove species in the siteXspp that are not in phylogeny
siteXspp<-siteXspp[,colnames(siteXspp) %in% tree$tip.label]


#Define Beta Function
matpsim <- function(phyl, com) # make sure nodes are labelled and that com and phyl species match
{
  
  require(phylobase)  # detail information for all phylo branches
  new <- phylo4(phyl)
  dat <- as.data.frame(print(new))
  allbr <- dat$edge.length
  names(allbr) <- getEdge(new)
  spp <- colnames(com)
  
  require(foreach)
  require(doSNOW)
  
  cl <- makeSOCKCluster() # create parellel clusters
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
  
  spp_br <- spp_br[,-(ncol(com)+1)] # removes root
  spp_br <- spp_br[,!colSums(spp_br) %in% c(0,nrow(spp_br))]  # here take out all common branches instead
  
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
  
  cl <- makeSOCKCluster() # create parellel clusters
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
  #psim <- as.dist(psim)
  print("its over")
  return(psim)
}

#compute
betaS<-matpsim(tree,siteXspp)

#Write to file
write.csv(betaS,paste(droppath,"BetaSim.csv",sep=""))
