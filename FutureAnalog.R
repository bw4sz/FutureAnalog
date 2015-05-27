# FutureAnalog.R ---------------------------------------------------------------

#This code goes through the results from AlphaMapping.R to determine the 
# number of analog hummingbird assemblages in Ecuador under future climate scenarios.
packages <- c("vegan", "picante", "reshape", "reshape2", "analogue", "doSNOW", "ape", "cluster", 
         "RColorBrewer", "raster", "ggplot2", "phylobase", "rgdal")

for(p in packages) {
  if (!p %in% installed.packages()) {
    install.packages(p)
  }
  require(p, character.only = TRUE)
}

#Load in source functions
source("AlphaMappingFunctions.R")
source("BenHolttraitDiversity.R")

# Set output folder
# set the cell size for the analysis - **DECISION**
cell_size = 0.1 

# create folders to output the models to
output_folder = "../FutureAnalog_output" 
out_path <- paste(output_folder, cell_size, sep = "/")

# Step 1) Load results from AlphaMapping.R -------------------------------------
load(paste(out_path, "siteXspps.rda", sep = "/"))
current <- siteXspps[[1]]
future <- siteXspps[2:4]

#Remove NAs from siteXspps so we can do the following analyses. Some species do
#not occur in Ecuador, so they should be removed from analysis here.
na.test <-  function (x) {
  w <- apply(x, 2, function(x)all(is.na(x)))
  if (any(w)) {
    fails <- names(which(w))
    print(paste("All NA in columns", paste(names(which(w)), collapse=", ")))
    return(fails)
  }
}

fails <- na.test(current)
current <- current[,!colnames(current) %in% fails]

future <- lapply(future, function(x){
  fails <- na.test(x)
  x[,!colnames(x) %in% fails]
})

# Step 2) Within time beta diversity -------------------------------------------
# Step 2a) Find within time SPECIES BETA DIVERSITY -----------------------------
within.current.dist <- vegdist(current, "bray")  # **DECISION** why Bray Curtis?
within.current <- as.matrix(within.current.dist)

within.future <- lapply(future, function(x){
  dist <- vegdist(x, "bray")
  as.matrix(dist)
})

#--------------- Find within PHYLO BETA DIVERSITY
#For phylobeta, there needs to be more than 2 species for a rooted tree
phylo.current <- current[,colnames(current) %in% trx$tip.label]
phylo.current <- phylo.current[!apply(phylo.current,1,sum)<=2,]   

phylo.future <- lapply(future, function(x){
  matched <- x[,colnames(x) %in% trx$tip.label]
  matched[!apply(matched,1,sum)<=2,] 
})


#Within current phylobetadiversity
system.time(holt.try <- matpsim(phyl=trx, com=phylo.current, clust=3))  
#turn beta measures into a matrix   
within.current.phylo <- as.matrix(holt.try)

within.future.phylo <- lapply(phylo.future, function(x){
  holt.try <- matpsim(phyl=trx, com=x, clust=3)
  as.matrix(holt.try)
})


#---------------- Find within FUNC BETA DIVERSITY

#----- Within current functional beta diversity
Func.current <- current[,colnames(current) %in% colnames(fco)]
Func.current <- Func.current[!apply(Func.current,1,sum)<=2,]  

####MNNTD method for integrating trait beta, used in the DimDiv script    
#   MNNTD = Mean nearest neighbor taxon distance                          TODO: recommend parallelizing all the func code in this script - takes HOURS to run on the fast desktop
#   Holt et al. 2012. An update of Wallace's zoogeographic regions of the world. Science.
sp.list <- lapply(rownames(Func.current), function(k){
  x <- Func.current[k,]
  names(x[which(x==1)])
})

names(sp.list) <- rownames(Func.current)
dists <- as.matrix(fco)
rownames(dists) <- rownames(fco)
colnames(dists) <- rownames(fco)

sgtraitMNTD <- sapply(rownames(Func.current), function(i){                    
  A <- i   #set iterator, takes a long time to run
  print(i)
  out <- lapply(rownames(Func.current)[1:(which(rownames(Func.current) == i))], function(B) {
    MNND(A, B, sp.list=sp.list, dists=dists)
    })
  names(out) <- rownames(Func.current)[1:(which(rownames(Func.current) == i))]
  return(out)
})

names(sgtraitMNTD) <- rownames(Func.current)  #rownames are site ID numbers
melt.MNTD <- melt(sgtraitMNTD)
colnames(melt.MNTD) <- c("MNTD","To","From")

# Turn new results into a matrix
within.current.func <- cast(melt.MNTD,To ~ From, value="MNTD")
rownames(within.current.func) <- within.current.func[,1]
within.current.func <- within.current.func[,-1]

within.current.func[lower.tri(within.current.func)] <- t(within.current.func[upper.tri(within.current.func)])


#----- Within future functional beta          
Func.future <- lapply(future, function(x){
  matched <- x[,colnames(x) %in% colnames(fco)]
  matched[!apply(matched,1,sum)<=2,]   
})

sp.list <- lapply(Func.future,function(x){
  a <- lapply(rownames(x),function(k){
    x <- x[k,]
    names(x[which(x==1)])
  })
  names(a) <- rownames(x)
})

##Turn cophenetic distance to matrix
dists <- as.matrix(fco)
rownames(dists) <- rownames(fco)
colnames(dists) <- rownames(fco)

within.future.func <- lapply(Func.future,function(x){
  sgtraitMNTD <- sapply(rownames(x),function(i){
    A<-i    #set iterator
    out<-lapply(rownames(x)[1:(which(rownames(x) == i))], function(B) {
      MNND(A, B, sp.list=sp.list, dists=dists)
      })
    names(out) <- rownames(x)[1:(which(rownames(x) == i))]
    return(out)
  })
  
  names(sgtraitMNTD) <- rownames(x)
  melt.MNTD <- melt(sgtraitMNTD)
  colnames(melt.MNTD) <- c("MNTD","To","From")
  
  #needs to be casted back into a matrix, see reshape2::dcast., name it betatime func
  within.future.func <- cast(melt.MNTD,To ~ From, value="MNTD")
  rownames(within.future.func) <- within.future.func[,1]
  within.future.func <- within.future.func[,-1]
  within.future.func[lower.tri(within.future.func)] <- t(within.future.func[upper.tri(within.future.func)])
  
  return(within.future.func)
})


###############################################################
#     BETWEEN TIME (compare current witih each future scenario)
###############################################################

#---------------- Find between TAXONOMIC BETA DIVERSITY
beta.time.taxa <- lapply(future, function(x){
  analogue::distance(current, x, "bray")
})


#---------------- Find between PHYLO BETA DIVERSITY
beta.time.phylo <- lapply(phylo.future, function(x){
  beta.time.phylo <- as.matrix(matpsim.pairwise(phyl=trx, com.x=phylo.current, com.y=x, clust=7))
})


#---------------- Find between FUNC BETA DIVERSITY
Beta.time.func <- lapply(Func.future, function(x){
  
  sp.list_current <- lapply(rownames(Func.current), function(k){
    g <- Func.current[k,]
    names(g[which(g==1)])
  })
  
  names(sp.list_current) <- rownames(Func.current)
  
  sp.list_future <- lapply(rownames(x), function(k){
    g <- x[k,]
    names(g[which(g==1)])
  })
  
  names(sp.list_future) <- rownames(x)
  
  #Get distances from the cophenetic matrix?
  dists <- as.matrix(fco)
  
  rownames(dists) <- rownames(fco)
  colnames(dists) <- rownames(fco)
  
  
  Beta.time.func <- lapply(rownames(x), function(fu){
    lapply(rownames(Func.current), function(cur){
      MNND_fc(fu, cur, sp.list_current, sp.list_future, dists)
    })
    
  })
  
  melt.MNTD <- melt(Beta.time.func)
  
  colnames(melt.MNTD) <- c("MNTD","To","From")
  
  #needs to be casted back into a matrix, see reshape2::dcast., name it betatime func
  beta.time.func <- dcast(melt.MNTD,To ~ From, value.var="MNTD")
  rownames(beta.time.func) <- beta.time.func[,1]
  beta.time.func <- beta.time.func[,-1]
  
  
  rownames(beta.time.func) <- rownames(Func.current)
  colnames(beta.time.func) <- rownames(x)
  return(beta.time.func)
})


#############################################     #TODO MONDAY: check dfs and start here
#     ANALOG ANALYSIS
#############################################
#TODO: Wrap this code in a function so we can test sensitivity of results to the threshold 
#      Set threshold for 5%, 10%, 50% and 100% - can present alternate results in appendices

#Set an arbitrary threshold      
arb.thresh <- 0.20

#################################
#PART I
#CURRENT COMMUNITIES THAT DO NOT HAVE ANALOGS IN FUTURE
#How many current communities have do not have analogs in the future?
#These are akin to communities which will disappear, "Disappearing"

#---------------- TAXONOMIC NON-ANALOGS - DISAPPEARING COMMUNITIES
#For each of the current communities how many future communities fall below the threshold (e.g., are analogous; 0 = similar, 1 = different)
current_to_future.analog <- lapply(beta.time.taxa, function(j){
  n.analogs <- sapply(rownames(j), function(x){
    sum(j[rownames(j) %in% x,] <= arb.thresh) #counts the number of assemblages with beta div values less than arb.thresh
  })
  current_to_future.analog <- data.frame(rownames(j), n.analogs)
  colnames(current_to_future.analog) <- c("cell.number", "numberofanalogs")  
  return(current_to_future.analog)
})

c_f_tax <- lapply(current_to_future.analog, function(x){
  fanalog <- cellVis(cell=x$cell.number, value=x$numberofanalogs)
  hist(x$numberofanalogs)
  writeRaster(fanalog, "NumberofFutureAnalogs_Taxon_ARB.tif", overwrite=TRUE)   
  return(fanalog)
})


#---------------- PHYLO NON-ANALOGS - DISAPPEARING COMMUNITIES
future.analog.phylo <- lapply(beta.time.phylo, function(j){
  n.analogs.phylo <- sapply(rownames(j), function(x){
    sum(j[rownames(j) %in% x,] <= arb.thresh)
  })
  
  #Create a dataframe of the number of analogs and the cellnumber
  future.analog.phylo <- data.frame(rownames(j), n.analogs.phylo)
  colnames(future.analog.phylo) <- c("cell.number", "numberofanalogs")
  return(future.analog.phylo)
})


#Visualize!
c_f_phylo <- lapply(future.analog.phylo, function(x){
  fanalog.phylo <- cellVis(cell=x$cell.number, value=x$numberofanalogs)
  hist(x$numberofanalogs)
  writeRaster(fanalog.phylo, "NumberofFutureAnalogs_Phylo_ARB.tif", overwrite=T)  
  return(fanalog.phylo)
})


#---------------- FUNC NON-ANALOGS - DISAPPEARING COMMUNITIES
future.analog.func <- lapply(Beta.time.func, function(j){
  n.analogs.func <- sapply(rownames(j), function(x){
    sum(j[rownames(j) %in% x,] <= arb.thresh)
  })
  
  future.analog.func <- data.frame(rownames(j), n.analogs.func)
  colnames(future.analog.func) <- c("cell.number", "numberofanalogs")
  return(future.analog.func)
})

c_f_func <- lapply(future.analog.func, function(x){
  fanalog.Func <- cellVis(cell=x$cell.number, value=x$numberofanalogs)
  writeRaster(fanalog.Func, "NumberofFutureAnalogs_Func_ARB.tif", overwrite=T)  
  return(fanalog.Func)
})


#############Visualize all three analog axes together
# Columns are the emissions scenarios, 
# Rows are taxonomic, phylogenetic and functional for CURRENT TO FUTURE analogs (dissapearing)
c_f <- stack(c(c_f_tax,c_f_phylo,c_f_func))

blues <- colorRampPalette(brewer.pal(9,"Blues"))(100) 
plot(c_f, col=blues)


#####################################################
#PART II
#FUTURE COMMUNITIES THAT DO NOT HAVE ANALOGS IN CURRENT
#How many future communities do not have analogs in the current time?
#These are akin to communities that are novel, "Novel"
######################################################

#---------------- TAXONOMIC NON-ANALOGS - NOVEL COMMUNITIES
#For each of the future communities how many future communities are analogous (< arb.thresh). Things that are different are "novel"
future_to_current.analog <- lapply(beta.time.taxa, function(j){
  n.analogs <- sapply(colnames(j), function(x){
    sum(j[,colnames(j) %in% x] <= arb.thresh)
  })
  
 future_to_current.analog <- data.frame(colnames(j), n.analogs)
  colnames(future_to_current.analog) <- c("cell.number", "numberofanalogs")  
  return(future_to_current.analog)
})

f_c_tax <- lapply(future_to_current.analog, function(x){
  fanalog <- cellVis(cell=x$cell.number, value=x$numberofanalogs)
  hist(x$numberofanalogs)
  writeRaster(fanalog, "NumberofCurrentAnalogs_Taxon_ARB.tif", overwrite=T)  
  return(fanalog)
})

#Data check, these should be different
RdBu <- colorRampPalette(brewer.pal(7,"RdBu"))(100)
plot(stack(f_c_tax) - stack(c_f_tax), col=RdBu)


#---------------- PHYLO NON-ANALOGS - NOVEL COMMUNITIES
future_to_current.phylo <- lapply(beta.time.phylo, function(j){
  n.analogs.phylo <- sapply(colnames(j), function(x){
    sum(j[,colnames(j) %in% x] <= arb.thresh)
  })
  
  #Create a dataframe of the number of analogs and the cellnumber
  future.analog.phylo <- data.frame(colnames(j), n.analogs.phylo)
  colnames(future.analog.phylo) <- c("cell.number", "numberofanalogs")
  return(future.analog.phylo)
})


#Visualize!
f_c_phylo <- lapply(future_to_current.phylo, function(x){
  fanalog.phylo <- cellVis(cell=x$cell.number, value=x$numberofanalogs)
  hist(x$numberofanalogs)
  writeRaster(fanalog.phylo, "NumberofCurrentAnalogs_Phylo_ARB.tif", overwrite=T)  
  return(fanalog.phylo)
})

plot(stack(f_c_phylo) - stack(c_f_phylo), col=RdBu)


#---------------- FUNC NON-ANALOGS - NOVEL COMMUNITIES
future_to_current.func <- lapply(Beta.time.func, function(j){
  n.analogs.func <- sapply(colnames(j), function(x){
    sum(j[,colnames(j) %in% x] <= arb.thresh)
  })
  
  future.analog.func <- data.frame(colnames(j), n.analogs.func)
  colnames(future.analog.func) <- c("cell.number", "numberofanalogs")
  return(future.analog.func)
})

f_c_func <- lapply(future_to_current.func, function(x){
  fanalog.Func <- cellVis(cell=x$cell.number, value=x$numberofanalogs)
  writeRaster(fanalog.Func, "NumberofCurrentAnalogs_Func_ARB.tif", overwrite=T)  
  return(fanalog.Func)
})


#############Visualize all three NON-ANALOGS together
f_c <- stack(c(f_c_tax,f_c_phylo,f_c_func))

#plot novel and disappearing assemblages
plot(f_c, col=rev(blues))  #novel
plot(f_c, col=rev(reds))   #disappearing

# #plot both as a panel              TODO: this needs to be improved - loop through emissions scenarios, plot in diff colors
# #Just try plotting one emission scenario across both disappearing and novel
# novel <- f_c[[c(1,4,7)]]
# disappear <- c_f[[c(1,4,7)]]
# 
# #This could be named correctly using 
# firstplot <- stack(novel,disappear)
# names(firstplot) <- c(paste("Novel",c("Tax","Phylo","Func")), paste("Disappearing",c("Tax","Phylo","Func")))
# cols = c(blues, blues, blues, reds, reds, reds)
# 
# blues <- colorRampPalette(brewer.pal(9,"Blues"))(100)
# reds <- colorRampPalette(brewer.pal(9,"Reds"))(100)
# plot(novel, col=rev(blues))
# plot(disappear, col=rev(reds))


#####################
#FIXME: CLEANED UNTIL HERE 8/25/2014 
#TODO: determine what the object names refer to - add loops? Finish code.
#TODO: The next sections appears to run some correlation tests and test the sensitivity of the analysis for 
#      the arbitrary threshold.  It also saves several output plots. Determine what we need, streamline and test it.
#       save results and make sure it could be run on many more scenarios.

#The rest should be fairly straightforward, correlating the rasters from above, 
#the f_c raster is the number of future analogs of current assemblages
#The c_f is the number of current analogs of future assemblages


##TODO Create a wrapper than takes in the above code as a function of the GCM layer. 
#    Maybe specify the gcm layer as a folder, ie each gcm is made up of three files, which are in their own folder, 
#     and it takes the folder name as input. 

#TODO: Make code that outputs main tables and figures for manuscript. Mean effects across all GCMs, 
#      plotted by scenario and dimension of biodiv (novel and disappearing; correlation table)

save.image("FutureAnalog.rData")
#s amnat
