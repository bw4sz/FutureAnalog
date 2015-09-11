# FutureAnalog.R ---------------------------------------------------------------

#This code goes through the results from AlphaMapping.R to determine the 
# number of analog hummingbird assemblages in Ecuador under future climate scenarios.
# PART I CALCULATE BETWEEN TIME BETA-DIVERSITY MEASURES ------------------------
runBetaDiv <- function(out_path, cell_size, clust = 7){

  # Step 1) Bring in Phylogenetic Data -------------------------------------------
  trx<-read.tree("InputData/hum294.tre")
  new<-str_extract(trx$tip.label,"(\\w+).(\\w+)")
  trx<-drop.tip(trx,trx$tip.label[duplicated(new)])
  trx$tip.label<-str_extract(trx$tip.label,"(\\w+).(\\w+)")

  # Step 2) Bring in trait data --------------------------------------------------
  traits <- getTraitData()
  
  # Step 3) Bring in niche models ------------------------------------------------
  all.niche <- list.files(out_path, pattern="ensemble.gri",full.name=T,recursive=T)
  
  # Clip to Extent and shape of desired countries (Ecuador for now)
  ec<-readOGR("InputData", "EcuadorCut")
  r<-raster(extent(ec))
  
  # Match cell size above from the SDM_SP function
  res(r) <- cell_size
  ec.r <- rasterize(ec,r)
  
  niche.crop <- lapply(all.niche,function(x){
    r <- crop(raster(x),extent(ec.r))
    filnam <- paste(strsplit(x,".gri$")[[1]][1],"crop",sep="")
    writeRaster(r,filnam,overwrite=TRUE)
  })
  
  # select only well fitting species (currently set as all where TSS is over 0.5
  # for all models and ROC is over 7.5 for all models)
  model_eval<-list.files(out_path, full.name=TRUE,recursive=T,pattern="Eval.csv")
  model_eval<-rbind_all(lapply(model_eval, 
                               function(x) read.csv(x, stringsAsFactors = FALSE)))
  colnames(model_eval)[1:2] <- c("Stat", "Species")

    # get the crop files
  niche.crops <- list.files(out_path,pattern="crop.gri",full.name=T,recursive=T)
  
  # create a blank raster object of the correct size and extent to have for
  # projecting the cell values
  blank <- raster(niche.crops[[1]])
  
  # Step 4) Get current niches for comparing to ----------------------------------
  current <- niche.crops[grep("current", niche.crops, value=FALSE)]
  current <- tableFromRaster(current, threshold = 0.05)
  
  #Remove NAs from current so we can do the following analyses. Some species do
  #not occur in Ecuador, so they should be removed from analysis here.
  fails <- na.test(current)
  current <- current[,!colnames(current) %in% fails]
  
  current.phylo <- current[,colnames(current) %in% trx$tip.label]
  current.phylo <- current.phylo[!rowSums(current.phylo)<=2,]  
  
  current.func <- current[,colnames(current) %in% colnames(fco)]
  current.func <- current.func[!apply(current.func,1,sum)<=2,]  
  
  # list to output results
  NonAnalogRasters <- list()
  
  # Get climate models in use
  clim.mods <- list.files("../worldclim_data/projections_2070/")
  
  for(mod in clim.mods){
    niche <- niche.crops[grep(mod,niche.crops,value=FALSE)]
    
    #Create siteXspp table from input rasters, function is from
    #AlphaMappingFunctions.R, sourced at the top.
    siteXspps <- tableFromRaster(niche, threshold = 0.05)
    
    #Remove NAs from siteXspps 
    fails <- na.test(siteXspps)
    siteXspps <- siteXspps[,!colnames(siteXspps) %in% fails]
    
    # Step 1) TAXONOMIC BETA DIVERSITY ---------------------------------------------
    tsor <- analogue::distance(current, siteXspps, "bray")
    
    tsim <- proxy::dist(current, siteXspps, "Simpson")
    
    tnes <- tsor - tsim
    
    # Step 2) PHYLO BETA DIVERSITY -------------------------------------------------
    # For phylobeta, there needs to be more than 2 species for a rooted tree
    phylo.dat <- siteXspps[,colnames(siteXspps) %in% trx$tip.label]
    phylo.dat <- phylo.dat[!rowSums(phylo.dat)<=2,]   
    
    strt <- Sys.time()
    pbeta <- matpsim.pairwise(phyl = trx, 
                                        com.x = current.phylo, 
                                        com.y = phylo.dat)
    
    psor <- pbeta$psor
    psim <- pbeta$psim
    pnes <- pbeta$pnes
    
    Sys.time()-strt
    
    # Step 3) FUNC BETA DIVERSITY ---------------------------------------------------
    func.dat <- siteXspps[,colnames(siteXspps) %in% colnames(fco)]
    func.dat <- func.dat[!apply(func.dat,1,sum)<=2,]  
    
    sp.list_current <- lapply(rownames(current.func), function(k){
      g <- current.func[k,]
      names(g[which(g==1)])
    })
    
    names(sp.list_current) <- rownames(current.func)
    
    sp.list_future <- lapply(rownames(func.dat), function(k){
      g <- func.dat[k,]
      names(g[which(g==1)])
    }) 
    
    names(sp.list_future) <- rownames(func.dat)

    site.pairs <- expand.grid(rownames(current.func), rownames(func.dat))
    betadiv <- apply(site.pairs, 1, function(x) {
      trait.betadiv(x[1], x[2], current.func, func.dat, "nes",
                    sp.list_current, sp.list_future, traits)
    })
    
    betadiv <- t(betadiv)
    betadiv <- cbind(site.pairs, betadiv)
    
    tsor <- spread(betadiv[-c(3,4)], Var2, funct.beta.sor)
    rownames(tsor) <- tsor[,1]
    tsor <- tsor[-1]
    
    tsim <- spread(betadiv[-c(4,5)], Var2, funct.beta.sim)
    rownames(tsim) <- tsim[,1]
    tsim <- tsim[-1]
    
    tnes <- spread(betadiv[-c(3,5)], Var2, funct.beta.sne)
    rownames(tnes) <- tnes[,1]
    tnes <- tnes[-1]
    
    res <- list(tsor=tsor, tsim=tsim, tnes=tnes, 
                psor=psor, psim=psim, pnes=pnes, 
                fsor=fsor, fsim=fsim, fnes=fnes)
    
    save(res, file = paste0(out_path, "/beta_diversity_", mod, ".rda"))
  }
}

# PART II: ANALOG ANALYSIS ---------------------------------------------------
runAnalogAnalysis <- function(arbthresh, out_path) {
  # get list of results from beta diversity analysis
  betadiv.files <- list.files(out_path, pattern = "beta_diversity", full.name = TRUE)
  if(!dir.exists(paste(out_path, arbthresh, sep = "/"))) dir.create(paste(out_path, arbthresh, sep="/"))

  for(f in betadiv.files){
    load(f)
    mod <- substr(f, nchar(f)-11, nchar(f)-4)
    
    # Step 1) CURRENT COMMUNITIES WITHOUT ANALOGS IN FUTURE (Disappearing) -------
    
    # How many current communities have do not have analogs in the future?
    # These are akin to communities which will disappear, "Disappearing"
    
    # For each of the current communities how many future communities fall below the
    # threshold (e.g., are analogous; 0 = similar, 1 = different)
    c_f_tax <- fnCurrent2Future(res$beta.time.taxa, arbthresh)
    c_f_phylo <- fnCurrent2Future(res$beta.time.phylo, arbthresh)
    c_f_func <- fnCurrent2Future(res$beta.time.func, arbthresh)
    
    # Create output raster stack (Disappearing) 
    c_f <- stack(c(c_f_tax,c_f_phylo,c_f_func))
    names(c_f) <- c("Taxonomic", "Phylogenetic", "Functional")
    
    # Step 2) FUTURE COMMUNITIES WITHOUT ANALOGS IN CURRENT (Novel) --------------
    
    #How many future communities do not have analogs in the current time?
    #These are akin to communities that are novel, "Novel"
    f_c_tax <- fnFuture2Current(res$beta.time.taxa, arbthresh)
    f_c_phylo <- fnFuture2Current(res$beta.time.phylo, arbthresh)
    f_c_func <- fnFuture2Current(res$beta.time.func, arbthresh)
    
    # Visualize all three NON-ANALOGS together --------------------------
    f_c <- stack(c(f_c_tax,f_c_phylo,f_c_func))
    names(f_c) <- c("Taxonomic", "Phylogenetic", "Functional")
    
    # output the raster stacks
    results <- stack(c_f, f_c)
    names(results) <- c(paste("Novel",c("Tax","Phylo","Func")), 
                        paste("Disappearing",c("Tax","Phylo","Func")))
    
    save(results, file=paste0(out_path, "/", arbthresh, "/NonAnalogRasters_", mod, ".rda"))
  }
}

# Misc functions needed for the above functions
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











## *** CODE CHECKED TO HERE *** ################################################



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

#save.image("FutureAnalog.rData")
#s amnat
