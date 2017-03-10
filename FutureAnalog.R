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
  morph <- read.csv("InputData/MorphologyShort.csv", na.strings="9999")
  
  #just get males & the 3 traits of interest
  mon <- filter(morph, Sex == "Macho") %>%
    select(SpID, ExpC, Peso, AlCdo) %>%
    group_by(SpID) %>%
    summarise_each(funs(mean(., na.rm = TRUE))) %>%
    filter(complete.cases(.))
  
  mon <- data.frame(mon)
  colnames(mon) <- c("Species","Bill","Mass","WingChord")
  rownames(mon) <- gsub(" ",".",mon$Species)
  mon <- mon[,-1]
  
  #principal component traits and get euclidean distance matrix
  means <- apply(mon, 2, mean)
  
  Bill <- (mon$Bill - means["Bill"])/sd(mon$Bill)
  Mass <- (mon$Mass - means["Mass"])/sd(mon$Mass)
  WingChord <- (mon$WingChord - means["WingChord"])/sd(mon$WingChord)
  
  z.scores <- data.frame(Bill, Mass, WingChord)
  rownames(z.scores) <- rownames(mon)
  
  fco <- as.matrix(dist(z.scores, method = "euclidean"))
  
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
  
  # get the crop files
  niche.crops <- list.files(out_path,pattern="crop.gri",full.name=T,recursive=T)

  
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
  
  # Step 5) calculate current null model (current-to-current betadiversity)
  # taxonomic beta-diversity
  beta.time.taxa.cnull <- analogue::distance(current, current, "bray")
  
  # phylogenetic beta-diversity
  beta.time.phylo.cnull <- matpsim.pairwise(phyl = trx, 
                                      com.x = current.phylo, 
                                      com.y = current.phylo)
  
  # functional beta-diversity 
  sp.list <- lapply(rownames(current.func), function(k){
    g <- current.func[k,]
    names(g[which(g==1)])
  })
  
  names(sp.list) <- rownames(current.func)
  
  #Get distances from the cophenetic matrix?
  dists <- as.matrix(fco)
  rownames(dists) <- rownames(fco)
  colnames(dists) <- rownames(fco)
  
  cl <- makeCluster(clust) # make another cluster
  registerDoSNOW(cl)
  
  beta.time.func.cnull <- foreach(a=rownames(current.func), .export = c("current.func", "MNND_fc"), .combine="cbind") %dopar%{
    sapply(rownames(current.func), function(b){
      MNND_fc(a, b, sp.list, sp.list, dists)
    })}
  
  stopCluster(cl)
  
  rownames(beta.time.func.cnull) <- rownames(current.func)
  colnames(beta.time.func.cnull) <- rownames(current.func)

  res <- list(beta.time.taxa.cnull, beta.time.phylo.cnull, beta.time.func.cnull)
  names(res) <- c("beta.time.taxa.cnull", "beta.time.phylo.cnull", "beta.time.func.cnull")
  save(res, file = paste0(out_path, "/beta_diversity_cnull.rda"))
  
  # Get climate models in use
  clim.mods <- list.files("../worldclim_data/projections_2070/")
  #for now
  #clim.mods <- list.files("sppXsite", full.names = TRUE)
  #clim.mods_not <- list.files("sppXsite", full.names = TRUE, pattern = "current")
  #clim.mods <- setdiff(clim.mods, clim.mods_not)
  
  for(mod in clim.mods){
    niche <- niche.crops[grep(mod,niche.crops,value=FALSE)]
    
    #Create siteXspp table from input rasters, function is from
    #AlphaMappingFunctions.R, sourced at the top.
    siteXspps <- tableFromRaster(niche, threshold = 0.05)
    
    #Remove NAs from siteXspps 
    fails <- na.test(siteXspps)
    siteXspps <- siteXspps[,!colnames(siteXspps) %in% fails]
    
    # Step 1) TAXONOMIC BETA DIVERSITY ---------------------------------------------
    beta.time.taxa <- analogue::distance(current, siteXspps, "bray")
    
    # future null 
    beta.time.taxa.fnull <- analogue::distance(siteXspps, siteXspps, "bray")
    # Step 2) PHYLO BETA DIVERSITY -------------------------------------------------
    # For phylobeta, there needs to be more than 2 species for a rooted tree
    phylo.dat <- siteXspps[,colnames(siteXspps) %in% trx$tip.label]
    phylo.dat <- phylo.dat[!rowSums(phylo.dat)<=2,]   
    
    strt <- Sys.time()
    beta.time.phylo <- matpsim.pairwise(phyl = trx, 
                                        com.x = current.phylo, 
                                        com.y = phylo.dat)
    
    # future null
    beta.time.phylo.fnull <- matpsim.pairwise(phyl = trx, 
                                        com.x = phylo.dat, 
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
    
    sp.list_future <- lapply(rownames(func.dat), function(k){
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
    
    beta.time.func <- foreach(fu=rownames(func.dat), .export = c("MNND_fc"), .combine="cbind") %dopar%{
      sapply(rownames(current.func), function(cur){
        MNND_fc(fu, cur, sp.list_current, sp.list_future, dists)
      })}
    
    
    stopCluster(cl)
    
    rownames(beta.time.func) <- rownames(current.func)
    colnames(beta.time.func) <- rownames(func.dat)
    
    cl <- makeCluster(clust) # make another cluster
    registerDoSNOW(cl)
    
    # future null
    beta.time.func.fnull <- foreach(a=rownames(func.dat), .export = c("func.dat", "MNND_fc"), .combine="cbind") %dopar%{
      sapply(rownames(func.dat), function(b){
        MNND_fc(a, b, sp.list_future, sp.list_future, dists)
      })}
    
    stopCluster(cl)
    
    rownames(beta.time.func.fnull) <- rownames(func.dat)
    colnames(beta.time.func.fnull) <- rownames(func.dat)
    
    res <- list(beta.time.taxa, beta.time.phylo, beta.time.func)
    names(res) <- c("beta.time.taxa", "beta.time.phylo", "beta.time.func")
    save(res, file = paste0(out_path, "/beta_diversity_", mod, ".rda"))
    
    res_fnull <- list(beta.time.taxa.fnull, beta.time.phylo.fnull, beta.time.func.fnull)
    names(res_fnull) <- c("beta.time.taxa.fnull", "beta.time.phylo.fnull", "beta.time.func.fnull")
    save(res_fnull, file = paste0(out_path, "/beta_diversity_fnull_", mod, ".rda"))
  }
}

# PART II: ANALOG ANALYSIS ---------------------------------------------------
runAnalogAnalysis <- function(arbthresh, out_path) {
  # get the crop files
  niche.crops <- list.files(out_path,pattern="crop.gri",full.name=T,recursive=T)
  
  
  # create a blank raster object of the correct size and extent to have for
  # projecting the cell values
  blank <- raster(niche.crops[[1]])
  
  # get list of results from beta diversity analysis
  betadiv.files <- list.files(out_path, pattern = "fnull", full.name = TRUE)
  if(!dir.exists(paste(out_path, arbthresh, sep = "/"))) dir.create(paste(out_path, arbthresh, sep="/"))
  
  # load the current null expectation
  load(paste0(out_path, "/beta_diversity_cnull.rda"))
  res_cnull <- res
  
  for(f in betadiv.files){
    load(f)
    f <- gsub("fnull_", "", f)
    load(f)
    mod <- substr(f, nchar(f)-11, nchar(f)-4)
    
    # Step 1) CURRENT COMMUNITIES WITHOUT ANALOGS IN FUTURE (Disappearing) -------
    
    # How many current communities have do not have analogs in the future?
    # These are akin to communities which will disappear, "Disappearing"
    
    # For each of the current communities how many future communities fall below the
    # threshold (e.g., are analogous; 0 = similar, 1 = different)
    c_f_tax <- fnCurrent2Future(res$beta.time.taxa, arbthresh, blank)
    c_f_phylo <- fnCurrent2Future(res$beta.time.phylo, arbthresh, blank)
    c_f_func <- fnCurrent2Future(res$beta.time.func, arbthresh, blank)
    
    # Create output raster stack (Disappearing) 
    c_f <- stack(c(c_f_tax,c_f_phylo,c_f_func))
    names(c_f) <- c("Taxonomic", "Phylogenetic", "Functional")
    
    # calculate the null expectation
    c_c_tax <- fnCurrent2Future(res_cnull$beta.time.taxa.cnull, arbthresh, blank)
    c_c_phylo <- fnCurrent2Future(res_cnull$beta.time.phylo.cnull, arbthresh, blank)
    c_c_func <- fnCurrent2Future(res_cnull$beta.time.func.cnull, arbthresh, blank)
    
    # calculate corrected number of analogs
    c_f_tax_null <- c_f_tax / c_c_tax
    c_f_phylo_null <- c_f_phylo / c_c_phylo
    c_f_func_null <- c_f_func / c_c_func
    
    # stack
    c_f_null <- stack(c(c_f_tax_null, c_f_phylo_null, c_f_func_null))
    names(c_f_null) <- c("Taxonomic", "Phylogenetic", "Functional")
    
    # Step 2) FUTURE COMMUNITIES WITHOUT ANALOGS IN CURRENT (Novel) --------------
    
    #How many future communities do not have analogs in the current time?
    #These are akin to communities that are novel, "Novel"
    f_c_tax <- fnFuture2Current(res$beta.time.taxa, arbthresh, blank)
    f_c_phylo <- fnFuture2Current(res$beta.time.phylo, arbthresh, blank)
    f_c_func <- fnFuture2Current(res$beta.time.func, arbthresh, blank)
    
    # Visualize all three NON-ANALOGS together --------------------------
    f_c <- stack(c(f_c_tax,f_c_phylo,f_c_func))
    names(f_c) <- c("Taxonomic", "Phylogenetic", "Functional")
    
    # calculate null expectation
    f_f_tax <- fnFuture2Current(res_fnull$beta.time.taxa.fnull, arbthresh, blank)
    f_f_phylo <- fnFuture2Current(res_fnull$beta.time.phylo.fnull, arbthresh, blank)
    f_f_func <- fnFuture2Current(res_fnull$beta.time.func.fnull, arbthresh, blank)
    
    # calculate corrected number of analogs
    f_c_tax_null <- f_c_tax / f_f_tax
    f_c_phylo_null <- f_c_phylo / f_f_phylo
    f_c_func_null <- f_c_func / f_f_func
    
    # stack
    f_c_null <- stack(c(f_c_tax_null, f_c_phylo_null, f_c_func_null))
    names(f_c_null) <- c("Taxonomic", "Phylogenetic", "Functional")
    
    # output the raster stacks
    results <- stack(c_f, f_c)
    names(results) <- c(paste("Novel",c("Tax","Phylo","Func")), 
                        paste("Disappearing",c("Tax","Phylo","Func")))
    
    save(results, file=paste0(out_path, "/", arbthresh, "/NonAnalogRasters_", mod, ".rda"))
    
    # output null raster stacks
    results <- stack(c_f_null, f_c_null)
    names(results) <- c(paste("Novel",c("Tax","Phylo","Func")), 
                        paste("Disappearing",c("Tax","Phylo","Func")))
    
    save(results, file=paste0(out_path, "/", arbthresh, "/NonAnalogRastersNull_", mod, ".rda"))
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
