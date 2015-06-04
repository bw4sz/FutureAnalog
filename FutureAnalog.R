# FutureAnalog.R ---------------------------------------------------------------

#This code goes through the results from AlphaMapping.R to determine the 
# number of analog hummingbird assemblages in Ecuador under future climate scenarios.
packages <- c("vegan", "picante", "analogue", "doSNOW", "ape", "cluster", 
         "RColorBrewer", "raster", "ggplot2", "phylobase", "rgdal", "tidyr", "stringr")

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
# set the cell size for the analysis
cell_size = 0.1 

# create folders to output the models to
output_folder = "../FutureAnalog_output" 
out_path <- paste(output_folder, cell_size, sep = "/")

# PART I: Bring in completed model, phylo and trait data -----------------------
# Step 1) Bring in Phylogenetic Data -------------------------------------------
trx<-read.nexus("InputData/ColombiaPhylogenyUM.tre")
spnames<-read.table("InputData/SpNameTree.txt" , sep = "\t", header = TRUE)

# Replace tip.label with Spnames, replace the tiplabels with periods, which is
# the biomod default Cophenetic distance is the distance between all pairs of
# species- measure of relatedness
trx$tip.label <- gsub("_",".",as.character(spnames$SpName))
co<-cophenetic(trx)

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

niche.crop <- lapply(all.niche,function(x){
  r <- crop(raster(x),extent(ec.r))
  filnam <- paste(strsplit(x,".gri")[[1]][1],"crop",sep="")
  writeRaster(r,filnam,overwrite=TRUE)
})

# get the crop files
niche.crops <- list.files(out_path,pattern="crop.gri",full.name=T,recursive=T)

# Step 4) Get current niches for comparing to ----------------------------------
current <- niche.crops[grep("current", niche.crops, value=FALSE)]
current <- tableFromRaster(current, threshold = 0.05)

#Remove NAs from current so we can do the following analyses. Some species do
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

current.phylo <- current[,colnames(current) %in% trx$tip.label]
current.phylo <- current.phylo[!rowSums(current.phylo)<=2,]  

current.func <- current[,colnames(current) %in% colnames(fco)]
current.func <- current.func[!apply(current.func,1,sum)<=2,]  

#### LOOP WILL START HERE, AT PRESENT TESTING ON:
mod <- "mc26bi70"

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
phylo.dat <- phylo.dist[!rowSums(phylo.dist)<=2,]   

beta.time.phylo <- matpsim.pairwise(phyl = trx, 
                                    com.x = current.phylo, 
                                    com.y = phylo.dat)

# Step 3) FUNC BETA DIVERSITY ---------------------------------------------------
func.dat <- siteXspps[,colnames(siteXspps) %in% colnames(fco)]
func.dat <- func.dat[!apply(func.dat,1,sum)<=2,]  

sp.list_current <- lapply(rownames(current.func), function(k){
    g <- current.func[k,]
    names(g[which(g==1)])
    })
  
names(sp.list_current) <- rownames(current.func)
  
sp.list_future <- lapply(phylo.dat, function(k){
    g <- phylo.dat[k,]
    names(g[which(g==1)])
  })
  
names(sp.list_future) <- rownames(phylo.dat)
  
#Get distances from the cophenetic matrix?
dists <- as.matrix(fco)
rownames(dists) <- rownames(fco)
colnames(dists) <- rownames(fco)
  
  
beta.time.func <- sapply(rownames(func.dat), function(fu){
    sapply(rownames(current.func), function(cur){
      MNND_fc(fu, cur, sp.list_current, sp.list_future, dists)
    })
  })
  


# PART II: ANALOG ANALYSIS -----------------------------------------------------

# TODO: Wrap this code in a function so we can test sensitivity of results to the threshold 
#      Set threshold for 5%, 10%, 50% and 100% - can present alternate results in appendices

# Set an arbitrary threshold      
arb.thresh <- 0.20

# Step 1) CURRENT COMMUNITIES WITHOUT ANALOGS IN FUTURE (Disappearing) ---------

# How many current communities have do not have analogs in the future?
# These are akin to communities which will disappear, "Disappearing"

# Step 1a) TAXONOMIC NON-ANALOGS - DISAPPEARING COMMUNITIES --------------------

# For each of the current communities how many future communities fall below the
# threshold (e.g., are analogous; 0 = similar, 1 = different)
current_to_future.analog <- lapply(beta.time.taxa, function(j){
  n.analogs <- sapply(rownames(j), function(x){
    sum(j[rownames(j) %in% x,] <= arb.thresh) #counts the number of assemblages with beta div values less than arb.thresh
  })
  current_to_future.analog <- data.frame(rownames(j), n.analogs)
  colnames(current_to_future.analog) <- c("cell.number", "numberofanalogs")  
  return(current_to_future.analog)
})

load(paste(out_path, "blank_raster.rda", sep="/")) # blank raster of correct size & extent
c_f_tax <- lapply(current_to_future.analog, function(x){
  fanalog <- cellVis(cell=x$cell.number, value=x$numberofanalogs)
  hist(x$numberofanalogs)
  writeRaster(fanalog, paste(out_path, "NumberofFutureAnalogs_Taxon_ARB.tif", 
                             sep="/"), overwrite=TRUE)   
  return(fanalog)
})

# Step 1b) PHYLO NON-ANALOGS - DISAPPEARING COMMUNITIES ------------------------
future.analog.phylo <- lapply(beta.time.phylo, function(j){
  n.analogs.phylo <- sapply(rownames(j), function(x){
    sum(j[rownames(j) %in% x,] <= arb.thresh)
  })
  
  #Create a dataframe of the number of analogs and the cellnumber
  future.analog.phylo <- data.frame(rownames(j), n.analogs.phylo)
  colnames(future.analog.phylo) <- c("cell.number", "numberofanalogs")
  return(future.analog.phylo)
})


# Visualize!
c_f_phylo <- lapply(future.analog.phylo, function(x){
  fanalog.phylo <- cellVis(cell=x$cell.number, value=x$numberofanalogs)
  hist(x$numberofanalogs)
  writeRaster(fanalog.phylo, paste(out_path, "NumberofFutureAnalogs_Phylo_ARB.tif",
                                   sep="/"), overwrite=T)  
  return(fanalog.phylo)
})


# Step 1c) FUNC NON-ANALOGS - DISAPPEARING COMMUNITIES -------------------------
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
  writeRaster(fanalog.Func, paste(out_path, "NumberofFutureAnalogs_Func_ARB.tif",
                                  sep="/"), overwrite=T)  
  return(fanalog.Func)
})


# Step 1d) Visualize all three analog axes together ----------------------------
# Columns are the emissions scenarios, Rows are taxonomic, phylogenetic and
# functional for CURRENT TO FUTURE analogs (dissapearing)
c_f <- stack(c(c_f_tax,c_f_phylo,c_f_func))

blues <- colorRampPalette(brewer.pal(9,"Blues"))(100) 
plot(c_f, col=blues)

# Step 2) FUTURE COMMUNITIES WITHOUT ANALOGS IN CURRENT (Novel) ----------------

#How many future communities do not have analogs in the current time?
#These are akin to communities that are novel, "Novel"

#Step 2a) TAXONOMIC NON-ANALOGS - NOVEL COMMUNITIES --------------------------- 
#For each of the future communities how many future communities are analogous 
# (< arb.thresh). Things that are different are "novel"
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
  writeRaster(fanalog, paste(out_path, "NumberofCurrentAnalogs_Taxon_ARB.tif", sep="/")
              , overwrite=T)  
  return(fanalog)
})

#Data check, these should be different
RdBu <- colorRampPalette(brewer.pal(7,"RdBu"))(100)
plot(stack(f_c_tax) - stack(c_f_tax), col=RdBu)


# Step 2b) PHYLO NON-ANALOGS - NOVEL COMMUNITIES -------------------------------
future_to_current.phylo <- lapply(beta.time.phylo, function(j){
  n.analogs.phylo <- sapply(colnames(j), function(x){
    sum(j[,colnames(j) %in% x] <= arb.thresh)
  })
  
  #Create a dataframe of the number of analogs and the cellnumber
  future.analog.phylo <- data.frame(colnames(j), n.analogs.phylo)
  colnames(future.analog.phylo) <- c("cell.number", "numberofanalogs")
  return(future.analog.phylo)
})


# Visualize!
f_c_phylo <- lapply(future_to_current.phylo, function(x){
  fanalog.phylo <- cellVis(cell=x$cell.number, value=x$numberofanalogs)
  hist(x$numberofanalogs)
  writeRaster(fanalog.phylo, paste(out_path, "NumberofCurrentAnalogs_Phylo_ARB.tif",
                                   sep="/"), overwrite=T)  
  return(fanalog.phylo)
})

plot(stack(f_c_phylo) - stack(c_f_phylo), col=RdBu)


# Step 2c) FUNC NON-ANALOGS - NOVEL COMMUNITIES --------------------------------
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
  writeRaster(fanalog.Func, paste(out_path, "NumberofCurrentAnalogs_Func_ARB.tif",
                                  sep = "/"), overwrite=T)  
  return(fanalog.Func)
})


# Step 2d) Visualize all three NON-ANALOGS together ----------------------------
f_c <- stack(c(f_c_tax,f_c_phylo,f_c_func))
reds <- colorRampPalette(brewer.pal(9,"Reds"))(100) 

#plot novel and disappearing assemblages
plot(c_f, col=rev(blues))  #novel
plot(f_c, col=rev(reds))   #disappearing








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

save.image("FutureAnalog.rData")
#s amnat
