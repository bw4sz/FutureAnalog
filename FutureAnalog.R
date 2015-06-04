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
source("FutureAnalogFunctions.R")

# Set output folder
# set the cell size for the analysis
cell_size = 0.1 

# create folders to output the models to
output_folder = "../FutureAnalog_output" 
out_path <- paste(output_folder, cell_size, sep = "/")

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

# create a blank raster object of the correct size and extent to have for
# projecting the cell values
blank <- raster(niche.crops[[1]])

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

# list to output results
NonAnalogRasters <- list()

# Need to write code to create mod.list - will probably be listing the files in
# the worldclim directory...

for(mod in mod.list){
  # PART I BETWEEN TIME BETA-DIVERSITY MEASURES --------------------------------
  # see FutureAnalogFunctions.R for how these are calculated
  fnBetaDiv(mod)

  # PART II: ANALOG ANALYSIS ---------------------------------------------------
  
  # This code has been set up so that additional lines with different thresholds
  # can be added in for sensitivity analysis
  
  # Step 1) CURRENT COMMUNITIES WITHOUT ANALOGS IN FUTURE (Disappearing) -------
  
  # How many current communities have do not have analogs in the future?
  # These are akin to communities which will disappear, "Disappearing"
  
  # For each of the current communities how many future communities fall below the
  # threshold (e.g., are analogous; 0 = similar, 1 = different)
  c_f_tax <- fnCurrent2Future(beta.time.taxa, "Taxon")
  c_f_phylo <- fnCurrent2Future(beta.time.phylo, "Phylo")
  c_f_func <- fnCurrent2Future(beta.time.func, "Func")
  
  # Create output raster stack (Disappearing) 
  c_f <- stack(c(c_f_tax,c_f_phylo,c_f_func))
  names(c_f) <- c("Taxonomic", "Phylogenetic", "Functional")
  
  # Step 2) FUTURE COMMUNITIES WITHOUT ANALOGS IN CURRENT (Novel) --------------
  
  #How many future communities do not have analogs in the current time?
  #These are akin to communities that are novel, "Novel"
  f_c_tax <- fnFuture2Current(beta.time.taxa, arbthresh)
  f_c_phylo <- fnFuture2Current(beta.time.phylo, arbthresh)
  f_c_func <- fnFuture2Current(beta.time.func, arbthresh)
  
  # Visualize all three NON-ANALOGS together --------------------------
  f_c <- stack(c(f_c_tax,f_c_phylo,f_c_func))
  names(f_c) <- c("Taxonomic", "Phylogenetic", "Functional")

  # output the raster stacks
  results <- stack(c_f, f_c)
  names(results) <- c(paste("Novel",c("Tax","Phylo","Func")), 
                      paste("Disappearing",c("Tax","Phylo","Func")))
  NonAnalogRasters[[mod]] <- results
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

save.image("FutureAnalog.rData")
#s amnat
