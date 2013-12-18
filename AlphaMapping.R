##Alpha Mapping
#Ben Weinstein Stony Brook University

##Goals: To develop parallel computing methods to turn niche models into assemblages and perform phylogenetic and functional betadiversity metrics
#This is the first script, and calls a source function to compute non-analog assemblages

#See Stralberg, D., D. Jongsomjit, C. a Howell, M. a Snyder, J. D. Alexander, J. a Wiens, and T. L. Root. 2009. Re-shuffling of species with climate disruption: a no-analog future for California birds? PloS One 4:e6825.

#October 16th 2012 - Ben Weinstein. Stony Brook University

#require packages for alpha analysis
require(RColorBrewer)
require(biomod2)
require(maptools)
require(ggplot2)
require(reshape)
require(raster)
require(rgdal)
require(doSNOW)
require(ape)
require(stringr)
require(FD)
require(picante)
require(parallel)
require(ecodist)
require(vegan)
require(MASS)


#Set dropbox and github paths
# Ben's
# droppath<-"C:\\Users\\Ben\\Dropbox\\"
# gitpath<-"C:\\Users\\Ben\\Documents\\FutureAnalog\\"

# Anusha's paths
droppath <- "C:\\Users/Anusha/Documents/Dropbox/Hummingbirds/Lab paper 1 Predicted vs observed assemblages/"
gitpath <- "C:\\Users/Anusha/Documents/GitHub/FutureAnalog/"
setwd(droppath)
#######################################################################################################################
#Please note all paths must be changed, we are switching over to Github workflow, credit sarah for the push (no pun...)
#######################################################################################################################
##

#load workspace if needed on reset
load("C:\\Users\\Ben\\Dropbox\\Lab paper 1 Predicted vs observed assemblages\\AlphaMapping.rData")

#source in all the Alpha Mapping functions
source(paste(gitpath,"AlphaMappingFunctions.R",sep=""))

#Bring in Phylogenetic Data
trx<-read.nexus(paste(gitpath,"InputData/ColombiaPhylogenyUM.tre",sep=""))
spnames<-read.table(paste(gitpath,"InputData/SpNameTree.txt",sep="") , sep = "\t", header = TRUE)

#Replace tip.label with Spnames#
#replace the tiplabels with periods, which is the biomod default
# Cophenetic distance is the distance between all pairs of species- measure of relatedness
trx$tip.label<-gsub("_",".",as.character(spnames$SpName))
co<-cophenetic(trx)

#Bring in morphology
###Bring in trait data
morph <- read.csv(paste(gitpath,"InputData/MorphologyShort.csv",sep=""),na.strings="9999")

#just get males & the 3 traits of interest
morph.male<-morph[morph$Sex=="Macho",c("SpID","ExpC","Peso","AlCdo")]
morph.complete<-morph.male[complete.cases(morph.male),]

#aggregate for species
agg.morph<-aggregate(morph.complete,list(morph.complete$SpID),mean)
mon<-agg.morph[,-2]
colnames(mon)<-c("Species","Bill","Mass","WingChord")
rownames(mon)<-gsub(" ",".",mon[,1])
mon<-mon[,-1]

#principal component traits and get euclidean distance matrix
means <- apply(mon, 2, mean)

Bill <- mon$Bill - means["Bill"]/sd(mon$Bill)
Mass <- mon$Mass - means["Mass"]/sd(mon$Mass)
WingChord <- (mon$WingChord - means["WingChord"])/sd(mon$WingChord)

z.scores <- data.frame(Bill, Mass, WingChord)
rownames(z.scores) <- rownames(mon)

fco <- as.matrix(dist(z.scores, method = "euclidean"))

##################################################
#Niche models for each species need to be run
#Beware Biomod Outputs are HUGE, i'd suggest starting with a couple species and saving locally
#########################################################################################

#source SDM function
source(paste(gitpath,"SDM.r",sep=""))

#Function computes ensemble niche models with default methods (Models=GLM,GBM,MAXENT) and paramters (background draws, keep ROC > .75)
#The SDM function could be tweaked to take in any set of parameters, but given the enormous number of options, start simple. 
#current arguements cell size (in degrees, ie 1 degree = 112km^2 at the equator) and output directory (If it doesn't exist, script will create it)

########################################################
#Before you run the function, go through this checklist
########################################################
#1.Set the Extent (Line 95)
#The extent was the bounding box we used for another paper and could be set a variable, or changed internally if you like

#2. Run for how many species
#i'd suggest just running on a few species at first, go to line 159 (x=1:length(spec)) and change length(spec) to the total number of species desired, rather than all species

#3 Set the climate scenarios (line 58 for current, 91 for future)

#There really isn't an elegant way of coding which climate scenerios you want, so please pay close attention to lines 
#The climate scenerios are very large, and need to be held locally

#In general the script is meant as a helpful wrapper, but attention needs to be paid to alot of the manual details in the setup of BIOMOD, there are simply too many options and dependencies for me to take a all encompassing function

####################################
#Perform Niche Models
####################################
#Define these variables outside the function so they can be used below.
# Cell size is in degrees. 1 degree = 112km
cell_size=.75
#output_folder="C:/Users/Ben/Desktop/Testmod"
output_folder = "C:/Users/Anusha/Desktop/Testmod"
SDM_SP(cell_size,output_folder)

##############################
#Bring in Completed Model Data
##############################

####Bring in niche models, from the output directory specified above.
#get all the niche model data
## TO DO------------ package update changed naming structure ("TotalConsensus.gri"?)
niche<-list.files(output_folder,pattern="TotalConsensus_EMbyROC.gri",full.name=T,recursive=T)

#split into current and future
#Get current models
current_niche<-niche[grep("current",niche,value=FALSE)]

#Get future models, for now its just one, check the SDM.R script
MICROC2070rcp26_niche<-niche[grep("MICROC2070rcp26",niche,value=FALSE)]
#MICROC2070rcp85_niche<-niche[grep("MICROC2070rcp85",niche,value=FALSE)]
#MICROC2070rcp45_niche<-niche[grep("MICROC2070rcp45",niche,value=FALSE)]

#create list of input rasters
#input.niche<-list(current_niche,MICROC2070rcp26_niche,MICROC2070rcp45_niche,MICROC2070rcp85_niche)

input.niche<-list(current_niche,MICROC2070rcp26_niche)
#names(input.niche)<-c("current","MICROC2070rcp26","MICROC2070rcp45", "MICROC2070rcp85")
names(input.niche)<-c("current","MICROC2070rcp26")


###############################################
#Clip to Extent and shape of desired countries (Ecuador for now)
###############################################

ec<-readShapePoly(paste(gitpath,"Inputdata/EcuadorCut.shp",sep=""))
r<-raster(extent(ec))

#Match cell size above from the SDM_SP function
res(r)<-cell_size
plot(ec.r<-rasterize(ec,r))
niche.crop<-lapply(niche,function(x){
  r<-crop(raster(x),extent(ec.r))
filnam<-paste(strsplit(x,".gri")[[1]][1],"crop",sep="")
  writeRaster(r,filnam,overwrite=TRUE)
  })

#get the crop files
niche.crops<-list.files(paste(output_folder,cell_size,sep="/"),pattern="TotalConsensus_EMbyROCcrop.gri",full.name=T,recursive=T)

#Get current models
current_niche<-niche.crops[grep("current",niche.crops,value=FALSE)]

#Get future models, for now its just
MICROC2070rcp26_niche<-niche.crops[grep("MICROC2070rcp26",niche.crops,value=FALSE)]
#MICROC2070rcp85_niche<-niche.crops[grep("MICROC2070rcp85",niche.crops,value=FALSE)]
#MICROC2070rcp45_niche<-niche.crops[grep("MICROC2070rcp45",niche.crops,value=FALSE)]

#create list of input rasters
#input.niche<-list(current_niche,MICROC2070rcp26_niche,MICROC2070rcp45_niche,MICROC2070rcp85_niche)
#names(input.niche)<-c("current","MICROC2070rcp26","MICROC2070rcp45", "MICROC2070rcp85")

input.niche<-list(current_niche,MICROC2070rcp26_niche)
names(input.niche)<-c("current","MICROC2070rcp26")

#Create siteXspp table from input rasters, function is from AlphaMappingFunctions.R, sourced at the top. 
siteXspps<-lapply(input.niche,tableFromRaster,threshold=700)

####################################################
#Niche Models Completed!
####################################################

#################################################
#Alpha Cell Statistics - find the taxonomic, phylogenetic and trait diversity at each cell
#################################################
#Functions are sourced from AlphaMappingFunctions.R

###create a blank raster object of the correct size and extent to have for projecting the cell values
blank<-raster(niche.crops[[1]])

#################################
#########################################
#######Phylogenetic Alpha Diversity (MPD)
#########################################

#Remove communities with less than 1 species in a row
#just get where diversity > 1, there is no phylgoenetic diversity or functional diversity of species with richness = 1

MPDs<-lapply(siteXspps,AlphaPhylo)

########################################
########## Trait Alpha Tree Diversity (MFD)
########################################

MFDs<-lapply(siteXspps,AlphaFunc)

##########################
#Visualize Mapping Metrics
##########################
#set to the number of climate scenerios.
par(mfrow=c(2,2))

cellVisuals<-function(inp.name){

  #Taxonomic richness
  richness<-cellVis(cells=rownames(siteXspps[names(siteXspps) %in% inp.name][[1]]),value=apply(siteXspps[names(siteXspps) %in% inp.name][[1]],1,sum))

  #Phylogenetic richness

  MPD.vis<-cellVis(cells=MPDs[names(MPDs) %in% inp.name][[1]]$Cell,value=MPDs[names(MPDs) %in% inp.name][[1]]$MPD)

  #Func Tree
  MFD.vis<-cellVis(MFDs[names(MFDs) %in% inp.name][[1]]$Cell,MFDs[names(MFDs) %in% inp.name][[1]]$MFD)

  out<-list(richness,MPD.vis,MFD.vis)
  names(out)<-c("Richness","Phylogenetic","Trait")
  return(out)
}

cell.Rasters<-lapply(names(siteXspps),cellVisuals)
names(cell.Rasters)<-names(siteXspps)

################################################
#Calculate differences among climate projections
################################################

## Assumes that it calls the first current in the list. Might want to change [[1]] to current
current<-cell.Rasters[[1]]

plot(stack(current))

#Compute Differences

diff.raster<-lapply(2:length(cell.Rasters),function(x){
  out<-stack(current[[1]]-cell.Rasters[[x]][[1]],
            current[[2]]-cell.Rasters[[x]][[2]],
            current[[3]]-cell.Rasters[[x]][[3]])
  names(out)<-c("Richness","Phylo","Func")
  return(out)
})

names(diff.raster)<-names(siteXspps[-1])

par(mfrow=c(3,2))
plot(diff.raster[[1]],col=brewer.pal(7,"RdBu"))
#plot(diff.raster[[2]],col=brewer.pal(7,"RdBu"))

#####################
#Correlate rasters
al<-lapply(1:length(diff.raster),function(x){
  within.cor<-cor(values(diff.raster[[x]]),use="complete.obs")
  within.cor<-melt(within.cor)
  a<-qplot(data=within.cor,x=X1,y=X2,fill=value,geom="tile") + xlab("") + ylab("")+ scale_fill_continuous(low="blue",high="red") + geom_text(aes(label=round(value,2)))
  return(a)
})

#Write difference in alpha out to file.

#############
#This needs to be inspected, does this work for multiple climate scenerios, need to test?
#############
#Write alpha rasters to file
lapply(1:length(cell.Rasters),function(x){
  writeRaster(stack(cell.Rasters[[x]]),paste(paste(gitpath,"Figures/",sep=""),names(cell.Rasters)[x],sep=""),overwrite=TRUE,bylayer=TRUE,suffix='names')
})

#Write difference raster to file
lapply(1:length(diff.raster),function(x){
  writeRaster(diff.raster[[x]],bylayer=TRUE,paste(gitpath,"Figures/AlphaChange.tif",sep=""),overwrite=TRUE,suffix=names(diff.raster[[x]]))
})

save.image(paste(droppath,"NASA_Anusha\\AlphaMapping.rData",sep=""))
