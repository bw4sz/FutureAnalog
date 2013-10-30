##Alpha Mapping
#Ben Weinstein Stony Brook University

##Goals: To develop parallel computing methods to turn niche models into assemblages and perform phylogenetic and functional betadiversity metrics
#This is the first script, and calls a source function to compute non-analog assemblages

#See Stralberg, D., D. Jongsomjit, C. a Howell, M. a Snyder, J. D. Alexander, J. a Wiens, and T. L. Root. 2009. Re-shuffling of species with climate disruption: a no-analog future for California birds? PloS One 4:e6825.

#October 16th 2012 - Ben Weinstein. Stony Brook University

#require packages
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

##
#######################################################################################################################
#Please note all paths must be changed, we are switching over to Github workflow, credit sarah for the push (no pun...)
#######################################################################################################################
##

#load workspace if needed on reset
#load("C:\\Users\\Ben\\Dropbox\\Lab paper 1 Predicted vs observed assemblages\\AlphaMapping.rData")

#source in all the Alpha Mapping functions
source("C:\\Users\\Ben\\Dropbox\\Scripts/Dimensions/AlphaMappingFunctions.R")

#Bring in Phylogenetic Data
trx<-read.nexus("C:\\Users\\Ben\\Dropbox\\Shared Ben and Catherine\\DimDivEntire\\\\Files for Analysis\\ColombiaPhylogenyUM.tre")
spnames<-read.table(file="C:\\Users\\Ben\\Dropbox\\Shared Ben and Catherine\\DimDivEntire\\Files for Analysis\\SpNameTree.txt" , sep = "\t", header = TRUE)

#Replace tip.label with Spnames#
#replace the tiplabels with periods, which is the biomod default
trx$tip.label<-gsub("_",".",as.character(spnames$SpName))
co<-cophenetic(trx)

#Bring in morphology
###Bring in trait data
morph <- read.csv("C:\\Users\\Ben\\Dropbox\\Lab paper 1 Predicted vs observed assemblages\\MorphologyShort.csv",na.strings="9999")

#just get males
morph.male<-morph[morph$Sex=="Macho",c("SpID","ExpC","Peso","AlCdo")]
morph.complete<-morph.male[complete.cases(morph.male),]

#aggregate for species
agg.morph<-aggregate(morph.complete,list(morph.complete$SpID),mean)
mon<-agg.morph[,-2]
colnames(mon)<-c("Species","Bill","Mass","WingChord")
rownames(mon)<-gsub(" ",".",mon[,1])
mon<-mon[,-1]

#principal component traits and get euclidean distance matrix
fco<-as.matrix(dist(prcomp(mon)$x))

#bring in traits
morph <- read.csv("C:\\Users\\Ben\\Dropbox\\Lab paper 1 Predicted vs observed assemblages\\MorphologyShort.csv",na.strings="9999")

#just get males
morph.male<-morph[morph$Sex=="Macho",c("SpID","ExpC","Peso","AlCdo")]
morph.complete<-morph.male[complete.cases(morph.male),]

#aggregate for species
agg.morph<-aggregate(morph.complete,list(morph.complete$SpID),mean)
mon<-agg.morph[,-2]
colnames(mon)<-c("Species","Bill","Mass","Wing Chord")


####Bring in niche models
  setwd("D:/Niche_Models/")

##############################
#Bring in Model Data
##############################

#get all the niche model data
niche<-list.files(getwd(),pattern="TotalConsensus_EMbyROC.gri",full.name=T,recursive=T)

#split into current and future
#Get current models
current_niche<-niche[grep("current",niche,value=FALSE)]

#Get future models, for now its just
MICROC2070rcp26_niche<-niche[grep("MICROC2070rcp26",niche,value=FALSE)]
MICROC2070rcp85_niche<-niche[grep("MICROC2070rcp85",niche,value=FALSE)]
MICROC2070rcp45_niche<-niche[grep("MICROC2070rcp45",niche,value=FALSE)]

#create list of input rasters
input.niche<-list(current_niche,MICROC2070rcp26_niche,MICROC2070rcp45_niche,MICROC2070rcp85_niche)

names(input.niche)<-c("current","MICROC2070rcp26","MICROC2070rcp45", "MICROC2070rcp85")


##############################
#Optional Clip by extent, skip to line 101
##############################
ec<-readShapePoly("C:/Users/Ben/Dropbox/Shared Ben and Catherine/FutureAnalog/EcuadorCut.shp")
r<-raster(extent(ec))
res(r)<-.1
plot(ec.r<-rasterize(ec,r))
niche.crop<-lapply(niche,function(x){
  r<-crop(raster(x),extent(ec.r))
filnam<-paste(strsplit(x,".gri")[[1]][1],"crop",sep="")
  writeRaster(r,filnam,overwrite=TRUE)
  })

#get the crop files
niche.crops<-list.files(getwd(),pattern="TotalConsensus_EMbyROCcrop.gri",full.name=T,recursive=T)

#Get current models
current_niche<-niche.crops[grep("current",niche.crops,value=FALSE)]

#Get future models, for now its just
MICROC2070rcp26_niche<-niche.crops[grep("MICROC2070rcp26",niche.crops,value=FALSE)]
MICROC2070rcp85_niche<-niche.crops[grep("MICROC2070rcp85",niche.crops,value=FALSE)]
MICROC2070rcp45_niche<-niche.crops[grep("MICROC2070rcp45",niche.crops,value=FALSE)]

#create list of input rasters
input.niche<-list(current_niche,MICROC2070rcp26_niche,MICROC2070rcp45_niche,MICROC2070rcp85_niche)

names(input.niche)<-c("current","MICROC2070rcp26","MICROC2070rcp45", "MICROC2070rcp85")
#################################################### Start here if not cropping!

#Create siteXspp table from input rasters
siteXspps<-lapply(input.niche,tableFromRaster,threshold=700)

#Get future models, for now its just
MICROC2070rcp26_niche<-niche[grep("MICROC2070rcp26",niche,value=FALSE)]
MICROC2070rcp85_niche<-niche[grep("MICROC2070rcp85",niche,value=FALSE)]
MICROC2070rcp45_niche<-niche[grep("MICROC2070rcp45",niche,value=FALSE)]

#create list of input rasters
input.niche<-list(current_niche,MICROC2070rcp26_niche,MICROC2070rcp45_niche,MICROC2070rcp85_niche)

names(input.niche)<-c("current","MICROC2070rcp26","MICROC2070rcp45", "MICROC2070rcp85")


#################################################
#Alpha Cell Statistics - fid the tasxonomic, phylgoenetic and trait diversity at each cell
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

MFDs<-lapply(siteXspps,AlphaFunc.tree)

###################################################################
#Functional Dispersion Metric Lalibert?, E., and P. Legendre (2010)
##################################################################

FDs<-lapply(siteXspps,AlphaFunc.FD,traits=mon)

##########################
#Visualize Mapping Metrics
##########################
par(mfrow=c(2,3))

cellVisuals<-function(inp.name){

  #Test with richness
richness<-cellVis(cells=rownames(siteXspps[names(siteXspps) %in% inp.name][[1]]),value=apply(siteXspps[names(siteXspps) %in% inp.name][[1]],1,sum))

#Phylogenetic

MPD.vis<-cellVis(cells=MPDs[names(MPDs) %in% inp.name][[1]]$Cell,value=MPDs[names(MPDs) %in% inp.name][[1]]$MPD)

#Func Tree
MFD.vis<-cellVis(MFDs[names(MFDs) %in% inp.name][[1]]$Cell,MFDs[names(MFDs) %in% inp.name][[1]]$MFD)
return(list(richness,MPD.vis,MFD.vis))}

cell.Rasters<-lapply(names(siteXspps),cellVisuals)

##FD Functional Divergence
#apply visualization to all columns
FD_metrics<-lapply(FDs,function(x){
  apply(x,2,cellVis,cells=rownames(x))
})
  
  
#create a stack of metrics

for ( x in 1:length(cell.Rasters)){
  cell.Rasters[[x]][[4]]<-stack(FD_metrics[[x]][-c(1,2,4,8)])
}

#Normalize the metrics
FD.norm<-FD.stack/cellStats(FD.stack,"max")

plot(FD.norm)

################################################
#Calculate differences among climate projections
################################################

#
current<-cell.Rasters[[1]]

diff.raster<-lapply(2:length(cell.Rasters),function(x){
  out<-list(current[[1]]-cell.Rasters[[x]][[1]],
            current[[2]]-cell.Rasters[[x]][[2]],
            current[[3]]-cell.Rasters[[x]][[3]],
current[[4]]-cell.Rasters[[x]][[4]])
  names(out)<-c("Richness","Phylo","Func")
  return(stack(out))}
  )

names(diff.raster)<-names(siteXspps[-1])

par(mfrow=c(3,2))
plot(diff.raster[[1]],col=brewer.pal(7,"RdBu"))
plot(diff.raster[[2]],col=brewer.pal(7,"RdBu"))

#####################
#Correlate rasters
al<-lapply(1:length(diff.raster),function(x){
within.cor<-cor(values(diff.raster[[x]]),use="complete.obs")
within.cor<-melt(within.cor)
a<-qplot(data=within.cor,x=X1,y=X2,fill=value,geom="tile") + xlab("") + ylab("")+ scale_fill_continuous(low="blue",high="red") + geom_text(aes(label=round(value,2)))
return(a)
})

#Write difference raster to file
writeRaster(diff.raster[[2]],bylayer=TRUE,"C:\\Users\\Ben\\Dropbox\\Shared Ben and Catherine\\FutureAnalog\\Alpha\\AlphaChange.tif",overwrite=TRUE,suffix=names(diff.raster[[2]]))

save.image("C:\\Users\\Ben\\Dropbox\\Lab paper 1 Predicted vs observed assemblages\\AlphaMapping.rData")
