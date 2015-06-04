# AlphaMapping.R ---------------------------------------------------------------
# Step 1) Bring in completed model data ----------------------------------------

# Bring in niche models, from the output directory specified above.
# get all the niche model data
niche <- list.files(out_path, pattern="ensemble.gri",full.name=T,recursive=T)

#split into current and future
#Get current models
current_niche <- niche[grep("current",niche,value=FALSE)]

#Get future models (emissions scenarios),check the SDM.R script
MICROC2070rcp26_niche<-niche[grep("MICROC2070rcp26",niche,value=FALSE)]
MICROC2070rcp85_niche<-niche[grep("MICROC2070rcp85",niche,value=FALSE)]
MICROC2070rcp45_niche<-niche[grep("MICROC2070rcp45",niche,value=FALSE)]

# create list of input rasters
input.niche<-list(current_niche, MICROC2070rcp26_niche, MICROC2070rcp45_niche,
                  MICROC2070rcp85_niche)
names(input.niche)<-c("current","MICROC2070rcp26","MICROC2070rcp45", "MICROC2070rcp85")


# Clip to Extent and shape of desired countries (Ecuador for now)
ec<-readOGR("InputData", "EcuadorCut")
r<-raster(extent(ec))

#Match cell size above from the SDM_SP function
res(r) <- cell_size
plot(ec.r <- rasterize(ec,r))

niche.crop <- lapply(niche,function(x){
  r <- crop(raster(x),extent(ec.r))
  filnam <- paste(strsplit(x,".gri")[[1]][1],"crop",sep="")
  writeRaster(r,filnam,overwrite=TRUE)
})

#get the crop files
niche.crops <- list.files(out_path,pattern="crop.gri",full.name=T,recursive=T)

#Get current models
current_niche <- niche.crops[grep("current",niche.crops,value=FALSE)]

#Get future models, for now its just
MICROC2070rcp26_niche <- niche.crops[grep("MICROC2070rcp26",niche.crops,value=FALSE)]
MICROC2070rcp45_niche <- niche.crops[grep("MICROC2070rcp45",niche.crops,value=FALSE)]
MICROC2070rcp85_niche <- niche.crops[grep("MICROC2070rcp85",niche.crops,value=FALSE)]

#create list of input rasters
input.niche <- list(current_niche, MICROC2070rcp26_niche, MICROC2070rcp45_niche, 
                    MICROC2070rcp85_niche)
names(input.niche) <- c("current","MICROC2070rcp26","MICROC2070rcp45", "MICROC2070rcp85")

#Create siteXspp table from input rasters, function is from AlphaMappingFunctions.R, sourced at the top. 
siteXspps <- lapply(input.niche, tableFromRaster, threshold=0.05) # **DECISION**
save(siteXspps, file = paste(out_path, "siteXspps.rda", sep = "/"))
     
# Step 2) Alpha Cell Statistics ------------------------------------------------
# find the taxonomic, phylgoenetic and trait diversity at each cell

# Bring in Phylogenetic Data
trx<-read.nexus("InputData/ColombiaPhylogenyUM.tre")
spnames<-read.table("InputData/SpNameTree.txt" , sep = "\t", header = TRUE)

# Replace tip.label with Spnames, replace the tiplabels with periods, which is
# the biomod default Cophenetic distance is the distance between all pairs of
# species- measure of relatedness
trx$tip.label <- gsub("_",".",as.character(spnames$SpName))
co<-cophenetic(trx)
save(trx, file=paste(out_path, "trx.rda", sep = "/"))

# Bring in morphology
# Bring in trait data
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
save(fco, file=paste(out_path, "fco.rda", sep = "/"))
# create a blank raster object of the correct size and extent to have for
# projecting the cell values
blank <- raster(niche.crops[[1]])
save(blank, file = paste(out_path, "blank_raster.rda", sep="/"))
# Step 2a) Phylogenetic Alpha Diversity (MPD) ----------------------------------

#Remove communities with less than 1 species in a row
#just get where diversity > 1, there is no phylogenetic diversity or functional
#diversity of species with richness = 1
MPDs <- lapply(siteXspps,AlphaPhylo) 

# Step 2b) Trait Alpha Tree Diversity (MFD) ------------------------------------

MFDs <- lapply(siteXspps,AlphaFunc) 

#Visualize Mapping Metrics
#set to the number of climate scenarios.
par(mfrow=c(4,3))

cell.Rasters <- lapply(names(siteXspps),cellVisuals) 
names(cell.Rasters) <- names(siteXspps)

# Step 3) Calculate differences among climate projections ---------------------

pdf(file=paste(out_path, "current_alpha.pdf", sep = "/"), height=4.2, width=7)
current<-cell.Rasters[[1]]
blues <- colorRampPalette(brewer.pal(9,"Blues"))(100)
plot(stack(current), col=blues)
dev.off()

#Compute Differences

diff.raster<-lapply(2:length(cell.Rasters),function(x){
  out<-stack(current[[1]]-cell.Rasters[[x]][[1]],
             current[[2]]-cell.Rasters[[x]][[2]],
             current[[3]]-cell.Rasters[[x]][[3]])
  names(out)<-c("Richness","Phylo","Func")
  return(out)}
)

names(diff.raster)<-names(siteXspps[-1])

#plot the differences between current and future alpha diversity for each of the scenarios (Richness, Phylogenetic, and Functional)
pdf(file=paste(out_path, "compare_alpha.pdf", sep="/"), height=8.5, width=11)
#par(mfrow=c(3,3))
cols <- colorRampPalette(brewer.pal(7,"RdBu"))(100)
plot(diff.raster[[1]], col=cols)
plot(diff.raster[[2]], col=cols)
plot(diff.raster[[3]], col=cols)
dev.off()

#I'm not entirely sure of the purpose of this section - need to find out.
#Leaving as is for now 

# Correlate rasters
al <- lapply(1:length(diff.raster), function(x){
  within.cor <- cor(values(diff.raster[[x]]), use="complete.obs")
  #within.cor <- data.frame(id = rownames(within.cor), within.cor)
  within.cor <- melt(within.cor)
  a <- qplot(data=within.cor, x=Var1, y=Var2, fill=value, geom="tile") + xlab("") + ylab("") + 
    scale_fill_continuous(low="white", high="red") + geom_text(aes(label=round(value,2))) + 
    theme(text=element_text(size=20)) + ggtitle(names(diff.raster[x]))
  return(a)
})

#Write difference in alpha out to file. 
#############
#This needs to be inspected, does this work for multiple climate scenarios, need to test?
#############
#Write alpha rasters to file
lapply(1:length(cell.Rasters),function(x){
  writeRaster(stack(cell.Rasters[[x]]), "/Figures/", names(cell.Rasters)[x],sep=""),
  overwrite=TRUE,bylayer=TRUE,suffix='names')
})

#Write difference raster to file
lapply(1:length(diff.raster),function(x){
  writeRaster(diff.raster[[x]], bylayer=TRUE, "Figures/AlphaChange.tif",
              overwrite=TRUE,suffix=names(diff.raster[[x]]))
})

# save the workspace - this gets picked back up by FutureAnalog.R
save.image(paste(out_path,"/AlphaMapping.rData",sep=""))
