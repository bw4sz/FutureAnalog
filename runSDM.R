# runSDM.R: --------------------------------------------------------------------

# Code to input species and climate data and then run the SDMs
# (ensemble and projections).

# load required packages (installing if not already done)
packages <- c("raster", "dplyr", "tidyr", "doSNOW", "stringr")

for(p in packages) {
  if (!p %in% installed.packages()) {
    install.packages(p)
  }
  require(p, character.only = TRUE)
}

# bring in functions (AlphaMappingFunctions and fnSDM)
source("AlphaMappingFunctions.R")
source("fnSDM.R")

# set the cell size for the analysis - **DECISION**
cell_size = 0.1 

# create folders to output the models to
output_folder = "../FutureAnalog_output" 
out_path <- paste(output_folder, cell_size, sep = "/")
if(!dir.exists(output_folder)) dir.create(output_folder)
if(!dir.exists(out_path)) dir.create(out_path)
if(!dir.exists(paste(out_path, "logs", sep = "/"))) 
  dir.create(paste(out_path, "logs", sep = "/"))

# Step 1) Load presence data ###################################################

# Lets go get the presence data on hummingbird distributions
PA <- read.csv("InputData/MASTER_POINTLOCALITYarcmap_review.csv")

# Just take the columns you want. 
PAdat <- select(PA, RECORD_ID, SPECIES, COUNTRY, LOCALITY, LATDECDEG, 
                LONGDECDEG, Decision, SpatialCheck, MapDecision)

#Just get the clean localities
gooddata <- c("ok", "Ok", "OK", "Y") #note blanks and REJECT data are excluded
loc_clean <- filter(PAdat, SpatialCheck=="Y", MapDecision %in% gooddata)

# Step 2) Load climate data #################################################### 

# The Paths to the climate layers must be changed. The layers are too large to
# hang out on dropbox and github (40gb) Unzip the files to a local directory and
# change the paths.

# Import environmental data from worldclim, three variables
# Bio1 = annual mean temp, Bio12 = annual precip, Bio15 = precip seasonality
myExpl <- c("../worldclim_data/bio1-9_30s_bil/bio_1.bil",
            "../worldclim_data/bio10-19_30s_bil/bio_12.bil",
            "../worldclim_data/bio10-19_30s_bil/bio_15.bil")

myExpl <- stack(myExpl)

# get extent of the presence data
extPoint <- SpatialPoints(loc_clean[,c("LONGDECDEG","LATDECDEG")])

# crop by this extent
myExpl <- crop(myExpl,extPoint)

# set cell size
fact <- cell_size/res(myExpl) # the "factor" to aggregate by

# Set cell size to ~ cell_size degree
myExpl <- aggregate(myExpl,fact)

# Step 3) Climate Scenarios and Futute Climate ################################# 

# Bring in future climate layers (change here to add in more scenarios - all
# below code will need adapting to include new scenarios) 
# Modelname_year_emmissionscenario
MICROC_2070_rcp26 <- stack("../worldclim_data/mc26bi70/mc26bi701.tif",
                           "../worldclim_data/mc26bi70/mc26bi7012.tif",
                           "../worldclim_data/mc26bi70/mc26bi7015.tif")
MICROC_2070_rcp85 <- stack("../worldclim_data/mc85bi70/mc85bi701.tif",
                           "../worldclim_data/mc85bi70/mc85bi7012.tif",
                           "../worldclim_data/mc85bi70/mc85bi7015.tif")
MICROC_2070_rcp45 <- stack("../worldclim_data/mc45bi70/mc45bi701.tif",
                           "../worldclim_data/mc45bi70/mc45bi7012.tif",
                           "../worldclim_data/mc45bi70/mc45bi7015.tif")

# make the names consistent for all 
names(MICROC_2070_rcp26) <- names(myExpl)
names(MICROC_2070_rcp85) <- names(myExpl)
names(MICROC_2070_rcp45) <- names(myExpl)

# Step 4) Set the Extent to project *into*. ####################################

# Presence points are still taken from everywhere Avoid projecting into areas 
# where sample size is really low **DECISION**
exte<-extent(c(-81.13411,-68.92061,-5.532386,11.94902))

# Cut by the extent. Crop by this layer 
myExpl.crop<-stack(crop(myExpl,exte))
MICROC_2070_rcp26.c<-stack(crop(MICROC_2070_rcp26,exte))
MICROC_2070_rcp85.c<-stack(crop(MICROC_2070_rcp85,exte))
MICROC_2070_rcp45.c<-stack(crop(MICROC_2070_rcp45,exte))

# set resolution for the future layers equivalent to cell size
if (!res(MICROC_2070_rcp26)[[1]] == cell_size){
  fact1<-cell_size/res(MICROC_2070_rcp26.c)
  fact2<-cell_size/res(MICROC_2070_rcp85.c)
  fact3<-cell_size/res(MICROC_2070_rcp45.c)
  
  # Set cell size to ~ cell_size degree
  MICROC_2070_rcp26.c<-stack(aggregate(MICROC_2070_rcp26.c,fact1))
  MICROC_2070_rcp85.c<-stack(aggregate(MICROC_2070_rcp85.c,fact2))
  MICROC_2070_rcp45.c<-stack(aggregate(MICROC_2070_rcp45.c,fact3))
}

# create a list of all env to project into
projEnv<-list(myExpl.crop,MICROC_2070_rcp26.c,MICROC_2070_rcp45.c,MICROC_2070_rcp85.c)
names(projEnv)<-c("current","MICROC2070rcp26","MICROC2070rcp45","MICROC2070rcp85")

# Step 5) List the species to run models for (needs updating) ##################

# See code from the original SDM.R file for info on how it was done

# for now just list all of the species with enough points
spec <- table(loc_clean$SPECIES)
spec <- names(spec[which(spec >= 10)])

# Step 6) Run SDM_SP ###########################################################
setwd(out_path)
for(x in 4:10) {
  SDM_SP(spec[x], loc_clean, myExpl, projEnv, out_path)
}
setwd("../../FutureAnalog")

# Step 7) create evaluation plots ######################################################
#Get the model evaluation from file
model_eval<-list.files(out_path, full.name=TRUE,recursive=T,pattern="Eval.csv")
model_eval<-rbind_all(lapply(model_eval, 
                             function(x) read.csv(x, stringsAsFactors = FALSE)))
colnames(model_eval)[1:2] <- c("Stat", "Species")
model_eval  <- gather(data.frame(model_eval), Model, value, -Stat, -Species)
colnames(model_eval)<-c("Stat", "Species", "Model", "Value")

# heatmap of model evaluations
ggplot(model_eval, aes(x=Species,y=Model,fill=Value)) + 
  geom_tile() + 
  facet_wrap(~Stat, ncol = 1) +
  scale_fill_gradient("ROC",limits=c(0,1),low="blue",high="red",na.value="white") + 
  theme_classic() +
  theme(axis.text.x=element_text(angle=-90))
ggsave(paste(out_path, "ModelEvaluations.jpeg", sep = "/"), 
       dpi=600, height = 6, width=11)

#Plot correlation of ROC and TSS scores
model_compare <- spread(model_eval, Stat, Value)
ggplot(model_compare, aes(TSS, ROC)) + 
  geom_point() + 
  stat_smooth(method="lm") + 
  theme_classic() + 
  theme(text=element_text(size=20))
lm1=lm(ROC~TSS, data=model_compare)
ggsave(paste(out_path, "ModelComparison_ROC-TSS.jpeg", sep = "/"), 
       dpi=600, height = 6, width=11)

model_thresh<-sapply(seq(.5,.95,.05),function(x){
  table(model_eval$Value > x,model_eval$Model)["TRUE",]
})

# Model thresholds plot - number of species where ROC/TSS scores are above
# varying thresholds
model_thresh <- data.frame(Model = rownames(model_thresh), model_thresh)
colnames(model_thresh)[-1] <- seq(.5,.95,.05)
model_thresh <- gather(model_thresh, variable, value, -Model)

names(model_thresh)<-c("Model","ROC_Threshold","Number_of_Species")
ggplot(model_thresh,aes(x=ROC_Threshold,y=Number_of_Species,col=Model)) + 
  geom_line() + 
  geom_point() + 
  geom_text(aes(label=Number_of_Species),vjust=4,size=5) + 
  xlab("Model Threshold") + 
  ylab("Number of species included") +
  theme_classic() + theme(text=element_text(size=20))
ggsave(paste(out_path, "ModelThresholding.jpeg", sep = "/"), 
       dpi=600, height=8, width=8)

#Get the variable importance from file
varI <- list.files(paste(out_path, sep = "/"), 
                   full.name=TRUE,recursive=T,pattern="VarImportance.csv")
varI<-rbind_all(lapply(varI, function(x) read.csv(x, stringsAsFactors = FALSE)))
varI <- data.frame(varI[,-1])

#Melt variable for plotting
mvar <- gather(varI, variable, value, -X1, -spec)
colnames(mvar)<-c("Bioclim","Species","Model","value")

#Plot variable importance across all models
ggplot(mvar, aes(x=Species,y=Bioclim,fill=value)) + 
  geom_tile() + 
  scale_fill_gradient(limits=c(0,1),low="blue",high="red",na.value="white") + 
  theme_classic() + 
  theme(axis.text.x=element_text(angle=-90)) + 
  facet_grid(Model ~ .)

ggsave(paste(out_path, "VariableImportance.jpeg", sep = "/"), 
       dpi=600, height = 6, width=11)

# Step 8) Bring in completed model data ########################################

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
siteXspps <- lapply(input.niche, tableFromRaster, threshold=0.05)
