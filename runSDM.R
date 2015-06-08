# runSDM.R: --------------------------------------------------------------------

# Code to input species and climate data and then run the SDMs
# (ensemble and projections).

# load required packages (installing if not already done)
packages <- c("raster", "dplyr", "tidyr", "biomod2")

for(p in packages) {
  if (!p %in% installed.packages()) {
    install.packages(p)
  }
  require(p, character.only = TRUE)
}

# bring in functions (AlphaMappingFunctions and fnSDM)
#source("AlphaMappingFunctions.R")
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

# Step 1) Load presence data ---------------------------------------------------

# Lets go get the presence data on hummingbird distributions
PA <- read.csv("InputData/MASTER_POINTLOCALITYarcmap_review.csv")

# Just take the columns you want. 
PAdat <- select(PA, RECORD_ID, SPECIES, COUNTRY, LOCALITY, LATDECDEG, 
                LONGDECDEG, Decision, SpatialCheck, MapDecision)

#Just get the clean localities
gooddata <- c("ok", "Ok", "OK", "Y") #note blanks and REJECT data are excluded
loc_clean <- filter(PAdat, SpatialCheck=="Y", MapDecision %in% gooddata)

# Step 2) Load climate data ---------------------------------------------------- 

# The Paths to the climate layers must be changed. The layers are too large to 
# hang out on dropbox and github (40gb) Unzip the files to a local directory and
# change the paths. 

# Laura's note (change to above) - if the climate layers are in a folder called
# "worldclim_data" in the same level as the project folder, only the references
# to the individual climate data sets need changing here

# Import environmental data from worldclim, three variables Bio1 = annual mean
# temp, Bio12 = annual precip, Bio15 = precip seasonality
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

# Presence points are still taken from everywhere Avoid projecting into areas 
# where sample size is really low **DECISION**
exte<-extent(c(-81.13411,-68.92061,-5.532386,11.94902))

# Cut by the extent. Crop by this layer 
myExpl.crop<-stack(crop(myExpl,exte))

# Step 3) List the species to run models for -----------------------------------
spec <- table(loc_clean$SPECIES)
spec <- names(spec[which(spec >= 10)])

# Step 4) Run SDM_SP and project current ---------------------------------------
setwd(out_path)
for(x in 1:length(spec)) {
  if(!file.exists(gsub(" ", ".", spec[x]))) {
    SDM_SP(spec[x], loc_clean, myExpl)
  } 
  if(!file.exists(paste0(gsub(" ", ".", spec[x]), "/proj_current"))) {
    bio_project(spec[x], myExpl.crop, "current")
  }
}
setwd("../../FutureAnalog")

# Step 5) create evaluation plots ----------------------------------------------
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
