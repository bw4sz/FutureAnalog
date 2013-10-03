#########Biomod2 - Ben Weinstein, Stony Brook University 10/11/2012

#Install packages - only needs to be done the first run
#install.packages("biomod2", repos = "http://R-Forge.R-project.org", dependencies = TRUE)
#install.packages("gam")


#Wrap this into a function to be called from another script

SDM_SP<-function(cell_size,output_folder){

#If you have already installed, let's start here.
#Call the packages we are going to need in this tutorial
require(biomod2)
require(maptools)
require(ggplot2)
require(reshape)
require(raster)
require(rgdal)
require(doSNOW)
require(stringr)

#set a working directory, where do we want to save files
#save locally for now
setwd("output_folder")
dir.create(paste(getwd(),cell_size,sep="/"))
setwd(paste(getwd(),cell_size,sep="/"))
dir.create("logs")

#To perform the biomod, you must have three pieces of data
#1) Presence Absence Matrix
#2) Input localities in a 2 column matrix
#3) Environmental Variables - unclear whether this should be masked or not

#################################

#1) presence absence data matrix for the desired species
#For the graham lab, or users of the dropbox, just change the user below to your name, or the whatever the dropbox path is:"

#Lets go get the presence data on hummingbird distributions
PA<-read.csv("C:\\Users\\Ben\\Dropbox\\Lab paper 1 Predicted vs observed assemblages\\MASTER_POINTLOCALITYarcmap_review.csv")
#Just take the columns you want. 
PAdat<-PA[,colnames (PA) %in% c("RECORD_ID","SPECIES","COUNTRY","LOCALITY","LATDECDEG","LONGDECDEG","Decision","SpatialCheck","MapDecision")]

PAdat<-PAdat[!PAdat$LONGDECDEG==-6,]
#We are going to use the layers that juan used for his Amnat paper. Only on laptop
#myExpl <- c("C:\\Users\\Ben\\Documents\\GIS\\bio_5m_bil/bio3.bil", "C:\\Users\\Ben\\Documents\\GIS\\bio_5m_bil/bio4.bil","C:\\Users\\Ben\\Documents\\GIS\\bio_5m_bil/bio7.bil","C:\\Users\\Ben\\Documents\\GIS\\bio_5m_bil/bio11.bil","C:\\Users\\Ben\\Documents\\GIS\\bio_5m_bil/bio12.bil" )

##############################
#The Paths to the climate layers must be changed. The layers are too large to hang out on dropbox and github (40gb)
#Unzip the files to a local directory and change the paths.

#myExpl <- c("D:\\worldclim_bio1-9_30s_bil\\bio_1.bil",
            #"D:\\worldclim_bio1-9_30s_bil\\bio_2.bil",
            ##"D:\\worldclim_bio1-9_30s_bil\\bio_4.bil",
            #"D:\\worldclim_bio1-9_30s_bil\\bio_12.bil",
            #"D:\\worldclim_bio1-9_30s_bil\\bio_15.bil")

#Paths must be changed to local directory! unzip
#Import environmental data from worldclim, three variables
myExpl <- c("D:\\worldclim_bio1-9_30s_bil\\bio_1.bil",
            "D:\\worldclim_bio1-9_30s_bil\\bio_12.bil",
            "D:\\worldclim_bio1-9_30s_bil\\bio_15.bil")

myExpl<-stack(myExpl)


#Just get the clean localities
loc_clean<-PAdat[PAdat$Decision=="OK"|PAdat$MapDecision=="OK"|PAdat$MapDecision=="",]
loc_clean<-PAdat[PAdat$SpatialCheck=="Y",]
extPoint<-SpatialPoints(loc_clean[,c("LONGDECDEG","LATDECDEG")])
#exte<-extent(c(-81.13411,-68.92061,-5.532386,11.94902))

#Crop by this layer, 
myExpl<-crop(myExpl,extPoint)
res(myExpl)

#Set Cell size
####################################
#set layer sie

fact<-cell_size/res(myExpl)
####################################

#Set cell size to ~ .1 degree
myExpl<-aggregate(myExpl,fact)

#Bring in future climate layers (need to uncomment here when ready, anusha, and subsequent lines)
#Start with one layer
#MICROC_2070_rcp26<-stack("D:\\Future GCM Layers\\MICROC 2070\\MICROCrcp26\\biovars.grd")[[c(1,12,15)]]
#MICROC_2070_rcp85<-stack("D:\\Future GCM Layers\\MICROC 2070\\MICROCrcp85\\biovars.grd")[[c(1,12,15)]]
#MICROC_2070_rcp45<-stack("D:\\Future GCM Layers\\MICROC 2070\\MICROCrcp45\\biovars.grd")[[c(1,12,15)]]


exte<-extent(c(-81.13411,-68.92061,-5.532386,11.94902))

#Cut by the extent

#Crop by this layer, 
myExpl.crop<-stack(crop(myExpl,exte))
#MICROC_2070_rcp26.c<-stack(crop(MICROC_2070_rcp26,exte))
#MICROC_2070_rcp85.c<-stack(crop(MICROC_2070_rcp85,exte))
#MICROC_2070_rcp45.c<-stack(crop(MICROC_2070_rcp85,exte))

#names(MICROC_2070_rcp26.c)<-names(myExpl)
#names(MICROC_2070_rcp85.c)<-names(myExpl)
#names(MICROC_2070_rcp45.c)<-names(myExpl)

#create a list of all env to project into
#projEnv<-list(myExpl.crop,MICROC_2070_rcp26.c,MICROC_2070_rcp45.c,MICROC_2070_rcp85.c)
projEnv<-list(myExpl.crop)

#names(projEnv)<-c("current","MICROC2070rcp26","MICROC2070rcp45","MICROC2070rcp85")
names(projEnv)<-c("current")

#######################
#Which species to run
#######################

#How many records per species?
rec<-table(loc_clean$SPECIES)

#let's grab some species that have been checked 
spec<-names(rec[which(rec >= 10)])

#remove any species that have already been run
#name the list with the correct species names from file
#get all the niche model data
niche<-list.files(getwd(),pattern="TotalConsensus_EMbyROC.gri",full.name=T,recursive=T)

#split into current and future
#Get current models
completed_models<-lapply(names(projEnv), function(x){
run_mod<-niche[grep(x,niche,value=FALSE)]
if(length(run_mod)==0) return(NA)
completed<-str_match(run_mod,pattern="Models/(\\w+.\\w+)")[,2]
})

names(completed_models)<-names(projEnv)

#Which species have been run in which GCM models?
modelXspp<-sapply(spec,function(x){
  sapply(completed_models, function(y){
    gsub(" ",".",x) %in% y
  })
})

#rownames(modelXspp)<-names(projEnv)

#Only run a species if it is NOT run in all models
if(length(projEnv)==1){spec<-spec[!(modelXspp*1)==length(projEnv)]}
if(!length(projEnv)==1){spec<-spec[!apply(modelXspp,2,sum)==length(projEnv)]}


paste("Species to be modeled",spec,sep=": ")

cl<-makeCluster(8,"SOCK")
registerDoSNOW(cl)
system.time(niche_loop<-foreach(x=1:length(spec),.packages=c("reshape","biomod2"),.verbose=T,.errorhandling="pass") %dopar% {
  sink(paste("logs/",paste(spec[[x]],".txt",sep=""),sep=""))
  
  #remove sites that have no valid records
  #For the moment, only get the clean records from Decision, or the cleaned map localities. 
  print(paste("Start Time is",Sys.time()))
  ###############Step 1) Get presence records for species
  PA_species<-loc_clean[loc_clean$SPECIES %in% spec[x],]
  
  #How many records do we have?
  records<-nrow(PA_species)
  
  #Have the species names as a vector that will helpful later
  spname<-levels(factor((PA_species$SPECIES)))
  
  #get unique localities
  pts<-aggregate(PA_species,list(PA_species$LOCALITY),FUN=mean)
  
  #The format required is x,y,presence columns as their "response"
  #response<-cbind(PA_species[,c("LONGDECDEG","LATDECDEG")],1)
  #colnames(response)<-c("x","y",spname)
  
  #########################The above gets you presence only records
  #To get presence/psuedoabscence  
  #Get the presence absence matrix
    loc_matrix<-table(loc_clean$LOCALITY,loc_clean$SPECIES)
 
  #Select the species you'd like
  sp_matrix<-as.data.frame(melt(loc_matrix[,spec[x]]))
  unique.loc<-unique(loc_clean[,c("LOCALITY","LONGDECDEG","LATDECDEG")])
  p_a<-merge(sp_matrix,unique.loc,by.x="row.names",by.y="LOCALITY")
  p_a<-p_a[!duplicated(p_a),]
  
  #name the columns
  colnames(p_a)<-c("Locality","Response","LONGDECDEG","LATDECDEG")
  
  p_a[p_a$Response > 1,"Response"]<-1
  p_a<-p_a[!p_a$Locality=="",]
  p_a<-aggregate(p_a,list(p_a$Locality),FUN=mean)
  
  #the 0 are not true absences, they are NA's psuedoabsences
  #p_a[p_a$Response == 0,"Response"]<-NA
  
  #Always get the right number of psuedo absences, there are 2115 sites
  number_pSA<-nrow(p_a) - length(pts)
  
  #Format the data
  
  myBiomodData <- BIOMOD_FormatingData(resp.var = p_a[,"Response"],
                                       expl.var = stack(myExpl),
                                       resp.xy = p_a[,c("LONGDECDEG","LATDECDEG")],
                                       resp.name = gsub(" ","_",spec[x]),
  )
  
  plot(myBiomodData)
  #Define modeling options
  myBiomodOption <- BIOMOD_ModelingOptions(    
    MAXENT = list( path_to_maxent.jar = "D:/Niche_Models/maxent.jar",
                   maximumiterations = 200,
                   visible = TRUE,
                   linear = TRUE,
                   quadratic = TRUE,
                   product = TRUE,
                   threshold = TRUE,
                   hinge = TRUE,
                   lq2lqptthreshold = 80,
                   l2lqthreshold = 10,
                   hingethreshold = 15,
                   beta_threshold = -1,
                   beta_categorical = -1,
                   beta_lqp = -1,
                   beta_hinge = -1,
                   defaultprevalence = 0.5)
  )
   
  #Give current project a name, so we can go get the files later
  projnam<-'current'
  myBiomodModelOut<-BIOMOD_Modeling( myBiomodData, 
                                     models = c("GLM","GBM","MAXENT"), 
                                     models.options = myBiomodOption, 
                                     NbRunEval=2, 
                                     DataSplit=80, 
                                     Yweights=NULL, 
                                     VarImport=3, 
                                     models.eval.meth = c('ROC'),
                                     SaveObj = TRUE )
  
  # get all models evaluation                                     
  myBiomodModelEval <- getModelsEvaluations(myBiomodModelOut)
  
  
  # print the dimnames of this object
  dimnames(myBiomodModelEval)
  
  # let's print the ROC and TSS scores of all selected models, get the mean value for all the combined runs.
  stat<-myBiomodModelEval["ROC","Testing.data",,"Full",]
  
  #need to write this to file
  filename<-paste(paste(getwd(),gsub(" ",".",spec[x]),sep="/"),"ModelEval.csv",sep="/")
  write.csv(cbind(spec[x],stat),filename)
  
  #Let's look at variable importance
  m.var<-melt(getModelsVarImport(myBiomodModelOut)[,,"Full",])
  c.var<-cast(m.var,X1~X2)
  
  #Write variable importance to file
  filename<-paste(paste(getwd(),gsub(" ",".",spec[x]),sep="/"),"VarImportance.csv",sep="/")
  write.csv(cbind(c.var,spec[x]),filename)
  
  #Ensemble model outputs
  
  # get evaluation scores??, as is see it i want the data for the full runs of the ensemble mean value
  #ens_stat<-melt(getEMeval(myBiomodEM))
  #ens_stat[ens_stat$x3=="em.mean",]
  
  # projection over the globe under current conditions  
  #The clamping issue is a big one here!

    #Ensemble model
  myBiomodEM <- BIOMOD_EnsembleModeling( 
    modeling.output = myBiomodModelOut,
    chosen.models = 'all',
    em.by='all',
    eval.metric = c('ROC'),
    eval.metric.quality.threshold = c(0.75),
    prob.mean = T,
    prob.cv = T,
    prob.ci = T,
    prob.ci.alpha = 0.05,
    prob.median = T,
    committee.averaging = T,
    prob.mean.weight = T,
    prob.mean.weight.decay = 'proportional' )
  
  #################
  #Project SDM into env projections
  #All the gcm and current worldclim layer (first) are put together in a list
  #Loop through this list, only run if it has not been run before
  #name it corrently
bio_project<-function(GCM,nam){
  paste("Running Env", nam)
  if(!gsub(" ",".",spec[x]) %in% completed_models[[nam]]){
  myBiomodProjection<- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = GCM,
    proj.name = nam,
    selected.models = 'all',
    binary.meth = 'ROC',
    compress = 'xz',
    clamping.mask = T)
  
  #Ensemble projection
  EnsBas<-BIOMOD_EnsembleForecasting(projection.output = myBiomodProjection, EM.output = myBiomodEM)
  }}
  
  
  mapply(bio_project,projEnv,names(projEnv))
  ##################################
  
  #end file output
  print(paste("End Time is",system.time()))
  sink()
  return(stat)
})
stopCluster(cl)

print("ModelsComplete")

############################Get model evaluation stats
###########################

#Get the model evaluation from file
model_eval<-list.files(full.name=TRUE,recursive=T,pattern="Eval.csv")
model_eval<-rbind.fill(lapply(model_eval,read.csv))
colnames(model_eval)<-c("Model","Species","Stat")
model_eval<-melt(model_eval,id.var=c("Model","Species","Stat"))
#model_eval<-cast(model_eval,Species~Model)


#remove NA's?
ggplot(model_eval, aes(x=Species,y=Model,fill=Stat)) + geom_tile() + scale_fill_gradient("ROC",limits=c(0,1),low="blue",high="red",na.value="white") + opts(axis.text.x=theme_text(angle=-90))
ggsave("ModelEvaluations.jpeg")

ggplot(model_eval, aes(x=Species,y=Model,fill=Stat)) + geom_tile() + scale_fill_gradient("ROC",limits=c(0,1),low="blue",high="red",na.value="white") + opts(axis.text.x=theme_text(angle=-90))

model_thresh<-sapply(seq(.5,.95,.05),function(x){
  table(model_eval$Stat > x,model_eval$Model)["TRUE",]
})

colnames(model_thresh)<-seq(.5,.95,.05)
model_thresh<-melt(model_thresh)

names(model_thresh)<-c("Model","ROC_Threshold","Number_of_Species")
ggplot(model_thresh,aes(x=ROC_Threshold,y=Number_of_Species,col=Model)) + geom_line() + geom_point() + geom_text(aes(label=Number_of_Species),vjust=4,size=5)
ggsave("ModelThresholding.jpeg",dpi=300,height=8,width=8)
#Get the variable importance from file
varI<-list.files(full.name=TRUE,recursive=T,pattern="VarImportance.csv")
varI<-rbind.fill(lapply(varI,read.csv))
varI<-varI[,-1]

#Melt variable for plotting
mvar<-melt(varI)
colnames(mvar)<-c("Bioclim","Species","Model","value")

#Plot variable importance across all models
ggplot(mvar, aes(x=Species,y=Bioclim,fill=value)) + geom_tile() + scale_fill_gradient(limits=c(0,1),low="blue",high="red",na.value="white") + opts(axis.text.x=theme_text(angle=-90)) + facet_grid(Model ~ .)

ggsave("VariableImportance.jpeg")
}

print("SDM Function Defined")