runSppXsite <- function(out_path){
  # select only well fitting species (currently set as all where TSS is over 0.5
  # for all models and ROC is over 7.5 for all models)
  model_eval<-list.files(out_path, full.name=TRUE,recursive=T,pattern="Eval.csv")
  model_eval<-rbind_all(lapply(model_eval, 
                               function(x) read.csv(x, stringsAsFactors = FALSE)))
  colnames(model_eval)[1:2] <- c("Stat", "Species")
  
  # change here for sensitivity analyis
  model_eval$TEST <- (apply(model_eval, 1, function(x) min(x[3:5], na.rm=TRUE) < 0.5) 
                      & model_eval$Stat == "TSS") |
    (apply(model_eval, 1, function(x) min(x[3:5], na.rm=TRUE) < 0.75) 
     & model_eval$Stat == "ROC")
  
  well.fitting.models <- subset(model_eval, !TEST)
  well.fitting.species <- unique(well.fitting.models$Species)
  well.fitting.species <- gsub(" ", ".", well.fitting.species)
  
  # get the crop files
  niche.crops <- list.files(out_path,pattern="crop.gri",full.name=T,recursive=T)
  
  files <- c()
  for(s in well.fitting.species){
    files <- c(files,grep(s, niche.crops))
  }
  
  niche.crops <- niche.crops[files]
  
  # create folder to store results
  if(!dir.exists(paste(out_path, "sppXsite", sep="/"))) 
    dir.create(paste(out_path, "sppXsite", sep="/"))
  
  # get list of environmental layers
  env.list <- list.files("../worldclim_data/projections_2070/")
  env.list <- c("current", env.list)
  
  for(env in env.list){
    f <- niche.crops[grep(env, niche.crops, value=FALSE)]
    sppXsite <- tableFromRaster(f, threshold = 0.05)
    sppXsite <- data.frame(sppXsite)
    sppXsite$cell <- rownames(sppXsite)
    xy <- as.data.frame(raster(f[[1]]), xy=TRUE)
    xy$cell <- rownames(xy)
    sppXsite <- merge(sppXsite, xy)
    save(sppXsite, file=paste0(out_path, "/sppXsite/", env, ".rda"))
    rm(sppXsite)
  }
}

getTraitData <- function() {
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
  z.scores$species <- rownames(mon)
  return(z.scores)
}