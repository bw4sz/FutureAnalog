# GDM.R ------------------------------------------------------------------------
# to explore using generalized dissimilarity models as an alternate (additional)
# approach for the betadiversity analysis

fitGDM <- function() {
  # This function loads the species and current climate data and fits the GDM -
  # will probably need cleaning a bit
  
  # load and clean species data --------------------------------------------------
  # Lets go get the presence data on hummingbird distributions
  PA <- read.csv("InputData/MASTER_POINTLOCALITYarcmap_review.csv")
  
  # Just take the columns you want. 
  PAdat <- select(PA, RECORD_ID, SPECIES, COUNTRY, LOCALITY, LATDECDEG, 
                  LONGDECDEG, Decision, SpatialCheck, MapDecision)
  
  # Just get the clean localities
  gooddata <- c("ok", "Ok", "OK", "Y") #note blanks and REJECT data are excluded
  loc_clean <- filter(PAdat, SpatialCheck=="Y", MapDecision %in% gooddata)
  
  # select only well fitting species (currently set as all where TSS is over 0.5
  # for all models and ROC is over 7.5 for all models)
  model_eval<-list.files("../FutureAnalog_output/5_arcmins/", full.name=TRUE,recursive=T,pattern="Eval.csv")
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
  
  loc_clean <- filter(loc_clean, SPECIES %in% well.fitting.species)
  
  # get into site by species format
  sppXsite <- select(loc_clean, RECORD_ID, LATDECDEG, LONGDECDEG, SPECIES) %>%
    mutate(pres = 1) %>%
    spread(SPECIES, pres)
  sppXsite[is.na(sppXsite)] <- 0
  
  # load climate data - predictors -----------------------------------------------
  # Import environmental data from worldclim, three variables Bio1 = annual mean
  # temp, Bio12 = annual precip, Bio15 = precip seasonality
  myExpl <- c("../worldclim_data/bio_5m_bil/bio1.bil",
              "../worldclim_data/bio_5m_bil/bio12.bil",
              "../worldclim_data/bio_5m_bil/bio15.bil")
  
  myExpl <- stack(myExpl)
  
  # get extent of the presence data
  extPoint <- SpatialPoints(loc_clean[,c("LONGDECDEG","LATDECDEG")])
  
  # crop by this extent
  myExpl <- crop(myExpl,extPoint)
  ec <- readOGR("InputData", "EcuadorCut")
  cur <- crop(myExpl, ec)
  cur <- rasterToPoints(cur)
  cur <- data.frame(cur)
  # get data info correct format for bioFormat 4 (2 wouldn't work :( )) ----------
  # get sppXsite into raster format
  sppXsite.r <- rasterize(select(sppXsite, LONGDECDEG, LATDECDEG), myExpl[[1]],
                          select(sppXsite, -RECORD_ID, -LATDECDEG, -LONGDECDEG), 
                          fun = max)
  datastack <- stack(myExpl, sppXsite.r)
  
  # all data into table with lat long
  dat <- rasterToPoints(datastack)
  # remove rows without species
  dat <- subset(dat, rowSums(dat[,6:ncol(dat)]) > 0)
  # need to get rid of any rows with no environmental data (this is just one
  # offshore record) create 
  dat <- dat[which(complete.cases(dat)),]
  dat <- data.frame(dat)
  
  # site-pair table format -------------------------------------------------------
  bioData <- data.frame(site = rownames(dat), long = dat$x, lat = dat$y, dat[,6:ncol(dat)])
  envData <- data.frame(site = rownames(dat), long = dat$x, lat = dat$y, dat[,3:5])
  site.pair <- formatsitepair(bioData, 1, siteColumn="site", XColumn="long", YColumn="lat",
                              predData=envData, weightType = "richness")
  
  # fit GDM ----------------------------------------------------------------------
  gdm.model <- gdm(site.pair)
  gdm.out <- list(gdm.model = gdm.model, cur = cur)
  return(gdm.out)
}


createSitePairs <- function(cur, fu) {
  ec <- readOGR("InputData", "EcuadorCut")
  
  fu.dat <- c(paste0("../worldclim_data/projections_2070/", fu, "/" , fu, "1.tif"),
              paste0("../worldclim_data/projections_2070/", fu, "/" , fu, "12.tif"),
              paste0("../worldclim_data/projections_2070/", fu, "/" , fu, "15.tif"))
  
  fu.dat <- stack(fu.dat)
  fu.dat <- crop(fu.dat, ec)
  fu.dat <- rasterToPoints(fu.dat)
  fu.dat <- data.frame(fu.dat)
  names(fu.dat) <- names(cur)
  fu.dat$fuID <- 1:nrow(fu.dat)
  cur$curID <- 1:nrow(cur)
  site.pairs <- expand.grid(cur$curID, fu.dat$fuID)
  site.pairs <- merge(site.pairs, cur, by.x="Var1", by.y="curID")
  site.pairs <- merge(site.pairs, fu.dat, by.x="Var2", by.y="fuID")
  site.pairs <- data.frame(distance = rep(0, nrow(site.pairs)), 
                           weights = rep(0, nrow(site.pairs)), 
                           xCoord.S1 = site.pairs$x.x,
                           yCoord.S1 = site.pairs$y.x,
                           xCoord.S2 = site.pairs$x.y,
                           yCoord.S2 = site.pairs$y.y,
                           s1.bio1 = site.pairs$bio1.x, 
                           s1.bio12 = site.pairs$bio12.x, 
                           s1.bio15 = site.pairs$bio15.x,
                           s2.bio1 = site.pairs$bio1.y, 
                           s2.bio12 = site.pairs$bio12.y, 
                           s2.bio15 = site.pairs$bio15.y, 
                           cur = site.pairs$Var1,
                           fu = site.pairs$Var2)
  return(site.pairs)
}

predictGDM <- function(gdm.mod, site.pairs){
  predDiss <- predict(gdm.mod, site.pairs)
  predDiss <- cbind(site.pairs, predDiss)
  betadiv <- spread(predDiss)
}
