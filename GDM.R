# GDM.R ------------------------------------------------------------------------
# to explore using generalized dissimilarity models as an alternate (additional)
# approach for the betadiversity analysis

# load required packages (installing if not already done) ----------------------
packages <- c("gdm", "dplyr", "tidyr")

for(p in packages) {
  if (!p %in% installed.packages()) {
    install.packages(p)
  }
  require(p, character.only = TRUE)
}

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
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              #                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   need to get rid of any rows with no environmental data (this is just one
# offshore record) create 
dat <- dat[which(complete.cases(dat)),]
dat <- data.frame(dat)
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
# site-pair table format -------------------------------------------------------
bioData <- data.frame(site = rownames(dat), long = dat$x, lat = dat$y, dat[,6:ncol(dat)])
envData <- data.frame(site = rownames(dat), long = dat$x, lat = dat$y, dat[,3:5])
site.pair2 <- formatsitepair(bioData, 1, siteColumn="site", XColumn="long", YColumn="lat",
                             predData=envData, weightType = "richness")

# fit GDM ----------------------------------------------------------------------
gdm.model <- gdm(site.pair2)
summary.gdm(gdm.model)

png("gdmplot.png")
plot.gdm(gdm.model, plot.layout = c(3,2))
dev.off()

