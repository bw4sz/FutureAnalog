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

# cell size - currently analyses are being done at 0.1 degrees, but might do 
# sensitivity analysis
cell_size = 0.1

# load and clean species data --------------------------------------------------
# Lets go get the presence data on hummingbird distributions
PA <- read.csv("InputData/MASTER_POINTLOCALITYarcmap_review.csv")

# Just take the columns you want. 
PAdat <- select(PA, RECORD_ID, SPECIES, COUNTRY, LOCALITY, LATDECDEG, 
                LONGDECDEG, Decision, SpatialCheck, MapDecision)

# Just get the clean localities
gooddata <- c("ok", "Ok", "OK", "Y") #note blanks and REJECT data are excluded
loc_clean <- filter(PAdat, SpatialCheck=="Y", MapDecision %in% gooddata)

# get into site by species format
sppXsite <- select(loc_clean, RECORD_ID, LATDECDEG, LONGDECDEG, SPECIES) %>%
                   mutate(pres = 1) %>%
                   spread(SPECIES, pres)
sppXsite[is.na(sppXsite)] <- 0

# load climate data - predictors -----------------------------------------------
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
                             predData=envData)

# fit GDM ----------------------------------------------------------------------
gdm.model <- gdm(site.pair)
summary.gdm(gdm.model)

predDiss <- predict(gdm.model, site.pair)

plot(site.pair$distance, predDiss)
