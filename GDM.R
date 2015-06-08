# GDM.R -------------------------------------------------------------------------------
# script to explore using generalized dissimilarity models as an alternate (additional)
# approach for the betadiversity analysis

# load required packages (installing if not already done) ----------------------
packages <- c("gdm", "dplyr")

for(p in packages) {
  if (!p %in% installed.packages()) {
    install.packages(p)
  }
  require(p, character.only = TRUE)
}

# cell size - currently analyses are being done at 0.1 degrees, but might do sensitivity analysis
cell_size = 0.1

# load and clean species data 
# Lets go get the presence data on hummingbird distributions
PA <- read.csv("InputData/MASTER_POINTLOCALITYarcmap_review.csv")

# Just take the columns you want. 
PAdat <- select(PA, RECORD_ID, SPECIES, COUNTRY, LOCALITY, LATDECDEG, 
                LONGDECDEG, Decision, SpatialCheck, MapDecision)

# Just get the clean localities
gooddata <- c("ok", "Ok", "OK", "Y") #note blanks and REJECT data are excluded
loc_clean <- filter(PAdat, SpatialCheck=="Y", MapDecision %in% gooddata)

# site by species matrix (from locations data, at cell_size resolution)

# load climate data - predictors

# create site-pair table format

# fit GDM

