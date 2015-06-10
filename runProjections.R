# runProjections.R -------------------------------------------------------------
# This code will run the projections for any GCM / Species combination for which
# the projections have not already been created.

# load required packages (installing if not already done) ----------------------
packages <- c("raster", "dplyr", "biomod2")

for(p in packages) {
  if (!p %in% installed.packages()) {
    install.packages(p)
  }
  require(p, character.only = TRUE)
}

# load fnSDM.R
source("fnSDM.R")


# set cell size
cell_size = 0.0833333333
cell = "5_arcmins"

output_folder = "../FutureAnalog_output" 
out_path <- paste(output_folder, cell, sep = "/")

# Step 1) List the future climate data available and collate variables ---------
# Create list of climate models
clim.mods <- list.files("../worldclim_data/projections_2070/")

for (mod in clim.mods) {
  # get the required climate variables
  mod.var1 <- list.files(paste("../worldclim_data/projections_2070/", mod, sep="/"),
                         pattern = "701.tif$", full.names = TRUE)
  
  mod.var12 <- list.files(paste("../worldclim_data/projections_2070/", mod, sep="/"),
                         pattern = "7012.tif$", full.names = TRUE)
  
  mod.var15 <- list.files(paste("../worldclim_data/projections_2070/", mod, sep="/"),
                         pattern = "7015.tif$", full.names = TRUE)
  
  mod.stack <- stack(mod.var1, mod.var12, mod.var15)
  
  # make the names consistent
  names(mod.stack) <- c("bio1", "bio12", "bio15")
  
  # Step 2) Set the Extent to project *into*. ----------------------------------
  
  # Presence points are still taken from everywhere Avoid projecting into areas 
  # where sample size is really low **DECISION**
  exte<-extent(c(-81.13411,-68.92061,-5.532386,11.94902))
  
  # Cut by the extent. Crop by this layer 
  mod.stack.c <- stack(crop(mod.stack, exte))
  
  # set resolution for the future layers equivalent to cell size
  fact <- cell_size/res(mod.stack.c)
 if(round(fact)[1] > 1) {
    # Set cell size to ~ cell_size degree
    mod.stack.c <- stack(aggregate(mod.stack.c, fact))
  }
  
  # Step 3) Get a list of the species which have SDMs already ------------------
  setwd(out_path)
  spec <- list.files(".")
  spec <- spec[file.info(spec)$isdir]
  spec <- spec[!grepl("logs",spec)] # don't want the logs folder included
  
  # Step 4) Create projections into future climate -----------------------------
  for (sp in spec){
    # project into future climate (bio_project() defined in fnSDM.R)
    if(!file.exists(paste0(sp, "/proj_", mod))) bio_project(sp, mod.stack.c, mod)
  }
  setwd("../../FutureAnalog/")
}
