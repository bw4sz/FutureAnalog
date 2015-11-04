# Load required packages
packages <- c("vegan", "picante", "analogue", "doSNOW", "ape", "cluster", 
              "RColorBrewer", "raster", "ggplot2", "phylobase", "rgdal", 
              "tidyr", "stringr", "dplyr", "biomod2", "rasterVis", "grid", 
              "devtools", "broom", "gridExtra", "proxy", "geometry", "rcdd")

for(p in packages) {
  if (!p %in% installed.packages()) {
    install.packages(p)
  }
  require(p, character.only = TRUE)
}

# Load in source functions
source("AlphaMappingFunctions.R")
source("BenHolttraitDiversity.R")
source("FutureAnalogFunctions.R")
source("runSDM.R")
source("fnSDM.R")
source("runProjections.R")
source("FutureAnalog.R")
source("TraitMappingFunctions.R")
source("BetaSorFunc.R")

# variables for sensitivity analysis
cell_size = 0.0833333333 # cell size numerical (degrees)
# must be a multiple of 5 arcmins or will need alternative data source
cell = "5_arcmins" # cell size string (text for folder name)

# create folders to output the models to
output_folder = "../FutureAnalog_output" 
out_path <- paste(output_folder, cell, sep = "/")
if(!dir.exists(output_folder)) dir.create(output_folder)
if(!dir.exists(out_path)) dir.create(out_path)
proj_folder <- getwd()

# Run SDMs
runSDM(cell_size, out_path, proj_folder)

# Project into future climate (will go into the folder with the climate
# projections and run species projections for all of these)
runProjections(cell_size, out_path, proj_folder)

# Create matrices of between time betadiversity and save them to the output
# folder (one per climate scenario)
runBetaDiv(out_path, cell_size)

# Run analog analysis (creates raster stacks of the number of novel and 
# disappearing communities - taxa, phylo, functional - for each climate model) 
# At the moment it's for 5%, 10%, 20% and 50% community similarity as the
# threshold for analog
runAnalogAnalysis(0.05, out_path)
runAnalogAnalysis(0.1, out_path)
runAnalogAnalysis(0.2, out_path)
runAnalogAnalysis(0.5, out_path)

# Create the plots for the MS
source("runFAPlots.R")

