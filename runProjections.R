# runProjections.R -------------------------------------------------------------
# This code will run the projections for any GCM put as input. 
# Step 3) Climate Scenarios and Futute Climate #################################
# TO DO: Once code is working fully, update the climate scenarios loading to
# collate all climate scenarios in the folder - this will involve searching the
# "worldclim_data" folder, creating a list, and 'listifying' all of the climate
# layer bits and pieces. Should be straightforward and will make automation
# easier. Also, think about removing climate scenarios which have already been
# run from the list?

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





for (i in 1:length(projEnv)){
  bio_project(projEnv[[i]], names(projEnv[i]))
}
