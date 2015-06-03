# runProjections.R -------------------------------------------------------------
# This code will run the projections for any GCM put as input. 

for (i in 1:length(projEnv)){
  bio_project(projEnv[[i]], names(projEnv[i]))
}
