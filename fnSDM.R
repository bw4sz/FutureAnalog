# Code to run the SDMs ---------------------------------------------------------


# Input: loc_clean, spec, myExpl

# Output: Maxent, GLM, GBM and ensemble models. 

# Packages:
SDM_SP <- function(spec, loc_clean, myExpl) {
  packages <- c("biomod2", "dplyr")
  
  for(p in packages) {
    if (!p %in% installed.packages()) {
      install.packages(p)
    }
    require(p, character.only = TRUE)
  }
  print(paste0("Modelling for ", spec))
  strt <- Sys.time()
  # Step 1) Get presence/pseudo-absence records for species --------------------
  
  # get the presence data for the species, select necessary columns
  p_species <- filter(loc_clean, SPECIES == spec) %>%
    select(LOCALITY, LONGDECDEG, LATDECDEG) %>%
    mutate(PRES = 1)
  p_species <- unique(p_species) # remove duplicate records
  
  # get unique locations
  unique.loc<-unique(loc_clean[,c("LOCALITY","LONGDECDEG","LATDECDEG")])
  
  # join presence data onto the unique locations keeping all locations (pseudoabsences will
  # automatically be NA)
  p_a <- merge(unique.loc, p_species, all.x = TRUE)
  
  #name the columns
  colnames(p_a)<-c("Locality", "LONGDECDEG", "LATDECDEG", "Response")
  
  # Make sure all locations have the same lat long (created as the mean of all lat
  # longs for that location) - 
  # NB. I'M REALLY NOT SURE ABOUT THIS STEP - NEED TO CHECK.
  #p_a<-group_by(p_a, Locality) %>%
  #  summarise_each("mean")
  
  #we want 2000 psuedoabsences, randomly pick 2000 NA rows. 
  if(length(which(is.na(p_a$Response))) < 2000){
    psuedos<-which(is.na(p_a$Response))
  }
  
  if(length(which(is.na(p_a$Response))) > 2000){
    psuedos<-sample(which(is.na(p_a$Response)),2000)
  }
  
  pres_pts<-which(!is.na(p_a$Response))
  
  #Bring together the presence and pseudoabsence data (biomod does not like tbl_df
  #format that's output by dplyr)
  p_a <- data.frame(p_a[c(pres_pts,psuedos),])
  
  #Format the data
  myBiomodData <- BIOMOD_FormatingData(resp.var = p_a[,"Response"],
                                       expl.var = stack(myExpl),
                                       resp.xy = p_a[,c("LONGDECDEG","LATDECDEG")],
                                       resp.name = gsub(" ","_",spec),
  )
  
  #Define modeling options - **DECISION**
  myBiomodOption <- BIOMOD_ModelingOptions(    
    MAXENT = list( path_to_maxent.jar = "../../FutureAnalog/maxent.jar",
                   maximumiterations = 200,
                   visible = FALSE,
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
  projnam <- 'current'
  
  # Individual model outputs -----------------------------------------------------
  # **DECISION**
  myBiomodModelOut<-BIOMOD_Modeling( myBiomodData, 
                                     models = c("GLM","GBM","MAXENT"), 
                                     models.options = myBiomodOption, 
                                     NbRunEval=2, 
                                     DataSplit=80, 
                                     Yweights=NULL, 
                                     VarImport=3, 
                                     models.eval.meth = c('ROC',"TSS"),
                                     SaveObj = TRUE )
  
  # get all models evaluation                                     
  myBiomodModelEval <- get_evaluations(myBiomodModelOut)
  
  # get ROC and TSS scores of all selected models (mean value for all the combined
  # runs) and write to file
  stat <- myBiomodModelEval[c("ROC","TSS"), "Testing.data", ,"Full",]
  filename <- paste(gsub(" ",".",spec), "ModelEval.csv", sep="/")
  write.csv(cbind(spec,stat),filename)
  
  # get variable importance and write to file
  m.var <- melt(get_variables_importance(myBiomodModelOut)[,,"Full",])
  c.var <- cast(m.var,X1~X2)
  filename <- paste(gsub(" ",".",spec),"VarImportance.csv",sep="/")
  write.csv(cbind(c.var,spec),filename)
  
  # Ensemble model outputs -------------------------------------------------------
  # projection over the globe under current conditions  
  #The clamping issue is a big one here!
  
  #Ensemble model - **DECISION**
  # only uses models that did not fail in previous step (myBiomodModelOut)
  myBiomodEM <- BIOMOD_EnsembleModeling( 
    modeling.output = myBiomodModelOut,
    chosen.models =  get_built_models(myBiomodModelOut),
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
  
  # save modelling outputs for use in env projections
  save(myBiomodModelOut, file = paste(gsub(" ", ".", spec), "myBiomodModelOut.rda", sep="/"))
  save(myBiomodEM, file = paste(gsub(" ", ".", spec), "myBiomodEM.rda", sep="/"))

  print(Sys.time() - strt)
}

# Project SDM into env projections ---------------------------------------------
# Input: spec, GCM, nam

# Output: projections for species (spec) into projected climate (GCM)
bio_project<-function(spec, GCM, nam){
  load(paste(gsub(" ", ".", spec), "myBiomodModelOut.rda", sep="/"))
  load(paste(gsub(" ", ".", spec), "myBiomodEM.rda", sep="/"))
  paste("Running Env", nam)
  # **DECISION** 
  myBiomodProjection <- BIOMOD_Projection(
    modeling.output = myBiomodModelOut,
    new.env = GCM,
    proj.name = nam,
    selected.models = 'all',
    binary.meth = c('ROC'),
    compress = 'xz',
    clamping.mask = T)
  
  #Ensemble projection
  EnsBas<-BIOMOD_EnsembleForecasting(projection.output = myBiomodProjection, EM.output = myBiomodEM)
}