# list of required climate models
clim.mods <- c(
  "cc26bi70", "cc45bi70", "cc85bi70", "cn26bi70", "cn45bi70", "cn85bi70", 
  "gs26bi70", "gs45bi70", "gs85bi70", "he26bi70", "he45bi70", "he85bi70", 
  "ip26bi70", "ip45bi70", "ip85bi70", "mp26bi70", "mp45bi70", "mp85bi70",
  "mc26bi70", "mc45bi70", "mc85bi70"
)

# download and unzip climate model one by one
for(mod in clim.mods){
  download.file(paste0("http://biogeo.ucdavis.edu/data/climate/cmip5/5m/", mod, ".zip"),
                paste0("../worldclim_data/projections_2070/", mod, ".zip"))
  unzip(paste0("../worldclim_data/projections_2070/", mod, ".zip"),
        exdir = paste0("../worldclim_data/projections_2070/", mod))
}