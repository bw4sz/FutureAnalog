# get list of sppXsite matrices for each environmental scenario
sppXsite.files <- list.files("sppXsite", full.names = TRUE)

# get the number/proportion of cells occupied for each species in each climate scenario
sXs.summary <- list()
for(sppXsite.file in sppXsite.files){
  load(sppXsite.file)
  
  # get climate model name
  clim.mod <- str_match(sppXsite.file,pattern=paste("sppXsite/(\\w+.\\w+).rda",sep="/"))[,2]
  
  # get rid of cell information
  sppXsite <- sppXsite[,2:(ncol(sppXsite)-3)]
  fails <- na.test(sppXsite)
  sppXsite <- sppXsite[,!colnames(sppXsite) %in% fails]
  sXs.summary[[clim.mod]] <- data.frame(species=colnames(sppXsite), no.cells=colSums(sppXsite), 
                            perc.occ=colSums(sppXsite)/nrow(sppXsite), clim.mod=clim.mod)
}

sXs.summary <- do.call("rbind", sXs.summary)

# load trait data
traits <- getTraitData()

# remove species from summary where we do not have trait data and vice versa
sXs.summary <- filter(sXs.summary, species %in% rownames(traits))
traits <- filter(traits, rownames(traits) %in% sXs.summary$species)

# separate out the current from future scenario results
sXs.current <- filter(sXs.summary, clim.mod=="current")
sXs.future <- filter(sXs.summary, clim.mod!="current")

# separate the GCM and RCP in the scenario identifier
sXs.future$GCM <- substr(sXs.future$clim.mod, 1, 2)
sXs.future$RCP <- substr(sXs.future$clim.mod, 3, 4)

# calculate the difference in proportion of cells occupied between the current
# and all future scenarios
sXs.future$diff <- sXs.current$perc.occ - sXs.future$perc.occ

# create data frames to be used for groups, point size and facets in the biplot
sXs.RCP <- group_by(sXs.future, RCP, species) %>%
  summarise(diff=mean(diff)) %>%
  mutate(groups=ifelse(diff < 0, "More cells occupied in future environment", "More cells occupied in current environment"),
         diff=abs(diff)*100,
         scenario=RCP)

sXs.GCM <- group_by(sXs.future, GCM, species) %>%
  summarise(diff=mean(diff)) %>%
  mutate(groups=ifelse(diff < 0, "More cells occupied in future environment", "More cells occupied in current environment"),
         diff=abs(diff)*100,
         scenario=GCM)

# PCA on the species traits
pca <- prcomp(traits)

# plots
ggbiplot2(pca, add.dat = sXs.RCP) + theme_classic() +
  labs(size="Absolute difference in % cells occupied", color="")
ggsave("Figures/traitsbyRCP.pdf", width=17, height=9)

ggbiplot2(pca, add.dat = sXs.GCM) + theme_classic() +
  labs(size="Absolute difference in % cells occupied", color="") + 
  theme(legend.position = c(1, 0), legend.justification = c(1, 0)) 
ggsave("Figures/traitsbyGCM.pdf", width=17, height=9)
