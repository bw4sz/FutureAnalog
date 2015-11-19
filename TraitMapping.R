# get list of sppXsite matrices for each environmental scenario
sppXsite.files <- list.files("sppXsite", full.names = TRUE)

# load trait data
traits <- getTraitData()

# load hillshade data for plotting
load("InputData/srtm_5arcmin.rda")
slope = terrain(elev, opt='slope')
aspect = terrain(elev, opt='aspect')
hill = hillShade(slope, aspect, 40, 270)
hdf <- rasterToPoints(hill)
hdf <- data.frame(hdf)

# ecuador boundary - get into format for ggplot
ec <- readOGR("InputData", "EcuadorCut")
ec@data$id = rownames(ec@data)
ec.df = fortify(ec, region="id")


# APPROACH #1 ------------------------------------------------------------------
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

# APPROACH #2 ------------------------------------------------------------------
# get current species list and calculate mean traits for each site
current <- list.files("sppXsite", pattern="current", full.names = TRUE)
currentTraits <- getSiteTraitValues(current, traits)

# get future species lists and calculate mean traits for each site and each model
future <- sppXsite.files[which(sppXsite.files != current)]
futureTraits <- list()
for(f in future) {
  futureTraits[[f]] <- getSiteTraitValues(f, traits)
}

futureTraits <- do.call("rbind", futureTraits)

# get the mean trait values grouped by RCP
futureTraits <- mutate(futureTraits, mod=substr(mod, 3, 4)) %>%
  group_by(x, y, mod) %>%
  summarise(Bill=mean(Bill, na.rm=TRUE), 
            Mass=mean(Mass, na.rm=TRUE), 
            WingChord=mean(WingChord, na.rm=TRUE))

# bind together and get the values for the first two principal components
allTraits <- rbind(currentTraits, futureTraits)
pca <- princomp(allTraits[,4:6])
allTraits$PCA1 <- pca$scores[,1]
allTraits$PCA2 <- pca$scores[,2]

# create subsets for current and the three RCPs
currentTraits <- subset(allTraits, mod=="current")
RCP26Traits <- subset(allTraits, mod=="26")
RCP45Traits <- subset(allTraits, mod=="45")
RCP85Traits <- subset(allTraits, mod=="85")

# calculate the diff between current and each future scenario for PCA
RCP26PCADiff <- data.frame(mod="RCP 2.6", x=currentTraits$x, y=currentTraits$y, currentTraits[,7:8] - RCP26Traits[,7:8])
RCP45PCADiff <- data.frame(mod="RCP 4.5", x=currentTraits$x, y=currentTraits$y, currentTraits[,7:8] - RCP45Traits[,7:8])
RCP85PCADiff <- data.frame(mod="RCP 8.5", x=currentTraits$x, y=currentTraits$y, currentTraits[,7:8] - RCP85Traits[,7:8])

# calculate the diff between current and each future scenario for actual trait value
RCP26TraitDiff <- data.frame(mod="RCP 2.6", x=currentTraits$x, y=currentTraits$y, currentTraits[,4:6] - RCP26Traits[,4:6])
RCP45TraitDiff <- data.frame(mod="RCP 4.5", x=currentTraits$x, y=currentTraits$y, currentTraits[,4:6] - RCP45Traits[,4:6])
RCP85TraitDiff <- data.frame(mod="RCP 8.5", x=currentTraits$x, y=currentTraits$y, currentTraits[,4:6] - RCP85Traits[,4:6])

# get into dataframe for plotting
PCADiff <- rbind(RCP26PCADiff, RCP45PCADiff, RCP85PCADiff)
PCADiff <- gather(PCADiff, key, value, -x, -y, -mod)

TraitDiff <- rbind(RCP26TraitDiff, RCP45TraitDiff, RCP85TraitDiff)
TraitDiff <- gather(TraitDiff, key, value, -x, -y, -mod)

# PCA difference values plotted
PCAPlot <- ggplot(NULL, aes(x, y)) + 
  geom_raster(data = PCADiff, aes(fill=value)) +
  geom_raster(data = hdf, aes(alpha=layer)) +
  geom_path(data = ec.df, aes(x=long, y=lat)) +
  scale_fill_gradient2(name="Current value - Future value") +
  guides(fill = guide_colorbar()) +
  scale_alpha(range = c(0, 0.5), guide = "none") +
  facet_grid(variable ~ mod) + 
  scale_x_continuous(name=expression(paste("Longitude (", degree, ")"))) + 
  scale_y_continuous(name=expression(paste("Latitude (", degree, ")"))) +
  coord_equal() + theme_classic(base_size=15) + 
  theme(strip.background = element_blank(), panel.margin = unit(2, "lines"))

# biplot for PCA so that the axes make sense
require(ggbiplot)
traitBiplot <- ggbiplot(pca, alpha=0.005) + theme_classic()

png("Figures/trait_pca_changes.png", width=1311, height=515)
grid.arrange(PCAPlot, traitBiplot, nrow=1, widths=c(2, 1))
dev.off()

# Actual trait value differences plotted
ggplot(NULL, aes(x, y)) + 
  geom_raster(data = TraitDiff, aes(fill=value)) +
  geom_raster(data = hdf, aes(alpha=layer)) +
  geom_path(data = ec.df, aes(x=long, y=lat)) +
  scale_fill_gradient2(name="Current value - Future value") +
  guides(fill = guide_colorbar()) +
  scale_alpha(range = c(0, 0.5), guide = "none") +
  facet_grid(variable ~ mod) + 
  scale_x_continuous(name=expression(paste("Longitude (", degree, ")"))) + 
  scale_y_continuous(name=expression(paste("Latitude (", degree, ")"))) +
  coord_equal() + theme_classic(base_size=15) + 
  theme(strip.background = element_blank(), panel.margin = unit(2, "lines"))

ggsave("Figures/trait_changes.png", width=9, height=9)

# PLOTTING THE HYPERVOLUMES
require(hypervolume)
strt <- system.time()
# get file names
files <- list.files("sppXsite", pattern="85", full.names = TRUE)
cur <- list.files("sppXsite", pattern="current", full.names = TRUE)
files <- c(files, cur)

# get trait data
traits <- getTraitData()

# calculate hypervolume
hv_res <- list()

for(f in files) {
  hv_res[[f]] <- calcHV(f)
}

hv_res <- do.call("cbind",hv_res)
future <- rowMeans(hv_res[,which(colnames(hv_res)!=cur)])
current <- hv_res[,cur]
load(cur)
sppXsite <- sppXsite[1:10,]
xy <- sppXsite[,c("x", "y")]
trait_vol <- data.frame(xy, current, future)
print(paste0("hypervolume calculations took ", system.time()-strt))
      