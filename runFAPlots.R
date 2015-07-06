rasterToDataFrame <- function(out_path){
  
  out.rasters <- list.files(out_path, pattern = "NonAnalogRasters", full.names = TRUE, recursive = TRUE)
  
  out <- lapply(out.rasters, function(x) {
    load(x)
    
    arbthresh <- str_match(x,pattern=paste(cell,"(\\w+.\\w+)/NonAnalog",sep="/"))[,2]
    GCM <- substr(str_match(x,pattern=paste("Rasters_(\\w+.\\w+)bi70",sep="/"))[,2], 1, 2)
    RCP <- substr(str_match(x,pattern=paste("Rasters_(\\w+.\\w+)bi70",sep="/"))[,2], 3, 4)
    
    dat <- as.data.frame(results, xy = TRUE) %>%
      gather(key, value, -x, -y, na.rm = TRUE) %>%
      mutate(arbthresh = arbthresh, GCM = GCM, RCP = RCP) %>%
      separate(variable, into = c("comm.type", "measure"), sep = "\\.")
  }
  )
  out <- do.call("rbind", out)
}

dat <- rasterToDataFrame(out_path)

dat$measure <- factor(dat$measure, levels=c("Tax", "Phylo", "Func"),
                      labels=c("Taxonomic", "Phylogenetic", "Functional"))

dat$RCP <- factor(dat$RCP, levels=c(26, 45, 85), 
                  labels=c("RCP 2.6", "RCP 4.5", "RCP 8.5"))

dat$GCM <- factor(dat$GCM, levels=unique(dat$GCM),
                  labels=c("CCSM4", "CNRM-CM5", "GISS-E2-R", "HadGEM2-ES", 
                           "IPSL-CM5A-LR", "MIROC5", "MPI-ESM-LR"))


# Sensitivity analysis plots by threshold
thresh_sa.novel <- filter(dat, comm.type=="Novel") %>% 
  group_by(x, y, measure, arbthresh) %>%
  summarise(NoOfAnalogs = mean(value))

ggplot(thresh_sa.novel, aes(x, y, fill = NoOfAnalogs)) + 
  scale_fill_gradient(low = "blue", high = "white", name="Number of analog\ncommunities") +
  labs(x="Longitude", y="Latitude") +
  facet_grid(measure ~ arbthresh) + 
  geom_raster() + coord_equal() + theme_classic(base_size=15) + 
  theme(strip.background = element_blank(), panel.margin = unit(2, "lines"))

ggsave("Figures/Threshold_SA_Novel.png", width = 17, height = 9)

thresh_sa.diss <- filter(dat, comm.type=="Disappearing") %>% 
  group_by(x, y, measure, arbthresh) %>%
  summarise(NoOfAnalogs = mean(value))

ggplot(thresh_sa.diss, aes(x, y, fill = NoOfAnalogs)) + 
  scale_fill_gradient(low = "red", high = "white", name="Number of analog\ncommunities") +
  facet_grid(measure ~ arbthresh) + 
  labs(x="Longitude", y="Latitude") +
  geom_raster() + coord_equal() + theme_classic(base_size=15) + 
  theme(strip.background = element_blank(), panel.margin = unit(2, "lines"))

ggsave("Figures/Threshold_SA_Disappearing.png", width = 17, height = 9)

# Novel average of GCMs (res for each RCP shown)
novel.20.rcp <- subset(dat, comm.type=="Novel", arbthresh = 0.2) %>%
  group_by(x, y, measure, RCP) %>%
  summarise(NoOfAnalogs = mean(value))

ggplot(novel.20.rcp, aes(x, y, fill = NoOfAnalogs)) + 
  scale_fill_gradient(low = "blue", high = "white", name="Number of analog\ncommunities") +
  facet_grid(measure ~ RCP) + 
  labs(x="Longitude", y="Latitude") +
  geom_raster() + coord_equal() + theme_classic(base_size=15) + 
  theme(strip.background = element_blank(), panel.margin = unit(2, "lines"))

ggsave("Figures/Novel_by_RCP_20perc_thres.png", width = 17, height = 9)

novel.20.gcm <- subset(dat, comm.type=="Novel", arbthresh = 0.2) %>%
  group_by(x, y, measure, GCM) %>%
  summarise(NoOfAnalogs = mean(value))

ggplot(novel.20.gcm, aes(x, y, fill = NoOfAnalogs)) + 
  scale_fill_gradient(low = "blue", high = "white", name="Number of analog\ncommunities") +
  facet_grid(measure ~ GCM) + 
  labs(x="Longitude", y="Latitude") +
  geom_raster() + coord_equal() + theme_classic(base_size=15) + 
  theme(strip.background = element_blank(), panel.margin = unit(2, "lines"))

ggsave("Figures/Novel_by_GCM_20perc_thres.png", width = 18, height = 9)

# Disappearing average of GCMs (res for each RCP shown)
diss.20.rcp <- subset(dat, comm.type=="Disappearing", arbthresh = 0.2) %>%
  group_by(x, y, measure, RCP) %>%
  summarise(NoOfAnalogs = mean(value))

ggplot(diss.20.rcp, aes(x, y, fill = NoOfAnalogs)) + 
  scale_fill_gradient(low = "red", high = "white", name="Number of analog\ncommunities") +
  facet_grid(measure ~ RCP) + 
  labs(x="Longitude", y="Latitude") +
  geom_raster() + coord_equal() + theme_classic(base_size=15) + 
  theme(strip.background = element_blank(), panel.margin = unit(2, "lines"))

# Disappearing average of RCPs (res for each GCM shown)
ggsave("Figures/Disappearing_by_RCP_20perc_thres.png", width = 17, height = 9)

diss.20.gcm <- subset(dat, comm.type=="Disappearing", arbthresh = 0.2) %>%
  group_by(x, y, measure, GCM) %>%
  summarise(NoOfAnalogs = mean(value))

ggplot(diss.20.gcm, aes(x, y, fill = NoOfAnalogs)) + 
  scale_fill_gradient(low = "red", high = "white", name="Number of analog\ncommunities") +
  facet_grid(measure ~ GCM) + 
  labs(x="Longitude", y="Latitude") +
  geom_raster() + coord_equal() + theme_classic(base_size=15) + 
  theme(strip.background = element_blank(), panel.margin = unit(2, "lines"))

ggsave("Figures/Disappearing_by_GCM_20perc_thres.png", width = 18, height = 9)

