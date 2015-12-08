# GET DATA FOR THE CURRENT/FUTURE COMPARISON -----------------------------------


rasterToDataFrame <- function(out_path){ 
  
  out.rasters <- list.files(out_path, pattern = "NonAnalogRasters_", full.names = TRUE, recursive = TRUE)
  
  out <- lapply(out.rasters, function(x) { 
    load(x)
    
    results <- stack(get("results"), elev)
    
    arbthresh <- str_match(x,pattern=paste(cell,"(\\w+.\\w+)/NonAnalog",sep="/"))[,2] 
    GCM <- substr(str_match(x,pattern="NonAnalogRasters_(\\w+.\\w+)bi70")[,2], 1, 2) 
    RCP <- substr(str_match(x,pattern="NonAnalogRasters_(\\w+.\\w+)bi70")[,2], 3, 4)  
    
    dat <- as.data.frame(results, xy = TRUE) %>%
      gather(key, value, -x, -y, -output_srtm, na.rm = TRUE) %>%
      mutate(arbthresh = arbthresh, GCM = GCM, RCP = RCP) %>%
      separate(key, into = c("comm.type", "measure"), sep = "\\.")
  } 
  ) 
  out <- do.call("rbind", out)
}

# hillshade data
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

# get results data
dat <- rasterToDataFrame(out_path)

dat$measure <- factor(dat$measure, levels=c("Tax", "Phylo", "Func"),
                      labels=c("Taxonomic", "Phylogenetic", "Functional"))


rcp.list <- unique(dat$RCP)
thresh.list <- unique(dat$arbthresh)
for(rcp in rcp.list) {
  for(thresh in thresh.list) {
    
    # data for the main plots will use the most severe RCP and analog threshold of 20%
    dat.NvsD <- filter(dat, RCP==rcp, arbthresh == thresh) %>%
      group_by(x, y, output_srtm, measure, comm.type) %>%
      summarise(NoOfAnalogs = mean(value)) %>%
      spread(comm.type, NoOfAnalogs) %>%
      mutate(NminusD=Novel-Disappearing)
    
    plot.NvsD <- ggplot(NULL, aes(x, y)) + 
      geom_raster(data = dat.NvsD, aes(fill=NminusD)) +
      geom_raster(data = hdf, aes(alpha=layer)) +
      geom_path(data = ec.df, aes(x=long, y=lat)) +
      scale_fill_gradient2(name="# Novel analogs \n- # Disappearing \nanalogs") +
      guides(fill = guide_colorbar()) +
      facet_wrap(~measure, nrow=1) +
      scale_alpha(range = c(0, 0.5), guide = "none") +
      scale_x_continuous(name=expression(paste("Longitude (", degree, ")"))) + 
      scale_y_continuous(name=expression(paste("Latitude (", degree, ")"))) +
      coord_equal() + 
      theme(strip.background = element_blank(), panel.margin = unit(2, "lines"))
    
    dat.PvsF <- filter(dat, RCP==rcp, arbthresh==thresh) %>%
      group_by(x, y, output_srtm, measure, comm.type) %>%
      summarise(NoOfAnalogs=mean(value)) %>%
      spread(measure, NoOfAnalogs) %>%
      mutate(FminusP=Functional-Phylogenetic)
    
    plot.PvsF <- ggplot(NULL, aes(x, y)) + 
      geom_raster(data = dat.PvsF, aes(fill=FminusP)) +
      geom_raster(data = hdf, aes(alpha=layer)) +
      geom_path(data = ec.df, aes(x=long, y=lat)) +
      scale_fill_gradient2(name="# Functional analogs \n- # Phylogenetic analogs") +
      guides(fill = guide_colorbar()) +
      facet_wrap(~comm.type, nrow=1) +
      scale_alpha(range = c(0, 0.5), guide = "none") +
      scale_x_continuous(name=expression(paste("Longitude (", degree, ")"))) + 
      scale_y_continuous(name=expression(paste("Latitude (", degree, ")"))) +
      coord_equal() + 
      theme(strip.background = element_blank(), panel.margin = unit(2, "lines"))
    
    dat.novel <- filter(dat, RCP==rcp, arbthresh == thresh, comm.type=="Novel") %>%
      group_by(x, y, output_srtm, measure) %>%
      summarise(NoOfAnalogs = mean(value))
    
    plot.novel <- ggplot(NULL, aes(x, y)) + 
      geom_raster(data = dat.novel, aes(fill=NoOfAnalogs)) +
      geom_raster(data = hdf, aes(alpha=layer)) +
      geom_path(data = ec.df, aes(x=long, y=lat)) +
      scale_fill_gradient(low = "blue", high = "white", name="Number of analog\ncommunities") +
      guides(fill = guide_colorbar()) +
      scale_alpha(range = c(0, 0.5), guide = "none") +
      facet_wrap(~ measure, nrow=1) + 
      scale_x_continuous(name=expression(paste("Longitude (", degree, ")"))) + 
      scale_y_continuous(name=expression(paste("Latitude (", degree, ")"))) +
      coord_equal() + 
      theme(strip.background = element_blank(), panel.margin = unit(2, "lines"))
    
    # data for the main plots will use the most severe RCP and analog threshold of 20%
    dat.dis <- filter(dat, RCP==rcp, arbthresh == thresh, comm.type=="Disappearing") %>%
      group_by(x, y, output_srtm, measure) %>%
      summarise(NoOfAnalogs = mean(value))
    
    plot.dis <- ggplot(NULL, aes(x, y)) + 
      geom_raster(data = dat.dis, aes(fill=NoOfAnalogs)) +
      geom_raster(data = hdf, aes(alpha=layer)) +
      geom_path(data = ec.df, aes(x=long, y=lat)) +
      scale_fill_gradient(low = "red", high = "white", name="Number of analog\ncommunities") +
      guides(fill = guide_colorbar()) +
      scale_alpha(range = c(0, 0.5), guide = "none") +
      facet_wrap(~ measure, nrow=1) + 
      scale_x_continuous(name=expression(paste("Longitude (", degree, ")"))) + 
      scale_y_continuous(name=expression(paste("Latitude (", degree, ")"))) +
      coord_equal() + 
      theme(strip.background = element_blank(), panel.margin = unit(2, "lines"))
    
    gam.novel <- ggplot(dat.novel, aes(x=output_srtm, y=round(NoOfAnalogs, 0))) + geom_point(alpha=0.01) + 
      facet_wrap(~ measure) + 
      geom_smooth(method="gam",formula = y~s(x, k=20), colour="blue") +
      labs(x=expression("Elevation (m)"), y=expression("Number of analog\ncommunities")) +
      theme(strip.background=element_blank())
    
    gam.dis <- ggplot(dat.dis, aes(x=output_srtm, y=round(NoOfAnalogs, 0))) + geom_point(alpha=0.01) + 
      facet_wrap(~ measure) + 
      geom_smooth(method="gam",formula = y~s(x, k=20), colour="red") +
      labs(x=expression("Elevation (m)"), y=expression("Number of analog\ncommunities")) +
      theme(strip.background=element_blank())
    
    output.plot <- plot_grid(plot.novel, gam.novel, plot.dis, gam.dis, plot.NvsD, plot.PvsF, labels=c("A", "B", "C", "D", "E", "F"), ncol=2, align="h")
    save_plot(paste0("Figures/Main_Results_", rcp, "_", thresh, ".png"), output.plot, ncol=2, nrow=2, base_aspect_ratio = 1.3, base_width = 8.75, base_height = 4.66)
  }
}

# plot of the variance to check uncertainty
novel.20.sd <- filter(dat, comm.type=="Novel", arbthresh == 0.2, RCP=="85") %>%
  group_by(x, y, output_srtm, measure) %>%
  summarise(NoOfAnalogs = sd(value)/mean(value))

cvplot <- ggplot(NULL, aes(x, y)) + 
  geom_raster(data = novel.20.sd, aes(fill=NoOfAnalogs)) +
  geom_raster(data = hdf, aes(alpha=layer)) +
  geom_path(data = ec.df, aes(x=long, y=lat)) +
  scale_fill_gradient(low = "white", high = "blue", name="CV of # of analog\ncommunities") +
  guides(fill = guide_colorbar()) +
  scale_alpha(range = c(0, 0.5), guide = "none") +
  facet_wrap(~measure) + 
  scale_x_continuous(name=expression(paste("Longitude (", degree, ")"))) + 
  scale_y_continuous(name=expression(paste("Latitude (", degree, ")"))) +
  coord_equal() +
  theme(strip.background = element_blank(), panel.margin = unit(2, "lines"))

save_plot("Figures/Novel_by_RCP85_20perc_thres_cv.png", cvplot, base_aspect_ratio = 1.3, base_width = 8.75, base_height = 4.66)

# SENSITIVITY ANALYSIS FRIEDMAN TEST -------------------------------------------
# friedman test for each group
# input.dat <- spread(dat, arbthresh, value)
# res <- list()
# for(com in unique(input.dat$comm.type)) {
#   for(mes in unique(input.dat$measure)) {
#     dat <- filter(input.dat, comm.type==com, measure==mes)
#     f <- friedman.test(as.matrix(dat[8:11]))
#     res[[paste0(com, mes)]] <- cbind(tidy(f), measure=mes, comm.type=com)
#   }
# }
# 
# res <- do.call("rbind", res)

# Additional plots, not for the MS ---------------------------------------------
# CORRELATION OF SPECIES RICHNESS AND NUMBER OF ANALOGS ------------------------
# get the species counts for each cell
load("sppXsite/current.rda")
sp_richness <- data.frame(x=sppXsite$x, y=sppXsite$y, 
                          sp_richness=rowSums(sppXsite[,2:(ncol(sppXsite)-3)], na.rm=TRUE))

novel.20.rcp.sprich <- merge(novel.20.rcp, sp_richness)
cor.sprich <- group_by(novel.20.rcp.sprich, RCP, measure) %>%
  summarise(cor(sp_richness, NoOfAnalogs, method="spearman"))
novel.20.rcp.sprich <- merge(novel.20.rcp.sprich, cor.sprich)
names(novel.20.rcp.sprich)[8] <- "cor"
ggplot(novel.20.rcp.sprich, aes(x=sp_richness, y=NoOfAnalogs)) + geom_point() +
  facet_grid(RCP~measure) +
  geom_text(aes(label=paste("rho ==", round(cor, 2))), x=3, y=2000, parse=TRUE) +
  theme_classic()
ggsave("Figures/Novel_rich_analog_cor.png", width = 18, height = 9)

diss.20.rcp.sprich <- merge(diss.20.rcp, sp_richness)
cor.sprich <- group_by(diss.20.rcp.sprich, RCP, measure) %>%
  summarise(cor(sp_richness, NoOfAnalogs, method="spearman"))
diss.20.rcp.sprich <- merge(diss.20.rcp.sprich, cor.sprich)
names(diss.20.rcp.sprich)[8] <- "cor"
ggplot(diss.20.rcp.sprich, aes(x=sp_richness, y=NoOfAnalogs)) + geom_point() +
  facet_grid(RCP~measure) +
  geom_text(aes(label=paste("rho ==", round(cor, 2))), x=3, y=2000, parse=TRUE) +
  theme_classic()
ggsave("Figures/Diss_rich_analog_cor.png", width = 18, height = 9)

# UNCERTAINTY IN CLIMATE VARIABLES ---------------------------------------------
ann_mean_temp <- list.files("../worldclim_data/projections_2070/",
                            pattern = "701.tif$", full.names = TRUE, recursive = TRUE)
ann_mean_temp <- stack(ann_mean_temp)
ann_mean_temp <- stack(crop(ann_mean_temp, elev))

ann_prec <- list.files("../worldclim_data/projections_2070/",
                       pattern = "7012.tif$", full.names = TRUE, recursive = TRUE)
ann_prec <- stack(ann_prec)
ann_prec <- stack(crop(ann_prec, elev))

prec_seasonality <- list.files("../worldclim_data/projections_2070/",
                               pattern = "7015.tif$", full.names = TRUE, recursive = TRUE)
prec_seasonality <- stack(prec_seasonality)
prec_seasonality <- stack(crop(prec_seasonality, elev))


ann_mean_temp <- as.data.frame(ann_mean_temp, xy = TRUE) %>%
  gather(key, value, -x, -y, na.rm = TRUE) %>%
  mutate(var="ann_mean_temp", 
         GCM=substr(variable, 1, 2),
         RCP=substr(variable, 3, 4)) %>%
  group_by(x, y, var, RCP) %>%
  summarise(CV=sd(value)/mean(value))

ann_prec <- as.data.frame(ann_prec, xy = TRUE) %>%
  gather(key, value, -x, -y, na.rm = TRUE) %>%
  mutate(var="ann_prec", 
         GCM=substr(variable, 1, 2),
         RCP=substr(variable, 3, 4)) %>%
  group_by(x, y, var, RCP) %>%
  summarise(CV=sd(value)/mean(value))

prec_seasonality <- as.data.frame(prec_seasonality, xy = TRUE) %>%
  gather(key, value, -x, -y, na.rm = TRUE) %>%
  mutate(var="prec_seasonality", 
         GCM=substr(variable, 1, 2),
         RCP=substr(variable, 3, 4)) %>%
  group_by(x, y, var, RCP) %>%
  summarise(CV=sd(value)/mean(value))

climate.CV <- rbind(ann_mean_temp, ann_prec, prec_seasonality)

ggplot(NULL, aes(x, y)) + 
  geom_raster(data = climate.CV, aes(fill=CV)) +
  geom_raster(data = hdf, aes(alpha=layer)) +
  geom_path(data = ec.df, aes(x=long, y=lat)) +
  scale_fill_gradient(low = "white", high = "red", name="CV for climate variable") +
  guides(fill = guide_colorbar()) +
  scale_alpha(range = c(0, 0.5), guide = "none") +
  facet_grid(var ~ RCP) + 
  scale_x_continuous(name=expression(paste("Longitude (", degree, ")"))) + 
  scale_y_continuous(name=expression(paste("Latitude (", degree, ")"))) +
  coord_equal() + theme_classic(base_size=15) + 
  theme(strip.background = element_blank(), panel.margin = unit(2, "lines"))

ggsave("Figures/climate_CV.png", width=9, height=9)

ann_temp.p <- ggplot(NULL, aes(x, y)) + 
  geom_raster(data = ann_mean_temp, aes(fill=CV)) +
  geom_raster(data = hdf, aes(alpha=layer)) +
  geom_path(data = ec.df, aes(x=long, y=lat)) +
  scale_fill_gradient(low = "white", high = "red", name="CV mean annual temperature") +
  guides(fill = guide_colorbar()) +
  scale_alpha(range = c(0, 0.5), guide = "none") +
  facet_grid(. ~ RCP) + 
  scale_x_continuous(name=expression(paste("Longitude (", degree, ")"))) + 
  scale_y_continuous(name=expression(paste("Latitude (", degree, ")"))) +
  coord_equal() + theme_classic(base_size=15) + 
  theme(strip.background = element_blank(), panel.margin = unit(2, "lines"))

ann_prec.p <- ggplot(NULL, aes(x, y)) + 
  geom_raster(data = ann_prec, aes(fill=SD)) +
  geom_raster(data = hdf, aes(alpha=layer)) +
  geom_path(data = ec.df, aes(x=long, y=lat)) +
  scale_fill_gradient(low = "white", high = "red", name="SD annual precipitation") +
  guides(fill = guide_colorbar()) +
  scale_alpha(range = c(0, 0.5), guide = "none") +
  facet_grid(. ~ RCP) + 
  scale_x_continuous(name=expression(paste("Longitude (", degree, ")"))) + 
  scale_y_continuous(name=expression(paste("Latitude (", degree, ")"))) +
  coord_equal() + theme_classic(base_size=15) + 
  theme(strip.background = element_blank(), panel.margin = unit(2, "lines"))

prec_seasonality.p <- ggplot(NULL, aes(x, y)) + 
  geom_raster(data = prec_seasonality, aes(fill=SD)) +
  geom_raster(data = hdf, aes(alpha=layer)) +
  geom_path(data = ec.df, aes(x=long, y=lat)) +
  scale_fill_gradient(low = "white", high = "red", name="SD precipitation seasonality") +
  guides(fill = guide_colorbar()) +
  scale_alpha(range = c(0, 0.5), guide = "none") +
  facet_grid(. ~ RCP) + 
  scale_x_continuous(name=expression(paste("Longitude (", degree, ")"))) + 
  scale_y_continuous(name=expression(paste("Latitude (", degree, ")"))) +
  coord_equal() + theme_classic(base_size=15) + 
  theme(strip.background = element_blank(), panel.margin = unit(2, "lines"))

png("Figures/climate_variable_sd.png", width = 864, height = 864)
grid.arrange(ann_temp.p, ann_prec.p, prec_seasonality.p)
dev.off()



