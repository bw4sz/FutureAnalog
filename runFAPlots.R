# GET DATA FOR THE CURRENT/FUTURE COMPARISON -----------------------------------


rasterToDataFrame <- function(out_path){ 
  
  out.rasters <- list.files(out_path, pattern = "NonAnalogRasters_", full.names = TRUE, recursive = TRUE)
  
  out <- lapply(out.rasters, function(x) { 
    load(x)
    
    results <- stack(get("results"), elev)
    results <- mask(results, ec)
    arbthresh <- str_match(x,pattern=paste(cell,"(\\w+.\\w+)/NonAnalog",sep="/"))[,2] 
    GCM <- substr(str_match(x,pattern="NonAnalogRasters_(\\w+.\\w+)bi70")[,2], 1, 2) 
    RCP <- substr(str_match(x,pattern="NonAnalogRasters_(\\w+.\\w+)bi70")[,2], 3, 4)  
    
    dat <- rasterToPoints(results) %>%
      data.frame() %>%
      gather(key, value, -x, -y, -output_srtm, na.rm = TRUE) %>%
      mutate(arbthresh = arbthresh, GCM = GCM, RCP = RCP) %>%
      separate(key, into = c("comm.type", "measure"), sep = "\\.")
  } 
  ) 
  out <- do.call("rbind", out)
}

# hillshade data
ec <- readOGR("InputData", "EcuadorCut")
load("InputData/srtm_5arcmin.rda")
hill <- raster("MAPA_VEGETACION_MAE_FINAL_OFICIAL/hillshade_cropped.tif")
hdf <- as.data.frame(hill, xy=TRUE)
hdf <- data.frame(hdf)
hdf <- hdf[complete.cases(hdf),]
colnames(hdf) <- c("x", "y", "layer")

# Get results data ----
dat <- rasterToDataFrame(out_path)

dat$measure <- factor(dat$measure, levels=c("Tax", "Phylo", "Func"),
                      labels=c("Taxonomic", "Phylogenetic", "Functional"))


rcp.list <- unique(dat$RCP)
thresh.list <- unique(dat$arbthresh)

# Figure 2 is a conceptual diagram and created elsewhere ----

# Figure 3 (and supp mat): main results at each RCP and threshold ----
for(rcp in rcp.list) {
  for(thresh in thresh.list) {
    
    # data for the main plots will use the most severe RCP and analog threshold of 20%
    dat.NvsD <- filter(dat, RCP==rcp, arbthresh == thresh) %>%
      group_by(x, y, output_srtm, measure, comm.type) %>%
      summarise(NoOfAnalogs = mean(value)) %>%
      spread(comm.type, NoOfAnalogs) %>%
      mutate(NminusD=Novel-Disappearing)
    
    minval <- min(dat.NvsD$NminusD) + 0.3*min(dat.NvsD$NminusD)
    maxval <- max(dat.NvsD$NminusD) - 0.3*max(dat.NvsD$NminusD)
    maxval <- max(abs(minval), maxval)
    
    plot.NvsD <- ggplot(NULL, aes(x, y)) + 
      geom_raster(data = dat.NvsD, aes(fill=NminusD)) +
      scale_fill_gradient2(name="Difference in\nno. analog\ncommunities", breaks=c(-maxval, 0, maxval), 
                           labels=c("Future", "No difference", "Current"),
                           limits=c(-maxval, maxval)) +
      facet_wrap(~measure, nrow=1) +
      geom_raster(data = hdf, aes(alpha=layer)) +
      scale_alpha(range = c(0, 0.5), guide = "none") +
      scale_x_continuous(name="") + 
      scale_y_continuous(name="") +
      coord_equal() + 
      theme(strip.background = element_blank(), panel.margin = unit(2, "lines"), 
            axis.ticks=element_blank(), axis.text=element_blank(), axis.line=element_blank())
    
    dat.PvsF <- filter(dat, RCP==rcp, arbthresh==thresh) %>%
      group_by(x, y, output_srtm, measure, comm.type) %>%
      summarise(NoOfAnalogs=mean(value)) %>%
      spread(measure, NoOfAnalogs) %>%
      mutate(FminusP=Functional-Phylogenetic)
    
    dat.PvsF$comm.type <- factor(dat.PvsF$comm.type, levels=c("Novel", "Disappearing"), labels=c("Current", "Future"))
    minval <- min(dat.PvsF$FminusP) + 0.3*min(dat.PvsF$FminusP)
    maxval <- max(dat.PvsF$FminusP) - 0.3*max(dat.PvsF$FminusP)
    maxval <- max(abs(minval), maxval)
    
    plot.PvsF <- ggplot(NULL, aes(x, y)) + 
      geom_raster(data = dat.PvsF, aes(fill=FminusP)) +
      scale_fill_gradient2(name="No. analog\ncommunities", breaks=c(-maxval, 0, maxval), 
                           labels=c("Phylogenetic", "No difference", "Functional"),
                           limits=c(-maxval, maxval)) +
      facet_wrap(~comm.type, nrow=1) +
      geom_raster(data = hdf, aes(alpha=layer)) +
      scale_alpha(range = c(0, 0.5), guide = "none") +
      scale_x_continuous(name="") + 
      scale_y_continuous(name="") +
      coord_equal() + 
      theme(strip.background = element_blank(), panel.margin = unit(2, "lines"), 
            axis.ticks=element_blank(), axis.text=element_blank(), axis.line=element_blank())
    
    
    dat.novel <- filter(dat, RCP==rcp, arbthresh == thresh, comm.type=="Novel") %>%
      group_by(x, y, output_srtm, measure) %>%
      summarise(NoOfAnalogs = mean(value))
    
    dat.novel.zero <- filter(dat.novel, NoOfAnalogs==0)
    
    plot.novel <- ggplot(NULL) + 
      geom_raster(data = dat.novel, aes(x = x, y = y, fill=NoOfAnalogs)) +
      geom_raster(data = dat.novel.zero, color = "black", aes(x = x, y = y)) +
      scale_fill_gradient2(name="No. current\nanalogs") +
      facet_wrap(~measure, nrow=1) +
      geom_raster(data = hdf, aes(x = x, y = y, alpha=layer)) +
      scale_alpha(range = c(0, 0.5), guide = "none") +
      scale_x_continuous(name="") + 
      scale_y_continuous(name="") +
      coord_equal() + 
      theme(strip.background = element_blank(), panel.margin = unit(2, "lines"), 
            axis.ticks=element_blank(), axis.text=element_blank(), axis.line=element_blank())
    
    # data for the main plots will use the most severe RCP and analog threshold of 20%
    dat.dis <- filter(dat, RCP==rcp, arbthresh == thresh, comm.type=="Disappearing") %>%
      group_by(x, y, output_srtm, measure) %>%
      summarise(NoOfAnalogs = mean(value))
    
    dat.dis.zero <- filter(dat.dis, NoOfAnalogs==0)
    plot.dis <- ggplot(NULL, aes(x, y)) + 
      geom_raster(data = dat.dis, aes(fill=NoOfAnalogs)) +
      geom_raster(data = dat.dis.zero, color = "black") +
      scale_fill_gradient2(low="white", high="red", name="No. future\nanalogs") +
      facet_wrap(~measure, nrow=1) +
      geom_raster(data = hdf, aes(alpha=layer)) +
      scale_alpha(range = c(0, 0.5), guide = "none") +
      scale_x_continuous(name="") + 
      scale_y_continuous(name="") +
      coord_equal() + 
      theme(strip.background = element_blank(), panel.margin = unit(2, "lines"), 
            axis.ticks=element_blank(), axis.text=element_blank(), axis.line=element_blank())
    
    gam.novel <- ggplot(dat.novel, aes(x=output_srtm, y=round(NoOfAnalogs, 0))) + geom_point(alpha=0.01) + 
      facet_wrap(~ measure) + 
      geom_smooth(method="gam",formula = y~s(x, k=20), colour="blue") +
      labs(x=expression("Elevation (m)"), y=expression("No. current analogs")) +
      ylim(0, 3000) +
      theme(strip.background=element_blank())
    
    gam.dis <- ggplot(dat.dis, aes(x=output_srtm, y=round(NoOfAnalogs, 0))) + geom_point(alpha=0.01) + 
      facet_wrap(~ measure) + 
      geom_smooth(method="gam",formula = y~s(x, k=20), colour="red") +
      labs(x=expression("Elevation (m)"), y=expression("No. future analogs")) +
      ylim(0, 3000) +
      theme(strip.background=element_blank())
    
    output.plot <- plot_grid(plot.novel, gam.novel, plot.dis, gam.dis, plot.NvsD, plot.PvsF, labels=c("A", "B", "C", "D", "E", "F"), ncol=2, align="h")
    save_plot(paste0("Figures/Main_Results_", rcp, "_", thresh, ".png"), output.plot, ncol=2, nrow=2, base_aspect_ratio = 1.3, base_width = 8.75, base_height = 4.66)
  }
}

# plot of the variance to check uncertainty
cvplot.dat <- filter(dat, arbthresh == 0.2, RCP=="85") %>%
  group_by(x, y, output_srtm, measure, comm.type) %>%
  summarise(NoOfAnalogs = sd(value)/mean(value))

cvplot.dat$comm.type <- factor(cvplot.dat$comm.type, levels=c("Novel", "Disappearing"), labels=c("Current", "Future"))

cvplot <- ggplot(NULL, aes(x, y)) + 
  geom_raster(data = cvplot.dat, aes(fill=NoOfAnalogs)) +
  geom_raster(data = hdf, aes(alpha=layer)) +
  scale_fill_gradient(low = "white", high = "darkgreen", name="CV of no. of analog\ncommunities") +
  guides(fill = guide_colorbar()) +
  scale_alpha(range = c(0, 0.5), guide = "none") +
  facet_grid(comm.type~measure) + 
  scale_x_continuous(name="") + 
  scale_y_continuous(name="") +
  coord_equal() + 
  theme(strip.background = element_blank(), panel.margin = unit(2, "lines"), 
        axis.ticks=element_blank(), axis.text=element_blank(), axis.line=element_blank())

save_plot("Figures/RCP85_20perc_thres_cv.png", cvplot, base_aspect_ratio = 1.3, base_width = 8.75, base_height = 4.66)

# UNCERTAINTY IN CLIMATE VARIABLES ---------------------------------------------
ann_mean_temp <- list.files("../worldclim_data/projections_2070/",
                            pattern = "701.tif$", full.names = TRUE, recursive = TRUE)
ann_mean_temp <- stack(ann_mean_temp)
ann_mean_temp <- mask(ann_mean_temp, ec)

ann_prec <- list.files("../worldclim_data/projections_2070/",
                       pattern = "7012.tif$", full.names = TRUE, recursive = TRUE)
ann_prec <- stack(ann_prec)
ann_prec <- mask(ann_prec, ec)

prec_seasonality <- list.files("../worldclim_data/projections_2070/",
                               pattern = "7015.tif$", full.names = TRUE, recursive = TRUE)
prec_seasonality <- stack(prec_seasonality)
prec_seasonality <- mask(prec_seasonality, ec)


ann_mean_temp <- as.data.frame(ann_mean_temp, xy = TRUE) %>%
  gather(key, value, -x, -y, na.rm = TRUE) %>%
  mutate(var="ann_mean_temp", 
         GCM=substr(key, 1, 2),
         RCP=substr(key, 3, 4)) %>%
  filter(RCP=="85") %>%
  group_by(x, y, var) %>%
  summarise(CV=sd(value)/mean(value))

ann_prec <- as.data.frame(ann_prec, xy = TRUE) %>%
  gather(key, value, -x, -y, na.rm = TRUE) %>%
  mutate(var="ann_prec", 
         GCM=substr(key, 1, 2),
         RCP=substr(key, 3, 4)) %>%
  filter(RCP=='85') %>%
  group_by(x, y, var) %>%
  summarise(CV=sd(value)/mean(value))

prec_seasonality <- as.data.frame(prec_seasonality, xy = TRUE) %>%
  gather(key, value, -x, -y, na.rm = TRUE) %>%
  mutate(var="prec_seasonality", 
         GCM=substr(key, 1, 2),
         RCP=substr(key, 3, 4)) %>%
  filter(RCP=='85') %>%
  group_by(x, y, var) %>%
  summarise(CV=sd(value)/mean(value))

climate.CV <- rbind(ann_mean_temp, ann_prec, prec_seasonality)

ggplot(NULL, aes(x, y)) + 
  geom_raster(data = climate.CV, aes(fill=CV)) +
  geom_raster(data = hdf, aes(alpha=layer)) +
  scale_fill_gradient(low = "white", high = "darkgreen", name="CV for climate variable") +
  guides(fill = guide_colorbar()) +
  scale_alpha(range = c(0, 0.5), guide = "none") +
  facet_wrap(~var) + 
  scale_x_continuous(name=expression(paste("Longitude (", degree, ")"))) + 
  scale_y_continuous(name=expression(paste("Latitude (", degree, ")"))) +
  coord_equal() + theme_classic(base_size=15) + 
  theme(strip.background = element_blank(), panel.margin = unit(2, "lines"))

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            ggsave("Figures/climate_CV.png", width=9, height=3)

ann_temp.p <- ggplot(NULL, aes(x, y)) + 
  geom_raster(data = ann_mean_temp, aes(fill=CV)) +
  geom_raster(data = hdf, aes(alpha=layer)) +
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



