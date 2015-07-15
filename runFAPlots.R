# GET DATA FOR THE CURRENT/FUTURE COMPARISON -----------------------------------
rasterToDataFrame <- function(out_path){
  load("InputData/srtm_5arcmin.rda")
  out.rasters <- list.files(out_path, pattern = "NonAnalogRasters", full.names = TRUE, recursive = TRUE)
  
  out <- lapply(out.rasters, function(x) {
    load(x)
    results <- stack(results, elev)
    arbthresh <- str_match(x,pattern=paste(cell,"(\\w+.\\w+)/NonAnalog",sep="/"))[,2]
    GCM <- substr(str_match(x,pattern=paste("Rasters_(\\w+.\\w+)bi70",sep="/"))[,2], 1, 2)
    RCP <- substr(str_match(x,pattern=paste("Rasters_(\\w+.\\w+)bi70",sep="/"))[,2], 3, 4)
    
    dat <- as.data.frame(results, xy = TRUE) %>%
      gather(key, value, -x, -y, -output_srtm, na.rm = TRUE) %>%
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
                  labels=c("RCP2.6", "RCP4.5", "RCP8.5"))

dat$GCM <- factor(dat$GCM, levels=unique(dat$GCM),
                  labels=c("CCSM4", "CNRM-CM5", "GISS-E2-R", "HadGEM2-ES", 
                           "IPSL-CM5A-LR", "MIROC5", "MPI-ESM-LR"))

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

# FIGURE 1 WILL BE A CONCEPTUAL DIAGRAM TYPE THING -----------------------------

# FIGURE 2 NOVEL COMMUNITIES BY RCP --------------------------------------------
novel.20.rcp <- filter(dat, comm.type=="Novel", arbthresh == 0.2) %>%
  group_by(x, y, output_srtm, measure, RCP) %>%
  summarise(NoOfAnalogs = mean(value))

novel.20.rcp.summary <- group_by(novel.20.rcp, measure, RCP) %>%
  summarise(min.analog=min(NoOfAnalogs),
            max.analog=max(NoOfAnalogs),
            mean.analog=mean(NoOfAnalogs),
            sd.analog=sd(NoOfAnalogs),
            no.analogs=sum(NoOfAnalogs==0))

ggplot(NULL, aes(x, y)) + 
  geom_raster(data = novel.20.rcp, aes(fill=NoOfAnalogs)) +
  geom_raster(data = hdf, aes(alpha=layer)) +
  geom_path(data = ec.df, aes(x=long, y=lat)) +
  scale_fill_gradient(low = "blue", high = "white", name="Number of analog\ncommunities") +
  guides(fill = guide_colorbar()) +
  scale_alpha(range = c(0, 0.5), guide = "none") +
  facet_grid(measure ~ RCP) + 
  scale_x_continuous(name=expression(paste("Longitude (", degree, ")"))) + 
  scale_y_continuous(name=expression(paste("Latitude (", degree, ")"))) +
  coord_equal() + theme_classic(base_size=15) + 
  theme(strip.background = element_blank(), panel.margin = unit(2, "lines"))

ggsave("Figures/Novel_by_RCP_20perc_thres.png", width = 9, height = 9)

# 2b Analog ~ Elevation GAMs
ggplot(novel.20.rcp, aes(x=output_srtm, y=round(NoOfAnalogs, 0))) + geom_point(alpha=0.01) + 
  facet_grid(measure~RCP) + 
  geom_smooth(method="gam",formula = y~s(x, k=20)) +
  labs(x="Elevation", y="Number of analog communities") +
  theme_classic(base_size=15) + theme(strip.background=element_blank())

ggsave("Figures/Novel_elevation_gam_RCP_20perc.png", width = 9, height = 9)

# FIGURE 3 DISAPPEARING COMMUNITIES BY RCP -------------------------------------
diss.20.rcp <- filter(dat, comm.type=="Disappearing", arbthresh == 0.2) %>%
  group_by(x, y, output_srtm, measure, RCP) %>%
  summarise(NoOfAnalogs = mean(value))

diss.20.rcp.summary <- group_by(diss.20.rcp, measure, RCP) %>%
  summarise(min.analog=min(NoOfAnalogs),
            max.analog=max(NoOfAnalogs),
            mean.analog=mean(NoOfAnalogs),
            sd.analog=sd(NoOfAnalogs),
            no.analogs=sum(NoOfAnalogs==0))

ggplot(NULL, aes(x, y)) + 
  geom_raster(data = diss.20.rcp, aes(fill=NoOfAnalogs)) +
  geom_raster(data = hdf, aes(alpha=layer)) +
  geom_path(data = ec.df, aes(x=long, y=lat)) +
  scale_fill_gradient(low = "red", high = "white", name="Number of analog\ncommunities") +
  guides(fill = guide_colorbar()) +
  scale_alpha(range = c(0, 0.5), guide = "none") +
  facet_grid(measure ~ RCP) + 
  scale_x_continuous(name=expression(paste("Longitude (", degree, ")"))) + 
  scale_y_continuous(name=expression(paste("Latitude (", degree, ")"))) +
  coord_equal() + theme_classic(base_size=15) + 
  theme(strip.background = element_blank(), panel.margin = unit(2, "lines"))

ggsave("Figures/Disappearing_by_RCP_20perc_thres.png", width = 9, height = 9)

# 3b Analog ~ Elevation GAMs
ggplot(diss.20.rcp, aes(x=output_srtm, y=round(NoOfAnalogs, 0))) + geom_point(alpha=0.01) + 
  facet_grid(measure~RCP) + 
  geom_smooth(method="gam",formula = y~s(x, k=20), colour="red") +
  labs(x="Elevation", y="Number of analog communities") +
  theme_classic(base_size=15) + theme(strip.background=element_blank())

ggsave("Figures/Diss_elevation_gam_RCP_20perc.png", width = 9, height = 9)

# FIGURE 4 BOXPLOTS OF DIFFERENCES IN NUMBER OF ANALOGUES BETWEEN SCENARIOS-----

# function to create all pairwise differences
# pair.diff <- function(input){
#   nm1 <- outer(colnames(input), colnames(input), paste, sep=" - ")
#   indx1 <-  which(lower.tri(nm1, diag=TRUE))
#   res <- outer(1:ncol(input), 1:ncol(input), 
#                function(x,y) input[,x]-input[,y])
#   colnames(res) <- nm1
#   res <- res[-indx1]
#   return(res)
# }
# 
# input.dat <- filter(dat, arbthresh == 0.2) %>%
#   group_by(x, y, measure, comm.type, RCP) %>%
#   summarise(NoOfAnalogs = mean(value)) %>%
#   spread(RCP, NoOfAnalogs)
# 
# input <- input.dat[,5:7]
# input.diff <- pair.diff(input)
# input.diff <- cbind(comm.type=input.dat$comm.type, measure = input.dat$measure, input.diff)
# input.diff <- gather(input.diff, var, value, -comm.type, -measure)
# input.diff$comm.type <- factor(input.diff$comm.type, levels=c("Novel", "Disappearing"))
# input.diff$variable <- factor(input.diff$variable, levels=c("RCP2.6 - RCP4.5",
#                                                             "RCP4.5 - RCP8.5",
#                                                             "RCP2.6 - RCP8.5"))
# ggplot(input.diff, aes(x=variable, y = value)) + 
#   geom_boxplot() + 
#   labs(x="", y="Difference in number of analog communities") +
#   geom_hline(colour="red") + 
#   facet_grid(measure~comm.type) + 
#   theme_classic() + 
#   theme(strip.background=element_blank())
# 
# # for GCM
# input.dat <- filter(dat, arbthresh == 0.2) %>%
#   group_by(x, y, measure, comm.type, GCM) %>%
#   summarise(NoOfAnalogs = mean(value)) %>%
#   spread(GCM, NoOfAnalogs)
# 
# input <- input.dat[,5:11]
# input.diff <- pair.diff(input)
# input.diff <- cbind(comm.type=input.dat$comm.type, measure = input.dat$measure, input.diff)
# input.diff <- group_by(input.diff, comm.type, measure) %>%
#   summarise_each(funs(mean(., na.rm=TRUE)))
# 
# 
# input.diff <- gather(input.diff, var, value, -comm.type, -measure)
# input.diff$comm.type <- factor(input.diff$comm.type, levels=c("Novel", "Disappearing"))
# 
# ggplot(input.diff, aes(x=variable, y = value)) + 
#   geom_boxplot() + 
#   labs(x="", y="Difference in number of analog communities") +
#   geom_hline(colour="red") + 
#   facet_grid(measure~comm.type) + 
#   theme_classic() + 
#   theme(strip.background=element_blank()) + 
#   theme(axis.text.x=element_text(angle=-90))
# 
# input <- as.matrix(input.dat[,5:11])
# input.diff <- list(data.frame(var = "RCP2.6 - RCP4.5", value = input[,1] - input[,2]),
#                    data.frame(var = "RCP4.5 - RCP8.5", value = input[,2] - input[,3]),
#                    data.frame(var = "RCP2.6 - RCP8.5", value = input[,1] - input[,3]))
# input.diff <- do.call("rbind", input.diff)
# input.diff <- cbind(comm.type=input.dat$comm.type, measure = input.dat$measure, input.diff)
# input.diff$comm.type <- factor(input.diff$comm.type, levels=c("Novel", "Disappearing"))
# ggplot(input.diff, aes(x=var, y = value)) + 
#   geom_boxplot() + 
#   labs(x="", y="Difference in number of analog communities") +
#   geom_hline(colour="red") + 
#   facet_grid(measure~comm.type) + 
#   theme_classic() + 
#   theme(strip.background=element_blank())
# 
# # friedman test for each group
# res <- list()
# for(com in unique(input.dat$comm.type)) {
#   for(mes in unique(input.dat$measure)) {
#     dat <- filter(input.dat, comm.type==com, measure==mes)
#     f <- friedman.test(as.matrix(dat[5:11]))
#     res[[paste0(com, mes)]] <- cbind(tidy(f), measure=mes, comm.type=com)
#   }
# }
# 
# res <- do.call("rbind", res)



# SUPP. MAT. FIGURE X (THRESHOLD SENSITIVITY NOVEL) ----------------------------
thresh_sa.novel <- filter(dat, comm.type=="Novel") %>% 
  group_by(x, y, measure, arbthresh) %>%
  summarise(NoOfAnalogs = mean(value))

ggplot(NULL, aes(x, y)) + 
  geom_raster(data = thresh_sa.novel, aes(fill=NoOfAnalogs)) +
  geom_raster(data = hdf, aes(alpha=layer)) +
  geom_path(data = ec.df, aes(x=long, y=lat)) +
  scale_fill_gradient(low = "blue", high = "white", name="Number of analog\ncommunities") +
  guides(fill = guide_colorbar()) +
  scale_alpha(range = c(0, 0.5), guide = "none") +
  facet_grid(measure ~ arbthresh) + 
  scale_x_continuous(name=expression(paste("Longitude (", degree, ")"))) + 
  scale_y_continuous(name=expression(paste("Latitude (", degree, ")"))) +
  coord_equal() + theme_classic(base_size=15) + 
  theme(strip.background = element_blank(), panel.margin = unit(2, "lines"))

ggsave("Figures/Threshold_SA_Novel.png", width = 17, height = 9)

# SUPP. MAP. FIGURE X (THRESHOLD SENSITIVITY DISAPPEARING) ---------------------
thresh_sa.diss <- filter(dat, comm.type=="Disappearing") %>% 
  group_by(x, y, measure, arbthresh) %>%
  summarise(NoOfAnalogs = mean(value))

ggplot(NULL, aes(x, y)) + 
  geom_raster(data = thresh_sa.diss, aes(fill=NoOfAnalogs)) +
  geom_raster(data = hdf, aes(alpha=layer)) +
  geom_path(data = ec.df, aes(x=long, y=lat)) +
  scale_fill_gradient(low = "red", high = "white", name="Number of analog\ncommunities") +
  guides(fill = guide_colorbar()) +
  scale_alpha(range = c(0, 0.5), guide = "none") +
  facet_grid(measure ~ arbthresh) + 
  scale_x_continuous(name=expression(paste("Longitude (", degree, ")"))) + 
  scale_y_continuous(name=expression(paste("Latitude (", degree, ")"))) +
  coord_equal() + theme_classic(base_size=15) + 
  theme(strip.background = element_blank(), panel.margin = unit(2, "lines"))


ggsave("Figures/Threshold_SA_Disappearing.png", width = 17, height = 9)

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
# SUPP. MAT. FIGURE X (NOVEL COMMUNITIES BY GCM) -------------------------------
novel.20.gcm <- filter(dat, comm.type=="Novel", arbthresh == 0.2) %>%
  group_by(x, y, measure, GCM) %>%
  summarise(NoOfAnalogs = mean(value))

novel.20.gcm.summary <- group_by(novel.20.gcm, measure, GCM) %>%
  summarise(min.analog=min(NoOfAnalogs),
            max.analog=max(NoOfAnalogs),
            mean.analog=mean(NoOfAnalogs),
            sd.analog=sd(NoOfAnalogs),
            no.analogs=sum(NoOfAnalogs==0))

ggplot(NULL, aes(x, y)) + 
  geom_raster(data = novel.20.gcm, aes(fill=NoOfAnalogs)) +
  geom_raster(data = hdf, aes(alpha=layer)) +
  geom_path(data = ec.df, aes(x=long, y=lat)) +
  scale_fill_gradient(low = "blue", high = "white", name="Number of analog\ncommunities") +
  guides(fill = guide_colorbar()) +
  scale_alpha(range = c(0, 0.5), guide = "none") +
  facet_grid(measure ~ GCM) + 
  scale_x_continuous(name=expression(paste("Longitude (", degree, ")"))) + 
  scale_y_continuous(name=expression(paste("Latitude (", degree, ")"))) +
  coord_equal() + theme_classic(base_size=15) + 
  theme(strip.background = element_blank(), panel.margin = unit(2, "lines"))

ggsave("Figures/Novel_by_GCM_20perc_thres.png", width = 18, height = 9)

# SUPP. MAT. FIGURE X (DISAPPEARING COMMUNITIES BY GCM) ------------------------
diss.20.gcm <- filter(dat, comm.type=="Disappearing", arbthresh == 0.2) %>%
  group_by(x, y, measure, GCM) %>%
  summarise(NoOfAnalogs = mean(value))

diss.20.gcm.summary <- group_by(diss.20.gcm, measure, GCM) %>%
  summarise(min.analog=min(NoOfAnalogs),
            max.analog=max(NoOfAnalogs),
            mean.analog=mean(NoOfAnalogs),
            sd.analog=sd(NoOfAnalogs),
            no.analogs=sum(NoOfAnalogs==0))

ggplot(NULL, aes(x, y)) + 
  geom_raster(data = diss.20.gcm, aes(fill=NoOfAnalogs)) +
  geom_raster(data = hdf, aes(alpha=layer)) +
  geom_path(data = ec.df, aes(x=long, y=lat)) +
  scale_fill_gradient(low = "red", high = "white", name="Number of analog\ncommunities") +
  guides(fill = guide_colorbar()) +
  scale_alpha(range = c(0, 0.5), guide = "none") +
  facet_grid(measure ~ GCM) + 
  scale_x_continuous(name=expression(paste("Longitude (", degree, ")"))) + 
  scale_y_continuous(name=expression(paste("Latitude (", degree, ")"))) +
  coord_equal() + theme_classic(base_size=15) + 
  theme(strip.background = element_blank(), panel.margin = unit(2, "lines"))

ggsave("Figures/Disappearing_by_GCM_20perc_thres.png", width = 18, height = 9)
