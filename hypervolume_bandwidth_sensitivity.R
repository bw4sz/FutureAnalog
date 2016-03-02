#load results for current scenario
cur <- list.files("sppXsite", pattern="current", full.names = TRUE)

# get trait data
traits <- getTraitData()
traits.pca <- prcomp(traits)$x # traits too correlated to go straight into hypervolume

load(cur)
sppXsite <- sppXsite[,2:(ncol(sppXsite)-3)]

comm.hv <- sapply(rownames(sppXsite), function(k){
  g <- sppXsite[k,]
  sp <- names(g[which(g==1)])
  traits.sub <- subset(traits.pca, rownames(traits) %in% sp)
  cv.strt <- Sys.time()
  cv.bw <- estimate_bandwidth(traits.pca, method="cross-validation")
  cv.time <- Sys.time() - bw.strt
  cv.hv <- get_volume(hypervolume(traits.sub, bandwidth = cv.bw))
  pi.strt <- Sys.time()
  pi.bw <- estimate_bandwidth(traits.pca, method="plug-in")
  pi.time <- Sys.time() - bw.strt
  pi.hv <- get_volume(hypervolume(traits.sub, bandwidth = pi.bw))
  sl.strt <- Sys.time()
  sl.bw <- estimate_bandwidth(traits.pca, method="silverman")
  sl.time <- Sys.time() - bw.strt
  sl.hv <- get_volume(hypervolume(traits.sub, bandwidth = sl.bw))
  no.species <- nrow(traits.sub)
  out <- c(cv.time, cv.hv, pi.time, pi.hv, sl.time, sl.hv, no.species)
  return(out)
}) 

comm.hv <- t(comm.hv)
colnames(comm.hv) <- c("cv.time", "cv.hv", "pi.time", "pi.hv", "sl.time", "sl.hv", "no.species")

sl.cv.diff <- comm.hv[,2] - comm.hv[,6]
sl.pi.diff <- comm.hv[,4] - comm.hv[,6]
cv.pi.diff <- comm.hv[,4] - comm.hv[,2]

diffs <- data.frame(sl.cv=sl.cv.diff, sl.pi=sl.pi.diff, cv.pi=cv.pi.diff)
ggpairs(comm.hv, columns=c(2,4,6))
hist(sl.cv.diff)

out <-comm.hv[,c(2,4,6,7)]
save(out, file="hypervolume_sensitivity_analysis.rda")
