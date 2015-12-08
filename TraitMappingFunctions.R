runSppXsite <- function(out_path){
  # select only well fitting species (currently set as all where TSS is over 0.5
  # for all models and ROC is over 7.5 for all models)
  model_eval<-list.files(out_path, full.name=TRUE,recursive=T,pattern="Eval.csv")
  model_eval<-rbind_all(lapply(model_eval, 
                               function(x) read.csv(x, stringsAsFactors = FALSE)))
  colnames(model_eval)[1:2] <- c("Stat", "Species")
  
  # change here for sensitivity analyis
  model_eval$TEST <- (apply(model_eval, 1, function(x) min(x[3:5], na.rm=TRUE) < 0.5) 
                      & model_eval$Stat == "TSS") |
    (apply(model_eval, 1, function(x) min(x[3:5], na.rm=TRUE) < 0.75) 
     & model_eval$Stat == "ROC")
  
  well.fitting.models <- subset(model_eval, !TEST)
  well.fitting.species <- unique(well.fitting.models$Species)
  well.fitting.species <- gsub(" ", ".", well.fitting.species)
  
  # get the crop files
  niche.crops <- list.files(out_path,pattern="crop.gri",full.name=T,recursive=T)
  
  files <- c()
  for(s in well.fitting.species){
    files <- c(files,grep(s, niche.crops))
  }
  
  niche.crops <- niche.crops[files]
  
  # create folder to store results
  if(!dir.exists(paste(out_path, "sppXsite", sep="/"))) 
    dir.create(paste(out_path, "sppXsite", sep="/"))
  
  # get list of environmental layers
  env.list <- list.files("../worldclim_data/projections_2070/")
  env.list <- c("current", env.list)
  
  for(env in env.list){
    f <- niche.crops[grep(env, niche.crops, value=FALSE)]
    sppXsite <- tableFromRaster(f, threshold = 0.05)
    sppXsite <- data.frame(sppXsite)
    sppXsite$cell <- rownames(sppXsite)
    xy <- as.data.frame(raster(f[[1]]), xy=TRUE)
    xy$cell <- rownames(xy)
    sppXsite <- merge(sppXsite, xy)
    save(sppXsite, file=paste0(out_path, "/sppXsite/", env, ".rda"))
    rm(sppXsite)
  }
}

getTraitData <- function() {
  morph <- read.csv("InputData/MorphologyShort.csv", na.strings="9999")
  
  #just get males & the 3 traits of interest
  mon <- filter(morph, Sex == "Macho") %>%
    select(SpID, ExpC, Peso, AlCdo, TarsL) %>%
    group_by(SpID) %>%
    summarise_each(funs(mean(., na.rm = TRUE))) %>%
    filter(complete.cases(.))
  
  mon <- data.frame(mon)
  colnames(mon) <- c("Species","Bill","Mass","WingChord")
  rownames(mon) <- gsub(" ",".",mon$Species)
  mon <- mon[,-1]
  
  #principal component traits and get euclidean distance matrix
  means <- apply(mon, 2, mean)
  
  Bill <- (mon$Bill - means["Bill"])/sd(mon$Bill)
  Mass <- (mon$Mass - means["Mass"])/sd(mon$Mass)
  WingChord <- (mon$WingChord - means["WingChord"])/sd(mon$WingChord)
  
  z.scores <- data.frame(Bill, Mass, WingChord)
  rownames(z.scores) <- rownames(mon)
  return(z.scores)
}

ggbiplot2 <- function (pcobj, add.dat, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                       obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
                       ellipse = FALSE, ellipse.prob = 0.68, labels = NULL, labels.size = 3, 
                       alpha = 1, pt.size, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
                       varname.size = 3, varname.adjust = 1.5, varname.abbrev = FALSE, 
                       ...) 
{
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  
  nobs.factor <- sqrt(nrow(pcobj$x) - 1)
  d <- pcobj$sdev
  u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
  v <- pcobj$rotation
  
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, 
                              FUN = "*"))
  v <- sweep(v, 2, d^var.scale, FUN = "*")
  df.v <- as.data.frame(v[, choices])
  names(df.u) <- c("xvar", "yvar")
  names(df.v) <- names(df.u)
  
  df.u <- df.u * nobs.factor
  
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  if (obs.scale == 0) {
    u.axis.labs <- paste("standardized PC", choices, sep = "")
  } else {
    u.axis.labs <- paste("PC", choices, sep = "")
  }
  u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)", 
                                            100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  
  # get the additional data onto the PCA scores dataframe
  df.u <- cbind(df.u, add.dat)
  
  # create the dataframe with the variable axes
  df.v$varname <- rownames(v)
  df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
  
  # create ellipse
  theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
  circle <- cbind(cos(theta), sin(theta))
  ell <- ddply(df.u, .(scenario, groups), function(x) {
    if (nrow(x) <= 2) {
      return(NULL)
    }
    sigma <- var(cbind(x$xvar, x$yvar))
    mu <- c(weighted.mean(x$xvar, x$diff), weighted.mean(x$yvar, x$diff))
    ed <- sqrt(qchisq(ellipse.prob, df = 2))
    data.frame(sweep(circle %*% chol(sigma) * ed, 2, 
                     mu, FUN = "+"), groups = x$groups[1], scenario = x$scenario[1])
  })
  names(ell)[1:2] <- c("xvar", "yvar")
  
  
  # plot it all
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + 
    #scale_x_continuous(limits=c(0, 2.5)) + 
    #scale_y_continuous(limits=c(-2.5, 2)) + 
    xlab(u.axis.labs[1]) + 
    ylab(u.axis.labs[2]) +
    coord_equal()
  
  
  g <- g + geom_point(aes(size = diff, color = groups), alpha = alpha)
  
  g <- g + facet_wrap(~scenario, nrow=floor(length(unique(add.dat$scenario))/2))
  
  g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  
  g <- g + geom_segment(data = df.v, 
                        aes(x = 0, y = 0, xend = xvar, yend = yvar), 
                        arrow = arrow(length = unit(1/2, "picas")), color = muted("red"))
  
  g <- g + geom_text(data = df.v, 
                     aes(label = varname, x = xvar, y = yvar, angle = angle, hjust = hjust), 
                     color = "darkred", size = varname.size)
  
  
  
  return(g)
}

getSiteTraitValues <- function(f, traits) {
  fName <- str_match(f,pattern=paste("sppXsite/(\\w+.\\w+).rda",sep="/"))[,2]
  
  load(f)
  siteXspp <- t(sppXsite[,2:(ncol(sppXsite) - 3)])
  out <- matrix(nrow=ncol(siteXspp),ncol=ncol(traits))
  
  for (x in 1:ncol(siteXspp)){
    #subset site
    site <- siteXspp[,x]
    #get the trait matrix for species that are present
    site.trait<-traits[rownames(traits) %in% names(site[site==1]),]
    #mean position
    out[x,] <- colMeans(site.trait)
    
  }
  colnames(out) <- colnames(traits)
  out <- data.frame(x=sppXsite$x, y=sppXsite$y, mod=fName, out)
  return(out)
}

# function to calculate the hypervolume for each community
calcHV <- function(f, traits) {
  load(f)
  xy <- sppXsite[,c("x", "y")]
  sppXsite <- sppXsite[,2:(ncol(sppXsite)-3)]
  comm.hv <- sapply(rownames(sppXsite), function(k){
    g <- sppXsite[k,]
    sp <- names(g[which(g==1)])
    traits.sub <- subset(traits, rownames(traits) %in% sp)
    bw <- estimate_bandwidth(traits.sub, method="silverman")
    trait.hv <- ifelse(nrow(traits.sub) < 2, 0, get_volume(hypervolume(traits.sub, bandwidth = bw)))
  })
  comm.hv <- unlist(comm.hv)
  return(data.frame(xy, comm.hv=comm.hv))
}
