#This code goes through the results from AlphaMapping.R to determine the 
# number of analog hummingbird assemblages in Ecuador under future climate scenarios.
require(vegan)
require(picante)
require(reshape)
require(reshape2)
require(analogue)
require(doSNOW)
require(ape)
require(cluster)
require(RColorBrewer)
require(raster)
require(ggplot2)
require(phylobase)

#Sarah's
droppath <- "C:\\Users\\sarah\\Dropbox\\Hummingbirds\\NASA_Anusha\\" #Is this where we want the results to go?
gitpath <- "C:\\Users\\sarah\\Documents\\GitHub\\FutureAnalog\\"
output_folder <- "C:\\Users\\sarah\\Desktop\\Testmod"
rdata <- paste(output_folder, "\\AlphaMapping.RData", sep="")

#Load in data
load(rdata)

#Load in source functions
source(paste(gitpath,"\\AlphaMappingFunctions.R",sep=""))
source(paste(gitpath, "BenHolttraitDiversity.R", sep=""))

setwd(paste(droppath,"FutureAnalog",sep=""))


#If running the code with full dataset, for the full analysis
current <- siteXspps[[1]]
future <- siteXspps[2:4]

#Remove NAs from siteXspps so we can do the following analyses   
# Some species do not occur in Ecuador, so they should be removed from analysis here.
na.test <-  function (x) {
  w <- apply(x, 2, function(x)all(is.na(x)))
  if (any(w)) {
    fails <- names(which(w))
    print(paste("All NA in columns", paste(names(which(w)), collapse=", ")))
    return(fails)
  }
}

fails <- na.test(current[,])
current <- current[,!colnames(current) %in% fails]

future <- lapply(future, function(x){
  fails <- na.test(x)
  x[,!colnames(x) %in% fails]
})


#---------------- Find within SPECIES BETA DIVERSITY
within.current.dist <- vegdist(current, "bray")  
within.current <- as.matrix(within.current.dist)

within.future <- lapply(future, function(x){
  dist <- vegdist(x, "bray")
  as.matrix(dist)
})

#--------------- Find within PHYLO BETA DIVERSITY
#For phylobeta, there needs to be more than 2 species for a rooted tree
phylo.current <- current[,colnames(current) %in% trx$tip.label]
phylo.current <- phylo.current[!apply(phylo.current,1,sum)<=2,]   

phylo.future <- lapply(future, function(x){
  matched <- x[,colnames(x) %in% trx$tip.label]
  matched[!apply(matched,1,sum)<=2,] 
})


#Within current phylobetadiversity
system.time(holt.try <- matpsim(phyl=trx, com=phylo.current, clust=3))  
#turn beta measures into a matrix   
within.current.phylo <- as.matrix(holt.try)

within.future.phylo <- lapply(phylo.future, function(x){
  holt.try <- matpsim(phyl=trx, com=x, clust=3)
  as.matrix(holt.try)
})


#---------------- Find within FUNC BETA DIVERSITY

#----- Within current functional beta diversity
Func.current <- current[,colnames(current) %in% colnames(fco)]
Func.current <- Func.current[!apply(Func.current,1,sum)<=2,]  

####MNNTD method for integrating trait beta, used in the DimDiv script    TODO:  needs to be checked (@BenWeinstein)
#   MNNTD = Mean nearest neighbor taxon distance
#   Holt et al. 2012. An update of Wallace's zoogeographic regions of the world. Science.
sp.list<-lapply(rownames(Func.current),function(k){
  x<-Func.current[k,]
  names(x[which(x==1)])
})

names(sp.list) <- rownames(Func.current)
dists <- as.matrix(fco)
rownames(dists) <- rownames(fco)
colnames(dists) <- rownames(fco)

sgtraitMNTD <- sapply(rownames(Func.current),function(i){                    #TODO: why does this take so long!?
  A<-i   #set iterator
  out<-lapply(rownames(Func.current)[1:(which(rownames(Func.current) == i))], function(B) {
    MNND(A,B,sp.list=sp.list,dists=dists)
    })
  names(out) <- rownames(Func.current)[1:(which(rownames(Func.current) == i))]
  return(out)
})

names(sgtraitMNTD) <- rownames(Func.current)  #rownames are site ID numbers
melt.MNTD <- melt(sgtraitMNTD)
colnames(melt.MNTD) <- c("MNTD","To","From")

# Turn new results into a matrix
within.current.func <- cast(melt.MNTD,To ~ From, value="MNTD")
rownames(within.current.func) <- within.current.func[,1]
within.current.func <- within.current.func[,-1]

within.current.func[lower.tri(within.current.func)] <- t(within.current.func[upper.tri(within.current.func)])


#----- Within future functional beta          TODO: TEST that this works
Func.future <- lapply(future, function(x){
  matched <- x[,colnames(x) %in% colnames(fco)]
  matched[!apply(matched,1,sum)<=2,]   
})

sp.list <- lapply(Func.future,function(x){
  a <- lapply(rownames(x),function(k){
    x <- x[k,]
    names(x[which(x==1)])
  })
  names(a) <- rownames(x)
})

##Turn cophenetic distance to matrix
dists <- as.matrix(fco)
rownames(dists) <- rownames(fco)
colnames(dists) <- rownames(fco)

within.future.func <- lapply(Func.future,function(x){
  sgtraitMNTD <- sapply(rownames(x),function(i){
    A<-i    #set iterator
    out<-lapply(rownames(x)[1:(which(rownames(x) == i))], function(B) {
      MNND(A, B, sp.list=sp.list, dists=dists)
      })
    names(out) <- rownames(x)[1:(which(rownames(x) == i))]
    return(out)
  })
  
  names(sgtraitMNTD) <- rownames(x)
  melt.MNTD <- melt(sgtraitMNTD)
  colnames(melt.MNTD) <- c("MNTD","To","From")
  
  #needs to be casted back into a matrix, see reshape2::dcast., name it betatime func
  within.future.func <- cast(melt.MNTD,To ~ From, value="MNTD")
  rownames(within.future.func) <- beta.time.func[,1]
  within.future.func <- within.future.func[,-1]
  within.future.func[lower.tri(within.future.func)] <- t(within.future.func[upper.tri(within.future.func)])
  
  return(within.future.func)
})


###############################################################
#BETWEEN TIME
###############################################################

###########################
#Between time taxonomic betadiversity
###########################

beta.time<-lapply(future,function(x){
  analogue::distance(current,x,"bray")
})


###########################
#Between time phylogenetic
###########################
beta.time.phylo<-lapply(phylo.future,function(x){
  beta.time.phylo<-as.matrix(matpsim.pairwise(phyl=trx,com.x=phylo.current,com.y=x,clust=7))
})


###########################
#Between time functional
###########################

#Repeat steps above for within time trait, but replacing Func.current with Func.future
#
Beta.time.func<-lapply(Func.future,function(x){
  
  sp.list_current<-lapply(rownames(Func.current),function(k){
    g<-Func.current[k,]
    names(g[which(g==1)])
  })
  
  names(sp.list_current)<-rownames(Func.current)
  
  sp.list_future<-lapply(rownames(x),function(k){
    g<-x[k,]
    names(g[which(g==1)])
  })
  
  names(sp.list_future)<-rownames(x)
  
  #Get distances from the cophenetic matrix?
  dists <- as.matrix(fco)
  
  rownames(dists) <- rownames(fco)
  colnames(dists) <- rownames(fco)
  
  
  Beta.time.func<-lapply(rownames(x),function(fu){
    lapply(rownames(Func.current),function(cur){
      MNND_fc(fu,cur,sp.list_current,sp.list_future,dists)
    })
    
  })
  
  melt.MNTD<-melt(Beta.time.func)
  
  colnames(melt.MNTD)<-c("MNTD","To","From")
  
  #needs to be casted back into a matrix, see reshape2::dcast., name it betatime func
  beta.time.func<-dcast(melt.MNTD,To~From,value.var="MNTD")
  rownames(beta.time.func)<-beta.time.func[,1]
  beta.time.func<-beta.time.func[,-1]
  
  
  rownames(beta.time.func) <- rownames(Func.current)
  colnames(beta.time.func) <- rownames(x)
  return(beta.time.func)
})

#############################################
#             ANALOG ANALYSIS
#############################################
#TODO: Wrap this code in a function so we can test sensitivity of results to the threshold 
#      Set threshold for 5%, 10%, 50% and 100% - can present alternate results in appendices

#Set an arbitrary threshold      
arb.thresh<-.2

#################################
#PART I
#CURRENT ANALOGS IN FUTURE
#How many current communities have analogs in the future?
#These are akin to communities which will disappear, "Disappearing"
###################
#Taxonomic Analogs
###################

#For each of the current communities how many future communities are less different than the threshold
current_to_future.analog<-lapply(beta.time,function(j){
  n.analogs<-sapply(rownames(j), function(x){
    sum(j[rownames(j) %in% x,] <= arb.thresh)
  })
  current_to_future.analog<-data.frame(rownames(j),n.analogs)
  colnames(current_to_future.analog)<-c("cell.number","numberofanalogs")  
  return(current_to_future.analog)
})

c_f_tax<-lapply(current_to_future.analog,function(x){
  fanalog<-cellVis(cell=x$cell.number,value=x$numberofanalogs)
  hist(x$numberofanalogs)
  writeRaster(fanalog,"NumberofFutureAnalogs_Taxon_ARB.tif",overwrite=T)
})


###################
#Phylogenetic Analogs
####################

future.analog.phylo<-lapply(beta.time.phylo,function(j){
  n.analogs.phylo<-sapply(rownames(j), function(x){
    sum(j[rownames(j) %in% x,] <= arb.thresh)
  })
  
  #Create a dataframe of the number of analogs and the cellnumber
  future.analog.phylo<-data.frame(rownames(j),n.analogs.phylo)
  colnames(future.analog.phylo)<-c("cell.number","numberofanalogs")
  return(future.analog.phylo)
})


#Visualize!
c_f_phylo<-lapply(future.analog.phylo,function(x){
  fanalog.phylo<-cellVis(cell=x$cell.number,value=x$numberofanalogs)
  hist(x$numberofanalogs)
  writeRaster(fanalog.phylo,"NumberofFutureAnalogs_Phylo_ARB.tif",overwrite=T)
})

###################
#Functional Analogs
####################

future.analog.func<-lapply(Beta.time.func,function(j){
  n.analogs.func<-sapply(rownames(j), function(x){
    sum(j[rownames(j) %in% x,] <= arb.thresh)
  })
  
  future.analog.func<-data.frame(rownames(j),n.analogs.func)
  colnames(future.analog.func)<-c("cell.number","numberofanalogs")
  return(future.analog.func)
})

c_f_func<-lapply(future.analog.func,function(x){
  fanalog.Func<-cellVis(cell=x$cell.number,value=x$numberofanalogs)
  writeRaster(fanalog.Func,"NumberofFutureAnalogs_Func_ARB.tif",overwrite=T)
})


#############Visualize all three
c_f<-stack(c(c_f_tax,c_f_phylo,c_f_func))

plot(c_f)

###Across are the emissions scenerios, down are taxonomic, phylogenetic and functional for CURRENT TO FUTURE analogs (dissapearing)

#####################################################
#PART II
#FUTURE ANALOGS IN CURRENT - NON-analog communities ("Novel")
#How many current communities have analogs in the future?
######################################################

###################
#Taxonomic Analogs
####################

#For each of the future communities how many future communities are less different 5th current quantile
future_to_current.analog<-lapply(beta.time,function(j){
  n.analogs<-sapply(colnames(j), function(x){
    sum(j[,colnames(j) %in% x] <= arb.thresh)
  })
 future_to_current.analog<-data.frame(colnames(j),n.analogs)
  colnames(future_to_current.analog)<-c("cell.number","numberofanalogs")  
  return(future_to_current.analog)
})

f_c_tax<-lapply(future_to_current.analog,function(x){
  fanalog<-cellVis(cell=x$cell.number,value=x$numberofanalogs)
  hist(x$numberofanalogs)
  writeRaster(fanalog,"NumberofCurrentAnalogs_Taxon_ARB.tif",overwrite=T)
})

#Data check, these should be different
plot(stack(f_c_tax) - stack(c_f_tax))

###################
#Phylogenetic Analogs
####################

future_to_current.phylo<-lapply(beta.time.phylo,function(j){
  n.analogs.phylo<-sapply(colnames(j), function(x){
    sum(j[,colnames(j) %in% x] <= arb.thresh)
  })
  
  #Create a dataframe of the number of analogs and the cellnumber
  future.analog.phylo<-data.frame(colnames(j),n.analogs.phylo)
  colnames(future.analog.phylo)<-c("cell.number","numberofanalogs")
  return(future.analog.phylo)
})


#Visualize!
f_c_phylo<-lapply(future_to_current.phylo,function(x){
  fanalog.phylo<-cellVis(cell=x$cell.number,value=x$numberofanalogs)
  hist(x$numberofanalogs)
  writeRaster(fanalog.phylo,"NumberofCurrentAnalogs_Phylo_ARB.tif",overwrite=T)
})

plot(stack(f_c_phylo) - stack(c_f_phylo))


###################
#Functional Analogs
####################

future_to_current.func<-lapply(Beta.time.func,function(j){
  n.analogs.func<-sapply(colnames(j), function(x){
    sum(j[,colnames(j) %in% x] <= arb.thresh)
  })
  
  future.analog.func<-data.frame(colnames(j),n.analogs.func)
  colnames(future.analog.func)<-c("cell.number","numberofanalogs")
  return(future.analog.func)
})

f_c_func<-lapply(future_to_current.func,function(x){
  fanalog.Func<-cellVis(cell=x$cell.number,value=x$numberofanalogs)
  writeRaster(fanalog.Func,"NumberofCurrentAnalogs_Func_ARB.tif",overwrite=T)
})


#############Visualize all three
f_c<-stack(c(f_c_tax,f_c_phylo,f_c_func))

#plot novel assemblages
plot(f_c)


#plot both as a panel, this needs to be improved
#Just try plotting one emission scenerio across both disappearing and novel
novel<-f_c[[c(1,4,7)]]
disappear<-c_f[[c(1,4,7)]]

#This could be named correctly using 
firstplot<-stack(novel,disappear)
names(firstplot)<-c(paste("Novel",c("Tax","Phylo","Func")),paste("Disappearing",c("Tax","Phylo","Func")))
plot(firstplot)

#####################
#FIX ME: CLEANED UNTIL HERE 7/3/2014
#The rest should be fairly straightforward, correlating the rasters from above, the f_c raster is the number of future analogs of current assemblages
#The c_f is the number of current analogs of future assemblages


###########################
#Correlation among outputs, this currently onlymake sense for the beta diversity, the clusters are non-ordinal
###########################

#cluster.cor<-cor(values(clusters), use="complete.obs")
current.cor<-cor(values(current.ras), use="complete.obs")
future.cor<-cor(values(future.ras), use="complete.obs")

#make it a spare matrix
cluster.cor[upper.tri(cluster.cor)] <- NA
current.cor[upper.tri(current.cor)] <- NA
future.cor[upper.tri(future.cor)] <- NA

#Plot the correlations?

#m.current<-melt(current.cor)
#p<-ggplot(na.omit(m.current),aes(x=X1,y=X2,fill=as.numeric(value))) + geom_tile() + theme_bw()
#p<-p+coord_flip() + scale_x_reverse() + xlab("") + ylab("") + scale_x_discrete(labels=c("Taxonomic","Phylogenetic","Trait")) + scale_y_discrete(labels=c("Taxonomic","Phylogenetic","Trait"))
#p+scale_fill_continuous("Pearson Correlation",high="red",low="blue")

plot(current_arb<-all.raster[[c("NumberofFutureAnalogs_Taxon_ARB","NumberofFutureAnalogs_Phylo_ARB","NumberofFutureAnalogs_Func_ARB")]])

#standardize
current_standard<-current_arb/cellStats(current_arb,"max")

#corrlate with richness
ric<-raster(paste(gitpath, "Figures\\AlphaChange_Richness.tif", sep=""))
rc<-cor(values(current_arb[[1]]),values(ric),use="complete.obs")

#Pairwise plots of results, and glms
#Create giant dataframe
current_val<-data.frame(values(ric),values(current_arb))
colnames(current_val)<-c("Richness","Tax","Phylo","Trait")

ggplot(current_val,aes(Richness,Tax))+geom_point() + theme_bw()
ggplot(current_val,aes(Richness,Phylo))+geom_point() + theme_bw()
ggplot(current_val,aes(Richness,Trait))+geom_point() + theme_bw()

ggplot(current_val,aes(Phylo,Trait))+geom_point() + theme_bw()
ggplot(current_val,aes(Tax,Phylo))+geom_point() + theme_bw()
ggplot(current_val,aes(Tax,Trait))+geom_point() + theme_bw()


save.image("FutureAnalog.rData")

qplot(beta.time[1010,]) + theme_bw() + xlab("Community1010") + geom_vline(xintercept=quantile(beta.time[1010,],.2),col="red",linetype="dashed")

#######################
#Test number of analogs as a function of the threshold?  TODO: needs to be reviewed
########################

cl<-makeCluster(3,"SOCK")
registerDoSNOW(cl)

#Create a range of arb thresholds, output the number of analogs in each dimensions
novel.frame<-foreach(arb.thresh=seq(0,.5,.05)) %dopar% {
  require(reshape2)

  #Taxonomic
n.analogs<-sapply(colnames(beta.time), function(x){
  sum(beta.time[,colnames(beta.time) %in% x] <= arb.thresh)
})
future_to_current.analog<-data.frame(colnames(beta.time),n.analogs)
  colnames(future_to_current.analog)<-c("Cell","Taxonomic")

  #Phylogenetic
n.analogs.phylo<-sapply(colnames(beta.time.phylo), function(x){
  sum(beta.time.phylo[,colnames(beta.time.phylo) %in% x] <= arb.thresh)
})
  
future.analog.phylo<-data.frame(colnames(beta.time.phylo),n.analogs.phylo)
  colnames(future.analog.phylo)<-c("Cell","Phylogenetic")
  
  #Trait
n.analogs.func<-sapply(colnames(beta.time.func), function(x){
  sum(beta.time.func[,colnames(beta.time.func) %in% x] <= arb.thresh)
})
future.analog.func<-data.frame(rownames(Func.future),n.analogs.func)
  colnames(future.analog.func)<-c("Cell","Trait")
  
novel<-list(future_to_current.analog,future.analog.phylo,future.analog.func)
within.melt<-melt(novel)
  names(within.melt)[4]<-"del"
return(within.melt)}
stopCluster(cl)

names(novel.frame)<-seq(0,.5,.05)

#plot the number of analogs as a function of threshold, with each dimension as a series
m.novel<-melt(novel.frame)

m.novel<-m.novel[!m.novel$variable.1=="del",]

#pick 1000 random cells?
novel.samp<-m.novel[m.novel$Cell %in% sample(m.novel$Cell,1000),]

head(m.novel)
ggplot(m.novel,aes(L1,value)) + geom_smooth(aes(L1,value,group=1)) + facet_grid(.~variable) + theme_bw() + geom_point()
ggplot(novel.samp,aes(L1,value)) + geom_smooth(aes(L1,value,group=1)) + facet_grid(.~variable) + theme_bw() + geom_path(aes(group=Cell))

#Split out any arb value
novel.rasters<-lapply(1:length(novel.frame),function(x){
  fr<-novel.frame[[x]]
  fr.var<-split(fr,fr$variable)
  ras<-lapply(fr.var,function(y){
    cellVis(y$Cell,y$value)
  })
              names(ras)<-levels(fr$variable)
  return(stack(ras))
})
              
names(novel.rasters)<-names(novel.frame)

  s<-stack(novel.rasters[[2]],novel.rasters[[4]],novel.rasters[[6]],novel.rasters[[8]],novel.rasters[[10]])

names(s)<-sapply(c(2,4,6,8,10),function(x){
  paste(levels(fr$variable),names(novel.rasters)[[x]])
  })

plot(s,nc=3,zlim=c(0,1500))

save.image("FutureAnalog.rData")

##############
#######################
#Test number of disappearing analogs as a function of the threshold
########################

cl<-makeCluster(3,"SOCK")
registerDoSNOW(cl)

#Create a range of arb thresholds, output the number of analogs in each dimensions
disappear.frame<-foreach(arb.thresh=seq(0,.5,.05)) %dopar% {
  require(reshape2)
  
  #Taxonomic
  n.analogs<-sapply(rownames(beta.time), function(x){
    sum(beta.time[rownames(beta.time) %in% x,] <= arb.thresh)
  })
  future_to_current.analog<-data.frame(rownames(beta.time),n.analogs)
  colnames(future_to_current.analog)<-c("Cell","Taxonomic")
  
  #Phylogenetic
  n.analogs.phylo<-sapply(rownames(beta.time.phylo), function(x){
    sum(beta.time.phylo[rownames(beta.time.phylo) %in% x,] <= arb.thresh)
  })
  
  future.analog.phylo<-data.frame(rownames(beta.time.phylo),n.analogs.phylo)
  colnames(future.analog.phylo)<-c("Cell","Phylogenetic")
  
  #Trait
  n.analogs.func<-sapply(rownames(beta.time.func), function(x){
    sum(beta.time.func[rownames(beta.time.func) %in% x,] <= arb.thresh)
  })
  future.analog.func<-data.frame(rownames(beta.time.func),n.analogs.func)
  colnames(future.analog.func)<-c("Cell","Trait")
  
  novel<-list(future_to_current.analog,future.analog.phylo,future.analog.func)
  within.melt<-melt(novel)
  names(within.melt)[4]<-"del"
  return(within.melt)}
stopCluster(cl)

names(disappear.frame)<-seq(0,.5,.05)

fr<-disappear.frame[[1]]

#Split out any arb value
disappear.rasters<-lapply(1:length(disappear.frame),function(x){
  fr<-disappear.frame[[x]]
  fr.var<-split(fr,fr$variable)
  ras<-lapply(fr.var,function(y){
    cellVis(y$Cell,y$value)
  })
  names(ras)<-levels(fr$variable)
  return(stack(ras))
})

s.dis<-stack(disappear.rasters[[2]],disappear.rasters[[4]],disappear.rasters[[6]],disappear.rasters[[8]],disappear.rasters[[10]])

names(s.dis)<-sapply(c(2,4,6,8,10),function(x){
  paste(levels(fr$variable),names(disappear.rasters)[[x]])
})

plot(s.dis,nc=3,zlim=c(0,1500))

save.image("FutureAnalog.rData")
#s amnat
