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
droppath <- "C:\\Users\\sarah\\Dropbox\\Hummingbirds\\NASA_Anusha\\"
gitpath <- "C:\\Users\\sarah\\Documents\\GitHub\\FutureAnalog\\"
rdata <- paste(output_folder, "\\AlphaMapping.RData", sep="")

#Load in data
load(rdata)

#Load in source functions
source(paste(gitpath,"AlphaMappingFunctions.R",sep=""))

setwd(paste(droppath,"FutureAnalog",sep=""))

#If testing the script grab a much smaller chunk
#current <- siteXspps[[1]][sample(1:nrow(siteXspps[[1]]),1000),]
#future <- siteXspps[[3]][sample(1:nrow(siteXspps[[3]]),1000),]

#If running the code with full dataset, for the full analysis
current <- siteXspps[[1]]
future <- siteXspps[[3]]

#Find within betadiversity
within.current.dist <- vegdist(current, "bray")
within.current <- as.matrix(within.current.dist)

#Find within phylobetadiversity
#For phylobeta, there needs to be more than 2 species for a rooted tree
phylo.current <- current[,colnames(current) %in% trx$tip.label]
phylo.current <- phylo.current[!apply(phylo.current,1,sum)<=2,]

phylo.future <- future[,colnames(future) %in% trx$tip.label]
phylo.future <- phylo.future[!apply(phylo.future,1,sum)<=2,]

#Find within Func betadiversity
Func.current <- current[,colnames(current) %in% colnames(fco)]
Func.current <- Func.current[!apply(Func.current,1,sum)<=2,]

Func.future <- future[,colnames(future) %in% colnames(fco)]
Func.future <- Func.future[!apply(Func.future,1,sum)<=2,]

#Within current phylobetadiversity
system.time(holt.try<-matpsim(phyl=trx,com=phylo.current,clust=3))

#Within current func betadiversity
#system.time(holt.func<-matpsim(phyl=tree.func,com=Func.current,clust=7))

####MNNTD method for integrating trait beta, needs to be checked, used in the DimDiv script
# MNNTD = Mean nearest neighbor taxon distance, from Holt et al. 2012. 
# An update of Wallace's zoogeographic regions of the world. Science.

source(paste(gitpath, "BenHolttraitDiversity.R", sep=""))

#create sp.list
sp.list<-lapply(rownames(Func.current),function(k){
  x<-Func.current[k,]
  names(x[which(x==1)])
})

names(sp.list)<-rownames(Func.current)

dists <- as.matrix(fco)

rownames(dists) <- rownames(fco)
colnames(dists) <- rownames(fco)

sgtraitMNTD <- sapply(rownames(Func.current),function(i){
  
  #Iterator count
  #print(round(which(rownames(siteXspp_traits)==i)/nrow(siteXspp_traits),3))
  
  #set iterator
  A<-i
  
  #
  out<-lapply(rownames(Func.current)[1:(which(rownames(Func.current) == i))], function(B) {MNND(A,B,sp.list=sp.list,dists=dists)})
  names(out)<-rownames(Func.current)[1:(which(rownames(Func.current) == i))]
  return(out)
})

names(sgtraitMNTD) <- rownames(Func.current)
melt.MNTD<-melt(sgtraitMNTD)

colnames(melt.MNTD)<-c("MNTD","To","From")

#turn beta measures into a matrices
within.current.phylo<-as.matrix(holt.try)

#needs to cast into a matrix to fit old formatting
#############needs to be done#######################
#turn into a matrix
within.current.func<-dcast(melt.MNTD,To~From,value.var="MNTD")
rownames(within.current.func)<-within.current.func[,1]
within.current.func<-within.current.func[,-1]

#within.current.func<-as.matrix(holt.func)

##################################
#Quantile Delineation Approach - sensu Strahlberg et al. 2009 - Not Currently Using
##################################
#Find the 5th quantile for each community
#quant.5<-apply(within.current,1,function(x){
  #quantile(x,.95)})
#names(quant.5)<-rownames(within.current)

#quant.phylo.5<-apply(within.current.phylo,1,function(x){
  #quantile(x,.95)})
#names(quant.phylo.5)<-rownames(within.current.phylo)

#quant.func.5<-apply(within.current.func,1,function(x){
  #quantile(x,.95)})
#names(quant.func.5)<-rownames(within.current.func)

#Once the clusters have been set, we can remove the large within current filezs
#rm(within.current,within.current.dist,holt.func,holt.try,within.current.func,within.current.phylo)
#gc()

###########################
#Between time taxonomic betadiversity
###########################

beta.time<-analogue::distance(current,future,"bray")

#####################
#CURRENT IS ROWS
#FUTURE IS COLUMNS
#####################

#For phylobetadiversity
#Between time phylobetadiversity
beta.time.phylo<-as.matrix(matpsim.pairwise(phyl=trx,com.x=phylo.current,com.y=phylo.future,clust=8))

#Repeat steps above for within time trait, but replacing Func.current with Func.future
#create sp.list
sp.list<-lapply(rownames(Func.future),function(k){
  x<-Func.future[k,]
  names(x[which(x==1)])
})

names(sp.list)<-rownames(Func.future)

dists <- as.matrix(fco)

rownames(dists) <- rownames(fco)
colnames(dists) <- rownames(fco)

sgtraitMNTD <- sapply(rownames(Func.future),function(i){
  
  #Iterator count
  #print(round(which(rownames(siteXspp_traits)==i)/nrow(siteXspp_traits),3))
  
  #set iterator
  A<-i
  
  #
  out<-lapply(rownames(Func.future)[1:(which(rownames(Func.future) == i))], function(B) {MNND(A,B,sp.list=sp.list,dists=dists)})
  names(out)<-rownames(Func.future)[1:(which(rownames(Func.future) == i))]
  return(out)
})

names(sgtraitMNTD) <- rownames(Func.future)
melt.MNTD<-melt(sgtraitMNTD)

colnames(melt.MNTD)<-c("MNTD","To","From")

#needs to be casted back into a matrix, see reshape2::dcast., name it betatime func
beta.time.func<-dcast(melt.MNTD,To~From,value.var="MNTD")
rownames(beta.time.func)<-beta.time.func[,1]
beta.time.func<-beta.time.func[,-1]

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

#For each of the current communities how many future communities are less difference 5th current quantile
n.analogs<-sapply(rownames(beta.time), function(x){
  sum(beta.time[rownames(beta.time) %in% x,] <= arb.thresh)
})

#Create a dataframe of cell cell numbers and number of analogs
current_to_future.analog<-data.frame(rownames(beta.time),n.analogs)
colnames(current_to_future.analog)<-c("cell.number","numberofanalogs")

#Visualize raster of analogs
fanalog<-cellVis(cell=current_to_future.analog$cell.number,value=current_to_future.analog$numberofanalogs)
hist(current_to_future.analog$numberofanalogs)
writeRaster(fanalog,"NumberofFutureAnalogs_Taxon_ARB.tif",overwrite=T)

#If a am counting the number of cells whose betadiveristy is within the 5th quantile, 
#the threshold for non analog communities should be the number of sites that oc
#If there are 5000 cells, then 5000*.05=250
threshold<-nrow(current)*.005

sum(current_to_future.analog$numberofanalogs <= threshold)
plot(noAnalogTax<-fanalog < threshold)

#For each of the currents communities how many future communities are less difference 5th current quantile
#n.analogs<-sapply(rownames(beta.time), function(x){
#sum(beta.time[rownames(beta.time) %in% x,] >= quant.5[names(quant.5) %in% x])
#})

###################
#Phylogenetic Analogs
####################
n.analogs.phylo<-sapply(rownames(beta.time.phylo), function(x){
  sum(beta.time.phylo[rownames(beta.time.phylo) %in% x,] <= arb.thresh)
})


#Create a dataframe of the number of analogs and the cellnumber
future.analog.phylo<-data.frame(rownames(phylo.current),n.analogs.phylo)
colnames(future.analog.phylo)<-c("cell.number","numberofanalogs")

#Visualize!
fanalog.phylo<-cellVis(cell=future.analog.phylo$cell.number,value=future.analog.phylo$numberofanalogs)

#Write to file
writeRaster(fanalog.phylo,"NumberofFutureAnalogs_Phylo_ARB.tif",overwrite=T)
hist(future.analog.phylo$numberofanalogs)
#If I am counting the number of cells whose betadiveristy is within the 5th quantile, 
#the threshold for non analog communities should be the number of sites that oc
#If there are 5000 cells, then 5000*.05=250

threshold<-nrow(phylo.current)*.05
sum(future.analog.phylo$numberofanalogs <= threshold)
plot(noAnalogPhylo<-fanalog.phylo < threshold)

###################
#Functional Analogs
####################
n.analogs.func<-sapply(rownames(beta.time.func), function(x){
  sum(beta.time.func[rownames(beta.time.func) %in% x,] <= arb.thresh)
})

#n.analogs.Func<-sapply(rownames(beta.time.phylo), function(x){
# sum(beta.time.func[rownames(beta.time.func) %in% x,] <= quant.func.5[names(quant.func.5) %in% x])
#})

future.analog.func<-data.frame(rownames(beta.time.func),n.analogs.func)
colnames(future.analog.func)<-c("cell.number","numberofanalogs")

fanalog.Func<-cellVis(cell=future.analog.func$cell.number,value=future.analog.func$numberofanalogs)

hist(future.analog.func$numberofanalogs)
#If a am counting the number of cells whose betadiveristy is within the 5th quantile, the threshold for non analog communities should be the number of sites that oc
#If there are 5000 cells, then 5000*.05=250

threshold<-nrow(Func.current)*.05
sum(future.analog.func$numberofanalogs <= threshold)

plot(fanalog.Func < threshold)
writeRaster(fanalog.Func,"NumberofFutureAnalogs_Func_ARB.tif",overwrite=T)

#####################################################
#PART II
#FUTURE ANALOGS IN CURRENT - NON-analog communities ("Novel")
#How many current communities have analogs in the future?
######################################################

###################
#Taxonomic Analogs
####################

#For each of the future communities how many future communities are less different 5th current quantile
n.analogs<-sapply(colnames(beta.time), function(x){
  sum(beta.time[,colnames(beta.time) %in% x] <= arb.thresh)
})

#Create a dataframe of cell cell numbers and number of analogs
future_to_current.analog<-data.frame(colnames(beta.time),n.analogs)
colnames(future_to_current.analog)<-c("cell.number","numberofanalogs")

#Visualize raster of analogs
fanalog<-cellVis(cell=future_to_current.analog$cell.number,value=future_to_current.analog$numberofanalogs)
#hist(future_to_current.analog$numberofanalogs)

#If a am counting the number of cells whose betadiveristy is within the 5th quantile, the threshold for non analog communities should be the number of sites that oc
#If there are 5000 cells, then 5000*.05=250
#threshold<-nrow(current)*.005


#sum(future_to_current.analog$numberofanalogs <= threshold)
#plot(noAnalogTax<-fanalog < threshold)

writeRaster(fanalog,"NumberofCurrentAnalogs_Taxo.tif",overwrite=T)

#For each of the currents communities how many future communities are less difference 5th current quantile
#n.analogs<-sapply(rownames(beta.time), function(x){
#sum(beta.time[rownames(beta.time) %in% x,] >= quant.5[names(quant.5) %in% x])
#})

###################
#Phylogenetic Analogs
####################

n.analogs.phylo<-sapply(colnames(beta.time.phylo), function(x){
  sum(beta.time.phylo[,colnames(beta.time.phylo) %in% x] <= arb.thresh)
})

future.analog.phylo<-data.frame(colnames(beta.time.phylo),n.analogs.phylo)
colnames(future.analog.phylo)<-c("cell.number","numberofanalogs")

fanalog.phylo<-cellVis(cell=future.analog.phylo$cell.number,value=future.analog.phylo$numberofanalogs)
hist(future.analog.phylo$numberofanalogs)
#If a am counting the number of cells whose betadiveristy is within the 5th quantile, the threshold for non analog communities should be the number of sites that oc
#If there are 5000 cells, then 5000*.05=250

writeRaster(fanalog.phylo,"NumberofCurrentAnalogs_Phylo.tif",overwrite=T)
threshold<-nrow(phylo.current)*.05
sum(future.analog.phylo$numberofanalogs <= threshold)
plot(noAnalogPhylo<-fanalog.phylo < threshold)

###################
#functional Analogs
####################
n.analogs.func<-sapply(colnames(beta.time.func), function(x){
  sum(beta.time.func[,colnames(beta.time.func) %in% x] <= arb.thresh,na.rm=TRUE)
})

future.analog.func<-data.frame(rownames(Func.future),n.analogs.func)
colnames(future.analog.func)<-c("cell.number","numberofanalogs")

fanalog.func<-cellVis(cell=future.analog.func$cell.number,value=future.analog.func$numberofanalogs)
hist(future.analog.func$numberofanalogs)
#If a am counting the number of cells whose betadiveristy is within the 5th quantile, the threshold for non analog communities should be the number of sites that oc
#If there are 5000 cells, then 5000*.05=250

threshold<-nrow(Func.current)*.05
sum(future.analog.func$numberofanalogs <= threshold)
plot(fanalog.func < threshold)
writeRaster(fanalog.func,"NumberofCurrentAnalogs_Func.tif",overwrite=T)


####################
#Compare outputs
####################
comp<-list.files(pattern=".tif",full.names=T)
all.raster<-stack(comp)
plot(all.raster)

########################
#Split into clusters, current and future
########################

#####################################
#naming of clusters is wrong here??  TODO: What is wrong? FIXME?
###################################
#plot(clusters<-all.raster[[c("TaxonomicClusters","PhylogenticClusters","FunctionalClusters")]],col=rainbow(5))
plot(current.ras<-all.raster[[c("NumberofCurrentAnalogs_Taxo","NumberofCurrentAnalogs_Phylo","NumberofCurrentAnalogs_Func")]])
plot(future.ras<-all.raster[[c("NumberofFutureAnalogs_Taxon_ARB","NumberofFutureAnalogs_Phylo_ARB","NumberofFutureAnalogs_Func_ARB")]])

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
