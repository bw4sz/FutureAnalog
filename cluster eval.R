# pbsim = phylobetasim dist object
s.all <- sum(pbsim) # total phylobeta

# firstly hierarchical cluster dendrograms are produced using all the methods in Kreft & Jetz 2010

library(doSNOW)
METHODS <- as.vector(c("ward", "single", "complete", "average", "mcquitty", "median", "centroid"))


cl <- makeCluster(4)
registerDoSNOW(cl)

meths <- foreach(i = METHODS)%dopar%(hclust(as.dist(pbsim), method = i))

stopCluster(cl)

library(cluster)
date()
m.DIANA <- diana(pbsim, diss = T)
date()

# nj.  N.B. nj produces an unrooted cluster dendrogram and therefore there is no obvious way to cut the dendrogram and make cluster
# The following code produces an nj tree and arbitrarily roots it using an imaginary grid cell that has 100% turnover from all the others

library(phangorn)
n.am <- cbind(as.matrix(pbsim),1)
n.am <- rbind(n.am,1)
rownames(n.am)[nrow(n.am)] <- "OG"
colnames(n.am)[ncol(n.am)] <- "OG"
besttree <- nj(n.am)
besttree <- root(besttree, "OG")
besttree <- drop.tip(besttree,"OG")
besttree$edge.length <- 0.1 + besttree$edge.length - min(besttree$edge.length)


m.list <- append(meths, list(as.hclust(m.DIANA), as.hclust(chronopl(besttree, lambda = 0, eval.max = 1, iter.max = 1))))
names(m.list)  <- c("ward", "single", "complete", "average", "mcquitty", "median", "centroid", "DIANA", "nj")

#  all of these clustering schemes are then evaluated for all possible numbers of clusters
# function codes in "Beval.R"


h.eval <- Beval.meth(m.list,as.matrix(pbsim),s.all, maxclust = attr(pbsim, "Size")-1)

# the PAM non-hierarchical method is also evaluated

nh.eval <- Beval.pam(pbsim,as.matrix(pbsim),s.all, maxclust = attr(pbsim, "Size")-1)

#results

eval <- rbind(h.eval,nh.eval)
rownames(eval) <- c(names(m.list),"PAM")
colnames(eval) <- 2:(attr(pbsim, "Size")-1)

