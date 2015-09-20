# Function to calculate the functional betadiversity between each grid cell in
# two different matrices. Output is an RxC matrix where R is the number of sites
# in matrix 1 and C is the number of sites in matrix 2. The purpose here is to
# compare the current community composition to a projected future community
# composition for the end of calculating the number of analog cells.

# Based on functional.beta.pair from the betapart package by Andres Baselga,
# David Orme, Sebastien Villeger, Julien De Bortoli and Fabien Leprieur

functional.beta.c2f.pair <- function (cur, fu, traits, clust=20) 
{
  require(geometry)
  require(rcdd)
  CUR <- nrow(cur)
  FU <- nrow(fu)
  
  cur_FRi <- rep(NA, CUR)
  names(cur_FRi) <- row.names(cur)
  fu_FRi <- rep(NA, FU)
  names(fu_FRi) <- row.names(fu)
  
  coord_vert_i <- list()
  
  # this calculates the total functional richness based on 
  for (i in 1:CUR) {
    #current
    tr_i <- traits[which(cur[i, ] == 1), ]
    vert0 <- convhulln(tr_i, "Fx TO 'vert.txt'")
    vert1 <- scan("vert.txt", quiet = T)
    verti <- (vert1 + 1)[-1]
    coord_vert_i[[i]] <- tr_i[verti, ]
    cur_FRi[i] <- convhulln(tr_i[verti, ], "FA")$vol
  }
  
  for (i in 1:FU) {
    #future
    tr_i <- traits[which(fu[i, ] == 1), ]
    vert0 <- convhulln(tr_i, "Fx TO 'vert.txt'")
    vert1 <- scan("vert.txt", quiet = T)
    verti <- (vert1 + 1)[-1]
    coord_vert_i[[i]] <- tr_i[verti, ]
    fu_FRi[i] <- convhulln(tr_i[verti, ], "FA")$vol
  }
  
  #matrix in which to save intersection results
  vol_inter2_mat <- matrix(0, CUR, FU, dimnames = list(row.names(cur), 
                                                       row.names(fu)))
  
  #calculate the intersection between trait spaces
  unique.comm.pairs <- data.frame(set1=as.character(), set2=as.character(),
                                  stringsAsFactors = FALSE)
  interij <- as.numeric()
  
  vol_inter2_mat <- matrix(numeric(), nrow=CUR, ncol=FU, dimnames = c(rownames(CUR), rownames(FU)))
  
  for(i in 1:CUR) {
    for(j in 1:FU) {
      sp_i <- cur[i,]
      sp_i <- paste(names(sp_i[which(sp_i==1)]), collapse = ", ")
      sp_j <- fu[j,]
      sp_j <- paste(names(sp_j[which(sp_j==1)]), collapse = ", ")
      # this line finds a row which matches the current species community pair
      val <- which(apply(unique.comm.pairs, 1, function(x) all(x == c(sp_i, sp_j))))
      # if the above doesn't find a match, it searches for the pair reversed
      if(length(val)==0) val <- which(apply(unique.comm.pairs, 1, function(x) all(x == c(sp_j, sp_i)))) 
      
      if(length(val)==0) {
        unique.comm.pairs[nrow(unique.comm.pairs) + 1,] <- c(sp_i, sp_j)
        tr_i <- traits[which(cur[i, ] == 1), ]
        tr_j <- traits[which(fu[j, ] == 1), ]
        val <- length(interij) + 1
        interij[val] <- get_intersection(tr_i, tr_j)
        
      }
      vol_inter2_mat[i, j] <- interij[val]
    }
  }
  
  #use the above calculations to get the amount shared/not shared etc.
  shared <- vol_inter2_mat
  not.shared.cur <- apply(shared, 2, function(x) cur_FRi - x)
  not.shared.fu <- t(apply(shared, 1, function(x) fu_FRi - x))
  
  sum.not.shared <- not.shared.cur + not.shared.fu
  max.not.shared <- pmax(not.shared.cur, not.shared.fu)
  min.not.shared <- pmin(not.shared.cur, not.shared.fu)
  
  funct.beta.sim <- min.not.shared/(min.not.shared + shared)
  funct.beta.sne <- ((max.not.shared - min.not.shared)/((2 * shared) + sum.not.shared)) * (shared/(min.not.shared + shared))
  funct.beta.sor <- sum.not.shared/(2 * shared + sum.not.shared)
  functional.pairwise <- list(funct.beta.sim = funct.beta.sim, 
                              funct.beta.sne = funct.beta.sne, 
                              funct.beta.sor = funct.beta.sor)
  return(functional.pairwise)
}

#function to get the intersection of traitspace for two sites (required for above function)
get_intersection <- function(set1, set2) {
  set1rep <- d2q(cbind(0, cbind(1, set1)))
  set2rep <- d2q(cbind(0, cbind(1, set2)))
  polytope1 <- redundant(set1rep, representation = "V")$output
  polytope2 <- redundant(set2rep, representation = "V")$output
  H_chset1 <- scdd(polytope1, representation = "V")$output
  H_chset2 <- scdd(polytope2, representation = "V")$output
  H_inter <- rbind(H_chset1, H_chset2)
  V_inter <- scdd(H_inter, representation = "H")$output
  vert_1n2 <- q2d(V_inter[, -c(1, 2)])
  vol_inter <- 0
  if (is.matrix(vert_1n2) == T) 
    if (nrow(vert_1n2) > ncol(vert_1n2)) {
      vol_inter <- convhulln(vert_1n2, "FA")$vol
    }
  return(vol_inter)
}
