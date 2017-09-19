# ----------------
# Import libraries
# ----------------

library("NMF") #Used for NMF
library("vegan") #Used for rrarefy, and decostand
library("corrplot") #used to view dist matrix
library("fpc") #Used with plotcluster
library("clusteval") #Used for cluster_simularity function 

# ----------------
# Helper Functions
# ----------------

# go from x y # to dist matrix in no time flat 
nameDistToMatrix <- function(x){
  x.names <- sort(unique(c(x[[1]], x[[2]])))
  x.dist <- matrix(0, length(x.names), length(x.names))
  dimnames(x.dist) <- list(x.names, x.names)
  x.ind <- rbind(cbind(match(x[[1]], x.names), match(x[[2]], x.names)), cbind(match(x[[2]], x.names), match(x[[1]], x.names)))
  x.dist[x.ind] <- rep(x[[3]], 2)
  return(x.dist)
}

#creates a table of your favorite simulairties 
clusterSimularityTable <- function(clusteringsTable){
  simularityList = matrix(ncol=6)
  for( i in 1:ncol(clusteringsTable) ){
    for( n in i:ncol(clusteringsTable) ){
      colnamea <- colnames(clusteringsTable[i])
      colnameb <- colnames(clusteringsTable[n])
      jac <- cluster_similarity(clusteringsTable[[i]], clusteringsTable[[n]], similarity = c("jaccard"), method = "independence")
      rand <- cluster_similarity(clusteringsTable[[i]], clusteringsTable[[n]], similarity = c("rand"), method = "independence")
      wilcox <- wilcox.test(clusteringsTable[,i],clusteringsTable[,n])
      wilcox_p_value <- wilcox$p.value
      wilcox_statistic <- wilcox$statistic
      nameDistRow <- cbind(colnamea,colnameb,jac,rand,wilcox_p_value,wilcox_statistic )
      simularityList <- rbind(simularityList, nameDistRow)
    }
  }
  return(simularityList)
}

MantelTable <- function(mantelTable){
  mantelList = matrix(ncol=3)
  for( i in 1:length(mantelTable) ){
    for( n in i:length(mantelTable) ){
      colnamea <- names(mantelTable)[i]
      colnameb <- names(mantelTable)[n]
      if(!is.null(colnamea) && !is.null(colnameb)){
        mantel.info <- mantel(mantelTable[[i]], mantelTable[[n]])
        nameDistRow <- cbind(colnamea,colnameb,mantel.info$statistic)
        mantelList <- rbind(mantelList, nameDistRow)
      }
    }
  }
  return(mantelList)
}

# ----------------
# Prepare Some Color
# ----------------

warm1 <- c("#ff3333","#ff9933","#ffff33")

col1 <- colorRampPalette(c("#ff9933","white","#ff3333"))(200)
col2<- c("#cc33ff","#33ff33","#33ccff","#ff3333")
col3 <- colorRampPalette(c("#ff3333","white","#ff9933"))(200)
col4 <- colorRampPalette(c("#ff3333","white","#3399ff"))(200)

# ----------------
# Prepare Otu Data
# ----------------

otus <- read.csv("OTU-simple.txt", as.is=TRUE)

#prepare uneffected otu table
otus.plain <- otus[,3:61]

#prepare log transform data
otus.log <- decostand(otus.plain, method="log")

#preare normalized data
otus.norm <- decostand(otus.plain, method="normalize")

#prepare rareified data
otus.rarefy <- rrarefy(otus.plain, rowSums(otus.plain))

# ----------------
# Looking at heatmaps
# ----------------

otu.mantel.box <- list()

#png(filename = "distance_corrs.png", width = 1024, height = 1024)
#par( mfrow = c( 2, 2 ) )

otus.plain.dist <- as.dist(vegdist(otus.plain))
corrplot(as.matrix(otus.plain.dist), 
         main="Plain OTU",
         cex.main=2,
         mar=c(0,0,1,0),
         order="original", 
         method="color", 
         col=col4, 
         cl.lim=c(0,1), 
         is.corr=FALSE, 
         tl.col="white")
         #addgrid.col="white")

otu.mantel.box[[1]] <- otus.plain.dist

otus.log.clean.dist <- as.dist(vegdist(otus.log))
corrplot(as.matrix(otus.log.clean.dist), 
         main="Log OTU",
         cex.main=2,
         mar=c(0,0,1,0),
         order="original", 
         method="color", 
         col=col4, 
         cl.lim=c(0,1), 
         is.corr=FALSE, 
         tl.col="white")
#addgrid.col="white")

otu.mantel.box[[2]] <- otus.log.clean.dist

otus.norm.dist <- as.dist(vegdist(otus.norm))
corrplot(as.matrix(otus.norm.dist), 
         main="Normalized OTU",
         cex.main=2,
         mar=c(0,0,1,0),
         order="original", 
         method="color", 
         col=col4, 
         cl.lim=c(0,1), 
         is.corr=FALSE, 
         tl.col="white")
#addgrid.col="white")

otu.mantel.box[[3]] <- otus.norm.dist

otus.rarefy.dist <- as.dist(vegdist(otus.rarefy))
corrplot(as.matrix(otus.rarefy.dist), 
         main="Rarefied OTU",
         cex.main=2,
         mar=c(0,0,1,0),
         order="original", 
         method="color", 
         col=col4, 
         cl.lim=c(0,1), 
         is.corr=FALSE, 
         tl.col="white")
#addgrid.col="white")

#dev.off()

otu.mantel.box[[4]] <- otus.rarefy.dist

names(otu.mantel.box) <- c("otus.plain.dist","otus.log.dist","otus.norm.dist","otus.rarefy.dist")
otu.mantel.table <- MantelTable(otu.mantel.box)

otu.mantel.table <- otu.mantel.table[2:11,]
write.csv(otu.mantel.table, file="otu.mantel.table.csv")
otu.mantel.table.csv <- read.csv("otu.mantel.table.csv", as.is=TRUE)
otu.mantel.table.dist <- nameDistToMatrix(otu.mantel.table.csv[,2:4])

#png(filename = "mantel_test.png", width = 340, height = 340)
corrplot(otu.mantel.table.dist, 
         order="hclust", 
         method="color", 
         col=col3,
         cl.lim=c(0,1), 
         mar=c(0,0,1,0),
         is.corr=FALSE, 
         hclust.method="centroid", 
         tl.col="black", 
         addgrid.col="white", 
         addCoef.col="black", 
         type="upper")
#dev.off()

# ----------------
# Looking at PCA
# ----------------

#png(filename = "PCA_andSil.png", width = 512, height = 1024)
#par( mfrow = c( 4, 2 ) )
par( mfrow = c( 4, 2 ) )
otu.pca.box <- data.frame(1:69)

otus.plain.kmeans <- kmeans(otus.plain,3)
otus.plain.trans = preProcess(otus.plain,
                       method=c("BoxCox", "center", 
                                "scale", "pca"))
otus.plain.PC = predict(otus.plain.trans, otus.plain)
plotcluster(otus.plain.PC, otus.plain.kmeans$cluster, col=col2[otus.plain.kmeans$cluster], pch=16, method="dc")

otus.plain.dissE <- daisy(otus.plain.PC)
otus.plain.dE2   <- otus.plain.dissE^2
otus.plain.sk2   <- silhouette(otus.plain.kmeans$cluster, otus.plain.dE2)
plot(otus.plain.sk2, col=warm1, main="Kmeans PCA Plain OTUs", cex.axis=2, cex.names=2, cex.lab=0.01, do.n.k = FALSE, do.clus.stat = TRUE)

otu.pca.box$otus.plain <- otus.plain.kmeans$cluster


otus.log.kmeans <- kmeans(otus.log,3)
otus.log.trans = preProcess(otus.log,
                              method=c("BoxCox", "center", 
                                       "scale", "pca"))
otus.log.PC = predict(otus.log.trans, otus.log)
plotcluster(otus.log.PC, otus.log.kmeans$cluster, col=col2[otus.log.kmeans$cluster], pch=16, method="dc")

otus.log.dissE <- daisy(otus.log.PC)
otus.log.dE2   <- otus.log.dissE^2
otus.log.sk2   <- silhouette(otus.log.kmeans$cluster, otus.log.dE2)
plot(otus.log.sk2, col=warm1, main="Kmeans PCA Log OTUs", cex.axis=2, cex.names=2, cex.lab=0.01, do.n.k = FALSE, do.clus.stat = TRUE)

otu.pca.box$otus.log <- otus.log.kmeans$cluster


otus.norm.kmeans <- kmeans(otus.norm,3)
otus.norm.trans = preProcess(otus.norm,
                            method=c("BoxCox", "center", 
                                     "scale", "pca"))
otus.norm.PC = predict(otus.norm.trans, otus.norm)
plotcluster(otus.norm.PC, otus.norm.kmeans$cluster, col=col2[otus.norm.kmeans$cluster], pch=16, method="dc")

otus.norm.dissE <- daisy(otus.norm.PC)
otus.norm.dE2   <- otus.norm.dissE^2
otus.norm.sk2   <- silhouette(otus.norm.kmeans$cluster, otus.norm.dE2)
plot(otus.norm.sk2, col=warm1, main="KKmeans PCA Norm OTUs", cex.axis=2, cex.names=2, cex.lab=0.01, do.n.k = FALSE, do.clus.stat = TRUE)

otu.pca.box$otus.norm <- otus.norm.kmeans$cluster


otus.rarefy.kmeans <- kmeans(otus.rarefy,3)
otus.rarefy.trans = preProcess(otus.rarefy,
                             method=c("BoxCox", "center", 
                                      "scale", "pca"))
otus.rarefy.PC = predict(otus.rarefy.trans, otus.rarefy)
plotcluster(otus.rarefy.PC, otus.rarefy.kmeans$cluster, col=col2[otus.rarefy.kmeans$cluster], pch=16, method="dc")

otus.rarefy.dissE <- daisy(otus.rarefy.PC)
otus.rarefy.dE2   <- otus.rarefy.dissE^2
otus.rarefy.sk2   <- silhouette(otus.rarefy.kmeans$cluster, otus.rarefy.dE2)
plot(otus.rarefy.sk2, col=warm1, main="Kmeans PCA Rarefy OTUs", cex.axis=2, cex.names=2, cex.lab=0.01, do.n.k = FALSE, do.clus.stat = TRUE)

#dev.off()

otu.pca.box$otus.rarefy <- otus.rarefy.kmeans$cluster

otu.pca.box <- otu.pca.box[,2:5]

otu.pca.box.sim <- clusterSimularityTable(otu.pca.box)
otu.pca.box.sim <- otu.pca.box.sim[2:11,]
write.csv(otu.pca.box.sim, file="otu.pca.box.sim.csv")
otu.pca.box.csv <- read.csv("otu.pca.box.sim.csv", as.is=TRUE)
otu.pca.box.dist <- nameDistToMatrix(otu.pca.box.csv[,2:4])

#png(filename = "Jaccard_PCA_Corr.png", width = 340, height = 340)
corrplot(otu.pca.box.dist, 
         order="hclust",
         method="color", 
         col=col1, 
         cl.lim=c(0,1), 
         mar=c(0,0,1,0),
         is.corr=FALSE, 
         hclust.method="centroid", 
         tl.col="black", 
         addgrid.col="white", 
         addCoef.col="black", 
         type="upper")
#dev.off()

# ----------------
# Looking at NMF
# ----------------

png(filename = "NMF_coefmaps.png", width = 1024, height = 2048)
par( mfrow = c( 4, 2 ) )

otu.nmf.box <- data.frame(1:59)

otus.plain.nmf <- nmf(otus.plain,3)
otus.plain.nmf.sil <- silhouette(otus.plain.nmf)
plot(otus.plain.nmf.sil, col=warm1, main="NMF On Plain OTUs, K=3", cex.axis=2, cex.names=2, cex.lab=0.01,  do.n.k = FALSE, do.clus.stat = TRUE)
coefmap(otus.plain.nmf)

otu.nmf.box$otus.plain.nmf <- predict(otus.log.nmf)


otus.log.nmf <- nmf(otus.log,3)
otus.log.nmf.sil <- silhouette(otus.log.nmf)
plot(otus.log.nmf.sil, col=warm1, main="NMF On Log OTUs, K=3", cex.axis=2, cex.names=2, cex.lab=0.01,  do.n.k = FALSE, do.clus.stat = TRUE)
coefmap(otus.log.nmf)

otu.nmf.box$otus.log.nmf <- predict(otus.log.nmf)


otus.norm.nmf <- nmf(otus.norm,3)
otus.norm.nmf.sil <- silhouette(otus.norm.nmf)
plot(otus.norm.nmf.sil, col=warm1, main="NMF On Normalized OTUs, K=3", cex.axis=2, cex.names=2, cex.lab=0.01,  do.n.k = FALSE, do.clus.stat = TRUE)
coefmap(otus.norm.nmf)

otu.nmf.box$otus.norm.nmf <- predict(otus.norm.nmf)


otus.rarefy.nmf <- nmf(otus.rarefy,3)
otus.rarefy.nmf.sil <- silhouette(otus.rarefy.nmf)
plot(otus.rarefy.nmf.sil, col=warm1, main="NMF On Rarefied OTUs, K=3", cex.axis=2, cex.names=2, cex.lab=0.01,  do.n.k = FALSE, do.clus.stat = TRUE)
coefmap(otus.rarefy.nmf)

otu.nmf.box$otus.rarefy.nmf <- predict(otus.rarefy.nmf)

dev.off()

compare(list(otus.plain.nmf, otus.log.nmf,otus.norm.nmf,otus.rarefy.nmf))

otu.nmf.box <- otu.nmf.box[,2:5]
write.csv(otu.nmf.box, file="otu.nmf.box.csv")
otu.nmf.box.csv <- read.csv("otu.nmf.box.csv", as.is=TRUE)
otu.nmf.box.csv <- otu.nmf.box.csv[,2:5]

otu.nmf.box.sim <- clusterSimularityTable(otu.nmf.box.csv)
otu.nmf.box.sim <- otu.nmf.box.sim[2:11,]
write.csv(otu.nmf.box.sim, file="otu.nmf.box.sim.csv")
otu.nmf.box.csv <- read.csv("otu.nmf.box.sim.csv", as.is=TRUE)
otu.nmf.box.dist <- nameDistToMatrix(otu.nmf.box.csv[,2:4])

#png(filename = "NMF_Jaccard_simu.png", width = 340, height = 340)
corrplot(otu.nmf.box.dist,
         order="hclust", 
         method="color", 
         col=col1,
         cl.lim=c(0,1), 
         mar=c(0,0,1,0),
         is.corr=FALSE, 
         hclust.method="centroid", 
         tl.col="black", 
         addgrid.col="white", 
         addCoef.col="black", 
         type="upper")
#dev.off()

