
nCores <- 4
goIndividual <- 1
goCluster <- 0
goClassify <- 0
bed1Name <- c("par-100kbW-50kbS")
bed2Name <- c("par-100kbW_50kbS-randomized")

options(max.print=1000000)

library(ggbiplot)
library(ggplot2)
library(kSamples)
library(Matching)
library(parallel)
library(gplots)
library(ape)
library(data.table)

library(magrittr)
library(vegan)
library(pvclust)
library(mclust)
library(fpc)
library(class)
library(Hmisc)
library(ggdendro)
library(philentropy)
library(RColorBrewer)

library(tidyverse)
library(cluster)
library(factoextra)
library(Rtsne)
library(dendextend)
library(parallelDist)
library(abind)
#library(rgl)

b1File <- c("bed1Stats.tsv")
b2File <- c("bed2Stats.tsv")

outDir <- "Clustering/"
dir.create("Clustering")
dir.create("Clustering/pam")
dir.create("Clustering/intervalRepeatSimilarities")
dir.create("Clustering/mantel")

distType <- 0
hclustMethodInUse <- "ward.D2"
mdsDimensions <- 3
i<-2
bedNameConverts <- c(bed1Name, bed2Name)

###read in data
b1FullData <- as.data.frame(fread(b1File, header=FALSE, sep = "\t"))
b2FullData <- as.data.frame(fread(b2File, header=FALSE, sep = "\t"))

catLabels <- c("species", "class", "family")
repeatFileList <- c("species-repeatList.tsv", "class-repeatList.tsv", "family-repeatList.tsv")


regionList <- as.data.frame(fread("intervalLabels-c.tsv", header=FALSE, sep = "\t"))
regionList <- sapply(regionList, as.character)

rmRG <- grepl ("--bed1", regionList)
bed1Reg <- regionList[rmRG]
bed2Reg <- regionList[!rmRG]

pamSpam <- function(kPam, distMat){
  pam_fit <- pam(distMat, diss = TRUE, k = kPam)
  print(c("on k: ", kPam))
  result <- c(pam_fit$silinfo$avg.width, kPam)

  conOut <- paste(outDir, "/pam/", paste("pam_clustering-info-at-k", kPam, repCat, sep = "-"), ".txt", sep = "")
  sink(file = conOut)
  print("medoids")
  pam_fit$medoids
  print("cluster info")
  pam_fit$clusinfo
  sink()
  
  bed <- matrix(unlist(strsplit(gsub("[:-]", " ", names(pam_fit$clustering)), " +")), ncol = 4, byrow=TRUE)  
  bed <- unname(cbind(bed, as.matrix(pam_fit$clustering)))
  bedOut <- paste(outDir, "/pam/", paste("pam_clustering-at-k", kPam, repCat, sep = "-"), ".bed", sep = "")
  write.table(bed, file=bedOut, col.names=F, row.names=F, sep="\t", quote=FALSE)

  pdfFN <- paste(outDir, "/pam/", paste("pam_sillouhettePlot-k", kPam, repCat, sep = "-"), ".pdf", sep = "")
  pdf(pdfFN)
  plot(pam_fit)
  dev.off()
  
  return(pam_fit$silinfo$avg.width)
}

createBedPE <- function(preBedpe, repName){
	peOut <- c("", "", "", "", "", "", "", "")
	for(ipbpe in 1:dim(preBedpe)[1]){
		pe <- preBedpe[ipbpe,]
		if(!(pe[1] == pe[5] && pe[2] == pe[6] && pe[3] == pe[7])){
			peLabel <- paste("diffBed", pe[1], pe[5], sep="-")
			if(identical(pe[4], pe[8])){
				peLabel <- paste("sameBed", pe[1], pe[5], sep="-")
			}
			peLabel <- paste(peLabel, repName, sep="")
			outLine <- c(pe[1], pe[2], pe[3], pe[5], pe[6], pe[7], peLabel, pe[9])
			peOut <- rbind(peOut, outLine)
		}
	}
	peOut <- peOut[-1,]
	return(peOut)
}

bedSimWrites <- function(repNum, repIntComMat1, repIntComMat2){ 
	crMat <- rbind(repIntComMat1[,repNum,], repIntComMat2[,repNum,])
	#remove the extra column that abind adds
	crMat <- crMat[,1:length(statsTypes)]
	colnames(crMat) <- statsTypes
	crMat[is.na(crMat)] <- 0

	crdist <- getDist(crMat, distType)

	crBedPE <- melt(crdist)[complete.cases(melt(crdist)),]

	preBedpe <- paste(crBedPE[,1], crBedPE[,2], crBedPE[,3], sep=" ")
	preBedpe <- matrix(unlist(strsplit(gsub("[:-]", " ", preBedpe), " +")), ncol = 9, byrow=TRUE)

	bed <- createBedPE(preBedpe, reps[repNum])
	colnames(bed) <- c("#seg1", "start1", "stop1", "seg2", "start2", "stop2", "label", "score") #useful for circa
    bedFN <- paste(outDir, "/intervalRepeatSimilarities/", paste("bedScoresFor-", repCat, reps[repNum], sep = "-"), ".bed", sep = "")
    write.table(bed, file=bedFN, col.names=T, row.names=F, sep="\t", quote=FALSE)
}

getDist <- function(matData, choice){
	if(choice == 0){
		matDistGot <- distance(matData, method="gower", )
	}else if(choice == 1){
		matDistGot <- parDist(as.matrix(matData), method="bray", threads=nCores)
		matDistGot <- as.matrix(matDistGot)
	}
	colnames(matDistGot) <- rownames(matData)
	rownames(matDistGot) <- rownames(matData)
	matDistGot <- as.dist(matDistGot)
	
	return(matDistGot)
}

pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni'){
#modified from https://www.researchgate.net/post/How_can_I_do_PerMANOVA_pairwise_contrasts_in_R
	co = combn(unique(as.character(factors)),2)
	pairs = c()
	F.Model =c()
	R2 = c()
	p.value = c()

	for(elem in 1:ncol(co)){
		print(paste("on: ", elem, " ", paste(co[1,elem],'vs',co[2,elem]), sep=""))
		x1 = vegdist(data.matrix(x[factors %in% c(co[1,elem],co[2,elem]),]),method="gower")
		ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])], permutations=9999, parallel = nCores );
		pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
		F.Model =c(F.Model,ad$aov.tab[1,4]);
		R2 = c(R2,ad$aov.tab[1,5]);
		p.value = c(p.value,ad$aov.tab[1,6])
	}
	p.adjusted = p.adjust(p.value,method="bonferroni")
	sig = c(rep('',length(p.adjusted)))
	sig[p.adjusted <= 0.05] <-'.'
	sig[p.adjusted <= 0.01] <-'*'
	sig[p.adjusted <= 0.001] <-'**'
	sig[p.adjusted <= 0.0001] <-'***'

	pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
	print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
	return(pairw.res)

} 


###create plots for individual types of repeats
###remove empty elements

statsTypes <- c("avgDivRepSeq", "avgDelRepSeq", "avgInsRepSeq", "avgIdtRepSeq", "avgLenRepSeq", "sumRepLenSeq", "sumRepCountSeq", 
"sumRepLenNormRepTotalLen", "sumRepLenNormDNATotalLen", "sumRepLenNormRepIntervalLen", "sumRepLenNormDNAIntervalLen", 
"sumRepCountNormRepTotalLen", "sumRepCountNormDNATotalLen", "sumRepCountNormRepIntervalLen", "sumRepCountNormDNAIntervalLen", 
"sumRepCountNormRepTotalCount", "sumRepCountNormRepIntervalCount")

humanStats <- c("average divergence", "average deletion", "average insertion", "average identity", "average length", "sum of repeat length", "sum of rep count", 
"sum of repeat length / total length of repeats", "sum of repeat length / total length of DNA", "repeat length as fraction of repeats in interval", "repeat length as fraction of interval",
"number of repeats normalized to length of repeats", "number of repeats normalized to total length of DNA", "number of repeats normalized to length of repeats on interval", "number of repeats normalized to length of interval", 
"fraction of number of repeats", "fraction of number of repeats in interval")

for(i in 3:1){
repCat <- catLabels[i];
reps <- t(as.data.frame(fread(repeatFileList[i], header=FALSE, sep = "\t")))
repLen <- length(reps)
repsName <- reps[1:repLen]
repType <- catLabels[i]


mat1Intervals <- data.frame(matrix(NA, nrow = length(bed1Reg), ncol=0))
mat2Intervals <- data.frame(matrix(NA, nrow = length(bed2Reg), ncol=0))

for(sm in 1:length(statsTypes)){
	###load in and add labels to data
	b1FN <- paste(paste(statsTypes[sm], repCat, "b1", sep = "-"), ".tsv", sep = "")
	b1Mat <- as.data.frame(fread(b1FN, header=FALSE, sep = "\t"))
	names(b1Mat) <- reps
	row.names(b1Mat) <- bed1Reg
	colnames(b1Mat) <- paste(statsTypes[sm], colnames(b1Mat), sep = "_-_")
	if(sm == 1 || sm == 2 || sm == 3){
		b1Mat <- 100 - b1Mat
	}
	b1Mat[is.na(b1Mat)] <- 0
	mat1Intervals <- cbind(mat1Intervals, b1Mat)
	
	b2FN <- paste(paste(statsTypes[sm], repCat, "b2", sep = "-"), ".tsv", sep = "")
	b2Mat <- as.data.frame(fread(b2FN, header=FALSE, sep = "\t"))
	names(b2Mat) <- reps
	row.names(b2Mat) <- bed2Reg
	colnames(b2Mat) <- paste(statsTypes[sm], colnames(b2Mat), sep = "_-_")
	if(sm == 1 || sm == 2 || sm == 3){
		b2Mat <- 100 - b2Mat
	}
	b2Mat[is.na(b2Mat)] <- 0
	mat2Intervals <- cbind(mat2Intervals, b2Mat)
	
}
matCIntervals <- rbind(mat1Intervals, mat2Intervals)
comMatLabels <- c(rep("Bed1", length(bed1Reg)), rep("Bed2", length(bed2Reg)))
comMatColors <- c(rep("blue", length(bed1Reg)), rep("red", length(bed2Reg)))

mciDist <- getDist(matCIntervals, distType)

mciHclust <- hclust(mciDist, hclustMethodInUse)
mciDen <- as.dendrogram(mciHclust)

colorCodes <- c(Bed1 = "blue", Bed2 = "red")
leafColors <- colorCodes[comMatLabels][order.dendrogram(mciDen)]
writeLabels <- F
if(length(leafColors) < 200000){
	writeLabels <- T
}

mciDen <- mciDen %>% set("leaves_pch", 19) 
mciDen <- assign_values_to_leaves_nodePar(mciDen, leafColors, "col")
mciDen <- mciDen %>% set("labels_cex", 0.5)
        
pdfFN <- paste(outDir, "/", paste("dendrogram_allIntervals", repCat, sep = "-"), ".pdf", sep = "")
#cairo_pdf(file = pdfFN)
pdf(pdfFN, height = (length(leafColors)/8), width = (length(leafColors)/6))

mciGGplot <- ggplot(mciDen, horiz=T, labels=writeLabels) + expand_limits(y=c(-0.5,-0.5), x=c(-2,1))
mciGGplot <- mciGGplot + labs(title="Gower similarity of all regions\nusing ward.D2 to cluster") #+ ylab("Regions") + xlab("Distance")
mciGGplot <- mciGGplot + theme(plot.title = element_text(hjust = 0.5), axis.title.x=element_text(vjust=-2))
plot(mciGGplot)
#legend(x="bottomleft",legend= paste("blue =", bed1Name, " regions\nred =", bed2Name , " regions", sep = " "), fill="grey94")

dev.off()

intCMDS <- metaMDS(mciDist, autotransform=FALSE, k=mdsDimensions, trymax=10000, parallel=as.integer(nCores)) #plot=TRUE
     
pdfFN <- paste(outDir, "/", paste("metaMDS_stressplot_allIntervals", sep = "-"), ".pdf", sep = "")
#cairo_pdf(file = pdfFN)
pdf(pdfFN)
stressplot(intCMDS)
dev.off()

pdfFN <- paste(outDir, "/", paste("metaMDS_ordiplot_allIntervals", repCat, sep = "-"), ".pdf", sep = "")
pdf(pdfFN)
ordiplot(intCMDS)
ordihull(intCMDS, groups=comMatLabels, draw='lines', col="darkgrey", label=F)
for(icml in unique(comMatLabels)) {
  ordiellipse(intCMDS$point[grep(icml,comMatLabels),],draw="polygon",
  groups=comMatLabels[comMatLabels==icml],col=comMatColors[grep(icml,comMatLabels)],label=F, conf=0.95, alpha=30) 
} 
orditorp(intCMDS,display="sites",col=c(rep("blue", length(bed1Reg)), rep("red", length(bed2Reg))), labels=F, pch=1)
z <- kde2d(intCMDS$points[,1], intCMDS$points[,2], n=50)
k <- 10
my.cols <- rev(brewer.pal(k, "Spectral"))
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE, lwd=1)
dev.off()

cMat <- matCIntervals[,apply(matCIntervals, 2, var, na.rm=TRUE) != 0]
comMatLabels <- c(rep("Bed1", length(bed1Reg)), rep("Bed2", length(bed2Reg)))

if(0){
	comMatPCA <- prcomp(cMat, scale=T)
	pdfFN <- paste(outDir, "/", paste("ggbiplot_allIntervals", repCat, sep = "-"), ".pdf", sep = "")
	pdf(file = pdfFN)
	ggbp <- ggbiplot(comMatPCA, obs.scale = 1, var.scale = 1, var.axes=FALSE,
			groups = comMatLabels, ellipse = TRUE, circle = TRUE) +
			scale_color_discrete(name = "") +
			theme(legend.direction = "horizontal", legend.position = "top")
	print(ggbp)
	dev.off()
}

rglCMDS <- c("")
if(mdsDimensions > 2 && 0){
	library(rgl)
	NMDS.1 <- intCMDS$points[,1]
	NMDS.2 <- intCMDS$points[,2]
	NMDS.3 <- intCMDS$points[,3]

	#rgl.open(useNULL=TRUE)
	rgl.open(useNULL = rgl.useNULL())
	plot3d(NMDS.1, NMDS.2, NMDS.3, col=comMatColors)
	movie3d(spin3d(axis=c(0,1,1), rpm=5), dev=rgl.cur(), dir=paste(getwd(), outDir, sep="/"), duration=24, movie="rgl-combined-3D-MDS")
	#rgl.quit()
}else{
	rglCMDS <- cbind(rglCMDS, intCMDS) #save it just in case...
}

#plot(fviz_dist(mciUsed, gradient = list(low = "blue", mid = "white", high = "red", ordered=T)) + theme(text = element_text(size = 8)))
pdfFN <- paste(outDir, "/", paste("heatmap_regions_allIntervals", repCat, sep = "-"), ".pdf", sep = "")
#cairo_pdf(file = pdfFN)
pdf(pdfFN, height = (length(leafColors)/12 + 16), width = (length(leafColors)/12 + 16))

heatmap.2(as.matrix(mciDist), col=colorRampPalette(c("darkgreen", "white", "darkorchid"))(550), trace="none", Rowv=mciDen, Colv=rev(mciDen), margins=c(16,16))
dev.off()

pam.2 <- pam(mciDist, 2) 
conOut <- paste(outDir, "/", paste("combined-pam2clusters-reg-", repCat, sep = "-"), ".txt", sep = "")
pamClustRes <- pam.2$clustering
sink(file = conOut)
pam.2$clusinfo
cbind(names(pamClustRes), pamClustRes)
sink()

bed <- matrix(unlist(strsplit(gsub("[:-]", " ", names(pam.2$clustering)), " +")), ncol = 4, byrow=TRUE)  
bedNames <- as.numeric(gsub('^.{3}', '', bed[,4]))
matchingState <- bedNames == pam.2$clustering

bed <- unname(cbind(bed[,-4], paste(bedNameConverts[bedNames] ,"_cluster-", pam.2$clustering, sep=""), as.numeric(matchingState)))
bedOut <- paste(outDir, "/", paste("pam_clustering-k2_", repCat, sep = "-"), ".bed", sep = "")
write.table(bed, file=bedOut, col.names=F, row.names=F, sep="\t", quote=FALSE)

pdfFN <- paste(outDir, "/", paste("pam_clustering-k2_sillouhettePlot-2_regions", repCat, sep = "-"), ".pdf", sep = "")
pdf(pdfFN)
plot(pam.2)
dev.off()


silKRes <- mclapply(2:min(10, as.integer(length(matchingState)/10)), pamSpam, distMat = mciDist, mc.cores=nCores, mc.preschedule = FALSE)
silKData <- unlist(silKRes)

bestSilK <- which(max(silKData) == silKData) + 1

pdfFN <- paste(outDir, "/", paste("pam_sillouhetteScore-v-K-Best_regions_allIntervals", repCat, sep = "-"), ".pdf", sep = "")
pdf(pdfFN)
plot(1:length(silKRes), silKRes, xlab = "Number of clusters", ylab = "Silhouette Width", main=paste("Best k is: ", bestSilK, sep=""))
lines(1:length(silKRes), silKRes)
dev.off()


pam_fit <- pam(mciDist, bestSilK) 
tsne_obj <- Rtsne(mciDist, is_distance = TRUE)

tsne_data <- tsne_obj$Y %>% data.frame() %>% setNames(c("TSNE_1", "TSNE_2")) %>% mutate(clust = factor(pam_fit$clustering), name = rownames(matCIntervals))

tsne_data$bed <- c(rep(bed1Name, length(bed1Reg)), rep(bed2Name, length(bed2Reg)))
clustNames <- paste(rep("Cluster-", bestSilK), 1:bestSilK, sep="")

pdfFN <- paste(outDir, "/", paste("tsnePlot-BestK_regions_allIntervals", repCat, sep = "-"), ".pdf", sep = "")
pdf(pdfFN)

tsneGG <- ggplot(tsne_data, aes(x = TSNE_1, y = TSNE_2, color = bed)) + labs(title=paste("TSNE of intervals at the ", repCat, " level", sep=""))
tsneGG <- tsneGG + scale_colour_manual(name="Source Files", labels = c(bed1Name, bed2Name), values = c("blue", "red"))
tsneGG <- tsneGG + geom_point(aes(shape=clust))
tsneGG <- tsneGG + scale_shape_discrete(name="PAM Cluster")#, labels=clustNames)
tsneGG <- tsneGG + theme(legend.position="bottom")
tsneGG <- tsneGG + guides(shape = guide_legend(nrow=2,byrow=TRUE), colour = guide_legend(nrow=2,byrow=TRUE))
tsneGG <- tsneGG + theme(plot.title = element_text(hjust = 0.5))
plot(tsneGG)

dev.off()

out <- paste(outDir, "/", paste("tsneResults_regions", repCat, sep = "-"), ".tsv", sep = "")
write.table(tsne_data, file=out, col.names=T, row.names=F, sep="\t", quote=FALSE)


pdfFN <- paste(outDir, "/", paste("metaMDS_ordiplot_allIntervals_pam.2", repCat, sep = "-"), ".pdf", sep = "")
pdf(pdfFN)
ordiplot(intCMDS)
ordihull(intCMDS, groups=pam.2$clustering, draw='lines', col="darkgrey", label=F)

clustColors <- brewer.pal(max(length(unique(pam.2$clustering)),3), "Dark2")
for(icml in unique(pam.2$clustering)) {
  ordiellipse(intCMDS$point[which(icml == pam.2$clustering),],draw="polygon",
  groups=pam.2$clustering[pam.2$clustering==icml],col=clustColors[icml],label=F, conf=0.95, alpha=30) 
} 
orditorp(intCMDS,display="sites",col=c(rep("blue", length(bed1Reg)), rep("red", length(bed2Reg))), labels=c(F), pch=1)

z <- kde2d(intCMDS$points[,1], intCMDS$points[,2], n=50)
k <- 10
my.cols <- rev(brewer.pal(k, "Spectral"))
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE, lwd=1)
dev.off()


pdfFN <- paste(outDir, "/", paste("metaMDS_ordiplot_allIntervals_pam.best", repCat, sep = "-"), ".pdf", sep = "")
pdf(pdfFN)
ordiplot(intCMDS)
ordihull(intCMDS, groups=pam_fit$clustering, draw='lines', col="darkgrey", label=F)

clustColors <- brewer.pal(max(length(unique(pam_fit$clustering)),3), "Dark2")
for(icml in unique(pam_fit$clustering)) {
  ordiellipse(intCMDS$point[which(icml == pam_fit$clustering),],draw="polygon",
  groups=pam_fit$clustering[pam_fit$clustering==icml],col=clustColors[icml],label=F, conf=0.95, alpha=30) 
} 
orditorp(intCMDS,display="sites",col=c(rep("blue", length(bed1Reg)), rep("red", length(bed2Reg))), labels=c(F), pch=1)

z <- kde2d(intCMDS$points[,1], intCMDS$points[,2], n=50)
k <- 10
my.cols <- rev(brewer.pal(k, "Spectral"))
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE, lwd=1)
dev.off()



###################

mat1Repeats <- data.frame(matrix(NA, nrow = repLen, ncol=0))
mat2Repeats <- data.frame(matrix(NA, nrow = repLen, ncol=0))

mtRes <- c("", "", "")
for(sm in 1:length(statsTypes)){
	###load in and add labels to data
	b1FN <- paste(paste(statsTypes[sm], repCat, "b1", sep = "-"), ".tsv", sep = "")
	b1Mat <- as.data.frame(fread(b1FN, header=FALSE, sep = "\t"))
	names(b1Mat) <- reps
	row.names(b1Mat) <- bed1Reg
	b1Mat <- t(b1Mat)
	if(sm == 1 || sm == 2 || sm == 3){
		b1Mat <- 100 - b1Mat #might help with the NaN stuff
	}
	colnames(b1Mat) <- paste(statsTypes[sm], colnames(b1Mat), sep = "_-_")
	b1Mat[is.na(b1Mat)] <- 0
	mat1Repeats <- cbind(mat1Repeats, b1Mat)
	
	b2FN <- paste(paste(statsTypes[sm], repCat, "b2", sep = "-"), ".tsv", sep = "")
	b2Mat <- as.data.frame(fread(b2FN, header=FALSE, sep = "\t"))
	names(b2Mat) <- reps
	row.names(b2Mat) <- bed2Reg
	b2Mat <- t(b2Mat)
	if(sm == 1 || sm == 2 || sm == 3){
		b2Mat <- 100 - b2Mat #might help with the NaN stuff
	}
	colnames(b2Mat) <- paste(statsTypes[sm], colnames(b2Mat), sep = "_-_")
	b2Mat[is.na(b2Mat)] <- 0
	mat2Repeats <- cbind(mat2Repeats, b2Mat)

	if(1){#high of a chance for failure with bray
		m1dist <- getDist(b1Mat, distType)
		m2dist <- getDist(b2Mat, distType)
		mt12res <- mantel(m1dist, m2dist, method="pearson", permutations = 100000, parallel = nCores)
		
		mtRes <- rbind( mtRes, c(humanStats[sm], mt12res$signif, mt12res$statistic ))
		pdfFN <- paste(outDir, "/mantel/", paste("mantelDistancePlots-", statsTypes[sm], repCat, sep = "-"), ".pdf", sep = "")
		pdf(pdfFN)
		plot(m1dist, m2dist)
		dev.off()

	}
}

m1dist <- getDist(mat1Repeats, distType)
m2dist <- getDist(mat2Repeats, distType)

#mtRes <- mantel(m1dist, m2dist, method="pearson", permutations = 100000, parallel = nCores)
mt12res <- mantel(m1dist, m2dist, method="pearson", permutations = 100000, parallel = nCores)
mtRes <- rbind( mtRes, c("all", mt12res$signif, mt12res$statistic ))
mtResM <- matrix(unlist(mtRes), ncol = 3, byrow=TRUE)

mtRes <- rbind( mtRes, c(humanStats[sm], mt12res$signif, mt12res$statistic ))
pdfFN <- paste(outDir, "/mantel/", paste("mantelDistancePlots-all", repCat, sep = "-"), ".pdf", sep = "")
pdf(pdfFN)
plot(m1dist, m2dist)
dev.off()

out <- paste(outDir, "/", paste("mantel-tests", repCat, sep = "-"), ".tsv", sep = "")
write.table(mtRes, file=out, col.names=F, row.names=F, sep="\t", quote=FALSE)

m1Den <- as.dendrogram(hclust(m1dist, "ward.D2"))
m2Den <- as.dendrogram(hclust(m2dist, "ward.D2"))
dend_list <- dendlist(m1Den, m2Den)

pdfFN <- paste(outDir, "/", paste("tanglegram_repeats", repCat, sep = "-"), ".pdf", sep = "")
pdf(pdfFN)
tanglegram(m1Den, m2Den, margin_inner= 7, center = TRUE,
  highlight_distinct_edges = FALSE, # Turn-on dashed lines
  common_subtrees_color_lines = TRUE, # Turn-on line colors
  common_subtrees_color_branches = TRUE, # Color common branches 
  main = paste("entanglement =", round(entanglement(dend_list), 2)),
  main_left = bed1Name, main_right = bed2Name )#, xlim = c(7,0))
dev.off()


hmNumColors <- min(length(reps)*length(reps), 500)
pdfFN <- paste(outDir, "/", paste("heatmap_repeats_bed1", repCat, sep = "-"), ".pdf", sep = "")
#cairo_pdf(file = pdfFN)
pdf(pdfFN, height = (length(reps)/12 + 16), width = (length(reps)/12 + 16))
heatmap.2(as.matrix(m1dist), col=colorRampPalette(c("darkgreen", "white", "darkorchid"))(hmNumColors), trace="none", Rowv=m1Den, Colv=rev(m1Den), margins=c(16,16))
dev.off()

pdfFN <- paste(outDir, "/", paste("heatmap_repeats_bed2", repCat, sep = "-"), ".pdf", sep = "")
#cairo_pdf(file = pdfFN)
pdf(pdfFN, height = (length(reps)/12 + 16), width = (length(reps)/12 + 16))
heatmap.2(as.matrix(m2dist), col=colorRampPalette(c("darkgreen", "white", "darkorchid"))(hmNumColors), trace="none", Rowv=m2Den, Colv=rev(m2Den), margins=c(16,16))
dev.off()

m12Diff <- as.dist(as.matrix(m1dist) - as.matrix(m2dist))
m12Dend <- as.dendrogram(hclust(m12Diff, "ward.D2"))

pdfFN <- paste(outDir, "/", paste("heatmap_repeats_bed1Dists-bed2Dists", repCat, sep = "-"), ".pdf", sep = "")
#cairo_pdf(file = pdfFN)
pdf(pdfFN, height = (length(reps)/12 + 16), width = (length(reps)/12 + 16))
heatmap.2(as.matrix(m12Diff), col=colorRampPalette(c("darkgreen", "white", "darkorchid"))(hmNumColors), trace="none", Rowv=m12Dend, Colv=rev(m12Dend), margins=c(16,16))
dev.off()

###############
###############
###############
}
###############
###############
###############





cn1 <- matrix(unlist(strsplit(bed1Reg, "--")), ncol=2, byrow=TRUE)
mat1DF <- cbind(rep(reps, each=length(bed1Reg)), cn1[,1], cn1[,2], rep("", (length(reps) * length(bed1Reg))))
mat1DF <- cbind(paste(mat1DF[,1], mat1DF[,2], mat1DF[,3], sep="__"), mat1DF[,-4])
rownames(mat1DF) <- mat1DF[,1]
colnames(mat1DF) <- c("uid", "reps", "interval", "bed")

cn2 <- matrix(unlist(strsplit(bed2Reg, "--")), ncol=2, byrow=TRUE)
mat2DF <- cbind(rep(reps, each=length(bed2Reg)), cn2[,1], cn2[,2], rep("", (length(reps) * length(bed2Reg))))
mat2DF <- cbind(paste(mat2DF[,1], mat2DF[,2], mat2DF[,3], sep="__"), mat2DF[,-4])
rownames(mat2DF) <- mat2DF[,1]
colnames(mat2DF) <- c("uid", "reps", "interval", "bed")

for(sm in 1:length(statsTypes)){
	###load in and add labels to data
	b1FN <- paste(paste(statsTypes[sm], repCat, "b1", sep = "-"), ".tsv", sep = "")
	b1Mat <- as.data.frame(fread(b1FN, header=FALSE, sep = "\t"))
	m1Vec <- as.vector(as.matrix(b1Mat))
	names(m1Vec) <- rownames(mat1DF)
	cn1 <- colnames(mat1DF)

	m1Vec[is.na(m1Vec)] <- 0
	mat1DF <- cbind(m1Vec, mat1DF)
	colnames(mat1DF) <- c(statsTypes[sm],cn1)
	
	b2FN <- paste(paste(statsTypes[sm], repCat, "b2", sep = "-"), ".tsv", sep = "")
	b2Mat <- as.data.frame(fread(b2FN, header=FALSE, sep = "\t"))
	m2Vec <- as.vector(as.matrix(b2Mat))
	names(m2Vec) <- rownames(mat2DF)
	cn2 <- colnames(mat2DF)
	if(sm == 1 || sm == 2 || sm == 3){
		m2Vec <- 100 - m2Vec
	}
	m2Vec[is.na(m2Vec)] <- 0
	mat2DF <- cbind(m2Vec, mat2DF)
	colnames(mat2DF) <- c(statsTypes[sm],cn2)
}
matDF <- rbind(mat1DF, mat2DF)
matDF <- as.data.frame(matDF)


dataMat <- data.matrix(matDF[, 1:(dim(matDF)[2]-4)])




bedSources <- c(rep("bed1", dim(mat1Repeats)[2]), rep("bed2", dim(mat2Repeats)[2]))
comm<-as.data.frame(t(cbind(mat1Repeats, mat2Repeats)))
labels <- cbind(bedSources, statsTypes, humanStats)
rownames(labels) <- rownames(comm)
colnames(labels) <- c("bed", "parameter", "human")
labels <- as.data.frame(labels)


comDist <- parDist(as.matrix(comm), method=distMethodInUse, threads=nCores)
comDist <- as.matrix(comDist)

rownames(comDist) <- rownames(comm)
colnames(comDist) <- rownames(comm)
comDist <- as.dist(comDist)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

interaction.plot(labels$bed, labels$human, comm[,ip], col=cbbPalette, lwd=2, fixed=TRUE, ylab=paste("mean of ", colnames(comm)[ip], sep=""), xlab="Source Bed File", trace.label="Tracked Statistics")


















#################################
#################################

cn1 <- matrix(unlist(strsplit(bed1Reg, "--")), ncol=2, byrow=TRUE)
blank1 <- rep("", length(bed1Reg))
mat1DF <- cbind(bed1Reg, cn1[,1], cn1[,2])
mat1DF <- cbind(blank1, paste(mat1DF[,1], mat1DF[,2], mat1DF[,3], sep="__"), mat1DF)
rownames(mat1DF) <- mat1DF[,2]
colnames(mat1DF) <- c("uid", "interval", "bed", "metric")
mat1Base <- mat1DF[,c(2,3,4,5)]
mat1DF <- c("", "", "", "", "", "")

cn2 <- matrix(unlist(strsplit(bed2Reg, "--")), ncol=2, byrow=TRUE)
blank2 <- rep("", (length(reps) * length(bed2Reg)))
mat2DF <- cbind(rep(reps, each=length(bed2Reg)), cn2[,1], cn2[,2], rep("", (length(reps) * length(bed2Reg))))
mat2DF <- cbind(blank2, paste(mat2DF[,1], mat2DF[,2], mat2DF[,3], sep="__"), mat2DF)
rownames(mat2DF) <- mat2DF[,2]
colnames(mat2DF) <- c("value", "uid", "reps", "interval", "bed", "metric")
mat2Base <- mat2DF[,c(2,3,4,5)]
mat2DF <- c("", "", "", "", "", "")

for(sm in 1:length(statsTypes)){
	###load in and add labels to data
	b1FN <- paste(paste(statsTypes[sm], repCat, "b1", sep = "-"), ".tsv", sep = "")
	b1Mat <- as.data.frame(fread(b1FN, header=FALSE, sep = "\t"))
	if(sm == 1 || sm == 2 || sm == 3){
		b1Mat <- 100 - b1Mat
	}
	cs1Mat <- decostand(b1Mat, method="chi.square", na.rm = T)
	m1Vec <- as.vector(as.matrix(cs1Mat))
	names(m1Vec) <- rownames(mat1Base)
	cn1 <- c(colnames(mat1Base))
	
	m1Vec[is.na(m1Vec)] <- 0
	v1Mat <- cbind(m1Vec, statsTypes[sm])
	mat1DF <- rbind(mat1DF, cbind(v1Mat[,1], mat1Base, v1Mat[,2]))
	colnames(mat1DF) <- c("value", cn1, "metric")
	
	b2FN <- paste(paste(statsTypes[sm], repCat, "b2", sep = "-"), ".tsv", sep = "")
	b2Mat <- as.data.frame(fread(b2FN, header=FALSE, sep = "\t"))
	if(sm == 1 || sm == 2 || sm == 3){
		b2Mat <- 100 - b2Mat
	}
	cs2Mat <- decostand(b1Mat, method="chi.square", na.rm = T)
	m2Vec <- as.vector(as.matrix(cs2Mat))
	names(m2Vec) <- rownames(mat2Base)
	cn2 <- c(colnames(mat2Base))

	m2Vec[is.na(m2Vec)] <- 0
	v2Mat <- cbind(m2Vec, mat2DF)
	mat2DF <- rbind(mat2DF, cbind(v2Mat[,1], mat2Base, v2Mat[,2]))
	colnames(mat2DF) <- c("value", cn2, "metric")
}

matDF <- rbind(mat1DF[-1,], mat2DF[-1,])
matDF <- as.data.frame(matDF)
#################################
#################################
cn1 <- matrix(unlist(strsplit(bed1Reg, "--")), ncol=2, byrow=TRUE)
blank1 <- rep("", (length(reps) * length(bed1Reg)))
mat1DF <- cbind(rep(reps, each=length(bed1Reg)), cn1[,1], cn1[,2], blank1)
mat1DF <- cbind(blank1, paste(mat1DF[,1], mat1DF[,2], mat1DF[,3], sep="__"), mat1DF)
rownames(mat1DF) <- mat1DF[,2]
colnames(mat1DF) <- c("value", "uid", "reps", "interval", "bed", "metric")
mat1Base <- mat1DF[,c(2,3,4,5)]
mat1DF <- c("", "", "", "", "", "")

cn2 <- matrix(unlist(strsplit(bed2Reg, "--")), ncol=2, byrow=TRUE)
blank2 <- rep("", (length(reps) * length(bed2Reg)))
mat2DF <- cbind(rep(reps, each=length(bed2Reg)), cn2[,1], cn2[,2], rep("", (length(reps) * length(bed2Reg))))
mat2DF <- cbind(blank2, paste(mat2DF[,1], mat2DF[,2], mat2DF[,3], sep="__"), mat2DF)
rownames(mat2DF) <- mat2DF[,2]
colnames(mat2DF) <- c("value", "uid", "reps", "interval", "bed", "metric")
mat2Base <- mat2DF[,c(2,3,4,5)]
mat2DF <- c("", "", "", "", "", "")

for(sm in 1:length(statsTypes)){
	###load in and add labels to data
	b1FN <- paste(paste(statsTypes[sm], repCat, "b1", sep = "-"), ".tsv", sep = "")
	b1Mat <- as.data.frame(fread(b1FN, header=FALSE, sep = "\t"))
	if(sm == 1 || sm == 2 || sm == 3){
		b1Mat <- 100 - b1Mat
	}
	cs1Mat <- decostand(b1Mat, method="chi.square", na.rm = T)
	m1Vec <- as.vector(as.matrix(cs1Mat))
	names(m1Vec) <- rownames(mat1Base)
	cn1 <- c(colnames(mat1Base))
	
	m1Vec[is.na(m1Vec)] <- 0
	v1Mat <- cbind(m1Vec, statsTypes[sm])
	mat1DF <- rbind(mat1DF, cbind(v1Mat[,1], mat1Base, v1Mat[,2]))
	colnames(mat1DF) <- c("value", cn1, "metric")
	
	b2FN <- paste(paste(statsTypes[sm], repCat, "b2", sep = "-"), ".tsv", sep = "")
	b2Mat <- as.data.frame(fread(b2FN, header=FALSE, sep = "\t"))
	if(sm == 1 || sm == 2 || sm == 3){
		b2Mat <- 100 - b2Mat
	}
	cs2Mat <- decostand(b1Mat, method="chi.square", na.rm = T)
	m2Vec <- as.vector(as.matrix(cs2Mat))
	names(m2Vec) <- rownames(mat2Base)
	cn2 <- c(colnames(mat2Base))

	m2Vec[is.na(m2Vec)] <- 0
	v2Mat <- cbind(m2Vec, mat2DF)
	mat2DF <- rbind(mat2DF, cbind(v2Mat[,1], mat2Base, v2Mat[,2]))
	colnames(mat2DF) <- c("value", cn2, "metric")
}

matDF <- rbind(mat1DF[-1,], mat2DF[-1,])
matDF <- as.data.frame(matDF)
#################################

repIntComMat1 <- data.frame(matrix(NA, ncol = length(reps), nrow=length(bed1Reg)))
repIntComMat2 <- data.frame(matrix(NA, ncol = length(reps), nrow=length(bed2Reg)))

for(sm in 1:length(statsTypes)){
	###load in and add labels to data
	b1FN <- paste(paste(statsTypes[sm], repCat, "b1", sep = "-"), ".tsv", sep = "")
	b1Mat <- as.data.frame(fread(b1FN, header=FALSE, sep = "\t"))
	colnames(b1Mat) <- reps
	row.names(b1Mat) <- bed1Reg
	b1Mat <- as.matrix(b1Mat)
	
	if(sm == 1 || sm == 2 || sm == 3){
		b1Mat <- 100 - b1Mat
	}
	repIntComMat1 <- abind(repIntComMat1, b1Mat, along=3)
	
	b2FN <- paste(paste(statsTypes[sm], repCat, "b2", sep = "-"), ".tsv", sep = "")
	b2Mat <- as.data.frame(fread(b2FN, header=FALSE, sep = "\t"))
	colnames(b2Mat) <- reps
	row.names(b2Mat) <- bed2Reg
	b2Mat <- as.matrix(b2Mat)
	
	if(sm == 1 || sm == 2 || sm == 3){
		b2Mat <- 100 - b2Mat
	}
	repIntComMat2 <- abind(repIntComMat2, b2Mat, along=3)
}
#remove the extra bit generated from the creation of the 3d array
repIntComMat1 <- repIntComMat1[,,-1]
repIntComMat2 <- repIntComMat2[,,-1]

suppress <- mclapply(1:length(reps), bedSimWrites, repIntComMat1 = repIntComMat1, repIntComMat2 = repIntComMat2, mc.cores=nCores, mc.preschedule = FALSE)

repIntMat <- rbind(repIntComMat1, repIntComMat2)
repIntMatLabels <- c(rep("Bed1", length(rowLabels1)), rep("Bed2", length(rowLabels2)))
repIntMatColors <- c(rep("blue", length(rowLabels1)), rep("red", length(rowLabels2)))

if(switchDist){
	mcdist <- distance(repIntMat, method="gower")
}else{
	mcdist <- parDist(as.matrix(repIntMat), method=distMethodInUse, threads=nCores)
	mcdist <- as.matrix(mcdist)
}
rownames(mcdist) <- rownames(repIntMat)
colnames(mcdist) <- rownames(repIntMat)


mcdist <- as.dist(mcdist)




























###################################


mc1dist <- distance(matCIntervals, method="gower")
rownames(mc1dist) <- rownames(matCIntervals)
colnames(mc1dist) <- rownames(matCIntervals)
mc1dist <- as.dist(mc1dist)
m1Den <- as.dendrogram(hclust(mc1dist, "ward.D2"))

mc2dist <- distance(matCIntervals, method="sorensen")
rownames(mc2dist) <- rownames(matCIntervals)
colnames(mc2dist) <- rownames(matCIntervals)
mc2dist <- as.dist(mc2dist)
m2Den <- as.dendrogram(hclust(mc2dist, "ward.D2"))

dend_list <- dendlist(m1Den, m2Den)

pdfFN <- paste(outDir, "/", paste("tanglegram_repeats", repCat, sep = "-"), ".pdf", sep = "")
pdf(pdfFN)
tanglegram(m1Den, m2Den, margin_inner= 7, center = TRUE,
  highlight_distinct_edges = FALSE, # Turn-on dashed lines
  common_subtrees_color_lines = TRUE, # Turn-on line colors
  common_subtrees_color_branches = TRUE, # Color common branches 
  main = paste("entanglement =", round(entanglement(dend_list), 2)),
  main_left = bed1Name, main_right = bed2Name )#, xlim = c(7,0))
dev.off()






#####################



ggdendrogram(as.dendrogram(hclust(mciUsed, hclustMethodInUse)), rotate=TRUE)

repDist.bray <- vegdist(mat2Repeats, method = "bray")

meths <- getDistMethods()
res <- c()
metRes <- c()
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

# function to compute coefficient
ac <- function(x) {
  agnes(mciTesting, method = x)$ac
}
#meths <- c("sorensen", "soergel",  "gower")
for( j in 1:length(meths)){
	met <- meths[i]
	met <- "gower"
	if(i == 3){
		mciTesting <- distance(mat1Repeats, method=met, p=3)
	}else{
		mciTesting <- distance(mat1Repeats, method=met)
	}
	
	rownames(mciTesting) <- rownames(mat1Repeats)
	colnames(mciTesting) <- rownames(mat1Repeats)
	mciTesting <- as.dist(mciTesting)

	mciTree <- as.dendrogram(hclust(mciTesting, 'ward.D2'))
	ggd <- ggdendrogram(mciTree, rotate = TRUE)
	ggd <- ggd + labs(title = met)
	plot(ggd)
    Sys.sleep(3)
    
	#intCMDS <- metaMDS(mciTesting, autotransform=FALSE, k=2, trymax=500, parallel=4)
	#stressplot(intCMDS)
	#Sys.sleep(3)
	#res <- c(res, c(met, intCMDS$stress))
	# methods to assess
	
	plot(fviz_dist(mciTesting, gradient = list(low = "blue", mid = "white", high = "red", ordered=T)) + theme(text = element_text(size = 8)))

    Sys.sleep(3)
}


rClust <- as.dendrogram(hclust(dist(matCIntervals, method="minkowski"), 'ward.D2'))
cClust <- as.dendrogram(hclust(dist(t(matCIntervals), method="minkowski"), 'ward.D2'))

nbs<-1000
rc_result <- pvclust(t(matCIntervals), method.dist="minkowski", method.hclust="ward.D2", nboot=nbs, parallel=as.integer(nCores))
#out <- paste(outDir, "/", paste("bed2-pvclust-reg", statsTypes[l], repCat, sep = "-"), ".pdf", sep = "")
#pdf(out, height = (1*regTreeSize + 10), width = ((log(regTreeSize) * 200) + 10))
plot(rc_result)
pvrect(rc_result, alpha=0.95)


heatmap.2(as.matrix(data.prop.1),Rowv = rClust, Colv = cClust, col = scaleyellowred, RowSideColors = var1)

heatmap.2(as.dm, col = heat.colors(256), scale = coloring[cs], cexRow = 3, cexCol = 3, lhei = c(3,regTreeSize), lwid = c(2,repTreeSize), trace="none", Rowv = rClust, Colv = cClust)





matCRepeats <- cbind(mat1Repeats, mat2Repeats)
repDist.bray <- vegdist(matCRepeats, method = "bray")
plot(hclust(repDist.bray, method="ward.D2")


dmDist.bray <- vegdist(matC, method = "jaccard")
dmDistT.bray <- vegdist(t(dm), method = "jaccard")
dmDistT.bray.r <- replace(dmDistT.bray, is.na(dmDistT.bray), 0)
#col.clust <- hclust(dmDist.bray, "ward.D2")
#row.clust <- hclust(dmDistT.bray.r, "ward.D2")	
col.clust <- hclust(dmDist.bray, method = "average")
row.clust <- hclust(dmDistT.bray.r, method = "average")


###create the comvined matrix that contains all of the data
matC <- rbind(mat1, mat2)
comMatLabels <- combinedMat[,(repLen+1)]
noLabelMat <- combinedMat[,1:repLen]

b2MatCat <- rep("Bed2", length(bed2Reg))
	b2Mat$Classification <- b2MatCat


b1MatCat <- rep("Bed1", length(bed1Reg))
	b1Mat$Classification <- b1MatCat
