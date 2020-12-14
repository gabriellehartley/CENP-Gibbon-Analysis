setwd("repStatI-CHPO-85-MNASE_vs_CHPO-86-MNASE_500kReads")

nCores <- 2
goIndividual <- 1
goCluster <- 0
goClassify <- 0
goPies <- 1
goPVClust <- 1 #1
pvParaCore <- 1 #3
bed1Name <- c("SRR5723785-500k.bed")
bed2Name <- c("SRR5723786-500k.bed")


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


#following variables usesful for debugging
#i<-3
#l<-1
#st<-1

b1File <- c("bed1Stats.tsv")
b2File <- c("bed2Stats.tsv")
#read in data
b1FullData <- as.data.frame(fread(b1File, header=FALSE, sep = "\t"))
b2FullData <- as.data.frame(fread(b2File, header=FALSE, sep = "\t"))

catLabels <- c("species", "class", "family")
repeatFileList <- c("species-repeatList.tsv", "class-repeatList.tsv", "family-repeatList.tsv")

#bedDataColNames <- c("Repeat", "Species", "Family", "Class", "Region", "Region_Length", "Scaffold", "Divergence", "Deleted", "Inserted", "Length", "DNA_length_bed1", "DNA_length_bed2", "Rep_length_bed1", "Rep_length_bed2", "DNA_length_Total", "Rep_length_Total")
#ggplot(b1FullData, aes(x=Inserted,y=Length)) + geom_point() + geom_smooth() + facet_grid(Family ~ Region)

#columnLabels <- as.data.frame(fread("columnLabels.tsv", header=FALSE, sep = "\t"))

regionList <- as.data.frame(fread("intervalLabels-c.tsv", header=FALSE, sep = "\t"))
#regionList <- t(region List)
regionList <- sapply(regionList, as.character)
#get the regions that are overlapped, (search for the ::: in contig:::start-stop)
rmRG <- grepl ("--bed1", regionList)
bed1Reg <- regionList[rmRG]
bed2Reg <- regionList[!rmRG]
#regionList <- t(regionList)


individualApply <- function(curCombo, allCombos, statsTypes, statsNames, repLevel, b1, b2){
  combo <- allCombos[curCombo,]
  
  stat <- as.integer(paste(unlist(combo[2]), collapse=""))
  repName <- paste(unlist(combo[1]), collapse="")
  
  print(combo)
  
  rmb1 <- b1[which(b1[,1] == repName),stat+1]
  rmb2 <- b2[which(b2[,1] == repName),stat+1]
  
  if(length(rmb1) < 1 || length(rmb2) < 1){
    b1Len <- length(rmb1)
    b2Len <- length(rmb2)
    b1mean <- 0
    b2mean <- 0

    return(c(stat,repName,b1Len,b2Len,0,0,0,0))
  }
  if(identical(unique(rmb1), (unique(rmb2)))){
    b1Len <- length(rmb1)
    b2Len <- length(rmb2)
    return(c(stat,repName,b1Len,b2Len,1,1,1,1))
  }
  b1CDF <- ecdf(rmb1)
  b2CDF <- ecdf(rmb2)

  #plot combined
  adRes <- ad.test(rmb1, rmb2, method= "simulated", Nsim = 500, dist=TRUE)
  ksResG <- ks.boot(rmb1, rmb2, nboots=100, alternative = "g")
  ksResL <- ks.boot(rmb1, rmb2, nboots=100, alternative = "l")

  #consider looking into Cramer-von Mises, Wilcoxon-Mann-Whitney and ryan-joiner

  pdf(file = paste("individualRepeatTests/", repLevel, "-", statsTypes[stat],"/cdf/CDF-", statsTypes[stat], "-for-", repName, ".pdf", sep = ""))
  plot(b1CDF, col="blue", main= paste(statsNames[stat], " for ", repName, " (", repLevel, ") ", sep = ""), xlab=humanName)
  lines(b2CDF, col="red")
  legend(x="bottomright",legend= paste("blue =", bed1Name, " regions\nred =", bed2Name , " regions", sep = " "))
  dev.off()

  pdf(file = paste("individualRepeatTests/", repLevel, "-", statsTypes[stat],"/violin/violin-", statsTypes[stat], "-for-", repName, ".pdf", sep = ""))
  b1df <- as.data.frame(rmb1)
  b1df$type <- bed1Name
  colnames(b1df) <- c("stat", "bed")
  b2df <- as.data.frame(rmb2)
  b2df$type <- bed2Name
  colnames(b2df) <- c("stat", "bed")
  bcdf <- rbind(b1df, b2df)

  p <- ggplot(bcdf, aes(x=bed, y=stat, color=bed)) + ylab(humanName) + geom_violin()
  p <- p + stat_summary(fun.data="mean_sdl", geom="pointrange")
  p <- p+scale_color_manual(values=c("#0000FF", "#FF0000"))
  print(p)
  dev.off()
  

  b1Len <- length(rmb1)
  b2Len <- length(rmb2)
  #colnames(resultStats) <- c("repeat", "ROI counts", "background counts", "KS lower p-value", "KS greater p-value", "AD Cont. p-value", "AD Disc. p-value")
  results <- c(stat,repName,b1Len,b2Len,ksResL$ks.boot.pvalue,ksResG$ks.boot.pvalue,adRes$ad[1,3],adRes$ad[2,3])
  return(results)
}

clAnalysesApply <- function(l, statsTypes, repCat, reps, bed1Reg, bed2Reg, repLen){
  outDir <- paste("IntervalAnalyses/", statsTypes[l], sep="")
  dir.create(outDir)

  #load in and add labels to data
  b1FN <- paste(paste(statsTypes[l], repCat, "b1", sep = "-"), ".tsv", sep = "")
  b1Mat <- as.data.frame(fread(b1FN, header=FALSE, sep = "\t"))
  names(b1Mat) <- reps
  row.names(b1Mat) <- bed1Reg
  b1MatCat <- rep("Bed1", length(bed1Reg))
  b1Mat$Classification <- b1MatCat

  b2FN <- paste(paste(statsTypes[l], repCat, "b2", sep = "-"), ".tsv", sep = "")
  b2Mat <- as.data.frame(fread(b2FN, header=FALSE, sep = "\t"))
  names(b2Mat) <- reps
  row.names(b2Mat) <- bed2Reg
  b2MatCat <- rep("Bed2", length(bed2Reg))
  b2Mat$Classification <- b2MatCat

  #create the comvined matrix that contains all of the data
  combinedMat <- rbind(b1Mat,b2Mat)
  comMatLabels <- combinedMat[,(repLen+1)]
  noLabelMat <- combinedMat[,1:repLen]

  #remove cols without any variance for ggbiplot
  cMat <- noLabelMat[,apply(noLabelMat, 2, var, na.rm=TRUE) != 0]

  #run the MDS/PCA analyses
  if(goCluster){
    comMatPCA <- prcomp(cMat)

    #pca for unlabeled spots
    pdfFN <- paste(outDir, "/", paste(statsTypes[l], repCat, sep = "-"), ".pdf", sep = "")
    pdf(file = pdfFN)
    ggbp <- ggbiplot(comMatPCA, obs.scale = 1, var.scale = 1,
                     groups = comMatLabels, ellipse = TRUE, circle = TRUE) +
      scale_color_discrete(name = "") +
      theme(legend.direction = "horizontal", legend.position = "top")
    print(ggbp)
    dev.off()

    #pca for bed1 labeled
    rncm <- rownames(combinedMat)
    rmRG <- grepl ("--bed1", rncm)
    b1names <- rncm[rmRG]
    b2names <- rncm[!rmRG]
    b2l <- length(b2names)
    regName <- c(b1names, rep("", b2l))
    regName <- gsub("--bed1", "", regName)

    pdfFN <- paste(outDir, "/", paste(statsTypes[l], repCat, "bed1Named", sep = "-"), ".pdf", sep = "")
    pdf(file = pdfFN)
    ggbp <- ggbiplot(comMatPCA, obs.scale = 1, var.scale = 1,
                     groups = comMatLabels, ellipse = TRUE, circle = TRUE, labels = regName) +
      scale_color_discrete(name = "") +
      theme(legend.direction = "horizontal", legend.position = "top")
    print(ggbp)
    dev.off()

    #pca for bed2 labeled
    rncm <- rownames(combinedMat)
    rmRG <- grepl ("--bed1", rncm)
    b1names <- rncm[rmRG]
    b2names <- rncm[!rmRG]
    b1l <- length(b1names)
    regName <- c(rep("", b1l), b2names)
    regName <- gsub("--bed2", "", regName)

    pdfFN <- paste(outDir, "/", paste(statsTypes[l], repCat, "bed2Named", sep = "-"), ".pdf", sep = "")
    pdf(file = pdfFN)
    ggbp <- ggbiplot(comMatPCA, obs.scale = 1, var.scale = 1,
                     groups = comMatLabels, ellipse = TRUE, circle = TRUE, labels = regName) +
      scale_color_discrete(name = "") +
      theme(legend.direction = "horizontal", legend.position = "top")
    print(ggbp)
    dev.off()
  }

  coloring <- c("none","row", "column")

  #bed1
  #bed1
  #bed1
  #bed1

  dm <- b1Mat[,1:ncol(b1Mat)-1]
  dm <- data.matrix(dm)

  #dataframe used for multiple analyses
  df1 <- as.data.frame(b1Mat[,1:ncol(b1Mat)-1])

  repTreeSize <- as.integer(ncol(dm)/3)
  regTreeSize <- as.integer(nrow(dm)/2.5)
  nbs=100

  #bulk of bed1 processing
  if(nrow(dm) > 1 && ncol(dm) > 1){ #some sanity checks for overly small datasets

    if(goCluster){
      pdf(NULL)
      hm <- heatmap.2(dm)
      dev.off()
      hc <- as.hclust(hm$colDendrogram)
      t <- as.phylo(hc)
      treeOut <- paste(outDir, "/", paste("bed1", "repeats", statsTypes[l], repCat, sep = "-"), ".newick", sep = "")
      write.tree(phy=t, file=treeOut)

      if(0){ #this could be missing repeats in bed1 that are present in bed2, so must do another way below
        hc <- as.hclust(hm$rowDendrogram)
        t <- as.phylo(hc)
        treeOut <- paste(outDir, "/", paste("bed1", "regions", statsTypes[l], repCat, sep = "-"), ".newick", sep = "")
        write.tree(phy=t, file=treeOut)
      }
      #has to remove rows with 0 variance to get it to actually work....
      c_DF1 <- df1[apply(df1, 1, var) > 0,]
      tc_DF1 <- t(c_DF1)

      row_hclust <- tc_DF1 %>% as.matrix %>% t %>% vegdist("jaccard",binary=TRUE) %>% hclust(method="ward.D2")
      row.phylo <- as.phylo(row_hclust)
      treeOut <- paste(outDir, "/", paste("bed1", "regions", statsTypes[l], repCat, sep = "-"), ".newick", sep = "")
      write.tree(phy=row.phylo, file=treeOut)

      #pvclust
      c_DF1 <- df1[apply(df1, 1, var) > 0,]
      tc_DF1 <- t(c_DF1)
      rc_DF1 <- tc_DF1[apply(tc_DF1, 1, var) > 0,]
      trc_DF1 <- t(rc_DF1)
      if(goPVClust){
        options(max.print=1000000)
        rc_result <- pvclust(rc_DF1, nboot=nbs, parallel=as.integer(pvParaCore))
        out <- paste(outDir, "/", paste("bed1-pvclust-reg", statsTypes[l], repCat, sep = "-"), ".pdf", sep = "")
        pdf(out, height = (1*regTreeSize + 10), width = ((log(regTreeSize) * 200) + 10))
        plot(rc_result)
        pvrect(rc_result, alpha=0.95)
        dev.off()

        trc_result <- pvclust(trc_DF1, nboot=nbs, parallel=as.integer(pvParaCore))
        out <- paste(outDir, "/", paste("bed1-pvclust-rep", statsTypes[l], repCat, sep = "-"), ".pdf", sep = "")
        pdf(out, height = (1*repTreeSize + 10), width = ((log(repTreeSize) * 200) + 10))
        plot(trc_result)
        pvrect(trc_result, alpha=0.95)
        dev.off()

        conOut <- paste(outDir, "/", paste("bed1-pvclust-reg", statsTypes[l], repCat, sep = "-"), ".txt", sep = "")
        con<-textConnection(conOut, "w")
        sink(con)
        pvpick(rc_result)
        print(rc_result)
        sink()
        close(con)

        conOut <- paste(outDir, "/", paste("bed1-pvclust-rep", statsTypes[l], repCat, sep = "-"), ".txt", sep = "")
        con<-textConnection(conOut, "w")
        sink(con)
        pvpick(trc_result)
        print(trc_result)
        sink()
        close(con)
        options(max.print=1000)
      }

      #heatmap
      for( cs in 1:3){
        out <- paste(outDir, "/", paste("bed1", statsTypes[l], repCat, "coloring", coloring[cs], sep = "-"), ".pdf", sep = "")
        pdf(out, height = (35 + regTreeSize), width = (85 + repTreeSize))
        par(oma=c(20, 10, 10, 70)) # s, w, n, e
        heatmap.2(dm, col = heat.colors(256), scale = coloring[cs], cexRow = 3, cexCol = 3, lhei = c(3,regTreeSize), lwid = c(2,repTreeSize))
        #print(hm2graphic)
        dev.off()

        #bmpFN <- paste(outDir, "/", paste("bed1", statsTypes[l], repCat, "coloring", coloring[cs], sep = "-"), ".bmp", sep = "")
        #bmp(filename = bmpFN)
        #par(oma=c(20,3,3,60))
        #print(hm2graphic)
        #dev.off()

        out <- paste(outDir, "/", paste("bed1", statsTypes[l], repCat, "noTrace", "coloring", coloring[cs], sep = "-"), ".pdf", sep = "")
        pdf(out, height = (35+regTreeSize), width = (85+repTreeSize))
        par(oma=c(20, 10, 10, 70)) # s, w, n, e
        heatmap.2(dm, col = heat.colors(256), scale = coloring[cs], cexRow = 3, cexCol = 3, lhei = c(3,regTreeSize), lwid = c(2,repTreeSize), trace="none")
        #print(hm2graphic)
        dev.off()
        #bmpFN <- paste(outDir, "/", paste("bed1", statsTypes[l], repCat, "noTrace", "coloring", coloring[cs], sep = "-"), ".bmp", sep = "")
        #bmp(filename = bmpFN)
        #par(oma=c(20,3,3,60))
        #print(hm2graphic)
        #dev.off()
      }
    }

    #mclust for bed1 regions and repeats
    if(goClassify){
      #for the regions
      mc_df1 <- Mclust(df1, G=c(seq(1,9),seq(10,49,5), seq(50,149,10), seq(150,500,50)))
      mc_df1_class <- as.data.frame(mc_df1$classification)
      df1_temp <- data.frame(do.call("rbind", strsplit(row.names(mc_df1_class), "-+|:", perl = TRUE)))
      df1_temp["cluster"] <- mc_df1_class
      out <- paste(outDir, "/", paste("bed1-mclust-reg", statsTypes[l], repCat, sep = "-"), ".bed", sep = "")

      write.table(df1_temp, file=out, col.names=F, row.names=F, sep="\t", quote=FALSE)

      out <- paste(outDir, "/", paste("bed1-mclust-reg", statsTypes[l], repCat, sep = "-"), ".pdf", sep = "")
      pdf(out)
      plot(mc_df1, what="BIC")
      dev.off()

      conOut <- paste(outDir, "/", paste("bed1-mclust-reg-summary", statsTypes[l], repCat, sep = "-"), ".txt", sep = "")
      con <-textConnection(conOut, "w")
      con<-file(out, "w")
      sink(con)
      summary(mc_df1)
      sink()
      close(con)

      index <- runif(length(mc_df1$classification), min=1, max=100)
      indexTF <- index < 75
      df1_mc_train <- mc_df1$classification[indexTF]
      df1_mc_test <- mc_df1$classification[!indexTF]

      #f1_knn <- knn(df1[indexTF,], df1[!indexTF,], df1_mc_train)

      bestKreg <- 1
      numRight <- 0

  pseudoNoise <- runif(dim(df1)[1] * dim(df1)[2], min = 0, max = 0.0000000001)
  df1PN <- df1 + pseudoNoise
      for(tryK in seq(1,50, by = 2)){
        df1_knn <- knn(df1PN[indexTF,], df1PN[!indexTF,], df1_mc_train, k=tryK)
        er <- sum(df1_mc_test == df1_knn)
        if(er > numRight){
          if(tryK <= bestKreg){
             numRight <- er
             bestKreg <- tryK
          }
        }
      }

      conOut <- paste(outDir, "/", paste("bed1-knn-reg-summary", statsTypes[l], repCat, sep = "-"), ".txt", sep = "")
      con<-textConnection(conOut, "w")
      sink(con)
      print("Tried from 1 to 50 k nearest neighbors")
      print(c("best knn", bestKreg))
      print(c("Number correct", numRight))
      print(c("Number tested", sum(!indexTF)))
      print(c("Error rate", ((sum(!indexTF) - numRight) / sum(!indexTF))))
      sink()
      close(con)


      #mclust for bed1 repeats
      t_df1 <- t(df1)
      mc_tdf1 <- Mclust(t_df1, G=c(seq(1,19),seq(20,49,5), seq(50,100,10)))
      mc_tdf1_class <- as.data.frame(mc_tdf1$classification)
      out <- paste(outDir, "/", paste("bed1-mclust-rep", statsTypes[l], repCat, sep = "-"), ".tsv", sep = "")
      tdf1_temp <- cbind(row.names(mc_tdf1_class), mc_tdf1_class)

      write.table(tdf1_temp, file=out, col.names=F, row.names=F, sep="\t", quote=FALSE)

      out <- paste(outDir, "/", paste("bed1-mclust-rep", statsTypes[l], repCat, sep = "-"), ".pdf", sep = "")
      pdf(out)
      plot(mc_tdf1, what="BIC")
      dev.off()

      conOut <- paste(outDir, "/", paste("bed1-mclust-rep-summary", statsTypes[l], repCat, sep = "-"), ".txt", sep = "")
      con<-textConnection(conOut, "w")      
      sink(con)
      summary(mc_tdf1)
      sink()
      close(con)

      index <- runif(length(mc_tdf1$classification), min=1, max=100)
      indexTF <- index < 75
      tdf1_mc_train <- mc_tdf1$classification[indexTF]
      tdf1_mc_test <- mc_tdf1$classification[!indexTF]

      #f1_knn <- knn(df1[indexTF,], df1[!indexTF,], df1_mc_train)

      bestKrep <- 1
      numRight <- 0
  t_df1 <- t(df1PN)
      for(tryK in seq(1,50, by=2)){
        tdf1_knn <- knn(t_df1[indexTF,], t_df1[!indexTF,], tdf1_mc_train, k=tryK)
        er <- sum(tdf1_mc_test == tdf1_knn)
        if(er > numRight){
          if(tryK <= bestKrep){
            numRight <- er
            bestKrep <- tryK
          }
        }
      }

      conOut <- paste(outDir, "/", paste("bed1-knn-rep-summary", statsTypes[l], repCat, sep = "-"), ".txt", sep = "")
      con<-textConnection(conOut, "w")
      sink(con)
      print("Tried from 1 to 50 k nearest neighbors")
      print(c("best knn", bestKrep))
      print(c("Number correct", numRight))
      print(c("Number tested", sum(!indexTF)))
      print(c("Error rate", ((sum(!indexTF) - numRight) / sum(!indexTF))))
      sink()
      close(con)
    }
  }

  #bed2
  #bed2
  #bed2
  #bed2

  dm <- b2Mat[,1:ncol(b2Mat)-1]
  dm <- data.matrix(dm)

  df2 <- as.data.frame(b2Mat[,1:ncol(b2Mat)-1])

  repTreeSize <- as.integer(ncol(dm)/3)
  regTreeSize <- as.integer(nrow(dm)/2.5)

  #bulk of bed2 and combined processing
  if(nrow(dm) > 1 && ncol(dm) > 1){
    if(goCluster){
      pdf(NULL)
      hm <- heatmap.2(dm)
      hc <- as.hclust(hm$colDendrogram)
      t <- as.phylo(hc)
      treeOut <- paste(outDir, "/", paste("bed2", "repeats", statsTypes[l], repCat, sep = "-"), ".newick", sep = "")
      write.tree(phy=t, file=treeOut)
      dev.off()

      if(0){
        hc <- as.hclust(hm$rowDendrogram)
        t <- as.phylo(hc)
        treeOut <- paste(outDir, "/", paste("bed2", "regions", statsTypes[l], repCat, sep = "-"), ".newick", sep = "")
        write.tree(phy=t, file=treeOut)
      }

      c_DF2 <- df2[apply(df2, 1, var) > 0,]
      tc_DF2 <- t(c_DF2)

      row_hclust <- tc_DF2 %>% as.matrix %>% t %>% vegdist("jaccard",binary=TRUE) %>% hclust(method="ward.D2")
      row.phylo <- as.phylo(row_hclust)
      treeOut <- paste(outDir, "/", paste("bed2", "regions", statsTypes[l], repCat, sep = "-"), ".newick", sep = "")
      write.tree(phy=row.phylo, file=treeOut)

      #pvclust
      c_DF2 <- df2[apply(df2, 1, var) > 0,]
      tc_DF2 <- t(c_DF2)
      rc_DF2 <- tc_DF2[apply(tc_DF2, 1, var) > 0,]
      trc_DF2 <- t(rc_DF2)
      if(goPVClust){
        options(max.print=1000000)
        rc_result <- pvclust(rc_DF2, nboot=nbs, parallel=as.integer(pvParaCore))
        out <- paste(outDir, "/", paste("bed2-pvclust-reg", statsTypes[l], repCat, sep = "-"), ".pdf", sep = "")
        pdf(out, height = (1*regTreeSize + 10), width = ((log(regTreeSize) * 200) + 10))
        plot(rc_result)
        pvrect(rc_result, alpha=0.95)
        dev.off()

        trc_result <- pvclust(trc_DF2, nboot=nbs, parallel=as.integer(pvParaCore))
        out <- paste(outDir, "/", paste("bed2-pvclust-rep", statsTypes[l], repCat, sep = "-"), ".pdf", sep = "")
        pdf(out, height = (1*repTreeSize + 10), width = ((log(repTreeSize) * 200) + 10))
        plot(trc_result)
        pvrect(trc_result, alpha=0.95)
        dev.off()

        conOut <- paste(outDir, "/", paste("bed2-pvclust-reg", statsTypes[l], repCat, sep = "-"), ".txt", sep = "")
        con<-textConnection(conOut, "w")
        sink(con)
        pvpick(rc_result)
        print(rc_result)
        sink()
        close(con)

        outOut <- paste(outDir, "/", paste("bed2-pvclust-rep", statsTypes[l], repCat, sep = "-"), ".txt", sep = "")
        con<-textConnection(conOut, "w")
        sink(con)
        pvpick(trc_result)
        print(trc_result)
        sink()
        close(con)
        options(max.print=1000)
      }


      #heatmap stuff
      for( cs in 1:3){
        out <- paste(outDir, "/", paste("bed2", statsTypes[l], repCat, "coloring", coloring[cs], sep = "-"), ".pdf", sep = "")
        pdf(out, height = (35+regTreeSize), width = (85+repTreeSize))
        par(oma=c(20, 10, 10, 70)) # s, w, n, e
        heatmap.2(dm, col = heat.colors(256), scale = coloring[cs], cexRow = 3, cexCol = 3, lhei = c(3,regTreeSize), lwid = c(2,repTreeSize))
        #print(hm2graphic)
        dev.off()

        #bmpFN <- paste(outDir, "/", paste("bed1", statsTypes[l], repCat, "coloring", coloring[cs], sep = "-"), ".bmp", sep = "")
        #bmp(filename = bmpFN)
        #par(oma=c(20,3,3,60))
        #print(hm2graphic)
        #dev.off()

        out <- paste(outDir, "/", paste("bed2", statsTypes[l], repCat, "noTrace", "coloring", coloring[cs], sep = "-"), ".pdf", sep = "")
        pdf(out, height = (35+regTreeSize), width = (85+repTreeSize))
        par(oma=c(20, 10, 10, 70)) # s, w, n, e
        heatmap.2(dm, col = heat.colors(256), scale = coloring[cs], cexRow = 3, cexCol = 3, lhei = c(3,regTreeSize), lwid = c(2,repTreeSize), trace="none")
        #print(hm2graphic)
        dev.off()
        #bmpFN <- paste(outDir, "/", paste("bed2", statsTypes[l], repCat, "noTrace", "coloring", coloring[cs], sep = "-"), ".bmp", sep = "")
        #bmp(filename = bmpFN)
        #par(oma=c(20,3,3,60))
        #print(hm2graphic)
        #dev.off()
      }
    }

    if(goClassify){
      pseudoNoise <- runif(dim(df2)[1] * dim(df2)[2], min = 0, max = 0.0000000001)
      df2PN <- df2 + pseudoNoise
  
      df2_knn <- knn(df1PN, df2PN, mc_df1$classification, k=bestKreg)
      df2_knnClassified <- data.frame(do.call("rbind", strsplit(row.names(df2), "-+|:", perl = TRUE)))
      df2_knnClassified <- data.frame(df2_knnClassified, df2_knn)

      out <- paste(outDir, "/", paste("bed2-reg-knn-classified", statsTypes[l], repCat, sep = "-"), ".bed", sep = "")
      write.table(df2_knnClassified, file=out, col.names=F, row.names=F, sep="\t", quote=FALSE)
    }

    #combined
    #combined
    #combined
    #combined

    dm <- data.matrix(combinedMat)
    dm <- dm[,-ncol(dm)]
    sideColors <- c(rep("blue", nrow(b1Mat)), rep("red", nrow(b2Mat)))

    dfc <- as.data.frame(combinedMat)
    dfc <- dfc[,1:(dim(dfc)[2]-1)]

    repTreeSize <- as.integer(ncol(dm)/3)
    regTreeSize <- as.integer(nrow(dm)/2.5)
    if(goCluster){
      pdf(NULL)
      hm <- heatmap.2(dm)
      hc <- as.hclust(hm$colDendrogram)
      t <- as.phylo(hc)
      treeOut <- paste(outDir, "/", paste("combined", "repeats", statsTypes[l], repCat, sep = "-"), ".newick", sep = "")
      write.tree(phy=t, file=treeOut)
      dev.off()

      if(0){
        hc <- as.hclust(hm$rowDendrogram)
        t <- as.phylo(hc)
        treeOut <- paste(outDir, "/", paste("combined", "regions", statsTypes[l], repCat, sep = "-"), ".newick", sep = "")
        write.tree(phy=t, file=treeOut)
      }

      c_DFC <- dfc[apply(dfc, 1, var) > 0,]
      tc_DFC <- t(c_DFC)

      row_hclust <- tc_DFC %>% as.matrix %>% t %>% vegdist("jaccard",binary=TRUE) %>% hclust(method="ward.D2")
      row.phylo <- as.phylo(row_hclust)
      treeOut <- paste(outDir, "/", paste("combined", "regions", statsTypes[l], repCat, sep = "-"), ".newick", sep = "")
      write.tree(phy=row.phylo, file=treeOut)

      #pvclust
      c_DFC <- dfc[apply(dfc, 1, var) > 0,]
      tc_DFC <- t(c_DFC)
      rc_DFC <- tc_DFC[apply(tc_DFC, 1, var) > 0,]
      trc_DFC <- t(rc_DFC)
      if(goPVClust){
        options(max.print=1000000)
        rc_result <- pvclust(rc_DFC, nboot=1000, parallel=as.integer(pvParaCore))
        out <- paste(outDir, "/", paste("combined-pvclust-reg", statsTypes[l], repCat, sep = "-"), ".pdf", sep = "")
        pdf(out, height = (1*regTreeSize + 10), width = ((log(regTreeSize) * 200) + 10))
        plot(rc_result)
        pvrect(rc_result, alpha=0.95)
        dev.off()

        trc_result <- pvclust(trc_DFC, nboot=1000, parallel=as.integer(pvParaCore))
        out <- paste(outDir, "/", paste("combined-pvclust-rep", statsTypes[l], repCat, sep = "-"), ".pdf", sep = "")
        pdf(out, height = (1*repTreeSize + 10), width = ((log(repTreeSize) * 200) + 10))
        plot(trc_result)
        pvrect(trc_result, alpha=0.95)
        dev.off()

        conOut <- paste(outDir, "/", paste("combined-pvclust-reg", statsTypes[l], repCat, sep = "-"), ".txt", sep = "")
        con<-textConnection(conOut, "w")
        sink(con)
        pvpick(rc_result)
        print(rc_result)
        sink()
        close(con)

        conOut <- paste(outDir, "/", paste("combined-pvclust-rep", statsTypes[l], repCat, sep = "-"), ".txt", sep = "")
        con<-textConnection(conOut, "w")
        sink(con)
        pvpick(trc_result)
        print(trc_result)
        sink()
        close(con)
        options(max.print=1000)
      }

      #heatmap
      for( cs in 1:3){
        out <- paste(outDir, "/", paste("combined", statsTypes[l], repCat, "coloring", coloring[cs], sep = "-"), ".pdf", sep = "")
        pdf(out, height = (35+regTreeSize), width = (85+repTreeSize))
        par(oma=c(20, 10, 10, 70)) # s, w, n, e
        heatmap.2(dm, col = heat.colors(256), scale = coloring[cs], cexRow = 3, cexCol = 3, RowSideColors = sideColors, lhei = c(3,regTreeSize), lwid = c(2,repTreeSize))
        #print(hm2graphic)
        dev.off()

        #bmpFN <- paste(outDir, "/", paste("bed1", statsTypes[l], repCat, "coloring", coloring[cs], sep = "-"), ".bmp", sep = "")
        #bmp(filename = bmpFN)
        #par(oma=c(20,3,3,60))
        #print(hm2graphic)
        #dev.off()

        out <- paste(outDir, "/", paste("combined", statsTypes[l], repCat, "noTrace", "coloring", coloring[cs], sep = "-"), ".pdf", sep = "")
        pdf(out, height = (35+regTreeSize), width = (85+repTreeSize))
        par(oma=c(20, 10, 10, 70)) # s, w, n, e
        heatmap.2(dm, col = heat.colors(256), scale = coloring[cs], cexRow = 3, cexCol = 3, RowSideColors = sideColors, lhei = c(3,regTreeSize), lwid = c(2,repTreeSize), trace="none")
        #print(hm2graphic)
        dev.off()
        #bmpFN <- paste(outDir, "/", paste("bed2", statsTypes[l], repCat, "noTrace", "coloring", coloring[cs], sep = "-"), ".bmp", sep = "")
        #bmp(filename = bmpFN)
        #par(oma=c(20,3,3,60))
        #print(hm2graphic)
        #dev.off()
      }
    }

    if(goClassify){
      #mclust - predict number of clusters in another way
      fit <- Mclust(dfc, G=c(seq(1,9),seq(10,49,5), seq(50,149,10), seq(150,500,50)))
      conOut <- paste(outDir, "/", paste("combined-mclust-dfc", statsTypes[l], repCat, sep = "-"), ".txt", sep = "")
      con<-textConnection(conOut, "w")
      sink(con)
      summary(fit)
      sink()
      close(con)

      out <- paste(outDir, "/", paste("combined-mclust-reg", statsTypes[l], repCat, sep = "-"), ".txt", sep = "")
      write.csv(fit$classification, file=out)

      #kmeans
      pseudoNoise <- runif(dim(dfc)[1] * dim(dfc)[2], min = 0, max = 0.0000000001)
      dfcPN <- dfc + pseudoNoise
      km_dm <- kmeans(dfcPN,2)
      givenClust <- as.integer(gsub("Bed", "", combinedMat[,dim(combinedMat)[2]]))
      names(givenClust) <- rownames(combinedMat)
      mskiDist <- dist(dfcPN, method = "minkowski")
      cStats <- cluster.stats(mskiDist,km_dm$cluster, givenClust)

      conOut <- paste(outDir, "/", paste("combined-kmeans2-vs-supplied_minkowski-distance", statsTypes[l], repCat, sep = "-"), ".txt", sep = "")
      con<-textConnection(conOut, "w")
      sink(con)
      cStats
      sink()
      close(con)
    }
  }
}


#do the individual tests for each level of repeat
if(goIndividual){
  for ( i in 1:3){
    repType <- catLabels[i]

    #create plots for individual types of repeats
    #remove empty elements
    reps <- t(as.data.frame(fread(repeatFileList[i], header=FALSE, sep = "\t")))
    repLen <- length(reps)
    repsName <- reps[1:repLen]
    repeatNames <- reps
    #create CDF plot for each repeat
    dir.create("individualRepeatTests")
    #create CDF and scatter plots for divergence 3, del 4, ins 5, len 6
    #this is the portion of the script that takes the longest and does not need to be rerun after you
    #identify the repeats that are of interest
    statsTypes <- c("percent-divergence", "percent-deleted", "percent-inserted", "length", "lenRegRepNorm", "lenRegDNANorm", "lenRegRepNorm", "lenRegDNANorm", "percent-identity")
    statsNames <- c("percent divergence", "percent deleted", "percent inserted", "length", "fraction of interval repetitive DNA", "fraction of interval DNA", "fraction of interval repetitive DNA", "fraction of interval DNA", "percent identity")
    offset <- i + 1

    b1Build <- b1FullData[offset]
    b2Build <- b2FullData[offset]

    for( st in 1:6){#this is for the columns in the data files
      dir.create(paste("individualRepeatTests/", catLabels[i], "-", statsTypes[st], sep = ""))
      dir.create(paste("individualRepeatTests/", catLabels[i], "-",  statsTypes[st],"/violin", sep = ""))
      dir.create(paste("individualRepeatTests/", catLabels[i], "-",  statsTypes[st],"/cdf", sep = ""))
   }
   
   #see bed1Stats.tsv for organization of b1FullData
   for(st in 1:3){
      k <- st + 7
      #get data for each stats type, pay attention to the offset
      b1Build <- cbind(b1Build,b1FullData[k])
      b2Build <- cbind(b2Build,b2FullData[k])
      cdfStatName <- paste(catLabels[i], statsTypes[st], sep = "-")
      humanName <- statsNames[st]
      repLevel <- catLabels[i]
    }
    #for length
    b1Build <- cbind(b1Build,b1FullData[11])
    b2Build <- cbind(b2Build,b2FullData[11])
    #for normalized to rep length
    b1Build <- cbind(b1Build,(b1FullData[11]/b1FullData[18]))
    b2Build <- cbind(b2Build,(b2FullData[11]/b2FullData[18]))
    #for normalized to dna length
    b1Build <- cbind(b1Build,(b1FullData[11]/b1FullData[6]))
    b2Build <- cbind(b2Build,(b2FullData[11]/b2FullData[6]))
	#for normalized to total Rep length
	b1TotalRepLength <- sum(b1FullData[6])
	b2TotalRepLength <- sum(b2FullData[6])

    b1Build <- cbind(b1Build,(b1FullData[11]/b1TotalRepLength))
    b2Build <- cbind(b2Build,(b2FullData[11]/b2TotalRepLength))
    #for normalized to total DNA length
    b1TotalDNALength <- sum(b1FullData[18])
	b2TotalDNALength <- sum(b2FullData[18])
	
    b1Build <- cbind(b1Build,(b1FullData[11]/b1TotalDNALength))
    b2Build <- cbind(b2Build,(b2FullData[11]/b2TotalDNALength))
    
    #for total percent difference
    pidentity1 <- 100 - rowSums(b1FullData[,c(8,9,10)])
    pidentity2 <- 100 - rowSums(b2FullData[,c(8,9,10)])
    b1Build <- cbind(b1Build,pidentity1)
    b2Build <- cbind(b2Build,pidentity2)

    b1 <- b1Build
    b2 <- b2Build
    
    allCombos <- expand.grid(repeatNames, 1:11)
    
    #testSet <- sort(sample(1:dim(allCombos)[1], 20))
    #system.time(statsData <- mclapply(testSet, individualApply, allCombos = allCombos, statsTypes = statsTypes, statsNames = statsNames, repLevel = repLevel, b1=b1, b2=b2, mc.cores=nCores))
    system.time(statsData <- mclapply(1:dim(allCombos)[1], individualApply, allCombos = allCombos, statsTypes = statsTypes, statsNames = statsNames, repLevel = repLevel, b1=b1, b2=b2, mc.cores=nCores))

    resultColNames <- c("repeat", "Bed1 counts", "Bed2 counts", "KS boot p-value Lower", "KS p-value Greater", "AD Cont. p-value", "AD Disc. p-value") 
    statsData <- matrix(unlist(statsData), ncol = 8, byrow=TRUE)
    for( st in 1:6){#this is for the columns in the data files
      statsData[which(statsData[,1] == st),2:8]
      resultStats <- as.data.frame(statsData[which(statsData[,1] == st),2:8])
      colnames(resultStats) <- resultColNames
      write.csv(lapply(resultStats,as.character), file = paste("individualRepeatTests/summaryStats", "--", catLabels[i], "--", statsTypes[st], "--", repType, ".csv", sep = ""))
    }
  }
}


#do the clustering/classification for every type of stat on the defined intervals
if(goCluster || goClassify || goPies){
  for ( i in 1:3 ){
    repType <- catLabels[i]

    #create plots for individual types of repeats
    #remove empty elements
    reps <- t(as.data.frame(fread(repeatFileList[i], header=FALSE, sep = "\t")))
    repLen <- length(reps)
    repsName <- reps[1:repLen]

    statsTypes <- c("avg_divergence", "avg_deletion", "avg_insertion", "avg_length", "sum_rep_length", "sum_rep_count", "sum_rep_len_norm2totalRepLen", "sum_rep_len_norm2totalDNALen", "sum_rep_count_norm2totalRepLen", "sum_rep_count_norm2totalDNALen")
    humanStats <- c("avgerage divergence", "avgerage deletion", "avgerage insertion", "average length", "sum_rep_length", "sum_rep_count", "sum_rep_len_norm2totalRepLen", "sum_rep_len_norm2totalDNALen", "sum_rep_count_norm2totalRepLen", "sum_rep_count_norm2totalDNALen")
    dir.create("IntervalAnalyses")
    repCat <- catLabels[i];
    if(goCluster || goClassify){
      system.time(applyResults <- mclapply(1:length(statsTypes), clAnalysesApply, statsTypes = statsTypes, repCat = repCat, reps = reps, bed1Reg = bed1Reg, bed2Reg = bed2Reg, repLen = repLen, mc.cores=nCores))
    }

    #pie charts
    if(goPies){ #fast enough it doesn't need to be multithreaded
      #create the Pie charts and scatter plots
      statsTypes <- c("countRaw", "countNorm2Rep", "countNorm2DNA", "bpRaw", "bpNorm2Rep", "bpNorm2DNA")
      dir.create("Pies")

      for( l in 1:6){
        #load in and add labels to data
        pieFN <- paste(paste("pie", catLabels[i], statsTypes[l], sep = "-"), ".tsv", sep = "")
        pieData <- as.data.frame(fread(pieFN, header=FALSE, sep = "\t"))
        repName <- reps
        if(l == 4 || l == 6){
          repName <- c(repName, "Not_repetitive")
        }

        names(pieData) <- c("Repeats", "Bed1", "Bed2", "Merged")
        outDir <- "Pies"

        pdfFN <- paste(outDir, "/", paste(statsTypes[l], catLabels[i], "Bed1", sep = "-"), ".pdf", sep = "")
        pdf(file = pdfFN)
        p <- ggplot(pieData, aes(x=1, y=Bed1, fill=Repeats)) +
          ggtitle(paste(catLabels[i], statsTypes[l], sep = " ")) +
          coord_polar(theta="y")
        p <- p + geom_bar(stat="identity", color="black") +
          guides(fill=guide_legend(override.aes=list(colour=NA)))

        p <- p + theme(axis.ticks=element_blank(), axis.title=element_blank(), axis.text.y=element_blank())
        y.breaks <- cumsum(pieData$Bed1) - pieData$Bed1/2

        p <- p + theme(axis.text.x=element_text(color="black")) +
          scale_y_continuous(breaks=y.breaks, labels=pieData$Repeats) +
          theme(panel.background = element_rect(fill = "white"))
        #p <- p + theme(plot.margin=unit(c(3,3,3,3), "cm"))

        print(p)
        dev.off()

        pdfFN <- paste(outDir, "/", paste(statsTypes[l], catLabels[i], "Bed2", sep = "-"), ".pdf", sep = "")
        pdf(file = pdfFN)
        p <- ggplot(pieData, aes(x=1, y=Bed2, fill=Repeats)) +
          ggtitle(paste(catLabels[i], statsTypes[l], sep = " ")) +
          coord_polar(theta="y")
        p <- p + geom_bar(stat="identity", color="black") +
          guides(fill=guide_legend(override.aes=list(colour=NA)))

        p <- p + theme(axis.ticks=element_blank(), axis.title=element_blank(), axis.text.y=element_blank())
        y.breaks <- cumsum(pieData$Bed2) - pieData$Bed2/2

        p <- p + theme(axis.text.x=element_text(color="black")) +
          scale_y_continuous(breaks=y.breaks, labels=pieData$Repeats) +
          theme(panel.background = element_rect(fill = "white"))

        print(p)
        dev.off()

        if(l == 4 || l == 6){
          pieData <- pieData[-nrow(pieData),]
          pdfFN <- paste(outDir, "/", paste(statsTypes[l], catLabels[i], "Bed1-noDNA", sep = "-"), ".pdf", sep = "")
          pdf(file = pdfFN)
          p <- ggplot(pieData, aes(x=1, y=Bed1, fill=Repeats)) +
            ggtitle(paste(catLabels[i], statsTypes[l], sep = " ")) +
            coord_polar(theta="y")
          p <- p + geom_bar(stat="identity", color="black") +
            guides(fill=guide_legend(override.aes=list(colour=NA)))

          p <- p + theme(axis.ticks=element_blank(), axis.title=element_blank(), axis.text.y=element_blank())
          y.breaks <- cumsum(pieData$Bed1) - pieData$Bed1/2

          p <- p + theme(axis.text.x=element_text(color="black")) +
            scale_y_continuous(breaks=y.breaks, labels=pieData$Repeats) +
            theme(panel.background = element_rect(fill = "white"))

          print(p)
          dev.off()

          pdfFN <- paste(outDir, "/", paste(statsTypes[l], catLabels[i], "Bed2-noDNA", sep = "-"), ".pdf", sep = "")
          pdf(file = pdfFN)
          p <- ggplot(pieData, aes(x=1, y=Bed2, fill=Repeats)) +
            ggtitle(paste(catLabels[i], statsTypes[l], sep = " ")) +
            coord_polar(theta="y")
          p <- p + geom_bar(stat="identity", color="black") +
            guides(fill=guide_legend(override.aes=list(colour=NA)))

          p <- p + theme(axis.ticks=element_blank(), axis.title=element_blank(), axis.text.y=element_blank())
          y.breaks <- cumsum(pieData$Bed2) - pieData$Bed2/2

          p <- p + theme(axis.text.x=element_text(color="black")) +
            scale_y_continuous(breaks=y.breaks, labels=pieData$Repeats) +
            theme(panel.background = element_rect(fill = "white"))

          print(p)
          dev.off()
        }
      }
    }
  }
}

#savehistory(file="workspace-commands.Rhistory")
#save(file="workspace-all.RData")
save.image(file=workspace.RData")
