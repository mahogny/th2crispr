# 
# allde.mouse <- smerge(data.frame(
#   ens_mouse = ensconvert$ensembl_gene_id,
#   me=ensconvert$ensembl_gene_id %in% getdefromtable.mouse(tcmouse$de_early),
#   ml=ensconvert$ensembl_gene_id %in% getdefromtable.mouse(tcmouse$de_late), stringsAsFactors = FALSE), 
#   ortho_mouse_human_unique, all.x = TRUE)


#which are the 1-1 orthologs?

#### Convert human table into mouse equivalents. Only keep overlap
getmousehumanTPM <- function(){
  tpm1 <- tcmouse$mtpm
  tpm2 <- tchuman$mtpm
  thet <- ortho_mouse_human_unique[
    ortho_mouse_human_unique$ens_mouse %in% rownames(tpm1) &
      ortho_mouse_human_unique$ens_human %in% rownames(tpm2),]
  rownames(thet) <- thet$ens_mouse
  
  tpm1 <- tpm1[rownames(tpm1) %in% thet$ens_mouse,]
  tpm2 <- tpm2[rownames(tpm2) %in% thet$ens_human,]
  tpm2 <- tpm2[thet[rownames(tpm1),]$ens_human,]
  rownames(tpm2) <- rownames(tpm1)
  
  
  colnames(tpm1) <- str_replace(colnames(tpm1),"rep","")
  colnames(tpm1) <- str_replace(colnames(tpm1),"_salmon","")

  colnames(tpm2) <- str_replace(colnames(tpm2),"rep","")
  colnames(tpm2) <- str_replace(colnames(tpm2),"_salmon","")
  
  list(mouse=tpm1,human=tpm2)  
}
otpm <- getmousehumanTPM()

#Only keep genes that are sufficiently expressed
mintpm <- 0
keep <- apply(otpm$mouse,1,mean)>mintpm & apply(otpm$human,1,mean)>mintpm
otpm$mouse <- otpm$mouse[keep,]
otpm$human <- otpm$human[keep,]

# colnames(otpm$mouse)
# colnames(otpm$human)
# colnames(tcmouse$mtpm)
# cor(log(1+otpm$mouse[,1]), log(1+otpm$human[,1]))

ta <- apply(otpm$mouse,1,mean)
otpm$mouse <- otpm$mouse/ta
ta <- apply(otpm$human,1,mean)
otpm$human <- otpm$human/ta

plot(log(1+otpm$mouse[,"Th2_72h_1"]), log(1+otpm$human[,"Th2_72h_1"]))
plot(log(1+otpm$mouse[,"Th2_72h_1"]), log(1+otpm$human[,"Naive_2"]))


calccortimesmatrix <- function(){
  #otpm$mouse <- log()

  timesn <- c("Naive",sprintf("%sh",c("05","1","2","4","6","12","24","48","72")))
  times0 <- c("Naive",sprintf("Th0_%sh",c("05","1","2","4","6","12","24","48","72")))
  times2 <- c("Naive",sprintf("Th2_%sh",c("05","1","2","4","6","12","24","48","72")))
  allcor <- matrix(NA, ncol=length(times0),nrow=length(times0))
  for(i in 1:length(times0)){
    for(j in 1:length(times0)){
      thecol1 <- c(sprintf("%s_%s",times0[i],1:3), sprintf("%s_%s",times2[i],1:3))
      thecol2 <- c(sprintf("%s_%s",times0[j],1:3), sprintf("%s_%s",times2[j],1:3))
      cols1 <- which(colnames(otpm$mouse) %in% thecol1)
      cols2 <- which(colnames(otpm$human) %in% thecol2)

      cs <- c()
      for(k in cols1)
        for(m in cols2){
#          diff <- otpm$mouse[,k]/sum(otpm$mouse[,k])-otpm$human[,m]/sum(otpm$human[,m])
#          diff <- otpm$mouse[,k]-otpm$human[,m]
          #cs <- c(cs, diff*diff)
          cs <- c(cs, cor(otpm$mouse[,k], otpm$human[,m],method="spearman"))
        }
      allcor[i,j] <- mean(cs)#*1e8
    }
  }
  rownames(allcor) <- timesn
  colnames(allcor) <- timesn
  allcor
}
allcor <- calccortimesmatrix()


#only keep the expressed genes. above some TPM. or not really


#log-space comparison


#reply natalia

#Correlate 


p <- aracne_mousetc[aracne_mousetc$Target==toensid("Ccdc134"),]



############################################################################
############ Compare using DE, for every time point ########################
############################################################################


#### Load data
dat1 <- read.csv("out_tc/human/all_samples_human_genes_Salmon_counts.txt",sep="\t",row.names = "gene")
dat2 <- read.csv("out_tc/mouse/all_samples_mouse_genes_Salmon_counts.txt",sep="\t",row.names = "gene")
colnames(dat1) <- str_replace(colnames(dat1), "Thp","Naive")

#### Do DE for every timepoint
listM <- list()
listH <- list()
timesn <- c("Naive",sprintf("%sh",c("05","1","2","4","6","12","24","48","72")))
for(i in 2:length(timesn)){ #Cannot compare 1 obviously
  #Mouse
  x <- dat1[,grep(timesn[i],colnames(dat1))]
  ddsFullCountTable <- DESeqDataSetFromMatrix(
    countData = round(x),
    colData = data.frame(row.names = colnames(x), is2 = factor(colnames(x) %!in% colnames(x)[grep("Th2",colnames(x))])),
    design = ~ is2)
  ddsFullCountTable <- DESeq(ddsFullCountTable)
  res1 <- results(ddsFullCountTable)
  listM[[i]] <- res1  

  #Human
  x <- dat2[,grep(timesn[i],colnames(dat2))]
  ddsFullCountTable <- DESeqDataSetFromMatrix(
    countData = round(x),
    colData = data.frame(row.names = colnames(x), is2 = factor(colnames(x) %!in% colnames(x)[grep("Th2",colnames(x))])),
    design = ~ is2)
  ddsFullCountTable <- DESeq(ddsFullCountTable)
  res2 <- results(ddsFullCountTable)
  listH[[i]] <- res2
}

#### Plot comparison

#pdf("out_tc/mousehuman_fc_alltimes.pdf")
png("out_tc/mousehuman_fc_alltimes.png", w=2000,h=2000)
#ton <- 5 
ton <- length(timesn)
par(mfrow=c(ton-1,ton-1), mar=c(1,1,1,1))
#par(mfrow=c(length(timesn)-1,length(timesn)-1))
for(i in 2:ton){
  for(j in 2:ton){
    v1 <- listH[[i]]   #mixup, fix later TODO
    v2 <- listM[[j]]
    vc <- ortho_mouse_human_unique[
      ortho_mouse_human_unique$ens_mouse %in% rownames(v1) &
        ortho_mouse_human_unique$ens_human %in% rownames(v2),]
    rownames(vc) <- vc$ens_mouse
    
    #v1 <- v2[rownames(v1) %in% vc$ens_mouse,]
    #v2 <- v2[rownames(v2) %in% vc$ens_human,]
    
    v1 <- v1[vc$ens_mouse,]
    v2 <- v2[vc$ens_human,]
    
    v1 <- v1$log2FoldChange
    v2 <- v2$log2FoldChange
    v1[is.na(v1)] <- 0
    v2[is.na(v2)] <- 0
    
    plot(v1, v2,cex=0.1)
  }
}
dev.off()



allcorMH <- matrix(nrow = ton-1, ncol = ton-1)
for(i in 2:ton){
  for(j in 2:ton){
    v1 <- listH[[i]]   #mixup, fix later TODO
    v2 <- listM[[j]]
    vc <- ortho_mouse_human_unique[
      ortho_mouse_human_unique$ens_mouse %in% rownames(v1) &
        ortho_mouse_human_unique$ens_human %in% rownames(v2),]
    rownames(vc) <- vc$ens_mouse
    
    #v1 <- v2[rownames(v1) %in% vc$ens_mouse,]
    #v2 <- v2[rownames(v2) %in% vc$ens_human,]
    
    v1 <- v1[vc$ens_mouse,]
    v2 <- v2[vc$ens_human,]
    
    v1 <- v1$log2FoldChange
    v2 <- v2$log2FoldChange
    v1[is.na(v1)] <- 0
    v2[is.na(v2)] <- 0
    
    allcorMH[i-1,j-1] <- cor(v1,v2,method = "spearman")
    #plot(v1$log2FoldChange, v2$log2FoldChange,cex=0.1)
  }
}
write.table(allcorMH, file = "out_tc/mousehuman_corfc.csv", row.names = FALSE, sep="\t")
## row is human, col is mouse

#plot.new()
#plot.new()
# i<-2
# j<-3






