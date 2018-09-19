library(gplots)
library(reshape2)
library(scater)
library(stringr)
library(sqldf)
library(DESeq)
library(DESeq2)
library(Rtsne)
library(limma)
library(BiocParallel)

nworkers <- 16

options(MulticoreParam=quote(MulticoreParam(workers=nworkers)))

source("common_functions.R")


######################################################################
############# Read plate layout ######################################
######################################################################

## Note: we only provide the libraries that worked on arrayexpress. Use derived count tables
## instead of this mess

#Read layout, plate 1
oelay_f1 <- read.table("walkup final oe/cond_treatment_1.csv",   #index D
                       stringsAsFactors = FALSE, sep = ",")
oelay_m1 <- read.table("walkup final oe/cond_mouse_1.csv", 
                       stringsAsFactors = FALSE, sep = ",")
oelay1 <- NULL
for(i in 1:12){
  for(j in 1:8){
    oelay1 <- rbind(oelay1,
                    data.frame(
                      well=sprintf("%s%s",c("A","B","C","D","E","F","G","H")[j],i), 
                      mouse=oelay_m1[j,i], 
                      treatment=oelay_f1[j,i], 
                      stringsAsFactors = FALSE))
  }
}

#Read layout, plate 2
oelay_f2 <- read.table("walkup final oe/cond_treatment_2.csv",   #index B
                       stringsAsFactors = FALSE, sep = ",")
oelay_m2 <- read.table("walkup final oe/cond_mouse_2.csv", 
                       stringsAsFactors = FALSE, sep = ",")
oelay2 <- NULL
for(i in 1:12){
  for(j in 1:8){
    oelay2 <- rbind(oelay2,
                    data.frame(
                      well=sprintf("%s%s",c("A","B","C","D","E","F","G","H")[j],i), 
                      mouse=oelay_m2[j,i], 
                      treatment=oelay_f2[j,i], 
                      stringsAsFactors = FALSE))
  }
}



wellname <- matrix(nrow=8,ncol=12)
for(i in 1:12){
  for(j in 1:8){
    wellname[j,i]<-sprintf("%s%s",c("A","B","C","D","E","F","G","H")[j],i)
  }
}

####### Build up index numbers for rnaseq plates

rnaseqPlateIndex <- read.table("walkup final oe/rnaseqplateindex.csv",stringsAsFactors = FALSE, sep = ",")
rnaseqPlateIndex <- as.matrix(rnaseqPlateIndex)

mapLibindexWell <- rbind(
  data.frame(well=wellname[1:96], libraryid=rnaseqPlateIndex[1:96], plate=1),
  data.frame(well=wellname[1:96], libraryid=rnaseqPlateIndex[1:96]+96, plate=2),
  data.frame(well=wellname[1:96], libraryid=rnaseqPlateIndex[1:96]+96+96, plate=3),
  data.frame(well=wellname[1:96], libraryid=rnaseqPlateIndex[1:96]+96+96+96, plate=4))
mapLibindexWell$libraryname <- sprintf("samp_%s",mapLibindexWell$libraryid)

mapLibindexWell

oelay1$plate <- 4  #index set D
oelay2$plate <- 2  #index set B

oelay <- rbind(oelay1,oelay2)

cellcond <- sqldf("select * from mapLibindexWell natural left outer join oelay")
cellcond <- cellcond[order(cellcond$libraryid),]
cellcond


######################################################################
############# Read count data and make condition matrix ##############
######################################################################

dat1 <- read.csv("walkup final oe/counttable1.csv",stringsAsFactors = FALSE, sep = ",", row.names = 1)
dat2 <- read.csv("walkup final oe/counttable2.csv",stringsAsFactors = FALSE, sep = ",", row.names = 1)

dat <- dat1+dat2  #read.csv("walkup final oe/counttable.csv",stringsAsFactors = FALSE, sep = ",", row.names = 1)
rownames(dat) <- str_split_fixed(rownames(dat),"\\.",2)[,1]

mapGeneidTransid <- read.csv("walkup final oe/map_ensg_ensmust.csv", stringsAsFactors = FALSE)
rownames(mapGeneidTransid) <- mapGeneidTransid$transid
dat_geneid <- mapGeneidTransid[rownames(dat),]$ensemble_gene_id

#Sum on per-gene level. Drop unknown transcripts
dat_geneid[is.na(dat_geneid)] <- "NA"
sum(dat_geneid=="NA")
dat_pergene <- rowsum(dat, group=dat_geneid)
dat_pergene <- dat_pergene[rownames(dat_pergene)!="NA",]
dat <- dat_pergene


######################################################################
############# Read distribution across plate #########################
######################################################################


##Read distribution across plate. index B and D relevant
cellcond$total_reads <- apply(dat,2,sum)
outmat <- matrix(nrow=8, ncol=12)
for(i in 1:12){
  for(j in 1:8){
    outmat[j,i] <- cellcond$total_reads[cellcond$plate==4 & cellcond$well==wellname[j,i]]
  }
}
round(outmat)



######################################################################
############# QC #####################################################
######################################################################

#only include wells that may have some content
keep <- !is.na(cellcond$treatment) & cellcond$treatment!="none"
dat <- dat[,keep]
cellcond <- cellcond[keep,]

#Number of exonic reads
gene_count <- as.integer(colSums(dat[grep("ENSMUSG",rownames(dat)),]))
cellcond$gene_count <- gene_count 
#Number of detected genes
detected_genes <- as.integer(colSums(dat>0))
#Proportion of mitochondrial reads ------ this one is a bit suspicious...
mt_counts <- as.integer(colSums(na.rm = TRUE,dat[mt_genes$ensembl_gene_id,]))
mt_prop <- mt_counts / gene_count  



findbadcells <- function(forcells=1:ncol(dat),mt_cutoff=0.1){
  # Using scater to assign outliers
  libsize.drop <- isOutlier(gene_count[forcells], nmads=2, type="lower",log=TRUE)
  feature.drop <- isOutlier(detected_genes[forcells], nmads=2, type="lower", log=TRUE)
  mito.drop <- mt_prop[forcells]>mt_cutoff
  
  todrop <- libsize.drop | feature.drop | mito.drop
  print(dim(todrop))
  totc<-length(todrop)
  print(as.matrix(list(
    total=totc,
    kept=totc-length(which(todrop)),
    drop.by.size=length(which(libsize.drop)),
    drop.by.features=length(which(feature.drop)),
    drop.by.mito=length(which(mito.drop))
  )))
  
  todroptot <- rep(FALSE, ncol(dat))
  todroptot[(1:ncol(dat))[forcells][todrop]] <- TRUE
  
  par(mfrow=c(1,2))
  hist(log10(gene_count[forcells]), xlab="Mapped reads (log-scale)", main="",breaks=50, col="grey80", ylab="Number of cells")
  v<-max(gene_count[forcells][libsize.drop])
  abline(v=log10(v), col="red", lty=2, lwd=1.5)
  hist(detected_genes[forcells], xlab="Number of detected genes", main="",breaks=50, col="grey80", ylab="Number of cells")
  v<-max(detected_genes[forcells][feature.drop])
  abline(v=v, col="red", lty=2, lwd=1.5)
  # hist(mt_prop[forcells], xlab="mtDNA %", main="",breaks=50, col="grey80", ylab="Number of cells")
  # abline(v=mt_cutoff, col="red", lty=2, lwd=1.5)
  #dev.off()
  
  todroptot
}

#Remove crappy libraries
pdf("walkup final oe/qc_oescreen.pdf")
toremove <- findbadcells()
toremove <- detected_genes<1000  #override. only removes 3 samples. one is zfp397
dev.off()
keep <- !toremove   #wants to drop 11 in total. We can maybe save 4 of them if we want

cellcond[toremove,]


dat <- dat[,keep]
cellcond <- cellcond[keep,]

dim(dat)
dim(cellcond)


######################################################################
############# Write out for upload ###################################
######################################################################

genes_ae <- c("Lrrc40dec","Lrrc40","Pparg","Scara3","bfp","Ccdc134","Bhlhe40")
keep <- cellcond$treatment %in% genes_ae & cellcond$mouse %in% c("1","2","3")

write.csv(dat[,keep],"out/oe_cnt.csv")
write.csv(cellcond[keep,],"out/oe_cellcond.csv", row.names = FALSE)

######################################################################
############# Prepare for analysis ###################################
######################################################################

ncnt <- normalizeDeseqByAll(dat)


######################################################################
############# Cluster conditions #####################################
######################################################################


#Batch-correct normalized counts
batchcorr <- limma::removeBatchEffect(
  log2(1+ncnt),
  batch = cellcond$mouse,
  covariates = cellcond$depth)


keepgenes <- apply(ncnt,1,mean)>0.2
sum(keepgenes)



set.seed(0)
ncntred <- batchcorr[keepgenes,]
d <- stats::dist(t((ncntred)))
rtsne_out <- Rtsne(d,is_distance = TRUE, perplexity=30, verbose = TRUE, max_iter=5000,dims = 2)
keep <- cellcond$treatment!="FOOOOOOOO"

#keep <- cellcond$treatment!="FOOOOOOOO"
keep <- -grep("oe",cellcond$treatment)

pdf("walkup final oe/tsne.pdf")
plot(rtsne_out$Y[keep,], pch=16, main='',xlab="",ylab="",cex=0)  
text(rtsne_out$Y[keep,1], rtsne_out$Y[keep,2],
     labels=cellcond$treatment[keep],cex=1.0)#,   col=tabrora$colplate[keep])
dev.off()




######################################################################
############# Linear model to compare conditions ##################### 
######################################################################

cellcond$treatment <- factor(cellcond$treatment)
cellcond$mouse <- factor(cellcond$mouse)

########## Fit linear model
datred <- dat[apply(ncnt,1,mean)>0.3,]
nrow(datred)
save(datred, cellcond, file="oe_in.RData")    ### for running on cluster
load("oe_in.RData")
dds <- DESeqDataSetFromMatrix(countData = round(datred),
                              colData = cellcond,
                              design = ~treatment + mouse)
dds <- DESeq(dds,
             parallel = TRUE, 
             BPPARAM = MulticoreParam(workers=nworkers))


########## Test all conditions vs WT
deseq_one_vs_all <- function(basecond){
  xbp_de <- list()
  for(curcond in unique(cellcond$treatment)){
    if(curcond!=basecond){
      print(curcond)
      v <- results(dds,
                   contrast=c("treatment",curcond, basecond),    #numerator / denominator
                   parallel = TRUE, 
                   BPPARAM = MulticoreParam(workers=nworkers)) 
      v$ensembl_gene_id <- rownames(v)
      v <- as.data.frame(v[order(v$padj),])
      xbp_de[[curcond]] <- v
    }
  }
  xbp_de
}


########## Assemble list of tests into a matrix
de_list2matrix <- function(xbp_de){
  v<-xbp_de[[names(xbp_de)[1]]]
  xbp_fc    <- matrix(ncol=length(xbp_de), nrow=nrow(v))
  xbp_p     <- matrix(ncol=length(xbp_de), nrow=nrow(v))
  xbp_fcse  <- matrix(ncol=length(xbp_de), nrow=nrow(v))
  rownames(xbp_fc)   <- v$ensembl_gene_id
  rownames(xbp_p)    <- v$ensembl_gene_id
  rownames(xbp_fcse) <- v$ensembl_gene_id
  colnames(xbp_fc)   <- names(xbp_de)
  colnames(xbp_p)    <- names(xbp_de)
  colnames(xbp_fcse) <- names(xbp_de)
  #store lfSE too?
  for(i in 1:length(xbp_de)){
    curcond <- names(xbp_de)[i]
    v <- xbp_de[[curcond]]
    xbp_fc  [v$ensembl_gene_id,i] <- v$log2FoldChange
    xbp_p   [v$ensembl_gene_id,i] <- v$padj
    xbp_fcse[v$ensembl_gene_id,i] <- v$lfcSE
  }
  list(
    fc=as.data.frame(xbp_fc), 
    padj=as.data.frame(xbp_p), 
    fcse=as.data.frame(xbp_fcse))  
}
oe_de <- de_list2matrix(deseq_one_vs_all("bfp"))
save.image(file="oe_out.RData") 


########## Test all conditions vs mean
deseq_one_vs_mean <- function(){
  xbp_de <- list()
  for(curcond in unique(cellcond$treatment)){
    print(curcond)
    
    thecont <- rep(0, length(resultsNames(dds)))
    indgene <- resultsNames(dds)==sprintf("treatment%s",curcond)  
    thecont[indgene] <- 1
    v <- results(dds,
                 contrast=thecont,    
                 parallel = TRUE, 
                 BPPARAM = MulticoreParam(workers=nworkers)) 
    v$ensembl_gene_id <- rownames(v)
    v <- as.data.frame(v[order(v$padj),])
    xbp_de[[curcond]] <- v
  }
  xbp_de
}
oe_de_mean <- de_list2matrix(deseq_one_vs_mean())
save.image(file="oe_out_vsmean.RData") 










##################################################################################################################
##################################################################################################################
###################################### Interpret results #########################################################
##################################################################################################################
##################################################################################################################





#######################
## for MA plots
plot_sym_ma <- function(thegene, filteroef=TRUE, dosort=TRUE, filterChIP=NULL){
  #thegene <- "Gata3"
  if(filteroef){
    oefs <- grep("oe",colnames(oe_de_sym$padj))
    keepg <- apply(oe_de_sym$padj[,oefs]<1e-5,1,function(x) sum(x, na.rm = TRUE))<2
    print(sum(!keepg))
    #keepg <- apply(oe_de_sym$padj<1e-5,1,function(x) sum(x, na.rm = TRUE))<4
    oe_de_sym$fc   <- oe_de_sym$fc[keepg,]
    oe_de_sym$padj <- oe_de_sym$padj[keepg,]
  }
  
  tpm <- data.frame(
    ensembl_gene_id=names(tcmouse$max_mtpm),
    max_tpm=round(tcmouse$max_mtpm))
  
  d <- data.frame(
    ensembl_gene_id=oe_de_sym$padj$ensembl_gene_id,
    mgi_symbol=oe_de_sym$padj$mgi_symbol, 
    fc=oe_de_sym$fc[,thegene],
    padj=oe_de_sym$padj[,thegene])  
  
  if(!is.null(filterChIP)){
    d <- d[d$ensembl_gene_id %in% filterChIP$Nearest.Ensembl,]
  }
  
  plot(d$fc, log10(d$padj),cex=0)  #  , ylim=c(-3,0)
  text(d$fc, log10(d$padj), d$mgi_symbol)
  
  d$fc <- round(d$fc, 2)
  d <- merge(d, tpm)
  if(dosort){
    d<-d[order(d$padj),]
  }
  d
}


#######################
## Check if the OEF are overrepresented. Then do not trust the target
checkOErep <- function(genesym){
  theid <- ensconvert$ensembl_gene_id[ensconvert$mgi_symbol==genesym]
  oefs <- grep("oe",colnames(oe_de_sym$padj))
  wilcox.test(
    as.double(oe_de$padj[theid, oefs]),
    as.double(oe_de$padj[theid, -oefs])
  )$p.value
}

#######################
## Check effect on particular targets
listTargetSortedByPadj <- function(thegene){
  sort(oe_de$padj[ensconvert$ensembl_gene_id[ensconvert$mgi_symbol==thegene],])  
}







##############################
## Load the DE data again
load("walkup final oe/oe_out_ccdc134.RData") 
load("walkup final oe/oe_out.RData") 

# load("walkup final oe/oe_out_vsmean.RData") 
# oe_de <- oe_de_mean

write.csv(oe_de$fc[,c("Lrrc40","Pparg","Scara3","Ccdc134","Bhlhe40")],"out/de_fc.csv")
write.csv(oe_de$padj[,c("Lrrc40","Pparg","Scara3","Ccdc134","Bhlhe40")],"out/de_padj.csv")


oe_de$fc[,normalizesym(genes_ae)]
normalizesym(genes_ae)



#load("walkup final oe/oe_out_withmouse_notred.RData") 
#load("walkup final oe/oe_out_nomouse.RData") 


##############################
## Integrate Xbp1 data
load("../oe_xbp1.RData")                      #careful with up or down. check.
oe_de$fc$Xbp1   <- xbp_de$fc[rownames(oe_de$fc),]$oe
oe_de$padj$Xbp1 <- xbp_de$padj[rownames(oe_de$padj),]$oe


##############################
## Filter out genes which might be due to a faulty control
oe_de_sym <- oe_de
oe_de_sym$padj$ensembl_gene_id <- rownames(oe_de_sym$padj)
oe_de_sym$fc$ensembl_gene_id <- rownames(oe_de_sym$fc)
oe_de_sym$padj <- merge(oe_de_sym$padj, ensconvert)
oe_de_sym$fc <- merge(oe_de_sym$fc, ensconvert)



