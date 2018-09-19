#############################################################################################
###                                                                                       ###
###             Part of the paper ...                                                     ###
###             Author: Johan Henriksson (mahogny@areta.org)                              ###
###                                                                                       ###
###             This code preprocesses the KO RNAseq data: fold change etc                ###
###                                                                                       ###
#############################################################################################

library(stringr)
library(DESeq2)
library(DESeq)
library(biomaRt)
library(scater)
library(sqldf)
library(gplots)
library(RColorBrewer)
library(RUVSeq)
library(Rtsne)
library(gProfileR)


######################################################################
### Load data ########################################################
######################################################################

listnongenes <- c(
  "no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique","grna-numread",
  "__no_feature","__ambiguous","__too_low_aQual","__not_aligned","__alignment_not_unique")

# read normal gene counts
readcounttable20161209one <- function(f){
  dat <- read.csv(f,stringsAsFactors=FALSE)
  rownames(dat) <- dat[,1]
  dat <- dat[,-1]
  rownames(dat)[which(rownames(dat)=="CAS9")] <- "cas9"
  rownames(dat)[which(rownames(dat)=="grna-Ptma-F-1")] <- "grna-ptma1"
  rownames(dat)[which(rownames(dat)=="grna-Ptma-F-2")] <- "grna-ptma2"
  rownames(dat)[which(rownames(dat)=="grna-Ptma-F-3")] <- "grna-ptma3"
  rownames(dat)[which(rownames(dat)=="grna-Ptma-F-4")] <- "grna-ptma4"
  rownames(dat)[which(rownames(dat)=="grna-Ptma-F-5")] <- "grna-ptma5"
  dat
}
readcounttable20170211 <- function(){
  dat <- readcounttable20161209one("ko/raw20170211.csv")
  colnames(dat) <- sprintf("%s",1:ncol(dat))
  dat
}
dat <- readcounttable20170211()

### Regularized log
log1 <- function(x) as.double(log(1+x))


### Compute normalized counts
mynormalize <- function(dat){
  dat2 <- dat
  dat666 <- dat[grep("ENS",rownames(dat)),] 
  sf <- estimateSizeFactorsForMatrix(dat666)
  sfnew <- rep(1,ncol(dat))
  sfnew[(1:ncol(dat))]<-sf
  for(i in 1:ncol(dat)){
    dat2[,i] <- dat2[,i]/sfnew[i]
  }
  dat2
}
ncount <- mynormalize(dat)
ncount <- as.matrix(ncount)




### Construct plate design
getplatedesign20170211 <- function(){
  
  numplates <- ncol(dat)/96 
  
  thecols <- rep(1:12,8)
  thecols <- rep(thecols,numplates)
  therows <- as.double(sapply(1:8,function(i) rep(i,12)))
  therows <- rep(therows,numplates)
  theplate <- as.double(sapply(1:numplates,function(i) rep(i,96)))
  isgood <- rep(TRUE,ncol(dat))
  
  genelay <- as.matrix(read.csv("ko/layouttcellfull.csv",stringsAsFactors=FALSE))
  #  genelay <- cbind(genelay,genelay,genelay,rep("f",8),rep("f",8),rep("f",8),rep("f",8)) #not really true
  ko <- rep("",length(thecols))
  mouse <- theplate
  for(i in 1:length(thecols))  {
    if(theplate[i] %in% c(1,2,3)){
      ko[i] <- genelay[therows[i],thecols[i]]
    }
  }
  
  cellcondition <- data.frame(
    col = thecols,
    row = therows,
    plate = theplate,
    isgood=isgood,
    ko=ko,
    mouse=mouse,
    stringsAsFactors=FALSE
  )
  rownames(cellcondition) <- colnames(dat)
  cellcondition$probe <- sprintf("%s_%s_%s",cellcondition$ko,cellcondition$col, cellcondition$row)
  cellcondition$probe[cellcondition$ko=="Thy1control"] <- "Thy1control"
  cellcondition
}
cellcondition <- getplatedesign20170211()
#write.csv(cellcondition, "ko/cellcondition.csv", row.names = FALSE)


######################################################################
### Study the impact of BFP ##########################################
######################################################################

### this is the key thing
badbfp <- c(9,  17,  28,  56,  64,  68,  79,  96, 105, 117, 152, 160, 175, 201, 209, 220, 256, 271, 288)

plot(sort(ncount["BFP",]))



###################################################################################################
########### Formulate linear model ################################################################  
###################################################################################################

### Subset data
usewell <- cellcondition$isgood &    !((1:nrow(cellcondition)) %in% badbfp)   #barely makes a difference to remove these
cellcondition_red <- cellcondition[usewell,]
dat_red <- as.matrix(dat[,usewell])
levbfp <- log(1+ncount["BFP",usewell])  #need to use the normalized one! log because the other genes are log too?
levdepth <- log(1+apply(dat_red,2,sum))

#Make KO part of design matrix.
#note: Thy1control should not get a variable
if(FALSE){  #KO-level
} else { #Probe-level
  un <- setdiff(unique(cellcondition_red$ko), "Thy1control")
  cellcondition_red$ko=cellcondition_red$probe
}
un <- setdiff(unique(cellcondition_red$ko), "Thy1control")
ddsColKo <- matrix(0,nrow=nrow(cellcondition_red), ncol=length(un))
for(i in 1:length(un)){
  w <- cellcondition_red$ko==un[i]
  ddsColKo[w,i] <- 1
  #ddsColKo[w,i] <- levbfp[w] #this is the best but by a very fine margin    #used this on probes, good result
  #  ddsColKo[w,i] <- exp(levbfp[w]) #this is the worst
}
colnames(ddsColKo) <- sprintf("ko%s",un)
ddsColKo <- ddsColKo[,apply(ddsColKo,2,sum)>0]  #Remove non-existent probes/genes

makebfpn <- function(m){
  x <- levbfp
  x[cellcondition_red$mouse==m] <- max(x[cellcondition_red$mouse==m]) - x[cellcondition_red$mouse==m]
  x[cellcondition_red$mouse!=m] <- 0
  x
}

makedepthn <- function(m){
  x <- levdepth
  x[cellcondition_red$mouse==m] <- max(x[cellcondition_red$mouse==m]) - x[cellcondition_red$mouse==m]
  x[cellcondition_red$mouse!=m] <- 0
  x
}

##### Assembly the total condition matrix
ddsCol <- cbind(data.frame(
  #bfp=levbfp, #d[usewell],
  # bfp1=makebfpn(1),
  # bfp2=makebfpn(2),
  # bfp3=makebfpn(3),
  #  bfp=max(levbfp) - levbfp, #bfp=0 means best condition
  #  depth=log(1+apply(dat_red,2,sum)),
  # d1=makedepthn(1),
  # d2=makedepthn(2),
  # d3=makedepthn(3),
  ko=relevel(factor(cellcondition_red$ko),"Thy1control"),
  mouse=factor(cellcondition_red$mouse)
),ddsColKo)
#further kick out genes to speed things up. In particular, do NOT keep BFP in here because it is in the cond matrix
dat_red <- as.matrix(dat_red[rownames(dat_red) %in% expressedGenes10Id & rownames(dat_red)!="BFP",])

## Write data to disk to solve on cluster
outdir <- "ko/test9"
write.table(dat_red, sprintf("%s/dat.csv", outdir), row.names = TRUE)
write.table(ddsCol,  sprintf("%s/col.csv", outdir), row.names = TRUE)

#Build formula: linear combination of all columns
pastevecsep <- function(sdata, sep) do.call(paste, c(as.list(sdata), sep=sep))
f <- as.formula(paste("~", pastevecsep(colnames(ddsCol)," + ")))





############################################################################################################################
##########  Analyze DESeq2 result ##########################################################################################
############################################################################################################################

######
load("ko/with_bfpmouse/de1a.Rdata")
load("ko/with_bfpmouse/de2a.Rdata")
load("ko/degrna_withmouse.Rdata")
load("ko/degene_withmouse.Rdata")

#Get smallest P-value from each condition
v<-round(sort(log(apply(newkoPadj,2,function(x) min(x,na.rm = TRUE)))))
v<-v[order(names(v))]
v

# stopgosym(
#   togenesym2(na.exclude(rownames(newkoPadj)[newkoPadj[,"koThy1_12_1"]<1e-3])),
#   togenesym2(rownames(newkoPadj))
# )

#Rules:
## Kick out one replicate if difference >5
## Remove anything == 0
#or: check that p-value better in combined. otherwise remove the bad one



if(FALSE){
  toclust <- newkoFC
  toclust <- toclust[apply(!is.na(toclust),1,all),]
  d <- stats::dist(t((toclust)))
  set.seed(0) 
  rtsne_out <- Rtsne(d,is_distance = TRUE, perplexity=20, verbose = TRUE, max_iter=5000,dims = 2)
  plot(rtsne_out$Y, pch=16, main='',xlab="",ylab="",cex=0)  # col=rcol[toycellcondition$isgood], 
  text(rtsne_out$Y[,1], rtsne_out$Y[,2], 
       labels=colnames(newkoFC),cex=0.5) 
  heatmap(cor(toclust))
}





# calcreprodfor <- function(gene){
#   toclustcor[grep(gene,colnames(toclustcor)),grep(gene,colnames(toclustcor))]
# }
# calcreprodfor("koXbp1") #something
# calcreprodfor("koGata3") #poor
# calcreprodfor("koIl4") #something
# calcreprodfor("koIl13")
# calcreprodfor("koThy1") #great
# calcreprodfor("koZc3h12a") #still great


sort(repcor)

##########################################################################################
##### Assess reproducibility with and without BFP term ###################################
##########################################################################################

##################################
# Function: calculate reproducibility between gRNAs
calcrepcor <- function(){
  toclust <- newkoFC
  toclust <- toclust[apply(!is.na(toclust),1,all),]
  toclustcor <- cor(toclust)
  repcor <- c()
  un <- totest <- unique(unlist(lapply(colnames(toclustcor),function(s) unlist(str_split(s,"_"))[[1]])))
  for(i in 1:length(un)){
    tt <- grep(un[i], colnames(toclustcor))
    if(length(tt)>1){
      repcor <- c(repcor, toclustcor[tt[1],tt[2]])
    } else {
      repcor <- c(repcor, NA)
    }
  }
  names(repcor) <- un
  plot(density(bw=0.01,na.exclude(repcor)),xlim=c(-1,1))
  lines(density(bw=0.01,as.double(toclustcor)),col="blue")
}


load("ko/with_bfpmouse/de1a.Rdata")
load("ko/with_bfpmouse/de2a.Rdata")
repcorMouseBfp <- calcrepcor()

load("ko/degrna_withmouse.Rdata")
#load("ko/degene_withmouse.Rdata")
repcorMouse <- calcrepcor()
vGene <- v

if(FALSE){
  plot(repcorMouse, repcorMouseBfp,cex=0)
  text(repcorMouse, repcorMouseBfp, labels = names(repcorMouse))
  lines(c(-1,1),c(-1,1))
}

if(FALSE){
  ### Correlation vs pvalues in the replicates
  compCorPval <- na.omit(smerge(data.frame(n=names(repcorMouse), repcor=repcorMouse),
                        data.frame(n=names(vGene), v=vGene)))
  plot(compCorPval$repcor, compCorPval$v, cex=0)
  text(compCorPval$repcor, compCorPval$v, labels=compCorPval$n)
}

#How many DE genes in each condition
sort(apply(newkoPadj<1e-2,2,function(x) sum(x,na.rm = TRUE)))
#How many genes with high FC in each condition
sort(apply(abs(newkoFC)>log(1.5,2),2,function(x) sum(x,na.rm = TRUE)))

compnumde <- function(){
  x<-data.frame(
    byp=apply(newkoPadj<1e-2,2,function(x) sum(x,na.rm = TRUE)),
    byfc=apply(abs(newkoFC)>log(1.5,2),2,function(x) sum(x,na.rm = TRUE)))
  qtextscatter(log(1+x$byp), log(1+x$byfc), rownames(x))
  x
}


##########################################################################################
#### Detect poor gRNAs using pvalues gene vs pvalues one grna ############################
##########################################################################################

load("ko/degrna_withmouse.Rdata")
vguide<-round(sort(log(apply(newkoPadj,2,function(x) min(x,na.rm = TRUE)))))
load("ko/degene_withmouse.Rdata")
vgene<-round(sort(log(apply(newkoPadj,2,function(x) min(x,na.rm = TRUE)))))

w<-0
for(i in 1:length(vgene)){
  pguides <- vguide[grep(names(vgene)[i],names(vguide))]
  if(any(pguides<vgene[i])){
    print(names(vgene)[i])
    print(pguides)
    print(as.double(vgene[i]))
    w <- w+1
  }
}
w 
length(vgene)


##########################################################################################
### Remove bad gRNAs and show new result #################################################
##########################################################################################

## Based on above pvalues, we decided to remove these:
listBadGrna <- c(
  #"koBcl11b_1_6","koBcl11b_1_7","koBcl11b_1_5",  #Bcl11b looks fine without removing these
  "koFam26f_5_1",
  "koPgk1_10_1",
  "koNrd1_9_3",
  "koHerc6_6_2",
  "koIl13_7_1",
  "koGata3_5_8",
  "koIl12a_6_8")

load("ko/degrna_withmouse.Rdata")
grna_newkoFC <- newkoFC
grna_newkoPadj <- newkoPadj
grna_newkoFC <- grna_newkoFC[,!(colnames(grna_newkoFC) %in% listBadGrna)]
grna_newkoPadj <- grna_newkoPadj[,!(colnames(grna_newkoPadj) %in% listBadGrna)]
load("ko/degene_withmouse.Rdata")
compnumde()   #before fixing

### Update pvalue/FC for the genes where there is now only one replicate left (=use pval/fc of that replicate)
recalcgene <- unique(unlist(lapply(listBadGrna,function(s) unlist(str_split(s,"_"))[[1]])))
for(i in 1:length(recalcgene)){
  j <- which(colnames(newkoPadj)==recalcgene[i])
  k <- grep(recalcgene[i],colnames(grna_newkoFC))
  newkoPadj[,j] <- grna_newkoPadj[,k]
  newkoFC[,j] <- grna_newkoFC[,k]
}
compnumde()   #after fixing

### Wherever there is an NA-value, set as "non-hit"
for(i in 1:ncol(newkoPadj)){
  newkoPadj[is.na(newkoPadj[,i]),i] <- 1
  newkoFC[is.na(newkoFC[,i]),i] <- 0
}

### Store updated KO results
write.table(newkoPadj,"out_ko/newko_pval.csv",quote = FALSE, sep="\t")
write.table(newkoFC,"out_ko/newko_fc.csv",quote = FALSE, sep="\t")

write.table(grna_newkoPadj,"out_ko/grna_newko_pval.csv",quote = FALSE, sep="\t")
write.table(grna_newkoFC,"out_ko/grna_newko_fc.csv",quote = FALSE, sep="\t")



############################################################################################################################
################# Calculate GO for each KO #################################################################################
############################################################################################################################



calcnewkogo <- function(n=400){
  listfg <- list()
  listbg <- list()
  for(i in 1:ncol(newkoFC)){
    print(i)
    
    #sort by FC, then by pval. so in case of tie, ordered by fc
    v <- na.omit(data.frame(n=rownames(newkoFC), pval=newkoPadj[,i], fc=newkoFC[,i], stringsAsFactors = FALSE))
    v <- v[order(abs(v$fc),decreasing = TRUE),]    
    v <- v[order(v$pval, decreasing = FALSE),]
    
    listfg[[i]] <- togenesym2(v$n[1:n])
    listbg[[i]] <- togenesym2(rownames(dat_red))
  }
  
  newkoGO <- getgomatrix(listfg, listbg, listnames=colnames(newkoFC), cutoff=1e-1, useTopgo=TRUE)
  colnames(newkoGO) <- colnames(newkoFC)
  newkoGO
}
   
newkoGO300 <- calcnewkogo(300)
newkoGO <- calcnewkogo(400)
newkoGO1000 <- calcnewkogo(1000)
#write.csv(newkoGO,"out_ko/newko_go.csv")

#write.csv(newkoGO300,"out_ko/newko_go300.csv")
#write.csv(newkoGO,"out_ko/newko_go.csv")

plot.new()
pdf("out_ko/newko_go.pdf",width = 10, height = 10)
plotgomatrix(newkoGO)#, keepimmuno = TRUE)
plotgomatrix(newkoGO1000, keepimmuno = FALSE, cutoff = -7)#, keepimmuno = TRUE)
dev.off()
plot.new()
plotgomatrix(newkoGO, keepimmuno = TRUE)
  #dev.off()

listweakko <- c("koIl4","koIl12a","koXbp1","koPgk1","koCxcr7",
                "koBcl11b","koF2rl1","koCrls1","koEtv2","koCd180",
                "koRyr1","koOas1f","koNhedc2","koLag3","koCcdc134",
                "koGata3","koMed8","koZc3h12a","koMapkapk3","koCcna2")
plotgomatrix(newkoGO[,colnames(newkoGO) %in% listweakko])





######################################################################
### KO similarity to Th_x - classification ###########################
######################################################################

##Normalization of strength between the KO since we have an FC cutoff
newkoFCnorm <- newkoFC
for(i in 1:ncol(newkoFC)){
  newkoFCnorm[,i] <- newkoFCnorm[,i]/sd(newkoFCnorm[,i])
}


### by using DE fold changes
theall<-c()
allsd<-c()
numgenetested <- c()
for(j in 1:ncol(thep)){
  print(sprintf("Thx: %s",j))
  #Only consider DE genes according to ThExpress
  totest <- !is.na(thep[,j]) & thep[,j]<1e-2
  totest <- intersect(rownames(thefc)[totest],rownames(newkoPadj))

  cmin<-function(x) sum(na.omit(x)==-1)
  cplus<-function(x) sum(na.omit(x)==1)
  kov <- c()
  for(i in 1:ncol(outpadj)){
    #Further only consider DE genes according to the KO
    usefc <- newkoFCnorm
    totest3 <- intersect(
      totest,
      rownames(usefc)[abs(usefc[,i])>log(2,2)])     #inc/dec 10%  is good... here more

    numgenetested <- c(numgenetested,length(totest3))
    s1<-sign(thefc[totest3,j])
    s2<-sign(usefc[totest3,i])
    a<-sum(s1*s2)
    
    # unfortunately this is not a great bootstrap distribution so not attempting any fancy corrections
    # boot <- rep(NA,1000)
    # for(k in 1:1000){
    #   s2b<-sample(s2,length(s2))   #sign(usefc[totest3,i])
    #   boot[i] <- sum(s1*s2b)
    # }
    
    kov <- c(kov, a/length(totest)) #-mean(boot)
  }
  if(j==1)  
    theall<-kov
  else
    theall<-rbind(theall,kov)
}
colnames(theall) <- colnames(newkoPadj)
rownames(theall) <- sprintf("Th2 vs %s",colnames(thep))
hist(numgenetested,breaks=20)

thenorm <- theall
for(i in 1:nrow(thenorm)){
  thenorm[i,] <- thenorm[i,]/sd(thenorm[i,])
}

colnames(theall) <- sapply(colnames(theall),function(s) str_sub(s,3))
colnames(thenorm) <- colnames(theall)

pdf("out_ko/newko_thtype_norm_new.pdf",width = 10, height = 2.5)
thecol <- colorRampPalette(c("#dd0152", "white","#2719dd"))(n = 299) #red to blue
heatmap.2(thenorm[-1,],
          trace="none",
          #density.info="none", 
          dendrogram = "none",
          cexRow = 0.8,
          col=thecol,
          Rowv=FALSE)#,
#          hclustfun = function(x) hclust(x,method="ward.D2"))
dev.off()

# pdf("out_ko/newde_thtype_notnorm.pdf",width = 10, height = 2.5)
# thecol <- colorRampPalette(c("#8e0152", "white","#276419"))(n = 299) #red to green
# heatmap.2(theall,
#           trace="none",
#           #density.info="none", 
#           dendrogram = "none",
#           cexRow = 0.8,
#           col=thecol,
#           Rowv=FALSE)#,
# #          hclustfun = function(x) hclust(x,method="ward.D2"))
# dev.off()



# http://www.sciencedirect.com/science/article/pii/S1074761300806251
# this is maybe THE IL13 paper, from andrew.
# Il13 is fairly high in Th1 too. 
# Il13 persists as anti-Th2 no matter the parameters
# Il13 needed in vitro, but not in vivo. note that here it is KOed during mid-development


# BIG QUESTION: is IL4 and Il13 regulated by the same genes? this we should write about
# onepermutation <- function(x) sample(x,length(x))
plot(log10(1+sgenescorer2_matrix[,"Il13"]),log10(1+sgenescorer2_matrix[,"Il4"]))

rownames(sgenescorer2_matrix)[sgenescorer2_matrix[,"Il13"]<200 & sgenescorer2_matrix[,"Il4"]<100]
#Mpv17l2  Ryr1 Nfatc1

#block rbc matur: ccdc134
#DMDD check online
#Camsap3
#got ccdc134 sections  E14.5  
#antonella ag10@ .... 

getaracnerelatedto <- function(net, genes){
  net[net$Regulator %in% genes | net$Target %in% genes,]
}
getaracnerelatedto(aracne_ko, c("Ryr1","Mpv17l2","Nfatc1"))  #  Mpv17l2  is DE!

mtpm[toensid2("Mpv17l2"),]
plot(av_mtpm["ENSMUSG00000035559",],type="l")   #should highlight  Mpv17l2
which(togenesym2(rownames(mtpm)) == "Mpv17l2")

ensconvert[ensconvert$mgi_symbol=="Mpv17l2",]
toensid2("Mpv17l2")

# plot(log10(1+sgenescorer2_matrix[,"Il13"]),onepermutation(log10(1+sgenescorer2_matrix[,"Il4"])))
# 
# plot(log10(1+sgenescorer2_matrix[,"Il13"]),log10(1+sgenescorer2_matrix[,"Gata3"]))
# plot(log10(1+sgenescorer2_matrix[,"Il13"]),onepermutation(log10(1+sgenescorer2_matrix[,"Gata3"])))
# 
#plot(log10(1+sgenescorer2_matrix[,"Il13"]),log10(1+sgenescorer2_matrix[,"Irf4"]))

#?sample



######################################################################
### KO similarity to Th_x - pinwheel #################################
######################################################################


pdf("out_ko/ko_plast_pinwheel_neutral_new.pdf")
fortsne <- t(thenorm)[,-1]
colnames(fortsne) <- c("Th1","Th17","iTreg","nTreg")
npin <- ncol(fortsne)+1  #last one is Th2
pinx <- cos((1:npin)*2*pi/npin - 5*pi/5/2)
piny <- sin((1:npin)*2*pi/npin - 5*pi/5/2)
plot(c(-2,2),c(-2,2),cex=0)
#lines(c(pinx,pinx[1]),c(piny,piny[1]), col="gray")
for(i in 1:npin){ #nrow(fortsne)
  lines(c(0,pinx[i]),c(0,piny[i]), col="gray")
}
text(pinx*1.1,piny*1.1, labels = c(colnames(fortsne),"Th2"),cex=0.8)
text(0,0, labels = "Neutral")
for(i in 1:npin){ 
  pinx[i] <- pinx[i] - pinx[npin]
  piny[i] <- piny[i] - piny[npin]
}
for(i in 1:nrow(fortsne)){
  #sc<-4
  sc<-0.16  #for the normalized
  thecol <- "black"
  if(rownames(fortsne)[i] %in% c("Cd180","Cd200","Cxcr7","F2rl1","Ifngr1","Lag3","Orm3","Scara3","Thy1"))
    thecol <- "#c80000ff"
  if(rownames(fortsne)[i] %in% c("Bcl11b","Bhlhe40","Ern1","Etv2","Gata3","Stat6","Tbx21","Xbp1","Zc3h12a"))
    thecol <- "#0000c8ff"
  if(rownames(fortsne)[i] %in% c("Abcg4","Crls1","Nhedc2","Pgk1","Pxk","Slc5a1","Slc25a3"))
    thecol <- "#b900c8ff"
  if(rownames(fortsne)[i] %in% c("Ccdc134","Il2","Il4","Il12a","Il13"))
    thecol <- "#00c800ff"
  v <- c(fortsne[i,],0)
  
  text(sc*sum(v*pinx), sc*sum(v*piny), rownames(fortsne)[i],cex=0.5, col=thecol)
}
dev.off()


######################################################################
### KO activation/differentiation ####################################
######################################################################


####### First perform the DE

#### Calculate differentiation DE: Th0/Th2 at 72h
datTcDiff <- read.csv("out_tc/mouse/all_samples_mouse_genes_Salmon_counts.txt",sep="\t",row.names = "gene")
datTcDiff <- datTcDiff[,grep("_72h",colnames(datTcDiff))]
colnames(datTcDiff)
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = round(datTcDiff),
  colData = data.frame(row.names = colnames(datTcDiff), is2 = factor(colnames(datTcDiff) %!in% colnames(datTcDiff)[grep("Th2",colnames(datTcDiff))])),
  design = ~ is2)
ddsFullCountTable <- DESeq(ddsFullCountTable)
resTcDiff <- results(ddsFullCountTable)

#### Calculate activity DE: Th0 0h vs Th2 72h   #Changed from Th2
datTcAct <- read.csv("out_tc/mouse/all_samples_mouse_genes_Salmon_counts.txt",sep="\t",row.names = "gene")
datTcAct <- datTcAct[,c(grep("Naive_",colnames(datTcAct)), grep("Th0_72h_",colnames(datTcAct)))]
ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = round(datTcAct),
  colData = data.frame(row.names = colnames(datTcAct), isact = factor(colnames(datTcAct) %!in% colnames(datTcAct)[grep("Naive",colnames(datTcAct))])),
  design = ~ isact)
ddsFullCountTable <- DESeq(ddsFullCountTable)
resTcAct <- results(ddsFullCountTable)


resTcActDiffP  <- data.frame(row.names = rownames(resTcAct), 
                             act=resTcAct$padj, 
                             diff=resTcDiff$padj)
resTcActDiffFC <- data.frame(row.names = rownames(resTcAct), 
                             act=resTcAct$log2FoldChange, 
                             diff=resTcDiff$log2FoldChange)

plot(resTcActDiffFC$act, resTcActDiffFC$diff)

### by using DE fold changes
compareDEwithRef <- function(thep, thefc){
  theall<-c()
  allsd<-c()
  for(j in 1:ncol(thep)){
    numgenetested <- c()
    print(sprintf("ref: %s",j))
    #Only consider DE genes according to reference
    totest <- !is.na(thep[,j]) & thep[,j]<1e-10
    totest <- intersect(rownames(thefc)[totest],rownames(newkoPadj))

    cmin<-function(x) sum(na.omit(x)==-1)
    cplus<-function(x) sum(na.omit(x)==1)
    kov <- c()
    for(i in 1:ncol(outpadj)){  ################### respval vs outpadj???? TODO
      #Further only consider DE genes according to our dataset. Use better cut-off
      usefc <- newkoFCnorm
      totest3 <- intersect(
        totest,
        rownames(usefc)[abs(usefc[,i])>log(2,2)])       #was 0.5

      numgenetested <- c(numgenetested,length(totest3))
      s1<-sign(thefc[totest3,j])
      s2<-sign(usefc[totest3,i])
      a<-sum(s1*s2)

      kov <- c(kov, sum(a))
    }
    print(numgenetested)
    if(j==1)  
      theall<-kov
    else
      theall<-rbind(theall,kov)
  }
  colnames(theall) <- colnames(respval)
  rownames(theall) <- sprintf("Th2 vs %s",colnames(thep))
  hist(numgenetested,breaks=20)
  
  thenorm <- theall
  for(i in 1:nrow(thenorm)){
    thenorm[i,] <- thenorm[i,]/sd(thenorm[i,])
  }
  
  list(
    theall=theall,
    thenorm=thenorm)
  #  theall#list(theall=theall)  
}

###Could generalize this function to plasticity analysis too!
compkoTcAcDiff <- compareDEwithRef(resTcActDiffP, resTcActDiffFC)


listReallyLowReprod <- c("Cxcr7","Il13","Pgk1","Ern1","Nrd1","Gata3","Fam26f","Il12a","Herc6")
listFairlyLowReprod <- c("Cd180","Il4","Mta3","Pxk","Lrrc40","Oas1f","Pxk")
listReallyGoodReprod <- c("Ifngr1","Hsph1","Ccdc134","Cd200","Bcl11b","Ccna2",
                          "Tbx21","Abcg4","Thy1","F2rl1","Fam32a","Fam71b",
                          #This last row more skeptical
                          "Bhlhe40","Orm3","Stat6","Slc5a1","Cntrob","Slc25a3","B230219D22Rik","Scara3")
vis <- compkoTcAcDiff$theall
#vis <- compkoTcAcDiff$theall[,colnames(compkoTcAcDiff$theall) %!in% listReallyLowReprod]

#colnames(theall) <- sprintf("%s KO",colnames(theall))

thecol <- rep("black",ncol(vis))
thecol[colnames(vis) %in% c("Cd180","Cd200","Cxcr7","F2rl1","Ifngr1","Lag3","Orm3","Scara3","Thy1")]    <- "#c80000ff"
thecol[colnames(vis) %in% c("Bcl11b","Bhlhe40","Ern1","Etv2","Gata3","Stat6","Tbx21","Xbp1","Zc3h12a")] <- "#0000c8ff"
thecol[colnames(vis) %in% c("Abcg4","Crls1","Nhedc2","Pgk1","Pxk","Slc5a1","Slc25a3")]                  <- "#b900c8ff"
thecol[colnames(vis) %in% c("Ccdc134","Il2","Il4","Il12a","Il13")]                                      <- "#00c800ff"

pdf("out_ko/ko_actdiff_new.pdf")
plot(vis[2,], vis[1,],cex=0,
     ylab="Activation related", xlab="Differentiation related")
lines(minmax(vis[2,]),c(0,0),col="gray")
lines(c(0,0),minmax(vis[1,]),col="gray")
#vis <- compkoTcAcDiff$theall[,colnames(compkoTcAcDiff$theall) %!in% listReallyGoodReprod]
#col <- rep("black",ncol(vis))
#col[colnames(vis) %!in% listReallyGoodReprod] <- "gray"
text(vis[2,], vis[1,],col=thecol,
     labels = colnames(vis))
dev.off()






######################################################################
### Store data in the format data.teichlab.org #######################
######################################################################

######### Sample description #################
res_samplemeta_fc <- data.frame(sample=sprintf("ko_%s_fc",colnames(respval)))
res_samplemeta_fc$target_ensembl<-colnames(respval)
res_samplemeta_fc$organism<-"mouse"
res_samplemeta_fc$'Cell Type'<-"Th2"
res_samplemeta_fc$hours<-120
res_samplemeta_fc$method<-"CRISPRKO"
res_samplemeta_fc$target_mgi_symbol<-res_samplemeta$sample

res_samplemeta_pval <- data.frame(sample=sprintf("ko_%s_pval",colnames(respval)))
res_samplemeta_pval$target_ensembl<-colnames(respval)
res_samplemeta_pval$organism<-"mouse"
res_samplemeta_pval$'Cell Type'<-"Th2"
res_samplemeta_pval$hours<-120
res_samplemeta_pval$method<-"CRISPRKO"
res_samplemeta_pval$target_mgi_symbol<-res_samplemeta$sample


write.csv(res_samplemeta_fc,"out_teichlab/th2crispr_ko_samplemeta_fc.csv",row.names = FALSE, quote = FALSE)
write.csv(res_samplemeta_pval,"out_teichlab/th2crispr_ko_samplemeta_pval.csv",row.names = FALSE, quote = FALSE)

## Which genes to store
badbifurlo <- c("ENSMUSG00000097058",toensid(c("4930556m19rik","Gm20796","Mical3","Mtcp1","Supt20")))
keep <- (rownames(respval) %!in% badbifurlo) #& apply(abs(resfc),1,min)>0.05

## KO pval
temp <- respval[keep,]
colnames(temp) <- sprintf("ko_%s_pval",colnames(temp))
temp <- temp[rownames(temp) %!in% badbifurlo,]
write.csv(temp,"out_teichlab/th2crispr_ko_pval.csv",row.names = TRUE, quote = FALSE)

## KO FC
temp <- resfc[keep,]
colnames(temp) <- sprintf("ko_%s_fc",colnames(temp))
temp <- temp[rownames(temp) %!in% badbifurlo,]
write.csv(temp,"out_teichlab/th2crispr_ko_fc.csv",row.names = TRUE, quote = FALSE)






######### Sample description #################

out <- t(compkoTcAcDiff$theall)

quickens <- function(temp){
  te <- ensconvert
  te <- te[te$mgi_symbol %in% rownames(temp),]
  te <- te[isUnique(te$ensembl_gene_id),]
  te <- te[isUnique(te$mgi_symbol),]
  rownames(te) <- te$mgi_symbol
  temp <- temp[rownames(temp) %in% te$mgi_symbol,]
  rownames(temp) <- te[rownames(temp),]$ensembl_gene_id
  temp
}  
out <- quickens(out)  #losing 2 genes here
colnames(out) <- c("th2crispr_ko_act","th2crispr_ko_diff")

samplemeta <- data.frame(sample=colnames(out))
samplemeta$organism<-"mouse"
samplemeta$'Cell Type'<-"Th2"
write.csv(samplemeta,"out_teichlab/th2crispr_ko_diffact_samplemeta.csv",row.names = FALSE, quote = FALSE)

write.csv(out, sprintf("out_teichlab/th2crispr_ko_diffact_data.csv",wm), row.names = TRUE, quote = FALSE)



######################################################################
############### Volcano of act/diff genes ############################
######################################################################


x <- as.data.frame(resTcAct)
x$ensembl_gene_id <- rownames(x)
x <- merge(x, ensconvert)
x <- x[grep("Il",x$mgi_symbol),]

plot(x$log2FoldChange, log10(x$padj),cex=0)#, ylim=c(0,-50))
text(x$log2FoldChange, log10(x$padj),labels = x$mgi_symbol)

x <- x[order(x$padj),]
x

#il2 il2ra  Cdk4   Ybx1   (eif5a)  (rps27l)    Ranbp1   Cdca7   Rpl7l1   Plk2    Lif   Cdkn1a   Cacna1i

#Focusing on act diff: Il12rb1, Il7r, Il6st, Il6ra, Il16, Il17ra,   (il10 is 2e-2)


x[order(x$mgi_symbol),]

x <- as.data.frame(resTcDiff)
x$ensembl_gene_id <- rownames(x)
x <- merge(x, ensconvert)
x <- x[grep("Il",x$mgi_symbol),]

plot(x$log2FoldChange, log10(x$padj),cex=0)#, ylim=c(-40,0), xlim=c(-12,12))
text(x$log2FoldChange, log10(x$padj),labels = x$mgi_symbol)


#Ms4a4b
#Il4
#Ccr4
#Ccr6
#Batf
#Il9r
#Il18r
#Il13, Il24
#Il10
#Il12rb1