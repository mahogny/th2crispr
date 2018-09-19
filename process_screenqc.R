#############################################################################################
###                                                                                       ###
###             Part of the paper ...                                                     ###
###             Author: Johan Henriksson (mahogny@areta.org)                              ###
###                                                                                       ###
###             This code ...                                   ###
###                                                                                       ###
#############################################################################################


# source("common_geneid.R")
# source("common_gotest.R")
# source("common_readrnaatacchip.R")
# 

################################################################################
############ Gather screen pos and neg lists ###################################
################################################################################

fname_pneg <- "out_screenqc/screen_pneg.csv"
fname_ppos <- "out_screenqc/screen_ppos.csv"
if(!file.exists(fname_pneg)){
  files <- c("s9_stg", "s8a_stl", "sx2_stl", "s10_foxp3", "sx4_irf4",
             "first_il4","s8b_xbp1","s11_il13", "allneg","sc1_il13",
             "sc2a_irf4","sc2b_xbp1","sc3_gata3","gatasurv1","gatasurv2") 
  snegp<-c()
  sposp<-c()
  for(curf in 1:length(files)){
    print(files[curf])
    scr <- read.csv(sprintf("screen/mageck/%s.gene_summary.txt",files[curf]),sep="\t",stringsAsFactors=FALSE)
    tn<-scr[,c("id","neg.p.value")]
    tp<-scr[,c("id","pos.p.value")]
    colnames(tn)<-c("id",files[curf])
    colnames(tp)<-c("id",files[curf])
    if(curf==1){
      snegp<-tn
      sposp<-tp
    } else {
      snegp<-merge(snegp,tn)
      sposp<-merge(sposp,tp)
    }
  }
  rownames(snegp)<-snegp$id
  rownames(sposp)<-sposp$id
  snegp<-snegp[,-1]
  sposp<-sposp[,-1]
  snegp<-snegp[order(rownames(snegp)),]
  sposp<-sposp[order(rownames(sposp)),]
  
  ####### write to file, supplementary
  write.csv(snegp, fname_pneg)
  write.csv(sposp, fname_ppos)
} else {
  ######### read screening results ##########
  snegp<-read.csv(fname_pneg,row.names = 1)
  sposp<-read.csv(fname_ppos,row.names = 1)
}
screenispos <- sposp>snegp



################################################################################
############ Score screens #####################################################
################################################################################

## Get rank of each gene in the screens
srank<-sposp[,keepscreens]
for(i in 1:ncol(srank)){
  srank[,i] <- pmin(rank(srank[,i]),rank(snegp[,keepscreens][,i]))
}

## Get pval of each gene in the screens
spval<-sposp[,keepscreens]
for(i in 1:ncol(spval)){
  spval[,i] <- pmin(spval[,i],snegp[,keepscreens][,i])
}


## Take the lower rank as the consensus score
getscreenconsensus <- function(x) sort(x)[2]
getscreenconsensus_m <- function(x) apply(x,1,getscreenconsensus)

sgenescorer       <- apply(srank,1,getscreenconsensus)
sgenescorer_il4   <- apply(srank[,screens_il4],  1,getscreenconsensus)
sgenescorer_il13  <- apply(srank[,screens_il13], 1,getscreenconsensus)
sgenescorer_irf4  <- apply(srank[,screens_irf4], 1,getscreenconsensus)
sgenescorer_xbp1  <- apply(srank[,screens_xbp1], 1,getscreenconsensus)
sgenescorer_gata3 <- apply(srank[,screens_gata3],1,getscreenconsensus)
sgenescorer_matrix <- data.frame(
  Il4=sgenescorer_il4, 
  Il13=sgenescorer_il13, 
  Irf4=sgenescorer_irf4, 
  Xbp1=sgenescorer_xbp1, 
  Gata3=sgenescorer_gata3)
sgenescorer_total <- apply(sgenescorer_matrix,1,min)


## Take the pval rank as the consensus score
sgenescorep       <- apply(spval,1,getscreenconsensus)
sgenescorep_il4   <- apply(spval[,screens_il4],  1,getscreenconsensus)
sgenescorep_il13  <- apply(spval[,screens_il13], 1,getscreenconsensus)
sgenescorep_irf4  <- apply(spval[,screens_irf4], 1,getscreenconsensus)
sgenescorep_xbp1  <- apply(spval[,screens_xbp1], 1,getscreenconsensus)
sgenescorep_gata3 <- apply(spval[,screens_gata3],1,getscreenconsensus)
sgenescorep_matrix <- data.frame(
  il4  =sgenescorep_il4, 
  il13 =sgenescorep_il13, 
  irf4 =sgenescorep_irf4, 
  xbp1 =sgenescorep_xbp1, 
  gata3=sgenescorep_gata3)
sgenescorep_total <- apply(sgenescorep_matrix,1,min)




### Alternative scoring method which better handles if there is only one hit
getscreenconsensus2 <- function(x) gm_mean(sort(x)[1:2])
getscreenconsensus2_m <- function(x) apply(x,1,getscreenconsensus2)

sgenescorer2_il4   <- apply(srank[,screens_il4],  1,getscreenconsensus2)
sgenescorer2_il13  <- apply(srank[,screens_il13], 1,getscreenconsensus2)
sgenescorer2_irf4  <- apply(srank[,screens_irf4], 1,getscreenconsensus2)
sgenescorer2_xbp1  <- apply(srank[,screens_xbp1], 1,getscreenconsensus2)
sgenescorer2_gata3 <- apply(srank[,screens_gata3],1,getscreenconsensus2)
sgenescorer2_matrix <- data.frame(
  Il4=sgenescorer2_il4, 
  Il13=sgenescorer2_il13, 
  Irf4=sgenescorer2_irf4, 
  Xbp1=sgenescorer2_xbp1, 
  Gata3=sgenescorer2_gata3)
sgenescorer2_matrix <- round(sgenescorer2_matrix)

#This is used for the network building later - be careful in QC
for(i in list_screen_genes){
  sgenescorer2_matrix[i,i]<-1
}

# sgenescorer2_matrix["Gata3","Gata3"]<-1
# sgenescorer2_matrix["Il4",  "Il4"]  <-1
# sgenescorer2_matrix["Il13", "Il13"] <-1
# sgenescorer2_matrix["Xbp1", "Xbp1"] <-1
# sgenescorer2_matrix["Irf4", "Irf4"] <-1


################################################################################
############ go analysis (+miRNA) for all screens ##############################
################################################################################

#Perform GO analysis for each screen
fname_screengo <- "out_screenqc/screens_go.csv"
if(!file.exists(fname_screengo)){
  outgo <- list()
  for(curf in 1:ncol(sposp)){
    #curf<-1
    print(curf)
#    ishit <- snegp[,curf]<1e-1 | sposp[,curf]<1e-1   #why 5 at end?   was 1e-2 below
    ishit <- snegp[,curf]<quantile(snegp[,1],0.15) | sposp[,curf]<quantile(sposp[,1],0.15)   #why 5 at end?   was 1e-2 below
    bg <- rownames(sposp)[which(rownames(sposp) %in% expressedGenes)]
    oneset <- rownames(sposp)[which(ishit & rownames(sposp) %in% expressedGenes)]
    x <- stopgosym(oneset,bg,cutoff=1e-1)
    outgo[[curf]] <- x
  }
  names(outgo) <- colnames(sposp)
  
  ### Combine GO analyses into a matrix instead
  allcat <- c()
  for(i in 1:length(outgo))
    allcat <- union(allcat, outgo[[i]]$term.name)
  screengo <- matrix(nrow=length(allcat),ncol=ncol(sposp))
  rownames(screengo) <- allcat
  colnames(screengo) <- colnames(sposp)  
  for(i in 1:length(outgo))
    screengo[outgo[[i]]$term.name,i] <- outgo[[i]]$p.value
  screengo[is.na(screengo)] <- 1

  ### Save data  
  write.csv(screengo,fname_screengo)  
} else {
  screengo <- read.csv(fname_screengo, row.names = 1)
}



### Read list of interesting categories
intcat <- unique(read.csv("screen_interestinggocat.txt",header = FALSE,stringsAsFactors = FALSE)[,1])
intcat_mi <- intcat[grep("MI",intcat)]
intcat <- intcat[-grep("MI",intcat)]



########## Plot GO heatmaps
screengo <- screengo[order(apply(log(screengo),1,function(x) sum(na.omit(x)))),]
heatmap.2(log(screengo[intcat,keepscreens]),
          hclustfun=function(d) hclust(d, method="ward.D2"),  
          trace="none",density.info="none", 
          margins=c(8,20),
          col=colorRampPalette(c("red", "yellow","white","blue","pink"))(n = 300))
m2 <- binarize(log(screengo[intcat,keepscreens])< -15)
m3 <- binarize(log(screengo[intcat,keepscreens])< -8)
m2 <- log(screengo[intcat,keepscreens])[apply(m2,1,sum)>0 | apply(m3,1,sum)>1,]
colnames(m2) <- keepscreens_ren
png("out_screenqc/heatmap go screens normalgene.png",w=1500,h=2000)
heatmap.2(m2,
          hclustfun=function(d) hclust(d, method="ward.D2"),  
          trace="none",density.info="none", 
          margins=c(15,50),
          cexRow = 2,
          cexCol = 3,
          col=colorRampPalette(c("#7a0177","#c51b8a","#f768a1","#fbb4b9","#feebe2",
                                 "#feebe2","blue","blue","blue","blue"))(n = 33))
dev.off()
pdf("out_screenqc/heatmap go screens normalgene.pdf",w=6)
heatmap.2(m2,
          hclustfun=function(d) hclust(d, method="ward.D2"),  
          trace="none",density.info="none", 
          margins=c(15,15),
          cexRow = 0.45,
          cexCol = 0.5,
          col=colorRampPalette(c("#7a0177","#c51b8a","#f768a1","#fbb4b9","#feebe2",
                                 "#feebe2","blue","blue","blue","blue"))(n = 33))
dev.off()

heatmap.2(log(screengo[intcat_mi,keepscreens]),
          hclustfun=function(d) hclust(d, method="ward.D2"),  
          trace="none",density.info="none", 
          margins=c(8,20),
          col=colorRampPalette(c("red", "yellow","white","blue","pink"))(n = 300))
m2 <- log(screengo[intcat_mi,keepscreens])< -15
m2[m2]<-1
m2[!m2]<-0
m2 <- log(screengo[intcat_mi,keepscreens])[apply(m2,1,sum)>0,]

rownames(m2) <- str_replace(rownames(m2),"MI:","")
colnames(m2) <- str_replace(colnames(m2),"_stg"," Gata3")
colnames(m2) <- str_replace(colnames(m2),"_gata3"," Gata3")
colnames(m2) <- str_replace(colnames(m2),"_xbp1"," Xbp1")
colnames(m2) <- str_replace(colnames(m2),"_il13"," IL13")
colnames(m2) <- str_replace(colnames(m2),"_irf4"," Irf4")
colnames(m2) <- str_replace(colnames(m2),"_stl","  IL4")
colnames(m2) <- str_replace(colnames(m2),"_il4","  IL4")

pdf("out_screenqc/heatmap go screens mir.pdf")
heatmap.2(m2,
          hclustfun=function(d) hclust(d, method="ward.D2"),      # change key manually to show p-value
          trace="none",density.info="none",  
          cexCol = 0.8,
          cexRow = 0.8,
          #          colsep = 1:ncol(m2),
          #          rowsep = 1:nrow(m2),
          margins=c(8,20),
          col=colorRampPalette(c("#7a0177","#c51b8a","#f768a1","#fbb4b9","#feebe2",
                                 "#feebe2","blue","blue","blue","blue"))(n = 33))
dev.off()

pdf("~/Desktop/newscreengo.pdf")
screengo_red <- log(screengo[,keepscreens])
colnames(screengo_red) <- keepscreens_ren2
screengo_red <- screengo_red[-grep("Factor",rownames(screengo_red)),]
head(sort(screengo_red[,"Xbp1 b"]))
screengo_red <- screengo_red[apply(screengo_red,1,function(x) sort(x)[2])< 0,]
rownames(screengo_red)
screengo_red
heatmap.2(screengo_red,
          hclustfun=function(d) hclust(d, method="ward.D2"),  
          trace="none",density.info="none", 
          margins=c(8,20),
          col=colorRampPalette(c("red", "yellow","white","blue","pink"))(n = 300))
dev.off()




################################################################################
############ The big overview heatmap of all screens ###########################
################################################################################

########## Pick hits for each screen separately
cutoff <- 280
overrideinclude <- c(colnames(respval),"Slc9b2","Ern1","Ackr3")#,,"Cxcr4","Ire1", "Fam71b","Cxcr7","Ackr3","Ern1")
#overrideinclude %in% rownames(srank)
#overrideinclude[!(overrideinclude %in% rownames(r))]
keepgenes <- unique(c(overrideinclude,
  names(which((sgenescorer_il4)<cutoff*0.7)),
  names(which((sgenescorer_il13)<cutoff)),
  names(which((sgenescorer_irf4)<cutoff)),
  names(which((sgenescorer_xbp1)<cutoff)),
  names(which((sgenescorer_gata3)<cutoff*1.2))))
length(keepgenes)
r2<-srank[keepgenes,] 
colnames(r2) <- keepscreens_ren2

### only consider expressed genes
r3<-r2[rownames(r2) %in% c(overrideinclude,expressedGenes),]
pdf("out_screenqc/all.pdf",width=3,height=16)
min((as.matrix(r3)))
max((as.matrix(r3)))
png("out_screenqc/all.png",width=2000,height=2000)
heatmap.2(log(1+as.matrix(r3)),
          trace="none",
          dendrogram = "none", 
          # lhei = c(0.02,1),
          # lwid = c(0.02,1),
          # Colv=FALSE,
          # cexCol=1.2,
          # cexRow=1.2,
          col=colorRampPalette(c("#dd1c77","#c994c7", "#e7e1ef","white"))(n = 300))
dev.off()


############## what is the overlap between the screens?
cutoff2 <- 1000 
expressedGenes <- th[apply(th[,-1]>10,1,sum)>0,1]
vd <- cbind(
  sgenescorer_il4  < sort(sgenescorer_il4  [rownames(r) %in% expressedGenes])[cutoff2] & rownames(r) %in% expressedGenes,
  sgenescorer_il13 < sort(sgenescorer_il13 [rownames(r) %in% expressedGenes])[cutoff2] & rownames(r) %in% expressedGenes,
  sgenescorer_irf4 < sort(sgenescorer_irf4 [rownames(r) %in% expressedGenes])[cutoff2] & rownames(r) %in% expressedGenes,
  sgenescorer_xbp1 < sort(sgenescorer_xbp1 [rownames(r) %in% expressedGenes])[cutoff2] & rownames(r) %in% expressedGenes,
  sgenescorer_gata3< sort(sgenescorer_gata3[rownames(r) %in% expressedGenes])[cutoff2] & rownames(r) %in% expressedGenes
)
colnames(vd) <- c("IL4","IL13","IRF4","XBP1","GATA3")
vc <- vennCounts(vd)
vennDiagram(vc,cex=c(1.5,1.5,1.5))
## Correlation heatmap
oc <- cor(vd,method="spearman")
for(i in 1:ncol(oc))
  oc[i,i]<-0
heatmap(oc)
##least correlated: Gata3/IL3
##most correlated:  Xbp1/Il13 & irf4/IL4 & xbp1/il4
heatmap.2(oc,
  trace="none",
  density.info="none",
  #Rowv=FALSE,
  dendrogram = "none",
  #Colv="NA",
  #dendrogram="rows",
  col=colorRampPalette(c("white","blue"))(n = 300))

plot(hclust(dist(t(vd))))
#500: Il13/xbp1 & il4/irf4. gata3 is an outlier
#1000: same as 500. 

#Could also do a heatmap based on correlations.





### new score: max of the min scores!
hist(as.double(as.matrix(log10(w))),breaks=30)



####What about for the expressed TFs? ... some ATAC peak overlap cut-off too?
hist(apply(r2,1,function(x) sort(x)[6]))
r2<-r[rownames(r) %in% motif_explevel$mgi_symbol[motif_explevel$maxexp>20],]  
r2<-r2[apply(r2,1,function(x) sort(x)[6])<4000,]
pdf("out_screenqc/overview_atactf.pdf",width=3,height=10)
heatmap.2(log(1+as.matrix(r2)),
          trace="none",
          dendrogram = "none", 
          lhei = c(0.02,1),
          lwid = c(0.02,1),
          Colv=FALSE,
          cexCol=0.8,
          col=colorRampPalette(c("#dd1c77","#c994c7", "#e7e1ef","white"))(n = 300))
dev.off()
##Pou6f1, Yy1, Nfatc1, Fli1, Bhlhe40, Foxo3, Hinfp, Zfp523 are generally up to things. xbp1 hit for il4 (as known)
##and: Mnt, Smad3, 
#Bhlhe40 & etv6 & ewsr1 seem to attack Irf4?
#The best screen suggest several

#Fam71b close to il4
#Zbtb
#ddit3,ormdl2,olfr*** close to stat6

#Zfp523 close to abcg1 and notch3
#Lag3, Vmn* close to Cd4
#Thy1, Abcg4 & Hinfp nearby each other!!!   and Bcl9l
#Pou2f2 close to Vmn* and Vmn1r181

#these were also in the good screen:
# Vmn1r195           61.5
# Vmn1r232           61.5
# Vmn2r102           61.5
# Vmn2r104           61.5

head(r[order(r[,"sc3 Gata3"]),][,"sc3 Gata3",drop=FALSE],n=200)



################################################################################
############ Overview heatmap of the KO genes ##################################
################################################################################

########## Pick hits for each screen separately
kon <- colnames(respval)
kon[kon=="Cxcr7"]<-"Ackr3"
kon[kon=="Nhedc2"]<-"Slc9b2"

keepgenes <- kon
r2<-srank[kon,] 
colnames(r2) <- keepscreens_ren2

### only consider expressed genes
pdf("out_ko/screenscores.pdf",width=2.5,height=8)
heatmap.2(log(1+as.matrix(r2)),
          trace="none",
          dendrogram = "none", 
          lhei = c(0.02,1),
          lwid = c(0.02,1),
          Colv=FALSE,
          cexCol=1.2,
          cexRow=1.2,
          col=colorRampPalette(c("#dd1c77","#c994c7", "#e7e1ef","white"))(n = 300))
dev.off()



################################################################################
##################### plot an example screening vs FC ##########################
################################################################################


getscreenPFC <- function(sname){
  p<-sposp[,sname]
  p[p>snegp[,sname]] <- -snegp[p>snegp[,sname],sname]
  p <- data.frame(id=rownames(sposp),p=p,stringsAsFactors=FALSE)
  fc <- data.frame(id=th$id,fc = log10(th$t72/th$t0), tpm=apply(th[,-1],1,min), stringsAsFactors=FALSE)
  p <- merge(p,fc)
  p <- p[p$tpm>15,] #gata3 is 24. stat6 is 104
  p$p <- sign(p$p)*log10(abs(p$p)) #convert p values into log10
  p
}

plotPFC <- function(p){
  plab <- p$p>2.3 | p$p < -2.1
  plot(p$p[!plab],p$fc[!plab],col="blue",pch=19,cex=0.1,ylim=c(-1,1),xlim=c(-6,6),ylab="Fold change",xlab="Log10 p-value")
  ###draw gray line cross hair
  pcol <- rep("black",nrow(p))
  pcol[p$id %in% c("Gata3","Stat6")]<-"red"
  #pcol[p$id=="Stat6"]<-"red"
  psize <- rep(0.5,nrow(p))
  psize[which(pcol=="red")] <- 0.8
  text(p$p[plab],p$fc[plab],labels=p$id[plab],col=pcol[plab],cex=psize[plab])
}

pdf("out_screenqc/ex_sc3_gata3.pdf")
plotPFC(getscreenPFC("sc3_gata3"))
dev.off()



################################################################################
##################### plot an example screening vs Th20 pval ###################
################################################################################

getscreenPFC_pval <- function(sname){
  sname <- "sc3_gata3"
  p<-sposp[,sname]
  p[p>snegp[,sname]] <- -snegp[p>snegp[,sname],sname]
  p <- data.frame(id=rownames(sposp),p=p,stringsAsFactors=FALSE)
  fc <- data.frame(id=tcmouse$de_early$mgi_symbol,fc = tcmouse$de_early$pval, stringsAsFactors=FALSE)
  fc <- sqldf("select id, min(fc) as fc from fc group by id")
  p <- merge(p,fc)
  expgenes <- ensconvert$mgi_symbol[ensconvert$ensembl_gene_id %in% rownames(tcmouse$av_mtpm)[tcmouse$av_mtpm>15]]
#  print(expgenes)
#  print(head(p))
  p <- p[p$id %in% expgenes,]
  p$fc <- log10(p$fc) 
  #print(head(p))
  #print(6666)
#  p <- p[p$tpm>15,] #gata3 is 24. stat6 is 104
  p$p <- sign(p$p)*log10(abs(p$p)) #convert p values into log10
  p
}

plotdottext_ex <- function(x,y,labels,cex,col=rep("black",length(x)),donew=FALSE){
  
  xnorm <- x/(max(x)-min(x))  #ratio
  ynorm <- y/(max(y)-min(y))
  pd <- data.frame(x=xnorm,y=ynorm)
  pd <- as.matrix(dist((pd)))
  
  #666
  #Local density
  nn <- apply(pd<0.05,1,sum)
  asp <- nn>3
  #Make exceptions
  cd <- apply(pd,1,function(v) sort(v)[2])
  asp[cd>0.1] <- FALSE
  
  asp[col!="black"] <- FALSE
  
  if(donew)
    plot(x,y,cex=0)
  points(x[ asp], y[ asp],pch = 20,cex=0.5,col="gray")
  text(  x[!asp], y[!asp],labels = labels[!asp],cex=cex,col=col[!asp])
  #  text(  x, y,labels = nn,cex=cex,col=col)
}


plotPFC_pval <- function(p){
  plot(p$p, p$fc,cex=0,
       xlab="Screen hit Log10 p-value",
       ylab="DE Log10 p-value")
  col <- rep("black",length(p$p))
  col[p$id %in% c("Gata3","Stat6")] <- "red"
  plotdottext_ex(p$p, p$fc, p$id, cex=1, col=col)
}

pdf("out_screenqc/ex_sc3_gata3_th02.pdf")
plotPFC_pval(getscreenPFC_pval("sc3_gata3"))
dev.off()


# plab <- p$p>2.7 | p$p < -2.6 | p$fc < -21
# plot(p$p,p$fc,cex=0,xlim=c(-6,6),ylab="DE Log10 p-value",xlab="Screen hit Log10 p-value")
# points(p$p[!plab],p$fc[!plab],col="blue",pch=19,cex=0.1)
# ###draw gray line cross hair
# pcol <- rep("black",nrow(p))
# pcol[p$id %in% c("Gata3","Stat6")]<-"red"
# #pcol[p$id=="Stat6"]<-"red"
# psize <- rep(0.5,nrow(p))
# psize[which(pcol=="red")] <- 0.8
# text(p$p[plab],p$fc[plab],labels=p$id[plab],col=pcol[plab],cex=psize[plab])




################################################################################
#################### DE genes vs hit scores ####################################
################################################################################

plot_de_screen <- function(sde,sele,sell,fearly,limit=-1,tcol=NULL,plott=NULL){
  #Set color of labels
  if(is.null(tcol)){
    tcol <- rep("",nrow(sde))
    tcol[sele] <- "red"         #Only early
    tcol[sell] <- "#00DD00"       #Only late (green)
    tcol[sele & sell] <- "blue" #Both early and late
  }
  if(is.null(plott)){
    plott <- tcol!=""
  }
  tcol[tcol==""] <- "black"
  print(sum(plott))
  print(tcol[plott])  
  
  shifty <- function(y){
    s <- y>limit
    y[s] <- y[s] + (y[s]-limit)*2
    y
  }

  ############ Plot early, p-value based
  pdf(fearly,width=14,height=7)
  par(mfrow=c(1,2))
  plot(log10(sde$sp), log10(sde$de_early),cex=0,xlab="Screen p-value",ylab="Early DE p-value",ylim=c(-27,2))
  points(log10(sde$sp)[!plott], log10(sde$de_early)[!plott],cex=0.2,pch=19, col=tcol[!plott])
  text(log10(sde$sp)[plott], shifty(log10(sde$de_early)[plott]), 
       labels = sde$mgi_symbol[plott], cex=tsize, col=tcol[plott])

  ############ Plot late, p-value based
  plot(log10(sde$sp), log10(sde$de_late),cex=0,xlab="Screen p-value",ylab="Late DE p-value",ylim=c(-27,2))
  points(log10(sde$sp)[!plott], log10(sde$de_late)[!plott],cex=0.2,pch=19)
  text(log10(sde$sp)[plott], shifty(log10(sde$de_late)[plott]), 
       labels = sde$mgi_symbol[plott], cex=tsize, col=tcol[plott])
  dev.off()
}



### Combine screen score and DE score
de1 <- na.omit(sqldf("select distinct `ext_gene` as mgi_symbol, ens_gene as id, min(qval) as p from de_early group by mgi_symbol"))
de2 <- na.omit(sqldf("select distinct `ext_gene` as mgi_symbol, ens_gene as id, min(qval) as p from de_late  group by mgi_symbol"))
sde <- data.frame(mgi_symbol=names(sgenescorep_total),sp=as.double(sgenescorep_total), sr=as.double(sgenescorer_total), stringsAsFactors = FALSE)
sde <- smerge(sde, data.frame(mgi_symbol=de1$mgi_symbol, de_early=de1$p, stringsAsFactors = FALSE), all.x = TRUE)
sde <- smerge(sde, data.frame(mgi_symbol=de2$mgi_symbol, de_late =de2$p, stringsAsFactors = FALSE), all.x = TRUE)
sde$de_early[is.na(sde$de_early)] <- 1
sde$de_late [is.na(sde$de_late)] <- 1

#TODO: replot with the 1 added

### Suitable for all genes
sde2 <- sde[sde$mgi_symbol %in% expressedGenes,]
sele <- log10(sde2$sp)*10 + log10(sde2$de_early) < -27 | log10(sde2$sp)< -2.2 | sde2$mgi_symbol=="Stat6"
sell <- log10(sde2$sp)*10 + log10(sde2$de_late)  < -27 | log10(sde2$sp)< -2.2 | sde2$mgi_symbol=="Stat6"
sum(sele)
tsize=0.7
plot_de_screen(sde2,sele,sell,
               "out_screenqc/compare_de_screen.pdf")
#rnf19b cool. http://www.informatics.jax.org/diseasePortal/popup?isPhenotype=true&markerID=MGI:1922484&header=immune%20system

### Suitable for membrane proteins
sde2 <- sde[sde$id %in% c(list_protatlas_membrane) & sde$mgi_symbol %in% expressedGenes,]
sele <- log10(sde2$sp)*10 + log10(sde2$de_early) < -27 + 6
sell <- log10(sde2$sp)*10 + log10(sde2$de_late)  < -27 + 6
sum(sele)
tsize=0.7
plot_de_screen(sde2,sele,sell,
               "out_screen_qc/compare_de_screen_memb.pdf")

### Suitable for secreted proteins
sde2 <- sde[sde$id %in% c(list_protatlas_secreted) & sde$mgi_symbol %in% expressedGenes,]
sele <- log10(sde2$sp)*10 + log10(sde2$de_early) < -27 + 6
sell <- log10(sde2$sp)*10 + log10(sde2$de_late)  < -27 + 6
sum(sele)
tsize=0.7
plot_de_screen(sde2,sele,sell,
               "out_screenqc/compare_de_screen_secreted.pdf")

### Suitable for TFs only
sde2 <- sde[sde$id %in% c(list_tf,list_co,list_cm) & sde$mgi_symbol %in% expressedGenes,]
sele <- log10(sde2$sp)*10 + log10(sde2$de_early) < -27 + 8
sell <- log10(sde2$sp)*10 + log10(sde2$de_late)  < -27 + 8
sum(sell)
tsize=0.7
plot_de_screen(sde2,sele,sell,
               "out_screenqc/compare_de_screen_tfs.pdf")


### Suitable for expressed TFs seen in MARA (NOTE: different different cut-off)
sde2 <- sde[sde$id %in% toensid(motif_explevel$mgi_symbol[motif_explevel$maxexp>20]),]
sele <- log10(sde2$sp)*10 + log10(sde2$de_early) < -27 + 13
sell <- log10(sde2$sp)*10 + log10(sde2$de_late)  < -27 + 13
sum(sell)
tsize=0.8
plot_de_screen(sde2,sele,sell,
               "out_screenqc/compare_de_screen_atacmotif.pdf")




################################################################################
#################### Merge all types of scores for TFs #########################
################################################################################

## Connect motif gene and ATAC score
motif_ael <- merge(merge(merge(
  data.frame(jasparname=genemotif$tf, stringsAsFactors = FALSE), 
  map_jaspar_namegenesym),ensconvert, stringsAsFactors = FALSE),
  data.frame(jasparname=names(score_ael),ael=as.double(score_ael)), stringsAsFactors = FALSE)
motif_ael <- unique(motif_ael)
motif_ael <- sqldf("select jasparname, mgi_symbol, avg(ael) as mael from motif_ael group by mgi_symbol")

## Merge late DE, ATAC early/late and screen score
sael <- data.frame(mgi_symbol=names(sgenescorep_total),sp=as.double(sgenescorep_total), sr=as.double(sgenescorer_total), stringsAsFactors = FALSE)
sael <- smerge(sael,motif_ael)
sael <- smerge(sael, data.frame(mgi_symbol=de2$mgi_symbol, de_late=de2$p, stringsAsFactors = FALSE), all.x = TRUE)
sael$de_late[is.na(sael$de_late)] <- 1
sael


#### Plot it with screen pval vs ATAC early/late. then color by DE score
pdf("out_screenqc/compare_ael_screen_tfs.pdf")
plott <- rep(FALSE,nrow(sael))
plott[sael$sp < 10^-1.8 | abs(sael$mael)>0.5 | sael$mgi_symbol %in% c("Nfil3","Irf4") ]<-TRUE
plot(log10(sael$sp), sael$mael,cex=0,xlab="Screen p-value",ylab="ATAC early/late")#,ylim=c(-27,2))
points(log10(sael$sp)[!plott], sael$mael[!plott],cex=0.2,pch=19,col="blue")
text(log10(sael$sp)[plott], sael$mael[plott], 
     labels = sael$mgi_symbol[plott], cex=0.8, col="red")
dev.off()

#### Plot it with screen pval vs DE pval
sde2 <- sde[sde$mgi_symbol %in% c("Etv2",motif_explevel$mgi_symbol[motif_explevel$maxexp>0.2]),]
sde2 <- smerge(sde2, motif_ael)
elc <- 0.9
sele <- log10(sde2$sp)*10 + log10(sde2$de_early) < -27 + 7 | abs(sde2$mael)>elc
sell <- log10(sde2$sp)*10 + log10(sde2$de_late)  < -27 + 7 | abs(sde2$mael)>elc
tcol <- rep("",nrow(sde2))
tcol <- col_from_ael(sde2$mael)
plott <- sell | sele
sum(sell)
tsize<-0.8
plot_de_screen(sde2,sele,sell,
               "out_screenqc/compare_de_screen_ael_atacmotif.pdf",tcol=tcol,plott=plott,limit = 1)



#### Plot it with screen pval vs ATAC early/late. then color by DE score (lowest early/late)
pdf("out_screenqc/compare_ael_screen_tfs2.pdf")
des <- scale01(-log10(pmin(sde2$de_early,sde2$de_late)))^0.3
#hist(des)
tcol <- rgb(des,0,0)
plott <- rep(FALSE,nrow(sael))
#plott[des>0.75]<-TRUE
plott[sde2$sp < 10^-1.8 | abs(sde2$mael)>0.5 | sde2$mgi_symbol %in% c("Nfil3","Irf4","Fli1") | des>0.8]<-TRUE
tcol[!plott] <-"gray"
plot(log10(sde2$sp), sde2$mael,cex=0,xlab="Screen p-value",ylab="ATAC early/late")
points(log10(sde2$sp)[!plott], sde2$mael[!plott],cex=0.2,pch=19,col=tcol[!plott])
text(log10(sde2$sp)[plott], 
     sde2$mael[plott], 
     labels = sde2$mgi_symbol[plott], cex=0.8, col=tcol[plott])
dev.off()

#Gata3: only one both differentially expressed, a screen hit and "atac early/late". a bit late
#Npas2/Atf3 might come next (a bit early) but pretty crappy
#Fli1 might be more interesting then

#TODO: an online interactive plot with cut-offs for DE/screen score
#screen score, which method to use too. can precalculate and upload
