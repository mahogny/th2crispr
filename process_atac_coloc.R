#############################################################################################
###                                                                                       ###
###             Part of the paper ...                                                     ###
###             Author: Johan Henriksson (mahogny@areta.org)                              ###
###                                                                                       ###
###             This code ...                                   ###
###                                                                                       ###
#############################################################################################


source("common_geneid.R")
source("common_gotest.R")
source("common_readrnaatacchip.R")

#load("ATAC_data.RData")
#load("RNA_data.RData") 

#chipbatf <- read.csv("../chip/batf_72h.csv",sep="\t",stringsAsFactors = FALSE)
#chipirf4 <- read.csv("../chip/irf4_72h.csv",sep="\t",stringsAsFactors = FALSE)
###chipbatf$Nearest.Ensembl  this could go into promoter activation analysis!

###coloc - gene or peak level? trivial to integrate into 

#infoTable_ATAC$mouse$
#no info on peak size. so definitely need to do coloc on a per-gene level
#which peaks are different Th0/Th2? where is the Th0 data? can we say that the corresponding TFs are 
#important for Th2/Th0 differentiation and run a test?
###TODO: for time course analysis, use the table 1/0. also, any of them working together?
#if we have TF vs gene/time, can we separate out different pathways?
###########whee. position in peak file is relative to that snippet. not global!!
### sites which a TF target - are there TFs for which the genes come up together?
# can categorize TFs as broad in time and in space. 2d plot. can also cluster genes based on
#their space/time pattern. genes x time as feature matrix
#xi: previous plot on regions the ATAC peaks are in. add the plot? can add a time axis make a figure
#like I did for the dogseq paper?
# plot(apply(infoTable_ATAC$mouse[,3:8],2,sum),type="l",ylim=c(0,50e3)) #why going down later? stabilizing?
# head(infoTable_ATAC$mouse)
# unique(infoTable_ATAC$mouse$feature)
#   
# ### Read mapping gene -> motifID
# genemotif <- read.csv("JASPAR2016_MA_PB_C2H2_nonred.meme.names",header=FALSE,sep=" ",stringsAsFactors = FALSE)[,c(2,3)]
# colnames(genemotif) <- c("motifid","tf")
# normalizesym <- function(s) paste(str_sub(s,1,1),str_to_lower(str_sub(s,2)),sep="")
# for(i in 1:nrow(genemotif))
#   genemotif$tf[i] <- normalizesym(genemotif$tf[i])
# 
# 
# ### Read motifs detected in peaks
# peakpat <- read.csv("fimo.txt",stringsAsFactors = FALSE,sep="\t")
# head(peakpat)
# colnames(peakpat)[1]<-"motifid"
# colnames(peakpat)[2]<-"peak"
# peakpat$peakid <- str_to_lower(peakpat$peakid)
# 
# unique(peakpat$motifid)
# x <- merge(peakpat,genemotif)
# y<-x[,c(2,5,7,8,10)]
# hist(log10(x$p.value[x$tf=="Batf3"]))




#de_early <- read.csv("early_DE_Th0Th2_genes.txt",sep="\t",stringsAsFactors = FALSE)
#de_late <- read.csv("late_DE_Th0Th2_genes.txt",sep="\t",stringsAsFactors = FALSE)



#######################################################################################
######## when are TFs activated? check which regulate DE genes ########################
#######################################################################################

#Only consider expressed TFs
testtf <- intersect(genemotif$tf,expressed_atacTF_50)
## it gets quite boring like this. original tested all TFs

## Function to perform bootstrap test of one TF driving DE genes
f2 <- function(settot, seta,setb){
  #settot <- #union(seta,setb)
  m<-matrix(data=c(length(settot), length(setdiff(seta,setb)),
                   length(setdiff(setb,seta)), length(intersect(seta,setb))),ncol = 2)
  p<-m[1,2]/(m[1,1]+m[1,2])
  q<-m[2,1]/(m[1,1]+m[2,1])
  exp2 <- p*q*sum(as.double(m))
  c(fisher.test(m)$p.value,m[2,2]>exp2)
}


### Test promoters of early genes
ep <- c()
epa <- c()
for(i in 1:length(testtf)){
  print(i)
  oneft <- f2(
    unique(de_early$ens_gene),
    unique(de_early[de_early$qval<1e-2,]$ens_gene),
    infoTable_ATAC$mouse[infoTable_ATAC$mouse$peak %in% x$peak[x$tf==testtf[i]],]$TSS_ensg)
  ep <- c(ep,oneft[1]) 
  epa <- c(epa,oneft[2]) 
}
names(ep) <- testtf
names(epa) <- testtf



### Test promoters of late genes
lp <- c()
lpa <- c()
for(i in 1:length(testtf)){
  print(i)
  oneft <- f2(
    unique(de_late$ens_gene),
    unique(de_late[de_late$qval<1e-2,]$ens_gene),
    infoTable_ATAC$mouse[infoTable_ATAC$mouse$peak %in% x$peak[x$tf==testtf[i]],]$TSS_ensg)
  lp <- c(lp,oneft[1]) 
  lpa <- c(lpa,oneft[2]) 
}
names(lp) <- testtf
names(lpa) <- testtf


hist(log10(1e-50+ep),breaks=20)
hist(log10(ep[ep!=0]),breaks=20) #-e300 is cool. -320 is cool

elp <- merge(
  data.frame(name=names(ep),p.early=ep, e.over=epa),
  data.frame(name=names(lp),p.late=lp, p.over=lpa))
elp<-elp[order(pmax(elp$p.early,elp$p.late)),]
elp[elp$p.early< 1e-300 | elp$p.late< 1e-300,]
#so all of these over expected

#elp[(elp$p.early< 1e-300 | elp$p.late< 1e-300) & !elp$p.over,]

earlytf <- sort(names(sort(ep[ep< 1e-300])))
latetf <- sort(names(sort(lp[lp< 1e-300])))
unique(earlytf)
unique(latetf)
setdiff(earlytf,latetf)
onlylatetf <- unique(setdiff(latetf,earlytf)) #strictly only more coming on
intersect(earlytf,latetf)

write.csv(collapsegenenames(as.character(unique(elp$name[elp$p.early<1e-300]))),
          "collapsed_atac_regearly.red.csv")
write.csv(collapsegenenames(as.character(unique(elp$name[elp$p.late<1e-300 & elp$p.early>1e-300]))),
          "collapsed_atac_reglate.red.csv")


################### genes DE early or late? 
#library(limma)
allgene<-unique(c(de_early$ext_gene, de_late$ext_gene))
vd <- data.frame(
  early=allgene %in% de_early$ext_gene[de_early$qval<1e-2],
  late=allgene %in% de_late$ext_gene[de_late$qval<1e-2]
)
vc <- vennCounts(vd)
vennDiagram(vc,cex=c(2,2,2))
pastewspace <- function(sdata) do.call(paste, c(as.list(sdata), sep=" "))


#stat6. is stat6 a tf???

#Bcl6 comes on late

#Make a new peak table as a heatmap over genes? maybe no need

#Zbtb14 is here too. 
#Ikzf3 early. Etv6. Fosl2 also, from Th17 paper





#######################################################################################
########################## co-loc correlation between TFs #############################
#######################################################################################

########### Create input data for java program
write.csv(x[,c(2,10,3,4,5)],"out_coloc/out_peaktf.csv",quote = FALSE)
#### TODO still called x?  put in common_

########### Read java program output
# tfcoloc <- read.csv("tf_coloc.jaccard.csv",stringsAsFactors = FALSE)
rownames(tfcoloc)<-colnames(tfcoloc)
for(i in 1:nrow(tfcoloc))
  tfcoloc[i,i]<-0
#tfcoloc[1:30,1:30]
tfcoloc <- as.matrix(tfcoloc)
tfcoloc[is.infinite(tfcoloc)]<-0


#mi <- apply(tfcoloc,1,mean)>0.05
mi <- apply(tfcoloc>0.5,1,any)
mc <- tfcoloc[mi,mi]
heatmap.2(mc,#tfcoloc,
          trace="none",
          density.info="none", 
          #Colv="NA",
          #dendrogram="rows",
          col=colorRampPalette(c("red", "green","blue"))(n = 300))
#hist(as.double(tfcoloc),breaks=50)
#as.double(tfcoloc)


### read fisher test input
fa <- read.csv("out_coloc/tf_coloc.a.csv",stringsAsFactors = FALSE)
fb <- read.csv("out_coloc/tf_coloc.b.csv",stringsAsFactors = FALSE)
fc <- read.csv("out_coloc/tf_coloc.c.csv",stringsAsFactors = FALSE)
fd <- read.csv("out_coloc/tf_coloc.d.csv",stringsAsFactors = FALSE)


####Calculate fisher statistic
fname_ftexpd <- "out_coloc/tf_ftexpd.csv"
if(!file.exists(fname_ftexpd)){
  ft <- matrix(nrow=nrow(fa), ncol=nrow(fa))
  rownames(ft)<-colnames(fa)
  colnames(ft)<-colnames(fa)
  for(i in 1:nrow(fa)){
    print(i)  
    for(j in 1:nrow(fa)){
      oneft <- fisher.test(matrix(data = c(fa[i,j],fb[i,j],fc[i,j],fd[i,j]),nrow = 2))
      ft[i,j] <- oneft$p.value
    }
  }
  ft_p <- fb/(fa+fb)
  ft_q <- fc/(fa+fc)
  ft_exp_d <- ft_p*ft_q*(fa+fb+fc+fd)
  
  write.csv(ft_exp_d,fname_ftexpd)
  write.csv(ft,"out_coloc/tf_ft.csv")
} else {
  #### Read fisher stats
  ft_exp_d <- read.csv(fname_ftexpd,row.names = 1)
  ft <- read.csv("out_coloc/tf_ft.csv",row.names = 1)
  rownames(ft_exp_d) <- colnames(ft_exp_d)
  ft <- as.matrix(ft)
}

#### Only keep TFs that are expressed. OBSOLETE AND ANNOYING
# keeptf <- rep(TRUE,ncol(ft))#colnames(ft) %in% expressed_atacTF_20
# keeptf <- colnames(ft) %in% expressed_atacTF_20
# #keeptf <- colnames(ft) %in% expressed_atacTF_50
# ft_exp_d <- ft_exp_d[keeptf,keeptf]
# ft <- as.matrix(ft[keeptf, keeptf])
# fd <- as.matrix(fd[keeptf, keeptf])

##### Rescale correlations to log, and make anti-correlation positive on the log scale
for(i in 1:nrow(ft))
  ft[i,i]<-1
lft <- log10(ft)
lft[is.infinite(lft)]<- -50
invi <- fd<ft_exp_d              
lft[invi] <- -lft[invi]

#### Only keep TFs that are expressed
keeptf <- rep(TRUE,ncol(ft))#colnames(ft) %in% expressed_atacTF_20
keeptf <- colnames(ft) %in% expressed_atacTF_20
lft_red <- lft[keeptf,keeptf]
##TODO: use the red form below, if any

########### find representative genes
collapsegenenames <- function(n){
  for(cg in c("Znf","Tbx","Fox","Zbtb","Znp","Olig","Zscan","Tp","Hox","Atf","Vax","Alx","Evx",
              "Gbx","Sp","Etv","Pou","Jun","Lmx","Nfat","Tcf","Vsx","Klf","Elf","Lmx","Zfp","Tgif","Msx",
              "Fos","Creb","Tfap","Lhx","Tead","Maf","Yy","Pknox","Scrt","Stat","Irf","Meox","Phox","Neuro",
              "Gsx","Pax","E2f","Egr","Zic","Elk")){
    i<-grep(cg,n)
    if(length(i)>1){
      nn <- do.call(paste, c(as.list(
        c(n[i[1]],str_sub(n[i[2:length(i)]], str_length(cg)+1))
      ), sep="/"))
      n <- c(n[-i],nn)
    }
  }
  n
}

numclust<-30
#mi <- 1:nrow(lft)#  apply(abs(lft) >10,1,any) #throw away the least interesting ones at least
mi <- apply(abs(lft) >20,1,function(x)sum(x)>3) #throw away the least interesting ones at least
mi <- 1:nrow(lft)
mi <- (rownames(lft) %in% earlytf | rownames(lft) %in% latetf)
#nrow(lft)
sum(mi)
set.seed(123) 
km <- kmeans(lft[mi,mi],numclust)
pasten <- function(sdata) do.call(paste, c(as.list(sdata), sep=" "))
repgenes <- c()
for(i in 1:numclust){
  n <- names(which(km$cluster==i))
  rn <- pasten(collapsegenenames(n))
  repgenes <- c(repgenes,rn)
}
if(FALSE){
  set.seed(123)
  rtsne_out <- Rtsne(as.matrix(lft[mi,mi]), pca = FALSE, verbose = TRUE,perplexity = 100)
  cols <- palette(rainbow(numclust))
  plot(rtsne_out$Y, asp = 1, pch = 20,# col = "blue",
       cex = 0.75, cex.axis = 1.25, cex.lab = 1.25, cex.main = 1.5,
       col = cols[km$cluster],
       xlab = "t-SNE dimension 1", ylab = "t-SNE dimension 2",
       main = "2D t-SNE projection")
}

#### Ideally construct an average correlation
avcor <- matrix(nrow=numclust,ncol=numclust)
for(i in 1:numclust)
  for(j in 1:numclust)# set.seed(123) 

    avcor[i,j] <- mean(as.double(lft[km$cluster==i, km$cluster==j]))
colnames(avcor) <- repgenes
rownames(avcor) <- repgenes

heatmap.2(avcor,
          trace="none",
          density.info="none", 
          hclustfun=function(d) hclust(d, method="ward.D2"),  
          margins=c(15,15),
          col=colorRampPalette(c("red", "green","blue"))(n = 300))



########### plot individual fisher p values
mi <- apply(abs(lft) >200,1,any)
#mi <- apply(abs(lft) >100,1,any)
mi <- apply(abs(lft) >10,1,any)
#mi <- 1:nrow(lft)
isznf <- rownames(lft) %in% rownames(lft)[grep("Znf",rownames(lft))]
isztf <- rownames(lft) %in% rownames(lft)[grep("Zt",rownames(lft))]
iszb <- rownames(lft) %in% rownames(lft)[grep("Zb",rownames(lft))]
#mi <- apply(abs(lft) >200,1,function(x) sum(x)>10)  & !isznf
mi <- rownames(lft) %in% onlylatetf & !isztf
mi <- (rownames(lft) %in% earlytf | rownames(lft) %in% latetf) & !isztf & !isznf & !iszb

mi <- (rownames(lft) %in% earlytf | rownames(lft) %in% latetf) & apply(lft >10,1,any) & !(rownames(lft) %in% rownames(lft)[grep("Z",rownames(lft))])
mi <- (rownames(lft) %in% earlytf | rownames(lft) %in% latetf) & apply(lft >10,1,any)# & !(rownames(lft) %in% rownames(lft)[grep("Z",rownames(lft))])
#mi <- rownames(lft) %in% repgenes & apply(abs(lft) >100,1,function(x) sum(x)>5)

#mi <- apply(abs(lft) >100,1,function(x)sum(x)>3) #throw away the least interesting ones at least
#pdf("heatmap coloc corr earlylate.pdf")
############ Supplementary figure really
png("heatmap coloc corr earlylate.png",width=5000, height = 5000)
mi <- (rownames(lft) %in% earlytf | rownames(lft) %in% latetf) & apply(lft >10,1,any) & !(rownames(lft) %in% rownames(lft)[grep("Zn",rownames(lft))])
mi <- rep(TRUE,nrow(lft))
plot.new()
pdf("heatmap coloc corr earlylate.red.pdf",width=20, height = 20)
heatmap.2(as.matrix(lft[mi,mi]),
          trace="none",
          density.info="none", 
          hclustfun=function(d) hclust(d, method="ward.D2"),   #ward.D2?
          col=colorRampPalette(c("red", "green","blue"))(n = 300), 
          margins=c(15,15),
          cex=0.5,
          dendrogram = "none",
          keysize = 0.5,
          cexRow = 1, cexCol = 1)
dev.off()


png("heatmap coloc corr earlylate reduced.png",width=3000, height = 3000)
#mi <- (rownames(lft) %in% earlytf | rownames(lft) %in% latetf) & apply(lft >10,1,any) & !(rownames(lft) %in% rownames(lft)[grep("Zn",rownames(lft))])
mi <- rownames(lft) %in% c("Ctcf","Elf5","Etv6","Elf4","Fos","Fosl1","Junb","Spic","Spi1",
                           "Irf2","Irf7","Irf8","Prdm6","Foxp3","Foxo1","Zic2_2","Zfp740_1","Plag1",
                           "Elk1","Elk4","Klf13","Etv5",
                           "Insm","Yy1","Nr2c2","Ikzf3","Rest","Zbtb7c","Thap1",
                           "Elk4","Nrf","Hinfp","E2f4","Klf10","Klf12","Klf14","Glis1","Klf4",
                           "Zbtb14","Zbtb26","Zbtb6",
                           "Batf3","Batf..jun")  ##these are not really there
#what about batf?
heatmap.2(lft[mi,mi],
          trace="none",
          density.info="none", 
          hclustfun=function(d) hclust(d, method="ward.D2"),   #ward.D2?
          col=colorRampPalette(c("red", "green","blue"))(n = 300), 
          margins=c(20,20),
          cex=3,
          dendrogram = "none",
          keysize = 0.5,
          cexRow = 5, cexCol = 5)
dev.off()


#which genes correlate similarly? correlate by correlation
hist(as.double(cor(lft,method="spearman")))


colnames(lft)[grep("Batf",colnames(lft))]
colnames(lft)[grep("Ben",colnames(lft))]





#Jundm2_2/Fos/Fosl2/Junb  are anticorr with Hinfp/Zfp161_1/Tcfl5
#Batf:Jun anticorrelated with E2f4 and other things

#log10(x$p.value[x$tf=="T"])
#head(genemotif)

###hm. for each motif pairs, check overlap. as a function of time???

###check for a given gene the spread in time for its genes





################################################################################################
###### Plot ATAC peak categories over time
################################################################################################

gethstat <- function(f) read.csv(f,sep="\t",stringsAsFactors=FALSE,header=TRUE)[1:12,]
htimes <- c(0,2,4,24,48,72)
hstat <- c()
for(i in 1:length(htimes)){
  f <- sprintf("hstat_%sh",htimes[i])
  if(i==1)
    hstat <- gethstat(f)[,c(1,2)]
  else
    hstat <- cbind(hstat,gethstat(f)[,2])
}
colnames(hstat) <- c("n",htimes)
rownames(hstat) <- hstat[,1]
hstat <- hstat[,-1]
hstat <- as.matrix(hstat)
class(hstat)<-"numeric"

mulfac <- c(100,5000,100,50,1000,
       10,1,1,1,10,
       1,3000)
for(i in 1:nrow(hstat))
  hstat[i,] <- mulfac[i]*hstat[i,]
# for(i in 1:ncol(hstat))
#   hstat[,i] <- hstat[,i]/sum(hstat[,i])
plot(hstat[1,],type="l",xaxt = 'n',ylim=c(0,0.6))
for(i in 2:nrow(hstat))
  lines(1:ncol(hstat),hstat[i,],type="l",xaxt = 'n')
colnames(hstat) <- sprintf("%sh",htimes)
axis(side=1,at=1:ncol(hstat),labels=sprintf("%sh",htimes))
pdf("out_atacvst/atac_cat_rescale.pdf",width = 3)
barplot(hstat[-11,], col = brewer.pal(11,"Set3"), legend.text = TRUE)
dev.off()


############################################################################
#### Find potential interactors of Irf4 ####################################
############################################################################


### Compute DE vs motif-by-atac interaction
findDEatacpartners <- function(jasparname){
  irfi <- sort(lft[,jasparname])
  tc <- data.frame(cooc = as.double(irfi), jasparname=names(irfi))
  irfi <- smerge(tc,motif_explevel)
  sqldf("select distinct max(cooc) as cooc, mgi_symbol, maxexp from irfi group by mgi_symbol")
}

### Plot DE vs motif-by-atac interaction
plotDEatacpartners <- function(irfi, showaboveexp=2.4, showabovecooc=5.2, highlight){
  plott<-abs(irfi$cooc)>showabovecooc | log10(1+irfi$maxexp)>showaboveexp | irfi$mgi_symbol %in% highlight
  plot(log10(1+irfi$maxexp), irfi$cooc,cex=1*binarize(!plott),pch=19,col="gray",
       ylab="IRF4 interaction index",xlab="Log10 TPM")
  text(log10(1+irfi$maxexp)[plott], irfi$cooc[plott],labels = irfi$mgi_symbol[plott], cex=0.7)
}
colnames(deatac_irf)


pdf("out_coloc/irf4_potpart.pdf")
deatac_irf <-findDEatacpartners("Irf4_1")
deatac_irf <- deatac_irf[-grep("Irf",deatac_irf$mgi_symbol),]
plotDEatacpartners(deatac_irf, highlight=c("Spib","Stat6","Nfatc1","Nfatc2","Spi1","Pou6f1"))
dev.off()


## I don't think this is good enough to be used yet. definitely need to integrate position information
pdf("out_coloc/etv2_potpart.pdf")
deatac_etv2 <-findDEatacpartners("Etv2")
deatac_etv2 <- deatac_etv2[-grep("Irf",deatac_etv2$mgi_symbol),]  #Irf1 similar
deatac_etv2 <- deatac_etv2[!(deatac_etv2$mgi_symbol %in% c("Elf3","Elf5","Etv5","Ehf")),]  #very similar
deatac_etv2 <- deatac_etv2[!(deatac_etv2$mgi_symbol %in% c("Spic","Spi1")),]  #quite similar
plotDEatacpartners(deatac_etv2, highlight=c("Spib","Stat6","Nfatc1","Nfatc2","Spi1","Pou6f1"))
dev.off()


### should cite
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2920730/
#about irf4 being an integrator than a regulator





####
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2920730/
### This review suggests some modes of interaction for Th2,0,17
  
### Know to interact with Stat6, Nfatc1/2 and Batf
#so also Runx1? Stat1,3,4  foxp1  listed here. all candidates looking at motif



####### Foxo1 on Irf4 and IL9 promoter in Th9. claimed a master regulator of Th9
#http://www.jimmunol.org/content/196/1_Supplement/133.23
##

##############
#http://www.nature.com/onc/journal/v19/n21/full/1203485a.html
#Bcl6 and Stat6 bind to the same sites and compete
#

##############
# http://www.cell.com/cell-reports/pdf/S2211-1247(16)31099-3.pdf   macrophages
## Stat6 & Stat3 -> Ire1a?  
# induced, and part of UPR: NFE2, ATF6, Rfx-family members, HOXA5/B5, NR1H4, ATF2, ARNT, and XBP1
# But in my view, large amount of new protein synthesis -> UPR should be activated
# These guys also do MARA



