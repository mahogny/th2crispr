#Typical size of an irf4/batf peak is 600bp


###### for gata3, might want 2 out 4, as it is shakier


### New take at analyzing; what is the Batf/Irf4 overlap? is there a stat6/gata3 overlap?
#stat6 looks particularly difficult
#gata3 closer
####################################################################
##### gata + batf + irf4 + stat6 ###################################
####################################################################
colnames(dchiptot)
#dchiptot <- readallchiptot(fname="gbi")
ic <- grep("_72h",colnames(dchiptot))
#ip <- which(apply(dchiptot[,ic]>0,1,any))
vd <- dchiptot[,ic]
colnames(vd) <- c("Gata3","Batf","Irf4","Stat6")
vc <- vennCounts(vd)
pdf("out_chip/gbis_venn.pdf")
vennDiagram(vc,cex=c(1.5,1.5,1.5))
dev.off()

ipall <- apply(vd>0,1,sum)==4
sort(dchiptot[ipall,]$Gene.Name)
dchiptot[ipall,c("Gene.Name","Start","Distance.to.TSS")]



#72h only: all agree on
#Mir101c - yes. but actually in intron of Gm26870. other Gm nearby. this is a lincRNA. blast shows many copies
#Lrrc4c - actually peak is on Gm10800. 1.3M away from gene. should probably have a cut-off!
#redo MARA with cut-offs?


####################################################################
##### gata + batf + irf4 + xbp1 ####################################
####################################################################
dgbix <- readallchiptot(fname="chip3")
ic <- grep("_72h",colnames(dgbix))
vd <- dgbix[,ic]
colnames(vd) <- c("Gata3","Batf","Irf4","Xbp1")
vc <- vennCounts(vd)
pdf("out_chip/gbix_venn.pdf")
vennDiagram(vc,cex=c(1.5,1.5,1.5))
dev.off()

#xbp1 & gata3 fairly separate. gata3 would be a bit more pro-irf4 than pro-batf, if anything
corjaccard(vd)

genes_gbix <- dgbix[apply(vd>0,1,sum)>=4,c("Gene.Name","Chr","Start","Distance.to.TSS")]
genes_gbix

#Venn just gata3 & xbp1
vd <- dgbix[,c("Gata3","Xbp1")]
vc <- vennCounts(vd)
vennDiagram(vc,cex=c(1.5,1.5,1.5))

cs <-data.frame(
  mgi_symbol=dgbix$Gene.Name,
  score_gata3=apply(dgbix[,grep("Gata3",colnames(dgbix))],1,sum),
  score_xbp1=dgbix$Xbp1_72h)
cs <- cs[cs$score_gata3>1 & cs$score_xbp1>0,]
cs <- cs[order(cs$score_gata3),]
cs

#Rora, Il2rb is a target with score 3
#Bach2 score 1.5

####################################################################
##### gata + batf + irf4 + merged stat6 ############################
####################################################################
ic <- c(grep("_72h",colnames(dgbis)))
vd <- as.data.frame(binarize(dgbis[,ic]>0))
colnames(vd) <- c("Gata3","Batf","Irf4","Stat6")
vc <- vennCounts(vd)
pdf("out_chip/gbism_venn.pdf")
vennDiagram(vc,cex=c(1.5,1.5,1.5))
dev.off()

corjaccard(vd) #Stat6 is closer to Irf4 than others

genes_gbis <- dgbis[apply(vd,1,sum)>=4,c("Gene.Name","Start","Distance.to.TSS")]
#(B3gnt2,) Mir7001, Mir101c, Rn45s

ipall <- vd$Stat6m_72h & apply(vd,2,sum)<=1 & abs(dgbis$Distance.to.TSS)<50e3
genes_gbis <- dgbis[ipall,c("Gene.Name","Chr","Start","Distance.to.TSS")]
genes_gbis
#Ccr2? maybe, with gata3 & stat6 in the same peak

ipall <- vd$Stat6m_72h & apply(vd,2,sum)<=1 & abs(dgbis$Distance.to.TSS)<50e3
genes_gbis <- dgbis[ipall,c("Gene.Name","Chr","Start","Distance.to.TSS")]
genes_gbis
#Oas1b might be real!



################# Compare with screen score
ipall <- apply(dgbis[,ic],1,sum)>=3 & abs(dgbis$Distance.to.TSS)<100e3
#sort(dgbi[ipall,]$Gene.Name)
cs <- cbind(
  sgenescorer2_matrix[dchiptot[ipall,c("Gene.Name")],],
  #sgenescorer_matrix[dchiptot[ipall,c("Gene.Name")],],
  dgbis[ipall,c("Gene.Name","Chr","Start","Distance.to.TSS")])
cs <- na.omit(cs[apply(cs[,1:5]<800,1,any),])
cs <- cs[order(apply(cs[,1:5],1,function(x)sort(x)[2])),]
unique(cs$Gene.Name)
head(cs,n=100)


# [1] "Rnpep"         "Atoh8"         "Rybp"          "Ddit4"         "Vdac1"        
# [6] "Krt222"        "Dclk1"         "Cwc27"         "Fyn"           "Vwa5a"        
# [11] "1190005I06Rik" "Cdyl"          "Smim12"        "Chrnb4"        "Cd44"         
# [16] "Xbp1"          "Elovl5"        "Gabrb3"        "Stt3a"         "Myo5c"        
# [21] "Itgb1"         "Cytip"         "Cxcr4"         "Cd82"     


####################################################################
##### gata + batf + irf4 ###########################################
####################################################################
#For stat6: could take 72 & 48 together. or consensus of all

#colnames(dchiptot)

ic <- grep("_72h",colnames(dgbi))
#ip <- which(apply(dgbi[,ic]>0,1,any))
vd <- binarize(dgbi[,ic]>0)
colnames(vd) <- c("Gata3 72h","Batf 72h","Irf4 72h")
vc <- vennCounts(vd)
pdf("out_chip/gbi_venn.pdf")
vennDiagram(vc,cex=c(1.5,1.5,1.5))
dev.off()

vecjaccard(vd[,1], vd[,2])
vecjaccard(vd[,1], vd[,3])  
vecjaccard(vd[,3], vd[,2])  #no collaboration gata3  vs batf/irf4


# testPeakOverlap <- function(seta, setb){
#   totalsites <- length(seta)
#   n <- length(intersect(mygeneset, genesinterm))
#   m <- length(setdiff(mygeneset, genesinterm))
#   tm <- matrix(c(
#     n, length(genesinterm)-n,
#     m, length(backgroundset)-m
#   ),nrow=2)
#   fisher.test(tm)$p.value
# }


genes_gbi <- sort(unique(dgbi[apply(vd,1,sum)>=3,]$Gene.Name))
genes_gbi

sort(intersect(expressedGenes,genes_gbi))

#DE early & have a peak with GBI
sort(intersect(de_early$ext_gene[de_early$qval<1e-2],genes_gbi))  
#DE late  & have a peak with GBI. This list is longer
sort(intersect(de_late$ext_gene[de_late$qval<1e-2],genes_gbi))  


#differentially expressed?
#IL4, Il20, 

################# Compare with screen score
ipall <- apply(dgbi[,ic],1,sum)>=3 & abs(dgbi$Distance.to.TSS)<100e3
#sort(dgbi[ipall,]$Gene.Name)
cs <- cbind(
    sgenescorer2_matrix[dchiptot[ipall,c("Gene.Name")],],
  #sgenescorer_matrix[dchiptot[ipall,c("Gene.Name")],],
  dgbi[ipall,c("Gene.Name","Chr","Start","Distance.to.TSS")])
cs <- na.omit(cs[apply(cs[,1:5]<800,1,any),])
cs <- cs[order(apply(cs[,1:5],1,function(x)sort(x)[2])),]
sort(unique(cs$Gene.Name))
head(cs,n=100)





####################################################################
##### Gata3 timecourse analysis ####################################
####################################################################

timesgata3 <- c(0, 24, 48, 72,   0, 24, 48, 72)
repgata3   <- c(1,  1,  1,  1,   2,  2,  2,  2)
cntgata3 <- read.table("chip/chipcov_gata3.csv",stringsAsFactors = FALSE)
colnames(cntgata3)[1:4] <- c("Chr","Start","End","Peakid") 
cntgata3.ann <- read.csv("chip/chipcov_gata3.csv.ann",sep="\t",stringsAsFactors = FALSE)
colnames(cntgata3.ann)[1]<-"Peakid"

#Reorder to match up rows. not the same order in files!!!
cntgata3 <- cntgata3[order(cntgata3$Peakid),]
cntgata3.ann <- cntgata3.ann[order(cntgata3.ann$Peakid),]
#all(cntgata3$Peakid==cntgata3.ann$Peakid)


#size factor normalization. not using the controls in the end
cntgata3.av <- cntgata3[,c(5+c(1,2,3))] + cntgata3[,c(9+c(1,2,3))]
for(i in 1:ncol(cntgata3.av)){
  cntgata3.av[,i] <- cntgata3.av[,i] / sum(cntgata3.av[,i])
}
#normalize each peak against initial level
cntgata3.sum <- cntgata3.av[,1]
cntgata3.ann$totalcount <- apply(cntgata3.av,1,sum)
for(i in 1:ncol(cntgata3.av)){
  cntgata3.av[,i] <- cntgata3.av[,i] / cntgata3.sum
}
cntgata3.ann$inc <- cntgata3.av[,3]
cntgata3.ann.noinf <- cntgata3.ann[!is.infinite(cntgata3.ann$inc),]
cntgata3.ann.noinf$incorder <- order(cntgata3.ann.noinf$inc)/nrow(cntgata3.ann.noinf)

plot(1:3,c(0,0,0),ylim=c(-2,2),cex=0)
for(i in 1:nrow(cntgata3.av)){
  lines(1:3,log10(cntgata3.av[i,]))
}

#order(c(1,2,3))

#the most decreasing
go_gata3_dec <- stopgosym(
  genelist = unique(cntgata3.ann$Gene.Name[order(cntgata3.av[,3])[1:1000]]),
  bg = unique(cntgata3.ann$Gene.Name))

#the most increasing
go_gata3_inc <- stopgosym(
  genelist = unique(cntgata3.ann$Gene.Name[order(decreasing = TRUE,cntgata3.av[,3])[1:1000]]),
  bg = unique(cntgata3.ann$Gene.Name))

#the most increasing, no infinity
stopgosym(
  genelist = unique(cntgata3.ann.noinf$Gene.Name[order(decreasing = TRUE,cntgata3.ann.noinf$inc)[1:1000]]),
  bg = unique(cntgata3.ann$Gene.Name))

#peaks overall
go_gata3_all <- stopgosym(
  genelist = unique(cntgata3.ann$Gene.Name),
  bg = unique(ensconvert$mgi_symbol))

####### where are the interesting peaks? cases
cntgata3.ann.noinf[which(cntgata3.ann.noinf$Gene.Name=="Il4"),]
cntgata3.ann.noinf[which(cntgata3.ann.noinf$Gene.Name=="Il13"),]
cntgata3.ann.noinf[which(cntgata3.ann.noinf$Gene.Name=="Tbx21"),]  #### one is 0.017
cntgata3.ann.noinf[which(cntgata3.ann.noinf$Gene.Name=="Gata3"),]

####### where are the interesting peaks? rank
x<-cntgata3.ann.noinf[order(cntgata3.ann.noinf$inc,decreasing = TRUE),]
x<-x[x$totalcount>200,]
nrow(x)
x[1:100,c(2:4,9,10,16,20,21)]

sgenescorer2_matrix["Col10a1",]

cntgata3.ann.noinf[order(cntgata3.ann.noinf$inc,decreasing = FALSE),][1:100,c(2:4,9,10,16,20,21)]  #most decreasing
# cntgata3.ann.noinf[order(cntgata3.ann.noinf$inc,decreasing = TRUE), ][1:200,c(2:4,9,10,16,20,21)]

# increasing Fam134a,   Gpr132,   1700054M17Rik

# which(cntgata3.ann$Gene.Name=="Il13")
#chr11 53617069 53618098

### Volcano style

qtextscatter(
  log10(cntgata3.ann.noinf$inc), 
  log10(1+cntgata3.ann.noinf$totalcount), 
  cntgata3.ann.noinf$Gene.Name,
  cex=0.5)

### Store data for supplementary

write.table(cntgata3.ann, "chip/out_gata3/dynamics.csv", row.names = FALSE, sep="\t")
write.table(go_gata3_inc, "chip/out_gata3/go_inc.csv",   row.names = FALSE, sep="\t")
write.table(go_gata3_dec, "chip/out_gata3/go_dec.csv",   row.names = FALSE, sep="\t")
write.table(go_gata3_all, "chip/out_gata3/go_all.csv",   row.names = FALSE, sep="\t")

