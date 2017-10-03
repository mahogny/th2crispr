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



