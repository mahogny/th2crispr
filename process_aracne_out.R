
ncount2 <- dat
for(i in 1:ncol(ncount2)){
  ncount2[,i] <- ncount2[,i]/sum(ncount2[,i])*1e6
}



rownames(fpkm_thymoma)

fpkm_thymoma_av <- apply(matrix2genesym(fpkm_thymoma, ensconvert_ = human_ensconvert, removeclashing = TRUE),1,mean)
ncount_mean <- apply(matrix2genesym(ncount2, ensconvert_ = ensconvert, removeclashing = TRUE),1,mean)
#names(ncount_mean) <- togenesym2(names(ncount_mean))
compare_tcga <- smerge(
  data.frame(mgi_symbol=names(ncount_mean), ko=as.double(ncount_mean), stringsAsFactors = FALSE),
  data.frame(mgi_symbol=names(fpkm_thymoma_av), tcga=fpkm_thymoma_av, stringsAsFactors = FALSE), all = TRUE)

#Il2ra = Cd25. Cd28 has two ENSids
compare_tcga[compare_tcga$mgi_symbol %in% 
               c("Gata3","Tbx21","Irf4","Il4","Il13","Nfil3","Stat6","Stat5a","Stat5b",
                 "Cd4","Cd3e","Il2ra","Cd8"),]
#ensconvert$ensembl_gene_id[ensconvert$mgi_symbol=="Cd28"]
#rownames(ncount)


hist((as.double(fpkm_thymoma[toensid2("Gata3", ensconvert_ = human_ensconvert),]))) #20
hist((as.double(fpkm_thymoma[toensid2("Tbx21", ensconvert_ = human_ensconvert),]))) #1
hist((as.double(fpkm_thymoma[toensid2("Batf", ensconvert_ = human_ensconvert),])),breaks=50) #10
hist((as.double(fpkm_thymoma[toensid2("Irf4", ensconvert_ = human_ensconvert),])),breaks=50) #2
hist((as.double(fpkm_thymoma[toensid2("Il4", ensconvert_ = human_ensconvert),])),breaks=50) #0.2
hist((as.double(fpkm_thymoma[toensid2("Il13", ensconvert_ = human_ensconvert),])),breaks=50) #0.05



fpkm_thymoma_av[toensid2("Gata3", ensconvert_ = human_ensconvert)]

head(fpkm_thymoma)
rownames(fpkm_thymoma)

# head(fpkm_thymoma)
# sum(apply(fpkm_thymoma,1,mean)>1)
# #hist(log(1+apply(fpkm_thymoma,1,mean)),breaks=50)
# 
# filestart <- "tcga_thymoma"
# write.table(fpkm_thymoma,sprintf("aracne/%s.exp",filestart),quote = FALSE,row.names = FALSE,sep="\t")
# 
# list_tf_red <- togenesym2(list_tf)
# list_tf_red <- data.frame(tf=intersect(list_tf_red,rownames(fpkm_thymoma)),stringsAsFactors = FALSE)
# write.table(list_tf_red,sprintf("aracne/%s.tf",filestart),quote = FALSE, row.names = FALSE, col.names = FALSE)
# ## Full list of genes
# write.table(rownames(fpkm_thymoma),sprintf("aracne/%s.genes",filestart),quote = FALSE, row.names = FALSE, col.names = FALSE)



####### Distance analysis, genes to calculate
# getdistfor <- unique(c("Irf4","Xbp1","Batf","Gata3","Tbx21","Il13","Il4",colnames(respval),togenesym(list_tf)))
# write.table(getdistfor,"aracne/listgenes.csv",row.names = FALSE, col.names = FALSE, quote = FALSE)




###############################################################################################
###### Distance analysis
###############################################################################################


#aracne_mousetc <- read.aracne_net("aracne/allgene/mousetc.out/network.txt")
#aracne_humantc <- read.aracne_net("aracne/allgene/humantc.out/network.txt")  
#aracne_ko <- read.aracne_net("aracne/allgene/ko.out/network.txt")




aracne_humantc$Regulator <- normalizesym(aracne_humantc$Regulator)  ##Maybe resave it
aracne_humantc$Target <- normalizesym(aracne_humantc$Target)



adist_ko <- as.matrix(read.csv("aracne/ko.out/network.dist",sep="\t",stringsAsFactors = FALSE, row.names = 1))
#adist_ko[1:3,1:3]
table(as.double(as.matrix(adist_ko)))

plotdistfromgene <- function(forgene){
  plot(density(bw=0.1, as.double(adist_ko[forgene,])),ylim=c(0,5))
  th <- rownames(sgenescorer2_matrix)[order(sgenescorer2_matrix[,forgene])] [1:n]
  lines(density(bw=0.1, as.double(adist_ko[forgene,colnames(adist_ko) %in% th])),col="blue")
}

###############################################################################################
###### Distance analysis
###############################################################################################

n<-50
plotdistfromgene("Xbp1")
plotdistfromgene("Gata3") 
plotdistfromgene("Il13")
plotdistfromgene("Irf4")
#plotdistfromgene("Il4")
#### Xbp1 and Gata3 more central than il13 and irf4. irf4 very distant. hits not closer than average

sort(adist_ko["Gata3",])
#sort(adist_ko["Stat6",])

adist_ko["Gata3","Stat6"]

"Stat6" %in% colnames(adist_ko)
"Nfil3" %in% colnames(adist_ko)
#mtpm[toensid("Stat6"),]


###############################################################################################
###### Neightbour analysis - KO
###############################################################################################


checkrelated <- function(aracne, gene){
#  v<-aracne[aracne$Target==gene | aracne$Regulator==gene,]
  v1<-aracne[aracne$Target==gene,]
  v2<-aracne[aracne$Regulator==gene,]
  v1<-cbind(v1$Regulator,v1)
  v2<-cbind(v2$Target,v2)
  colnames(v1)[1]<-"mgi_symbol"
  colnames(v2)[1]<-"mgi_symbol"
  v<-rbind(v1,v2) 
  v <- v[order(v$MI, decreasing = TRUE),]
  v <- v[order(v$pvalue),]
  v$rank <- 1:nrow(v)
  v
}

v<-checkrelated(aracne_ko, "Gata3")
#aracne_ko[aracne_ko$Regulator=="Gata3",]
#interesting high: Junb Nfatc3 Ikzf3
#aracne_ko[aracne_ko$Target=="Gata3",]
v

aracne_ko[aracne_ko$Regulator=="Irf4",]
#interesting high: Yy1 Ormdl3
aracne_ko[aracne_ko$Target=="Irf4",]
#Yy1

aracne_ko[aracne_ko$Regulator=="Batf",]
#Leprot superhigh 

aracne_ko[aracne_ko$Regulator=="Batf3",]
#Tnfrsf8 highest. Fli1 third. Rorc really high too. 
aracne_ko[aracne_ko$Target=="Batf3",]

aracne_ko[aracne_ko$Regulator=="Xbp1",]
aracne_ko[aracne_ko$Target=="Xbp1",]

checkrelated(aracne_ko, "Stat6")


################ Which genes are the most connected?
nn<-sqldf("select `Regulator`, count(`Regulator`) as c from aracne_ko group by `Regulator` order by c")
hist(nn$c,breaks=50) #over 400 is cool
rownames(nn) <- nn$Regulator
nn

# 456          Tbx6 588  down in many activated types
# 457         Creb1 638
# 458         Nr4a1 653
# 459        Zbtb32 677 up in all but Th2
# 460           Yy1 683 ****
# 461        Zfp943 689  high in all
# 462 2810021J22Rik 697  high in all T cells
# 463         Klf13 869

mm<-sqldf("select `Regulator`, count(`Target`) as c from aracne_ko group by `Target` order by c")
hist(mm$c,breaks=50) #over 40/50 is cool
tail(mm,n=50)
# 6177     Gata3  55
# 6178      Rere  55
# 6179    Zbtb43  56
# 6180     Nfkb1  56
# 6181    Pou2f1  57 high in mara?
# 6182    Tada2a  58
# 6183    Zfp943  58
# 6184    Zbtb24  66
# 6185     Pias1  66 inhibits Stat1 https://en.wikipedia.org/wiki/PIAS1
# 6186   Zkscan1  66
# 6187    Nfatc2  71 high in all
# 6188    Setdb2  74 stable h
# 6189    Zfp445  81 stable h
# 6190       Tef  91 stable h
# 6191    Zfp664  94 stable h
# 6192    Zfp275 109 stable high


sgenescorer2_matrix["Yy1",] #supercool
sgenescorer2_matrix["Yy2",] #gata3
sgenescorer2_matrix["Klf13",] #coolish irf4
sgenescorer2_matrix["Pias1",] #has potential for il13
sgenescorer2_matrix["Zkscan1",] #crap
#Zfp275 Zfp664 Tef Zfp445 Setdb2 Nfatc2 Zkscan1  uncool 

t.test(
  nn$c[nn$Regulator %in% curatedpos$mgi_symbol],
  nn$c)  #pos a bit higher than average. not significant


mean(curatedpos$mgi_symbol %in% nn$Regulator)
mean(curatedneg$mgi_symbol %in% nn$Regulator) #excluded lowly expressed genes, not fair
#nn$
  

## Compare overall #neighbours, screen score
compn2s <- mergebyrow(
  data.frame(row.names = rownames(sgenescorer2_matrix), ss=apply(sgenescorer2_matrix, 1, min)),
  nn  
)    
plot(compn2s$ss, compn2s$c)
cor.test(compn2s$ss, compn2s$c, method = "spearman") #p=0.05, rho=0.09. unfortunately -> worse screen hit <-> more neighbours


##For Gata3, compare MI vs screen score

plotMIvsS <- function(gene){
  v<-checkrelated(aracne_ko, gene)
  w<-data.frame(mgi_symbol=rownames(sgenescorer2_matrix), ss=sgenescorer2_matrix[,gene])
  v<-smerge(v,w, all.y = TRUE)
  v$MI[is.na(v$MI)] <- 0
  plot(v$MI, v$ss)  
  print(cor.test(v$MI,v$ss,method = "spearman"))
}
plotMIvsS("Gata3") #12% pval, right sign
plotMIvsS("Irf4") # right sign
plotMIvsS("Xbp1") #17% pval, wrong sign
plotMIvsS("Il13") #wrong sign
#plot(density(v$))
#Il4 does not work


###########################
# so which TF does have most of the hits? a TF + downstream genes is called a "regulon" in K wording. and it is just
# a fisher test. test all TFs. generalize the GO GSEA code?








###############################################################################################
###############################################################################################
#### Do stuff ######################################################
###############################################################################################
###############################################################################################




doaracneDPI(aracne_tc_atac_de)  #kills 2%
  nrow(aracne_tc_atac_de)
doaracneDPI(aracne_tc_atac)  #kills 10% interactions
write_aracne_net(doaracneDPI(aracne_tc_atac), "aracne/out_tc_dpi.txt")
write_aracne_net(aracne_tc_atac, "aracne/out_tc_nodpi.txt")

# removearacnesingle <- function(net){
#   cnt <- countaracneneigh(net)
#   k <- names(which(cnt>1))
#   net <- net[net$Regulator %in% k & net$Target %in% k,]
#   print(dim(net))
#   net
# }

aracne_ko_atac <- merge_aracneATACchip(aracne_ko)
#keepgene <- c("Irf4","Xbp1","Gata3","Il4","Il13","Batf","Batf3")
#aracne_ko_atac_red <- aracne_ko_atac[(aracne_ko_atac$pvalue < 1e-15 & aracne_ko_atac$MI>0.15) | aracne_ko_atac$Regulator %in% keepgene | aracne_ko_atac$Target %in% keepgene,] 
aracne_ko_atac_red <- cutneighN(reducearacneDE(aracne_ko_atac),n=4)
# #aracne_ko_atac_red <- aracne_ko_atac[aracne_ko_atac$pvalue < 1e-15 & aracne_ko_atac$MI>0.15,] 
# nrow(aracne_ko)
# nrow(aracne_ko_atac)
# nrow(aracne_ko_atac_red)
#filterAracneATAC
write_aracne_net(aracne_ko_atac_red, "aracne/merge_koATAC.txt")

# "Irf4" %in% sort(unique(c(aracne_ko_atac$Regulator, aracne_ko_atac$Target)))
# "Batf" %in% sort(unique(c(aracne_ko_atac$Regulator, aracne_ko_atac$Target)))
# "Batf3" %in% sort(unique(c(aracne_ko_atac$Regulator, aracne_ko_atac$Target)))
# "Ern1" %in% sort(unique(c(aracne_ko_atac$Regulator, aracne_ko_atac$Target)))
#Lost Ern1


# sort(unique(c(aracne_ko_atac$Regulator, aracne_ko_atac$Target)))
# sort(unique(c(aracne_ko_atac_red$Regulator, aracne_ko_atac_red$Target)))

######################
######################   Save merged networks
######################

#TODO: ensure all genes at least connected by N edges

### Mouse TC
aracne_tc_atac <- merge_aracneATACchip(aracne_tc)
aracne_tc_atac_red <- filterAracneATAC(aracne_tc_atac, 1e-15, 0) 
nrow(aracne_tc)
nrow(aracne_tc_atac)
nrow(aracne_tc_atac_red)
write_aracne_net(aracne_tc_atac_red, "aracne/merge_tcATAC.txt")

### Human TC
aracne_humantc_atac <- merge_aracneATACchip(aracne_humantc, mapMotifGene_ = hg_mapMotifGene, ensconvert_ = human_ensconvert)
aracne_humantc_atac_red <- filterAracneATAC(aracne_humantc_atac, 1e-15, 0) 
write_aracne_net(aracne_humantc_atac_red, "aracne/merge_humantc.txt")

### Mouse and human TC merged
aracne_mergetc_atac <- intersectARACNe(aracne_humantc_atac,aracne_tc_atac)
write_aracne_net(aracne_mergetc_atac, "aracne/merge_mergetc.txt")

cutneighN(aracne_mergetc_atac)

aracne_mergetc_atac_red <- cutneighN(aracne_mergetc_atac, n=6)
write_aracne_net(aracne_mergetc_atac_red, "aracne/merge_mergetc_red.txt")

### Another way of checking #partner vs score
v <- table(c(aracne_mergetc_atac$Regulator,aracne_mergetc_atac$Target))
v<-data.frame(mgi_symbol=names(v), nn=as.double(v), ss=sgenescorer2_matrix[names(v),"Il4"])
v<-na.omit(v)
plot(v$ss, log10(v$nn)+runif(nrow(v))*0.1)
cor.test(v$ss, v$nn, method = "spearman")  #-4.5% correlation for Gata3. +4.5 for Ir4. -1.3% for xbp1
#sgenescorer2_matrix[names(v),"Il4"]





##############################
### Only keep TFs which are some kind of hit
reducearacnehit <- function(gene,n=1500,out=gene,net=aracne_mergetc_atac){
  keeptf <- names(which(apply(sgenescorer2_matrix[,gene,drop=FALSE]<n,1,any)))
  net <- net[net$Regulator %in% keeptf & net$Target %in% keeptf,]
  print(dim(net))
  write_aracne_net(net, sprintf("aracne/merge_mergetc_hit_%s.txt",out))
}
reducearacnehit("Il4",n=500)
reducearacnehit("Il13")
reducearacnehit("Irf4")
reducearacnehit("Xbp1")
reducearacnehit("Gata3")
reducearacnehit("Il4",n=500,out="mouseIl4",net=aracne_tc_atac)

#### DE
aracne_tc_atac_de <- reducearacneDE(aracne_tc_atac, qval = 1e-3)
write_aracne_net(aracne_tc_atac_de, "aracne/out_mouseDE.txt")

aracne_kotc_atac_de <- reducearacneDE(mergeARACne(aracne_tc_atac, aracne_ko_atac), qval = 1e-3)
write_aracne_net(aracne_kotc_atac_de, "aracne/out_mouse_kotcDE.txt")


aracne_mh_kotc_atac_de <- doaracneDPI(reducearacneDE(mergeARACne(aracne_mergetc_atac, aracne_ko_atac), qval = 1e-3))
#write_aracne_net(cutneighN(aracne_mh_kotc_atac_de), "aracne/out_mousehuman_kotcDE.txt")
write_aracne_net(aracne_mh_kotc_atac_de, "aracne/out_mousehuman_kotcDE.txt")

######## TODO: For TFs which there are no motifs, could instead rely on MI and pick top 10%?


un <- unique(c(aracne_mh_kotc_atac_de$Regulator,aracne_mh_kotc_atac_de$Target))
intersect(un,rownames(sgenescorer2_matrix)[sgenescorer2_matrix[,"Il4"]<300])
intersect(un,rownames(sgenescorer2_matrix)[sgenescorer2_matrix[,"Il13"]<500])
intersect(un,rownames(sgenescorer2_matrix)[sgenescorer2_matrix[,"Gata3"]<500])
intersect(un,rownames(sgenescorer2_matrix)[sgenescorer2_matrix[,"Xbp1"]<500])
intersect(un,rownames(sgenescorer2_matrix)[sgenescorer2_matrix[,"Irf4"]<500])
#### good number!!

#### Could one force our KOed genes, IL4,13, to be included?
# Bcl11b has no motif
# aracne_ko[aracne_ko$Target=="Bcl11b" | aracne_ko$Regulator=="Bcl11b",]   #Many connections but no motif. can we infer?
# aracne_ko[aracne_ko$Target=="Il13" | aracne_ko$Regulator=="Il13",]
# "Il13" %in% c(aracne_ko_atac$Target, aracne_ko_atac$Regulator)  #ATAC kicks IL13 out currently, but is in net. IL4 not in network. change MI cutoff?

write.table(
  as.data.frame(countaracneneigh(aracne_mh_kotc_atac_de)), 
  "aracne/out_mousehuman_kotcDE.nodes.txt",sep="\t", quote = FALSE, row.names = FALSE)

# entire paper on Bcl11b in T cells. http://www.pnas.org/content/114/23/5800.full  
#### but makes sense to include it in the network?

write_aracne_net(doaracneDPI(aracne_tc_atac_de), "aracne/out_mouseDE_dpi.txt")




#Dph6 in Th*. effector T cells. reg T cells



hist(log(as.double(1+aracne_ko_atac$pvalue)),breaks=100)

aracne_ko_atac[aracne_ko_atac$Regulator=="Gata3",]
aracne_ko_atac[aracne_ko_atac$Regulator=="Tbx21",]
aracne_ko_atac[aracne_ko_atac$Regulator=="Irf4",]
aracne_ko_atac[aracne_ko_atac$Regulator=="Batf",]
aracne_ko_atac[aracne_ko_atac$Regulator=="Yy1",]

aracne_tc_atac[aracne_tc_atac$Regulator=="Gata3",]
aracne_tc_atac[aracne_tc_atac$Regulator=="Tbx21",]
aracne_tc_atac[aracne_tc_atac$Regulator=="Irf4",]
aracne_tc_atac[aracne_tc_atac$Regulator=="Batf",]
aracne_tc_atac[aracne_tc_atac$Regulator=="Yy1",]

aracne_tc_atac[aracne_tc_atac$Regulator=="Stat6",]



############
### Combine 3 KO ARACNE networks. By running them separately, can avoid effects due to different mouse backgrounds
##### 

# Average of MI, log-average of p-value? 
# If not found in one, MI=0? intersect them all, min & max?

aracne_ko1 <- read.aracne_net("aracne/ko1.out/network.txt")
aracne_ko2 <- read.aracne_net("aracne/ko2.out/network.txt")
aracne_ko3 <- read.aracne_net("aracne/ko3.out/network.txt")
nrow(aracne_ko1)
nrow(aracne_ko2)
nrow(aracne_ko3)  #Plate #3 also suggested to be the most informative!!!

aracne_ko123 <- rbind(
 read.aracne_net("aracne/ko1.out/network.txt"),
 read.aracne_net("aracne/ko2.out/network.txt"),
  read.aracne_net("aracne/ko3.out/network.txt")
)
aracne_ko123c <- sqldf("select distinct `Regulator`,`Target`, max(`MI`), min(pvalue) from aracne_ko123 group by `Regulator`,`Target` ")
#aracne_ko123c <- sqldf("select *, count(*) as c from aracne_ko123 group by `Regulator`,`Target` ")
#hist(aracne_ko123c$c)









#Sanger talk on ALL. Il4. Utx-/Y males do not develop AML
#Utx+Uty same pheno as homoz phemales
#
# Utx- -> ETS program upreg
# siEWSR1-FLI1 phusion. 
# they got Tbx21 chipseq in Utx+/- progenitor cells
# Utx- => supression of gata program (all)
#interacts with Brg1 and Chd4
sgenescorer2_matrix["Kdm6a",] #=Utx
sgenescorer2_matrix["Smarca4",]
sgenescorer2_matrix["Chd4",]


#sgenescorer2_matrix["Bcl11b",]   #comes fairly high in all screens :o






### 40% of the edges are covered by ATAC. the other 60, minimize
mean(aracne_ko$Regulator %in% map_jaspar_namegenesym$mgi_symbol | aracne_ko$Target %in% map_jaspar_namegenesym$mgi_symbol)
mean(aracne_tc$Regulator %in% map_jaspar_namegenesym$mgi_symbol | aracne_tc$Target %in% map_jaspar_namegenesym$mgi_symbol)

#Add edge info: ATAC-validated

#50kb cutoff: 70k edges in atac_ko

#7712/(76765*0.39)  #ATAC keeps 25% of affected edges with 10kb cutoff



sum(aracne_ko[aracne_ko$Regulator %in% map_jaspar_namegenesym$mgi_symbol | aracne_ko$Target %in% map_jaspar_namegenesym$mgi_symbol,]$MI>0.1)



###############################################################################################
###############################################################################################
#### serious stuff ######################################################
###############################################################################################
###############################################################################################


### would be good to keep track of evidence
### alternative: merge everything with evidence. no DPI. then intersect when making final plot. but: huge!
# aracne_ko_atac_dpi <- doaracneDPI(aracne_ko_atac)
# aracne_ko_atac$isdirect<-TRUE
# aracne_ko_nontf <- cutNonAtac(aracne_ko)
# aracne_ko_nontf$isdirect<-FALSE
# test <- rbind(aracne_ko_atac,aracne_ko_nontf)



### TODO: remove redundant edges should prioritize directed over non-directed


### Mouse TC
aracne_mousetc_atac    <- merge_aracneATACchip(aracne_tc)
aracne_mousetc_nonatac <- cutNonAtac(aracne_tc)
aracne_mousetc_union <- directARACNe(removerededge(rbind(
  aracne_mousetc_atac,
  aracne_mousetc_nonatac)))

aracne2gml("aracne/tc_mouse_atac.gml", aracne_mousetc_union)

### Human TCz
aracne_humantc_atac    <- merge_aracneATACchip(aracne_humantc)
aracne_humantc_nonatac <- cutNonAtac(aracne_humantc)
aracne_humantc_union <- directARACNe(removerededge(rbind(
  aracne_humantc_atac,
  aracne_humantc_nonatac)))

aracne2gml("aracne/tc_human_atac.gml", aracne_humantc_union)

### Conserved TC human and mouse
intersectalltc <- intersectARACNe(
  aracne_mousetc_union,
  aracne_humantc_union
  )

### KO for mouse
aracne_mouseko_atac    <- merge_aracneATACchip(aracne_ko)
aracne_mouseko_nonatac <- cutNonAtac(aracne_ko)
aracne_mouseko_union   <- unionARACNE(
  aracne_mouseko_atac,
  aracne_mouseko_nonatac)

### KO and TC together
aracne_mergekotc   <- unionARACNE(
  intersectalltc,
  aracne_mouseko_union)


aracne_mergekotcde <- reducearacneDE(aracne_mergekotc, qval = 1e-10, thecut = 1e-10)

aracne_mergekotcdeInt <- reducearacneinteresting(aracne_mergekotcde)

aracne2gml("aracne/out_aracne_mergekotcdeInt.gml",aracne_mergekotcdeInt)


#####################################################################
#####################################################################
#####################################################################
#####################################################################

#test3 <- aracne_mergekotcde[aracne_mergekotcde$Regulator %in% interesting & test2$Target %in% interesting, ]




test2 <- directARACNe(removerededge(reducearacneDE(test, qval = 1e-10, thecut = 1e-10)))
print(nrow(test2))
#write_aracne_net(test2,"aracne/test.txt")
interesting <- c(togenesym2(list_tf),
                 ensconvert$mgi_symbol[startsWith(ensconvert$mgi_symbol,"Il")],
                 ensconvert$mgi_symbol[startsWith(ensconvert$mgi_symbol,"Cxc")],
                 ensconvert$mgi_symbol[startsWith(ensconvert$mgi_symbol,"Ccr")])
test3 <- test2[test2$Regulator %in% interesting & test2$Target %in% interesting, ]
nrow(test3)
aracne2gml("aracne/test.gml", test3)


test3[test3$Target=="Hic1",]


"Il13" %in% test$Regulator
"Il13" %in% test2$Regulator



list_tf

#"Il13" %in% 
de_late[de_late$mgi_symbol=="Il13",]

################# story: non-DE hits are *maintenance* - but not part of the program. a plot of DE genes vs hits? we got it already!
### bcl11b is not DE but impact, known. clearly maintenance.
#story again: what IS the characteristics of the hits?
  



#From the Viper paper, one could consider these as TFs:
# GO:0003700, 'transcription factor activity', 
# or as GO:0004677, 'DNA binding',
# and GO:0030528, 'transcription regulator activity', or as
# GO:0004677 and GO: 
#   0045449, 'regulation of transcription'), 969 transcriptional cofactors (a manually curated list, not overlapping with the transcription factor list, built upon genes annotated as GO:0003712, 'transcription cofactor activity', 
# or GO:0030528 or GO:0045449) or 3,370 signaling pathway related genes (annotated in GO Biological Process database as 
# GO:0007165 'signal transduction' and in GO cellular component database as GO:0005622, 'intracellular',
# or GO:0005886, 'plasma membrane'



############## For each regulon, check if the genes are in by a KS-test.
### 

aracne_ko_gsea
aracne_ko_gsea[aracne_ko_gsea$mgi_symbol=="Gata3",]

temp <-smerge(aracne_ko_gsea, data.frame(mgi_symbol=rownames(sgenescorer2_matrix),ss=sgenescorer2_matrix[,"Xbp1"]))
temp[order(temp$pt),]
#Zfp526, Arid5b, Rest, Wiz, Cdc5l: effect on Xbp1 and downstream also effect. arid5b maybe DE
#Stat2, Zfp422, Smad3, Pias1: effect on Il13 and downstream also effect. stat2 zfp422 bit DE, smad3 DE
#Zfp422, Thap1, Cic, Atf3, (Pou2f1): Gata3 ...   Zfp422, Cic a bit DE. Atf3 very DE.  
#Zfp526, (Zfp422), Arid3b, Foxo3  : Il4     Zfp422 a bit DE, maybe foxo3
#Zfp408 ! : Irf4.    only one. but decent expression. fairly stable   

aracne_mousetc_gsea
aracne_mousetc_gsea[aracne_mousetc_gsea$mgi_symbol=="Gata3",]
#Stat2 has a direct effect on Il13 as well as enriched genes beneath
#Zfp526 effect on Il4 and Xbp1

#### which regulons most overlap. do they correlate the same way
#### the TCGA network => network (this is cancer data)


sgenescorer2_matrix[c("Stat1","Stat2","Zfp512","Hivep1","Zfp526","Tox4"),]

aracne_ko_gsea
aracne_ko_gsea[aracne_ko_gsea$mgi_symbol=="Gata3",]


head(aracne_ko_gsea,     n=20)   #after dealing with the NA, gata3 went down. also, need to use lateDE list
head(aracne_mousetc_gsea,n=20)
head(aracne_humantc_gsea,n=20)
head(aracne_rawitc_gsea, n=20)



aracne_rawittc <- intersectARACNe(
  aracne_mousetc,
  aracne_humantc
)
aracne_kotcrawint <- intersectARACNe(
  aracne_mousetc,
  aracne_ko
)

aracne_humanrawint <- intersectARACNe(
  aracne_humantc,
  aracne_thymoma
)


aracne2gml("aracne/tc_rawintersect.gml",aracne_rawittc)
aracne2gml("aracne/kotc_rawintersect.gml",aracne_kotcrawint)
aracne2gml("aracne/human_rawintersect.gml",aracne_humanrawint)

#delist <- unique(names(which(getDElistfromBothEarlyLate()< -2)))
aracne_humanrawint_de <- reducearacneDE(aracne_humanrawint, qval = 1e-12)
nrow(aracne_humanrawint_de)
aracne2gml("aracne/human_rawintersect_de.gml",aracne_humanrawint_de)


###################

aracne_humanrawint_red <- intersectARACNe(
  aracne_humantc[aracne_humantc$pvalue<1e-30,],
  aracne_thymoma[aracne_thymoma$pvalue<1e-30,]
)
aracne_humanrawint_red_de <- reducearacneDE(aracne_humanrawint_red, qval = 1e-12)
nrow(aracne_humanrawint_red_de)
aracne_humanrawint_red_de <- merge_aracneATACchip(aracne_humanrawint_red_de)
#aracne_humanrawint_red_de <- cutneighN(aracne_humanrawint_red_de, n=20)
aracne2gml("aracne/aracne_humanrawint_de_red.gml",aracne_humanrawint_red_de)

#aracne_humanrawint_de_red <- aracne_humanrawint_de[aracne_humanrawint_de$pvalue<1e-20,]
# nrow(aracne_humanrawint_de_red)
# aracne_humanrawint_de_red <- cutneighN(aracne_humanrawint_de_red,n=2)
# nrow(aracne_humanrawint_de_red)
# #hist(aracne_humanrawint_de$MI)
# #nrow(aracne_humanrawint_de_red)
# aracne2gml("aracne/aracne_humanrawint_de_red.gml",aracne_humanrawint_de_red)

aracne_interesting <- list_tf   #maybe everything with a hit too?

ensc
#grep("Il",)





#### NEW: correlation test between all genes

log_mtpm <- log10(1+mtpm)
timestamp()
foo <- for(i in 1:10000){#nrow(mtpm)){
  cor(as.double(log_mtpm[1,]),(as.double(log_mtpm[i,])))
}
timestamp()




dim(mtpm)


#council tax ref 82288

###############################################################################################
#### Calculate correlations between all the genes  ############################################
###############################################################################################

#13407 genes in mouse

13407*2000*8/1e6  #200mb per correlation file

24*2000/3600  #13h to test all genes




#matrix2genesym(mtpm, ensconvert, removeclashing = TRUE)


saveTPMforCorr <- function(docor, x, ensconvert_, mintpm=2){
  nexpressedGenesId <- rownames(x)[apply(x>mintpm,1,any)]  #NOTE: log. not log10
  print(dim(x))
  x <- x[rownames(x) %in% nexpressedGenesId,]
  print(dim(x))
  x <- matrix2genesym(x, ensconvert_, removeclashing = TRUE)
  print(dim(x))
  #x <- log10(1+x)
  x <- t(x)
  m <- x[order(rownames(x)),]
  save(m, file=sprintf("%s/m.RData",docor))
}

saveTPMforCorr("aracne/corr/mousetc",mtpm,       ensconvert,       mintpm=5)
saveTPMforCorr("aracne/corr/humantc",human_mtpm, human_ensconvert, mintpm=5)
saveTPMforCorr("aracne/corr/ko",     ncount,     ensconvert,       mintpm=50)



################################

#docor <- "testcor"
#saveTPMforCorr(docor,mtpm)


runCorrBootstrap <- function(docor, tcrand=FALSE){
  load(sprintf("%s/m.RData",docor))
  numboot <- 60
  for(i in 1:numboot){
    print(i)
    if(tcrand){
      mboot <- m
      for(j in 1:ncol(m)){
        mboot[,j] <- m[c(
          sample(1:3),sample(4:6),sample(7:9),sample(10:12),sample(13:15),sample(16:18),
          sample(19:21),sample(22:24),sample(25:27),sample(28:30),sample(31:33),sample(34:36),
          sample(37:39),sample(40:42),sample(43:45),sample(46:48),sample(49:51),sample(52:54),sample(55:57)),j]    
      }
      oc<-cor(mboot,method="spearman")
      save(oc, file=sprintf("%s/oc_%s.RData",docor,unclass(Sys.time())))
    } else {
      #Jack-knife
      j<-round(runif(1,min=1,max=ncol(m)))
      fn<-sprintf("%s/oc_%s.RData",docor,j)
      if(!file.exists(fn)){
        mboot <- m[,-j]
        oc<-cor(mboot,method="spearman")
        if(!file.exists(fn)){
          save(oc, file=fn)
        }
      }
    }
  }
}

############################

library(reshape)
library(stringr)

doboot <- list.files(docor)
doboot <- doboot[grep("oc",doboot)]

load(sprintf("%s/m.RData",docor))

#The centered value
meancor2<-cor(m,method="spearman")

#Get average
load(sprintf("%s/%s",docor,doboot[1]))
meancor <- oc*0
i<-1
for(f in doboot){
  print(i)
  load(sprintf("%s/%s",docor,f))
  meancor <- meancor + oc
  i<-i+1
}
meancor <- meancor / length(doboot)

#Get variance
varcor <- oc*0
i<-1
for(f in doboot){
  print(i)
  load(sprintf("%s/%s",docor,f))
  varcor <- varcor + (oc-meancor)^2
  i<-i+1
}
varcor <- varcor / length(doboot)

#Calculate Zscore
zcor <- meancor2/sqrt(varcor)
#zcor <- (meancor2-mean(as.double(meancor2)))/sqrt(varcor)
qnorm(0.99)
sigcor <- abs(zcor)>qnorm(1-1e-10)  & abs(meancor2)>0.7
mean(as.double( sigcor ),na.rm = TRUE)

hist(as.double(zcor),breaks=50)  
#Extract relevant correlations

out <- cbind(melt(meancor2),melt(zcor)[,3])
colnames(out) <- c("Regulator","Target","corr","pvalue")
fout <- out[out$pvalue>qnorm(1-1e-10) & abs(out$corr)>0.8 & out$Regulator!=out$Target,]
write.table(fout, sprintf("%s/network.csv",docor),quote = FALSE)



#### No p-values
rownames(meancor2) <- normalizesym(rownames(meancor2))
colnames(meancor2) <- normalizesym(colnames(meancor2))
out <- cbind(melt(meancor2))
colnames(out) <- c("Regulator","Target","corr")



fout <- out[abs(out$corr)>0.8 & out$Regulator!=out$Target,]   #For the KO -> 30M?
nrow(fout)

quantile(abs(out$corr),0.9)

fout <- out[abs(out$corr)>0.3 & out$Regulator!=out$Target,]   #For the mouse TC  -> 10M edges? human TC -> 3M edges
nrow(fout)


write.table(fout, sprintf("%s/network.csv",docor),quote = FALSE,sep="\t",row.names = FALSE)


save.image()



#nexpressedGenes   <- unique(ensconvert$mgi_symbol[ensconvert$ensembl_gene_id %in% expressedGenesId])



hist(as.double(log(1+m[20,])),breaks=40)

#############################




library(reshape)
library(stringr)

calccorrM <- function(docor){
  load(sprintf("%s/m.RData",docor))
  
  meancor2<-cor(m,method="spearman")
  
  #### No p-values
  rownames(meancor2) <- normalizesym(rownames(meancor2))
  colnames(meancor2) <- normalizesym(colnames(meancor2))
  out <- cbind(melt(meancor2))
  colnames(out) <- c("Regulator","Target","corr")
  
  fout <- out[abs(out$corr)>quantile(abs(out$corr),0.95) & out$Regulator!=out$Target,]   
  print(nrow(fout))
  write.table(fout, sprintf("%s/network.csv",docor),quote = FALSE,sep="\t",row.names = FALSE)
  
  #TODO: keep only correlations TF/cyto - any gene
}

calccorrM("mousetc")
calccorrM("humantc")
calccorrM("ko")



###############################################################################################
#### Check MI network using corr  ############################################################
###############################################################################################

library(reshape2)
temp <- aracne_mousetc[,1:3]  #set MI to 1?
temp <- acast(temp, Regulator~Target, value.var="MI", fill=0)
thec <- cor(t(temp),method="spearman")
for(i in 1:nrow(thec))
  thec[i,i] <- 0

#For KO
keep <- apply(thec>quantile(thec,0.998) | 
                thec<quantile(thec,0.002) | 
                rownames(thec) %in% c("Gata3","Batf","Il4","Il13","Ir4","Xbp1","Tbx21","Fli1"),1,any)
#For TC
keep <- apply(thec>quantile(thec,0.9997) | 
                thec<quantile(thec,0.0003) | 
                rownames(thec) %in% c("Gata3","Batf","Il4","Il13","Ir4","Xbp1","Tbx21","Fli1","Yy1","Yy2"),1,any)

##Yy1 / 2 not a regulator? Yy1 is

thec2 <- thec[keep, keep]
heatmap.2(
  trace="none",
  thec2,
  col=colorRampPalette(c("blue", "white", "red"))(n = 300)
  )




###############################################################################################
#### Check corr network using corr  ############################################################
###############################################################################################


read.corr_net <- function(file){
  x <- read.csv(file,sep="\t",stringsAsFactors = FALSE)
  #x <- x[order(x$MI, decreasing = TRUE),]
  #x <- x[order(x$pvalue),]
  x
}


aracne_mousetc <- read.corr_net("aracne/corr/mousetc.csv")
aracne_humantc <- read.corr_net("aracne/corr/humantc.csv")
aracne_ko <- read.corr_net("aracne/corr/ko.csv")



corrinteresting <- c(togenesym2(list_tf), "Il4","Il13") 

library(reshape2)
temp <- aracne_mousetc[,1:3]  
#temp <- aracne_ko[,1:3]  
temp <- temp[temp$Regulator %in% corrinteresting | temp$Target %in% corrinteresting,]
temp <- acast(temp, Regulator~Target, value.var="corr", fill=0)
qval <- 1e-2
gde <- unique(c(de_early$mgi_symbol[de_early$qval<qval],de_late$mgi_symbol[de_late$qval<qval]))
temp <- temp[rownames(temp) %in% gde,
             colnames(temp) %in% corrinteresting]


       

#dim(temp)
thec <- cor(temp, method="spearman")  #t needed?
for(i in 1:nrow(thec))
  thec[i,i] <- 0
thec[is.na(thec)]<-0 #not sure why this is needed with few genes

#For TC
isexp <- rownames(thec) %in% togenesym2(rownames(mtpm)[apply(mtpm>10,1,any)])
keep <- isexp & (apply(thec < -0.4 | thec > 0.4,1,any) | 
                   rownames(thec) %in% c("Gata3","Batf","Il4","Il13","Ir4","Xbp1","Tbx21","Fli1"))
# keep <- isexp & (apply(thec < -0.9 | thec > 0.9,1,any) | 
#                 rownames(thec) %in% c("Gata3","Batf","Il4","Il13","Ir4","Xbp1","Tbx21","Fli1"))
sum(keep)

#For KO
isexp <- rownames(thec) %in% togenesym2(rownames(mtpm)[apply(mtpm>20,1,any)])
keep <- isexp & (apply(thec < -0.8 | thec > 0.8,1,any) | 
                   rownames(thec) %in% c("Gata3","Batf","Il4","Il13","Ir4","Xbp1","Tbx21","Fli1"))
sum(keep)

##Yy1 / 2 not a regulator? Yy1 is
thec2 <- thec[keep, keep]
plot.new()
heatmap.2(
  trace="none",
  thec2,
  col=colorRampPalette(c("blue", "white", "red"))(n = 300)
)










########################## Actually used ###############################


aracne_ko13_int_red <- merge_aracneATACchip(aracne_ko13_int, atac.mouse$mapPeakGene, ensconvert)
aracne_ko13_int_red
#aracne_ko13_int$pvalue
aracne2gml(aracne_ko13_int, file = "aracne/aracne_ko13_int.gml")
aracne2gml(aracne_ko13_int_red, file = "aracne/aracne_ko13_int_red.gml")
#aracne_ko13_int_red

aracne_ko123_int <- intersectARACNe(aracne_ko13_int, aracne_ko2) 
#aracne2gml(aracne_ko123_int, file = "aracne/aracne_ko123_int.hitsonly.gml",keephitonly = TRUE)

aracne_ko123c_red <- filterAracneATAC(aracne_ko123c, 1e-15, MI=0.39)
#write_aracne_net(aracne_ko123c_red, "aracne/merge_ko123.txt")
########################################################################




aracne_difftc_int <- intersectARACNe(aracne_difftc.mouse, aracne_difftc.human)
aracne2gml(merge_aracneATACchip(aracne_difftc.mouse, atac.mouse$mapPeakGene, ensconvert),
           file = "aracne/aracne_difftc.gml",keephitonly = TRUE)

aracne2gml(intersectARACNe(aracne_tc, aracne_humantc),
           file = "aracne/a.gml",keephitonly = TRUE)

### in this network, blue and green seem to seggregate
totest <- intersectARACNe(aracne_tc, aracne_humantc)







aracne2gml(merge_aracneATACchip(aracne_tc, atac.mouse$mapPeakGene, ensconvert),
           file = "aracne/test.gml",keephitonly = TRUE)


#Any point checking DE? built from Th2 after all
nummention <- sort(table(c(totest$Regulator,totest$Target)))
nummention

#cutting neighbours twice should be enough to remove all boring nodes
v <- undirectARACNe(cutneighN(cutneighN(filterAracneHit(totest,c("Il4","Il13","Irf4","Gata3")))))  #Il13
w <- sqldf("select distinct `Regulator`,count(`Target`) as c from v group by `Regulator` order by c desc")
head(w,n=10)
sgenescorer2_matrix["Mbd1",]


edges <- filterAracneHit(edges)
edges <- cutneighN(edges)





