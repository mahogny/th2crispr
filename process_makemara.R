#############################################################################################
###                                                                                       ###
###             Part of the paper ...                                                     ###
###             Author: Johan Henriksson (mahogny@areta.org)                              ###
###                                                                                       ###
###             This code ...                                   ###
###                                                                                       ###
#############################################################################################

library(reshape2)
library(Rtsne)
library(gplots)
library(RColorBrewer)
library(stringr)
library(sqldf)

# source("common_geneid.R")
# source("common_gotest.R")
# source("common_readrnaatacchip.R")




################################################################################
########## Generate MARA input files for time course ###########################
################################################################################


## Function to write the two MARA input files
writemara <- function(ensconvert_,  genetfcount,dir="mara_redde",mtpm,de_genes,minsites=20){
  
  #mtpm <- matrix2genesym(mtpm, ensconvert_,removeclashing=TRUE)
#  print(head(mtpm))
  d <- getmarasitecountmatrix(genetfcount)

  fdir <- sprintf("out_mara/%s",dir)
  if(!file.exists(fdir))
    dir.create(fdir)
    
  #Genes for which we will fit TFs. The annoying "tf" row disappears here
  #mara_genes <- intersect(rownames(d),rownames(mtpm)) #
  mara_genes <- intersect(intersect(rownames(d),rownames(mtpm)),de_genes)
#  print(head(d))
 # print(mara_genes)
#  print(mara_genes)
  
  # print(mara_genes)  
  # print(tail(mara_genes))
  #print(mara_genes[2000:2100])
  #mara_genes <- mara_genes[1:2050]
  
  #Require sufficient binding sites to include a gene
  # mara_genes <- mara_genes[apply(d[mara_genes,],1,function(x) sum(as.double(x)))>minsites]
  #print(mara_genes)
  
  dsite <- d[mara_genes,]
  dsignal <- log(1+mtpm[mara_genes,])
#  print(apply(dsite,1,function(x) sum(as.double(x)))[2000:2100])
#  print(apply(dsite,1,function(x) sum(as.double(x))))
 
  #Any TFs with too few sites to be fit?
  d <- d[,apply(dsite,2,function(x) sum(as.double(x)))>minsites]
#  print(apply(dsite,2,function(x) sum(as.double(x))))
  
  # #unknowns: #samples * #motifs.
  print(sprintf("# unknowns %s",ncol(mtpm)*ncol(d)))
  # #equations: #genes * #samples
  print(sprintf("# equations %s",ncol(mtpm)*length(mara_genes)))
  
  
  # print(dim(dsite))
  # print(dim(dsignal))
#  print(all(apply(dsite,2,function(x) sum(as.double(x)))>0))
  

  #Check for colinearity  
  # Nc <- cor(dsite)
  # for(i in 1:nrow(Nc))
  #   Nc[i,i]<-0
  # print(tail(sort(apply(Nc,1,max)),n=40))

  #TODO very aggressive!
  #dsite <- dsite[,apply(Nc,1,max)<0.95 ]    
  
  write.table(dsite,sprintf("out_mara/%s/mara_sitecount.csv",dir),sep="\t",quote = FALSE)
  write.table(dsignal,sprintf("out_mara/%s/mara_signal.csv",dir),sep="\t",quote = FALSE)
}

mtpm2 <- log(1+matrix2genesym(mtpm, ensconvert,removeclashing=TRUE))
mtpm2["Stat6",]

#log(1+mtpm["ENSMUSG00000002147",])

#allde_conserved_anytime


#tcmouse$mtpm


##################################################
########### finally ????? ########################
##################################################

tcmouse$all_genes <- rownames(tcmouse$mtpm)
tchuman$all_genes <- rownames(tchuman$mtpm)

writemara(ensconvert,       atac.mouse$cons_tfc,
          dir="mouse_alltf_cDEgenes_alltime_th20",tcmouse$mtpm,allde$ens_mouse[allde$conserved_loose])
writemara(human_ensconvert, atac.human$cons_tfc,
          dir="human_alltf_cDEgenes_alltime_th20",tchuman$mtpm,allde$ens_human[allde$conserved_loose])

writemara(ensconvert,       atac.mouse$cons_tfc,
          dir="mouse_alltf_allgenes_alltime_th20",tcmouse$mtpm,tcmouse$all_genes)
writemara(human_ensconvert, atac.human$cons_tfc,
          dir="human_alltf_allgenes_alltime_th20",tchuman$mtpm,tchuman$all_genes)

writemara(ensconvert,       atac.mouse$noncons_tfc,
          dir="mouse_alltf_allgenes_alltime_th20_nonconsP",tcmouse$mtpm,tcmouse$all_genes)
writemara(human_ensconvert, atac.human$noncons_tfc,
          dir="human_alltf_allgenes_alltime_th20_nonconsP",tchuman$mtpm,tchuman$all_genes)

writemara(ensconvert,       atac.mouse$noncons_tfc,
          dir="mouse_alltf_cDEgenes_alltime_th20_nonconsP",tcmouse$mtpm,allde$ens_mouse[allde$conserved_loose])
writemara(human_ensconvert, atac.human$noncons_tfc,
          dir="human_alltf_cDEgenes_alltime_th20_nonconsP",tchuman$mtpm,allde$ens_human[allde$conserved_loose])


writekomara <- function(){
  ncounts_corr <- normCounts(set2)
  red <- ncounts_corr
  colnames(red) <- sprintf("%s_%s",1:ncol(red),colnames(red))
  red_allgenes <- rownames(red)
  red_de <- allde$ens_mouse[allde$conserved_loose]
  writemara(ensconvert,     atac.mouse$noncons_tfc,
            dir="ko_alltf_allgenes_nonconsP",red,red_allgenes)  
  writemara(ensconvert,     atac.mouse$noncons_tfc,
            dir="ko_alltf_cDEgenes_nonconsP",red,red_de)  
  writemara(ensconvert,     atac.mouse$cons_tfc,
            dir="ko_alltf_allgenes",red,red_allgenes)  
  writemara(ensconvert,     atac.mouse$cons_tfc,
            dir="ko_alltf_cDEgenes",red,red_de)  
}
writekomara()



################################################################################
########################### overlaps peaks vs DE ###############################
################################################################################

#should include in supplementary

qtestfish <- function(df){
  v <- vennCounts(df)
  fm <- matrix(nrow=2,ncol=2)
  for(i in 1:nrow(v)){
    fm[v[i,1]+1, v[i,2]+1] <- v[i,3]
  }
  fisher.test(fm)
}

testTFinDEoverrep <- function(tfc){
  pm <- getmarasitecountmatrix(tfc)
  outratio <- rep(NA,ncol(pm))
  outp <- outratio
  outnum <- outp
  
  for(i in 1:ncol(pm)){
    if(i%%50==0)
      print(i)
    v <- qtestfish(data.frame(
      rownames(mtpm) %in% rownames(pm)[pm[,i]>0],
      rownames(mtpm) %in% allde$ens_mouse[allde$anytime_mouse]))
    outnum[i] <- sum(rownames(mtpm) %in% rownames(pm)[pm[,i]>0])
    outratio[i] <- v$estimate
    outp[i] <- v$p.value
  }
  v<-data.frame(p.value=outp, ods=outratio, num=outnum, jasparname=colnames(pm))
  v <- v[order(v$ods, decreasing = TRUE),]
  v$odsrank <- (1:nrow(v))/nrow(v)
  v
}
#Problem: Ridiculous p-values because of the large number of genes.
v<-testTFinDEoverrep(atac.mouse$cons_tfc)
#Yy1 0.11. chip_batf 0.098. sp*. fli1 0.064. zscan22 0.05. etv6 0.12   ctcf, sp*, zscan22 < 0.05
#Gata3: 0.0075 Pou6f1: 0.009 Runx3 0.01  Nfatc1/3: 0.025  Runx1: 0.05   Irf4: 0.06  Batf::jun 0.09
#Runx3 5x more sites than gata3, still same ratio
v
v[grep("Stat",v$jasparname),]

v2<-testTFinDEoverrep(atac.mouse$noncons_tfc)
v2
#Sp1,2,4: <0.01 ctcf: 0.036 fli1: 0.05  irf4: 0.08  chip_batf: 0.1  
v2[grep("Stat",v2$jasparname),]   #Stat6 same here, 15%. conserved

#Most conserved regulators
v3 <- smerge(sqldf("select jasparname, odsrank as cons from v"),sqldf("select jasparname, odsrank as noncons from v2"))
plotdottext(v3$cons, v3$noncons, labels=v3$jasparname, cex = 1, donew=TRUE)
v4 <- intersect(
  v$jasparname[v$odsrank<0.20],
  v2$jasparname[v2$odsrank<0.20])
v4
sgenescorer2_matrix[c("Runx1","Fosl1","Irf1","Sp2","Sp1","Klf15","Klf14","Sp4","Patz1","Irf4","Batf","Sp2","Stat1","Stat2","Stat6"),]
# [1] "Runx1"        "chip_Irf4"    "Stat1::stat2" "Jun(var.2)"   "Junb"         "Fosl1"        "chip_Batf"    "Irf1"         "Sp2"         
# [10] "Sp1"          "Fos::jun"     "Klf15"        "Stat6"        "Klf14"        "Sp4"          "Patz1"        "Klf1"        
#Sp1,Klf15 -> Il4   klf15 not really expressed
#Irf4,Stat2 -> il13   weakly


