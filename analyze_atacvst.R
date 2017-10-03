#############################################################################################
###                                                                                       ###
###             Part of the paper ...                                                     ###
###             Author: Johan Henriksson (mahogny@areta.org)                              ###
###                                                                                       ###
###             This code plots ATAC and ChIPseq stats over time and space                ###
###                                                                                       ###
#############################################################################################





################################################################################
############ Plot ATAC over time, version used in the paper in the end #########
################################################################################


#### This function works better. palette() is a weird function
plotatacvst.2 <- function(thelev,pgenes,ymax=10000, reorder=TRUE){
  thelev["Ctcf",] <- thelev["Ctcf",]/4
  rownames(thelev)[rownames(thelev)=="Ctcf"] <- "Ctcf/4"
  thelev <- thelev[pgenes,]
  print(thelev)
  if(reorder){
    no <- order(thelev[,6],decreasing = TRUE)
    print(no)
    pgenes <- pgenes[no]
    thelev <- thelev[no,]
  }
  print(thelev)
  print(pgenes)
  pcol <- brewer.pal(length(pgenes),name = "Set3")
  #  pcol <- palette(brewer.pal(length(pgenes),name = "Set3"))
  print(pcol)
  tp <- c("0h","2h","4h","24h","48h","72h")
  lwd <- 4
  plot(c(0,0),c(0,0),cex=0,xaxt = "n",
       ylab="Peak count",xlab="",xlim=c(1,6),ylim=c(0,ymax),lwd=lwd) 
  for(i in 1:length(pgenes)){
    lines(as.matrix(thelev)[i,],type="l",col=pcol[i],lwd=lwd) 
  }
  legend(2,ymax,pgenes,cex=1.2,fill = pcol,y.intersp=0.65,box.lwd = 0)
  axis(1, at = 1:length(tp), labels=tp)
}

head(rownames(tfatall_normtime2),n=100) #negative
tail(rownames(tfatall_normtime2),n=10)

v <- levatac.mouse.norm
v$crank <- 1:nrow(v) #from low to high

v<-smerge(
  map_jaspar_namegenesym,
  data.frame(jasparname=rownames(levatac.mouse.norm), crank=1:nrow(levatac.mouse.norm)))
v <- smerge(v, ensconvert)
v <- v[order(v$crank),]
v <- sgenescorer2_matrix_ensmouse[v$ensembl_gene_id,]
head(v,n=50)
tail(v,n=50)


pdf("out_atacvst/new_exampleTF.pdf")
plotatacvst.2(levatac.mouse.norm, c("Foxo1","Gata3_1","Batf::jun","Yy1","Stat6","E2f7","Runx1"),ymax=3)
dev.off()


##################################################################
############ Plot ATAC over time, individual #####################
##################################################################

plotatacvst1 <- function(gene){
#  png(sprintf("atacvst/gene %s.png",gene),w=400)
  barplot(as.matrix(tfattall)[gene,],ylab=gene)
 # dev.off()
}

plotatacvst1("Gmeb2") #Gmeb2 is surprisingly high early on


plot(tfattall[,2],tfattall[,6],cex=0.2,col="white")
text(tfattall[,2],tfattall[,6],rownames(tfattall),cex=0.6,col="red",xlab="Early activity",ylab="Late activity")


rlog <- function(x) log(x+1)
actearly <- tfattall[,2]
actlate <- apply(tfattall[,4:6],1,mean)
explate <- mean(actlate/actearly)*actearly

png("atacvst/early vs late.png",w=1200,height = 800)
plot(rlog(actearly),rlog(actlate),col="white")
text(rlog(actearly),rlog(actlate),rownames(tfattall),cex=0.8,col="red",xlab="Early activity",ylab="Late activity")
dev.off()

png("atacvst/early vs late regressed.png",w=1200,height = 800)
plot(rlog(actearly),rlog(actlate-explate),col="white")
text(rlog(actearly),rlog(actlate-explate),rownames(tfattall),cex=0.8,col="red",xlab="Early activity",ylab="Late activity")
dev.off()
#Top side: Fosl2, Batf, Jun, Nfe2
#Bottom side: 




##################################################################
#### Perform tSNE on the ATAC motif peaks (using time only) ######
##################################################################

set.seed(0)
tfattall_red <- tfattall[rownames(tfattall) %in% expressed_atacTF_80,]
rtsne_out <- Rtsne(tfattall_red, pca = FALSE, verbose = TRUE,perplexity = 4, max_iter = 5000)

pdf("out_atacvst/atac_tsne_time.pdf")
plot(rtsne_out$Y, pch = 20,
     cex = 0, cex.axis = 1.25, cex.lab = 1.25, cex.main = 1.5,
     xlab = "t-SNE dimension 1", ylab = "t-SNE dimension 2",
     main = "TF binding similarity over time and space")

thecol <- col_from_ael(score_ael[rownames(tfattall_red)])

rtsne_jitter <- 2
text(rtsne_out$Y[,1]+runif(nrow(rtsne_out$Y))*rtsne_jitter,
     rtsne_out$Y[,2]+runif(nrow(rtsne_out$Y))*rtsne_jitter,rownames(tfattall_red),cex=0.6,col=thecol)
dev.off()

csvtsne <- data.frame(x=rtsne_out$Y[,1],y=rtsne_out$Y[,2],tf=rownames(tfattall_red),time=wt, color=thecol)
csvtsne <- merge(csvtsne,genemotif)
write.csv(csvtsne,"out_atacvst/atac_tsne_time.csv",row.names = FALSE,quote = FALSE)




#######################################################################
#### Perform tSNE on the ATAC motif peaks (on time & gene level) ######
#######################################################################

## Merge time and space
d1 <- getmarasitecountmatrix(genetfcount_early)
d2 <- getmarasitecountmatrix(genetfcount_late)
d12 <- rbind(d1[-(1:2),],d2[-(1:2),])
d12 <- d12[,colnames(d12) %in% c(expressed_atacTF_50, "Etv2")]
fortsne <- t(d12)
class(fortsne)<-"double"

## Run tSNE
set.seed(0)
rtsne_out <- Rtsne(fortsne, pca = FALSE, verbose = TRUE,perplexity = 10, max_iter = 5000)

pdf("out_atacvst/atac_tsne_timespace.pdf")
set.seed(0)
rtsne_jitter <- 1
rtsne_outY <- cbind(
  rtsne_out$Y[,1]+runif(nrow(rtsne_out$Y))*rtsne_jitter,
  rtsne_out$Y[,2]+runif(nrow(rtsne_out$Y))*rtsne_jitter)
plot(rtsne_outY, pch = 20,
     cex = 0, cex.axis = 1.25, cex.lab = 1.25, cex.main = 1.5,
     xlab = "t-SNE dimension 1", ylab = "t-SNE dimension 2",
     main = "TF binding similarity over time and space")
thecol <- col_from_ael(score_ael[rownames(fortsne)])
#wt <- calc_score_ael(tfattall[rownames(fortsne),])
#thecol <- col_from_wt(wt)
#thecol <- kmcol[atackm$cluster[rownames(fortsne)]]  #For k-means color
text(rtsne_outY[,1],
     rtsne_outY[,2],rownames(fortsne),cex=0.6,col=thecol)
dev.off()

csvtsne <- data.frame(x=rtsne_out$Y[,1],y=rtsne_out$Y[,2],tf=rownames(fortsne),time=wt, color=thecol)
csvtsne <- merge(csvtsne,genemotif)
write.csv(csvtsne,"out_atacvst/atac_tsne_timespace.csv",row.names = FALSE,quote = FALSE)




#################################################################
############### Plot chip peak count  ###########################
#################################################################
plotChipPeakCount <- function(){
  pcol <- brewer.pal(4,name = "Set2")
  pname <- c("Batf","Irf4","Gata3","Stat6")
  plot(x=c(5,8,9,10),apply(dstat6[,4:7],2,sum),xaxt="n",xlab="",type="l",lwd=4,ylab="",
       xlim=c(1,10),ylim=c(0,9000),col=pcol[4],main="ChIP merged peak count")
  lines(x=c(9.5,10),apply(dbatf[,rep(7,2)],2,sum),lwd=4,lty=1,col=pcol[1])
  lines(x=c(9.5,10),apply(dirf4[,rep(7,2)],2,sum),lwd=4,lty=1,col=pcol[2])
  lines(x=c(8,9,10),apply(dgata3[,5:7],2,sum),lwd=4,col=pcol[3])
  axis(1, at=c(1,4,5,8,9,10), labels=c("Naive","2h","4h","24h","48h","72h"))
  legend(1,7500,pname,cex=1.2,fill = pcol,y.intersp=0.65,box.lwd = 0)
}
pdf("mara_summary/chip merged peak count.pdf",width = 6,height=5)
plotChipPeakCount()
dev.off()








#################################################################
###### ATAC peak categories over time ###########################
#################################################################

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

#plot(apply(hstat,2,sum),type="l",ylim=c(0,40e3))

mulfac <- c(100,5000,100,50,1000,
            10,1,1,1,10,
            1,3000)
for(i in 1:nrow(hstat))
  hstat[i,] <- mulfac[i]*hstat[i,]
plot(hstat[1,],type="l",xaxt = 'n',ylim=c(0,0.6))
for(i in 2:nrow(hstat))
  lines(1:ncol(hstat),hstat[i,],type="l",xaxt = 'n')
colnames(hstat) <- sprintf("%sh",htimes)
axis(side=1,at=1:ncol(hstat),labels=sprintf("%sh",htimes))
pdf("out_atacvst/atac_cat_rescale.pdf",width = 3)
barplot(hstat[-11,], col = brewer.pal(11,"Set3"), legend.text = TRUE)
dev.off()






##############################################################################################
######## Compare the ATAC predicted binding sites with the KOs. New version ##################
##############################################################################################


#### Perform one comparison KO/ATAC
compkoatac <- function(gene,jaspar,tss=10e3,peakcutoff=1,tpmcut=10){
  #Calculate num timepoints
  mapPeakInfo_red <- mapPeakInfo[mapPeakInfo$tf==jaspar,]
  mapMotifGene_red <- merge(mapPeakInfo_red,mapPeakGene_unfiltered)
  colnames(mapMotifGene_red)[6] <- "ensembl_gene_id"
  colnames(mapMotifGene_red)[7:12] <- sprintf("tp%s",1:6)
  #Apply distance cutoff
  mapMotifGene_red <- mapMotifGene_red[abs(mapMotifGene_red$TSS_distance)<tss,]
  #Count peaks  
  numtp <- apply(mapMotifGene_red[,sprintf("tp%s",1:6)],1,sum)
  
  #Merge: expression level. differential expression. num timepoints
  jsym <- map_jaspar_namegenesym[map_jaspar_namegenesym$mgi_symbol==gene,]$jasparname
  compkoatac<-data.frame(ensembl_gene_id=rownames(respval), de.p=respval[,gene], de.fc=resfc[,gene])
  compkoatac<-smerge(compkoatac,ensconvert)
  compkoatac <- smerge(compkoatac,
                       data.frame(ensembl_gene_id=mapMotifGene_red$ensembl_gene_id, tfcount=numtp),all.x=TRUE)
  compkoatac$tfcount[is.na(compkoatac$tfcount)]<-0

  #Only consider expressed genes
#  tpmcut<-10
  eId <- names(av_mtpm[apply(cbind(av_mtpm,av_mtpm0),1,max)>tpmcut,1])
  compkoatac <- compkoatac[compkoatac$ensembl_gene_id %in% eId,]
  
  t1<-t.test(
    log10(compkoatac$de.p[compkoatac$tfcount==0]),
    log10(compkoatac$de.p[compkoatac$tfcount>=peakcutoff])) #should be more negative

  list(
    fits=c(t1$estimate[2]<t1$estimate[1]),
    pval=c(t1$p.value,t2$p.value),
    diff=c(t1$estimate[2]-t1$estimate[1])  #should be positive
  )
}
#compkoatac("Tbx21","Tbx21",10e3) 
# compkoatac("Tbx21","Tbx21",50e3) 
# compkoatac("Tbx21","Tbx21",5000e3) 

#### Calculate a series of KO/ATAC comparisons as a function of TPM cutoff
calc_compkoatacvstpm <- function(gene,jaspar=NULL,peakcutoff,tpmcut=10,alltss=c(1000,2000,3000,4000, 5000*(1:20))){
  if(is.null(jaspar))
    jaspar=map_jaspar_namegenesym$jasparname[map_jaspar_namegenesym$mgi_symbol=="Irf4"]
  st.diff<-matrix(nrow=0,ncol=2)
  st.p<-matrix(nrow=0,ncol=2)
  for(tss in alltss){
    print(tss)
    cko <- compkoatac(gene, jaspar, tss=tss, peakcutoff = peakcutoff, tpmcut = tpmcut)
    st.diff <- c(st.diff,cko$diff)
    st.p <- rbind(st.p,cko$pval)
  }
  list(tss=alltss,diff=st.diff,p=st.p)  
}


#### Plot a number of KO/ATAC comparisons
plot_compkoatacvstpm <- function(){
  test_koatacvstpm_genes <- c("Tbx21","Stat6")
  test_koatacvstpm_col <- c("red","blue")
  alltss <- c(1000,2000,3000,4000, 5000*(1:40))
  for(i in 1:length(test_koatacvstpm_genes)){
    koatacvstpm <- calc_compkoatacvstpm(test_koatacvstpm_genes[i],peakcutoff = 1,tpmcut=1)
    koatacvstpm$diff <- koatacvstpm$diff-tail(koatacvstpm$diff,n=1) #Asymptote as baseline
    if(i==1)
      plot(koatacvstpm$tss,koatacvstpm$diff,type="l",col=test_koatacvstpm_col[i],
           ylim=c(0,0.05),xlab="< TSS distance",ylab="Gain P-value")
    else
      lines(koatacvstpm$tss,koatacvstpm$diff,col=test_koatacvstpm_col[i])
  }
  legend("topright",legend = test_koatacvstpm_genes, fill=test_koatacvstpm_col)
  lines(c(0,max(alltss)),c(0,0),col="gray")
}
pdf("out_atacvst/comp_ko_atac_vs_tpm.pdf")
plot_compkoatacvstpm()
dev.off()

