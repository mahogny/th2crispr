

##########################################
## Load all chipseq annotation in a folder
allchipanno <- list()
for(f in list.files("ext_chip2/", pattern = "*.anno")){
  print(f)
  anno_x <- read.csv(sprintf("ext_chip2/%s",f),sep="\t", stringsAsFactors = FALSE)
  fname  <- str_split_fixed(f,".bed",2)[1]
  g <- regexpr("_[^_]*$", fname)
  fname <- str_sub(fname,1,g[1]-1)
  allchipanno[[fname]] <- anno_x
}
sort(names(allchipanno))

## Extract the gene name
allchipanno_gene <- str_split_fixed(str_split_fixed(names(allchipanno),"-",2)[,1],"_",2)[,1]
allchipanno_gene <- sapply(allchipanno_gene, normalizesym)
names(allchipanno_gene)<-NULL


##########################################
## Make on large peak count table of the chipseq data
overlapAnnot <- function(annotlist){
  allg <- c()
  for(v in annotlist){
    allg <- unique(c(allg,v$Nearest.Ensembl))
  }
  out <- data.frame(row.names=allg, stringsAsFactors = FALSE)
  for(n in names(annotlist)){
    x <- data.frame(gname=annotlist[[n]]$Nearest.Ensembl, stringsAsFactors = FALSE)
    y <- sqldf("select gname, count(gname) as c from x group by gname")
    out[,n] <- 0
    out[y$gname,n] <- y$c
  }
  out
}



##########################################
## Write 0s on the matrix diagonal
zerodiag <- function(thecor){
  for(i in 1:nrow(thecor)){
    thecor[i,i] <- 0
  }
  thecor  
}


dat <- overlapAnnot(allchipanno)
dat_npeak <- apply(dat,2,sum)
dat <- dat[,dat_npeak>1000]

dat_npeak
# dat <- overlapAnnot(list(
#   Bhlhe40=anno_bhlhe40, 
#   Pparg  =anno_pparg, 
#   
#   Gata3  =anno_gata3, 
#   Batf   =anno_batf, 
#   
#   #anno_rora, 
#   #anno_zfp207,
#   Yy1=anno_yy1,
#   Fli1_1=anno_fli1_1,
#   #  Fli1_2=anno_fli1_2,   #could merge the two. here hiding one
#   Irf4=anno_irf4,
#   Rorc=anno_rorc_th0,
#   Xbp1=anno_xbp1,
#   Foxp3=anno_foxp3,
#   #Bmi1_1=anno_bmi1_1,   #misbehaves
#   Bmi1_2=anno_bmi1_2,
#   Cdk9=anno_cdk9,
#   #Ets1=anno_ets1,   #very different
#   
#   #weird polycomb stuff
#   #  Ezh2_1=anno_ezh2_1,
#   #  Ezh2_2=anno_ezh2_2,
#   #  Men1_1=anno_men1_1,
#   #  Men1_2=anno_men1_2,
#   
#   #Smarca4_1=anno_smarca4_1,   #weird thing, and inconsistent, rather not include
#   #Smarca4_2=anno_smarca4_2,
#   Spns1=anno_spns1,
#   Stat6=anno_stat6
# ))



#Add batf and gata3 here? 

#vennDiagram(vennCounts(dat))



######### Binarize
binarizematrix <- function(dat){
  dat[dat>0] <- 1
  dat
}


# ######### tSNE on ChIPseq  -- continuous
# set.seed(0)
# d <- stats::dist(t((dat)))
# rtsne_out <- Rtsne(d,is_distance = TRUE, perplexity=20, verbose = TRUE, max_iter=15000,dims = 2)
# pdf("out/chipseq_tsne_cont.pdf")
# plot(rtsne_out$Y[,], pch=16, main='',xlab="",ylab="",cex=0)  
# text(rtsne_out$Y[,1], rtsne_out$Y[,2], labels=str_split_fixed(colnames(dat),"_",2)[,1],cex=1)
# dev.off()
# #pdf("out/chipseq_tsne_binary_long.pdf")
# plot(rtsne_out$Y[,], pch=16, main='',xlab="",ylab="",cex=0)  
# text(rtsne_out$Y[,1], rtsne_out$Y[,2], labels=colnames(dat),cex=1)
# #dev.off()



######### tSNE on ChIPseq  -- binary
datb <- binarizematrix(dat)
dat_npeak <- apply(dat,2,sum)
datb <- removeBatchEffect(datb, covariates = dat_npeak)
d <- stats::dist(t((datb)))
set.seed(111)
rtsne_out <- Rtsne(d,is_distance = TRUE, perplexity=10, verbose = TRUE, max_iter=30000,dims = 2)
pdf("out_chip/chipseq_tsne_binary.pdf")
simpname <- str_split_fixed(colnames(dat),"_",2)[,1]
simpname <- sapply(simpname, str_to_upper)
plot(rtsne_out$Y[,], pch=16, main='',xlab="",ylab="",cex=0)  
text(rtsne_out$Y[,1], rtsne_out$Y[,2], labels=simpname,cex=1, 
     col = rgb(dat_npeak/max(dat_npeak),0,0))
text(rtsne_out$Y[,1], rtsne_out$Y[,2], labels=simpname,cex=1, 
     col = rgb(dat_npeak/max(dat_npeak),0,0))

dev.off()
pdf("out_chip/chipseq_tsne_binary_long.pdf")
plot(rtsne_out$Y[,], pch=16, main='',xlab="",ylab="",cex=0)  
text(rtsne_out$Y[,1], rtsne_out$Y[,2], labels=colnames(dat),cex=1)
dev.off()





######### tSNE on ChIPseq  -- binary
datb <- binarizematrix(dat)
dat_npeak <- apply(dat,2,sum)
datb <- removeBatchEffect(datb, covariates = dat_npeak)
d <- stats::dist(t((datb)))
set.seed(222)
rtsne_out <- Rtsne(d,is_distance = TRUE, perplexity=10, verbose = TRUE, max_iter=30000,dims = 2)
pdf("out_chip/chipseq_tsne_binary.pdf")
simpname <- str_split_fixed(colnames(dat),"_",2)[,1]
simpname <- sapply(simpname, str_to_upper)
plot(rtsne_out$Y[,], pch=16, main='',xlab="",ylab="",cex=0)  
text(rtsne_out$Y[,1], rtsne_out$Y[,2], labels=simpname,cex=1, 
     col = rgb(dat_npeak/max(dat_npeak),0,0))
dev.off()
pdf("out_chip/chipseq_tsne_binary_long.pdf")
plot(rtsne_out$Y[,], pch=16, main='',xlab="",ylab="",cex=0)  
text(rtsne_out$Y[,1], rtsne_out$Y[,2], labels=colnames(dat),cex=1)
dev.off()



# 
# ######### Correlation based comparison, #peaks
# pdf("out/cor_global_cont.pdf")
# thecor <- cor(dat)
# heatmap(thecor,scale = "none") # Foxp3, Pparg, Gata3 closest, then Stat6.    Batf, Irf4, Xbp1, Bhlhe40
# dev.off()
# 
# heatmap(thecor^10,scale = "none") 
# heatmap(cor(thecor),scale = "none") 
# hist(as.double(thecor))
# 
# 
######### Correlation based comparison, peak-or-not
pdf("out/cor_global_binary.pdf")
thecor <- cor(binarizematrix(dat)) ####### stat6 + xp1.   Irf4, Batf, Pparg, Gata3
heatmap(thecor,scale = "none")
dev.off()

# ### how similar are Irf4 and Pparg overexpression?
# 
# 
# ######### "Local" correlation, continuous
# pdf("out/cor_local_cont.pdf")
# thecor <- matrix(nrow=ncol(dat), ncol=ncol(dat))
# for(i in 1:ncol(dat)){
#   for(j in 1:ncol(dat)){
#     keep <- dat[,i]>0 | dat[,j]>0
#     thecor[i,j] <- cor(dat[keep,i], dat[keep,j])
#   }
# }
# colnames(thecor) <- colnames(dat)
# rownames(thecor) <- colnames(dat)
# heatmap(thecor,scale = "none")    #stat6 and foxp3.  nearby xbp1.   pparg linked to all. batf/irf4/pparg own cluster
# #heatmap(zerodiag(thecor),scale = "none")
# dev.off()
# 
######### "Local" correlation, binary
pdf("out/cor_local_binary.pdf")
thecor <- matrix(nrow=ncol(dat), ncol=ncol(dat))
bdat <- binarizematrix(dat)
for(i in 1:ncol(dat)){
  for(j in 1:ncol(dat)){
    keep <- bdat[,i]>0 | bdat[,j]>0
    thecor[i,j] <- cor(bdat[keep,i], bdat[keep,j])
  }
}
colnames(thecor) <- colnames(dat)
rownames(thecor) <- colnames(dat)
heatmap(thecor,scale = "none")    #stat6 and foxp3.  nearby xbp1.   pparg linked to all. batf/irf4/pparg own cluster
#heatmap(zerodiag(thecor),scale = "none")
dev.off()
# 
# 
# ######### Jaccard index, digital version
# pdf("out/jaccard_binary.pdf")
# thecor <- matrix(nrow=ncol(dat), ncol=ncol(dat))
# for(i in 1:ncol(dat)){
#   for(j in 1:ncol(dat)){
#     thecor[i,j] <- sum(dat[,i] & dat[,j]) / sum(dat[,i] | dat[,j])
#   }
# }
# colnames(thecor) <- colnames(dat)
# rownames(thecor) <- colnames(dat)
# #heatmap(thecor,scale = "none")
# heatmap(zerodiag(thecor),scale = "none")
# dev.off()
# 
# 
# ######### Jaccard index, continuous version
# pdf("out/jaccard_con.pdf")
# thecor <- matrix(nrow=ncol(dat), ncol=ncol(dat))
# for(i in 1:ncol(dat)){
#   for(j in 1:ncol(dat)){
#     #    thecor[i,j] <- sum(pmin(dat[,i], dat[,j])) / nrow(dat)
#     thecor[i,j] <- sum(pmin(dat[,i], dat[,j])) / sum(pmax(dat[,i], dat[,j]))
#   }
# }
# colnames(thecor) <- colnames(dat)
# rownames(thecor) <- colnames(dat)
# #heatmap(thecor,scale = "none")
# heatmap(zerodiag(thecor),scale = "none")
# dev.off()
# 
# 
# #yy1 and bhlhe40 seem to belong together
# #rorc and spns1 definitely hangs together with bmi and men1
# #men1_1 seem to be an outlier
# #irf4 batf pparg belong together, with xbp1 being an extension
# #stat6 binds with everything but yy1?
# 
# 
# ###### alternative normalization. general mess, don't use
# thecor <- matrix(nrow=ncol(dat), ncol=ncol(dat))
# for(i in 1:ncol(dat)){
#   for(j in 1:ncol(dat)){
#     thecor[i,j] <- mean(dat[,i]*dat[,j]) / (sum(dat[,i]*dat[,i])*sum(dat[,j]*dat[,j]))
#   }
# }
# colnames(thecor) <- colnames(dat)
# rownames(thecor) <- colnames(dat)
# heatmap(thecor,scale = "none")
# heatmap(log10(thecor),scale = "none")
# heatmap(zerodiag(thecor),scale = "none")
# 
# 
# 
# sort(oe_de$padj[ensconvert$ensembl_gene_id[ensconvert$mgi_symbol=="Lrig1"],])  #Scara3 and Lrrc40 ... nope
# 
# 
# 
# 
# 
# 

##############################
############################## activation/differentiation
##############################



dat_act  <- dat[rownames(dat) %in% rownames(resTcActDiffP)[!is.na(resTcActDiffP$act) & resTcActDiffP$act<1e-20],]
dat_diff <- dat[rownames(dat) %in% rownames(resTcActDiffP)[!is.na(resTcActDiffP$diff) & resTcActDiffP$diff<1e-5],]
dim(dat_act)
dim(dat_diff)

dat_act  <- dat[rownames(dat) %in% rownames(resTcActDiffP)[!is.na(resTcActDiffP$act) & resTcActDiffP$act<1e-50],]
dat_diff <- dat[rownames(dat) %in% rownames(resTcActDiffP)[!is.na(resTcActDiffP$diff) & resTcActDiffP$diff<1e-7],]
dim(dat_act)
dim(dat_diff)

dat_act  <- dat[rownames(dat) %in% rownames(resTcActDiffP)[!is.na(resTcActDiffP$act) & resTcActDiffP$act<1e-80],]
dat_diff <- dat[rownames(dat) %in% rownames(resTcActDiffP)[!is.na(resTcActDiffP$diff) & resTcActDiffP$diff<1e-10],]
dim(dat_act)
dim(dat_diff)

dat_diffact  <- dat[rownames(dat) %in% c(
  rownames(resTcActDiffP)[!is.na(resTcActDiffP$act) & resTcActDiffP$act<1e-50],
  rownames(resTcActDiffP)[!is.na(resTcActDiffP$diff) & resTcActDiffP$diff<1e-7]
),]


sort(apply(dat_act,2,sum))
sort(apply(dat_diff,2,sum))

######### By overlap

diffact_chip_overlap <- data.frame(
  name=colnames(dat_act),
  act=apply(dat_act,2,sum),
  diff=apply(dat_diff,2,sum),
  stringsAsFactors = FALSE
)
# diffact_chip_overlap <- data.frame(
#   name=colnames(dat_act),
#   act=log10(1+apply(dat_act,2,sum)),
#   diff=log10(1+apply(dat_diff,2,sum)),
#   stringsAsFactors = FALSE
# )
plot(diffact_chip_overlap$diff, 
     diffact_chip_overlap$act, 
     pch=16, main='',xlab="",ylab="",cex=0)  
text(diffact_chip_overlap$diff, 
     diffact_chip_overlap$act, 
     labels=diffact_chip_overlap$name,
     cex=1)


plot(diffact_chip_overlap$diff/diffact_chip_overlap$act, 
     diffact_chip_overlap$act, 
     pch=16, main='',xlab="",ylab="",cex=0)  
text(diffact_chip_overlap$diff/diffact_chip_overlap$act, 
     diffact_chip_overlap$act, 
     labels=diffact_chip_overlap$name,
     cex=1)


sort(apply(dat_diffact,2,sum))

# for(i in 1:ncol(dat_act)){
#   dat_act[,i] <- dat_act[,i] / sum(dat_act[,i])
# }


######### Activation
dat_act <- dat_act[,apply(dat_act,2,sum)>10]
dim(dat_act)
set.seed(0)
d <- stats::dist(t((binarizematrix(dat_act))))
rtsne_out <- Rtsne(d,is_distance = TRUE, perplexity=12, verbose = TRUE, max_iter=20000,dims = 2)
#pdf("out/chipseq_tsne_binary.pdf")
plot(rtsne_out$Y[,], pch=16, main='',xlab="",ylab="",cex=0)  
text(rtsne_out$Y[,1], rtsne_out$Y[,2], labels=str_split_fixed(colnames(dat_act),"_",2)[,1],cex=1)
#dev.off()
#pdf("out/chipseq_tsne_binary_long.pdf")
plot(rtsne_out$Y[,], pch=16, main='',xlab="",ylab="",cex=0)  
text(rtsne_out$Y[,1], rtsne_out$Y[,2], labels=colnames(dat),cex=1)
#dev.off()



######### Differentiation
dat_diff <- dat_diff[,apply(dat_diff,2,sum)>10]
dim(dat_diff)
set.seed(0)
d <- stats::dist(t((binarizematrix(dat_diff))))
rtsne_out <- Rtsne(d,is_distance = TRUE, perplexity=10, verbose = TRUE, max_iter=20000,dims = 2)
#pdf("out/chipseq_tsne_binary.pdf")
plot(rtsne_out$Y[,], pch=16, main='',xlab="",ylab="",cex=0)  
text(rtsne_out$Y[,1], rtsne_out$Y[,2], labels=str_split_fixed(colnames(dat_diff),"_",2)[,1],cex=1)
#dev.off()
#pdf("out/chipseq_tsne_binary_long.pdf")
plot(rtsne_out$Y[,], pch=16, main='',xlab="",ylab="",cex=0)  
text(rtsne_out$Y[,1], rtsne_out$Y[,2], labels=colnames(dat),cex=1)
#dev.off()


######### Differentiation / activation 
set.seed(0)
d <- stats::dist(t((binarizematrix(dat_diffact))))
rtsne_out <- Rtsne(d,is_distance = TRUE, perplexity=10, verbose = TRUE, max_iter=20000,dims = 2)
#pdf("out/chipseq_tsne_binary.pdf")
plot(rtsne_out$Y[,], pch=16, main='',xlab="",ylab="",cex=0)  
text(rtsne_out$Y[,1], rtsne_out$Y[,2], labels=str_split_fixed(colnames(dat_diffact),"_",2)[,1],cex=1)
#dev.off()
#pdf("out/chipseq_tsne_binary_long.pdf")
plot(rtsne_out$Y[,], pch=16, main='',xlab="",ylab="",cex=0)  
text(rtsne_out$Y[,1], rtsne_out$Y[,2], labels=colnames(dat),cex=1)
#dev.off()




######### Differentiation / activation   -- reduced
dat_diffact_red <- dat_diffact[,apply(dat_diffact,2,sum)>100]
set.seed(0)
d <- stats::dist(t((binarizematrix(dat_diffact_red))))
rtsne_out <- Rtsne(d,is_distance = TRUE, perplexity=10, verbose = TRUE, max_iter=20000,dims = 2)
#pdf("out/chipseq_tsne_binary.pdf")
plot(rtsne_out$Y[,], pch=16, main='',xlab="",ylab="",cex=0)  
text(rtsne_out$Y[,1], rtsne_out$Y[,2], labels=str_split_fixed(colnames(dat_diffact_red),"_",2)[,1],cex=1)
#dev.off()
#pdf("out/chipseq_tsne_binary_long.pdf")
plot(rtsne_out$Y[,], pch=16, main='',xlab="",ylab="",cex=0)  
text(rtsne_out$Y[,1], rtsne_out$Y[,2], labels=colnames(dat),cex=1)
#dev.off()
