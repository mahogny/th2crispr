
######################################################################
############# Compare with ChIPseq data 3#############################
######################################################################


##########################################
## Load all chipseq annotation in a folder
allchipanno <- list()
for(f in list.files("../ext_chip2/", pattern = "*.anno")){
  print(f)
  anno_x <- read.csv(sprintf("../ext_chip2/%s",f),sep="\t", stringsAsFactors = FALSE)
  # fname  <- str_split_fixed(f,".bed",2)[1]
  # g <- regexpr("_[^_]*$", fname)
  # fname <- str_sub(fname,1,g[1]-1)
  allchipanno[[fname]] <- anno_x
}
allchipanno_fullname <- names(allchipanno)

shorten_anno_name <- function(fname){
  fname  <- str_split_fixed(f,".bed",2)[1]
  g <- regexpr("_[^_]*$", fname)
  fname <- str_sub(fname,1,g[1]-1)
  fname
}

sapply(allchipanno_fullname, shorten_anno_name)


sort(names(allchipanno))


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

######### Binarize
binarizematrix <- function(dat){
  dat[dat>0] <- 1
  dat
}


######### tSNE on ChIPseq  -- continuous
set.seed(0)
d <- stats::dist(t((dat)))
rtsne_out <- Rtsne(d,is_distance = TRUE, perplexity=20, verbose = TRUE, max_iter=15000,dims = 2)
pdf("out_oe/chipseq_tsne_cont.pdf")
plot(rtsne_out$Y[,], pch=16, main='',xlab="",ylab="",cex=0)  
text(rtsne_out$Y[,1], rtsne_out$Y[,2], labels=str_split_fixed(colnames(dat),"_",2)[,1],cex=1)
dev.off()
pdf("out_oe/chipseq_tsne_binary_long.pdf")
plot(rtsne_out$Y[,], pch=16, main='',xlab="",ylab="",cex=0)  
text(rtsne_out$Y[,1], rtsne_out$Y[,2], labels=colnames(dat),cex=1)
dev.off()


######### tSNE on ChIPseq  -- binary
set.seed(0)
d <- stats::dist(t((binarizematrix(dat))))
#d <- stats::dist(t((dat)))
rtsne_out <- Rtsne(d,is_distance = TRUE, perplexity=20, verbose = TRUE, max_iter=15000,dims = 2)
pdf("out_oe/chipseq_tsne_binary.pdf")
plot(rtsne_out$Y[,], pch=16, main='',xlab="",ylab="",cex=0)  
text(rtsne_out$Y[,1], rtsne_out$Y[,2], labels=str_split_fixed(colnames(dat),"_",2)[,1],cex=1)
dev.off()
pdf("out_oe/chipseq_tsne_binary_long.pdf")
plot(rtsne_out$Y[,], pch=16, main='',xlab="",ylab="",cex=0)  
text(rtsne_out$Y[,1], rtsne_out$Y[,2], labels=colnames(dat),cex=1)
dev.off()




######### Correlation based comparison, #peaks
pdf("out_oe/cor_global_cont.pdf")
thecor <- cor(dat)
heatmap(thecor,scale = "none") # Foxp3, Pparg, Gata3 closest, then Stat6.    Batf, Irf4, Xbp1, Bhlhe40
dev.off()

heatmap(thecor^10,scale = "none") 
heatmap(cor(thecor),scale = "none") 
hist(as.double(thecor))


######### Correlation based comparison, peak-or-not
pdf("out_oe/cor_global_binary.pdf")
thecor <- cor(binarizematrix(dat)) ####### stat6 + xp1.   Irf4, Batf, Pparg, Gata3
heatmap(thecor,scale = "none")
dev.off()

### how similar are Irf4 and Pparg overexpression?


######### "Local" correlation, continuous
pdf("out_oe/cor_local_cont.pdf")
thecor <- matrix(nrow=ncol(dat), ncol=ncol(dat))
for(i in 1:ncol(dat)){
  for(j in 1:ncol(dat)){
    keep <- dat[,i]>0 | dat[,j]>0
    thecor[i,j] <- cor(dat[keep,i], dat[keep,j])
  }
}
colnames(thecor) <- colnames(dat)
rownames(thecor) <- colnames(dat)
heatmap(thecor,scale = "none")    #stat6 and foxp3.  nearby xbp1.   pparg linked to all. batf/irf4/pparg own cluster
#heatmap(zerodiag(thecor),scale = "none")
dev.off()

######### "Local" correlation, binary
pdf("out_oe/cor_local_binary.pdf")
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


######### Jaccard index, digital version
pdf("out_oe/jaccard_binary.pdf")
thecor <- matrix(nrow=ncol(dat), ncol=ncol(dat))
for(i in 1:ncol(dat)){
  for(j in 1:ncol(dat)){
    thecor[i,j] <- sum(dat[,i] & dat[,j]) / sum(dat[,i] | dat[,j])
  }
}
colnames(thecor) <- colnames(dat)
rownames(thecor) <- colnames(dat)
#heatmap(thecor,scale = "none")
heatmap(zerodiag(thecor),scale = "none")
dev.off()


######### Jaccard index, continuous version
pdf("out_oe/jaccard_con.pdf")
thecor <- matrix(nrow=ncol(dat), ncol=ncol(dat))
for(i in 1:ncol(dat)){
  for(j in 1:ncol(dat)){
    #    thecor[i,j] <- sum(pmin(dat[,i], dat[,j])) / nrow(dat)
    thecor[i,j] <- sum(pmin(dat[,i], dat[,j])) / sum(pmax(dat[,i], dat[,j]))
  }
}
colnames(thecor) <- colnames(dat)
rownames(thecor) <- colnames(dat)
#heatmap(thecor,scale = "none")
heatmap(zerodiag(thecor),scale = "none")
dev.off()


###### alternative normalization. general mess, don't use
thecor <- matrix(nrow=ncol(dat), ncol=ncol(dat))
for(i in 1:ncol(dat)){
  for(j in 1:ncol(dat)){
    thecor[i,j] <- mean(dat[,i]*dat[,j]) / (sum(dat[,i]*dat[,i])*sum(dat[,j]*dat[,j]))
  }
}
colnames(thecor) <- colnames(dat)
rownames(thecor) <- colnames(dat)
heatmap(thecor,scale = "none")
heatmap(log10(thecor),scale = "none")
heatmap(zerodiag(thecor),scale = "none")
