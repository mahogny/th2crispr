###################################################################################################
################################################################################################### 
########## Functions for use on the cluster only  #################################################
###################################################################################################
################################################################################################### 
###################################################################################################


#### for one gene vs control
resFull <- results(ddsFull, contrast=c("ko","Il2","Thy1control"), parallel = TRUE, BPPARAM = MulticoreParam(workers=20))
resFull <- resFull[order(resFull$padj),]
x <- data.frame(fc=resFull$log2FoldChange, padj=resFull$padj, mgi_symbol=togenesym2(rownames(resFull)), stringsAsFactors = FALSE)
head( x , n=100)



##### for one gene vs mid-position
testvsmid <- function(var){
  of <- rep(0, length(resultsNames(ddsFull)))
  of[resultsNames(ddsFull)==var] <- 1
  resFull <- results(ddsFull, contrast=of, parallel = TRUE, BPPARAM = MulticoreParam(workers=20))
  resFull <- resFull[order(resFull$log2FoldChange),]
  resFull <- resFull[order(resFull$padj),]
  data.frame(fc=resFull$log2FoldChange, padj=resFull$padj, mgi_symbol=togenesym2(rownames(resFull)), stringsAsFactors = FALSE)
}
x <- testvsmid("koCcdc134_2_3")
x <- testvsmid("koGata3_5_7")
#x <- testvsmid("koEtv2")
#x <- testvsmid("koIl4")
#x <- testvsmid("koGata3")
x[1:100,]

#head( x , n=100)


##### for several probes vs mid-position
testprobesvsmid <- function(var){
  pid <- grep(var, resultsNames(ddsFull))
  of <- rep(0, length(resultsNames(ddsFull)))
  of[pid] <- 1/length(pid)
  resFull <- results(ddsFull, contrast=of, parallel = TRUE, BPPARAM = MulticoreParam(workers=20))
  resFull <- resFull[order(resFull$log2FoldChange),]
  resFull <- resFull[order(resFull$padj),]
  data.frame(fc=resFull$log2FoldChange, padj=resFull$padj, mgi_symbol=togenesym2(rownames(resFull)), stringsAsFactors = FALSE)
}
x <- testprobesvsmid("koCcdc134")
#x <- testprobesvsmid("koEtv2")
#x <- testprobesvsmid("koIl4")
#x <- testprobesvsmid("koGata3")
#x <- testprobesvsmid("koGata3")
x[1:100,]




##### for one gene vs average
testvsav <- function(var){
  of <- rep(0, length(resultsNames(ddsFull)))
  of[resultsNames(ddsFull)==var] <- 1
  tk <- grep("ko",resultsNames(ddsFull))
  of[tk] <- of[tk] - 1/length(tk)    
  resFull <- results(ddsFull, contrast=of, parallel = TRUE, BPPARAM = MulticoreParam(workers=20))
  resFull <- resFull[order(resFull$padj),]
  data.frame(fc=resFull$log2FoldChange, padj=resFull$padj, mgi_symbol=togenesym2(rownames(resFull)), stringsAsFactors = FALSE)
}
#x <- testvsav("koIl4")
x <- testvsav("koGata3_5_7")
x <- testvsav("koGata3_5_8")
x[1:100,]


testprobesvsav <- function(var){
  of <- rep(0, length(resultsNames(ddsFull)))
  pid <- grep(var, resultsNames(ddsFull))
  of[pid] <- 1/length(pid)
  tk <- grep("ko",resultsNames(ddsFull))
  of[tk] <- of[tk] - 1/length(tk)    
  resFull <- results(ddsFull, contrast=of, parallel = TRUE, BPPARAM = MulticoreParam(workers=20))
  resFull <- resFull[order(resFull$log2FoldChange),]
  resFull <- resFull[order(resFull$padj),]
  data.frame(fc=resFull$log2FoldChange, padj=resFull$padj, mgi_symbol=togenesym2(rownames(resFull)), stringsAsFactors = FALSE)
}
x <- testprobesvsav("koCcdc134")
x <- testprobesvsav("koGata3")
head( x , n=100)


##### for one bfp
of <- rep(0, length(resultsNames(ddsFull)))
of[resultsNames(ddsFull)=="bfp"] <- 1
resFull <- results(ddsFull, contrast=of, parallel = TRUE, BPPARAM = MulticoreParam(workers=20))
resFull <- resFull[order(resFull$padj),]
x <- data.frame(fc=resFull$log2FoldChange, padj=resFull$padj, mgi_symbol=togenesym2(rownames(resFull)), stringsAsFactors = FALSE)
head( x , n=100)


##### for several bfp
of <- rep(0, length(resultsNames(ddsFull)))
of[resultsNames(ddsFull)=="bfp1"] <- 1
resFull <- results(ddsFull, contrast=of, parallel = TRUE, BPPARAM = MulticoreParam(workers=20))
resFull <- resFull[order(resFull$padj),]
x <- data.frame(fc=resFull$log2FoldChange, padj=resFull$padj, mgi_symbol=togenesym2(rownames(resFull)), stringsAsFactors = FALSE)
head( x , n=100)



#x$padj<0.99 & (1:nrow(x))<200
stopgosym(x$mgi_symbol[1:400], x$mgi_symbol)[,c(2,6)]

#Adam8 is an artifact in the WT control?
which(x[,3]=="Adam8")
which(x[,3]=="Cd3d")




