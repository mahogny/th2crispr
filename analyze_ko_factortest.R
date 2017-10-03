# These functions were used to analyze different linear models. They are not needed to reproduce any data in the paper


###################################################################################################
########### Solve linear model. Quick testing only - used DESeq2 separately for final data ########
###################################################################################################

#Model: Y = X B
y.train <- t(log10(1+as.matrix(dat_red)))
x.train <- as.matrix(model.matrix(f, ddsCol ))

if(FALSE){
  ## For use with python solvers
  write.table(x.train, "../python/x.csv", row.names = FALSE, col.names = FALSE, sep=",")
  write.table(y.train, "../python/y.csv", row.names = FALSE, col.names = FALSE, sep=",")
  b.train <- read.table("../python/fit.csv",sep=",")
  colnames(b.train) <- colnames(x.train)
  rownames(b.train) <- rownames(dat_red)
}

if(FALSE){
  ## Compute with regular lm
  b.train <- t(lm.fit(x.train, y.train)$coefficients)  #to be compatible with python, which is transposed
  dim(b.train)
}

###############################################
########### QC: check sources of variation ##############
###############################################
#Y=XB.
calcpartvar <- function(){
  av <- c()
  for(i in 1:ncol(b.train)){
    av <- c(av,
            sd(as.double(x.train[,i,drop=FALSE] %*% t(b.train[,i,drop=FALSE]))))
  }
  names(av)<-colnames(b.train)
  round(sort(av)/sum(av), digits = 5)
}
v<-calcpartvar()
v
mean(v[grep("ko",names(v))])

mean(levdepth[cellcondition_red$mouse==1])
mean(levdepth[cellcondition_red$mouse==2])
mean(levdepth[cellcondition_red$mouse==3])


#colnames(b.train)
#t(b.train[1,,drop=FALSE])

head(b.train)
#head(x.train)
#sum(abs(b.train$bfp-median(b.train$bfp))>0.1)



###############################################
########### check top DE genes ################
###############################################
showtopfc <- function(x,col,useall=TRUE, doplot=TRUE){
  v <- x[,col]
  if(useall){
    v <- v - apply(b.train[,grep("ko",colnames(b.train))],1,mean)
  }
  v <- v-median(v)
  if(doplot)
    hist(v,breaks=100)
  d <- data.frame(mgi_symbol=togenesym2(rownames(x)), fc=v, stringsAsFactors = FALSE)
  d <- d[order(abs(d$fc),decreasing = TRUE),]
  #  togenesym2( rownames(x)[order(v,decreasing = TRUE)] )
  d
}
showtopfc(b.train,"bfp")[1:40,]
showtopfc(b.train,"koIfngr1")[1:40,]
showtopfc(b.train,"koIl4")[1:40,]
showtopfc(b.train,"koIl13")[1:40,]
showtopfc(b.train,"koGata3")[1:40,]
showtopfc(b.train,"koGata3_5_7")[1:40,]
showtopfc(b.train,"koNhedc2")[1:40,]
showtopfc(b.train,"koEtv2")[1:40,]
showtopfc(b.train,"koB230219D22Rik")[1:40,]


v<-showtopfc(b.train,"koXbp1")#[1:40,]
v <- v[v$mgi_symbol %in% togenesym2(c_alltime$gene[c_alltime$tf=="chip_Xbp1"]),]
v[1:100,]

v<-showtopfc(b.train,"koGata3")#[1:40,]
v <- v[v$mgi_symbol %in% togenesym2(c_alltime$gene[c_alltime$tf=="chip_Gata3"]),]
v[1:100,]   #Gata3 itself comes high

#apply(dat[grep("ERCC",rownames(dat)),],2,sum)  #no ERCCs anyway. but should remove these

###############################################################
############# where genes are different? ######################
###############################################################

### Normalize. first for some of these being off-center, then for averages
b.train.center <- b.train
for(i in 1:ncol(b.train)){
  b.train.center[,i] <- b.train.center[,i] - mean(b.train.center[,i])
}
# ov <- apply(b.train.center[,grep("ko",colnames(b.train))],1,mean)
# for(i in 1:ncol(b.train)){
#   b.train.center[,i] <- b.train.center[,i] - ov  #median(b.train.center[,i])
# }

v<-b.train.center[toensid("Il13"),]; round(v[order(abs(v),decreasing = TRUE)],digits = 4)
v<-b.train.center[toensid("Il4"),]; round(v[order(abs(v),decreasing = TRUE)],digits = 4)
v<-b.train.center[toensid("Gata3"),]; round(v[order(abs(v),decreasing = TRUE)],digits = 4)


sort(apply(b.train,2,sd))

#il13 heavily affected by depth(?). il4?

#############################################
############# check go ######################
#############################################
gettopfc <- function(col,n=100){
  v <- showtopfc(b.train,col,doplot = FALSE)
  v <- v[order(v$fc),]
  v2 <- c(v$mgi_symbol[1:(n/2)], reverse(v$mgi_symbol)[1:(n/2)])
  v2  
}
testtopgo <- function(col,n=100){
  v <- showtopfc(b.train,col)
  v <- v[order(v$fc),]
  v2 <- c(v$mgi_symbol[1:(n/2)], reverse(v$mgi_symbol)[1:(n/2)])
  
  stopgosym(
    v2,
    togenesym2(rownames(b.train))  
  )[1:60,c(2,6)]
}
gon <- 500
testtopgo("bfp1",gon)
testtopgo("bfp",gon)
testtopgo("koIfngr1",gon)
testtopgo("koIl4",gon)
testtopgo("koIl13",gon) #cyto prod!
testtopgo("koBcl11b",gon)
testtopgo("koCxcr7",gon) #reg cyto. MAPK. cytokine production.  many nice terms
testtopgo("koGata3",gon)
testtopgo("koXbp1",gon) #immune process
testtopgo("koBcl11b",gon) #cata. MAPK. nice
testtopgo("koB230219D22Rik",gon) #lot of aggregation. catabolic. 

testtopgo("koGata3_5_7",gon)


###############################
###### cluster individual probes?
##########################################

toclust <- b.train
# for(i in 1:ncol(toclust))
#   toclust[,i] <- toclust[,i]/sd(toclust[,i])  #might not be a good idea
d <- stats::dist(t((toclust)))
set.seed(0) 
rtsne_out <- Rtsne(d,is_distance = TRUE, perplexity=30, verbose = TRUE, max_iter=5000,dims = 2)
#rcol <- colorbyko()
#rcol <- colorbymouse()
#rcol <- colorbywt()
#plot(rtsne_out$Y, col=rcol[toycellcondition$isgood], pch=16, main='',xlab="",ylab="",cex=0.8)
plot(rtsne_out$Y, pch=16, main='',xlab="",ylab="",cex=0)  # col=rcol[toycellcondition$isgood], 
text(rtsne_out$Y[,1], rtsne_out$Y[,2], 
     #col=rcol[toycellcondition$isgood],
     labels=colnames(b.train),cex=0.5) 
heatmap(cor(toclust))


###############################
#### Consensus DE
###############################

getkoconcensus <- function(listko,useall=TRUE,n=500){
  v <- apply(b.train[,listko,drop=FALSE],1,mean)
  if(useall){
    v <- v - apply(b.train[,grep("ko",colnames(b.train))],1,mean)
  }
  hist(as.double(v),breaks=100)
  d <- data.frame(mgi_symbol=togenesym2(rownames(b.train)), fc=v, stringsAsFactors = FALSE)
  d <- d[order(d$fc,decreasing = TRUE),]
  d[c(1:(n/2),  nrow(d) - (1:(n/2)) +1 ),]
}

v<-getkoconcensus(c("koStat6_11_4","koStat6_11_5"))

v<-getkoconcensus(c("koXbp1_12_4","koXbp1_12_5"), useall = FALSE)
v2<-getkoconcensus(c("koXbp1_12_4","koXbp1_12_5"),n=1000, useall = FALSE)
v2[v2$mgi_symbol %in% dgbix$Gene.Name[dgbix$Xbp1_72h>0],]

v<-getkoconcensus(c("koErn1_4_4")) #one grna failed

v2<-getkoconcensus(c("koXbp1_12_4","koXbp1_12_5", "koErn1_4_4")) #if we believe ern1 & xbp1 has similar effect:
v

stopgosym(v$mgi_symbol, togenesym2(rownames(b.train)))



# v<-intersect(
#   gettopfc("koStat6_11_4",n=2000),
#   gettopfc("koStat6_11_5",n=2000))
# 
# v<-intersect(
#   gettopfc("koBcl11b_1_7",n=2000),
#   gettopfc("koBcl11b_1_8",n=2000))
# 
# v<-intersect(
#   gettopfc("koBhlhe40_2_1",n=2000),
#   gettopfc("koBhlhe40_2_2",n=2000))
# 
# v<-intersect(
#   gettopfc("koCcdc134_2_3",n=2000),  #this ignores sign - later problem
#   gettopfc("koCcdc134_2_4",n=2000))
# 
# stopgosym(v, togenesym2(rownames(b.train)))
#gettopfc("koCcdc134_2_3",n=100)

# hist(b.train$koGata3,breaks=100)
# hist(b.train$koIl13,breaks=100)
# hist(b.train$koIl4,breaks=100)


head(b.train)

usealpha <- 0.5
fit.elnet <- glmnet(tx, ty,  alpha=usealpha, family="gaussian")  
plot(fit.elnet, xvar="lambda")  #this only makes sense for a few param
#rownames(b.train) <- colnames(x.train)  

fit.elnet
#what is that extra number??? - it is intercept. needed? already in model. kick out?

head(coef(fit.elnet)[, 10])

### Cross-validation
cv.glmmod <- cv.glmnet(tx, ty, alpha=usealpha)   #Note. can make parallel. 
best.lambda <- cv.glmmod$lambda.min
best.lambda
plot(cv.glmmod)


b.train <- untransformB(y.train, coef(fit.elnet)[, 10])


dim(y.train)
dim(x.train)
length(coef(fit.elnet)[, 10])
dim(b.train)

length(ty)
dim(tx)


############################################################################################################################
################# What if we cannot trust DEseq pvalues because too uncertain, and only rely on FC? ########################
############################################################################################################################

plotresvolcano <- function(gene)
  plot(resfc[,gene],-log10(respval[,gene]))
plotresvolcano("Tbx21")
plotresvolcano("Gata3")
plotresvolcano("Il4")
plotresvolcano("Il13")

