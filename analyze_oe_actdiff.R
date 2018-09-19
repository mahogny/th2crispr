
######################################################################
### OE activation/differentiation ####################################
######################################################################


####### First perform the DE

##Normalization of strength between the OE since we have a FC cutoff
sdNewkoFC <- apply(newkoFC,1,sd)
oeFCnorm <- oe_de$fc

### by using DE fold changes
compareDEwithRef_oe <- function(thep, thefc){
  theall<-c()
  allsd<-c()
  for(j in 1:ncol(thep)){
    numgenetested <- c()
    print(sprintf("ref: %s",j))
    #Only consider DE genes according to reference
    totest <- !is.na(thep[,j]) & thep[,j]<1e-10     
    totest <- intersect(rownames(thefc)[totest],rownames(oe_de$padj)) 
    
    
    cmin<-function(x) sum(na.omit(x)==-1)
    cplus<-function(x) sum(na.omit(x)==1)
    kov <- c()
    for(i in 1:ncol(oeFCnorm)){ 
      #Further only consider DE genes according to our dataset. Use better cut-off
      usefc <- oeFCnorm
      totest3 <- intersect(
        totest,
        rownames(usefc)[abs(usefc[,i])>log(2,2)])   
      
      
      numgenetested <- c(numgenetested,length(totest3))
      s1<-sign(thefc[totest3,j])
      s2<-sign(usefc[totest3,i])
      a<-sum(s1*s2)
      
      kov <- c(kov, sum(a))
    }
    print(numgenetested)
    if(j==1)  
      theall<-kov
    else
      theall<-rbind(theall,kov)
  }
  colnames(theall) <- colnames(oeFCnorm)
  rownames(theall) <- sprintf("Th2 vs %s",colnames(thep))
  hist(numgenetested,breaks=20)
  
  thenorm <- theall
  for(i in 1:nrow(thenorm)){
    thenorm[i,] <- thenorm[i,]/sd(thenorm[i,])
  }
  
  list(
    theall=theall,
    thenorm=thenorm)
  #  theall#list(theall=theall)  
}




###Could generalize this function to plasticity analysis too!
compoeTcAcDiff <- compareDEwithRef_oe(resTcActDiffP, resTcActDiffFC)


keep <- c("Lrrc40","Pparg","Scara3","B230219D22Rik", "Ccdc134", "Bhlhe40")  

inprev <- c("Gata3", "Lrrc40","Ccdc134","Scara3","Bhlhe40")

vis <- -compoeTcAcDiff$theall[,keep]

thecol <- rep("black",ncol(vis))
thecol[colnames(vis) %in% inprev]    <- "#c80000ff"

pdf("out_oe/oe_actdiff.pdf")
plot(vis[2,], vis[1,],cex=0,
     ylab="Activation related", xlab="Differentiation related",
     xlim=c(-12,12),ylim=c(-10,10))
lines(minmax(vis[2,]),c(0,0),col="gray")
lines(c(0,0),minmax(vis[1,]),col="gray")
text(vis[2,], vis[1,],col=thecol,
     labels = colnames(vis))
dev.off()




######################################################################
### KO similarity to Th_x - classification ###########################
######################################################################

##Normalization of strength between the KO since we have an FC cutoff
oeFCnorm <- -oe_de$fc


### by using DE fold changes
theall<-c()
allsd<-c()
numgenetested <- c()
for(j in 1:ncol(thep)){
  print(sprintf("Thx: %s",j))
  #Only consider DE genes according to ThExpress
  totest <- !is.na(thep[,j]) & thep[,j]<1e-2              ################# e-2
  totest <- intersect(rownames(thefc)[totest],rownames(newkoPadj))
  
  cmin<-function(x) sum(na.omit(x)==-1)
  cplus<-function(x) sum(na.omit(x)==1)
  kov <- c()
  for(i in 1:ncol(oeFCnorm)){
    #Further only consider DE genes according to the KO
    usefc <- oeFCnorm
    totest3 <- intersect(
      totest,
      rownames(usefc)[abs(usefc[,i])>log(0.5,2)])                    #2,2
    
    numgenetested <- c(numgenetested,length(totest3))
    s1<-sign(thefc[totest3,j])
    s2<-sign(usefc[totest3,i])
    a<-sum(s1*s2)
    
    kov <- c(kov, a/length(totest)) 
  }
  if(j==1)  
    theall<-kov
  else
    theall<-rbind(theall,kov)
}
colnames(theall) <- colnames(oeFCnorm)
rownames(theall) <- sprintf("Th2 vs %s",colnames(thep))

thenorm <- theall
for(i in 1:nrow(thenorm)){
  thenorm[i,] <- thenorm[i,]/sd(thenorm[i,])
}
colnames(thenorm) <- colnames(theall)


thenorm <- thenorm[,-grep("oe",colnames(thenorm))]

pdf("out_oe/oe_thtype_norm_new.pdf",width = 10, height = 2.5)
thecol <- colorRampPalette(c("#dd0152", "white","#2719dd"))(n = 299) #red to blue
heatmap.2(thenorm[-1,],
          trace="none",
          #density.info="none", 
          dendrogram = "none",
          cexRow = 0.8,
          col=thecol,
          Rowv=FALSE)#,
dev.off()



######################################################################
### KO similarity to Th_x - pinwheel #################################
######################################################################


pdf("out_oe/oe_plast_pinwheel_neutral_new.pdf")
fortsne <- t(thenorm)[,-1]
colnames(fortsne) <- c("Th1","Th17","iTreg","nTreg")
npin <- ncol(fortsne)+1  #last one is Th2
pinx <- cos((1:npin)*2*pi/npin - 5*pi/5/2)
piny <- sin((1:npin)*2*pi/npin - 5*pi/5/2)
plot(c(-2,2),c(-2,2),cex=0)
#lines(c(pinx,pinx[1]),c(piny,piny[1]), col="gray")
for(i in 1:npin){ #nrow(fortsne)
  lines(c(0,pinx[i]),c(0,piny[i]), col="gray")
}
text(pinx*1.1,piny*1.1, labels = c(colnames(fortsne),"Th2"),cex=0.8)
text(0,0, labels = "Neutral")
for(i in 1:npin){ 
  pinx[i] <- pinx[i] - pinx[npin]
  piny[i] <- piny[i] - piny[npin]
}
for(i in 1:nrow(fortsne)){
  sc<-0.1  #for the normalized
  #  sc<-0.08  #for the normalized
  thecol <- "red"
  if(rownames(fortsne)[i] %in% c("Cd180","Cd200","Cxcr7","F2rl1","Ifngr1","Lag3","Orm3","Scara3","Thy1"))
    thecol <- "#c80000ff"
  if(rownames(fortsne)[i] %in% c("Bcl11b","Bhlhe40","Ern1","Etv2","Gata3","Stat6","Tbx21","Xbp1","Zc3h12a"))
    thecol <- "#0000c8ff"
  if(rownames(fortsne)[i] %in% c("Abcg4","Crls1","Nhedc2","Pgk1","Pxk","Slc5a1","Slc25a3"))
    thecol <- "#b900c8ff"
  if(rownames(fortsne)[i] %in% c("Ccdc134","Il2","Il4","Il12a","Il13"))
    thecol <- "#00c800ff"
  v <- c(fortsne[i,],0)
  
  text(sc*sum(v*pinx), sc*sum(v*piny), rownames(fortsne)[i],cex=0.5, col=thecol)
}
dev.off()




