sumcolpairs <- function(x){
  y <- matrix(0,ncol=ncol(x)/2,nrow=nrow(x))
  for(i in 1:ncol(y)){
    y[,i] <- x[,i*2-1]+x[,i*2]
  }
  y
}

### Read peak annotation
newatac_ann <- read.csv("atac/mouse/ATACall_peaks.red.ann.csv",sep="\t",stringsAsFactors = FALSE)
colnames(newatac_ann)[1] <- "peakid" #or something
head(newatac_ann)

### Read background counts and figure out average counts
newatac_inv <- read.table("atac/mouse/ATACall_peaks.inv.bed",sep="\t",stringsAsFactors = FALSE)
newatac_inv_sum <- sum((newatac_inv$V3-newatac_inv$V2))
newatac_bg <- read.table("atac/mouse/counts.f.bg.csv",sep="\t",stringsAsFactors = FALSE)
newatac_bg <- apply(sumcolpairs(newatac_bg[,-(1:3)]),2,sum)
newatac_bg_avgreads <- newatac_bg/newatac_inv_sum

### Read peak counts
newatac_peaks <- read.table("atac/mouse/counts.f.peaks.csv",sep="\t",stringsAsFactors = FALSE)
newatac_peaks <- cbind(newatac_peaks[,4,drop=FALSE],sumcolpairs(newatac_peaks[,-(1:6)]))
colnames(newatac_peaks) <- c("peakid","Naive","Th2_2h","Th2_4h","Th2_24h","Th2_48h","Th2_72h")  #consider other parts of the code
newatac_peaks <- smerge(newatac_peaks,newatac_ann)
#newatac_peaks <- newatac_peaks[order(newatac_peaks$Th2_4h,decreasing = TRUE),]
newatac_peaks <- newatac_peaks[order(newatac_peaks$Th2_72h,decreasing = TRUE),]
newatac_peaks <- newatac_peaks[order(newatac_peaks$Th2_4h/newatac_peaks$Th2_72h,decreasing = TRUE),]
newatac_peaks <- newatac_peaks[order(apply(newatac_peaks[,1+(1:6)],1,sum),decreasing = TRUE),]
newatac_peaks$rcrank <- 1:nrow(newatac_peaks)

#Normalize peaks by background and length
newatac_peaks_norm <- newatac_peaks
newatac_peakslen <- (newatac_peaks$End-newatac_peaks$Start)
for(i in 1:6){
  #Arbitrary unit to make it easier to think. now most peaks in 0-2. with up to 10
  newatac_peaks_norm[,i+1] <- 1e5*(newatac_peaks[,i+1]-newatac_bg_avgreads[i])/newatac_bg[i]/newatac_peakslen
}




#head(newatac_peaks$Detailed.Annotation,n=100)
#t(newatac_peaks[47,])

# plot(as.double(newatac_peaks[1,1+(1:6)]),type="l",ylim=c)
# for(i in 1:100){
#   lines(as.double(newatac_peaks[i,1+(1:6)]),type="l")
# }

#grep("Satellite",newatac_peaks$Detailed.Annotation)
#newatac_peaks <- newatac_peaks[-grep("Satellite",newatac_peaks$Detailed.Annotation),]




###### scaling law over time?
h<-hist(log(apply(newatac_peaks[,1+(1:6)],1,sum)),breaks=100)
plot(h$mids,log(h$density),type="l")

###### scaling law over time?  this with log-y
h<-hist(log(newatac_peaks_norm[,1+1]),breaks=100,plot = FALSE)
plot(h$mids,log(h$density),type="l",xlim=c(-13,0)) #naive
h<-hist(log(newatac_peaks_norm[,2+1]),breaks=100,plot = FALSE)
lines(h$mids,log(h$density),type="l",xlim=c(-13,0)) #2h
h<-hist(log(newatac_peaks_norm[,3+1]),breaks=100,plot = FALSE)
lines(h$mids,log(h$density),type="l",xlim=c(-13,0)) #4h
h<-hist(log(newatac_peaks_norm[,4+1]),breaks=100,plot = FALSE)  ##Suddenly all low-expressing peaks went up 4-fold
lines(h$mids,log(h$density),type="l",xlim=c(-13,0)) #24h
h<-hist(log(newatac_peaks_norm[,5+1]),breaks=100,plot = FALSE) #slight broadening but same trend
lines(h$mids,log(h$density),type="l",xlim=c(-13,0)) #48h
h<-hist(log(newatac_peaks_norm[,6+1]),breaks=100,plot = FALSE) #here everything went down except the mid-population
lines(h$mids,log(h$density),type="l",xlim=c(-13,0))


###### scaling law over time?  this w/o log-y
h<-hist(log(newatac_peaks_norm[,1+1]),breaks=100,plot = FALSE)
plot(h$mids,(h$density),type="l",xlim=c(-13,0)) #naive
h<-hist(log(newatac_peaks_norm[,2+1]),breaks=100,plot = FALSE)
lines(h$mids,(h$density),type="l",xlim=c(-13,0)) #2h
h<-hist(log(newatac_peaks_norm[,3+1]),breaks=100,plot = FALSE)
lines(h$mids,(h$density),type="l",xlim=c(-13,0)) #4h
h<-hist(log(newatac_peaks_norm[,4+1]),breaks=100,plot = FALSE)  ##Suddenly all low-expressing peaks went up 4-fold
lines(h$mids,(h$density),type="l",xlim=c(-13,0)) #24h
h<-hist(log(newatac_peaks_norm[,5+1]),breaks=100,plot = FALSE) #slight broadening but same trend
lines(h$mids,(h$density),type="l",xlim=c(-13,0)) #48h
h<-hist(log(newatac_peaks_norm[,6+1]),breaks=100,plot = FALSE) #here everything went down except the mid-population
lines(h$mids,(h$density),type="l",xlim=c(-13,0))




h<-hist(log(apply(newatac_peaks[,2+(1:6)],1,sum)),breaks=100)
plot(h$mids,log(h$density),type="l",xlim=c(-13,0))




plot(as.double(newatac_peaks_norm[1001,1+(1:6)]),type="l")
plot(apply(newatac_peaks_norm[1:500,1+(1:6)],2,mean),type="l")
plot(apply(newatac_peaks_norm[1:500,1+(1:6)],2,mean),type="l")

plot(as.double(newatac_peaks_norm[70000,1+(1:6)]),type="l")
plot(as.double(newatac_peaks_norm[1,1+(1:6)]),type="l",ylim=c(0,2))
plot(as.double(newatac_peaks_norm[1,1+(1:6)]),type="l",ylim=c(0,0.1))
for(i in 2:100){
  lines(as.double(newatac_peaks_norm[i,1+(1:6)]),type="l")
}
#### due to the unfolding, some regions become hyperaccessible. These are then taken down to normal levels. Which ones are affected?



#hist(newatac_peaks[newatac_peaks[,7]>1e-13,7])
# hist(newatac_peaks_norm[newatac_peaks_norm[,2]>0.1,2],breaks = 20)
# hist(newatac_peaks_norm[newatac_peaks_norm[,3]>0.1,3],breaks = 20)
# hist(newatac_peaks_norm[newatac_peaks_norm[,4]>0.1,4],breaks = 20)
# hist(newatac_peaks_norm[newatac_peaks_norm[,5]>0.1,5],breaks = 20)
# hist(newatac_peaks_norm[newatac_peaks_norm[,6]>0.1,6],breaks = 20)

###Now much more clear: eno
# plot( apply(newatac_peaks_norm[,1+(1:6)],2,function(x) quantile(x,0.9)),type="l",ylim=c(0,1e-3))
# lines(apply(newatac_peaks_norm[,1+(1:6)],2,function(x) quantile(x,0.7)),col="blue")
# lines(apply(newatac_peaks_norm[,1+(1:6)],2,function(x) quantile(x,0.5)),col="blue")
# lines(apply(newatac_peaks_norm[,1+(1:6)],2,function(x) quantile(x,0.3)),col="blue")
# lines(apply(newatac_peaks_norm[,1+(1:6)],2,function(x) quantile(x,0.1)),col="blue")


colnames(newatac_peaks)


newatac_peaks_sum <- apply(newatac_peaks[,-(1)],2,sum)
barplot(newatac_bg)
barplot(newatac_peaks_sum)
barplot(newatac_peaks_sum/newatac_bg)

### The peaks with most reads recondense. The peaks with fewer reads however do not.
#if we instead used the peaks with few reads as reference then the large ones would decrease even more!
newatac_peaks_norm_sum <- apply(newatac_peaks_norm[(1:1000),1+(1:6)],2,sum)
barplot(newatac_peaks_norm_sum)
newatac_peaks_norm_sum <- apply(newatac_peaks_norm[-(1:35000),1+(1:6)],2,sum)
barplot(newatac_peaks_norm_sum)
newatac_peaks_norm_sum <- apply(newatac_peaks_norm[-(33000:37000),1+(1:6)],2,sum)
barplot(newatac_peaks_norm_sum)
newatac_peaks_norm_sum <- apply(newatac_peaks_norm[-(1:50000),1+(1:6)],2,sum)
barplot(newatac_peaks_norm_sum)
newatac_peaks_norm_sum <- apply(newatac_peaks_norm[-(1:70000),1+(1:6)],2,sum)
barplot(newatac_peaks_norm_sum)

###Categorize the peaks around IL4. These are index 35000+, and average of these is decreasing a bit
newatac_peaks_norm[newatac_peaks_norm$peakid %in% c("ATACall_peak_11314","ATACall_peak_11312"),]

newatac_peaks


sqldf("select median(rcrank) from newatac_peaks_norm group by Chr")
sqldf("select avg(rcrank) from newatac_peaks_norm group by Chr")
#Chromesome 21 has a much lower mean and median rank! also the smallest chromosome. suggesting that this one is turn down in favour
#of others. https://en.wikipedia.org/wiki/Category:Genes_on_human_chromosome_21 
#Appears to have many of the immune genes
#Infgr2, Il10 receptor, runx1, Aire, 
#Chromosome 20 has the highest mean/median.
#Nfatc2. 
#Chromosome 8 and 17 among the lowest otherwise. 17 has the Ccl, Fli1, 

#### Read chromosome sizes
chrsize <- read.table("atac/mm10.chrom.sizes",sep="\t", stringsAsFactors = FALSE)
colnames(chrsize)<-c("Chr","Chrsize")
rownames(chrsize) <- chrsize$Chr
chrsize <- chrsize[order(chrsize$Chr),]
chrsize <- chrsize[-grep("_",chrsize$Chr),]
normalizebychrsize <- function(x){
  round(x/chrsize[names(x),]$Chrsize*1e8,digits = 3)
}

#For sorting by average size
normalizebychrsize(table(newatac_peaks_norm$Chr[1:100])) #chr9 has many annoying. nfil3 jak mtorc il33 here. 17 relatively many
normalizebychrsize(table(newatac_peaks_norm$Chr[1000:2000])) 

sort(table(newatac_peaks_norm$Chr[40000:80000])/table(newatac_peaks_norm$Chr[1000:20000]))
#table(newatac_peaks_norm$Chr[100:20000])

#Sorting by final size
normalizebychrsize(table(newatac_peaks_norm$Chr[1:100])) #Chromosome 9 has the biggest increases in openness. then 17, 7 in general
normalizebychrsize(table(newatac_peaks_norm$Chr[1:1000])) #9 second largest
normalizebychrsize(table(newatac_peaks_norm$Chr[1:10000])) #Y extremely small
normalizebychrsize(table(newatac_peaks_norm$Chr[1:30000])) #Y extremely small, X larger. otherwise even
normalizebychrsize(table(newatac_peaks_norm$Chr[1:70000])) #Y still very small and X not catching up. 11 most open
#1,9,16 supposed to have large chunks of heterochromatin https://en.wikipedia.org/wiki/Heterochromatin


#Easiest is to add up the replicates

head(newatac_peaks)





peakpat <- read.csv("out_tc/fimo.txt",stringsAsFactors = FALSE,sep="\t")
head(peakpat)
colnames(peakpat)[1]<-"motifid" ###there is the TF name there straight!
#colnames(peakpat)[2]<-"peak"
colnames(peakpat)[3]<-"peakrename"
hist(log(peakpat$p.value))


newatac_peaks_norm$peakrename <- sprintf("%s:%s-%s",newatac_peaks_norm$Chr, newatac_peaks_norm$Start, newatac_peaks_norm$End)

#chr7:25120342-25120810


