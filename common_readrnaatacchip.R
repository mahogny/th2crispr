#############################################################################################
###                                                                                       ###
###             Part of the paper ...                                                     ###
###             Author: Johan Henriksson (mahogny@areta.org)                              ###
###                                                                                       ###
###             This code ...                                   ###
###                                                                                       ###
#############################################################################################

#note: stat6 read depth is fine. but stat6 72h b is the wrong file!


library(Rtsne)
library(gplots)
library(RColorBrewer)
library(stringr)
library(sqldf)
library(reshape2)
library(limma)
library(GenomicRanges)
#library(BiocParallel)
#register(MulticoreParam(4))


## Cut-offs used. TODO move ATAC cutoff here
maxTSSdistCHIP <- 20e3


##################################################################
####### Common helper functions ##################################
##################################################################

######### Clean up memory
showmemuse <- function(){
  for (thing in ls()) {
    s <- object.size(
      get(thing)
    )
    if(s>20000000){
      cat(sprintf("%s\t%s\n", round(s/1000000),thing))  #in MB
    }
  }
}
#showmemuse()


## Not in-operator
'%!in%' <- function(x,y) !('%in%'(x,y))

## Rbind the elements of a list
rbindlist <- function(thelist){
  as.data.frame(data.table::rbindlist(thelist))
}

qtextscatter <- function(x,y,labels,cex=1){
  plot(x,y,cex=0, xlab=deparse(substitute(x)), ylab=deparse(substitute(y)))
  text(x,y,labels = labels,cex=cex)
}


detach_package <- function(pkg, character.only = FALSE){
  if(!character.only){
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while(search_item %in% search())
  {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}


bpsapply <- function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
  FUN <- match.fun(FUN)
  answer <- bplapply(X = X, FUN = FUN, ...)
  if (USE.NAMES && is.character(X) && is.null(names(answer))) 
    names(answer) <- X
  if (!identical(simplify, FALSE) && length(answer)) 
    simplify2array(answer, higher = (simplify == "array"))
  else answer
}


mergebyrow <- function(x,y){
  merge(x,y, by="row.names")
}

minmax <- function(x){
  c(min(x),max(x))
}
symrange <- function(x){
  c(-max(abs(x)),max(abs(x)))
}

normalizesym <- function(s) paste(str_to_upper(str_sub(s,1,1)),str_to_lower(str_sub(s,2)),sep="")

#########
## Function to calculate a "correlation matrix" of jaccard indices
corjaccard <- function(vd){
  v <- matrix(NA,ncol(vd),ncol(vd))
  colnames(v)<-colnames(vd)
  rownames(v)<-colnames(vd)
  for(i in 1:ncol(vd))
    for(j in 1:ncol(vd)){
      v[i,j] <- sum(vd[,i]>0 & vd[,j]>0) / sum(vd[,i]>0 | vd[,j]>0)
    }
  v
}

vecjaccard <- function(vi,vj){
  sum(vi>0 & vj>0) / sum(vi>0 | vj>0)
}

#########
### Merge columns with the same name by taking max value
mergecolmax <- function(x){
  cn <- unique(colnames(x))
  nx <- matrix(0,nrow = nrow(x), ncol=length(cn))
  colnames(nx) <- cn
  rownames(nx) <- rownames(x)
  for(i in 1:length(cn)){
    nx[,i] <- apply(x[,cn[i],drop=FALSE],1,max)
  }
  nx  
}

#########
### Return object but with class set to double
as.class.double <- function(x){
  class(x) <- "double"
  x
}

#########
### Function: geometric mean
gm_mean <- function(a){prod(a)^(1/length(a))}

#########
### Function: Scale values to a range to 0-1
scale01 <- function(x){
  x<-x-min(x)
  x<-x/max(x)
  x
}

#########
### Function: Turn a boolean matrix into a 1/0 matrix 
binarize <- function(m2){
  m2[m2]<-1
  m2[!m2]<-0
  m2
}


#########
### Safe merge: like merge() but first checks that there is at least one common column
smerge <- function(x, y, by = intersect(names(x), names(y)),
      by.x = by, by.y = by, all = FALSE, all.x = all, all.y = all,
      sort = TRUE, suffixes = c(".x",".y"),
      incomparables = NULL){
  if(length(by)==0){
    print(colnames(x))
    print(colnames(y))
    stop("No overlap between tables")
  } else {
    merge(x, y, by,
          by.x, by.y, all, all.x , all.y,
          sort , suffixes ,
          incomparables )
  }
}



keepscreens <- c("s8a_stl", "sx2_stl","first_il4",
                 "s11_il13","sc1_il13",
                 "sx4_irf4","sc2a_irf4",
                 "s8b_xbp1","sc2b_xbp1",
                 "sc3_gata3","s9_stg") 
keepscreens_ren2 <- c("IL4 a", "IL4 b", "IL4 c",   
                      "IL13 a","IL13 b",
                      "Irf4 a","Irf4 b",   
                      "Xbp1 a","Xbp1 b",    
                      "Gata3 a","Gata3 b") 

screens_il4 <- c("s8a_stl", "sx2_stl", "first_il4")
screens_il13 <- c("s11_il13","sc1_il13")
screens_irf4 <- c("sx4_irf4","sc2a_irf4")
screens_xbp1 <- c("s8b_xbp1","sc2b_xbp1")
screens_gata3 <- c("sc3_gata3","s9_stg") 

list_screen_genes <- c("Il4","Il13","Irf4","Xbp1","Gata3")
list_screens <- list(il4=screens_il4, il13=screens_il13, irf4=screens_irf4, xbp1=screens_xbp1, gata3=screens_gata3)



##################################################################
## Read ATAC motifs
atactf <- read.csv("out_motif/atactf.csv",stringsAsFactors=FALSE)
for(i in 1:nrow(atactf))
  atactf$motif[i] <- normalizesym(atactf$motif[i])
colnames(atactf)<-c("sym","p")
atactf$sym[atactf$sym=="Bhlh2b"] <- "Bhlhe40"
atactf$sym[atactf$sym=="Bhlh3b"] <- "Bhlhe41"
atactf$sym[atactf$sym=="Hinfp1"] <- "Hinfp"
#Tcfap2a -> ??? Ap2
#"Zbed1" -> ???


##################################################################
## Read TF name <-> TF ID
read.jaspar_namegenesym <- function(){
  #Read orthology map. Only consider unique mappings human->mouse
  map_ortho_humanmouse <- read.csv("map_ortho_humanmouse.csv",stringsAsFactors = FALSE)
  map_ortho_humanmouse <- map_ortho_humanmouse[!duplicated(map_ortho_humanmouse$human),]
  #map_ortho_humanmouse <- map_ortho_humanmouse[!duplicated(map_ortho_humanmouse$mouse),]
  rownames(map_ortho_humanmouse) <- str_to_lower(map_ortho_humanmouse$human)
  
  #Read multimap from jasparname to several genes involved.
  #Remap human gene names to mouse gene names
  map_jaspar_namegenesym <- read.csv("map_jasparname_sym.csv",stringsAsFactors = FALSE)
  altmap <- map_ortho_humanmouse[str_to_lower(map_jaspar_namegenesym$mgi_symbol),]$mouse
  ismouse <- map_jaspar_namegenesym$mgi_symbol %in% ensconvert$mgi_symbol
  map_jaspar_namegenesym$mgi_symbol[!ismouse] <- altmap[!ismouse]
  map_jaspar_namegenesym <- map_jaspar_namegenesym[!is.na(map_jaspar_namegenesym$mgi_symbol),]
  map_jaspar_namegenesym
}
map_jaspar_namegenesym <- read.jaspar_namegenesym()


###########################################################
############# read atac and rnaseq data ###################
###########################################################


##################################################################
### Read mapping TF name <-> motifID
genemotif <- read.csv("out_tc/JASPAR2016_MA_PB_C2H2_nonred.meme.names",header=FALSE,sep=" ",stringsAsFactors = FALSE)[,c(2,3)]
colnames(genemotif) <- c("motifid","jasparname")   #changed
normalizesym <- function(s) paste(str_sub(s,1,1),str_to_lower(str_sub(s,2)),sep="")
for(i in 1:nrow(genemotif))
  genemotif$tf[i] <- normalizesym(genemotif$tf[i])




##################################################################
### Read RNAseq time course
read.mtpm <- function(org, ensconvert_){
  # org <- "mouse"
  # ensconvert_=ensconvert
  
  mtpm <- read.csv(sprintf("out_tc/%s/tpm.txt",org),sep="\t",row.names = "gene")
  mtpm <- mtpm[,colnames(mtpm)!="row.names"]
  colnames(mtpm) <- str_replace_all(colnames(mtpm),"rep","")
  
  
  ### Calculate average TPM over time, Th2
  mtpm_times <- c(0,0.5,1,2,4,6,12,24,48,72)
  av_mtpm <- cbind(  #Th2
    apply(mtpm[,grep("Naive",colnames(mtpm))],1,mean),
    apply(mtpm[,grep("Th2_05h",colnames(mtpm))],1,mean),
    apply(mtpm[,grep("Th2_1h",colnames(mtpm))],1,mean),
    apply(mtpm[,grep("Th2_2h",colnames(mtpm))],1,mean),
    apply(mtpm[,grep("Th2_4h",colnames(mtpm))],1,mean),
    apply(mtpm[,grep("Th2_6h",colnames(mtpm))],1,mean),
    apply(mtpm[,grep("Th2_12h",colnames(mtpm))],1,mean),
    apply(mtpm[,grep("Th2_24h",colnames(mtpm))],1,mean),
    apply(mtpm[,grep("Th2_48h",colnames(mtpm))],1,mean),
    apply(mtpm[,grep("Th2_72h",colnames(mtpm))],1,mean)
  )
  colnames(av_mtpm) <- sprintf("%sh",mtpm_times)
  
  ### Calculate average TPM over time, Th0
  av_mtpm0 <- cbind(
    apply(mtpm[,grep("Naive",colnames(mtpm))],1,mean),
    apply(mtpm[,grep("Th0_05h",colnames(mtpm))],1,mean),
    apply(mtpm[,grep("Th0_1h",colnames(mtpm))],1,mean),
    apply(mtpm[,grep("Th0_2h",colnames(mtpm))],1,mean),
    apply(mtpm[,grep("Th0_4h",colnames(mtpm))],1,mean),
    apply(mtpm[,grep("Th0_6h",colnames(mtpm))],1,mean),
    apply(mtpm[,grep("Th0_12h",colnames(mtpm))],1,mean),
    apply(mtpm[,grep("Th0_24h",colnames(mtpm))],1,mean),
    apply(mtpm[,grep("Th0_48h",colnames(mtpm))],1,mean),
    apply(mtpm[,grep("Th0_72h",colnames(mtpm))],1,mean)
  )
#  print(head(av_mtpm0))
  max_mtpm <- apply(av_mtpm,1,max)
  
  mtpm_th0 <- mtpm[,c(
    grep("Naive",colnames(mtpm)),
    grep("Th0_",colnames(mtpm)))]
  mtpm_th2 <- mtpm[,c(
    grep("Naive",colnames(mtpm)),
    grep("Th2_",colnames(mtpm)))]
  
  mtpm_early <- mtpm[,c(
    grep("Naive",colnames(mtpm)),
    grep("_05h",colnames(mtpm)),
    grep("_1h",colnames(mtpm)),
    grep("_2h",colnames(mtpm)),
    grep("_4h",colnames(mtpm)))]
  mtpm_late <- mtpm[,c(
    grep("_6h",colnames(mtpm)),
    grep("_12h",colnames(mtpm)),
    grep("_24h",colnames(mtpm)),
    grep("_48h",colnames(mtpm)),
    grep("_72h",colnames(mtpm)))]
  
  expressedGenes10Id <- names(av_mtpm[apply((cbind(av_mtpm,av_mtpm0))>10,1,sum)>0,1])  
  expressedGenes10   <- unique(ensconvert$mgi_symbol[ensconvert_$ensembl_gene_id %in% expressedGenes10Id])

  #Differentially expressed genes
  de_early <- read.csv(sprintf("out_tc/%s/early_DE_Th0Th2_genes.txt",org),sep="\t",stringsAsFactors = FALSE)
  de_late <- read.csv(sprintf("out_tc/%s/late_DE_Th0Th2_genes.txt",org),sep="\t",stringsAsFactors = FALSE)
  colnames(de_early)[14]<-"ensembl_gene_id"
  colnames(de_early)[15]<-"mgi_symbol"
  colnames(de_late)[14]<-"ensembl_gene_id"
  colnames(de_late)[15]<-"mgi_symbol"

  #Get expresson levels for the TFs in particular
  motif_explevel <- smerge(smerge(smerge(
    data.frame(jasparname=genemotif$tf, stringsAsFactors = FALSE), 
    map_jaspar_namegenesym),ensconvert_),
    data.frame(ensembl_gene_id=names(max_mtpm), explevel=max_mtpm, stringsAsFactors = FALSE))
  motif_explevel <- sqldf("select jasparname, mgi_symbol, ensembl_gene_id, max(explevel) as maxexp from motif_explevel group by jasparname,mgi_symbol,ensembl_gene_id")
  

  list(
    mtpm=mtpm,
    av_mtpm=av_mtpm,
    av_mtpm0=av_mtpm0,
    max_mtpm=max_mtpm,
    mtpm_th0=mtpm_th0,
    mtpm_th2=mtpm_th2,
    mtpm_early=mtpm_early,
    mtpm_late=mtpm_late,
    motif_explevel=motif_explevel,
    #expressed_atacTF=expressed_atacTF,
    expressedGenes10Id=expressedGenes10Id,
    expressedGenes10=expressedGenes10,
    de_late=de_late,
    de_early=de_early,
    ensconvert=ensconvert_
  )
}

tcmouse <- read.mtpm("mouse", ensconvert_ = ensconvert)
tchuman <- read.mtpm("human", ensconvert_ = human_ensconvert)




expressedTFjaspar <- function(tc, tpm=1){
  unique(tc$motif_explevel$jasparname[tc$motif_explevel$maxexp>tpm])
}



######################################################################
### Check which genes are DE in both mouse and human #################
######################################################################

##################################################################
## DE for human and mouse in one big matrix
getconservedDE <- function(qval=5e-2){

  getdefromtable.mouse <- function(x){
    (unique(x$ensembl_gene_id[x$qval<qval]))  #should really normalize earlier!
  }
  getdefromtable.human <- function(x){
    (unique(x$ensembl_gene_id[x$qval<qval]))  #should really normalize earlier!
  }
  allde.mouse <- smerge(data.frame(
    ens_mouse = ensconvert$ensembl_gene_id,
    me=ensconvert$ensembl_gene_id %in% getdefromtable.mouse(tcmouse$de_early),
    ml=ensconvert$ensembl_gene_id %in% getdefromtable.mouse(tcmouse$de_late), stringsAsFactors = FALSE), 
    ortho_mouse_human_unique, all.x = TRUE)
    
  allde.human <- smerge(data.frame(
    ens_human = human_ensconvert$ensembl_gene_id,
    he=human_ensconvert$ensembl_gene_id %in% getdefromtable.mouse(tchuman$de_early),
    hl=human_ensconvert$ensembl_gene_id %in% getdefromtable.mouse(tchuman$de_late), stringsAsFactors = FALSE), 
    ortho_mouse_human_unique, all.x = TRUE)
  
  allde <- smerge(allde.mouse, allde.human, all=TRUE)
  allde$me[is.na(allde$me)] <- FALSE
  allde$ml[is.na(allde$ml)] <- FALSE
  allde$he[is.na(allde$he)] <- FALSE
  allde$hl[is.na(allde$hl)] <- FALSE

  allde <- smerge(allde, data.frame(
    ens_human=human_ensconvert$ensembl_gene_id, 
    sym_human=human_ensconvert$mgi_symbol, 
    stringsAsFactors = FALSE), all.x=TRUE)
  allde <- smerge(allde, data.frame(
    ens_mouse=ensconvert$ensembl_gene_id, 
    sym_mouse=ensconvert$mgi_symbol, 
    stringsAsFactors = FALSE), all.x=TRUE)

  repNAfalse <- function(x) {
    x[is.na(x)] <- FALSE
    x
  }
  allde$anytime_mouse <- repNAfalse(allde$me | allde$ml)
  allde$anytime_human <- repNAfalse(allde$he | allde$hl)
  allde$conserved_loose   <- repNAfalse((allde$me | allde$he) | (allde$ml | allde$hl))
  allde$conserved_alltime <- repNAfalse((allde$me & allde$he) & (allde$ml & allde$hl))
  allde$conserved_special <- repNAfalse(((allde$he | allde$hl) & allde$me)  | allde$ml)
  allde$conserved_anytime <- repNAfalse((allde$he & allde$me) | (allde$hl & allde$ml))
  
  allde                               
}
allde <- getconservedDE()

na.omit(allde[allde$sym_mouse %in% c("Il4","Fli1"),])


##################################################################
## Output data to put in Venn diagram (made manually)
if(FALSE){

  x <- getconservedDE(0.001)
  c( sum(x$me),        sum(x$me & x$ml),               sum(x$ml))
  c( sum(x$me & x$he), sum(x$ml & x$hl & x$me & x$he), sum(x$ml & x$hl))
  c( sum(x$he),        sum(x$he & x$hl),               sum(x$hl))

  sum(  c( sum(x$me & x$he), sum(x$ml & x$hl & x$me & x$he), sum(x$ml & x$hl)))
  
  x[x$conserved_alltime,]$sym_mouse   #Gata3, il2rb, mapkapk3 etc
  
  #How many % of the DE genes are in at least one of our screens?
  hitsinanyscreen <- names(which(sgenescorer2_matrix[,1]<500 | apply(sgenescorer2_matrix[,-1]<1000,1,any)))
  mean(unique(x$sym_mouse[x$anytime_mouse]) %in% hitsinanyscreen)
  
}

##################################################################
## 
printconservedde_all <- function(conservedde_all){
  out <- NULL
  for(i in 1:5){
    x<-sgenescorer2_matrix[conservedde_all,i]
    names(x) <- rownames(sgenescorer2_matrix[conservedde_all,])
    if(i==1)
      x <- x[x<500]
    else
      x <- x[x<1000]
    print(colnames(sgenescorer2_matrix)[i])
    print(sort(x))
    out <- rbind(out, data.frame(screen=colnames(sgenescorer2_matrix)[i], gene=names(sort(x))))
  }
  out
#  vc <- vennCounts(allde)
#  vennDiagram(vc,cex=c(1.5,1.5,1.5))
}
#printconservedde_all(allde_conserved_alltime)
#
x <- getconservedDE(0.001)
printconservedde_all(x$sym_mouse[x$conserved_anytime])



###########################################################
############# read chip data ##############################
###########################################################


##################################################################
##
readallchiptot <- function(chipgenes = c("Gata3","Batf","Irf4","Stat6","Stat6m","Xbp1"),fname="gbi"){
  foo <- read.csv(sprintf("chip/%s_total.csv",fname),sep="\t",stringsAsFactors = FALSE)
  times <- c(2, 4,24,48,72)
  tp <- c("peak","Naive",sprintf("Th2_%sh",times))
  out <- as.data.frame(matrix(ncol=0, nrow=nrow(foo),0))
  out$peak <- foo[,1]
  
  chiprep<-c("a","b")
  
  for(g in chipgenes)
    for(time in times){
#      print(g=="Stat6m" & time==72)
      if(g=="Stat6m" & time==72)
        f<-sprintf("out.%s_%s_peaks.narrowPeak",g,chiprep)
      else {
        f<-sprintf("out.%s_%sh_%s_peaks.narrowPeak",g,time,chiprep)
      }
 #     print(f)
      f<-intersect(f,colnames(foo))
      if(length(f)>0){
        print(f)
        rf <- sprintf("%s_%sh",g,time)
        w <- apply(foo[,f,drop=FALSE]!="",1,mean)
        out[,rf] <- w
#        out[w,rf] <- w
      }
    }
  ann <- read.csv(sprintf("chip/%s_ann.csv",fname),sep="\t",stringsAsFactors = FALSE)
  colnames(ann)[1] <- "peak"
  out <- merge(out,ann)
}

## gata + batf + irf4 + raw stat6
dchiptot<- readallchiptot(fname="chip")
#colnames(dchiptot)
## gata + batf + irf4 + xbp1
dgbix <- readallchiptot(fname="chip3")
## gata + batf + irf4 + merged stat6 
dgbis <- readallchiptot(fname="chip2")
## gata + batf + irf4
dgbi <- readallchiptot(fname="gbi")

dgata3<- readallchiptot(fname="Gata3")
dbatf <- readallchiptot(fname="Batf")
dirf4 <- readallchiptot(fname="Irf4")


##################################################################
## get if there is a chip peak for certain time points
chip_tp <- function(tf,out,tp){
  out <- out[abs(out$Distance.to.TSS)<maxTSSdistCHIP,]  ##TSS distance cut-off
  li <- apply(out[,tp]==1,1,any)
  v <- table(out$Nearest.Ensembl[li])
  v <- data.frame(
    jasparname=rep(tf,length(v)),
    TSS_ensg=names(v),
    cnt=as.double(v)
  )
  v
}
## get if there is a peak at any time
chip_alltime <- function(){
  tp <- 2:7
  rbind(
    #chip_tp("chip_Stat6",dstat6,tp),
    chip_tp("chip_Gata3",dgata3,tp),
    chip_tp("chip_Batf",dbatf,tp),
    chip_tp("chip_Irf4",dirf4,tp),
    chip_tp("chip_Xbp1",dxbp1,tp)
  )
}
c_alltime  <- chip_alltime()




##################################################################
#### Read and prepare ATAC peaks #################################
##################################################################

##################################################################
##
getnormATAC <- function(org, atac){

  sumcolpairs <- function(x){
    y <- matrix(0,ncol=ncol(x)/2,nrow=nrow(x))
    for(i in 1:ncol(y)){
      y[,i] <- x[,i*2-1]+x[,i*2]
    }
    y
  }

  newatac_ann <- read.csv(sprintf("atac/%s/ATACall_peaks.red.ann.csv",org),sep="\t",stringsAsFactors = FALSE)
  colnames(newatac_ann)[1] <- "peakid" #or something
  newatac_ann$Gene.Name <- normalizesym(newatac_ann$Gene.Name)

  ### Read background counts and figure out average counts
  newatac_inv <- read.table(sprintf("atac/%s/ATACall_peaks.inv.bed",org),sep="\t",stringsAsFactors = FALSE)
  newatac_inv_sum <- sum(as.numeric(newatac_inv$V3-newatac_inv$V2))
  newatac_bg <- read.table(sprintf("atac/%s/counts.f.bg.csv",org),sep="\t",stringsAsFactors = FALSE)
  newatac_bg <- apply(sumcolpairs(newatac_bg[,-(1:3)]),2,sum)
  newatac_bg_avgreads <- newatac_bg/newatac_inv_sum
  
  ### Read peak counts
  newatac_peaks <- read.table(sprintf("atac/%s/counts.f.peaks.csv",org),sep="\t",stringsAsFactors = FALSE)
  newatac_peaks <- cbind(newatac_peaks[,4,drop=FALSE],sumcolpairs(newatac_peaks[,-(1:6)]))
  colnames(newatac_peaks)
  colnames(newatac_peaks) <- c("peakid","Naive","Th2_2h","Th2_4h","Th2_24h","Th2_48h","Th2_72h")  #consider other parts of the code
  head(newatac_peaks)
  newatac_peaks <- smerge(newatac_peaks,newatac_ann)
  
  #Normalize peaks by background and length
  newatac_peaks_norm <- newatac_peaks
  newatac_peakslen <- (newatac_peaks$End-newatac_peaks$Start)
  for(i in 1:6){
    #Arbitrary unit to make it easier to think. now most peaks in 0-2. with up to 10
    newatac_peaks_norm[,i+1] <- 1e5*(newatac_peaks[,i+1]-newatac_bg_avgreads[i])/newatac_bg[i]/newatac_peakslen
  }
  print("Normalized ATAC peak counts by background/length")

  #Normalize over time
  newatac_peaks_norm_time <- newatac_peaks_norm
  ts <- newatac_peaks_norm_time[,1+2]    #normalize by second time point
  #print(head(ts))
  for(i in 1:6){
    newatac_peaks_norm_time[,i+1] <- newatac_peaks_norm_time[,i+1]/ts
  }
  print("Normalized ATAC peak counts over time")
  
  
  ### Peak -> Scaled size over time and distance to TSS
  mapPeakGene <- newatac_peaks_norm_time[,c("peakid","Nearest.Ensembl","Gene.Name","Naive","Th2_2h","Th2_4h","Th2_24h","Th2_48h","Th2_72h", "Distance.to.TSS")]
  colnames(mapPeakGene) <- c("peakid","TSS_ensg","gene","Naive","Th2_2h","Th2_4h","Th2_24h","Th2_48h","Th2_72h","TSS_distance")
  #mapPeakGene_unfiltered <- mapPeakGene
  mapPeakGene <- mapPeakGene[abs(mapPeakGene$TSS_distance)<30e3,]

  mapSiteGene <- smerge(mapPeakGene,atac$mapPeakInfo[,c("jasparname","peakid")])
  print("site->gene mapping done")  
  
  tfattall <- sqldf("select distinct jasparname, sum(`Naive`) as cnt1, sum(`Th2_2h`) as cnt2, sum(`Th2_4h`) as cnt3, sum(`Th2_24h`) as cnt4, sum(`Th2_48h`) as cnt5, sum(`Th2_72h`) as cnt6 from `mapSiteGene` group by jasparname")
  rownames(tfattall) <- tfattall$jasparname
  tfattall <- tfattall[,-1]
  tfattall
}



normlevatacTime <- function(tfatall){
  tfatall_normtime2 <- tfatall
  temp <- tfatall_normtime2[,2]
  for(i in 1:6){
    tfatall_normtime2[,i] <- tfatall_normtime2[,i]/temp
  }
  tfatall_normtime2[order(tfatall_normtime2[,6]),]
}

levatac.mouse <- getnormATAC("mouse", atac.mouse)
levatac.mouse.norm <- normlevatacTime(levatac.mouse)


### Store for website
write.csv(levatac.mouse.norm,"out_teichlab/th2crispr_mouse_TFchrom_data.csv",row.names = TRUE, quote = FALSE)



##################################################################
##
readnormATAC <- function(org){
  ### Read peak annotation
  newatac_ann <- read.csv(sprintf("atac/%s/ATACall_peaks.red.ann.csv",org),sep="\t",stringsAsFactors = FALSE)
  colnames(newatac_ann)[1] <- "peakid" #or something
  newatac_ann$Gene.Name <- normalizesym(newatac_ann$Gene.Name)
  head(newatac_ann)

  ### Peak -> global position of peak
  mapPeakPos <- newatac_ann[,c("Chr","Start","End","peakid","Nearest.Ensembl","Gene.Name","Distance.to.TSS")]
  colnames(mapPeakPos) <- c("Chr","peakstart","peakend","peakid","TSS_ensg","gene","TSS_distance")
  
  ### Peak -> global position of peak
  # mapPeakPos <- newatac_peaks_norm[,c("Chr","Start","End","peakid")]
  # colnames(mapPeakPos) <- c("Chr","peakstart","peakend","peakid")

  
  ### Peak -> local info about peak and TFs in it
  mapPeakInfo <- read.csv(sprintf("atac/%s/fimo.txt",org),stringsAsFactors = FALSE,sep="\t")   #not convinced this is right
  #head(mapPeakInfo)
  colnames(mapPeakInfo)[1]<-"motifid"
  colnames(mapPeakInfo)[2]<-"jasparname"
  colnames(mapPeakInfo)[3]<-"peakid"     
  
  #mapPeakInfo <- mapPeakInfo[mapPeakInfo$p.value<1e-5,]  #did not seem to filter before!
  mapPeakInfo$jasparname <- normalizesym(mapPeakInfo$jasparname)
  print("Got local coordinates of motifs")

  ### Figure out absolute position of motifs
  mapSiteInfo <- smerge(mapPeakInfo[,c("motifid","jasparname","peakid","start","stop","strand")], newatac_ann[,c("peakid","Chr","Start","Nearest.Ensembl","Gene.Name","Distance.to.TSS")])
#  mapMotifPos <- smerge(mapMotifPos, mapPeakInfo)
  mapSiteInfo$motifstart <- mapSiteInfo$start + mapSiteInfo$Start-1
  mapSiteInfo$motifend   <- mapSiteInfo$stop  + mapSiteInfo$Start-1   #I suspect this -1 is correct
  mapSiteInfo <- mapSiteInfo[,c("peakid","motifid","jasparname","strand","Chr","motifstart","motifend","Nearest.Ensembl","Gene.Name","Distance.to.TSS")]
  colnames(mapSiteInfo)[colnames(mapSiteInfo)=="Nearest.Ensembl"] <- "TSS_ensg"
  colnames(mapSiteInfo)[colnames(mapSiteInfo)=="Gene.Name"]       <- "gene"
  colnames(mapSiteInfo)[colnames(mapSiteInfo)=="Distance.to.TSS"] <- "TSS_distance"
  print("Got global coordinates of motifs")

  ### Write BED file with absolute coordinates of the sites
  abed <- mapSiteInfo[,c("Chr","motifstart","motifend","jasparname")]
  abed$motifstart <- format(abed$motifstart , scientific = FALSE)
  abed$motifend   <- format(abed$motifend , scientific = FALSE)
  write.table(abed,sprintf("atac/%s/sites.bed",org),row.names = FALSE,col.names=FALSE, quote = FALSE)
  print("Wrote TF site bed file")

  #Cache and return result  
  list(
    mapSiteInfo=mapSiteInfo,
    mapPeakInfo=mapPeakInfo
    )
}

##################################################################
##
writeBEDforATACsites <- function(org, mapSiteInfo, outf=sprintf("atac/%s/sites.bed",org)){
  abed <- mapSiteInfo[,c("Chr","motifstart","motifend","jasparname")]
  abed$motifstart <- format(abed$motifstart , scientific = FALSE)
  abed$motifend   <- format(abed$motifend , scientific = FALSE)
  write.table(abed,outf,row.names = FALSE,col.names=FALSE, quote = FALSE)
  print("Wrote TF site bed file")
}




############################
## Calculate # sites over time. Need to be rewritten
#colnames(newatac_peaks_norm)
mapsitelevelATAC <- function(d){
  ### Site -> Scaled size over time and distance to TSS
  mapMotifGene <- smerge(d$mapPeakInfo,d$mapPeakGene)

  ### TF -> Summed activity over sites at different times
  tfattall <- sqldf("select distinct jasparname, sum(`Naive`) as cnt1, sum(`Th2_2h`) as cnt2, sum(`Th2_4h`) as cnt3, sum(`Th2_24h`) as cnt4, sum(`Th2_48h`) as cnt5, sum(`Th2_72h`) as cnt6 from mapMotifGene group by jasparname")
  rownames(tfattall)<-tfattall$jasparname
  tfattall<-tfattall[,-1]
  colnames(tfattall)<-c("0h","2h","4h","24h","48h","72h")
  
  d$tfattall <- tfattall
  d
}


##################################################################
#### Putative binding site conservation ##########################
##################################################################

##################################################################
##
getConservedSites <- function(mapSiteInfo, flifted){
  #Turn peak info into a grange
  grPeakInfo<-makeGRangesFromDataFrame(data.frame(
    chr=sprintf("%s_%s",mapSiteInfo$Chr,mapSiteInfo$jasparname), 
    start  =mapSiteInfo$motifstart, 
    end    =mapSiteInfo$motifend,
    strand =mapSiteInfo$strand,
    peakid =mapSiteInfo$peakid # 1:nrow(x)
  ), keep.extra.columns=TRUE)
  
  
  #Get the lifted sequence and turn into a grange
  lifted <- read.table(flifted,sep="\t")
  colnames(lifted) <- c("chr","start","end","jasparname")
  lifted$chr <- sprintf("%s_%s", lifted$chr, lifted$jasparname)
  grLifted<-makeGRangesFromDataFrame(lifted)
  
  #See which TF sites are preserved
  grPeakInfo_int <- findOverlaps(grPeakInfo, grLifted, ignore.strand=TRUE) 

  mapSiteInfoConserved <- mapSiteInfo[unique(from(grPeakInfo_int)),]
  print(nrow(mapSiteInfoConserved)/nrow(mapSiteInfo))
  mapSiteInfoConserved
}





#Write new human BED file for all peaks
# writehumanBedATAC <- function(){
#   newatac_ann <- read.csv(sprintf("atac/human/ATACall_peaks.red.ann.csv"),sep="\t",stringsAsFactors = FALSE)
#   f <- function(y) format(y , scientific = FALSE)
#   
#   bedhumanatac <- data.frame(
#     chr=sprintf("chr%s",newatac_ann$Chr),
#     start=f(newatac_ann$Start), 
#     end=f(newatac_ann$End),
#     strand=newatac_ann$Strand,
#     stringsAsFactors = FALSE
#   )
#   write.table(bedhumanatac,sprintf("atac/lift/human.bed"),row.names = FALSE,col.names=FALSE, quote = FALSE)
# }
# writehumanBedATAC()




##################################################################
#### ATAC peak conservation ######################################
##################################################################

##################################################################
## Check conservation on peak level
getConservedPeaks <- function(org, flifted, istm=FALSE){
#  org <- "mouse"
  if(istm){
    ownpeak <- read.csv(sprintf("atac/%s/tm/ATACall_peaks.red.ann.csv",org),sep="\t",stringsAsFactors = FALSE)
  } else{
    ownpeak <- read.csv(sprintf("atac/%s/ATACall_peaks.red.ann.csv",org),sep="\t",stringsAsFactors = FALSE)
  }
 # flifted="atac/lift/lifted.peaks.human.bed"
  
  #mapPeakInfo <- atac.human$mapPeakInfo
  #flifted="atac/lift/lifted.sites.mouse.bed"
  
  #Turn peak info into a grange
  grOwnPeak<-makeGRangesFromDataFrame(data.frame(
    chr=ownpeak$Chr,
    start=ownpeak$Start, 
    end=ownpeak$End,
    strand=ownpeak$Strand
  ))
  
  #Get the lifted sequence and turn into a grange
  lifted <- read.table(flifted,sep="\t")
  colnames(lifted) <- c("chr","start","end","xxx")
  grLifted<-makeGRangesFromDataFrame(lifted)
  
  #See which peaks are preserved
  grPeakInfo_int <- findOverlaps(grOwnPeak, grLifted, ignore.strand=TRUE) 
  
  data.frame(
    nownTot=nrow(ownpeak),
    nownOverlap=length(unique(from(grPeakInfo_int))),
    nLifted=nrow(lifted))
}

##################################################################
## Plot how many peaks overlap
makeATACPeakOverlapPlot <- function(){
  #Note: no big difference with 0.2 and 0.6 cutoff in sequence conservation
  pdf("atac/lift/overlap.pdf",height = 3)
  if(FALSE){
    statPeakOverlap <- rbind(
      getConservedPeaks("mouse", "atac/lift.tm/lifted.peaks.human.bed", TRUE),
      getConservedPeaks("human", "atac/lift.tm/lifted.peaks.mouse.bed", TRUE))
  } else {
    statPeakOverlap <- rbind(
      getConservedPeaks("mouse", "atac/lift/lifted.peaks.human.bed"),
      getConservedPeaks("human", "atac/lift/lifted.peaks.mouse.bed"))
  }
  statPeakOverlap$nLifted <- rev(statPeakOverlap$nLifted)
  statPeakOverlap$notinother <- statPeakOverlap$nownTot - statPeakOverlap$nLifted
  statPeakOverlap$liftedbutnotoverlap <- statPeakOverlap$nLifted - rev(statPeakOverlap$nownOverlap)
  barplot(
    t(as.matrix(statPeakOverlap[,c("nownTot","liftedbutnotoverlap","nownOverlap")])),
    col=c(rgb(230,159,0,maxColorValue = 255),rgb(86,180,233,maxColorValue = 255),rgb(0,158,115,maxColorValue = 255)),
    horiz=TRUE,
    names.arg=c("Mouse","Human")
  )
  dev.off()
  
  ####TODO hmm... where are the missing peaks? are they closer to genes or anything?
  
}


#################################################################
#### Merge peaks and detected motifs in them ##################### Only for mouse right now
##################################################################

############
##### combine ATAC peak input files and count peaks per gene
calcgenetfcount <- function(mapMotifGene){
  zz<-mapMotifGene[,c("jasparname","TSS_ensg")]
 sqldf("select distinct jasparname, TSS_ensg, count(TSS_ensg) as cnt from zz group by jasparname, TSS_ensg")
}

##################################################################
## Get TF site count per gene
getmarasitecountmatrix <- function(genetfcount){
  #returns an annoying jasparname-row. but don't change, breaks code
  d <- dcast(genetfcount, jasparname~TSS_ensg, fill=0, value.var = "cnt")
  #colnames(d)[1:3]
  rownames(d) <- d[,1]
  d<-t(d)[-1,]   #removes jasparname row - might break some code!!!
  class(d) <- "double"
  d
}


##################################################################
## .... Use cached result if possible
fname_atac_mouse <- sprintf("atac/%s.RData","mouse")
fname_atac_human <- sprintf("atac/%s.RData","human")
if(file.exists(fname_atac_mouse)){
  atac.mouse <- readRDS(fname_atac_mouse)
  atac.human <- readRDS(fname_atac_human)
} else {
  
  atac.mouse <- readnormATAC("mouse")
  atac.human <- readnormATAC("human")

  object.size(atac.mouse)/1e6
  
  #Non-conserved TF-gene matrix
  #atac.mouse$noncons_tfc <- calcgenetfcount(atac.mouse$mapSiteInfo)
  atac.mouse$noncons_tfc <- rbind(c_alltime,calcgenetfcount(atac.mouse$mapSiteInfo))
  atac.human$noncons_tfc <- calcgenetfcount(atac.human$mapSiteInfo)  #TODO: should ideally have chipseq data here too?

  #object.size(atac.mouse$noncons_tfc)/1e6
  
  #Note: lifting, 0.2 vs 0.6: seems 25% more peaks are lifted over. but this has little improvement on the site overlap
  atac.mouse$mapSiteInfo <- getConservedSites(atac.mouse$mapSiteInfo, "atac/lift/lifted.sites.human.bed")
  #20% of mouse peaks left
  atac.human$mapSiteInfo <- getConservedSites(atac.human$mapSiteInfo, "atac/lift/lifted.sites.mouse.bed")
  #12% of human peaks left
  
  # need to rethink
  # atac.mouse <- mapsitelevelATAC(atac.mouse)  
  # atac.human <- mapsitelevelATAC(atac.human)  
  
  atac.mouse$cons_tfc <- rbind(c_alltime,calcgenetfcount(atac.mouse$mapSiteInfo))
  atac.human$cons_tfc <- calcgenetfcount(atac.human$mapSiteInfo)  #TODO chipseq?

  writeBEDforATACsites("mouse",atac.mouse$mapSiteInfo, "conservedsite.bed")
  
  
  #Cache  
  saveRDS(atac.mouse, fname_atac_mouse)
  saveRDS(atac.human, fname_atac_human)
}


##################################################################
### Conserved ATAC peaks #########################################
##################################################################

#Do they go up and down the same way? are the sizes similar?






##################################################################
### Extract absolute ATAC motif coordinates ######################
##################################################################


writeMotifBed <- function(motif=NULL){
  
  mapMotifPosBed <- unique(data.frame(
    chr=mapMotifPos$Chr,
    start=as.integer(mapMotifPos$motifstart),
    end=as.integer(mapMotifPos$motifend),
    name=mapMotifPos$jasparname,
    score=rep(1000,nrow(mapMotifPos)),
    strand=mapMotifPos$strand))
  list_interesting_tf
  
  if(is.null(motif)){
    mapMotifPosBed <- mapMotifPosBed[mapMotifPosBed$name %in% c(expressed_atacTF_50,"Etv2"),]
    write.table(mapMotifPosBed,sprintf("out_motif/motifs.ALL.bed"),col.names = FALSE, row.names = FALSE,quote = FALSE,sep="\t")
  } else {
    red <- mapMotifPosBed[mapMotifPosBed$name %in% motif,]
    write.table(red,sprintf("out_motif/motifs.red.bed"),col.names = FALSE, row.names = FALSE,quote = FALSE,sep="\t")
  }
}
# writeMotifBed()
# writeMotifBed(c("Yy1","Tbx21","Pou6f1","Pou2f2","Etv2","Etv6","E2f4","Runx1","Foxo1","Ctcf","Ewsr1-fli1","Nrf1","Spib","Spi1","Ikzf3",
#                 "Stat6","Stat4","Epas","Nfil3"))
# writeMotifBed(list_interesting_tf) #from the atac-screen combination




##################################################################
####### Score ATAC motifs as early/late ##########################
##################################################################

kmeans.atacT <- function(atac, tc){
  ## Perform k-kmeans on normalized trends.
  ## Groups are boring - from early to late.
  set.seed(0)
  forkm <- atac$tfattall
  for(i in 1:nrow(forkm)){
    forkm[i,] <- forkm[i,]/mean(forkm[i,])
  }
  atackm <- kmeans(forkm,5)
  kmcol <- brewer.pal(max(atackm$cluster),"Set1")
  
  ## Show the k-means groups
  plot(apply(forkm[atackm$cluster==1 & rownames(tfattall) %in% tc$expressed_atacTF,],2,mean),type="l",ylim=c(0,2),col=kmcol[1])
  for(i in 2:max(atackm$cluster))
    lines(apply(forkm[atackm$cluster==i & rownames(tfattall) %in% tc$expressed_atacTF,],2,mean),col=kmcol[i])
}


## Base color on when the motif is present
## This does the job as well as k-means. score=0 early. score=1 late
calc_score_ael <- function(tfattall_red){
  #wt <- apply(tfattall_red,1,function(x) sum(x*(1:6))/sum(x))
  wt <- 1-tfattall_red[,2]/tfattall_red[,5]
  wt <- (wt-min(wt))/(max(wt)-min(wt))
  wt <- wt - median(wt)
  wt[wt<0] <- wt[wt<0]/-min(wt)
  wt[wt>0] <- wt[wt>0]/max(wt)
  wt
}
col_from_ael <- function(wt){
  thecol <- rep("black",length(wt))
  rsc <- function(x) abs(x)^0.5
  r<-rsc(wt)

  x<-r[wt<0]
  thecol[wt<0]  <- rgb(x,0,0) 
  x<-r[wt>=0]
  thecol[wt>=0] <- rgb(0,0,x)
  
  thecol
  ##http://www.somersault1824.com/tips-for-designing-scientific-figures-for-color-blind-readers/
}

atac.mouse$score_ael <- calc_score_ael(atac.mouse$tfattall)
atac.human$score_ael <- calc_score_ael(atac.human$tfattall)




##################################################################
####### Curated list of genes ####################################
##################################################################

### Read list of transcription factors
list_cm <- rownames(read.csv("tflist/Mus_musculus_chromatin_remodeling_factors_gene_list.txt",sep="\t",stringsAsFactors = FALSE))
list_co <- rownames(read.csv("tflist/Mus_musculus_transcription_co-factors_gene_list.txt",sep="\t",stringsAsFactors = FALSE))
list_tf <- read.csv("tflist/Mus_musculus_transcription_factors_gene_list.txt",sep="\t",stringsAsFactors = FALSE)[,1]

### Read protein atlas protein annotation
protatlas <- read.csv("tflist/proteinatlas.tab",sep="\t",stringsAsFactors = FALSE)
#unique(protatlas$Subcellular.location)

list_protatlas_secreted <- ensconvert$ensembl_gene_id[
  str_to_lower(ensconvert$mgi_symbol) %in%
    str_to_lower(protatlas$Gene[c(
      grep("secreted",protatlas$Protein.class),
      grep("membraneNOPE",protatlas$Protein.class)
    )])]
list_protatlas_membrane <- ensconvert$ensembl_gene_id[
  str_to_lower(ensconvert$mgi_symbol) %in%
    str_to_lower(protatlas$Gene[c(
      grep("membrane",protatlas$Protein.class)
    )])]


##################################################################
####### ThExpress data ###########################################
##################################################################

d <-         read.csv("thexpress/Th2_vs_naive.txt",sep="\t",stringsAsFactors = FALSE)[,c(2,6),drop=FALSE]
d <- cbind(d,read.csv("thexpress/Th2_vs_Th1.txt",  sep="\t",stringsAsFactors = FALSE)[,c(2,6),drop=FALSE])
d <- cbind(d,read.csv("thexpress/Th2_vs_Th17.txt", sep="\t",stringsAsFactors = FALSE)[,c(2,6),drop=FALSE])
d <- cbind(d,read.csv("thexpress/Th2_vs_iTreg.txt",sep="\t",stringsAsFactors = FALSE)[,c(2,6),drop=FALSE])
d <- cbind(d,read.csv("thexpress/Th2_vs_nTreg.txt",sep="\t",stringsAsFactors = FALSE)[,c(2,6),drop=FALSE])
thep <- d[,c(2,4,6,8,10)]
thefc <- d[,c(1,3,5,7,9)]
colnames(thep)<-c("Naive/Th0","Th1","Th17","iTreg","nTreg")
colnames(thefc)<-c("Naive/Th0","Th1","Th17","iTreg","nTreg")








##################################################################
####### DMDD project data ########################################
##################################################################

dmdd <- read.csv("dmdd/dmdd_embryo_annotations_20170306.tsv",sep="\t",stringsAsFactors = FALSE)
dmdd <- dmdd[dmdd$MP.ID %in% c(
  "MP:0000690","MP:0000692","MP:0000694","MP:0000703","MP:0000705","MP:0000706",
  "MP:0001879","MP:0002368","MP:0010200","MP:0013970"),]
#"MP:0001914","MP:0001916","MP:0002633","MP:0002725","MP:0013970"
s <- sgenescorer2_matrix[rownames(sgenescorer2_matrix) %in% dmdd$Gene,"Xbp1",drop=FALSE]
s <- s[order(s[,1],decreasing = TRUE),,drop=FALSE]
s



