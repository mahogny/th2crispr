###############################################################################################
########### MARA reading functions ############################################################
###############################################################################################

###############################################
## 
gettfdestat <- function(){
  v1 <- rbind(
    smerge(ensconvert,tcmouse$de_early[,c("ensembl_gene_id","pval")]),  #combining like this a bit nasty
    smerge(ensconvert,tcmouse$de_late [,c("ensembl_gene_id","pval")]))
  v1$pval[is.na(v1$pval)] <- 1
  v1 <- smerge(v1, map_jaspar_namegenesym)
  tfde <- sqldf("select jasparname, min(pval) as pval from v1 group by jasparname")
  rownames(tfde) <- tfde$jasparname
  tfde  
}

###############################################
## Read TC MARA
readmara.tc <- function(dir){
  at    <- read.table(pipe(sprintf("tar -xOf %s/ismara_report.tar.gz ismara_report/activity_table",dir)), stringsAsFactors = FALSE)
  am    <- read.table(pipe(sprintf("tar -xOf %s/ismara_report.tar.gz ismara_report/active_matrices",dir)),stringsAsFactors = FALSE)
  delta <- read.table(pipe(sprintf("tar -xOf %s/ismara_report.tar.gz ismara_report/delta_table",dir)),    stringsAsFactors = FALSE)
  
  #Reorder scores in the same order as the other tables
  colnames(am)<-c("jasparname","maraz")
  am <- am[order(am$jasparname),]
  
  #Perform averaging
  doav <- function(thi,at){
    thtime <- c("Naive",sprintf("Th%s_%s",thi,c("05h","2h","4h","6h","12h","24h","48h","72h")))
    for(i in 1:length(thtime)){
      a <- apply(at[  #note. different naming in human
        rownames(at) %in% c(sprintf("%s_rep%s_salmon",thtime[i],1:3),  sprintf("%s_%s",thtime[i],1:3)  )
        ,],2,mean)
      if(i==1)
        avact <- a
      else
        avact <- cbind(avact, a)
    }
    colnames(avact) <- thtime
    avact <- avact[,which(!is.na(avact[1,]))]
    avact
  }

  list(table0=doav(0,at),
       table2=doav(2,at),
       score=am,
       scoregene=smerge(map_jaspar_namegenesym, am),  ######### todo: chip: make an exception here
       delta0=doav(0,delta),
       delta2=doav(2,delta)
       )
}



###############################################
## Read KO MARA
readmara.ko <- function(dir){
#  dir <- "out_mara/newmara_ko"
  
  at    <- read.table(pipe(sprintf("tar -xOf %s/ismara_report.tar.gz ismara_report/activity_table",dir)), stringsAsFactors = FALSE)
  am    <- read.table(pipe(sprintf("tar -xOf %s/ismara_report.tar.gz ismara_report/active_matrices",dir)),stringsAsFactors = FALSE)
  delta <- read.table(pipe(sprintf("tar -xOf %s/ismara_report.tar.gz ismara_report/delta_table",dir)),    stringsAsFactors = FALSE)
  
  #Reorder scores in the same order as the other tables
  colnames(am)<-c("jasparname","maraz")
  am <- am[order(am$jasparname),]
  
  #Perform averaging
  ko_name <- str_split_fixed(rownames(at),"_",2)[,2] 
  uko <- unique(ko_name)
  for(i in 1:length(uko)){
    a <- apply(at[grep(sprintf("_%s",uko[i]),rownames(at)),],2,mean)
    asd <- apply(at[grep(sprintf("_%s",uko[i]),rownames(at)),],2,sd)
    if(i==1){
      avact <- a
      sdact <- asd
    }
    else{
      avact <- cbind(avact, a)
      sdact <- cbind(sdact, asd)
    }
  }

  #Extract DE info
  tfde <- gettfdestat()
  pde <- tfde[am$jasparname,]$pval
  pde[is.na(pde)] <- 1
  pde <- log10(pde)
  
  ## Extract DE - needed? need an update
  # de_late <- tcmouse$de_late
  # dep <- sqldf("select min(pval) as de_pval, ensembl_gene_id from de_late group by ensembl_gene_id")
  # 
  # sg <- smerge(map_jaspar_namegenesym, am)
  # sg <- smerge(sg, dep,all.x=TRUE)
  
  colnames(avact) <- uko
  list(mean=avact,
       sd=sdact,
       score=am,
       scoregene=pde#sg  ######### todo: chip: make an exception here
  )
}


###############################################################################################
########### Reading the data ##################################################################
###############################################################################################


mara_human_alltf_allgenes_alltime_th20 <- readmara.tc("out_mara/human_alltf_allgenes_alltime_th20")
mara_mouse_alltf_allgenes_alltime_th20 <- readmara.tc("out_mara/mouse_alltf_allgenes_alltime_th20")

mara_human_alltf_cDEgenes_alltime_th20 <- readmara.tc("out_mara/human_alltf_cDEgenes_alltime_th20")
mara_mouse_alltf_cDEgenes_alltime_th20 <- readmara.tc("out_mara/mouse_alltf_cDEgenes_alltime_th20")

mara_human_alltf_allgenes_alltime_th20_nonconsP <- readmara.tc("out_mara/human_alltf_allgenes_alltime_th20_nonconsP")
mara_mouse_alltf_allgenes_alltime_th20_nonconsP <- readmara.tc("out_mara/mouse_alltf_allgenes_alltime_th20_nonconsP")

mara_human_alltf_cDEgenes_alltime_th20_nonconsP <- readmara.tc("out_mara/human_alltf_cDEgenes_alltime_th20_nonconsP")
mara_mouse_alltf_cDEgenes_alltime_th20_nonconsP <- readmara.tc("out_mara/mouse_alltf_cDEgenes_alltime_th20_nonconsP")

####

mara_ko_alltf_allgenes          <- readmara.ko("out_mara/ko_alltf_allgenes")
mara_ko_alltf_allgenes_nonconsP <- readmara.ko("out_mara/ko_alltf_allgenes_nonconsP")

mara_ko_alltf_cDEgenes          <- readmara.ko("out_mara/ko_alltf_cDEgenes")
mara_ko_alltf_cDEgenes_nonconsP <- readmara.ko("out_mara/ko_alltf_cDEgenes_nonconsP")


###############################################################################################
########### Plotting functions ################################################################
###############################################################################################


###############################################
## Coloring for mara overview plots
getmarapcol <- function(stat){
  pcol <- rep("black",nrow(stat))
  pcol[grep("Stat6",rownames(stat))]<-"red"
  cchip <- "blue"
  pcol[grep("chip_",rownames(stat))]<-"cchip"
  pcol[rownames(stat) %in% c("Xbp1","Gata3","Irf4","Batf")] <-cchip
  pcol  
}


###############################################
## Get global statistics about each motif from MARA file
getmaraoverallstat <- function(m, normalize=TRUE){
  
  maraoverridewithchip <- function(avact2){
    torem <- c("Xbp1","Gata3","Irf4","Batf")
    avact2 <- avact2[rownames(avact2) %!in% torem,]
    for(i in torem){
      rownames(avact2)[rownames(avact2)==sprintf("chip_%s",i)] <- i
    }
    avact2
  }
  
  avact2 <- m$table2
  avact0 <- m$table0
  avact2 <- maraoverridewithchip(avact2)
  avact0 <- maraoverridewithchip(avact0)
  #Here the score table stops matching
  newscore <- smerge(
    data.frame(jasparname=rownames(avact2), stringsAsFactors = FALSE),
    m$score, all.x = TRUE)

  #pscale <- apply(cbind(avact2,avact0,0),1,max)-apply(cbind(avact2,avact0,0),1,min)
  #if(!normalize)
    pscale <- 1
  pdiff <- (avact2[,ncol(avact2)]-avact2[,1])/pscale
  pinc <- apply((avact2-avact0),1,sum)/pscale
  pcol <- rep("black",nrow(avact2))
  pcol[grep("Stat6",rownames(avact2))]<-"red"
  pcol[grep("chip_",rownames(avact2))]<-"red"
  
  #  head(tcmouse$de_late[,c("ensembl_gene_id","qval")])
  tfde <- gettfdestat()
  
  pde <- tfde[names(pinc),]$pval
  pde[is.na(pde)] <- 1
  pde <- log10(pde)
  
  rownames(m$score) <- m$score$jasparname
  
  # print(length(pinc))
  # print(length(m$score$maraz))
  mstat <- data.frame(
    jaspar_name=names(pinc),
    diffth0=pinc,
    diffnaive=pdiff,
    pde=pde,
    pscale=pscale,
    maraz=newscore$maraz#[names(pinc),]$maraz
  )
  if(normalize){
    temp <- lm(diffth0 ~ diffnaive -1,data = mstat)
    mstat$diffth0 <- mstat$diffth0 - mstat$diffnaive*temp$coefficients
  }
  mstat
}
#mstat <- getmaraoverallstat(m, FALSE)





###############################################
### Th2 vs Th0 difference, overall 2d scatter plot. Naive/Th2 vs Th0/Th2
plotMARAoverallDiffvsDiff <- function(stat){
  plot(stat$diffth0,stat$diffnaive,cex=0,
       xlab="Th2/0 average activity difference",
       ylab="Th2 Naive/72h activity difference")
  lines(minmax(stat$diffth0),c(0,0),lty=3,col="gray")
  lines(c(0,0),minmax(stat$diffnaive),lty=3,col="gray")
  plotdottext(
    stat$diffth0,
    stat$diffnaive,
    labels = stat$jaspar_name,cex=.7,col=getmarapcol(stat))
}


###############################################
## Th2 vs Th0 difference, and Z-score
plotMARAoverallDiffvsZ <- function(stat){
  stat <- stat[!is.na(stat$maraz),]
  plot(stat$diffth0,stat$maraz,cex=0,
       xlab="Th2/0 average activity difference",
       ylab="MARA relevance score")
  lines(minmax(stat$diffth0),c(0,0),lty=3,col="gray")
  lines(c(0,0),minmax(stat$maraz),lty=3,col="gray")
  plotdottext(
    stat$diffth0,
    stat$maraz,
    labels = stat$jaspar_name,cex=.7,col=getmarapcol(stat))
}


###############################################
##
plotdottext <- function(x,y,labels,cex,col=rep("black",length(x)),donew=FALSE){
  xnorm <- x/(max(x)-min(x))  #ratio
  ynorm <- y/(max(y)-min(y))
  pd <- data.frame(x=xnorm,y=ynorm)
  pd <- as.matrix(dist((pd)))

  #666
  #Local density
  nn <- apply(pd<0.05,1,sum)
  asp <- nn>3
  #Make exceptions
  cd <- apply(pd,1,function(v) sort(v)[2])
  asp[cd>0.1] <- FALSE
  
  asp[labels %in% c("Yy1","Yy2")] <- FALSE
  asp[grep("Fli1",labels)] <- FALSE
  
  if(donew)
    plot(x,y,cex=0)
  points(x[ asp], y[ asp],pch = 20,cex=0.5,col="gray")
  text(  x[!asp], y[!asp],labels = labels[!asp],cex=cex,col=col)
#  text(  x, y,labels = nn,cex=cex,col=col)
}
# plotMARAoverallDiffvsDiff(mstat)
# plotMARAoverallDiffvsZ(mstat)
# plotMARAoverallDiffvsDE(mstat)

###############################################
### Th2 vs Th0 difference, overall 2d scatter plot. Naive/Th2 vs Th0/Th2
plotMARAoverallDiffvsDE <- function(stat){
  plot(stat$pde,stat$diffth0,cex=0,  ##need better names
       xlab="Th2/0 TF DE-score (Log10 p-value)",
       ylab="Th2/0 activity difference")
  lines(minmax(stat$pde),c(0,0),lty=3,col="gray")
  lines(c(0,0),minmax(stat$diffth0),lty=3,col="gray")
  plotdottext(
    stat$pde,
    stat$diffth0,
    labels = stat$jaspar_name,cex=.7,col=getmarapcol(stat))
}



###############################################################################################
########### Global plotting ###################################################################
###############################################################################################


m <- mara_human_alltf_allgenes_alltime_th20_nonconsP
m <- mara_mouse_alltf_allgenes_alltime_th20_nonconsP

m <- mara_human_alltf_allgenes_alltime_th20
m <- mara_mouse_alltf_allgenes_alltime_th20  #chips really stand out in this conserved one!

m <- mara_human_alltf_cDEgenes_alltime_th20
m <- mara_mouse_alltf_cDEgenes_alltime_th20  

############################
## Global statistics #######

#777
mstat <- getmaraoverallstat(m, FALSE)
plotMARAoverallDiffvsDiff(mstat)
plotMARAoverallDiffvsZ(mstat)
plotMARAoverallDiffvsDE(mstat)



doallMARAglobalplots <- function(){
  v <- list(
    human_alltf_allgenes_alltime_th20_nonconsP = mara_human_alltf_allgenes_alltime_th20_nonconsP,
    mouse_alltf_allgenes_alltime_th20_nonconsP = mara_mouse_alltf_allgenes_alltime_th20_nonconsP,
    human_alltf_allgenes_alltime_th20          = mara_human_alltf_allgenes_alltime_th20,
    mouse_alltf_allgenes_alltime_th20          = mara_mouse_alltf_allgenes_alltime_th20,  #chips really stand out in this conserved one!
    human_alltf_cDEgenes_alltime_th20          = mara_human_alltf_cDEgenes_alltime_th20,
    mouse_alltf_cDEgenes_alltime_th20          = mara_mouse_alltf_cDEgenes_alltime_th20  
  )
  for(i in 1:length(v)){
    print(i)
    mstat <- getmaraoverallstat(v[[i]], FALSE)
    
    pdf(sprintf("out_mara/%s/plot_diffvsdiff.pdf",names(v)[i]))
    plotMARAoverallDiffvsDiff(mstat)
    dev.off()
    
    pdf(sprintf("out_mara/%s/plot_diffvsZ.pdf",names(v)[i]))
    plotMARAoverallDiffvsZ(mstat)
    dev.off()
    
    pdf(sprintf("out_mara/%s/plot_diffvsDE.pdf",names(v)[i]))
    plotMARAoverallDiffvsDE(mstat)
    dev.off()
    
    mstat <- getmaraoverallstat(v[[i]], TRUE)
    
    pdf(sprintf("out_mara/%s/plot_diffvsdiff_norm.pdf",names(v)[i]))
    plotMARAoverallDiffvsDiff(mstat)
    dev.off()
    
    pdf(sprintf("out_mara/%s/plot_diffvsZ_norm.pdf",names(v)[i]))
    plotMARAoverallDiffvsZ(mstat)
    dev.off()
    
    pdf(sprintf("out_mara/%s/plot_diffvsDE_norm.pdf",names(v)[i]))
    plotMARAoverallDiffvsDE(mstat)
    dev.off()
    
    
  }  
  
}
doallMARAglobalplots()

#Tfap4: Gata3 
#Etv2: Irf4 hit. 
#Foxo1 is an Il13/Xbp1 hit. 
#Fos is Irf4/xbp1 hit (weak)
#E2f1: irf4
#Yy1: Il4, Il13, Xbp1, Gata3
#Nfkb2: Il4
#Yy2: Gata3   high Z
#E2f7^Il13

###############################################################################################
########### Individual plotting ###############################################################
###############################################################################################


###############################################
## Plot activity in Th2 and Th0 for a gene
plotMARAact2vs0 <- function(m, gene, tpmgene=gene, main=str_replace_all(deparse(substitute(gene)),"\"","")){
  i<-which(rownames(m$table2)==gene)
  tp <- c("0h",".5h","2h","4h","6h", "12h","24h","48h","72h")
  lwd <- 4
  if(is.matrix(m$table0)){
    plot(m$table2[i,],type="l",
         ylim=symrange(c(m$table0[i,],m$table2[i,],c(m$table0[i,],m$table2[i,]))), 
         col="black",ylab="",xlab="",
         main=main,lwd=lwd,xaxt = "n")
    lines(m$table0[i,],
          col=rgb(0, 158, 115, maxColorValue = 255),
          lwd=lwd)
    lines(c(0,9), c(0,0),
          col="gray",
          lwd=1)
  } else {
    plot(m$table2[i,],type="l",ylim=c(min(c(m$table2[i,])),max(c(m$table2[i,]))), 
         main=main,lwd=lwd)
  }
  axis(1, at = 1:length(tp), labels=tp)
  par(new=T)
  plot(tcmouse$av_mtpm[toensid(tpmgene),], lwd=2, lty=2,
       ylim=minmax(c(0, tcmouse$av_mtpm[toensid(tpmgene),])),
       col=rgb(213, 94, 0, maxColorValue = 255), axes=F, xlab=NA, ylab=NA, type="l")
  axis(side = 4)
  
  
  #atactp <- c("0h","2h","4h","24h","48h","72h")
  par(new=T)
  plot(c(1,2,3,7,8,9),
       levatac.mouse.norm[tpmgene,],ylim=c(0,3),type="l",axes=F, xlab=NA,ylab=NA,
        lwd=4, lty=3, col="blue")
}



#Spib is downregulated but not 0 in Th2. called Sfpi1 in Th-express. 
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1171190/  interaction with Irf4




################################################################# 
## Big PDF with all the panels!
if(TRUE){
  m <-  mara_mouse_alltf_allgenes_alltime_th20
  pdf("out_mara/multi.pdf",w=4*3.5, h=5*3.5)
  par(mfrow=c(5,4))
  plotMARAact2vs0(m,"Gata3_1","Gata3", main="Gata3")
  plotMARAact2vs0(m,"chip_Gata3","Gata3",main = "Gata3 (ChIP)")   #Show both to contrast?
  
  plotMARAact2vs0(m,"Xbp1")
  plotMARAact2vs0(m,"chip_Xbp1","Xbp1",  main = "Xbp1 (ChIP)")  #conserved, similar behavior
  
  plotMARAact2vs0(m,"Batf3")
  plotMARAact2vs0(m,"chip_Batf", "Batf", main = "Batf (ChIP)")
  
  plotMARAact2vs0(m,"Irf4_1", "Irf4",    main = "Irf4")
  plotMARAact2vs0(m,"chip_Irf4", "Irf4", main = "Irf4 (ChIP)")
  
  dev.off()
}
if(TRUE){
  m <-  mara_mouse_alltf_allgenes_alltime_th20
  #Tfap4: Gata3 
  #Etv2: Irf4 hit. 
  #Foxo1 is an Il13/Xbp1 hit. 
  #Fos is Irf4/xbp1 hit (weak)
  #E2f1: irf4
  #Yy1: Il4, Il13, Xbp1, Gata3
  #Nfkb2: Il4
  #Yy2: Gata3   high Z
    
  pdf("out_mara/multi2.pdf",w=4*3.5, h=5*3.5)
  par(mfrow=c(5,4))
  plotMARAact2vs0(m,"Tfap4")
  plotMARAact2vs0(m,"Etv2")
  plotMARAact2vs0(m,"Foxo1")
  plotMARAact2vs0(m,"Fos")
  plotMARAact2vs0(m,"E2f1")
  plotMARAact2vs0(m,"Yy1")
  plotMARAact2vs0(m,"E2f7")
  plotMARAact2vs0(m,"Nfkb2")
  
  plotMARAact2vs0(m,"Stat1")
  plotMARAact2vs0(m,"Stat3")
  plotMARAact2vs0(m,"Stat4")
  plotMARAact2vs0(m,"Stat5a..stat5b","Stat5a")
  plotMARAact2vs0(m,"Stat6")
  
  plotMARAact2vs0(m,"Fli1")
  
  dev.off()

  ######## Effect of Stats ###############
  
  #sgenescorer2_matrix[c("Stat1","Stat2","Stat3","Stat4","Stat5a","Stat5b","Stat6"),]
  #Stat6-> gata3  stat2 -> il13. that is all
  #ThE: only Stat1 different at end state
  #Th2/0 DE: Stat1 "Stat2"         "Stat3"         "Stat4"         "Stat5a"        "Stat5b"
  
  #https://www.ncbi.nlm.nih.gov/pubmed/21215659
  # important to cite: Stat3 needed for Th2 (and Th17). multiple Stats used!
}


# How is the overlap peaks vs DE genes? kind of did this before?
#tsne compares different TFs
#fisher test!



################################################################# 
## Do all of them for the website
if(TRUE){

  m <-  mara_mouse_alltf_allgenes_alltime_th20
  # temp_jaspar <- map_jaspar_namegenesym
  # rownames(temp_jaspar) <- temp_jaspar$jasparname
  listout <- list()
  ni <- 1
  for(i in 1:nrow(m$table0)){ 
    #why is .. batf going down?
    maraname <- rownames(m$table0)[i]
    jn <- maraname
    if(jn=="chip_Gata3")
      jn<-"Gata3_1"
    if(jn=="chip_Batf")
      jn<-"Batf3"
    if(jn=="chip_Irf4")
      jn<-"Irf4_1"
    if(jn=="chip_Xbp1")
      jn<-"Xbp1"
 #   print(jn)
    jn <- str_replace_all(jn,"\\.\\.","::")
#    print(jn)
    sym <- map_jaspar_namegenesym[map_jaspar_namegenesym$jasparname==jn,2]
    if(length(sym)>0){
      sym <- sym[1]
      
      if(toensid(sym) %in% rownames(tcmouse$av_mtpm)){
        listout[[i]]  <- data.frame(i=ni, sym=sym, maraname=maraname)
        svg(sprintf("out_mara/mouse_alltf_allgenes_alltime_th20/forsite/%s.svg",ni),width = 5, height = 5)
        plotMARAact2vs0(m,maraname, sym, main = sym)
        dev.off()
        ni <- ni + 1
      }
      
    } else {
      print(jn)
    }
  }

  listout <- rbindlist(listout)
  listout$fname <- sprintf("%s.svg",listout$i)    
  
  pastewspace <- function(sdata) do.call(paste, c(as.list(sdata), sep=" "))
  outhtml <- pastewspace(sprintf("<div style=\"position: float\"><h1>%s</h1><img src=\"%s\"/></div>",listout$sym, listout$fname  ))

  sprintf("<th>%s</th>",colnames(sgenescorer2_matrix))
  sprintf("<td>%s</td>",sgenescorer2_matrix[1,])
  #add screen score?
  
  write(outhtml,"out_mara/mouse_alltf_allgenes_alltime_th20/forsite/index.html")
}

