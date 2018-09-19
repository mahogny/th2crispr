###############################################################################################
###### Write ARACNE input files ###############################################################
###############################################################################################


################################################################
### Write an ARACNE input file
write_aracne_ens <- function(x, filestart, expressedGenes_=expressedGenes10Id, ensconvert_=ensconvert){
  x <- x[rownames(x) %in% expressedGenes_,] 
  x <- cbind(data.frame(gene=rownames(x)),x)
  write.table(x,sprintf("aracne/%s.exp",filestart),quote = FALSE,row.names = FALSE,sep="\t")
  
  #To be really sure we get all TFs. particularly for human!
  backup_tf <- ensconvert_$ensembl_gene_id[ensconvert_$mgi_symbol %in% togenesym(list_tf)]
  
  list_tf_red <- c(backup_tf,list_tf) # togenesym2(list_tf) #, ensconvert_ = ensconvert_  this is mouse 
  list_tf_red <- data.frame(tf=intersect(list_tf_red,rownames(x)),stringsAsFactors = FALSE)
  write.table(
    list_tf_red,sprintf("aracne/%s.tf",filestart),quote = FALSE, row.names = FALSE, col.names = FALSE)
  ## Full list of genes
  write.table(
    rownames(x),sprintf("aracne/%s.genes",filestart),quote = FALSE, row.names = FALSE, col.names = FALSE)
}


################################################################
## Calculate Th2 points with Th0 averages removed
calcth2sub0 <- function(tc){
  mtpm_th2sub0 <- tc$mtpm_th2
  times0 <- c("Naive",sprintf("Th0_%sh",c("05","1","2","4","6","12","24","48","72")))
  times2 <- c("Naive",sprintf("Th2_%sh",c("05","1","2","4","6","12","24","48","72")))
  for(i in 1:length(times0)){
    om <- apply(tc$mtpm[,startsWith(colnames(tc$mtpm),times0[i])],1,mean)
    for(j in which(startsWith(colnames(mtpm_th2sub0),times2[i]))){
      mtpm_th2sub0[,j] <- mtpm_th2sub0[,j] - om
    }
  }
  mtpm_th2sub0
}


write_aracne_humantc <- function(){
  write_aracne_ens(tchuman$mtpm,"tc.human",expressedGenes_ = human_expressedGenesId, ensconvert_ = tchuman$ensconvert)
}
write_aracne_mousetc <- function(){
  write_aracne_ens(tcmouse$mtpm,"tc.mouse")
}
write_aracne_tcga <- function(){
  fpkm_thymoma <- read.csv("tcga.thymoma/combine/tcga.csv", row.names = 1)
  e <- rownames(fpkm_thymoma)[apply(fpkm_thymoma,1,mean)>1] 
  print(length(e))
  write_aracne_ens(fpkm_thymoma,"tcga_thymoma",expressedGenes_ = e, ensconvert_ = human_ensconvert)
}
write_aracne_ko <- function(){
  colnames(ncounts_corr) <- as.character(phenoData(set2)$ko)
  write_aracne_ens(ncounts_corr,"ko")
}
write_aracne_ko1 <- function(m){
  x <- ncount
  colnames(x) <- cellcondition$ko
  x <- x[, cellcondition$isgood & cellcondition$mouse==m]
  write_aracne_ens(x,sprintf("ko%s",m))
}

for(i in 1:3)
  write_aracne_ko1(i)
write_aracne_ko()
write_aracne_mousetc()
write_aracne_humantc()   #very short list of TFs!
write_aracne_tcga()
write_aracne_ens(calcth2sub0(tcmouse),"difftc.mouse")
write_aracne_ens(calcth2sub0(tchuman),"difftc.human",expressedGenes_ = human_expressedGenesId, ensconvert_ = tchuman$ensconvert)



###############################################################################################
#### Functions to analyze the network #########################################################
###############################################################################################


###Filter out and keep non-TF interactions from ARACne. Keep MARA direction
#sum(!(aracne_ko$Regulator %in% utf) & !(aracne_ko$Target %in% utf))
#Turns out there are no such!!!

888
#####################
##### Function: Merge ARACNe and ATAC data. They output will be a directed network
filterAracneATACCHIP <- function(net, ensconvert_=ensconvert, tfc=atac.mouse$cons_tfc){

  #Should ideally do orthology mapping for the human TFs from map_jaspar

  tfc <- atac.mouse$cons_tfc[,c(1,2)]
  tfc$jasparname <- as.character(tfc$jasparname)
  tfc$TSS_ensg <- as.character(tfc$TSS_ensg)
  torem <- c("Xbp1","Gata3","Irf4","Batf")
  tfc <- tfc[tfc$jasparname %!in% torem,]
  for(i in torem){
    tfc$jasparname[tfc$jasparname==sprintf("chip_%s",i)] <- i
  }
  mapSiteInfo_ <- tfc
  
  print(nrow(mapSiteInfo_))
  
  
#  mapSiteInfo_ <- mapSiteInfo_[abs(mapSiteInfo_$TSS_distance) < 30e3,c("jasparname","TSS_ensg")]
  ensconvert_ <- ensconvert_[ensconvert_$mgi_symbol %in% map_jaspar_namegenesym$mgi_symbol,]
  tfl <- smerge(map_jaspar_namegenesym, ensconvert_)
  mapSiteInfo_ <- smerge(tfl, mapSiteInfo_)[,c("ensembl_gene_id","TSS_ensg")]  #TF is ensembl_gene_id, TSS_ensg is target

  print(nrow(mapSiteInfo_))
  print(head(mapSiteInfo_))
  
  ### Edges to retain as they are
  sub <- net[sprintf("%s_%s",net$Regulator,net$Target) %in% sprintf("%s_%s",mapSiteInfo_$ensembl_gene_id,mapSiteInfo_$TSS_ensg),]
  
  ### Edges to flip
  sub_rev <- net[sprintf("%s_%s",net$Target,net$Regulator) %in% sprintf("%s_%s",mapSiteInfo_$ensembl_gene_id,mapSiteInfo_$TSS_ensg),]
  sub_rev <- data.frame(Regulator=sub_rev$Target, Target=sub_rev$Regulator, MI=sub_rev$MI, pvalue=sub_rev$pvalue)
  
  rbind(sub,sub_rev)
}

#aracne_ko13_int_red <- filterAracneATACCHIP(aracne_ko13_int, atac.mouse$mapPeakGene, ensconvert)
#aracne_ko13_int_red


###################
# Function: Filter ARACNe by pvalue and MI cutoffs
filterAracneCutoff <- function(x, pvalue=1, MI=0){
  print(nrow(x))
  x <- x[aracne_ko_atac$pvalue < pvalue & aracne_ko_atac$MI>MI,] 
  print(nrow(x))
  x
}

##################
# Function: Intersect two networks - will not work between species!
intersectARACNe <- function(neta, netb){
  anyo <- c(sprintf("%s-%s",netb$Regulator, netb$Target),sprintf("%s-%s",netb$Target, netb$Regulator))
  netc <- neta[sprintf("%s-%s",neta$Regulator, neta$Target) %in% anyo,]
  print(dim(neta))
  print(dim(netb))
  print(dim(netc))
  netc
}


undirectARACNe <- function(net){
  rbind(
    net,
    data.frame(Target=net$Regulator, Regulator=net$Target, MI=net$MI, pvalue=net$pvalue)
  )
}
mergeARACne <- function(neta, netb){
  net <- rbind(neta,netb)
  #cheating a bit - MI and pval need not be from the same link, but most likely are
  sqldf("select distinct `Regulator`, `Target`, max(`MI`) as `MI`, min(pvalue) as pvalue from net group by `Regulator`, `Target`")
}

###################
# Function: Only keep DE genes in net
filterAracneDE <- function(net=aracne_mergetc_atac, alsokeep=c(), isMouse=TRUE, reqboth=FALSE){
  net$Regulator <- as.character(net$Regulator)
  net$Target <- as.character(net$Target)
  if(isMouse){
    the_de <- c(alsokeep, allde$ens_mouse[allde$conserved_loose])
  } else {
    the_de <- c(alsokeep, allde$ens_human[allde$conserved_loose])
  }

  if(reqboth)
    net <- net[net$Regulator %in% the_de & net$Target %in% the_de,] 
  else
    net <- net[net$Regulator %in% the_de | net$Target %in% the_de,] 
  #  print(dim(net))
  net
}



countaracneneigh <- function(net){
  table(c(net$Regulator, net$Target))
}

########################################
### Only keep TFs with N+ partners
filterAracneNeighN <- function(net,n=2, rep=1){
  net$Regulator <- as.character(net$Regulator)
  net$Target <- as.character(net$Target)
  keeptf <- names(which(table(c(net$Regulator,net$Target))>=n))
  net <- net[net$Regulator %in% keeptf & net$Target %in% keeptf,]
  print(sprintf("after neigh %s",nrow(net)))
  if(rep==1){
    net
  } else {
    filterAracneNeighN(net,n,rep-1)
  }
}  

########################################
### Perform DPI. Note that we now have a directed network
#i -> k  &   j -> k   &   i -> j
#then if MI suggests, can remove i->k
doaracneDPI <- function(net, tflist=unique(c(net$Regulator, net$Target))){
  net$Regulator <- as.character(net$Regulator)
  net$Target <- as.character(net$Target)
  onet <- net
  net$index <- 1:nrow(net)
  keep <- rep(TRUE, nrow(net))
  torem <- NULL
  for(i in 1:length(tflist)){
    tfi <- tflist[i]
    neti <- net[net$Regulator==tfi,]
    
    for(j in 1:nrow(neti)){
      tfj <- neti$Target[j]
      MI_ij <- neti$MI[j]
      netj <- net[net$Regulator==tfj  & net$Target %in% neti$Target,,drop=FALSE]
      
      if(nrow(netj)>0){
        compij <- smerge(
          data.frame(Target=netj$Target, MIj=netj$MI, indexj=netj$index),
          data.frame(Target=neti$Target, MIi=neti$MI, indexi=neti$index)
        )
        torem <- c(
          torem,
          compij$indexi[compij$MIj > compij$MIi]  #can kill i->j   but ignoring
          #          compij$indexi[compij$MIj > compij$MIi + MI_ij]  #can kill i->j
        )
        #        print(compij)
        #        print(c(tfi,tfj))
      }
    }
    
  }
  torem <- unique(torem)
  print(length(torem))
  print(nrow(onet))
  print(length(torem)/nrow(onet))
  onet[!(1:nrow(onet) %in% torem),]
  #net[keep,]
}


plotfilterAracneNeighN <- function(aracne_mergetc_atac){
  out <- c()
  for(n in 1:15){
    w <- filterAracneNeighN(aracne_mergetc_atac_red)
    score2mergetc <- sgenescorer2_matrix[sort(unique(c(w$Regulator, w$Target))),]
    out <- c(out,mean(apply(score2mergetc,1,min)<500, na.rm=TRUE))
  }
  out
  plot(out,type="l")
}






###########################
## Function: make network directed by removing redundant edges in opposite direction
directARACNe <- function(net){
  net[!(sprintf("%s-%s",net$Target, net$Regulator) %in% sprintf("%s-%s",net$Regulator, net$Target)),]
}


###########################
## Function: Cut down edges where there is no ATAC information. Do so by considering the distribution for each gene separately.
## Some genes are otherwise connected to almost everything and these take over
cutNonAtac <- function(net){
  
  cuthalfconns <- function(net){
    w <- sqldf("select distinct `Regulator`, median(`MI`) as a, count(`MI`) as b from net group by `Regulator`")
    w <- smerge(net,w)
    w[w$MI>w$a,]
  }
  
  net <- net[!(net$Regulator %in% map_jaspar_namegenesym$mgi_symbol | net$Target %in% map_jaspar_namegenesym$mgi_symbol),]
  print(nrow(net))  
  net <- undirectARACNe(net)
  print(nrow(net))  
  #keep 0.5^3 of the connections by removing 50% several times
  for(i in 1:3){
    net<-  cuthalfconns(net)
    net <- net[,c("Regulator","Target","MI","pvalue")]
    print(nrow(net))  
  }
  directARACNe(net)
}

###########################
## Function: Remove redundant edges by keeping the best ones
removerededge <- function(net){
  sqldf("select `Regulator`,`Target`,max(`MI`) as `MI`,min(pvalue) as pvalue from net group by `Regulator`,`Target`")
}


###########################
## Function: Keep interesting genes only
reducearacneinteresting <- function(net){
  interesting <- c(togenesym2(list_tf),
                   ensconvert$mgi_symbol[startsWith(ensconvert$mgi_symbol,"Il")],
                   ensconvert$mgi_symbol[startsWith(ensconvert$mgi_symbol,"Cxc")],
                   ensconvert$mgi_symbol[startsWith(ensconvert$mgi_symbol,"Ccr")])
  net <- net[net$Regulator %in% interesting & net$Target %in% interesting, ]
  print(nrow(net))
  net
}


###########################
## Function: Merge two networks by union
unionARACNE <- function(neta,netb){
  ### TODO: remove redundant edges should prioritize directed over non-directed
  directARACNe(removerededge(rbind(neta,netb)))
}

###########################
### Read net in ARACNe format
read.aracne_net <- function(file){
  x <- read.csv(file,sep="\t",stringsAsFactors = FALSE)
  x <- x[order(x$MI, decreasing = TRUE),]
  x <- x[order(x$pvalue),]
  x
}

################################################################
### Write an ARACNE network back to ARACNe format
write_aracne_net <- function(x, file){
  write.table(x, file,sep="\t", quote = FALSE, row.names = FALSE)
}



###############################################################################################
###### Read data ##############################################################################
###############################################################################################


aracne_difftc.mouse <- read.aracne_net("aracne/difftc.mouse.csv")
aracne_difftc.human <- read.aracne_net("aracne/difftc.human.csv")

#aracne_ko      <- read.aracne_net("aracne/ko.csv")
aracne_mousetc <- read.aracne_net("aracne/tc.mouse.csv")
aracne_humantc <- read.aracne_net("aracne/tc.human.csv")  
aracne_thymoma <- read.aracne_net("aracne/tcga_thymoma.csv")

aracne_ko123 <- rbind(
  read.aracne_net("aracne/ko1.csv"),
  read.aracne_net("aracne/ko2.csv"),
  read.aracne_net("aracne/ko3.csv")
)
aracne_ko123i <- intersectARACNe(
  intersectARACNe(read.aracne_net("aracne/ko1.csv"),read.aracne_net("aracne/ko2.csv")),
  read.aracne_net("aracne/ko3.csv"))

aracne_ko123c <- sqldf("select distinct `Regulator`,`Target`, max(`MI`) as MI, min(pvalue) as pvalue from aracne_ko123 group by `Regulator`,`Target` ")


###################
## Score matrix, indexed by ENSMUSG
sgenescorer2_matrix_ensmouse <- sgenescorer2_matrix
sgenescorer2_matrix_ensmouse$mgi_symbol <- rownames(sgenescorer2_matrix_ensmouse)
sgenescorer2_matrix_ensmouse <- smerge(sgenescorer2_matrix_ensmouse, ensconvert)
sgenescorer2_matrix_ensmouse <- sgenescorer2_matrix_ensmouse[
  sgenescorer2_matrix_ensmouse$ensembl_gene_id %in% rownames(tcmouse$mtpm) & 
    sgenescorer2_matrix_ensmouse$ensembl_gene_id!="ENSMUSG00000029268",]  #Ugt2a1 is annoying
rownames(sgenescorer2_matrix_ensmouse) <- sgenescorer2_matrix_ensmouse$ensembl_gene_id

###################
## Score matrix, indexed by ENSG
sgenescorer2_matrix_enshuman <- sgenescorer2_matrix_ensmouse
colnames(sgenescorer2_matrix_enshuman)[colnames(sgenescorer2_matrix_enshuman)=="ensembl_gene_id"] <- "ens_mouse"
sgenescorer2_matrix_enshuman <- smerge(sgenescorer2_matrix_enshuman, ortho_mouse_human)
sgenescorer2_matrix_enshuman <- sgenescorer2_matrix_enshuman[sgenescorer2_matrix_enshuman$ens_human %in% rownames(tchuman$mtpm),]
sgenescorer2_matrix_enshuman<-sqldf("select distinct avg(`Il4`) as `Il4`, avg(`Il13`) as `Il13`, avg(`Irf4`) as `Irf4`, avg(`Xbp1`) as `Xbp1`, avg(`Gata3`) as `Gata3`, ens_human from sgenescorer2_matrix_enshuman group by ens_human")
rownames(sgenescorer2_matrix_enshuman) <- sgenescorer2_matrix_enshuman$ens_human
#some 2000 genes affected. avg or min, not sure what is best
#sum(duplicated(sgenescorer2_matrix_enshuman$ens_human))
#colnames(sgenescorer2_matrix_enshuman)


###############################################################################################
###### Intersect networks #####################################################################
###############################################################################################

#############################
## Map network IDs from human to mouse. Need to be a 1-1 map
mapAracne2mouse <- function(net,ortho=ortho_mouse_human){
  m1 <- data.frame(Regulator=ortho$ens_human, Regulator_mouse=ortho$ens_mouse, stringsAsFactors = FALSE)
  m2 <- data.frame(Target   =ortho$ens_human, Target_mouse   =ortho$ens_mouse, stringsAsFactors = FALSE)
  net <- smerge(smerge(net,m1),m2)
  net <- net[,c("Regulator_mouse","Target_mouse","MI","pvalue")]
  colnames(net)[1:2] <- c("Target","Regulator")
  net
}
#############################
## Map network IDs from mouse to human Need to be a 1-1 map
mapAracne2human <- function(net,ortho=ortho_mouse_human){
  m1 <- data.frame(Regulator=ortho$ens_mouse, Regulator_human=ortho$ens_human, stringsAsFactors = FALSE)
  m2 <- data.frame(Target   =ortho$ens_mouse, Target_human   =ortho$ens_human, stringsAsFactors = FALSE)
  net <- smerge(smerge(net,m1),m2)
  net <- net[,c("Regulator_human","Target_human","MI","pvalue")]
  colnames(net)[1:2] <- c("Target","Regulator")
  net
}

# net <- mapAracne2mouse(aracne_tc)
# head(net)

### Conserved TC network
aracne_difftc.cons <- intersectARACNe(aracne_difftc.mouse,mapAracne2mouse(aracne_difftc.human))
aracne_tc.cons <- intersectARACNe(aracne_mousetc,mapAracne2mouse(aracne_humantc))
rm(aracne_humantc)
# [1] 1374128       4
# [1] 1698435       4
# [1] 8700    4


### ATAC on this
#aracne_tc.cons.atac <- filterAracneATACCHIP(aracne_tc.cons, atac.mouse$mapSiteInfo)
aracne_tc.mouse.atac <- filterAracneATACCHIP(aracne_mousetc, atac.mouse$mapSiteInfo)

nrow(aracne_tc.cons.atac)



nrow(mapAracne2mouse(aracne_difftc.human))
nrow(aracne_difftc.human)

nrow(ortho_mouse_human)


aracne_int_ko123tc <- intersectARACNe(aracne_mousetc, aracne_ko123c)



###############################################################################################
###### Store network graphs ###################################################################
###############################################################################################


###############################
## Function: store as GML
aracne2gml_colors <- c("Xbp1","Irf4","Il13","Il4","Gata3")
aracne2gml <- function(file, edges, nodes=NULL,keephitonly=FALSE,keepdeset=NULL,cutn=2, 
                       colorbyhit=aracne2gml_colors,cuthit=1000){ 
  psvec <- function(sdata) do.call(paste, c(as.list(sdata), sep=""))
  edges$Regulator <- as.character(edges$Regulator)
  edges$Target <- as.character(edges$Target)
  if(keephitonly){
    edges <- filterAracneHit(edges)
    edges <- filterAracneNeighN(edges,n = cutn)
  }
  
  ensconvert_ <- ensconvert
  #ensconvert_ <- ensconvert_[ensconvert_$ensembl_gene_id %in% c(edges$Regulator,edges$Target),]
  
  edges$Regulator <- togenesym(edges$Regulator, ensconvert_ = ensconvert_)
  edges$Target    <- togenesym(edges$Target,    ensconvert_ = ensconvert_)
  edges <- edges[edges$Target!="" & edges$Regulator!="",] #not quite sure how it gets in
  #print(tail(edges,n=100))
  isd <- duplicated(sprintf("%s_%s",edges$Regulator,edges$Target))
  edges <- edges[!isd,]
  
  #rownames(edges) <- NULL
  #rownames(edges) <- NULL
  
  if(!is.null(keepdeset)){
    hitlist <- c()
    for(i in list_screen_genes)
      hitlist <- c(hitlist, rownames(sgenescorer2_matrix)[sgenescorer2_matrix[,i]<cuthit])
    hitlist <- c(hitlist, togenesym(keepdeset, ensconvert_ = ensconvert_))
    edges <- edges[edges$Regulator %in% hitlist & edges$Target %in% hitlist,]    
  }
  
  if(is.null(nodes)){
    nodes<-data.frame(
      label=unique(c(edges$Regulator, edges$Target)),
      stringsAsFactors = FALSE)
  }
  #  print(head(nodes))
  nodes$id = 1:nrow(nodes)
  nodes$color = rep("#CCCCCC",nrow(nodes))
  
  uscol <- list()
  uscol[["Xbp1"]] <- "#00FFFF"
  uscol[["Irf4"]] <- "#FF55FF"
  uscol[["Il13"]] <- "#FFAAAA"
  uscol[["Il4"]] <- "#AAFFAA"
  uscol[["Gata3"]] <- "#AAAAFF"

  uscutoff <- list()
  for(i in list_screen_genes)
    uscutoff[[i]] <- 1000
  uscutoff[["Il4"]] <- 500

  for(i in names(uscol)){
    nodes$color[sgenescorer2_matrix[nodes$label,i ]<uscutoff[[i]]] <- uscol[[i]]
  }
  # nodes$color[sgenescorer2_matrix[nodes$label,"Xbp1" ]<1000] <- usc[["Xbp1"]]  #Xbp1  Cyan
  # nodes$color[sgenescorer2_matrix[nodes$label,"Irf4" ]<1000] <- "#"  #Irf4  Magenta
  # nodes$color[sgenescorer2_matrix[nodes$label,"Il13" ]<1000] <- "#FFAAAA"  #Il13  Red
  # nodes$color[sgenescorer2_matrix[nodes$label,"Il4"  ] <500] <- "#AAFFAA"  #Il4   Green
  # nodes$color[sgenescorer2_matrix[nodes$label,"Gata3"]<1000] <- "#AAAAFF"  #Gata3 Blue
  
  
  sec_node <- psvec(sprintf("node [\n  id %s\n  label \"%s\"\n graphics [ fill \"%s\" ]    ]\n", nodes$id, nodes$label, nodes$color))
  
  mapLabelId <- nodes
  rownames(mapLabelId) <- mapLabelId$label
  #  print(mapLabelId)
  # print(mapLabelId[51,1])
  print(nrow(edges))
  
  # print(nodes)
  id_regulator <- mapLabelId[edges$Regulator,]$id
  id_target    <- mapLabelId[edges$Target,]$id
  
  sec_edges <- psvec(sprintf(
    "edge [\n  source %s\n  target %s\n  MI %s\n  pvalue %s\n]\n",
    id_regulator, id_target, edges$MI, edges$pvalue))
  
  cat(
    paste("graph [  \ndirected 1\n", sec_node, "\n",sec_edges, "]", sep=""),
    file=file)
}
#aracne2gml("aracne/test.gml",aracne_int_ko123tc)

aracne2gml("aracne/f_tc.cons.gml",aracne_tc.cons,keephitonly = TRUE, keepdeset = allde$ens_mouse[allde$conserved_loose])

aracne2gml("aracne/f_tc.mouse.gml",aracne_mousetc,keephitonly = TRUE, keepdeset = allde$ens_mouse[allde$conserved_loose])

aracne2gml("aracne/f_tc.mouse.gml",aracne_mousetc,keephitonly = TRUE, keepdeset=c(666), cutn=8)

aracne2gml("aracne/f_tc.mouse.atac.de.gml",
           filterAracneDE(filterAracneATACCHIP(aracne_mousetc,ensconvert_ = ensconvert, tfc=atac.mouse$cons_tfc)), 
           keephitonly=TRUE, keepdeset=c(666))

aracne2gml("aracne/f_ko.mouse.gml",aracne_ko123c,keephitonly = TRUE, cutn=1)


filterAracneHit.loose <- function(edges, genes=list_screen_genes){
  hitlist <- c()
  for(i in genes){
    if(i=="Il4")
      hitlist <- c(hitlist, rownames(sgenescorer2_matrix)[sgenescorer2_matrix[,i]<500])
    else
      hitlist <- c(hitlist, rownames(sgenescorer2_matrix)[sgenescorer2_matrix[,i]<1000])
  }
  hitlist <- ensconvert$ensembl_gene_id[ensconvert$mgi_symbol %in% hitlist]
  edges <- edges[edges$Regulator %in% hitlist | edges$Target %in% hitlist,]    
  edges
}


aracne2gml("aracne/f_ko123i.mouse.gml",filterAracneHit.loose(aracne_ko123i),keephitonly = FALSE, cutn=1)





head(aracne_ko123i)
"Irf4" %in% togenesym(aracne_ko123i$Regulator)



filterAracneHit <- function(edges, genes=list_screen_genes,cutoff=500, isMouse=TRUE, reqboth=FALSE){
  hitlist <- c()
  if(isMouse){
    for(i in genes)
      hitlist <- c(hitlist, rownames(sgenescorer2_matrix_ensmouse)[sgenescorer2_matrix_ensmouse[,i]<cutoff])
  } else {
    for(i in genes)
      hitlist <- c(hitlist, rownames(sgenescorer2_matrix_enshuman)[sgenescorer2_matrix_enshuman[,i]<cutoff])
  }
  #  print(hitlist)
  if(reqboth)
    edges <- edges[edges$Regulator %in% hitlist & edges$Target %in% hitlist,]    
  else
    edges <- edges[edges$Regulator %in% hitlist | edges$Target %in% hitlist,]    
  edges
}


############### final networks in paper, regulons ####################
aracne2gml("aracne/f_tc_gata3.mouse.gml",
           filterAracneHit(genes = "Gata3",cutoff = 2000,reqboth = TRUE,    #barely used
             filterAracneATACCHIP(aracne_mousetc) 
             ),keepdeset = allde$ens_mouse[allde$conserved_loose],
           colorbyhit = c("Gata3",setdiff(aracne2gml_colors,"Gata3"))
          )
9999
# aracne2gml("aracne/f_tc_gata3.mouse.gml",
#            filterAracneHit(cutoff = 100,filterAracneHit(aracne_mousetc, genes = "Gata3",cutoff = 200)),keepdeset = allde$ens_mouse[allde$conserved_loose],
#            colorbyhit = c(aracne2gml_colors,"Gata3")
# )
# 
# aracne2gml("aracne/f_tc_gata3.mouse.gml",
#            filterAracneHit(reqboth = TRUE, cutoff = 1400,filterAracneHit(aracne_mousetc, genes = "Gata3",cutoff = 200)),keepdeset = allde$ens_mouse[allde$conserved_loose],
#            colorbyhit = c(aracne2gml_colors,"Gata3")
# )
# 
# aracne2gml("aracne/f_tc_gata3.mouse.gml",temp,
#            colorbyhit = c("Gata3",setdiff(aracne2gml_colors,"Gata3"))
#)

aracne2gml("aracne/f_tc_il4.mouse.gml",  ###this is what we kept
           filterAracneHit(filterAracneATACCHIP(aracne_mousetc), genes = "Il4",cutoff = 300),keepdeset = allde$ens_mouse[allde$conserved_loose],
           colorbyhit = c("Il4",setdiff(aracne2gml_colors,"Il4"))
)


aracne2gml("aracne/f_tc_gata3.cons.gml",
           filterAracneHit(aracne_tc.cons, genes = "Gata3",cutoff = 1000),keepdeset = allde$ens_mouse[allde$conserved_loose],
           colorbyhit = c(aracne2gml_colors,"Gata3")
)

#beeep. why is Spib not a hit, and yet colored as such?


aracne2gml("aracne/f_tc.mouse.atac.gml",aracne_tc.mouse.atac,keephitonly = TRUE, keepdeset=c(666))


aracne_mousetc[aracne_mousetc$Regulator==toensid("Gata3"),]
temp <- aracne_mousetc
temp <- filterAracneHit(aracne_mousetc, genes = "Gata3",cutoff = 500)
togenesym(temp[temp$Target==toensid("Gata3"),]$Regulator)
togenesym(temp[temp$Regulator==toensid("Gata3"),]$Target)

head(aracne_tc.cons)



net <- aracne_mousetc[aracne_mousetc$Regulator==toensid("Gata3") | aracne_mousetc$Target==toensid("Gata3"),]
aracne2gml("aracne/f_tc_regulon_gata3.mouse.gml",  
           filterAracneHit((net),cutoff = 1000),keepdeset = allde$ens_mouse[allde$conserved_loose],
           colorbyhit = c("Gata3",setdiff(aracne2gml_colors,"Gata3"))
)
#, genes = "Gata3"

net <- aracne_mousetc[aracne_mousetc$Regulator==toensid("Il4") | aracne_mousetc$Target==toensid("Il4"),]
aracne2gml("aracne/f_tc_regulon_il4.mouse.gml",  
           filterAracneHit((net), genes = "Il4",cutoff = 1000),keepdeset = allde$ens_mouse[allde$conserved_loose],
           colorbyhit = c("Il4",setdiff(aracne2gml_colors,"Il4"))
)


###used
net <- aracne_mousetc[aracne_mousetc$Regulator==toensid("Ube2m") | aracne_mousetc$Target==toensid("Ube2m"),]
aracne2gml("aracne/f_tc_regulon_ube2m.mouse.gml",  
           filterAracneHit(filterAracneATACCHIP(net), cutoff = 1000),keepdeset = allde$ens_mouse[allde$conserved_loose],
           colorbyhit = c(aracne2gml_colors)
)
thegene <- toensid("Ybx1")
net <- aracne_mousetc[aracne_mousetc$Regulator==thegene | aracne_mousetc$Target==thegene,]
aracne2gml("aracne/f_tc_regulon_ybx1.mouse.gml",  
           filterAracneHit(filterAracneATACCHIP(net), cutoff = 1000),keepdeset = allde$ens_mouse[allde$conserved_loose],
           colorbyhit = c(aracne2gml_colors)
)
thegene <- toensid("Nfkb2")   #shown network: removed non-hits
net <- aracne_mousetc[aracne_mousetc$Regulator==thegene | aracne_mousetc$Target==thegene,]
aracne2gml("aracne/f_tc_regulon_nfkb2.mouse.gml",  
           filterAracneHit(filterAracneATACCHIP(net), cutoff = 1000),keepdeset = allde$ens_mouse[allde$conserved_loose],
           colorbyhit = c(aracne2gml_colors)
)
thegene <- toensid("Rnpep")
net <- aracne_mousetc[aracne_mousetc$Regulator==thegene | aracne_mousetc$Target==thegene,]
aracne2gml("aracne/f_tc_regulon_rnpep.mouse.gml",  
           filterAracneHit(filterAracneATACCHIP(net), cutoff = 1000),keepdeset = allde$ens_mouse[allde$conserved_loose],
           colorbyhit = c(aracne2gml_colors)
)
thegene <- toensid("Sqstm1")
net <- aracne_mousetc[aracne_mousetc$Regulator==thegene | aracne_mousetc$Target==thegene,]
aracne2gml("aracne/f_tc_regulon_sqstm1.mouse.gml",  
           filterAracneHit(filterAracneATACCHIP(net), cutoff = 1000),keepdeset = allde$ens_mouse[allde$conserved_loose],
           colorbyhit = c(aracne2gml_colors)
)
thegene <- toensid("Irf8")
net <- aracne_mousetc[aracne_mousetc$Regulator==thegene | aracne_mousetc$Target==thegene,]
aracne2gml("aracne/f_tc_regulon_irf8.mouse.gml",  
           filterAracneHit(filterAracneATACCHIP(net), cutoff = 1000),keepdeset = allde$ens_mouse[allde$conserved_loose],
           colorbyhit = c(aracne2gml_colors)
)
thegene <- toensid("Ccdc134")   #shown network: removed non-hits
net <- aracne_mousetc[aracne_mousetc$Regulator==thegene | aracne_mousetc$Target==thegene,]
aracne2gml("aracne/f_tc_regulon_ccdc134.mouse.gml",  
           filterAracneHit((net), cutoff = 1000),keepdeset = allde$ens_mouse[allde$conserved_loose],
#           filterAracneHit(filterAracneATACCHIP(net), cutoff = 1000),keepdeset = allde$ens_mouse[allde$conserved_loose],
           colorbyhit = c(aracne2gml_colors)
)


net <- filterAracneHit(aracne_mousetc, cutoff = 500, reqboth = TRUE)
nn <- nametogenesym(sort(table(aracne_mousetc$Regulator)))
nn <- nn[order(names(nn))]
nn

###############################################################################################
#### Any genes iw  #################################################################
###############################################################################################

testRegulonScreenTarget <- function(net,gene,cutoff, ensconvert_=ensconvert){
  temp<-net
  intge <- toensid(  rownames(sgenescorer2_matrix)[sgenescorer2_matrix[,gene]<cutoff] , ensconvert_ = ensconvert_)
  temp2 <- undirectARACNe(temp[temp$Regulator %in% intge & temp$Target %in% intge,])
  temp3 <- undirectARACNe(temp[temp$Regulator %in% intge | temp$Target %in% intge,])
  v <- sqldf("select distinct `Regulator`, count(`Target`) as cde from temp2 group by `Regulator` ")
  v$mgi_symbol <- togenesym(v$Regulator, ensconvert_ = ensconvert_)
  v2 <- sqldf("select distinct `Regulator`, count(`Target`) as call from temp3 group by `Regulator` ")
  v <- smerge(v,v2)
  v$ratio <- v$cde/v$call
  
  allg <- unique(c(temp$Regulator,temp$Target))
  estp <- sum(allg %in% intge)/length(allg)

  v$pval <- NA
  for(i in 1:nrow(v)){
    v$pval[i] <- binom.test(v$cde[i],v$call[i],estp)$p.value
  }
  v$exp <- v$call*estp
  v$pval[v$exp>v$cde] <- 1
  v <- v[order(v$ratio,decreasing = TRUE),]
  v <- v[order(v$pval,decreasing = FALSE),]
  v  
}
v <- testRegulonScreenTarget(
  filterAracneDE(filterAracneATACCHIP(aracne_mousetc,atac.mouse$cons_tfc)),
  "Irf4",cutoff = 1000)
head(v[v$call>5,],n=20)

#Not working well
v <- testRegulonScreenTarget(
  filterAracneDE(filterAracneATACCHIP(aracne_tc.cons,atac.mouse$cons_tfc)),
  "Irf4",cutoff = 1000)
head(v[v$call>1,],n=20)


v <- testRegulonScreenTarget(
  ((aracne_tc.cons)),
  "Irf4",cutoff = 1000)
head(v[v$call>1,],n=20)
#Irf4. Ppp3r1

#Ppp3r1
888
temp <- filterAracneATACCHIP(aracne_humantc,tfc = atac.human$cons_tfc, ensconvert_ = human_ensconvert)
head(temp)  #this is really not working!
v <- testRegulonScreenTarget(
  filterAracneDE(temp, isMouse = FALSE),
  "Irf4",cutoff = 1000)
head(v[v$call>5,],n=20)



nrow(v)
#Ube2m -> Il4

# head(v,n=20)




###############################################################################################
#### Check network structure  #################################################################
###############################################################################################


###################
## Test if screen hits are in proximity - this version can also test between 2 screens
## Only works on mouse sgenescorer2_matrix_ensmouse
testproximityARACNe2 <- function(totest,thega="Il4",thegb="Gata3",cutoff=1000,nboot=200,scoremat=sgenescorer2_matrix_ensmouse){
  allg <- unique(c(totest$Regulator, totest$Target))
  ha <- scoremat[allg,thega]<cutoff
  hb <- scoremat[allg,thegb]<cutoff
  names(ha) <- allg
  names(hb) <- allg
  print(mean(is.na(ha)))
  ha[is.na(ha)]<-FALSE
  hb[is.na(hb)]<-FALSE
  scoref <- function(hx,hy){
    sum(hx[totest$Regulator] & hy[totest$Target]) + sum(hy[totest$Regulator] & hx[totest$Target])
  }
  b2 <- scoref(ha,hb)

  #Bootstrap over the same network but random assignments of which genes are hits
  bootb<-c()
  ha_orig <- ha
  hb_orig <- hb
  for(i in 1:nboot){
    # if(i%%100==0){
    #   print(i)
    # }
    ha <- sample(ha_orig)
    names(ha) <- allg
    if(thega==thegb){
      bootb <- c(bootb, scoref(ha,ha))
    } else {
      hb <- sample(hb_orig)
      names(hb) <- allg
      bootb <- c(bootb, scoref(ha,hb))
    }
  }
  mean(b2<bootb) #b is the number of neighbours. should be high. this score should be low
}
#testproximityARACNe <- function(totest,theg,cutoff) testproximityARACNe2(totest, theg, theg, cutoff)

###################
## Test if screen hits are in proximity: in this case, all combinations
testproximityARACNe2matrix <- function(){
  proximitymatrix <- matrix(nrow=length(list_screen_genes), ncol=length(list_screen_genes))
  colnames(proximitymatrix) <- list_screen_genes
  rownames(proximitymatrix) <- list_screen_genes
  for(i in 1:nrow(proximitymatrix)){
    for(j in i:nrow(proximitymatrix)){
      proximitymatrix[i,j] <- testproximityARACNe2(totest, list_screen_genes[i],list_screen_genes[j],  cutoff = 300)
      proximitymatrix[j,i] <- proximitymatrix[i,j]
    }
  }
  #heatmap(proximitymatrix)
  proximitymatrix
}

testproximityARACNe2list <- function(totest,cutoff=300,nboot=100,scoremat=sgenescorer2_matrix_ensmouse){
  out <- rep(0, length(list_screen_genes))
  names(out) <- list_screen_genes
  for(i in 1:nrow(proximitymatrix)){
    out[i] <- testproximityARACNe2(totest, list_screen_genes[i],  list_screen_genes[i], cutoff = cutoff, nboot=nboot, scoremat=scoremat)
  }
  out
}

testproximityARACNe2list(aracne_tc.cons,1000,nboot=1000)
#TODO

testproximityARACNe2list(aracne_ko123c,1000)
# Il4  Il13  Irf4  Xbp1 Gata3 
# 0.14  0.35  0.68  0.11  0.00 

testproximityARACNe2list(aracne_mousetc,1000,nboot=1000)
# Il4  Il13  Irf4  Xbp1 Gata3 
# 0.003 0.068 0.010 0.080 0.236 

testproximityARACNe2list(aracne_difftc.mouse,1000,nboot=1000)
# Il4  Il13  Irf4  Xbp1 Gata3 
# 0.034 0.078 0.305 0.054 0.007 

testproximityARACNe2list(aracne_thymoma,1000,nboot=1000,scoremat=sgenescorer2_matrix_enshuman)
# Il4  Il13  Irf4  Xbp1 Gata3 
# 1.000 1.000 1.000 0.993 0.996 






###############################################################################################
#### Check which regulons correlate with screens  #############################################
###############################################################################################


################################################################
## Given a named vector of genes and their scores, perform a t-test GSEA
performAracneGSEA <- function(net, checkde){
  checkde <- checkde[names(checkde) %in% net$Target]
  
  un <- unique(net$Regulator)
  unks <- rep(1, length(un))
  untp <- unks
  unt  <- unks
  na   <- unks
  unz  <- rep(0, length(un))
  for(i in 1:length(un)){
    if(i%%50==0)
      print(i)
    x <- checkde[names(checkde) %in% net$Target[net$Regulator==un[i]]]
    if(length(x)>3){
      #theks <- ks.test(x,checkde)
      #unks[i] <- theks$p.value
      thet    <- t.test(x,checkde, alternative = "greater")
      untp[i] <- thet$p.value
      unt[i]  <- thet$statistic
      unz[i]  <- (mean(x)-mean(checkde))/sd(x)      #if checkde rank, z should be positive
    }
    na[i] <- length(x)
  }
  en <- data.frame(ensembl_gene_id=un, pt=untp, z=unz, t=unt, na=na) #pKS=unks, 
  en <- smerge(en, rbind(ensconvert,human_ensconvert))
  en <- en[order(en$z,decreasing = FALSE),]
  en <- en[order(en$pt,decreasing = FALSE),]
  en
}


delistFromScreen <- function(theg, n=1000) {
  delist <- sgenescorer2_matrix_ensmouse[,theg]
  names(delist) <- rownames(sgenescorer2_matrix_ensmouse)
  delist[delist>n] <- n+sqrt(delist[delist>n]-n)
  delist
}

delist <- delistFromScreen("Gata3")
gsea_aracne_ko123c  <- performAracneGSEA(aracne_ko123c,  delist)   
gsea_aracne_tc.mouse <- performAracneGSEA(aracne_mousetc, delist)
gsea_aracne_tc.human <- performAracneGSEA(aracne_humantc, delist)
gsea_aracne_difftc.mouse <- performAracneGSEA(aracne_difftc.mouse, delist) #Il4: Zfp51..Nfkb2


head(gsea_aracne_ko123c,n=50)  
#Il4: Foxo4..
head(gsea_aracne_tc.mouse,n=50)
head(gsea_aracne_difftc.mouse,n=50)
#Gata3: Stat6 scores high


##########################
## Function: Plot how the GSEA compares for two networks (input: 2 gsea's)
plot2aracneGSEA <- function(neta,netb,usez=TRUE){
  neta <- neta[neta$na>50,]
  netb <- netb[netb$na>50,]
  if(usez){
    m <- smerge(
      data.frame(ensembl_gene_id=neta$ensembl_gene_id,  pa=neta$z),
      data.frame(ensembl_gene_id=netb$ensembl_gene_id,  pb=netb$z), all=TRUE)
    m$pa[is.na(m$pa)] <- 0
    m$pb[is.na(m$pb)] <- 0
    m <- m[m$pa<0 & m$pb<0,]
    qtextscatter(m$pa, m$pb, m$ensembl_gene_id)
  } else {
    m <- smerge(
      data.frame(mgi_symbol=neta$ensembl_gene_id,  pa=neta$pvalue),
      data.frame(mgi_symbol=netb$ensembl_gene_id,  pb=netb$pvalue), all=TRUE)
    m$pa[is.na(m$pa)] <- 1
    m$pb[is.na(m$pb)] <- 1
    qtextscatter(log10(m$pa), log10(m$pb), m$ensembl_gene_id)
  }
}

# plot2aracneGSEA(aracne_mousetc_gsea, aracne_humantc_gsea)
# plot2aracneGSEA(aracne_mousetc_gsea, aracne_ko_gsea)


################################################################
## Function: Compare ranks of TFs with the enrichment of genes downstream in a particular screen
testScreenVsOwnRegulon <- function(net, thescreen, rankcutoff=2000){
  if(TRUE){ #Use mageck
    delist <- delistFromScreen(thescreen, rankcutoff)
    # delist <- sgenescorer2_matrix_ensmouse[,thescreen]
    # names(delist) <- rownames(sgenescorer2_matrix_ensmouse)
    # delist <- smoothcutrank(delist,rankcutoff)
    allv <- data.frame(ensembl_gene_id=rownames(sgenescorer2_matrix_ensmouse),ss=sgenescorer2_matrix_ensmouse[,thescreen])
  } else {  #Use STAN
    # delist <- log10(newscreen_pval)[,str_to_lower(thescreen)]
    # names(delist) <- rownames(newscreen_pval)
    # allv <- data.frame(mgi_symbol=rownames(newscreen_pval),ss=newscreen_pval[,str_to_lower(thescreen)])
  }
  netgsea <- performAracneGSEA(net, delist)
  temp <-smerge(netgsea,allv)
  temp[order(temp$pt),]
}
testScreenVsOwnRegulonAll <- function(net, thescreen, rankcutoff=1000, ensconvert_=ensconvert){
  if(FALSE){
    #Only DE
    # gde <- unique(c(de_early$mgi_symbol[de_early$qval<qval],de_late$mgi_symbol[de_late$qval<qval]))
    # net <- net[net$Target %in% gde,]
  }  
  temp <- testScreenVsOwnRegulon(net, thescreen, rankcutoff)
  temp[temp$ensembl_gene_id==toensid(thescreen, ensconvert_=ensconvert_),,drop=FALSE]
}
testScreenVsOwnRegulonAll4net <- function(net, ensconvert_=ensconvert){
  #Here testing screen rank. want negative z
  rbindlist(lapply(list_screen_genes, function(x) testScreenVsOwnRegulonAll(net, x, ensconvert_=ensconvert_)))
}


svsownreg_ko123c  <- testScreenVsOwnRegulonAll4net(aracne_ko123c) 
svsownreg_tcmouse <- testScreenVsOwnRegulonAll4net(aracne_mousetc)
svsownreg_difftcmouse <- testScreenVsOwnRegulonAll4net(aracne_difftc.mouse)
svsownreg_tchuman <- testScreenVsOwnRegulonAll4net(aracne_thymoma, ensconvert_=human_ensconvert) 
#svsownreg_thymoma <- testScreenVsOwnRegulonAll4net(aracne_thymoma)  #same


################################################################
## Check if genes targeting X are enriched in screen X
ttest4upstreamaracne <- function(net, thescreen, mat=sgenescorer2_matrix_ensmouse, ensconvert_=ensconvert){
#  print(thescreen)
#  print(toensid(thescreen, ensconvert_))
  snet <- net[net$Target==toensid(thescreen, ensconvert_),]
#  print(head(snet))
  v <- smerge(
    net[net$Target==toensid(thescreen, ensconvert_),],
    data.frame(Regulator=rownames(mat),ss=mat[,thescreen])
  )  
#  print(head(v))
  if(nrow(v)>3){
    tt <- t.test(v$ss, mat[,thescreen])   #want first estimate to be lower
    if(tt$estimate[1]<tt$estimate[2])
      tt$p.value
    else
      -tt$p.value
  }
  else
    666
}
ttest4upstreamaracne4net <- function(net, mat=sgenescorer2_matrix_ensmouse, ensconvert_=ensconvert){
  v<-as.double((lapply(list_screen_genes, function(x) ttest4upstreamaracne(net, x, mat=mat, ensconvert_=ensconvert_))))
  names(v) <- list_screen_genes
  v
}


ttest4upstreamaracne4net(aracne_ko123c)
# Il4            Il13         Irf4           Xbp1         Gata3 
# -0.1297281     0.1084781    0.9515147      0.7894892    0.1027596 
ttest4upstreamaracne4net(aracne_mousetc)
# Il4            Il13         Irf4           Xbp1         Gata3 
# 0.99371199     0.00418827   0.15046471     666.00000000 666.00000000 

ttest4upstreamaracne4net(undirectARACNe(aracne_ko123c))
# Il4            Il13          Irf4          Xbp1          Gata3 
# -0.1297280542  0.1084780948  0.3768160511  0.1797998316  0.0006636917 
ttest4upstreamaracne4net(undirectARACNe(aracne_mousetc))
# Il4            Il13        Irf4            Xbp1          Gata3 
# 0.99371199     0.00418827  0.82353971      0.01073388    -0.78152771 
ttest4upstreamaracne4net(undirectARACNe(aracne_tc.cons))
# Il4        Il13        Irf4        Xbp1       Gata3 
# 666.0000000 666.0000000   0.7381014 666.0000000 666.0000000 

ttest4upstreamaracne4net(undirectARACNe(aracne_humantc),mat=sgenescorer2_matrix_enshuman, ensconvert_=human_ensconvert)



################################################################
## test specifically cor of the core TFs vs screen hits - all genes
testTFcorTarget <- function(x, gene, ensconvert_=ensconvert, mintpm=10){
  nexpressedGenesId <- rownames(x)[apply(x>mintpm,1,any)]  
  
  mtpm_red <- x[rownames(x) %in% nexpressedGenesId,]
  geneid <- toensid2(gene, ensconvert_=ensconvert_)
  genelev <- as.double(x[geneid,])
  #  print(genelev)
  
  v <- as.double(cor(genelev, t(mtpm_red),method = "spearman"))
  v[abs(v)<0.2] <- 0
  
  # s1 <- ensconvert_$ensembl_gene_id[ ensconvert_$mgi_symbol %in% rownames(sgenescorer2_matrix)[sgenescorer2_matrix[,gene]<1000] ]
  # s2 <- ensconvert_$ensembl_gene_id[ ensconvert_$mgi_symbol %in% rownames(sgenescorer2_matrix)[sgenescorer2_matrix[,gene]>5000] ]
  s1 <- rownames(sgenescorer2_matrix_ensmouse)[sgenescorer2_matrix_ensmouse[,gene]<100]  #requires mouse
  s2 <- rownames(sgenescorer2_matrix_ensmouse)[sgenescorer2_matrix_ensmouse[,gene]>2000]
  
  vn <- rownames(mtpm_red)
  i1 <- vn %in% s1
  i2 <- vn %in% s2
  t.test(
    abs(v[i1]),
    abs(v[i2])
  )
}
testTFcorTargetAll <- function(x, ensconvert_=ensconvert){
  alltestTFcorTarget <- lapply(list_screen_genes, function(gene) testTFcorTarget(x, gene, ensconvert_ = ensconvert_))
  data.frame(
    gene=list_screen_genes,
    meanCorrHits=sapply(alltestTFcorTarget, function(x) x$estimate[1]),
    meanCorrNonhits=sapply(alltestTFcorTarget, function(x) x$estimate[2]),
    dc_should_be_pos=sapply(alltestTFcorTarget, function(x) x$estimate[1]) - sapply(alltestTFcorTarget, function(x) x$estimate[2]),
    pvalue=sapply(alltestTFcorTarget, function(x) x$p.value)
  )
}
testTFcorTargetAll(tcmouse$mtpm)  #not impressive
testTFcorTargetAll(tchuman$mtpm, human_ensconvert)  #will not work
testTFcorTargetAll(ncount)
testTFcorTargetAll(calcth2sub0(tcmouse)) #interesting




###############################################################################################
#### Check which regulons seem interesting  ###################################################
###############################################################################################

###Of the main genes mentioned in http://www.nejm.org/doi/full/10.1056/NEJMoa1301689#t=article
#Cebpa mutated in TCGA cancer, hit-ish in IL4
#Cbfb fusion is another possibility. hit in Il13 and Xbp1
#Nup98 is almost an Il4 hit

######## Which TFs overlap?
net <- aracne_ko
un <- unique(net$Regulator)
tfo <- matrix(nrow=length(un), ncol=length(un))
tfr <- vector("list",length(un))
for(i in 1:length(un)){
  tfr[[i]] <- net$Target[net$Regulator==un[i]]
}
for(i in 1:length(un)){
  for(j in 1:length(un)){
    tfo[i,j] <- length(intersect(tfr[[i]], tfr[[j]]))
  }
}
colnames(tfo) <- un
rownames(tfo) <- un




########################## using the KO network ############################
############## These are interesting  based on de_late
#not Srf Carhsp1 Nr1h2 Setdb2
sgenescorer2_matrix["Smarcc2",] #meah
sgenescorer2_matrix["Zscan22",] #on the limit
sgenescorer2_matrix["Dmbx1",] #irf4
sgenescorer2_matrix["Foxo3",] #il4!
sgenescorer2_matrix["Dbp",] #gata3! il4!   partially DE itself
sgenescorer2_matrix["Gata3",] #gata3! xbp!
sgenescorer2_matrix["Zfp113",] #gata3!
### early
#sgenescorer2_matrix["Nfe2",]
#sgenescorer2_matrix["Usf3",]
sgenescorer2_matrix["Nrf1",] #il4!
sgenescorer2_matrix["Zfp408",] #irf!

sgenescorer2_matrix["Cbfb",] #good Z score. Il13 and Xbp1!
sgenescorer2_matrix["Zfp97",] 
sgenescorer2_matrix["Zfp955a",] #irf4! in commonish TC/KO
sgenescorer2_matrix["Stat4",] 




###############################################################################################
###### Plot interesting genes TC ##############################################################
###############################################################################################



tcmouse$mtpm_th2[toensid("Sqstm1"),]
  
plotTpmTC20 <- function(gname, mouse=TRUE, title){
  pcol <- c("#00FF00","#FF0000")
  tp <- c("0h","0.5h","1h","2h","4h","6h","12h","24h","48h","72h")
  lwd <- 2
  if(mouse){
    gid <- toensid(gname)
    usetc <- tcmouse
    lev2 <- as.double(usetc$av_mtpm[gid,])
    lev0 <- as.double(usetc$av_mtpm0[gid,])
  } else {
    gid <- toensid2(gname,ensconvert_ = human_ensconvert)
    usetc <- tchuman
    lev2 <- as.double(usetc$av_mtpm[gid,])
    lev0 <- as.double(usetc$av_mtpm0[gid,])
  }
  ymax <- max(c(lev0,lev2))
  plot(1:10,lev0,col=pcol[1],  xaxt = "n", type="l",
       ylab="TPM",xlab="",ylim=c(0,ymax),lwd=lwd, main=title) 
  lines(1:10,lev2,type="l",col=pcol[2],lwd=lwd) 
  legend(6,ymax,c("Th0","Th2"),cex=1.2,fill = pcol,y.intersp=0.65,box.lwd = 0)
  axis(1, at = 1:length(tp), labels=tp)
}

pdf("aracne/tc.pdf",width = length(genelist)*4, height=6)
genelist <- c("Sqstm1","Hk1","Rnpep","Irf8","Ube2m","Ybx1","Gata3","Il4","Tbx21","Xbp1","Irf4","Batf")
par(mfrow=c(2,length(genelist)))
for(formouse in c(TRUE,FALSE)){
  for(g in genelist){
    plotTpmTC20(g, mouse=formouse, g)
  }
}
dev.off()


plotTpmTC20("Sqstm1", mouse=formouse)
plotTpmTC20("Hk1",mouse=formouse)
plotTpmTC20("Rnpep",mouse=formouse)
plotTpmTC20("Irf8",mouse=formouse)   #in human fairly different. and mouse. different ways
plotTpmTC20("Ube2m",mouse=formouse)
plotTpmTC20("Xbp1",mouse=formouse)
plotTpmTC20("Ybx1",mouse=formouse)
plotTpmTC20("Gata3",mouse=formouse)
plotTpmTC20("Il4",mouse=formouse)


plotTpmTC20("Hk1",mouse=FALSE)
