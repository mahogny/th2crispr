library(gProfileR)
library(topGO)
library(GO.db)
mapGoTerm <- toTable(GOTERM) 

require(org.Mm.eg.db)  # generalize
#require(org.Hs.eg.db)  # generalize

annotation.mm <- list(
  mapping="org.Mm.eg.db",
  db=org.Mm.eg,
  dbSYMBOL=org.Mm.egSYMBOL,
  egGO=org.Mm.egGO,
  gene2ez=revmap(org.Mm.egSYMBOL)#annotation.dbSYMBOL)
)

# annotation.hs <- list(
#   mapping="org.Hs.eg.db",
#   db=org.Hs.eg,
#   dbSYMBOL=org.Hs.egSYMBOL,
#   egGO=org.Hs.egGO,
#   gene2ez=revmap(org.Hs.egSYMBOL)#annotation.dbSYMBOL)
# )


# source('http://www.bioconductor.org/biocLite.R') 
# biocLite(c('topGO', 'org.Mm.eg.db'))
#biocLite(c('biomaRt','topGO', 'org.Mm.eg.db'))
                                                            
######################################################################
### GO functions #####################################################
######################################################################

### Simplified use of topgo and gprofiler
stopgo <- function(all_data,DE_data,ID=c("genename","symbol","EnsemblID"),useTopgo=TRUE, keepCol=FALSE,cutoff=1e-4){
  if(useTopgo){
    relevant.genes <- rep(1,length(all_data))
    relevant.genes[all_data %in% DE_data] <- 0
    names(relevant.genes) <- all_data
    relevant.genes <- as.factor(relevant.genes)
    
    GOdata.BP <- new("topGOdata", ontology='BP', allGenes = relevant.genes,nodeSize=5,annot = annFUN.org,
                     mapping="org.Mm.eg.db", ID=ID, geneSel = function(p) p < 0.01)
    # apply Fisher's exact test with elimination mode:
    results <- runTest(GOdata.BP, algorithm = 'elim', statistic = 'fisher')
    x<-GenTable(GOdata.BP, results, topNodes = 60)#, pvalCutOff=0.1)   #was 40
    
    #### TODO: rename things
    colnames(x)[which(colnames(x)=="Term")] <- "term.name"
    colnames(x)[which(colnames(x)=="result1")] <- "p.value"
    x    
  } else {
    x<-gprofiler(
      DE_data,organism = "mmusculus", ordered_query = FALSE,
      correction_method = "fdr",
      custom_bg = all_data,
      hier_filtering = "moderate")
    x<-x[order(x$p.value),]
    x<-x[x$p.value < cutoff,]
    if(nofactor){
      x<-x[grep("Factor",x$term.name,invert = TRUE),]
    }
    if(!keepcol)
      x<-x[,c("p.value","term.name")]
    x
  }
}


stopgosym <- function(genelist,bg=unique(ensconvert$mgi_symbol),nofactor=FALSE, useTopgo=TRUE,cutoff=1e-4){
  stopgo(bg, genelist, "symbol", useTopgo = useTopgo, cutoff = cutoff) #or genename
}


combinegolist2matrix <- function(outgo, usenames){
  allcat <- c()
  for(i in 1:length(outgo))
    allcat <- union(allcat, outgo[[i]]$term.name)    #Also got GO.ID (in topgo)
  screengo <- matrix(1,nrow=length(allcat),ncol=length(outgo))
  rownames(screengo) <- allcat
  colnames(screengo) <- usenames #names(outgo)
  for(i in 1:length(outgo)){
    print(outgo[[i]])
    screengo[outgo[[i]]$term.name,i] <- as.double(outgo[[i]]$p.value)
    
  }
  #screengo[is.na(screengo)] <- 1
  screengo
}


########################
#### Calculate a matrix for each screen for each 
recalcgomatrix <- function(listfg, listbg, terms){
  listgo <- unique(mapGoTerm$go_id[mapGoTerm$Term %in% terms])
  #print(listgo)
  
  outgo <- matrix(NA,ncol=length(listfg), nrow=length(listgo))
  for(j in 1:length(listgo)){
    #print(listgo[j])
    genesinterm <- go2genes(listgo[j],with.descendants = TRUE, annotation = annotation.mm)
    for(i in 1:length(listfg)){
      #      fg <- rownames(newscreen_rank)[newscreen_rank[,i]<n]
      #     bg <- rownames(newscreen_rank)[newscreen_inc[,i]>0]
      outgo[j,i] <- testoneGO(listfg[[i]], listbg[[i]], genesinterm)
    }
  }
  colnames(outgo) <- names(listfg)
  rownames(outgo) <- listgo
  newterm <- NULL
  for(go in listgo)
    newterm <- c(newterm,mapGoTerm$Term[mapGoTerm$go_id==go][1])
  rownames(outgo) <- newterm
  outgo
}


#############
getgomatrix <- function(listfg, listbg, listnames=NULL, cutoff=1e-2, useTopgo=TRUE){
  outgo <- list()
  for(i in 1:length(listfg)){
    print(i)
    outgo[[i]] <- stopgosym(listfg[[i]],bg = listbg[[i]], nofactor = TRUE, cutoff=cutoff, useTopgo)
    print(outgo[[i]]) #entertainment
  }
  print(outgo)
  if(!is.null(listnames))
    listnames <- names(listfg)
  screengo <- combinegolist2matrix(outgo, listnames)
  colnames(screengo) <- listnames
  screengo  
  
  recalcgomatrix(listfg, listbg, rownames(screengo))
}



filtergomatriximmuno <- function(go){
  rn <- rownames(go)
  minp <- apply(go,1,min)
  keepgenes <- unique(c(
    #which(minp<1e-2),
    grep("T cell",rn),
    grep("B cell",rn),
    grep("immuno",rn),
    grep("MHC",rn),
    grep("TOR",rn),
    grep("inflam",rn),
    grep("immune",rn),
    grep("STAT",rn),
    grep("JAK",rn),
    grep("MAP",rn),
    grep("interleukin",rn),
    grep("platelet",rn),
    grep("vitamin",rn),  #Vitamin A enhances Th2  https://www.ncbi.nlm.nih.gov/pubmed/11970994
    grep("virus",rn),
    grep("cyt",rn),
    grep("MHC",rn),
    grep("calcium",rn)
    
    #grep("response",rn))
  ))
  go <- go[keepgenes,]
  minp <- apply(go,1,min)
  go <- go[minp<1,]
  go
}

plotgomatrix <- function(go, cutoff=-2.5,keepimmuno=TRUE, cutoffupper=-1000){
  go <- log10(go)
#  go <- go[grep("biological process",rownames(go),invert=TRUE),]
  go <- go[grep("biological",rownames(go),invert=TRUE),]
  go <- go[grep("intracellular",rownames(go),invert=TRUE),]
  go <- go[grep("MI",rownames(go),invert=TRUE),]
  go <- go[(apply(go,1,min))<cutoff,]
  go <- go[(apply(go,1,min))>=cutoffupper,]
  if(keepimmuno){
    go <- filtergomatriximmuno(go)
  }
  
  print(max(go))
  print(min(go))
  heatmap.2(
    trace="none",
    density.info = "none",
    go,
    margins = c(4,15),
    
    #Colv = FALSE, dendrogram="row",
    col=colorRampPalette(c("#7a0177","#c51b8a","#f768a1","#fbb4b9","#feebe2",
                           "#feebe2","blue","blue","blue","blue"))(n = 330))
  #  col=colorRampPalette(c("red", "yellow","white"))(n = 300))
}





#from here: https://gist.github.com/massyah/6029531

go2genes <- function(go,with.descendants=FALSE, annotation){
  genes=c()
  if(with.descendants){
    #Check if recursive
    goids<-go.descendants(go)
  }else{
    goids<-c(go)
  }
  ezgenes<-unique(unlist(mget(goids,revmap(annotation$egGO),ifnotfound=NA)))
  ezgenes<-ezgenes[!is.na(annotation.ezgenes)]
  symbols<-unique(unlist(mget(ezgenes,annotation$dbSYMBOL)))
  return(symbols)
}


gene2gos <- function(gene,flatten=T, annotation){
  ezgenes<-unique(unlist(mget(gene,annotation$gene2ez,ifnotfound=NA)))
  ezgenes<-ezgenes[!is.na(annotation$ezgenes)]
  if(flatten){
    gos<-unlist(mget(ezgenes,annotation$egGO),recursive=F,use.names=F)
    goids=unique(sapply(gos,"[[","GOID"))
    return(goids)
  }else{
    res<-melt(mget(annotation$ezgenes,annotation$egGO))
    res<-res[res$"L3"=="GOID",]
    res["gene"]=unlist(mget(res$"L1",annotation$dbSYMBOL))
    res["term"]=go2terms(as.character(res$L2))
    res<-res[,c("L1","gene","L2","term")]
    colnames(res)<-c("ezgene","gene","GO","term")
    rownames(res)<-NULL
    return(unique(res))
  }
}



#### Test one particular GO term, fisher stas
testoneGO <- function(mygeneset, backgroundset, genesinterm){
  n <- length(intersect(mygeneset, genesinterm))
  m <- length(setdiff(mygeneset, genesinterm))
  tm <- matrix(c(
    n, length(genesinterm)-n,
    m, length(backgroundset)-m
  ),nrow=2)
  fisher.test(tm)$p.value
}



###############
## Function: Get descendant GO categories
go.descendants <- function(go,onto=c("BP","MF","CC")){
  onto=match.arg(onto)
  if(onto=="BP"){
    onto=GOBPOFFSPRING
  }else if(onto=="MF"){
    onto=GOMFOFFSPRING
  }else if(onto=="CC"){
    onto=GOCCOFFSPRING
  }
  children<-unique(unlist(mget(go,onto,ifnotfound=NA)))
  children=union(children,go)
  return(children[complete.cases(children)])
}
