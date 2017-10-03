require(biomaRt)


######################################################################
### Read ensembl #####################################################
######################################################################
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
ensembl_genes <- as.data.frame(getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'gene_biotype'), mart=ensembl),stringsAsFactors=FALSE)
rownames(ensembl_genes) <- ensembl_genes$ensembl_gene_id
mt_genes <- ensembl_genes[which(ensembl_genes$gene_biotype=="Mt_rRNA" | ensembl_genes$gene_biotype=="Mt_tRNA"),]

human_ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

#Convert to gene symbols, or retain ID if no symbol
ensconvert <- as.data.frame(getBM(attributes=c('ensembl_gene_id', 'mgi_symbol'), mart=ensembl),stringsAsFactors=FALSE)

human_ensconvert <- as.data.frame(getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol'), mart=human_ensembl),stringsAsFactors=FALSE)
colnames(human_ensconvert)[2] <- "mgi_symbol"  #should use a different name. genesymbol
v<-sqldf("select *,count(mgi_symbol) as c from human_ensconvert group by ensembl_gene_id")
human_ensconvert <- human_ensconvert[human_ensconvert$ensembl_gene_id %in% v$ensembl_gene_id[v$c==1],] #there are some bastards including CCL3L3. this one better inserted by hand
human_ensconvert$mgi_symbol <- normalizesym(human_ensconvert$mgi_symbol)





nametogenesym<-function(x){
  names(x)<-togenesym(names(x))
  x
}

mtogenesym<-function(x){
  colnames(x)<-togenesym(colnames(x))
  rownames(x)<-togenesym(rownames(x))
  x
}






togenesym2 <- function(geneid, ensconvert_=ensconvert, dowarn=TRUE){
  ensconvert_ <- ensconvert_[ensconvert_$ensembl_gene_id %in% geneid,]
  ensconvert_ <- ensconvert_[order(ensconvert_$mgi_symbol),]
  ensconvert_ <- ensconvert_[order(ensconvert_$ensembl_gene_id),]
  ensconvert_ <- ensconvert_[!duplicated(ensconvert_$ensembl_gene_id, fromLast=TRUE),]
  rownames(ensconvert_) <- ensconvert_$ensembl_gene_id

  out <- ensconvert_[geneid,]$mgi_symbol
  out[is.na(out)] <- geneid[is.na(out)]
  out  
}
togenesym <- togenesym2
 
toensid2 <- function(geneid, ensconvert_=ensconvert, dowarn=TRUE){
  ensconvert_ <- ensconvert_[ensconvert_$mgi_symbol %in% geneid,]
  ensconvert_ <- ensconvert_[order(ensconvert_$mgi_symbol),]
  ensconvert_ <- ensconvert_[order(ensconvert_$ensembl_gene_id),]
  ensconvert_ <- ensconvert_[!duplicated(ensconvert_$mgi_symbol, fromLast=TRUE),]
  rownames(ensconvert_) <- ensconvert_$mgi_symbol
  
  out <- ensconvert_[geneid,]$ensembl_gene_id
  out[is.na(out)] <- geneid[is.na(out)]
  out  
}
toensid <- toensid2


########################################
## Orthology table human<->mouse
ortho_mouse_human <- read.csv("mouse_human_ortholog.csv",stringsAsFactors = FALSE)
colnames(ortho_mouse_human) <- c("ens_mouse","ens_human")
ortho_mouse_human <- ortho_mouse_human[ortho_mouse_human$ens_human!="" & ortho_mouse_human$ens_mouse!="",]
## Only 1-1 mappings
ortho_mouse_human_unique <- ortho_mouse_human[
  isUnique(ortho_mouse_human$ens_human) & 
  isUnique(ortho_mouse_human$ens_mouse),]
sum(duplicated(ortho_mouse_human_unique$ens_mouse)) 
sum(duplicated(ortho_mouse_human_unique$ens_human)) 


