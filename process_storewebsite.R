#A screen can be thought of as a sample, with FC and mageck p-value
#But further, a gene will have TPM and DE scores

######################################################################
### Store time-course data ###########################################
######################################################################

maketpmmeta <- function(mtpm, species){
  tph<-c("05h",sprintf("%sh",c(1,2,4,6,12,24,48,72)))
  tpd<-c(0.5,1,2,4,6,12,24,48,72)
  
  mtpm_samplemeta <- data.frame()
  for(thtype in c(0,2)){
    for(rep in 1:3){
      tph2 <- c(sprintf("Naive_%s_salmon",rep),sprintf("Th%s_%s_%s_salmon",thtype,tph,rep))
      tpd2 <- c(0,tpd)
      mtpm_samplemeta_part <- data.frame(sample=tph2, hours=tpd2, rep=rep, celltype=sprintf("Th%s",thtype),organism=species)
      if(nrow(mtpm_samplemeta)==0)
        mtpm_samplemeta <- mtpm_samplemeta_part
      else
        mtpm_samplemeta <- rbind(mtpm_samplemeta, mtpm_samplemeta_part)
    }
  }
  mtpm_samplemeta <- mtpm_samplemeta[mtpm_samplemeta$sample %in% colnames(mtpm),]
#  print(mtpm_samplemeta)
  mtpm_samplemeta$sample <- maketpmprefix(mtpm_samplemeta$sample,species)
  mtpm_samplemeta
}
maketpmprefix <- function(x, species){
  sprintf("th2crispr_%s_%s",species,x)
}
addtpmprefix <- function(mtpm, species){
  colnames(mtpm) <- maketpmprefix(colnames(mtpm),species)
  mtpm
}


write.csv(maketpmmeta(tcmouse$mtpm,"mouse"), "out_teichlab/th2crispr_mouse_tc_samplemeta.csv", row.names = FALSE, quote = FALSE)
write.csv(addtpmprefix(tcmouse$mtpm,"mouse"),"out_teichlab/th2crispr_mouse_tc_data.csv", row.names = TRUE, quote = FALSE)

write.csv(maketpmmeta(tchuman$mtpm,"human"), "out_teichlab/th2crispr_human_tc_samplemeta.csv", row.names = FALSE, quote = FALSE)
write.csv(addtpmprefix(tchuman$mtpm,"human"),"out_teichlab/th2crispr_human_tc_data.csv", row.names = TRUE, quote = FALSE)


######################################################################
### Store averaged time-course data ##################################
######################################################################

######### Sample description #################
write_teichlab_tcavg <- function(av_mtpm2, av_mtpm0, organism){
  colnames(av_mtpm0) <- sprintf("%s_tcavg_th0_%s",organism,colnames(av_mtpm2)) #on purpose!
  colnames(av_mtpm2) <- sprintf("%s_tcavg_th2_%s",organism,colnames(av_mtpm2))
  av_mtpm <- cbind(av_mtpm2, av_mtpm0)
  samplemeta <- data.frame(sample=colnames(av_mtpm))
  samplemeta$organism<-"mouse"
  samplemeta$'Cell Type'<-"Th2"  
  samplemeta$'Cell Type'[grep("th0",samplemeta$sample)] <- "Th0"
  samplemeta$method<-"RNAseq"
  samplemeta$hours <-c(0,0.5,1,2,4,6,12,24,48,72)

  #sgenescorer2_samplemeta$target_mgi_symbol<-res_samplemeta$sample
  write.csv(samplemeta,sprintf("out_teichlab/th2crispr_%s_tcavg_samplemeta.csv",organism),row.names = FALSE, quote = FALSE)
  write.csv(av_mtpm,   sprintf("out_teichlab/th2crispr_%s_tcavg_data.csv",organism), row.names = TRUE, quote = FALSE)
}

write_teichlab_tcavg(tcmouse$av_mtpm, tcmouse$av_mtpm0, "mouse")
write_teichlab_tcavg(tchuman$av_mtpm, tchuman$av_mtpm0, "human")




######################################################################
### Store screen data ################################################
######################################################################


##### the score. remove annoying rows
tempde <- rbind(tcmouse$de_early, tcmouse$de_late)
tempde <- sqldf("select distinct mgi_symbol, min(pval) as pval from tempde group by mgi_symbol")[-1,]
rownames(tempde) <- tempde$mgi_symbol

temp <- sgenescorer2_matrix
temp$pDE <- tempde[rownames(temp),]$pval
temp$pDE[is.na(temp$pDE)] <- 1
temp$pDE <- round(log10(temp$pDE),digits=3)
colnames(temp) <- sprintf("cr2_%s",colnames(temp))

te <- ensconvert
te <- te[te$mgi_symbol %in% rownames(temp),]
te <- te[isUnique(te$ensembl_gene_id),]
te <- te[isUnique(te$mgi_symbol),]
nasty <- c("ENSMUSG00000021223", "ENSMUSG00000021773", "ENSMUSG00000023845", "ENSMUSG00000025141", "ENSMUSG00000026614", 
           "ENSMUSG00000035620", "ENSMUSG00000036892", "ENSMUSG00000039126", "ENSMUSG00000040554", "ENSMUSG00000047363", 
           "ENSMUSG00000053624", "ENSMUSG00000053868", "ENSMUSG00000055891", "ENSMUSG00000060268", "ENSMUSG00000071686", 
           "ENSMUSG00000073608", "ENSMUSG00000090272" )
te <- te[te$ensembl_gene_id %!in% nasty,]
nasty <- toensid(c("Pcdha9","Vmn1r90","Pcdhgb6"))#,"Zdhhc8","Gm28898","Gm20849","Gm21045","Gm8011","Gm17055",
te <- te[te$ensembl_gene_id %!in% nasty,]
rownames(te) <- te$mgi_symbol
temp <- temp[rownames(temp) %in% te$mgi_symbol,]
rownames(temp) <- te[rownames(temp),]$ensembl_gene_id
write.csv(temp,"out_teichlab/th2crispr_screen_data.csv",row.names = TRUE, quote = FALSE)


######### Sample description #################
samplemeta <- data.frame(sample=colnames(temp))
#samplemeta <- data.frame(sample=sprintf("cr2_%s",colnames(temp)))
samplemeta$target_mgi_symbol<-colnames(temp)
samplemeta$organism<-"mouse"
samplemeta$'Cell Type'<-"Th2"
samplemeta$method<-"CRISPRSCREENRANK"   #technically not all columns
write.csv(samplemeta,"out_teichlab/th2crispr_screen_samplemeta.csv",row.names = FALSE, quote = FALSE)


######################## just DE data for human ###################

##### the score. remove annoying rows
temp <- rbind(tcmouse$de_early, tcmouse$de_late)
temp <- sqldf("select ensembl_gene_id, min(pval) as human_pDE from temp group by ensembl_gene_id")[-1,]
temp$human_pDE[is.na(temp$human_pDE)] <- 1
temp$human_pDE <- round(log10(temp$human_pDE), digits=3)

rownames(temp) <- temp$ensembl_gene_id
temp <- temp[,-1,drop=FALSE]
write.csv(temp,"out_teichlab/th2crispr_humanDE_data.csv",row.names = TRUE, quote = FALSE)


######### Sample description #################
samplemeta <- data.frame(sample=colnames(temp))
samplemeta$organism<-"human"
samplemeta$'Cell Type'<-"Th2"
write.csv(samplemeta,"out_teichlab/th2crispr_humanDE_samplemeta.csv",row.names = FALSE, quote = FALSE)




######################################################################
### sup: DE and TC RNAseq data #######################################
######################################################################

if(TRUE){
  org <- "mouse"
  
  mtpm <- read.csv(sprintf("out_tc/%s/tpm.txt",org),sep="\t",row.names = "gene")
  mtpm <- mtpm[,colnames(mtpm)!="row.names"]
  colnames(mtpm) <- str_replace_all(colnames(mtpm),"rep","")
  colnames(mtpm) <- str_replace_all(colnames(mtpm),"_salmon","")
  
  data.frame(gene)
  dde <- tcmouse$de_early[,c("ensembl_gene_id","pval")]
  dde <- sqldf("select ensembl_gene_id,min(pval) as pval from dde group by ensembl_gene_id")
  dde <- dde[!is.na(dde$ensembl_gene_id),]
  rownames(dde) <- dde$ensembl_gene_id
  mtpm$pval_early <- dde[rownames(mtpm),]$pval
  
  dde <- tcmouse$de_late[,c("ensembl_gene_id","pval")]
  dde <- sqldf("select ensembl_gene_id,min(pval) as pval from dde group by ensembl_gene_id")
  dde <- dde[!is.na(dde$ensembl_gene_id),]
  rownames(dde) <- dde$ensembl_gene_id
  mtpm$pval_late <- dde[rownames(mtpm),]$pval
  
  mtpm$pval_early[is.na(mtpm$pval_early)] <- 1
  mtpm$pval_late [is.na(mtpm$pval_late )] <- 1
  write.csv(mtpm, "out_teichlab/tpm_mouse.csv")  
}

if(TRUE){
  org <- "human"
  
  mtpm <- read.csv(sprintf("out_tc/%s/tpm.txt",org),sep="\t",row.names = "gene")
  mtpm <- mtpm[,colnames(mtpm)!="row.names"]
  colnames(mtpm) <- str_replace_all(colnames(mtpm),"rep","")
  colnames(mtpm) <- str_replace_all(colnames(mtpm),"_salmon","")
  
  data.frame(gene)
  dde <- tchuman$de_early[,c("ensembl_gene_id","pval")]
  dde <- sqldf("select ensembl_gene_id,min(pval) as pval from dde group by ensembl_gene_id")
  dde <- dde[!is.na(dde$ensembl_gene_id),]
  rownames(dde) <- dde$ensembl_gene_id
  mtpm$pval_early <- dde[rownames(mtpm),]$pval
  
  dde <- tchuman$de_late[,c("ensembl_gene_id","pval")]
  dde <- sqldf("select ensembl_gene_id,min(pval) as pval from dde group by ensembl_gene_id")
  dde <- dde[!is.na(dde$ensembl_gene_id),]
  rownames(dde) <- dde$ensembl_gene_id
  mtpm$pval_late <- dde[rownames(mtpm),]$pval
  
  mtpm$pval_early[is.na(mtpm$pval_early)] <- 1
  mtpm$pval_late [is.na(mtpm$pval_late )] <- 1
  write.csv(mtpm, "out_teichlab/tpm_human.csv")  
}

