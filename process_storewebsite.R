#A screen can be thought of as a sample, with FC and mageck p-value
#But further, a gene will have TPM and DE scores

######################################################################
### Store time-course data ###########################################
######################################################################

maketpmmeta <- function(){
  tph<-c("05h",sprintf("%sh",c(1,2,4,6,12,24,48,72)))
  tpd<-c(0.5,1,2,4,6,12,24,48,72)
  
  mtpm_samplemeta <- data.frame()
  for(thtype in c(0,2)){
    for(rep in 1:3){
      tph2 <- c("Naive",sprintf("Th%s_%s_rep%s_salmon",thtype,tph,rep))
      tpd2 <- c(0,tpd)
      mtpm_samplemeta_part <- data.frame(sample=tph2, hours=tpd2, rep=rep, celltype=sprintf("Th%s",thtype),species="mouse")
      if(nrow(mtpm_samplemeta)==0)
        mtpm_samplemeta <- mtpm_samplemeta_part
      else
        mtpm_samplemeta <- rbind(mtpm_samplemeta, mtpm_samplemeta_part)
    }
  }
  mtpm_samplemeta <- mtpm_samplemeta[mtpm_samplemeta$sample %in% colnames(mtpm),]
  mtpm_samplemeta
}


write.csv(maketpmmeta(),"out_teichlab/th2crispr_tc_samplemeta.csv", row.names = FALSE)
write.csv(mtpm,"out_teichlab/th2crispr_tc_tpm.csv", row.names = FALSE)

#Note: sample Naive is used twice! Because it is part of two datasets


######################################################################
### Store screen data ################################################
######################################################################



######### Sample description #################
sgenescorer2_samplemeta <- data.frame(sample=sprintf("cr2_%s",colnames(sgenescorer2_matrix)))
sgenescorer2_samplemeta$target_mgi_symbol<-colnames(sgenescorer2_matrix)
sgenescorer2_samplemeta$organism<-"mouse"
sgenescorer2_samplemeta$'Cell Type'<-"Th2"
sgenescorer2_samplemeta$method<-"CRISPRSCREENRANK"
#sgenescorer2_samplemeta$target_mgi_symbol<-res_samplemeta$sample
write.csv(sgenescorer2_samplemeta,"out_teichlab/th2crispr_screen_samplemeta.csv",row.names = FALSE, quote = FALSE)

##### the score. remove annoying rows
temp <- sgenescorer2_matrix
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
#                   "Gm23994","")) 
rownames(te) <- te$mgi_symbol
temp <- temp[rownames(temp) %in% te$mgi_symbol,]
rownames(temp) <- te[rownames(temp),]$ensembl_gene_id
colnames(temp) <- sprintf("cr2_%s",colnames(temp))
write.csv(temp,"out_teichlab/th2crispr_screen_data.csv",row.names = TRUE, quote = FALSE)

