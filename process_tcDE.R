library(sleuth)

sample_names_early_0=c(paste0("Naive_rep",seq(1,3),"_Th0"), paste0("Th0_05h_rep",seq(1,3)),
                       paste0("Th0_1h_rep",seq(1,3)), paste0("Th0_2h_rep",seq(1,3)),
                       paste0("Th0_4h_rep",seq(1,3)), paste0("Th0_6h_rep",seq(1,3)))
sample_names_early_2=c(paste0("Naive_rep",seq(1,3),"_Th2"), paste0("Th2_05h_rep",seq(1,3)),
                       paste0("Th2_1h_rep",seq(1,3)), paste0("Th2_2h_rep",seq(1,3)),
                       paste0("Th2_4h_rep",seq(1,3)), paste0("Th2_6h_rep",seq(1,3)))

sample_names_late_0=c(paste0("Naive_rep",seq(1,3),"_Th0"), paste0("Th0_6h_rep",seq(1,3)),
                      paste0("Th0_12h_rep",seq(1,3)), paste0("Th0_24h_rep",seq(1,3)),
                      paste0("Th0_48h_rep",seq(1,3)), paste0("Th0_72h_rep",seq(1,3)))
sample_names_late_2=c(paste0("Naive_rep",seq(1,3),"_Th2"), paste0("Th2_6h_rep",seq(1,3)),
                      paste0("Th2_12h_rep",seq(1,3)), paste0("Th2_24h_rep",seq(1,3)),
                      paste0("Th2_48h_rep",seq(1,3)), paste0("Th2_72h_rep",seq(1,3)))

timepoint_names_early_0=rep(c("Naive", "Th0_05h", "Th0_1h", "Th0_2h",
                              "Th0_4h", "Th0_6h"), each=3)
timepoint_names_late_0=rep(c("Naive", "Th0_6h", "Th0_12h", "Th0_24h",
                             "Th0_48h", "Th0_72h"), each=3)


timepoint_names_early_2=rep(c("Naive", "Th2_05h", "Th2_1h", "Th2_2h",
                              "Th2_4h", "Th2_6h"), each=3)
timepoint_names_late_2=rep(c("Naive", "Th2_6h", "Th2_12h", "Th2_24h",
                             "Th2_48h", "Th2_72h"), each=3)

cell_type=c(rep("Th0",18),rep("Th2",18))

timepoints_early=rep(c(0,0.5,1,2,4,6),each=3)
timepoints_late=rep(c(0,6,12,24,48,72),each=3)
names(sf_dirs_Th0_early) = sample_names_early_0
names(sf_dirs_Th0_late) = sample_names_late_0
names(sf_dirs_Th2_early) = sample_names_early_2


names(sf_dirs_Th2_late) = sample_names_late_2

data_summary_early=data.frame("sample"=c(sample_names_early_0, sample_names_early_2),
                              "timepoint"=c(timepoint_names_early_0, timepoint_names_early_2),
                              "time"=rep(timepoints_early,2), "cell"=cell_type,
                              "path"=c(sf_dirs_Th0_early, sf_dirs_Th2_early))
data_summary_late=data.frame("sample"=c(sample_names_late_0, sample_names_late_2),
                             "timepoint"=c(timepoint_names_late_0, timepoint_names_late_2),
                             "time"=rep(timepoints_late,2), "cell"=cell_type,
                             "path"=c(sf_dirs_Th0_late, sf_dirs_Th2_late))
data_summary_early[] <- lapply(data_summary_early, as.character)
data_summary_late[] <- lapply(data_summary_late, as.character)
data_summary_early[,"time"]=as.numeric(data_summary_early[,"time"])
data_summary_late[,"time"]=as.numeric(data_summary_late[,"time"])
data_summary_early$sample=ordered(data_summary_early$sample, levels=unique(data_summary_early$sample))
data_summary_late$sample=ordered(data_summary_late$sample, levels=unique(data_summary_late$sample))
data_summary_early$cell=factor(data_summary_early$cell, levels=unique(data_summary_early$cell))
data_summary_late$cell=factor(data_summary_late$cell, levels=unique(data_summary_late$cell))






## design matrix
### passed to sleuth. using ns() (B-splines basis matrix) to help
### Lower values for the df parameters will capture less “dynamics”,
### but will have more statistical power
### Maybe should also consider placing knots because of uneven tp
time_h_early <- data_summary_early$time
time_h_late <- data_summary_late$time
full_design_early <- model.matrix(formula(~ data_summary_early$cell * ns(time_h_early, df=5)+1))
full_design_late <- model.matrix(formula(~ data_summary_late$cell*ns(time_h_late, df=5)+1))
reduced_design_early <- model.matrix(formula(~ ns(time_h_early, df=5)+1))
reduced_design_late <- model.matrix(formula(~ ns(time_h_late, df=5)+1))
#reduced_design_early <- model.matrix(formula(~ data_summary_early$cell+ns(time_h_early, df=5)+0))
#reduced_design_late <- model.matrix(formula(~ data_summary_late$cell+ns(time_h_late, df=5)+0))



## associate external gene names to transcripts
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "mmusculus_gene_ensembl",
                         host = "www.ensembl.org")
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id",
                                     "ensembl_gene_id",
                                     "external_gene_name"),
                      mart = mart)
t2g <- dplyr::rename(t2g,
                     target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id,
                     ext_gene = external_gene_name)


# Sleuth analysis
## load data to object
so_early_2 <- sleuth_prep(data_summary_early, full_model = full_design_early, target_mapping = t2g)
so_late_2 <- sleuth_prep(data_summary_late, full_model = full_design_late, target_mapping = t2g)

## fit the defined model
so_early_2 <- sleuth_fit(so_early_2)
so_late_2 <- sleuth_fit(so_late_2)

## adding a null model for testing
so_early_2 <- sleuth_fit(so_early_2, formula = reduced_design_early, fit_name = "reduced")
so_late_2 <- sleuth_fit(so_late_2, formula = reduced_design_late, fit_name = "reduced")




so_early_2 <- sleuth_lrt(so_early_2, "reduced", "full")
lrt_results_early_2 <- sleuth_results(so_early_2, 'reduced:full', test_type = 'lrt')
table(lrt_results_early_2[,"qval"] < 0.01)
plot_qq(so_early_2, test = 'reduced:full', test_type = 'lrt', sig_level = 0.01)

so_late_2 <- sleuth_lrt(so_late_2, "reduced", "full")
lrt_results_late_2 <- sleuth_results(so_late_2, 'reduced:full', test_type = 'lrt')
table(lrt_results_late_2[,"qval"] < 0.01)
plot_qq(so_late_2, test = 'reduced:full', test_type = 'lrt', sig_level = 0.01)

??dplyr



# examine top genes
lrt_results_early_2 %>% head(n = 20) %>% dplyr::select(target_id, qval, ens_gene, ext_gene)
early_20=lrt_results_early_2 %>% head(n = 20) %>% dplyr::select(target_id, qval, ens_gene, ext_gene)
early_20=cbind(early_20$target_id, early_20$ext_gene)

lrt_results_late_2 %>% head(n = 20) %>% dplyr::select(target_id, qval, ens_gene, ext_gene)
late_20=lrt_results_late_2 %>% head(n = 20) %>% dplyr::select(target_id, qval, ens_gene, ext_gene)
late_20=cbind(late_20$target_id, late_20$ext_gene)





# Outputting tables
## DE genes
lrt_results_early_2$DE = complete.cases(lrt_results_early_2) & lrt_results_early_2$qval<.01
lrt_results_late_2$DE  = complete.cases(lrt_results_late_2) & lrt_results_late_2$qval<.01

write.table(lrt_results_early_2, "./saved_data/mouse/early_DE_Th0Th2_genes.txt", quote=F, sep="\t", row.names=F)
write.table(lrt_results_late_2,  "./saved_data/mouse/late_DE_Th0Th2_genes.txt", quote=F, sep="\t", row.names=F)




 ##DE human vs mouse

mtpm



nrow(cellcondition)
