######################################################################
############# Overall QC #############################################
######################################################################



overlapOEKO <- function(glist,koname){
  #ko* is with a linear model such that negative FC means it went down in KO
  #oe* is such that positive FC means it went up in OE
  kodata <- data.frame(
    ensembl_gene_id=rownames(newkoPadj),
    padj.ko=newkoPadj[,koname],
    fc.ko=newkoFC[,koname],
    fcnorm.ko=newkoFCnorm[,koname]
  )
  glist <- merge(glist, kodata)
  glist <- glist[glist$fc*glist$fcnorm.ko<0,]   #should be opposite?
  glist <- glist[order(glist$padj),]
  glist
}



######################################################################
############# Compare with ChIPseq data ##############################
######################################################################






###############
############### Pparg
###############

glist <- plot_sym_ma("Pparg",FALSE)
glist <- plot_sym_ma("Pparg",FALSE, filterChIP = anno_pparg)
glist

glist[glist$mgi_symbol %in% graph_genes,]  #Il5, Il4, Tnfrsf13b

listTargetSortedByPadj("Slc37a2")
listTargetSortedByPadj("Il5")     ### seems we can trust it


###############
############### lrrc40
###############

glist <- plot_sym_ma("Lrrc40", filteroef = FALSE)
glist

glist[glist$mgi_symbol %in% graph_genes,]  #Il5, Tnfrsf13b, Irf4, Il4, Stat6, Gata3

overlapOEKO(glist, "koLrrc40")
#OE: Lrrc40 up obviously, Igfbp4, Srgn up, Mis18a up, Cd86 up,  ... Cdk4 down, Ccr4 down, Ifngr1 up, Il4 up,   Tcrg*

glist1 <- plot_sym_ma("Lrrc40dec", filteroef = FALSE, dosort = FALSE)
glist2 <- plot_sym_ma("Lrrc40", filteroef = FALSE, dosort = FALSE)
glist2

glist[glist$ensembl_gene_id %in% rownames(newkoPadj)[order(newkoPadj$koLrrc40)][1:1000],]
#il12rb2   Rabgap1l   Lsp1   Fam107b
glist[glist$ensembl_gene_id %in% rownames(newkoPadj)[order(newkoPadj$koLrrc40)][1:100],]
#Galk1  Tpi1  Gatsl3   Cd19   Pou2af1  Il2rg  Ly6c1  



plot(glist1$fc, glist2$fc,cex=0)
text(glist1$fc, glist2$fc, labels = glist1$mgi_symbol)

glist_null <- plot_sym_ma("Foxp3", filteroef = FALSE, dosort = FALSE)
plot(glist1$fc, glist_null$fc,cex=0)
text(glist1$fc, glist_null$fc, labels = glist_null$mgi_symbol)
#Maybe here: Cd34, Thbs1,Tmcc2,Trdv5, Rgs8, Kynu, 

glist_null <- plot_sym_ma("Rorc", filteroef = FALSE, dosort = FALSE)
plot(glist1$fc, glist_null$fc,cex=0)
text(glist1$fc, glist_null$fc, labels = glist_null$mgi_symbol)
#This removes from above: 


plot(glist1$padj, glist2$padj)

glist1[glist1$fc*glist2$fc>0,]



listTargetSortedByPadj("Igfbp4")
listTargetSortedByPadj("Ikzf2")
listTargetSortedByPadj("Mis18a")         #comes up strongest here
listTargetSortedByPadj("Cd86")   #invalid






###############
############### Bhlhe40
###############

thegene <- "oe35_Bhlhe40"  

glist <- plot_sym_ma("Bhlhe40",FALSE)
glist <- plot_sym_ma("Bhlhe40",filterChIP = anno_bhlhe40)
glist

glist[glist$mgi_symbol %in% graph_genes,]  #Tnfrsf13b, Ir4

glist[glist$ensembl_gene_id %in% rownames(newkoPadj)[order(newkoPadj$koBhlhe40)][1:1000],]

overlapOEKO(glist, "koBhlhe40")


###############
############### Scara3
###############


glist <- plot_sym_ma("Scara3",FALSE)
glist







###############
############### Xbp1 (with chipseq from our lab)
###############

glist <- plot_sym_ma("Xbp1",FALSE)
glist <- glist[glist$ensembl_gene_id %in% anno_xbp1$Nearest.Ensembl,]
glist

xbp_fc[xbp_fc$mgi_symbol=="Gata3",]
xbp_fc[xbp_fc$mgi_symbol=="Xbp1",]
xbp_fc[xbp_fc$mgi_symbol=="Il4",]
xbp_fc[xbp_fc$mgi_symbol=="Il13",]
xbp_fc[xbp_fc$mgi_symbol=="Il5",]

sum(xbp_fc$inchip)
