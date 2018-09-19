#allchipanno$Batf_th2$Detailed.Annotation   # intron* Inter*

#sum(allchipanno$$Gene.Name=="Cccdc134")

# graph_genes <- unique(allchipanno_gene)

####################
## Final fig 6
graph_leafs <- c("Il4","Il13","Il5","Il10","Ifng","Tnfrsf13b",  "Tgfb1","Il2","Ccdc134","Lrrc40")       

####################
## Act sup fig
#graph_leafs <- c("Il2","Il2ra","Cdk4","Ybx1","Eif5a","Rps27l","Cdca7","Rpl7l1","Plk2","Lif","Cdkn1a")
graph_leafs <- c("Il2","Il2ra","Cdk4","Ybx1","Eif5a","Cdca7","Plk2","Lif")
#,"Cdkn1a","Rpl7l1","Rps27l",


# graph_leafs <- sapply(c("PCNA","CKS2","CDC28","NUSAP1","RRM2","ZWINT","PRC1","TFDP1","CCNA2","CCNB1","MELK","TPX2","BIRC5","NCAPG2",
# "RFWD3","TACC3","CDC2","KIAA0101","MCM2"),normalizesym)






graph_genes <- c("Batf","Bhlhe40","Irf4","Gata3","Stat6","Xbp1","Pparg",
                 graph_leafs)   



###Note: Stat6 -> Lrrc40, falls within an overlapping gene. thus not showing

#graph_genes <- unique(c(allchipanno_gene,graph_leafs))

graph_fromto <- matrix(0,length(graph_genes), length(graph_genes))
colnames(graph_fromto ) <- graph_genes
rownames(graph_fromto ) <- graph_genes
for(i in 1:length(allchipanno)){
  curchip <- allchipanno[[i]]
  
  curchip$len <- curchip$End-curchip$Start + 1
  
  curfrom <- allchipanno_gene[i]
  fromi <- which(graph_genes==curfrom)
  if(length(fromi)!=0){
    print(fromi)
    for(toj in 1:length(graph_genes)){
      curto <- graph_genes[toj]
      
      tssdist <- 20e3
      neartss <- abs(curchip$Distance.to.TSS)<tssdist | (graph_genes[fromi]=="Stat6" & graph_genes[toj]=="Lrrc40")
      
#      graph_fromto[fromi,toj] +
      graph_fromto[fromi,toj] <-  sum(curchip$len[curchip$Nearest.Ensembl==toensid(curto) & neartss])
      # graph_fromto[fromi,toj] <- graph_fromto[fromi,toj] + sum(curchip$Nearest.Ensembl==toensid(curto) & abs(curchip$Distance.to.TSS)<tssdist)
#      graph_fromto[fromi,toj] <- graph_fromto[fromi,toj] + sum(curchip$Peak.Score[curchip$Nearest.Ensembl==toensid(curto) & curchip$Distance.to.TSS<2e-4])
      # if(sum(curchip$Nearest.Ensembl==toensid(curto) & curchip$Distance.to.TSS<2e-4)>0){
      #   graph_fromto[fromi,toj] <- 1
      #   
      # }
    }
  }
}

#allchipanno$Batf_th2

mean(allchipanno$Batf_th2$Peak.Score)
mean(allchipanno$`Stat6-th2`$Peak.Score)


outf <- "/home/mahogny/Dropbox/applyPI/crpaper/followup/ucsc/network2.gv"
cat(file = outf,append = FALSE, "digraph G {\n")
cat(file = outf,append = TRUE, "{ rank = 8; Stat6; }\n")
for(g in graph_leafs){
  cat(file = outf, append = TRUE, sprintf("  { rank = sink; %s; }\n", g))
}
for(g in graph_genes){
  thecol <- "lightgray"   
  if(g %in% c("Bhlhe40","Pparg"))
    thecol <- "lightblue"
  if(g %in% graph_leafs)
    thecol <- "lightyellow"
  
  cat(file = outf, append = TRUE, sprintf("  %s [fillcolor=\"%s\", style=filled, fontsize=20]; \n", g, thecol))
}
for(fromi in 1:length(graph_genes)){
  for(toj in 1:length(graph_genes)){
    if(graph_fromto[fromi,toj]>0){
#      cat(file = outf, append = TRUE, sprintf("  %s -> %s;\n", graph_genes[fromi], graph_genes[toj]))
#      thew <- graph_fromto[fromi,toj]/500
#      thew <- log10(graph_fromto[fromi,toj])
      thew <- sqrt(graph_fromto[fromi,toj])/20
      cat(file = outf, append = TRUE, 
          sprintf("  %s -> %s [penwidth=%s];\n", graph_genes[fromi], graph_genes[toj],thew))
    }
  }
}
cat(file = outf,append = TRUE, "}")


curchip <- allchipanno$`Stat6-th2`
curchip[curchip$Gene.Name=="Lrrc40",]

graph_fromto


#Bhlhe40 and pparg has no width in the peaks???

