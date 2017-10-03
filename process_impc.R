impc <- read.csv("impc/ALL_genotype_phenotype.csv", stringsAsFactors = FALSE)
head(impc)

impc_terms <- unique(impc$mp_term_name)
impc_terms_imm <- impc_terms[unique(c(
  grep("T cell",impc_terms),
  grep("B cell",impc_terms),
  grep("NK cell",impc_terms),
  grep("Neutrophil",impc_terms),
  grep("Monocyte",impc_terms),
  grep("Eosinophil",impc_terms)))]
impc_terms_immT <- impc_terms[unique(c(
  grep("T cell",impc_terms)
  ))]
impc_terms_immT


#### Particular genes for T cells
impcT <- impc[impc$mp_term_name %in% impc_terms_immT & impc$p_value<1e-2 & impc$p_value!=0,c("mp_term_name","p_value","marker_symbol")]
impcT <- sort(unique(impcT$marker_symbol))
cat(intersect(togenesym2(list_tf), sort(unique(impcT))))
#TFs: Arid1b Elk4 Tcf7 Irf5 Mbd1 Tead3 Xbp1 Bach2 Zfp84 Zfp408 Zfp266
impcT
#Notable: Xbp1   Gimap6    Irf5   Arid1b   Tcf7   

sgenescorer2_matrix[intersect(togenesym2(list_tf), sort(unique(impcT))),]  
#Il13: Mbd1 really high
#Il4: Xbp1, Irf5 borderline. Tcf7
#Gata3: Arid1b decent
#Irf4: Zfp408 is reallyhigh.   Arid1b borderline
#Xbp1: *nothing*


#### Particular genes for any immune related cells
impcImm <- impc[impc$mp_term_name %in% impc_terms_imm & impc$p_value<1e-2 & impc$p_value!=0,c("mp_term_name","p_value","marker_symbol")]
impcImm <- sort(unique(impcImm$marker_symbol))
impcImm
cat(intersect(togenesym2(list_tf), sort(unique(impcImm))))
#Tfs: Arid1b Elk4 Tcf7   Pknox1   Irf5 Mbd1 Tead3 Xbp1 Bach2 Zfp84 Zfp408 Zfp266
sgenescorer2_matrix[intersect(togenesym2(list_tf), sort(unique(impcImm))),]  


#could actually do some kind of GO using this list! for all categories :)



#######################################
## Function: t-test of IMPC terms vs screens
ttestImpcHits <- function(theg){
  pterm <- rep(1, length(impc_terms))
  for(i in 1:length(impc_terms)){
    if(i%%100==0)
      print(i)
    term <- impc_terms[i]
    ppos <- na.omit(impc$p_value[impc$mp_term_name==term & impc$marker_symbol %in% rownames(sgenescorer2_matrix)[sgenescorer2_matrix[,theg]<1000]])
    pneg <- na.omit(impc$p_value[impc$mp_term_name==term & impc$marker_symbol %in% rownames(sgenescorer2_matrix)[sgenescorer2_matrix[,theg]>5000]])
    ppos <- ppos[ppos!=0]  #This and the NA omits are bad hacks
    pneg <- pneg[pneg!=0]  
    if(length(ppos)>3 & length(pneg)>3)
      pterm[i] <- t.test(log10(ppos), log10(pneg), alternative = "less")$p.value
  }
  impc_terms[pterm<1e-1]
}

#Gata3: dec eosinophil cell num      inc sodium level , chlor
#Il4:   plenty. leuk number. Cd4pos Cd25pos ab reg tcell
#Irf4:  not much interesting
#Il13:  not much int
#Xbp1:  likewise
