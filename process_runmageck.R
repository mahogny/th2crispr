############################################################################################
##                                                                                       ###
##             Part of the paper ...                                                     ###
##             Author: Johan Henriksson (mahogny@areta.org)                              ###
##                                                                                       ###
##             This code prepares a total count table for mageck.                        ###
##             Once run, use the run scripts in the mageck folder                        ###
##                                                                                       ###
############################################################################################


############################################
### Merge cnt20026.txt and cnt19691.txt ####
############################################

#should maybe sequence irf4pos, il13 0/1 deeper

############################################
# read all crispr data and output for mageck

dat <- read.csv("screen/cnt20026.txt",sep="\t")
dat <- cbind(dat[,1:2],  dat[,2+(1:19)]+dat[,2+19+(1:19)])  #sum 2 lanes together
colnames(dat)
ind <- c(10:19, 1:9) #undo wrong order
dat <- dat[,c(1,2,2+order(ind))]

propnames <- c(
  #il13 from low to high
  "s11_IL13_0",
  "s11_IL13_1",
  "s11_IL13_2",
  "s11_IL13_4",       ####3????
  "sx4_IRF4_pos",
  "sx4_IRF4_neg",
  "sx2_STL_0",
  "sx2_STL_1",
  "sx2_STL_2",
  "s8a_STL_pos",
  "s8a_STL_neg",
  "s10_FOXP3_pos",
  "s10_FOXP3_neg",
  "cd103_pos",
  "cd103_neg",
  "s9_STG_pos",
  "s9_STG_neg",
  "s8b_XBP1_pos",
  "s8b_XBP1_neg",
  #### first miseq
  "first_il4_neg",
  "lib",
  "first_il4_pos"
  )

##Add previous sequencing data
dat2 <- read.csv("screen/cnt19691.txt",sep="\t",stringsAsFactors=FALSE)
dat <- cbind(dat,dat2[,-(1:2)])
colnames(dat) <- c("sgRNA","Gene",propnames)

dat_firsttwo <- dat#mergecnt.txt

#write.table(dat, "mergecnt.txt",sep="\t",quote=FALSE, row.names = FALSE)


############################################
### Merge cnt20780.txt #####################
############################################





############################################
# read all crispr data and output for mageck

dat <- read.csv("screen/cnt20780.txt",sep="\t")
dat <- dat[,-18] #random crap

#Forgot this before
dat <- cbind(dat[,1:2],  dat[,2+(1:15)]+dat[,2+15+(1:15)])  #sum 2 lanes together

ind <- c(10:15,   1:9)
dat <- dat[,c(1,2,2+order(ind))]

#il13 from low to high

propnames <- c(
"sc2a_irf4_pos",
"sc2b_xbp1_pos",
"sc2b_xbp1_neg",
"sc1_il13_neg",
"lib_lib2",
"IL4_GFP_pos",
"s11_IL13_low",
"s11_IL13_mid-low",
"lib_origfull",
"lib_clone2",
"sx4_IRF4_pos2",
"sc3_gata3_pos",
"sc3_gata3_neg",
"sc1_il13_pos",
"sc2a_irf4_neg")
colnames(dat) <- c("sgRNA","Gene",propnames)

#TODO: sum sx4_IRF4_pos2 with first count?

dat <- cbind(dat,dat_firsttwo[,-(1:2)])

tc <- apply(dat[,-(1:2)],2,sum)
tc <- tc[order(names(tc))]
round(tc/1e6)


write.table(dat, "screen/grnacnt.txt",sep="\t",quote=FALSE, row.names = FALSE)



############################################
#### Experiment: Use the Mageck 0.5 MLE module, with this hack to combine biological replicates.
#### this is not used in paper which is based on mageck 0.4
############################################
####
#### Note that Magecks default way of combining replicates non-paired does NOT work because different
#### regrowths of the libraries have been used and thus the complexity is different


dat <- read.csv("screen/grnacnt.txt",sep="\t",stringsAsFactors = FALSE)

combinebiorep <- function(sname,mpos,mneg,subsample=TRUE,cutoffmax=0,cutoffavg=10,nullsf=1){ #max was 20
  combinebiorep2(sname,
                 mpos,mneg,subsample,simulateneg=FALSE,cutoffmax=cutoffmax,cutoffavg=cutoffavg,nullsf=nullsf)
  combinebiorep2(sprintf("neg_%s",sname),
                 mpos,mneg,subsample,simulateneg=TRUE, cutoffmax=cutoffmax,cutoffavg=cutoffavg,nullsf=nullsf)
}
combinebiorep2 <- function(sname,mpos,mneg,subsample=TRUE,simulateneg=FALSE,cutoffmax=0,cutoffavg=10,nullsf=1){  
  
  datpos <- dat[,mpos,drop=FALSE]
  datneg <- dat[,mneg,drop=FALSE]
  sizepos <- apply(datpos,2,sum)
  sizeneg <- apply(datneg,2,sum)
  print(sizepos)
  print(sizeneg)
  
  #subsample the libraries so pos & neg are the same level?
  #this wastes statistics but may avoid issuse with poorly sampled sgRNAs
  
  # ###### old wrong way
  # sizemin <- pmin(sizepos,sizeneg)
  # print(sizemin)
  # if(subsample){
  #   for(i in 1:length(mpos)){
  #     datpos[,i] <- floor(datpos[,i]*(sizepos[i]/sizemin[i]))  ## beep, error!
  #     datneg[,i] <- floor(datneg[,i]*(sizeneg[i]/sizemin[i])*nullsf)
  #   }
  # }
  
  msizepos <- min(sizepos)
  msizeneg <- min(sizeneg)   #possibly take the average instead?
  if(subsample){
    for(i in 1:length(mpos)){
      datpos[,i] <- floor(datpos[,i]*(msizepos/sizepos[i]))
      datneg[,i] <- floor(datneg[,i]*(msizeneg/sizeneg[i])*nullsf)
    }
  }
  print("---after resamp---")  
  sizepos <- apply(datpos,2,sum)
  sizeneg <- apply(datneg,2,sum)
  print(sizepos)
  print(sizeneg)
  
    
  newsgrna <- c()
  for(i in 1:length(mpos))
    newsgrna <- c(newsgrna, sprintf("%s_%s",i,dat$sgRNA))

  pos <- c()
  for(i in 1:length(mpos))
    pos <- c(pos, dat[,mpos[i]])

  neg <- c()
  for(i in 1:length(mneg))
    neg <- c(neg, dat[,mneg[i]])

  newdat <- data.frame(
    sgRNA=newsgrna,
    Gene=rep(dat$Gene, length(mpos)),
    pos=pos, neg=neg
  )
  
  if(simulateneg){
    newdat$neg <- rpois(nrow(newdat),newdat$pos)
  }
  
  #Remove low-complexity sgRNAs?
  avcnt <- (newdat$pos + newdat$neg)/2
  maxcnt <- pmax(newdat$pos,newdat$neg)
  keepgrna <- maxcnt>cutoffmax & avcnt>cutoffavg            ########## use the negative control to find a good cut-off!
  print(mean(keepgrna))
  newdat <- newdat[keepgrna,]
  
  ################ size factor normalization screws up this approach!!!
  ################ only works if subsampling ***ALL*** to the same level!!!
  
  #newdat
  write.table(newdat, sprintf("screen/combp_%s.txt",sname),sep="\t",quote=FALSE, row.names = FALSE)
}

### Gata3 
combinebiorep(
  "gata3a",
  c("s9_STG_pos","sc3_gata3_pos"),
  c("s9_STG_neg","sc3_gata3_neg"))
combinebiorep(
  "gata3a_inv",
  c("s9_STG_pos","sc3_gata3_neg"),
  c("s9_STG_neg","sc3_gata3_pos"))

combinebiorep(
  "gata3b",
  c("sc3_gata3_pos"),
  c("sc3_gata3_neg"))

combinebiorep(
  "gata3c",
  c("s9_STG_pos"),
  c("s9_STG_neg"))

combinebiorep(  ###this one is only for checking difference in library size
  "gata3d",
  cutoffmax = 20,
  cutoffavg = 0,
  nullsf = 2,
  c("sc3_gata3_pos"),
  c("sc3_gata3_neg"))


### Il4
combinebiorep("il4a",   #10M reads gone by subsampling
  c("s8a_STL_pos"),#"first_il4_pos"),  #sx2_STL_2 might have been mixed up pos/neg
  c("s8a_STL_neg"))#,"first_il4_neg"))

combinebiorep("il4b",
  c("sx2_STL_2"),#,"first_il4_pos"),  #sx2_STL_2 might have been mixed up pos/neg
  c("sx2_STL_1"))#,"first_il4_neg"))
write.table(newdat, "screen/combp_il4b.txt",sep="\t",quote=FALSE, row.names = FALSE)

combinebiorep("il4c",
  cutoffmax = 3,  #note: the subsampling is really nasty here! should consider eliminating, poor count number
  c("first_il4_pos"),  #sx2_STL_2 might have been mixed up pos/neg
  c("first_il4_neg"))




#sx2_STL_0        sx2_STL_1        sx2_STL_2

### Xbp1
combinebiorep("xbp1a",
  c("s8b_XBP1_pos","sc2b_xbp1_pos"),   #s8b fairly similar. sc2b_pos rather extreme outside
  c("s8b_XBP1_neg","sc2b_xbp1_neg"))

combinebiorep("xbp1b",
  c("sc2b_xbp1_pos"),   #s8b fairly similar. sc2b_pos rather extreme outside
  c("sc2b_xbp1_neg"))

combinebiorep("xbp1c",
  c("s8b_XBP1_pos"),   #s8b fairly similar. sc2b_pos rather extreme outside
  c("s8b_XBP1_neg"))



### Irf4
combinebiorep("irf4a",
  c("sc2a_irf4_pos","sx4_IRF4_pos2"),    #sc2a_irf4_neg good and clusters with sx4_irf4_neg!!
  c("sc2a_irf4_neg","sx4_IRF4_neg"))     
#note: sx4_irf4_pos2 is a resequencing of sx4_irf4_pos at higher depth
combinebiorep("irf4b",
  c("sc2a_irf4_pos"),    #sc2a_irf4_neg good and clusters with sx4_irf4_neg!!
  c("sc2a_irf4_neg"))     
combinebiorep("irf4c",
  c("sx4_IRF4_pos2"),    #sc2a_irf4_neg good and clusters with sx4_irf4_neg!!
  c("sx4_IRF4_neg"))     



# combinebiorep(
#   c("s11_IL13_0","sc1_il13_pos"),   ##_IL13_0 is crap!
#   c("s11_IL13_4","sc1_il13_neg"))
# combinebiorep(
#   c("s11_IL13_2","sc1_il13_pos"),   ##_IL13_0 is crap! replaced with _2
#   c("s11_IL13_4","sc1_il13_neg"))

### Il13
combinebiorep("il13a",  ##this does not seem great
  c("s11_IL13_4","sc1_il13_pos"),   ##_ is low to high?  and pos is high?
  c("s11_IL13_2","sc1_il13_neg"))

combinebiorep("il13b",
  c("sc1_il13_pos"),  
  c("sc1_il13_neg"))

combinebiorep("il13c",
  c("s11_IL13_4"),   ## only 2 and 4 are relevant
  c("s11_IL13_2"))

# 
# combinebiorep(
#   c("sx2_STL_2","s8a_STL_pos","first_il4_pos"),  #sx2_STL_2 might have been mixed up pos/neg
#   c("sx2_STL_1","s8a_STL_neg","first_il4_neg"))
# write.table(newdat, "screen/comb_il4.txt",sep="\t",quote=FALSE, row.names = FALSE)
# 
# #sx2_STL_0 - what is this? and what is 4?
# 
# #############
# 
# combinebiorep(
#   c("s9_STG_pos","sc3_gata3_pos"),
#   c("lib","lib_clone2"))
# write.table(newdat, "screen/comb_s_gata3.txt",sep="\t",quote=FALSE, row.names = FALSE)



combinebiorep("libtest",
              cutoffmax = 0,
              cutoffavg = 0,
              c("lib_origfull"),  
              c("lib_lib2"))




############################################
### Additional QC
############################################

sgrna_count <- read.csv("screen/grnacnt.txt",sep="\t",stringsAsFactors = FALSE)
rownames(sgrna_count) <- sgrna_count[,1]
sgrna_count_gene <- sgrna_count[,2]
sgrna_count <- sgrna_count[,-c(1:2)]
#head(sgrna_count)

#library(DESeq)
library(DESeq2)
mynormalizeGRNA <- function(dat){
  sf <- estimateSizeFactorsForMatrix(dat)
  for(i in 1:ncol(dat)){
    dat[,i] <- dat[,i]/sf[i]
  }
  print(sf)
  dat
}
tc <- apply(sgrna_count,2,sum)
toosmall <- tc<1e6
tc <- tc[order(names(tc))]
round(tc/1e6)

apply(sgrna_sub,2,sum)

#These are not used above right? check TODO. lib is really small. should maybe not use at all
#toosmall <- c("s11_IL13_1","IL4_GFP_pos","s11_IL13_0","s11_IL13_1","sx4_IRF4_pos")
sgrna_count <- sgrna_count[,!toosmall]

sgrna_ncount <- mynormalizeGRNA(sgrna_count)

#Subsampling
sgrna_sub <- sgrna_count
for(i in 1:ncol(sgrna_sub)){
  #if(sum(sgrna_sub[,i])>10e6){
    sgrna_sub[,i] <- floor(sgrna_sub[,i]*10e6/sum(sgrna_sub[,i]))
  #}
}

# #Dirty normalization
# sgrna_ncount <- sgrna_count
# for(i in 1:ncol(sgrna_count)){
#   sgrna_ncount[,i] <- sgrna_ncount[,i]/sum(sgrna_ncount[,i])
# }
hist(log10(as.double(apply(sgrna_ncount,1,sum))))
red <- sgrna_ncount
red <- sgrna_sub
#red <- red[apply(red,1,sum)>5e-5,]
#red <- red[apply(red,1,sum)>1e-6,]
#red <- sgrna_sub
heatmap.2(
  cor(log(1+red)),#,method="spearman"),
  trace = "none",
  margins = c(7,7),
  col=colorRampPalette(c("red", "yellow"))(n = 300)
  )

heatmap.2(
  cor(log(1+red)),
  margins = c(7,7),
  trace = "none"
)

red <- sgrna_ncount[,c(
  grep("sx2_STL",colnames(sgrna_ncount)),
  grep("lib",colnames(sgrna_ncount))
  )]
heatmap.2(
  cor(log(1e-50+red)),
  trace = "none"
)
####### Probably lib_clone2 was never really used. used lib2 instead


red <- sgrna_ncount[,c(
  grep("s8a_STL",colnames(sgrna_ncount)),
  grep("sx2_STL",colnames(sgrna_ncount))
)]
heatmap.2(
  cor(log(1e-50+red)),
  trace = "none",
  margins = c(15,15)
)



##### based on above, ignoring notes:
#fair to believe lib_origfull should replace lib whenever possible.
#or they should be summed together (believe this was a resequencing)

#libclone2 might have been used for sc1_
#sc1_ def seem fine to compare

#first_il4_pos has very low complexity and stands out

#lib_lib2 maybe used for sc3_ and sc2_?

#s8b_xbp1 with lib_origfull - highly confident that this worked

#S11_il13_mid.low clusters with foxp3.  but nothing that makes sense

#s8a_stl_neg stands out but high very high complexity

#s11 is 4 levels from a miseq. from low to high (manifest).
#resequenced and pulled out 2 and 4? s11_IL13_2/4 + s11_IL13_low/mid.low

#sc2a_irf4 pos OR sc2b_xbp1 pos, either is actually negative given clustering?

###############
### guaranteed good

#both irf4 negs cluster together: sc2a_irf4_neg   sx4_irf4_neg
## => one good irf4 all we need

#lib & lib_origfull are the same. the latter a resequencing. do not use lib or add together

#s8b_Xbp**** is good but likely few genes different since very close to library

#s11_il13 2 vs 4  good but likely few genes different


#sc3_gata3** is the best screen with antibody. it is at the outskirt of the main cluster
#but rather different from the other ones. weakly with the library
#sc3_gata3_pos is closest to the library  (hmmm. check by hand!!!)

#s9_stg_pos clusters marginally closer to sc3_gata3_neg 







######### weird comparisons!

#s11_il13_low vs s8a_stl_neg   13M vs 10M
#s10_foxp3_pos vs s11_il13_mid.low    8M vs 15M
#cd103_neg vs sx4_irf4_pos2   16M vs 14M
#s9_STG_neg vs s8b_XBP1_pos  17M vs 21M
#s8b_XBP1_neg vs s11_IL13_4
#sc1_il13 is already fine  -- il13b


combinebiorep2("weirda",
              c("s11_IL13_low"),  
              c("s8a_STL_neg"))
combinebiorep2("weirdb",
               c("s10_FOXP3_pos"),  
               c("s11_IL13_mid.low"))
combinebiorep2("weirdc",
               c("cd103_neg"),  
               c("sx4_IRF4_pos2"))
combinebiorep2("weirdd",  #poor scores. but gata3 comes out first even if compared this way
               c("s9_STG_neg"),     
               c("s8b_XBP1_pos"))
combinebiorep2("weirde",
               c("s8b_XBP1_neg"),  
               c("s11_IL13_4"))

#lib & lib_origfull closest. but lib2 correlates the best with these two



