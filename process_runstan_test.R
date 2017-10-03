trialgenes <- c("Gata3","Stat6","Nfatc1","Rora","Il4","Il13","Il2","Tbx21","Batf","Irf4","Xbp1","Fli1","Xbp1","Ern1","Ire1",
                "Etv6","Etv2","Pou6f1","Pou2f2","Yy1","Yy2","Rorc","Foxp3")









################ Run everything
grna_pos_screens <- c("sc3_gata3_pos","s9_STG_pos") #sc3 the best one. and by screeneff, it almost throws away s9_stg right now!
grna_neg_screens <- c("sc3_gata3_neg","s9_STG_neg")

grna_pos_screens <- c("s8a_STL_pos","sx2_STL_2")#  #sx2_STL_2 might have been mixed up pos/neg
grna_neg_screens <- c("s8a_STL_neg","sx2_STL_1")

grna_pos_screens <- c("s8a_STL_pos","sx2_STL_1")# THIS IS NOT HOW IT SHOULD BE??
grna_neg_screens <- c("s8a_STL_neg","sx2_STL_2")

grna_pos_screens <- c("s8b_XBP1_pos","sc2b_xbp1_pos")
grna_neg_screens <- c("s8b_XBP1_neg","sc2b_xbp1_neg")  #Claims effects smaller in second screen

grna_pos_screens <- c("s11_IL13_4","sc1_il13_pos")   #Rora and Fli1 comes out really high
grna_neg_screens <- c("s11_IL13_2","sc1_il13_neg")

grna_pos_screens <- c("sc2a_irf4_pos","sx4_IRF4_pos2")    #sc2a_irf4_neg good and clusters with sx4_irf4_neg!!
grna_pos_screens <- c("sc2a_irf4_neg","sx4_IRF4_neg")     


screendata <- collectscreen(grna_pos_screens,grna_neg_screens)
screendata_filtered <- filterscreen(screendata)
screendata_mod <- screen2stan(screendata_filtered)
screendata_mcmc <- runscreenMCMC(screendata_mod,num_iter = 10)
screendata_stat <- processMCMC(screendata_mod, screendata_mcmc)

screendata_stat$gene$mgi_symbol[1:100]


screendata_posterior <- as.array(screendata_mcmc)
grna.rhat <- rhat(screendata_mcmc)



pairs(screendata_posterior, pars = c("pos_sf[1]", "neg_sf[1]"))
pairs(screendata_posterior, pars = c("screeneff[1]", "screeneff[2]"))



color_scheme_set("red")
#mcmc_intervals(grna.posterior, pars = c("pos_sf[1]","pos_sf[2]"))
mcmc_areas(screendata_posterior, pars = c("pos_sf[1]","pos_sf[2]"))
mcmc_areas(screendata_posterior, pars = c("neg_sf[1]","neg_sf[2]"))
mcmc_areas(screendata_posterior, pars = c("screeneff[1]", "screeneff[2]"))
mcmc_areas(screendata_posterior, pars = c("geneeff[3635]"))


# hist(as.double(grna.posterior[,,"geneeff[1]"]))
# hist(as.double(grna.posterior[,,"geneeff[3635]"]),breaks=100)
# t.test(as.double(grna.posterior[,,"geneeff[3635]"]))$estimate



hist(log10(themle$gene$pval),n=50)




head(screendata_stat$gene,n=200)  
hist(log10(screendata_stat$gene$pval),n=50)

screendata_stat$gene[screendata_stat$gene$mgi_symbol %in% trialgenes,]

stopgosym(df$mgi_symbol[1:1000], df$mgi_symbol, nofactor=TRUE, useTopgo = TRUE)
stopgosym(df$mgi_symbol[1:1000], df$mgi_symbol, nofactor=TRUE,cutoff = 1e-1, useTopgo = FALSE)


t.test(as.double(grna.posterior[,,"geneeff[1]"]))

which(allscreen$genes=="Gata3")




head(as.matrix(grna.stan))

#summary(grna.stan)




####################################
########## MLE #####################
####################################

themle <- runscreenMLE(screendata_mod)
themle$screen
themle$gene[themle$gene$mgi_symbol %in% trialgenes,]

themle$gene[themle$gene$mgi_symbol %in% kon,]

#Gata3 only comes out when the screens are considered of very different weight

head(themle$gene,n=100)


#could the nofactor be a problem??
stopgosym(themle$gene$mgi_symbol[1:1000], themle$gene$mgi_symbol, nofactor=TRUE,cutoff = 1e-1, useTopgo = FALSE)
stopgosym(themle$gene$mgi_symbol[1:500],  themle$gene$mgi_symbol, nofactor=TRUE,cutoff = 1e-1, useTopgo = FALSE)

#topGO finds stuff gprofiler does not!
stopgosym(themle$gene$mgi_symbol[1:2000], themle$gene$mgi_symbol, nofactor=TRUE, useTopgo = TRUE)  
stopgosym(themle$gene$mgi_symbol[1:1000], themle$gene$mgi_symbol, nofactor=TRUE, useTopgo = TRUE)  
stopgosym(themle$gene$mgi_symbol[1:500],  themle$gene$mgi_symbol, nofactor=TRUE, useTopgo = TRUE) 

hist(themle$gene$eff,n=50)






#Might want to feed the average values as an initial guess
# for version 3
# getinitfunc.old <- function(grna.standata){
#   function(){
#     list(
#       pos_sf=grna.standata$pos_sfmean,  
#       neg_sf=grna.standata$neg_sfmean,  
#       screeneff = rnorm(grna.standata$N,mean=1,sd = 0.01),  #  rep(1,grna.standata$N),
#       geneeff   = rnorm(grna.standata$L,mean=0,sd = 0.001)   #rep(0,grna.standata$L)
#     )
#   }
# }

# screeninfo <- data.frame(pos=grna_pos_screens, neg=grna_neg_screens, eff=grna.opt$par$screeneff, 
#                          real_pos_sf=grna.opt$par$pos_sf, est_pos_sf=allscreen$pos_sfmean, 
#                          real_neg_sf=grna.opt$par$neg_sf, est_neg_sf=allscreen$neg_sfmean,
#                          stringsAsFactors = FALSE)

#  var_sf=0.005,
#  var_screeneff=0.1,






#### TODO store 
getgrnaforgene <- function(gene){
  i<-which(screendata_mod$genes==gene)
  keepgrna<-which(screendata_mod$geneindex==i)
  red <- data.frame(
    screeni=screendata_mod$screenindex,
    pos_ncount=screendata_mod$pos_count/screendata_stat$screen$real_pos_sf[screendata_mod$screenindex],
    neg_ncount=screendata_mod$neg_count/screendata_stat$screen$real_neg_sf[screendata_mod$screenindex]
  )
  red <- red[keepgrna,]
  red  
}


plotgrnaforgene <- function(gene){
  red <- getgrnaforgene(gene)
  scol<-c("red","green","blue") #add more...
  for(i in 1:nrow(red)){
    plot(
      c(0,1),
      c(red$pos_ncount[i],red$neg_ncount[i]),type="l",col=scol[red$screeni[i]])

  }
  #TODO
}



###############################################################################################
################# Post-process one MCMC  ######################################################
###############################################################################################


#load("screen/mcmc_gata3.dat")
load("screen.1/mcmc_il4.dat")
#load("screen/mcmc_il13.dat")
load("screen/mcmc_irf4.dat")

trialgenes <- c("Gata3","Stat6","Nfatc1","Rora","Il4","Il13","Il2","Tbx21","Batf","Irf4","Xbp1","Fli1","Xbp1","Ern1","Ire1",
                "Etv6","Etv2","Pou6f1","Pou2f2","Yy1","Yy2","Rorc","Foxp3","Ifng")
source("process_runstan.R")
screendata_stat<-processMCMC(screendata_mod, screendata_mcmc)
t(screendata_stat$screen)

screendata_stat$gene[screendata_stat$gene$mgi_symbol %in% trialgenes,]
head(screendata_stat$gene,n=50)

##TODO: should fix eff[1]=0. or eff[N]=0.  
#makes little difference. but would improve speed

#getgrnaforgene("Gata3")

getmostinteresting <- function(mgi_symbol, nr=2000){
  interesting <- mgi_symbol[(c(
    grep("Il",mgi_symbol[1:nr]),
    grep("Cd",mgi_symbol[1:nr]),
    grep("Gata",mgi_symbol[1:nr]),
    grep("Stat",mgi_symbol[1:nr]),
    grep("Cxc",mgi_symbol[1:nr]),
    grep("Ackr",mgi_symbol[1:nr]),
    grep("Ccr",mgi_symbol[1:nr]),
    grep("Stat",mgi_symbol[1:nr])))]
  #screendata_stat$gene[mgi_symbol %in% interesting,]
  mgi_symbol[mgi_symbol %in% interesting]
}
screendata_stat$gene[screendata_stat$gene$mgi_symbol %in% getmostinteresting(screendata_stat$gene$mgi_symbol),]


source("common_gotest.R")
stopgosym(screendata_stat$gene$mgi_symbol[1:1000],screendata_stat$gene$mgi_symbol, useTopgo = TRUE)
stopgosym(screendata_stat$gene$mgi_symbol[1:2000],screendata_stat$gene$mgi_symbol, useTopgo = TRUE) #works well on stan7, il4



stanc(fname_stan)



############## Potential distributions
hist(rbeta(1000,0.4,0.2),breaks=100)

hist(1/(1+exp(-1.5-4*runif(100000))),breaks = (0:100)/100)
hist(1/(1+exp(-0.5-4*runif(100000))),breaks = (0:100)/100)

hist(atan2(0.5,runif(1000,-10,10)))
hist(atan2(3,rnorm(10000)*40-30)/pi, breaks=50)
hist(atan2(8,200*(rnorm(10000)-0.5))/pi, breaks=50)



###############################################################################################
################# Post-process all MCMCs  #####################################################
###############################################################################################


source("process_runstan.R")

fname_screendir <- "screen"

### Load and summarize all the screens
screenresult <- list()
#batchmcmc <- c("gata3","gata3_inv","il4","xbp1","irf4","il13")   #,"il4all")
batchmcmc <- c("gata3","il4","xbp1","irf4","il13")   #,"il4all")
batchmcmc <- batchmcmc[file.exists(sprintf("%s/mcmc_%s.dat",fname_screendir,batchmcmc))]

for(ones in batchmcmc){
  print(ones)
  load(sprintf("%s/mcmc_%s.dat",fname_screendir,ones))
  screendata_stat<-processMCMC(screendata_mod, screendata_mcmc)
  screenresult[[ones]]<-screendata_stat
}


# require(org.Mm.eg.db)  # generalize
# annotation.db<-org.Mm.eg
# annotation.dbSYMBOL<-org.Mm.egSYMBOL
# annotation.egGO<-org.Mm.egGO

### Make a common matrix of genes
commongenes <- NULL
for(ones in screenresult)
  commongenes <- unique(c(commongenes,ones$gene$mgi_symbol))
newscreen_pval <- matrix(1,nrow=length(commongenes),ncol=length(batchmcmc))
rownames(newscreen_pval) <- commongenes
colnames(newscreen_pval) <- batchmcmc
newscreen_rank <- newscreen_pval*1000000
newscreen_fc  <- newscreen_pval*0
newscreen_inc <- newscreen_pval*0
for(i in 1:length(screenresult)){
  r <- screenresult[[i]]
  print(r)
  newscreen_pval[r$gene$mgi_symbol,i] <- r$gene$pval
  newscreen_rank[r$gene$mgi_symbol,i] <- r$gene$rank
  newscreen_fc  [r$gene$mgi_symbol,i] <- r$gene$fc10
  newscreen_inc [r$gene$mgi_symbol,i] <- 1
}
  
write.table(newscreen_pval,sprintf("%s/out_screen_pval.csv",fname_screendir))
write.table(newscreen_rank,sprintf("%s/out_screen_rank.csv",fname_screendir))
write.table(newscreen_fc,  sprintf("%s/out_screen_fc.csv",fname_screendir))
write.table(newscreen_inc, sprintf("%s/out_screen_inc.csv",fname_screendir))

###########


###############################################################################################
################# Check new screen data  #####################################################
###############################################################################################

fname_screendir <- "screen.0"

newscreen_pval <- read.table(sprintf("%s/out_screen_pval.csv",fname_screendir))
newscreen_rank <- read.table(sprintf("%s/out_screen_rank.csv",fname_screendir))
newscreen_fc   <- read.table(sprintf("%s/out_screen_fc.csv",fname_screendir))
newscreen_inc  <- read.table(sprintf("%s/out_screen_inc.csv",fname_screendir))

#library(WriteXLS)
# WriteXLS(sgenescorer2_matrix, "~/Dropbox/forjhuma/sgenescorer2.xls", row.names = TRUE)
# WriteXLS(newscreen_rank, "~/Dropbox/forjhuma/newscore_rank.xls", row.names = TRUE)
# WriteXLS(newscreen_pval, "~/Dropbox/forjhuma/newscore_pval.xls", row.names = TRUE)
# WriteXLS(newscreen_fc, "~/Dropbox/forjhuma/newscore_fc.xls", row.names = TRUE)

write.csv(sgenescorer2_matrix, "~/Dropbox/forjhuma/sgenescorer2.csv")
write.csv(newscreen_rank, "~/Dropbox/forjhuma/newscore_rank.csv")
write.csv(newscreen_pval, "~/Dropbox/forjhuma/newscore_pval.csv")
write.csv(newscreen_fc, "~/Dropbox/forjhuma/newscore_fc.csv")

############## Common GO analysis with topGO #############
calcscreentopgo <- function(n=1000){
  listfg <- list()
  listbg <- list()
  for(i in 1:ncol(newscreen_pval)){
    listfg[[colnames(newscreen_pval)[i]]] <- rownames(newscreen_rank)[newscreen_rank[,i]<n]
    listbg[[colnames(newscreen_pval)[i]]] <- rownames(newscreen_rank)[newscreen_inc[,i]>0]
    #Since some genes dropped out, should not use these in background!
  }
  getgomatrix(listfg, listbg, listnames=names(listfg), cutoff = 1e-1)
}
#newscreen_go100 <- calcscreentopgo(100)
newscreen_go300 <- calcscreentopgo(300)
newscreen_go500 <- calcscreentopgo(500)
newscreen_go1000 <- calcscreentopgo(1000)

colnames(newscreen_go300) <- c("Gata3","Il4","Xbp1","Irf4","Il13")
pdf("screen/stanwithP_go300.pdf",width = 5)
#png("screen/stanwithP_go300.png",width = 640, height=480)
plotgomatrix(newscreen_go300, keepimmuno = TRUE, cutoff = 0)    #this one looks best
dev.off()
plot.new()
plotgomatrix(newscreen_go300, keepimmuno = FALSE, cutoff = -1.5)    #this one looks best


plotgomatrix(newscreen_go500, keepimmuno = TRUE, cutoff = 0)

plotgomatrix(newscreen_go1000, keepimmuno = FALSE, cutoff = -1)
#plotgomatrix(newscreen_go1000[,c(1,3)], keepimmuno = FALSE, cutoff = -1)


# write.table(newscreen_go1000,sprintf("%s/out_screen_go1000.csv",fname_screendir))
# write.table(newscreen_go500,sprintf("%s/out_screen_go500.csv",fname_screendir))

newscreen_go1000 <- read.table(sprintf("%s/out_screen_go1000.csv",fname_screendir))

newscreen_go1000



# recalc_newscreen_go1000 <- recalcscreenggo(rownames(newscreen_go1000))
# 
# plotgomatrix(
#   as.matrix(recalc_newscreen_go1000)[,-2],keepimmuno = TRUE,cutoff = -0.1)
# plotgomatrix(
#   as.matrix(recalc_newscreen_go1000)[,-2],keepimmuno = FALSE,cutoff = -1.8)
# 


#round(filtergomatriximmuno(newscreen_ggo500[,-2]),digits = 3)


# trialgenes <- c("Gata3","Stat6","Nfatc1","Rora","Il4","Il13","Il2","Tbx21","Batf","Irf4","Xbp1","Fli1","Xbp1","Ern1","Ire1",
#                 "Etv6","Etv2","Pou6f1","Pou2f2","Yy1","Yy2","Rorc","Foxp3","Ifng")
# 
# newscreen_rank[rownames(newscreen_rank) %in% trialgenes,-2]
# newscreen_rank[rownames(newscreen_rank) %in% kon,-2]
# head(screendata_stat$gene,n=50)
# 
# 
# 
# rownames(newscreen_rank)[order(newscreen_rank[,"gata3"])][1:150]
# 
# getmostinteresting(rownames(newscreen_rank)[order(newscreen_rank[,"gata3"])],nr = 100)

gettopscoregenes <- function(newscreen_rank){
  topg <- matrix("",nrow=50,ncol=ncol(newscreen_rank))
  colnames(topg) <- colnames(newscreen_rank)
  for(i in 1:ncol(newscreen_rank)){
    topg[,i] <- rownames(newscreen_rank)[order(newscreen_rank[,i])][1:nrow(topg)]
  }
  topg
}
gettopscoregenes(newscreen_rank[,-2])


#####################################################
### Compare expression levels of top genes for different criteria
plotTPMdensityfor <- function(ensg, col="black",add=FALSE){
  d <- density( log(1+max_mtpm[names(max_mtpm) %in% ensg]), bw=0.5)
  if(add)
    lines(d,col=col)
  else
    plot(d,xlab="Log10 TPM")
}
forgene <- "Gata3"
genen <- 1000
plotTPMdensityfor(toensid(rownames(newscreen_rank)[newscreen_rank[,str_to_lower(forgene)]<genen]))
plotTPMdensityfor(add=TRUE,col="blue",toensid(rownames(sgenescorer_matrix)[sgenescorer_matrix[,forgene]<genen]))
plotTPMdensityfor(add=TRUE,col="red",toensid(rownames(sgenescorer2_matrix)[sgenescorer2_matrix[,forgene]<genen]))


#####################################################
### Compare de_late for top genes for different criteria
plotDEdensityfor <- function(ensg, col="black",add=FALSE){
  v <- log(de_late$pval[de_late$ensembl_gene_id %in% ensg])
  v[is.na(v)]<-1
  d <- density( v ) #, bw=0.5
  if(add)
    lines(d,col=col)
  else
    plot(d,xlab="Log10 TPM")
}
forgene <- "Xbp1"
genen <- 500
plotDEdensityfor(toensid(rownames(newscreen_rank)[newscreen_rank[,str_to_lower(forgene)]<genen]))
plotDEdensityfor(add=TRUE,col="blue",toensid(rownames(sgenescorer_matrix)[sgenescorer_matrix[,forgene]<genen]))
plotDEdensityfor(add=TRUE,col="red",toensid(rownames(sgenescorer2_matrix)[sgenescorer2_matrix[,forgene]<genen]))




toensid( rownames(newscreen_rank)[order(newscreen_rank[,"gata3"])][1:150]  )

toensid2( rownames(newscreen_rank)[order(newscreen_rank[,"gata3"])][1:150] )

#rownames(newscreen_rank)[order(newscreen_rank[,"gata3"])][1:150]
max_mtpm


#####################################################
### Overlap with mageck output

rownames(sgenescorer2_matrix)[sgenescorer2_matrix[,"Il13"]<940]

rownames(newscreen_rank)[newscreen_rank$il13<1000]

vd <- smerge(
  data.frame(mgi_symbol=rownames(newscreen_rank), sr=newscreen_rank$il13<300, stringsAsFactors = FALSE),
  data.frame(mgi_symbol=rownames(sgenescorer2_matrix), r2=sgenescorer2_matrix[,"Il13"]<300, stringsAsFactors = FALSE)
)
vd <- smerge(
  data.frame(mgi_symbol=rownames(newscreen_rank), sr=newscreen_rank$irf4<100, stringsAsFactors = FALSE),
  data.frame(mgi_symbol=rownames(sgenescorer2_matrix), r2=sgenescorer2_matrix[,"Irf4"]<250, stringsAsFactors = FALSE)
)
vd <- smerge(
  data.frame(mgi_symbol=rownames(newscreen_rank), sr=newscreen_rank$gata3<80, stringsAsFactors = FALSE),
  data.frame(mgi_symbol=rownames(sgenescorer2_matrix), r2=sgenescorer2_matrix[,"Gata3"]<200, stringsAsFactors = FALSE)
)
vd <- smerge(
  data.frame(mgi_symbol=rownames(newscreen_rank), sr=newscreen_rank$il4<50, stringsAsFactors = FALSE),
  data.frame(mgi_symbol=rownames(sgenescorer2_matrix), r2=sgenescorer2_matrix[,"Il4"]<80, stringsAsFactors = FALSE)
)
vennDiagram(vc <- vennCounts(vd[,-1]),cex=c(1.5,1.5,1.5))
vd$mgi_symbol[apply(vd[,-1],1,all)]
sum(apply(vd[,-1],1,all))/sum(apply(vd[,-1],1,any))

#Il13: Esrp2, Lag3, Phkb - Esrp2 goes down in mature T cells
#https://www.ncbi.nlm.nih.gov/pubmed/22037216  Esrp2

#Il4: "Kif26a" "Tmem47" "Zfp341" "Zic5"  

#Irf4: Slc5a1  , "Abcg4"    "Chd2"     "Naaladl1" "Otogl"    "Samd7"    "Scn8a"    "Slc25a3"  "Slc5a1"   "Zfp955a" 

#Gata3: Cacng4"   "Spink8"   "Stat6"    "Trappc12"
