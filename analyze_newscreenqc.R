


readfilternewscreen <- function(f,cutoff=3){
  d <- read.csv(f,sep="\t",stringsAsFactors = FALSE)
  d <- d[d$sgRNA>cutoff,] #not paired (yet)
  d <- d[d$Gene %in% expressedGenes,]
  colnames(d)[1]<-"mgi_symbol"
  d <- d[order(d$neg.wald.p.value),]
  d$rank <- 1:nrow(d)
  
  d <- smerge(d,ensconvert)
  d <- smerge(d,data.frame(ensembl_gene_id=names(max_mtpm),tpm=max_mtpm),all.x=TRUE)
  d <- d[order(d$neg.wald.p.value),]
  d
}

# ####### survival works better when not paired - at least for gata3!
# surv_irf4 <- readfilternewscreen("screen/new/mle_irf4_surv.gene_summary.txt")
# head(surv_irf4,n=10)
# 
# surv_gata3p <- readfilternewscreen("screen/new/mle_gata3_surv_p.gene_summary.txt")
# head(surv_gata3p,n=25)
# stopgosym(surv_gata3p$Gene[1:1000],bg = surv_gata3p$Gene, nofactor = TRUE)
# 
# surv_gata3 <- readfilternewscreen("screen/new/mle_gata3_surv.gene_summary.txt")
# head(surv_gata3,n=10)
# stopgosym(surv_gata3$Gene[1:2000],bg = surv_gata3$Gene, nofactor = TRUE)
# 
# intersect(surv_gata3$Gene[1:1000], surv_irf4$Gene[1:1000])

##############################################################
##############################################################
##############################################################

checkscreengo <- function(s){
  bg <- unique(sgrna_count_gene) #s$mgi_symbol
  stopgosym(s$mgi_symbol[1:3000],bg = bg, nofactor = TRUE, cutoff = 1e-3)
#  stopgosym(s$mgi_symbol[s$neg.wald.p.value<1e-8],bg = s$mgi_symbol, nofactor = TRUE)
}

####### Gata3
newscreen_gata3a <- readfilternewscreen("screen/new/mlep_gata3a.gene_summary.txt")
head(newscreen_gata3a,n=25)
go_gata3a <- checkscreengo(newscreen_gata3a)
go_gata3a

newscreen_gata3b <- readfilternewscreen("screen/new/mlep_gata3b.gene_summary.txt")
head(newscreen_gata3b,n=25)
go_gata3b <- checkscreengo(newscreen_gata3b)
go_gata3b

newscreen_gata3c <- readfilternewscreen("screen/new/mlep_gata3c.gene_summary.txt")
head(newscreen_gata3c,n=25)
go_gata3c <- checkscreengo(newscreen_gata3c)
go_gata3c


####### Negative control
newscreen_gata3neg <- readfilternewscreen("screen/new/mleneg_gata3a.gene_summary.txt")
head(newscreen_gata3neg,n=25)
go_gata3neg <- stopgosym(newscreen_gata3neg$mgi_symbol[1:2000],bg = newscreen_gata3neg$mgi_symbol, nofactor = TRUE)
go_gata3neg

####### Il4
newscreen_il4a <- readfilternewscreen("screen/new/mlep_il4a.gene_summary.txt")
head(newscreen_il4a,n=25)
go_il4a <- checkscreengo(newscreen_il4a)
go_il4a

newscreen_il4b <- readfilternewscreen("screen/new/mlep_il4b.gene_summary.txt")
head(newscreen_il4b,n=25)
go_il4b <- checkscreengo(newscreen_il4b)
go_il4b

newscreen_il4c <- readfilternewscreen("screen/new/mlep_il4c.gene_summary.txt")
head(newscreen_il4c,n=100)
go_il4c <- checkscreengo(newscreen_il4c)  #Coolest but weirdest GO terms!
go_il4c



####### Il13
newscreen_il13a <- readfilternewscreen("screen/new/mlep_il13a.gene_summary.txt")
head(newscreen_il13a,n=25)
go_il13a <- checkscreengo(newscreen_il13a)
go_il13a

newscreen_il13b <- readfilternewscreen("screen/new/mlep_il13b.gene_summary.txt")
head(newscreen_il13b,n=25)
go_il13b <- checkscreengo(newscreen_il13b)  #suspicious but good genes
go_il13b

newscreen_il13c <- readfilternewscreen("screen/new/mlep_il13c.gene_summary.txt")
head(newscreen_il13c,n=100)
go_il13c <- checkscreengo(newscreen_il13c)
go_il13c



####### Irf4
newscreen_irf4a <- readfilternewscreen("screen/new/mlep_irf4a.gene_summary.txt")
head(newscreen_irf4a,n=25)
go_irf4a <- checkscreengo(newscreen_irf4a)
go_irf4a

newscreen_irf4b <- readfilternewscreen("screen/new/mlep_irf4b.gene_summary.txt")
head(newscreen_irf4b,n=25)
go_irf4b <- checkscreengo(newscreen_irf4b)
go_irf4b

newscreen_irf4c <- readfilternewscreen("screen/new/mlep_irf4c.gene_summary.txt")
head(newscreen_irf4c,n=25)
go_irf4c <- checkscreengo(newscreen_irf4c)
go_irf4c


####### Xbp1
newscreen_xbp1a <- readfilternewscreen("screen/new/mlep_xbp1a.gene_summary.txt")   #il2 and il13 high!
head(newscreen_xbp1a,n=25)
go_xbp1a <- checkscreengo(newscreen_xbp1a)     #regulation of response to stress, signaling TORC1 sig. 
go_xbp1a

newscreen_xbp1b <- readfilternewscreen("screen/new/mlep_xbp1b.gene_summary.txt")   #il2 and il13 high!
head(newscreen_xbp1b,n=25)
go_xbp1b <- checkscreengo(newscreen_xbp1b)     #regulation of response to stress, signaling TORC1 sig. 
go_xbp1b  #### most interesting terms!

newscreen_xbp1c <- readfilternewscreen("screen/new/mlep_xbp1c.gene_summary.txt")   #il2 and il13 high!
head(newscreen_xbp1c,n=25)
go_xbp1c <- checkscreengo(newscreen_xbp1c)     #regulation of response to stress, signaling TORC1 sig. 
go_xbp1c ### also really cool immuno terms


###### library test: 
newscreen_libtest <- readfilternewscreen("screen/new/mlep_libtest.gene_summary.txt")
newscreen_libtest[100:150,]
go_libtest <- checkscreengo(newscreen_libtest)
go_libtest
#Shows up, do not trust (easily): Ccdc134, Nfatc1... that's it!
##quadruple-check now which libraries were sequenced when.
#then combine GO to get some trust
#GO: lymphocyte activation, symbios, these two terms are best avoided


##############################################################
##############################################################
##############################################################



newscreen_weirda <- readfilternewscreen("screen/new/mlep_weirda.gene_summary.txt")   
newscreen_weirdb <- readfilternewscreen("screen/new/mlep_weirdb.gene_summary.txt")   
newscreen_weirdc <- readfilternewscreen("screen/new/mlep_weirdc.gene_summary.txt")   
newscreen_weirdd <- readfilternewscreen("screen/new/mlep_weirdd.gene_summary.txt")   
newscreen_weirde <- readfilternewscreen("screen/new/mlep_weirde.gene_summary.txt")   

head(newscreen_weirdd,n=25)

go_xbp1c <- checkscreengo(newscreen_xbp1c)     
go_xbp1c ### also really cool immuno terms




##############################################################
##############################################################
##############################################################

## all crappy ways!
# #lib+lib2 vs origfull+clone2  --- Il27ra  irf8   Pou6f1  ccdc134  Trappc6a  
# newscreen_null <- readfilternewscreen("screen/new/mle_null.gene_summary.txt")
# head(newscreen_null,n=25)
# go_null <- stopgosym(newscreen_null$Gene[1:2000],bg = newscreen_null$Gene, nofactor = TRUE)
# go_null
# 
# #lib vs lib_origfull  --- Mlxip  Orm3    --- not that this is the resequencing!!!
# newscreen_null2 <- readfilternewscreen("screen/new/mle_null2.gene_summary.txt")
# head(newscreen_null2,n=25)
# go_null2 <- stopgosym(newscreen_null2$Gene[1:5000],bg = newscreen_null2$Gene, nofactor = TRUE)
# go_null2 <- stopgosym(newscreen_null2$Gene[1:2000],bg = newscreen_null2$Gene, nofactor = TRUE)
# go_null2
# 
# #lib_clone2 & lib_origfull  Ccdc134  Ube2m Cd3d
# newscreen_null3 <- readfilternewscreen("screen/new/mle_null3.gene_summary.txt")
# head(newscreen_null3,n=25)
# go_null2 <- stopgosym(newscreen_null2$Gene[1:2000],bg = newscreen_null2$Gene, nofactor = TRUE)
# go_null
# 
# #lib vs lib_clone2  -- Ccdc134 Cd3d
# newscreen_null4 <- readfilternewscreen("screen/new/mle_null4.gene_summary.txt")
# head(newscreen_null4,n=25)
# 
# #lib vs lib_lib2 -- Ccdc134  Irf8 Ube2m
# newscreen_null5 <- readfilternewscreen("screen/new/mle_null5.gene_summary.txt")
# head(newscreen_null5,n=25)



##############################################################
##############################################################
##############################################################


foo <- screengo[screengo$sc3_gata3<1e-5,"sc3_gata3",drop=FALSE]
foo <- foo[-grep("Factor",rownames(foo)),,drop=FALSE]
foo <- foo[order(foo[,1]),,drop=FALSE]
foo
sort(foo)
temp 

#temp$mgi_symbol <- toensid2(temp$Gene)
temp <- temp[temp$Gene %in% expressedGenes,]

head(max_mtpm)


foo <- rownames(sgenescorer_matrix)[order(sgenescorer_matrix$Gata3)]

stopgosym(newscreen_gata3$Gene[1:2000],nofactor = TRUE)
stopgosym(temp$Gene[1:872],bg = temp$Gene, nofactor = TRUE)
stopgosym(foo[1:872],bg = foo, nofactor = TRUE)
nrow(temp)

foo <- rownames(sgenescorer_matrix)[order(sgenescorer_matrix$Gata3)]
foo2 <- rownames(sgenescorer2_matrix)[order(sgenescorer2_matrix$Gata3)]
intersect(temp$Gene[1:2000], foo[1:2000])
intersect(temp$Gene[1:1000], foo2[1:1000]) #better overlap at least

intersect(newscreen_gata3p$Gene[1:1000], foo2[1:1000]) #better overlap at least
intersect(newscreen_gata3$Gene[1:1000], foo2[1:1000]) #paired has better overlap but less cool GO

intersect(newscreen_gata3p$Gene[1:1000], foo[1:1000]) #this method has less overlap than scorer2
intersect(newscreen_gata3$Gene[1:1000], foo[1:1000]) #again better overlap


intersect(intersect(temp$Gene[1:500], foo2[1:500]),togenesym2(list_tf))
intersect(temp$Gene[1:30], foo2[1:30]) 
#Mta3 Cic most important! Mta3 connected to Bcl6 and Nurd and Hif1   .. it is a Gata3-like zinc finger



#########################

intersect(newscreen_gata3$Gene[1:300], newscreen_il4$Gene[1:300])
intersect(newscreen_gata3$Gene[1:300], newscreen_null$Gene[1:300])

go_gata3





##############################################################
##############################################################
##############################################################

testgenesnewmageck <- c("Gata3","Stat6","Tbx21","Rorc","Yy1","Yy2","Cd3d","Il4","Il2","Il4ra","Il13","Nfatc1","Xbp1","Ern1","Atf4")

checktop <- function(s){
  s[s$mgi_symbol %in% testgenesnewmageck,]
}


checktop(newscreen_gata3a)
checktop(newscreen_gata3b) #gata3 & stat6 & il4ra
checktop(newscreen_gata3c) #Yy1(!) stat6 Cd3d.    but il4 comes from here (even if poor p-value)
#I think also this screen is 7 day(?). different genes might come out. il4 would be more important, and the receptor less. makes sense

checktop(newscreen_gata3neg)  #simulation of null-distribution. tbx21 and cd3d come out near rank 1000

checktop(newscreen_il13a) #gata3/il2/tbx21 high. there is hope. il13 crap
checktop(newscreen_il13b) #gata3 and tbx21 strongest here. il13 really crap. cd3d comes from this one
checktop(newscreen_il13c) #il13 itself is rank 1000. otherwise empty
#### Il13b is good but has some crap false positive. il13c is crap but maybe filters

checktop(newscreen_il4a) #il4ra,Gata3 ok
checktop(newscreen_il4b) #terrible except il4 and stat6 come out best rank. but >1000.
checktop(newscreen_il4c) #Rorc really high


checktop(newscreen_xbp1a) #il13,2,yy1,tbx21,rorc really high!
checktop(newscreen_xbp1b) #too many p=0 values, cannot say
checktop(newscreen_xbp1c) #fairly crap. gata3 & atf4 high

checktop(newscreen_irf4a) #only il13 high
checktop(newscreen_irf4b)  #il13,gata3,il4 on top but crappy rank. xbp1 quite high
checktop(newscreen_irf4c) #tbx21 and cd3d near top. and il13. 

head(newscreen_irf4b,n=25)

# 
# checktop(newscreen_irf4)
# checktop(newscreen_il4)
# checktop(newscreen_il13b)
# checktop(newscreen_gata3)
# checktop(newscreen_xbp1)
# 
# checktop(newscreen_null)
# checktop(newscreen_null2)
# checktop(newscreen_null3)
# checktop(newscreen_null4)
# 
# newscreen_null[newscreen_null$Gene %in% testgenesnewmageck,]
# newscreen_null2[newscreen_null2$Gene %in% testgenesnewmageck,]
# newscreen_null3[newscreen_null3$Gene %in% testgenesnewmageck,] #Cannot trust Cd3d
# newscreen_null4[newscreen_null4$Gene %in% testgenesnewmageck,] #Cannot trust Cd3d
# newscreen_null5[newscreen_null5$Gene %in% testgenesnewmageck,]


# newscreen_il13[grep("Cd3",newscreen_il13$Gene),]   
# newscreen_il13b[grep("Cd3",newscreen_il13b$Gene),] #Cd3d/g/(e)  top, does not seem an artefact
# newscreen_null[grep("Cd3",newscreen_null$Gene),]
# #### ISSUE: gprofiler does not use the bg-set. WTF? need to recode this
# newscreen_gata3[newscreen_gata3$Gene %in% testgenesnewmageck,]
# newscreen_gata3p[newscreen_gata3p$Gene %in% testgenesnewmageck,]

##############################################################
##############################################################
##############################################################
checknewdensitypos <- function(temp,xmin=-100,bw=.5){
  
  xmin<- -4
  bw <- .05
  plot(density(log10(temp$neg.p.value),bw=bw),xlim=c(xmin,0))
  lines(density(log10(temp$neg.p.value[temp$mgi_symbol %in% curatedpos$mgi_symbol]),bw=bw),col="blue")
  
  # plot(density(log10(temp$neg.wald.p.value),bw=bw),xlim=c(xmin,0))
  # lines(density(log10(temp$neg.wald.p.value[temp$mgi_symbol %in% curatedpos$mgi_symbol]),bw=bw),col="blue")
  
  print(sort(temp$rank[temp$mgi_symbol %in% kon]))
}
checknewzpos <- function(temp){
  plot(density((temp$neg.z),bw=1))
  lines(density((temp$neg.z[temp$mgi_symbol %in% curatedpos$mgi_symbol]),bw=1),col="blue")
  
  # plot(density((newscreen_gata3$neg.beta),bw=1))
  # lines(density((newscreen_gata3$neg.beta[temp$Gene %in% curatedpos$mgi_symbol]),bw=1),col="blue")
}

checknewdensitypos(newscreen_gata3a)  #perfect!
checknewdensitypos(xmin=-200,newscreen_gata3b)  #this one is good for certain. and it might be best alone comparing to previous hits
checknewdensitypos(xmin=-50,newscreen_gata3c)  #not convinced

plot(log10(newscreen_gata3a$neg.p.value),
     log10(newscreen_gata3a$neg.wald.p.value))



checknewdensitypos(newscreen_gata3neg) #perfect!



###checknewdensitypos(xmin=-320,newscreen_il4)  #what is this?
checknewdensitypos(xmin=-320,newscreen_il4a)  # most similar to current hit selection
checknewdensitypos(xmin=-120,newscreen_il4b) # really different from a and c.  a&b similar read number (b highest). but a & c same library. second most similar to previous
checknewdensitypos(xmin=-320,newscreen_il4c) # this one is really cool
###leaning toward a&b, skip c

checknewdensitypos(bw=0.1,newscreen_il13a)  #
checknewdensitypos(bw=0.1,newscreen_il13b)  # have more confidence here than c.  this one alone better than a. but annoying genes in list
checknewdensitypos(bw=0.1,newscreen_il13c)  # there might be something but weak!
##picked genes from both

checknewdensitypos(xmin=-150,newscreen_xbp1a) #
checknewdensitypos(xmin=-320,newscreen_xbp1b) #
checknewdensitypos(xmin=-40,bw=0.1,newscreen_xbp1c) #
##mainly picked genes from C but they are still in combined

checknewdensitypos(xmin=-200,newscreen_irf4a) #looking really good
checknewdensitypos(xmin=-200,newscreen_irf4b) #looking really good  
checknewdensitypos(xmin=-200,newscreen_irf4c) #looking really good - best of the 3?
#mainly picked from b. still valid after merge


checknewdensitypos(xmin=-200,newscreen_weirda)
checknewdensitypos(xmin=-200,newscreen_weirdb)
checknewdensitypos(xmin=-200,newscreen_weirdc)
checknewdensitypos(xmin=-200,newscreen_weirdd)##hmm
checknewdensitypos(xmin=-200,newscreen_weirde)



checknewzpos(newscreen_gata3a) #toward the right
checknewzpos(newscreen_il13a) #screwed up z
checknewzpos(newscreen_il13b) #toward the right
checknewzpos(newscreen_il13c) #toward the right
checknewzpos(newscreen_il4a)  #toward right, weird
checknewzpos(newscreen_il4b)  #to the right
checknewzpos(newscreen_il4c)  #to the left maybe? looks funky. should maybe not include
# checknewzpos(newscreen_irf4)  #toward the right and looking good
# 
checknewzpos(newscreen_xbp1a)  #messed up. maybe need to turn one xbp1 sample?
checknewzpos(newscreen_xbp1b)  #to the left and weird
checknewzpos(newscreen_xbp1c)  #weakly to the right
# 




compare2screen <- function(s1,s2,n=2000){
  f<-smerge(
    data.frame(mgi_symbol=s1$mgi_symbol[1:n], s1=s1$neg.beta[1:n]),
    data.frame(mgi_symbol=s2$mgi_symbol[1:n], s2=s2$neg.beta[1:n]))
  f<-smerge(
    data.frame(mgi_symbol=s1$mgi_symbol[1:n], s1=s1$neg.z[1:n]),
    data.frame(mgi_symbol=s2$mgi_symbol[1:n], s2=s2$neg.z[1:n]))
  xlim <- symrange(f$s1)
  plot(f$s1,f$s2,xlim=xlim)
}
compare2screen(newscreen_il4a,newscreen_il4b)  #best correlation, but should invert one of them
compare2screen(newscreen_il4a,newscreen_il4c)  #a,b maybe inverted
compare2screen(newscreen_il4b,newscreen_il4c)  #b,c disagree heavily
### maybe keep a and b only. s8a + sx2. drop first4. invert b

newscreen_il4b[newscreen_il4b$rank<1000 & newscreen_il4b$neg.z< -20,]  #Ube2m


compare2screen(newscreen_il13b,newscreen_il13c)  #terrible match. but can't say which to flip
compare2screen(newscreen_irf4b,newscreen_irf4c)  #terrible match. but can't say which to flip
compare2screen(newscreen_gata3b,newscreen_gata3c)  #terrible match. but can't say which to flip

compare2screen(n=8000,newscreen_xbp1b,newscreen_xbp1c)  #terrible match. and b only pulls out negative ones(?)




##############################################################
######### Check hits from simulated negatives ################
##############################################################

checknewdensitypos(readfilternewscreen("screen/new/mlep_neg_gata3b.gene_summary.txt"))
checknewdensitypos(readfilternewscreen("screen/new/mlep_neg_gata3c.gene_summary.txt"))

checknewdensitypos(readfilternewscreen("screen/new/mlep_neg_il4a.gene_summary.txt"))
checknewdensitypos(readfilternewscreen("screen/new/mlep_neg_il4b.gene_summary.txt"))
checknewdensitypos(readfilternewscreen("screen/new/mlep_neg_il4c.gene_summary.txt")) #crash

checknewdensitypos(readfilternewscreen("screen/new/mlep_neg_il13a.gene_summary.txt")) 
checknewdensitypos(readfilternewscreen("screen/new/mlep_neg_il13b.gene_summary.txt")) 
checknewdensitypos(readfilternewscreen("screen/new/mlep_neg_il13c.gene_summary.txt")) 

checknewdensitypos(readfilternewscreen("screen/new/mlep_neg_xbp1a.gene_summary.txt")) 
checknewdensitypos(readfilternewscreen("screen/new/mlep_neg_xbp1b.gene_summary.txt"))   #crash in pos
checknewdensitypos(readfilternewscreen("screen/new/mlep_neg_xbp1c.gene_summary.txt")) 

checknewdensitypos(readfilternewscreen("screen/new/mlep_neg_irf4a.gene_summary.txt")) 
checknewdensitypos(readfilternewscreen("screen/new/mlep_neg_irf4b.gene_summary.txt")) 
checknewdensitypos(readfilternewscreen("screen/new/mlep_neg_irf4c.gene_summary.txt")) 





##############################################################
##############################################################
##############################################################


getgo_listscreens <- function(listscreens) {
  bg <- list()
  fg <- list()
  for(curf in 1:length(listscreens)){
    bg[[curf]] <- unique(sgrna_count_gene)
    fg[[curf]] <- listscreens[[curf]]$mgi_symbol[1:2000]
  }
  names(fg) <- names(listscreens)
#  print(bg)
  getgomatrix(fg,bg)
}
# getgo_listscreens <- function(listscreens) {
#   bg <- unique(sgrna_count_gene) #read in runmageck
#   outgo <- list()
#   for(curf in 1:length(listscreens)){
#     #curf<-1nro
#     print(curf)
#     s<-listscreens[[curf]]
#     x<-stopgosym(s$mgi_symbol[1:2000],bg = bg, nofactor = TRUE, cutoff=1e-2)
#     outgo[[curf]] <- x
#   }
#   names(outgo) <- names(listscreens)
#   
#   ### Combine GO analyses into a matrix instead
#   allcat <- c()
#   for(i in 1:length(outgo))
#     allcat <- union(allcat, outgo[[i]]$term.name)
#   screengo <- matrix(nrow=length(allcat),ncol=length(outgo))
#   rownames(screengo) <- allcat
#   colnames(screengo) <- names(outgo)
#   for(i in 1:length(outgo))
#     screengo[outgo[[i]]$term.name,i] <- outgo[[i]]$p.value
#   screengo[is.na(screengo)] <- 1
#   screengo
# }

newscreengo <- getgo_listscreens(list(
  newscreen_il4a,
  newscreen_il4b#,
  # newscreen_il4c,
  # newscreen_il13a,
  # newscreen_il13b,
  # newscreen_il13c,
  # newscreen_irf4a,
  # newscreen_irf4b,
  # newscreen_irf4c,
  # newscreen_xbp1a,
  # newscreen_xbp1b,
  # newscreen_xbp1c,
  # newscreen_gata3a,
  # newscreen_gata3b,
  # newscreen_gata3c,
  # newscreen_libtest
))
#newscreengo
rownames(newscreengo)
colnames(newscreengo) <- c("il4a","il4b","il4c",
                           "il13a","il13b","il13c",
                           "irf4a","irf4b","irf4c",
                           "xbp1a","xbp1b","xbp1c",
                           "gata3a","gata3b","gata3c",
                           "libtest")
newscreengo_red <- newscreengo
newscreengo_red <- newscreengo_red[apply(newscreengo_red,1,function(x) sort(x)[2])<1,]
newscreengo_red <- newscreengo_red[-grep("intracellular",rownames(newscreengo_red)),]
rownames(newscreengo_red)
#newscreengo_red <- newscreengo_red[-c(1:8,16:19,39,12,),]
heatmap.2(
  trace="none",
  density.info = "none",
  log10(newscreengo_red),
  margins = c(4,15),
  Colv = FALSE,
  col=colorRampPalette(c("red", "yellow","white"))(n = 300)
)
