############################################################################################
##                                                                                       ###
##             Part of the paper ...                                                     ###
##             Author: Johan Henriksson (mahogny@areta.org)                              ###
##                                                                                       ###
##             This code prepares a total count table for mageck.                        ###
##                                                                                       ###
############################################################################################

#install.packages("bayesplot")
library(bayesplot)

library(DESeq2)
require(MASS)
library("rstan")
rstan_options(auto_write = TRUE)
sprintf("detected cores %s",parallel::detectCores())
fname_stan <- "callhit_SG.stan"

stan_use_cores<-parallel::detectCores()
stan_use_cores<-4



###############################################################################################
######################  STAN data preparation  ################################################
###############################################################################################

#####
## Read the sgrna counts for all the screens
sgrna_count <- read.csv("screen/grnacnt.txt",sep="\t",stringsAsFactors = FALSE)
rownames(sgrna_count) <- sgrna_count[,1]
colnames(sgrna_count)[2] <- "mgi_symbol"
sgrna_count <- sgrna_count[,-c(1)]

##########
## Function: Collect screens and prepare model
collectscreen <- function(grna_pos_screens,grna_neg_screens){
  ##Prepare empty list
  allscreen<-list()
  allscreen$pos_count   <- NULL
  allscreen$neg_count   <- NULL
  allscreen$avgrna      <- NULL
  allscreen$geneforgrna <- NULL
  allscreen$grna    <- NULL 
  allscreen$screenindex <- NULL
  allscreen$disp    <- NULL
  
  allscreen$pos_sfmean  <- NULL
  allscreen$neg_sfmean  <- NULL
  allscreen$pos_screens <- grna_pos_screens
  allscreen$neg_screens <- grna_neg_screens
  
  allscreen$numscreens <- length(grna_pos_screens)
  for(screeni in 1:allscreen$numscreens){
    onescreen_count <- sgrna_count[,c(grna_pos_screens[screeni], grna_neg_screens[screeni])]
    onescreen_geneforgrna <- sgrna_count$mgi_symbol 
    
    forde <- DESeqDataSetFromMatrix(onescreen_count,colData=data.frame(v=factor(c(1,1))),design=~1)
    forde <- estimateSizeFactors(forde)
    forde <- estimateDispersions(forde)

    onescreen_disp <- dispersions(forde)  #one per grna. will need to run deseq2 for each pair of experiments
    
    ##Get size factors and normalize counts
    onescreen_sf <- sizeFactors(forde)  #estimateSizeFactorsForMatrix(grna_count)
    onescreen_ncount <- onescreen_count
    for(i in 1:ncol(onescreen_count)){
      onescreen_ncount[,i] <- onescreen_ncount[,i]/onescreen_sf[i]
    }
    #Calculate average level for each sgRNA in this screen +/-
    onescreen_avgrna <- apply(onescreen_count,1,mean)
    
    
    ## Add to total screen list
    allscreen$pos_count <- c(allscreen$pos_count, onescreen_count[,1])
    allscreen$neg_count <- c(allscreen$neg_count, onescreen_count[,2])
    allscreen$pos_sfmean <- c(allscreen$pos_sfmean, onescreen_sf[1])
    allscreen$neg_sfmean <- c(allscreen$neg_sfmean, onescreen_sf[2])
    allscreen$disp <- c(allscreen$disp, onescreen_disp)
    #  allscreen$pos_disp <- c(allscreen$pos_disp, onescreen_disp)
    #  allscreen$neg_disp <- c(allscreen$neg_disp, onescreen_disp[,2])
    allscreen$avgrna <- c(allscreen$avgrna, onescreen_avgrna)
    allscreen$geneforgrna <- c(allscreen$geneforgrna, onescreen_geneforgrna)
    allscreen$grna <- c(allscreen$grna, rownames(sgrna_count))

    allscreen$screenindex <- c(allscreen$screenindex, rep(screeni, length(onescreen_avgrna)))
  }
  allscreen
}

##########
## Function: Filter out genes with too low counts or bad dispersions
filterscreen <- function(allscreen,minavgrna=300,mingrnapergene=4){  #was 4 grnas
  keep <- 
    !is.na(allscreen$pos_count) & !is.na(allscreen$neg_count) &
    !is.na(allscreen$avgrna) & !is.na(allscreen$disp) &
    allscreen$avgrna>minavgrna
  print(mean(keep))
  keep <- keep & allscreen$geneforgrna %in% names(table(allscreen$geneforgrna[keep])>=mingrnapergene)
  print(mean(keep))

  allscreen$pos_count   <- allscreen$pos_count[keep]
  allscreen$neg_count   <- allscreen$neg_count[keep]
  allscreen$disp        <- allscreen$disp[keep]
  allscreen$avgrna      <- allscreen$avgrna[keep]
  allscreen$geneforgrna <- allscreen$geneforgrna[keep]
  allscreen$grna        <- allscreen$grna[keep]
  allscreen$screenindex <- allscreen$screenindex[keep]

  allscreen  
}


##########
## Function: Turn the collected screens into data ready for stan
screen2stan <- function(allscreen){
  #Get the indices
  allscreen$genes <- unique(allscreen$geneforgrna)
  allscreen$geneindex<-unlist(lapply(allscreen$geneforgrna, function(x) which(allscreen$genes==x)))

  allscreen$grnas <- unique(allscreen$grna)
  allscreen$grnaindex <- unlist(lapply(allscreen$grna, function(x) which(allscreen$grnas==x)))
  
  ## Remaining data
  print(sprintf("remaining genes: %s",length(allscreen$genes)))
  print(sprintf("remaining sgRNA*experiments: %s",length(allscreen$geneindex)))

  list(
    N=allscreen$numscreens,
    M=length(allscreen$pos_count),
    L=length(allscreen$genes),
    P=length(allscreen$grnas),
  
    sd_geneeff=10,
    sd_sf=0.2,
    sd_screeneff=0.05,   #1 seems to give a bit too much freedom?
    sd_grnaeff=2.5,     #Probability of gRNA working
    
    #grnaeff_alpha=0.4,  #For the beta distribution version of grna efficiency
    #grnaeff_beta=0.2,
    
    pos_sfmean  =allscreen$pos_sfmean,  
    neg_sfmean  =allscreen$neg_sfmean,  
    
    pos_count   =allscreen$pos_count,
    neg_count   =allscreen$neg_count,
    
    avgrna      =allscreen$avgrna,
    disp        =allscreen$disp,
    
    grnaindex   =allscreen$grnaindex,
    geneindex   =allscreen$geneindex,
    screenindex =allscreen$screenindex,
    
    #Not used in STAN but for post-processing
    genes=allscreen$genes, #is this ok?
    grnas=allscreen$grnas,
    pos_screens=allscreen$pos_screens,
    neg_screens=allscreen$neg_screens
  )
}




###############################################################################################
######################  MCMC  #################################################################
###############################################################################################

getstaninitfunc <- function(s){
  function(){
    list(
      pos_sf=array(s$pos_sfmean*0), #since v4 this is a correction factor for the SF
      neg_sf=array(s$neg_sfmean*0),  
      screeneff = array(rnorm(s$N-1,mean=0,sd = 0.01)),
      geneeff   = array(rnorm(s$L,  mean=0,sd = 0.01)),
      #grnaeff   = array(runif(s$P, min=0.8,max=0.98))   #beta distribution version
      grnaeff   = array(rnorm(s$P,  mean=0,sd = 0.01))  #previous normal distribution version
    )
  }
}


##########
## Function: Run STAN MCMC
runscreenMCMC <- function(screendata_mod, num_iter=200, num_chains=stan_use_cores){
  options(mc.cores = stan_use_cores) #moved here - makes sense?
  screendata_stan<-stan(
    file=fname_stan,
    data = screendata_mod,
    verbose=TRUE,
    iter=num_iter,
    chains=num_chains,
    init = getstaninitfunc(screendata_mod),
    pars = c("pos_sf","neg_sf","screeneff","geneeff","grnaeff"))
}


##########
## Function: Extract information from an MCMC run
processMCMC <- function(screendata_mod, screendata_stan){
  
  grna.posterior <- as.array(screendata_stan) #expensive
  
  ##### For all the screens
  thena <- rep(NA,screendata_mod$N)
  screeninfo <- data.frame(
    #pos=screendata_mod$pos_screens, #missing info
    #neg=screendata_mod$neg_screens, 
    est_pos_sf=screendata_mod$pos_sfmean, 
    est_neg_sf=screendata_mod$neg_sfmean,
    
    real_pos_sf=thena,
    real_pos_sf.lower=thena,
    real_pos_sf.upper=thena,
    
    real_neg_sf=thena,
    real_neg_sf.lower=thena,
    real_neg_sf.upper=thena,
    
    eff=thena,
    eff.lower=thena,
    eff.upper=thena,
    stringsAsFactors = FALSE)

  gc <- function(x){
    list(
      estimate=mean(x),
      conf.int=quantile(x,c(0.2,0.8))
    )
  }  
  
  for(i in 1:screendata_mod$N){

    # with unfixed size - detect?    
    # tt <- gc(as.double(grna.posterior[,,sprintf("screeneff[%s]",i)]))
    # screeninfo$eff[i] <- exp(tt$estimate)
    # screeninfo$eff.lower[i] <- exp(tt$conf.int[1])
    # screeninfo$eff.upper[i] <- exp(tt$conf.int[2])
    
    
    tt <- gc(as.double(grna.posterior[,,sprintf("pos_sf[%s]",i)]))
    screeninfo$real_pos_sf[i]       <- exp(tt$estimate)   *screendata_mod$pos_sfmean[i]
    screeninfo$real_pos_sf.lower[i] <- exp(tt$conf.int[1])*screendata_mod$pos_sfmean[i]
    screeninfo$real_pos_sf.upper[i] <- exp(tt$conf.int[2])*screendata_mod$pos_sfmean[i]
    
    tt <- gc(as.double(grna.posterior[,,sprintf("neg_sf[%s]",i)]))
    screeninfo$real_neg_sf[i]       <- exp(tt$estimate)   *screendata_mod$neg_sfmean[i]
    screeninfo$real_neg_sf.lower[i] <- exp(tt$conf.int[1])*screendata_mod$neg_sfmean[i]
    screeninfo$real_neg_sf.upper[i] <- exp(tt$conf.int[2])*screendata_mod$neg_sfmean[i]
  }

  if(sprintf("screeneff[%s]",screendata_mod$N) %in% dimnames(grna.posterior)){
    for(i in 1:screendata_mod$N){
      #No fixed screen efficiencies
      tt <- gc(as.double(grna.posterior[,,sprintf("screeneff[%s]",i)]))
      screeninfo$eff[i] <- exp(tt$estimate)
      screeninfo$eff.lower[i] <- exp(tt$conf.int[1])
      screeninfo$eff.upper[i] <- exp(tt$conf.int[2])
    }
  } else {
    #With fixed efficiency of first screen
    screeninfo$eff[1] <- 1
    for(i in 2:screendata_mod$N){
      tt <- gc(as.double(grna.posterior[,,sprintf("screeneff[%s]",i-1)]))
      screeninfo$eff[i] <- exp(tt$estimate)
      screeninfo$eff.lower[i] <- exp(tt$conf.int[1])
      screeninfo$eff.upper[i] <- exp(tt$conf.int[2])
    }
  }

      
  ##### For every gene, extract P-value for each gene effect != 0, and get fold change
  grna.pval <- rep(0,length(screendata_mod$genes))
  grna.eff <- rep(0,length(screendata_mod$genes))
  for(i in 1:length(screendata_mod$genes)){ 
    #should probably not use t-test!
    v<-as.double(grna.posterior[,,sprintf("geneeff[%s]",i)])
    #tt <- t.test(v)
    p <- mean(v<0)
    grna.pval[i] <- min(p,1-p)#  tt$p.value
    grna.eff[i] <- mean(v)#tt$estimate
  }

  geneinfo <- data.frame(
    mgi_symbol=screendata_mod$genes, 
    pval=grna.pval, 
    eff=grna.eff, 
    fc10=log10(exp(2*grna.eff)), # fail!   need to define somehow else? # + mean(screeninfo$eff))), #[screendata_mod$screenindex])),    #log10(exp(2*grna.eff)),  
    stringsAsFactors = FALSE)
  geneinfo <- geneinfo[order(geneinfo$pval),]
  geneinfo$rank <- 1:nrow(geneinfo)
  
  
  
  list(
    gene=geneinfo,
    screen=screeninfo
  )
}


###############################################################################################
######################  MLE   #################################################################
###############################################################################################



##########
## Function: Calculate the MLE and estimate parameters
runscreenMLE <- function(screendata_mod, num_iter=4000){
  ## Run the optimization
  grna.opt <- optimizing(
    stan_model(file=fname_stan),  #topGO/annotationDbi might mess with select() here!
    data = screendata_mod,
    iter=num_iter,
    as_vector=FALSE,
    init = getstaninitfunc(screendata_mod))
  
  ## Calculate information about the genes - only fold change is possible here
  geneinfo <- data.frame(
    mgi_symbol=screendata_mod$genes,#allscreen$genes,rnorm(s$L,  mean=0,sd = 0.01)
    eff=grna.opt$par$geneeff,
    stringsAsFactors = FALSE)
  geneinfo <- geneinfo[order(abs(geneinfo$eff),decreasing = TRUE),]
  geneinfo$rank <- 1:nrow(geneinfo)
  
  ## Calculate information about the screens
  screeninfo <- data.frame(
    pos=screendata_mod$pos_screens, #missing right now 
    neg=screendata_mod$neg_screens, 
    eff=exp(c(0,grna.opt$par$screeneff)), 
    real_pos_sf=screendata_mod$pos_sfmean * exp(grna.opt$par$pos_sf),
    real_neg_sf=screendata_mod$neg_sfmean * exp(grna.opt$par$neg_sf),
    est_pos_sf=screendata_mod$pos_sfmean, 
    est_neg_sf=screendata_mod$neg_sfmean,
    stringsAsFactors = FALSE)

  ## Calculate info on grnas
  grnainfo <- data.frame(
    grna=screendata_mod$grnas,
    p.fine=grna.opt$par$grnaeff,
    stringsAsFactors = FALSE
  )
  
  list(
    gene=geneinfo,
    grna=grnainfo,
    screen=screeninfo)  
}


