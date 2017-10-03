source("process_runstan.R")

stan_use_cores<-12

toupload <- Sys.getenv("upload")



doallmcmc <- function(screenname, grna_pos_screens, grna_neg_screens){
  screendata <- collectscreen(grna_pos_screens,grna_neg_screens)
  screendata_filtered <- filterscreen(screendata, minavgrna = 200, mingrnapergene=3)
  screendata_mod <- screen2stan(screendata_filtered)
  screendata_mcmc <- runscreenMCMC(screendata_mod, num_iter = 1000)
  save(screendata_mod, screendata_mcmc, file=sprintf("screen/mcmc_%s.dat",screenname))
}

batchmcmc <- c("gata3","gata3_inv","il4","xbp1","irf4","il13","il4all")

if(toupload==1){
  doallmcmc(
    "gata3",
    c("s9_STG_pos","sc3_gata3_pos"),
    c("s9_STG_neg","sc3_gata3_neg"))
}
if(toupload==2){
  doallmcmc(
    "gata3_inv",
    c("s9_STG_pos","sc3_gata3_neg"),
    c("s9_STG_neg","sc3_gata3_pos"))
}
if(toupload==3){
  doallmcmc(
    "il4",
    c("s8a_STL_pos","sx2_STL_2"),#  #sx2_STL_2 might have been mixed up pos/neg
    c("s8a_STL_neg","sx2_STL_1"))
  
}
if(toupload==4){
  ### Xbp1
  doallmcmc(
    "xbp1",
    c("s8b_XBP1_pos","sc2b_xbp1_pos"),   #s8b fairly similar. sc2b_pos rather extreme outside
    c("s8b_XBP1_neg","sc2b_xbp1_neg"))
  
}
if(toupload==5){
  ### Irf4
  doallmcmc(
    "irf4",
    c("sc2a_irf4_pos","sx4_IRF4_pos2"),    #sc2a_irf4_neg good and clusters with sx4_irf4_neg!!
    c("sc2a_irf4_neg","sx4_IRF4_neg"))     
  
}
if(toupload==6){
  ### Il13
  doallmcmc(
    "il13",  ##this does not seem great
    c("s11_IL13_4","sc1_il13_pos"),   ##_ is low to high?  and pos is high?
    c("s11_IL13_2","sc1_il13_neg"))
}
if(toupload==7){
  doallmcmc(
    "il4all",
    c("s8a_STL_pos","sx2_STL_2","first_il4_pos"),   #  #sx2_STL_2 might have been mixed up pos/neg
    c("s8a_STL_neg","sx2_STL_1","first_il4_neg"))
  
}




#load("mcmc_gata3.dat")
#screendata_mcmc
###### Later on: load one at a time and extract parameters. Save table and download from the cluster

if(FALSE){
  
  
  load("screen/mcmc_gata3.dat")
  
  
  
  
}


