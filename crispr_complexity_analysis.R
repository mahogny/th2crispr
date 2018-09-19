
######################################################################
############# Original coverage analysis #############################
######################################################################

#coverage


covperc <- function(take,of) length(unique(sample(1:of,take,replace=TRUE)))/of
covperc(10000,80000)
sizes <- (0:50)*10000


g88 <- sapply(sizes,function(x) covperc(x,88000))
g40 <- sapply(sizes,function(x) covperc(x,40000))
g10 <- sapply(sizes,function(x) covperc(x,10000))
g5  <- sapply(sizes,function(x) covperc(x,5000))

sizes2 <- sizes/1000

png("coverage.png")
plot(sizes2,g88,type="l",col="red",xlab="#cells (k)", ylab="Coverage")
lines(sizes2,g40,type="l",col="blue")
lines(sizes2,g10,type="l",col="green")
lines(sizes2,g5,type="l",col="purple")
lines(sizes2,sizes*0+1,col="gray")
legend(x=400,y=0.8, legend=c("5k","10k","40k","88k"),col=c("purple","green","blue","red"))
dev.off()

### hmm... MOI=1 means on average on infection per cell. how does this translate back to
# prob for one virus into one cell?
#well. p_oneinf = 1/n^2  ???

mean(runif(88000,0,1)>1/(88000^2))








######################################################################
######### Simulation of T cell complexity using the EAT model ########
######################################################################


makelev<-function(d,lev){
  if(length(d)==0)
    return(data.frame(id=c(),stage=c()))
  else
    return(data.frame(id=d,stage=lev))
}

dodiv <- function(cells){
  
  p <- runif(n=nrow(cells),min=0,max=1)
  keep1 <- cells[which(cells$stage==1 & p>=0.24 & p<0.94),]
  div1 <- makelev(cells$id[which(cells$stage==1 & p>=0.94)],2)
  
  p <- runif(n=nrow(cells),min=0,max=1)
  keep2 <- cells[which(cells$stage==2 & p<0.67),]
  div2 <- makelev(cells$id[which(cells$stage==2 & p>=0.67 & p<0.98)],2)
  cont2 <- makelev(cells$id[which(cells$stage==2 & p>=0.98)],3)
  
  p <- runif(n=nrow(cells),min=0,max=1)
  keep3 <- cells[which(cells$stage==2 & p<0.22),]
  div3 <- makelev(cells$id[which(cells$stage==3 & p>=0.22)],3)
  
  cells <- rbind(keep1,div1,div1,
                 keep2,div2,div2,cont2,
                 keep3,div3,div3)
  return(cells)
}


cells <- data.frame(id=round(runif(1e5,min=0,max=88000)),stage=1)

numdiv <- round(5*24/14)
numdiv <- 8

#each 14h
for(i in 1:numdiv){
  cells <- dodiv(cells)
  print(length(which(cells$stage!=1))/nrow(cells))
  print(nrow(cells))
}


n <- as.numeric(table(cells$id))
hist(n)
N=100000
pmorethan2 <- length(which(abs(log2(n[runif(N,min=1,max=length(n))] / n[runif(N,min=1,max=length(n))]))>1))/N
pmorethan2


###### Starting from naive
#30M, p=0.1%   this still means 40k*0.001=40 genes
#10M cells, 88k targets, p more than twice higher, 5%
#5M, p=17%
#1M, p=54%


#hist(log(1+hist(cells$id,plot = FALSE)$counts),breaks=100)


length(unique(cells$id))



######################################################################
######### Chance of getting multiple infected cells ##################
######################################################################

lambdas <- (0:150)/100
prob_0     <- exp(-lambdas)
prob_1     <- lambdas*exp(-lambdas)
prob_2plus <- 1 - prob_0 - prob_1  #Not 0 and not 1
plot(lambdas,prob_1)
plot(lambdas,prob_1/(prob_1+prob_2plus),xlab="MOI",ylab="P[Observed cell is single-transduced]",type="l")



