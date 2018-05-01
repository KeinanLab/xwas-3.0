## This code is used to calculate the gene_based interaction tests
## use the methods described in the following paper:
## Gene-based testing of interactions in association studies of quantitative traits. PLoS Genet.
## Ma L, Clark AG, Keinan A.  9:e1003321 (2013)
## Changed by Yingjie GUO

args <- commandArgs(TRUE)
pvaluesFile = args[1]
corrFile = args[2]
gene1=args[3]
gene2=args[4]
outputFile = args[5]

library(mvtnorm)
library(corpcor)

if (file.info(corrFile)$size == 0) {
  print("Quitting gene_based_inter.R, your correlation file is empty!")
  quit(save = "no", status = 0, runLast = FALSE)
}

#################################################################################################
#################################################################################################
#################################################################################################
# function to get the truncated tail p-value
mytail=function(x) {
  temp=sort(x[x<myth])
  m=length(temp)
  if (m>0){
    ts=0
    for (i in 1:m) ts=ts+1-temp[i]*(length(x)+1)/i
    return(ts/length(x)) }
  else {return(-9999)}
}

# function to obtain the truncated prod p-value
myprod=function(x) {
  temp=sort(x[x<myth])
  m=length(temp)
  if (m>0){
    return(prod(temp))
  }
  else {
    return(1)
  }
}

# function to get the gene-based p-value
gbpvalue = function(reps,numsnps,co,pv) {
  # convert random statistics to p-values.
  # too large to handle whole matrix
  if (numsnps>100){
    obstailpv = mytail(pv)
    simpvtail = 0
    obsprodpv = myprod(pv)
    simpvprod = 0
    numiter = reps/1000
    for (i in 1:numiter){
      currsimpv = pchisq(rmvnorm(1000,mean=rep(0,numsnps),sigma=co)^2,df=1,lower.tail=FALSE)
      simpvtail=simpvtail+sum(apply(currsimpv,1,mytail)>=obstailpv)
      simpvprod=simpvprod+sum(apply(currsimpv,1,myprod)<=obsprodpv)
    }
    list(simpvtail/reps,simpvprod/reps)
  } else {
    pchisq( rmvnorm(reps,mean=rep(0,numsnps),sigma=co)^2,df=1,lower.tail=FALSE ) -> temp
    fractiontail = sum(apply(temp,1,mytail)>=mytail(pv))/reps
    fractionprod = sum(apply(temp,1,myprod)<=myprod(pv))/reps
    list(fractiontail,fractionprod)
  }
}

#################################################################################################
###############################END OF FUNCTIONS##################################################
#################################################################################################

#-----------------------------------------------------------------------------------------------#
## Run the gene_based test
reps=1000 # Maybe offer reps as a parameter? Defaul = 1000

myth=0.05

pv <- as.matrix(read.table(pvaluesFile,header=T))[,7]
pv <- as.numeric(pv)

cori=as.matrix(read.table(corrFile,header=FALSE,sep="\t"))
num.pv <- length(pv)
ptemp <- rep(NA,num.pv)

# calc pmin
minth=qnorm(min(pv)/2)
pmin=1-pmvnorm(lower=minth,upper=-minth,mean=rep(0, num.pv),corr=cori)

me=pmin/min(pv)
ptemp[1]=pmin
keff=1
myo=order(pv)
for (i in 2:(num.pv-1)) {
  keff=keff+sqrt(abs(1-max(abs(cori[myo[1:(i-1)],myo[i]]))^1.71))
  ptemp[i]=pv[myo[i]]*me/keff
}
ptemp[i+1]=pv[myo[i+1]]
pgate=min(ptemp)

# Following VEGAS with adaptive permutations
xx=gbpvalue(reps,num.pv,cori,pv)
if (xx[1] < 0.1 || xx[2] < 0.1 ) {
  reps=10000
  xx=gbpvalue(reps,num.pv,cori,pv)
  if(xx[1]<0.001||xx[2]<0.001){
    reps=1e6
    xx=gbpvalue(reps,num.pv,cori,pv)
  }
}

# write the result
phr = sprintf("\t%s\t%s\t%11.5e\t%11.5e\t%11.5e\t%11.5e", gene1,gene2,pmin,pgate,xx[1],xx[2])
write(phr,file=outputFile,append=TRUE)


