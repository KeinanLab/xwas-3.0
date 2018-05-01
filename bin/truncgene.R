## This code has been adapted from code released as part of the VEGAS software
## (http://www.cell.com/ajhg/abstract/S0002-9297(10)00312-5)
## to carry out the truncated tail strength and truncated product method
## for combining p-values 

args <- commandArgs(TRUE)
gene = args[1]
fname = args[2]
dataname = args[3]

library(corpcor)
library(mvtnorm)

set.seed(42)

###################################################################################################
######################################FUNCTIONS####################################################
###################################################################################################
# function to get the truncated tail p-value
mytail=function(x) {
	lx=length(x)
	temp=sort(x[x<0.05])
	m=length(temp)
	if (m>0){
		ts=0
		for (i in 1:m) ts=ts+1-temp[i]*(lx+1)/i
		return(ts/lx) 
	}
	else {return(0)}
}

# function to obtain the truncated prod p-value
myprod=function(x) {
	temp=sort(x[x<0.05])
	m=length(temp)
	if (m>0){
		return(prod(temp))
	}
	else {return(1)}
}


## function to get the gene-based p-value
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
###################################################################################################
######################################END OF FUNCTIONS#############################################
###################################################################################################

#-------------------------------------------------------------------------------------------------#
## Read in data and rune the gene-based tests. ###
reps=1000
# read in p-values for this gene
pv = as.numeric(as.matrix(read.table("temppv.txt"))[,2])
numsnps = length(pv)

fld = sprintf('%s_ld.ld',dataname)
matrix(scan(fld,quiet=T),nc=numsnps) -> co

#check that co is positive definite. Make diagonals 1.0001.
if(is.positive.definite(co)==F){
	co <- make.positive.definite(co)
	}
if(is.positive.definite(co)==F){
	matrix(scan(fld,quiet=T),nc=numsnps) -> co
	for(i in 1:numsnps){
		co[i,i] <- 1.0001
	}
}
if(is.positive.definite(co)==F){
	for(i in 1:numsnps){
		co[i,i] <- 1.001
	}
}
if(is.positive.definite(co)==F){
	for(i in 1:numsnps){
		co[i,i] <- 1.01
	}
}

# Following VEGAS with adaptive permutations
xx=gbpvalue(reps,numsnps,co,pv) 
if (xx[1] < 0.1 || xx[2] < 0.1 ) {
	reps=10000
	xx=gbpvalue(reps,numsnps,co,pv)
	if (xx[1] < 0.001 || xx[2] < 0.001 ) {
		reps = 1e6
		xx=gbpvalue(reps,numsnps,co,pv)
	}
}

# write the pvalue along with the gene
phr = sprintf("%s %d %f %f",gene,reps,xx[1],xx[2])
write(phr,file=fname,append=TRUE)
