#example:
# /broad/tools/bin/R CMD BATCH \
# profile="test.profile"
# method = "GMM"
# outDir = "out"
# verbose = T

options(echo = FALSE)
rm (list=ls())
library ("broadgap.canary")

profileFiles<-commandLineArg("profileFiles", default=NULL, verbose=T)

method<-commandLineArg("method", default="CANARY", verbose=T)
priorsFile<-commandLineArg("priorsFile", default=NULL, verbose=T)
offsetFile<-commandLineArg("offsetFile", default=NULL, verbose=T)
outDir<-commandLineArg("outDir", default=NULL, verbose=T)
verbose<-as.logical (commandLineArg("verbose", default=F, verbose=F))

#params
pseudopointFactor<-as.numeric (commandLineArg("pseudopointFactor", default=NULL, verbose=T))
simplePriorScaling<-as.logical(commandLineArg("simpleScalePriors", default=F, verbose=T))
af_weight<-as.numeric (commandLineArg("af_weight", default=NULL, verbose=T))
hwe_weight<-as.numeric (commandLineArg("hwe_weight", default=NULL, verbose=T))
closeness_weight<-as.numeric (commandLineArg("closeness_weight", default=NULL, verbose=T))
overlap_weight<-as.numeric (commandLineArg("overlap_weight", default=NULL, verbose=T))
globalpenalty_weight<-as.numeric (commandLineArg("globalpenalty_weight", default=NULL, verbose=T))
robust_means<-as.logical(commandLineArg("robust_means", default=NULL, verbose=T))
constrain_variance<-as.numeric(commandLineArg("constrain_variance", default=NULL, verbose=T))
bic_weight<-as.numeric (commandLineArg("bic_weight", default=NULL, verbose=T))

params<-list()
if (length(pseudopointFactor)!=0) params$pseudopointFactor<-pseudopointFactor
if (length(simplePriorScaling)!=0) params$simplePriorScaling<-simplePriorScaling
if (length(af_weight)!=0) params$af_weight<-af_weight
if (length(hwe_weight)!=0) params$hwe_weight<-hwe_weight
if (length(closeness_weight)!=0) params$closeness_weight<-closeness_weight
if (length(overlap_weight)!=0) params$overlap_weight<-overlap_weight
if (length(globalpenalty_weight)!=0) params$globalpenalty_weight<-globalpenalty_weight
if (length(robust_means)!=0) params$robust_means<-robust_means
if (length(constrain_variance)!=0) params$constrain_variance<-constrain_variance
if (length(bic_weight)!=0) params$bic_weight<-bic_weight

print (params)

if (file.exists(profileFiles)==F) cat ("profile file not found", profileFiles, "\n")
if (file.exists(priorsFile)==F) cat ('priors file not found', priorsFile, "\n")

if (length(profileFiles)==0) {
  print ("Must submit at least one profile file for clustering")  
}

#if (is.null(priorsFile)) {
#	r2<- cn_heuristics.callGenotypesByDF(r1$initialAssignments, r1$intensities, outDir=outDir, verbose=verbose)
#} else {

r1<-cnp_clustering.cluster(profileFiles, method=method, priorsFile=priorsFile, offsetFile=offsetFile, outDir=outDir, verbose=verbose, params=params) 

#suppress proc time output
q(runLast=FALSE)
