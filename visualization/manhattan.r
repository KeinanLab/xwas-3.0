library(qqman)
library(optparse)

# Parse command line arguments

option_list = list(
  make_option(c("-f", "--file"), action="store", default=NA, type="character",
              help="REQUIRED: Name of association results file with extension."),
  make_option(c("-c", "--chromosome"), action="store", default=NA, type="integer",
              help="Chromosome number to view. By default, all chromosomes are plotted."),
  make_option(c("-p", "--p-value"), action="store", default='P', type="character",
              help="Name of column containing p-value. Default is P."),
  make_option(c("-o", "--out"), action="store", default='manhattan', type="character",
              help="Name of output PNG. Default is manhattan.")
)

opt = parse_args(OptionParser(option_list=option_list))
if (is.na(opt$f)) {
  stop('Running this script requires an association results file. Run \'Rscript manhattan.r --help\' for more information.')
}

manfile <- paste(opt$o, ".png", sep="")

# Manhattan Plot
gwasResults <- read.table(gwasfile <- opt$f, header=T, stringsAsFactors=FALSE)
gwasResults <- gwasResults[!is.na(gwasResults[,opt$p]),] # filter NA p-vals

if (!is.na(opt$c)) {
  gwasResults <- gwasResults[gwasResults$CHR == opt$c,]
}

png(manfile, h=3, w=7, units="in", res=300)
par(mar=c(5,5,1,1)) # margins
manhattan(gwasResults, chr='CHR', bp='BP', p=opt$p, snp='SNP', suggestiveline = F, genomewideline = F, bty="n")
dev.off()

