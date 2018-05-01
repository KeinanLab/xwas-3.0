
library(ggplot2)
library(optparse)

suppressPackageStartupMessages(require(optparse))

option_list = list(
  make_option(c("-f", "--file"), action="store", default=NA, type="character",
              help="Name of results file with extension"),
  make_option(c("-u", "--unadjusted"), action="store", default='UNADJ', type="character",
              help="Name of column containing unadjusted p values"),
  make_option(c("-a", "--adjusted"), action="store", default='GC', type="character",
              help="Name of column containing adjusted p values"),
  make_option(c("-l", "--lambda"), action="store", default=NA, type="character",
              help="Genomic inflation factor (printed in XWAS log file)"),
  make_option(c("-n", "--name"), action="store", default='qqplot', type="character",
              help="Root name of output PNG files"),
  make_option(c("-d", "--destination"), action="store", default='./', type="character",
              help="Destination directory of all three PNGs")
)

opt = parse_args(OptionParser(option_list=option_list))

#Import data table that includes pvalues of interest
if (is.na(opt$f)) {
  stop('Running this script requires a p-value file. Run ./qqplot.R --help for more information')
}  

cat('Reading p-value file\n')

source_file = opt$f
data = read.table(source_file, sep="", header=TRUE)

#Create list of pvalues in column titled
column_unadjusted = opt$u
unadjusted = data[[column_unadjusted]]

column_adjusted = opt$a
adjusted = data[[column_adjusted]]

#Create and assign plot arguments f
xlab = expression(Expected~~-log[10](italic(p)))
ylab = expression(Observed~~-log[10](italic(p)))

unadjusted = sort(unadjusted,decreasing = F)
unadjusted = -log10(unadjusted)

adjusted = sort(adjusted, decreasing = F)
adjusted = -log10(adjusted)

expected = 1:length(unadjusted)
expected = expected / length(unadjusted)
expected = -log10(expected)


#making a function to save the plots


cat('Plotting unadjusted data\n')

df_unadjusted = data.frame(expected, unadjusted)

pdf(NULL)

plot_unadjusted = (ggplot(df_unadjusted, aes(x=expected, y=unadjusted))
                   + geom_point(aes(col = "unadjusted"))
                   + scale_color_manual(values = "black")
                   + geom_abline(intercept = 0, slope = 1, color = "red") 
                   + theme_classic()
                   + labs(x = xlab, y = ylab, colour = 'Legend'))

cat(paste('Saving plot of unadjusted data to ', opt$d, opt$n, '_1_unadjusted.png\n', sep = ''))

filename_unadjusted = paste(opt$d, opt$n, '_1_unadjusted.png', sep = "")
ggsave(filename_unadjusted, plot = plot_unadjusted)

#removing the unadjusted plot from memory (less space)
rm(plot_unadjusted)



cat('Plotting adjusted data\n')

df_adjusted = data.frame(expected, adjusted)

plot_adjusted = (ggplot(df_adjusted, aes(x=expected, y=adjusted))
                 + geom_point(aes(col="adjusted"))
                 + scale_color_manual(values = "blue")
                 + geom_abline(intercept = 0, slope = 1, color = "red") 
                 + theme_classic()
                 + labs(x = xlab, y = ylab, colour = 'Legend')
                 + if(!is.na(opt$l)) { 
                   geom_label(
                     aes( 
                       x = 0.5, 
                       y = max(unadjusted),
                       label = paste(expression(lambda),' = ', opt$l)
                     ), 
                     fill = "white")
                 }) 

cat(paste('Saving plot of adjusted data to ', opt$d, opt$n, '_2_adjusted.png
\n', sep = ''))

filename_adjusted = paste(opt$d, opt$n, '_2_adjusted.png', sep = "")
ggsave(filename_adjusted, plot = plot_adjusted)

rm(plot_adjusted)


cat('Plotting combined data for comparison\n')

df_combined = data.frame(expected, unadjusted, adjusted)

plot_compare = (ggplot(df_combined, aes(x=expected, y = value))
                + geom_point(aes(y = unadjusted, col = "unadjusted"))
                + geom_point(aes(y = adjusted, col = "adjusted")) 
                + scale_color_manual(values = c("blue","black"))
                + geom_abline(intercept = 0, slope = 1, color = "red") 
                + labs(x = xlab, y = ylab, colour = 'Legend')
                + theme_classic()
                + if(!is.na(opt$l)) { 
                  geom_label(
                    aes( 
                      x = 0.5, 
                      y = max(unadjusted),
                      label = paste(expression(lambda),' = ', opt$l)
                    ), 
                    fill = "white")
                }) 

cat(paste('Saving plot of combined data to ', opt$d, opt$n, '_3_compare.png\n', sep = ''))

filename_compare = paste(opt$d, opt$n, '_3_compare.png', sep = "")
ggsave(filename_compare, plot = plot_compare)
