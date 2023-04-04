#!/usr/bin/env Rscript 

suppressMessages(library(ggplot2))
suppressMessages(library(pROC))
suppressMessages(library(optparse))

# command line usage and options
options(error = function() traceback(2))
option_list <- list(
  make_option(c("-O", "--pdfname"), action="store", default="ROC_comparison", dest="pdf_name"
              , help="Specify the name of the pdf image [default: \"ROC_comparison\"]"),
  make_option(c("-d", "--direction"), action="store", default="auto", dest="direction"
              , help="Specify in which direction to make the comparison? “auto” (default): automatically define in which group the median is higher and take the direction accordingly. “>”: if the predictor values for the control group are higher than the values of the case group (controls > t >= cases). “<”: if the predictor values for the control group are lower or equal than the values of the case group (controls < t <= cases). You should set this explicity to “>” or “<” whenever you are resampling or randomizing the data, otherwise the curves will be biased towards higher AUC values. [default \"%default\"]"),
  make_option(c("-s", "--smoothing"), action="store", default=NULL, dest="smoothing"
              , help="Method for ROC curve smoothing. Available: binormal, density, fitdistr, logcondens, logcondens.smooth. [default: no smoothing]"),
  make_option(c("--negative_subsampling"), type="integer", default=0, dest="negative_subsampling"
              , help="Number of resampling [default: no subsampling]")
)
usage <- "%prog [options] outcomes_col_name predictor_col_name < TABLE

.META: stdin
	tab separated file

.META: stdout
	1	iteration
	2	sensitivities
	3	specificities
"
parser<-OptionParser(usage = usage, option_list=option_list)

# 1. Assign arguments to variables ----
arguments <- parse_args(parser, positional_arguments = c(1,Inf))
args <- arguments$args
opt<-arguments$options
# Read stdin 
stdin <- file("stdin")
open(stdin, blocking=TRUE)
z<-read.table(stdin, header=T, sep="\t")
# Take the first argument, the outcomes_col_name (the column name of the outcomes of the classification)
outcomes <- args[1]
z[,outcomes] <- as.factor(z[,outcomes])
outcomes_levels_ordered <- levels(z[,outcomes])[order(levels(z[,outcomes]))]
controls_level <- outcomes_levels_ordered[1]
cases_level <- outcomes_levels_ordered[2]

message("The first two levels of the \"outcomes\" variable are taken as \"control\" and \"cases\" respectively. The remaining levels are ignored.\nThe specified oredered outcomes are: ")
message(paste0("controls: ",controls_level))
message(paste0("cases: ",cases_level))
# Second positional argument is the predictor colname
predictors <- args[2]
# order outcomes column
z <- z[order(z[,outcomes]),]
# Coerce predictor values to 4 significant digits numbers
z[,predictors] <- sapply(z[,predictors],function(x){x <- round(as.numeric(x), digits = 4)})
# convert infinite values if any
z[,predictors][sapply(z[,predictors], is.infinite)] <- 1000


# 2. Define ROC obj ----
roc_formula <- as.formula(paste(collapse = "~", c(outcomes, predictors)))

z_cases <- z[z[,outcomes]==cases_level,]
z_controls <- z[z[,outcomes]==controls_level,]

sens_spec_df <- data.frame()
for (idx in seq(1:opt$negative_subsampling)) {
  z_sub <- rbind(z_cases, z_controls[sample(1:nrow(z_controls),nrow(z_cases)),])
  roc_obj <- suppressMessages(roc(roc_formula, data = z_sub, direction = opt$direction, na.rm = T))
  sens_spec_df <- rbind(sens_spec_df, cbind(idx, round(roc_obj$sensitivities,digits = 4), round(roc_obj$specificities,digits = 4)))
}
colnames(sens_spec_df) <- c("iteration","sensitivities","specificities")

# save dataframe to stdout
write.table(sens_spec_df, row.names =FALSE, quote = FALSE, file="", sep="\t")



