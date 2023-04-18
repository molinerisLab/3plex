#!/usr/bin/env Rscript 

suppressMessages(library("ggplot2"))
suppressMessages(library("pROC"))
suppressMessages(library("optparse"))

options(error = function() traceback(2))

option_list <- list(
	make_option(c("-O", "--pdfname"), action="store", default="ROC_comparison", dest="pdf_name", help="Specify the name of the pdf image [default: \"ROC_comparison\"]"),
  make_option(c("-d", "--direction"), action="store", default="auto", dest="direction", help="in which direction to make the comparison? “auto” (default): automatically define in which group the median is higher and take the direction accordingly. “>”: if the predictor values for the control group are higher than the values of the case group (controls > t >= cases). “<”: if the predictor values for the control group are lower or equal than the values of the case group (controls < t <= cases). You should set this explicity to “>” or “<” whenever you are resampling or randomizing the data, otherwise the curves will be biased towards higher AUC values. [default \"%default\"]")
)

# Esempio per passare stringa con più elementi da trasformare in un vettore di predittori
# eset2toptable

usage <- "%prog [options] outcomes_col_name [predictor_col1_name, predictor_col2_name, ...] < TABLE
.META: stdin
	tab separated file
.META: stdout
	1	predictor_1
	2	pedictor_2
	3	AUC_1
	4	AUC_2
	5	Pvalue_cfr_1_vs_2
"


parser<-OptionParser(usage = usage, option_list=option_list)


#### Assign arguments to objects
#### Define a function to take the first element of a list and remove it from it
shift_fn <- function(x) {
  if(length(x) == 0) {return(NA)}
  shiftret <- x[1]
  assign(as.character(substitute(x)), x[2:(length(x))], parent.frame())
  return(shiftret)
}

arguments <- parse_args(parser, positional_arguments = c(1,Inf))
args <- arguments$args
opt<-arguments$options

### stdin example: "/sto1/epigen/TPXcCRE/dataset/v8_ChIRP_neg_rand/TERC-cCRE.bed.tpx.raw_-L30-e20-l10-g70-froff.summary.clean.covered_frac.stability.custom_t_pot.neg_pos_rand"
stdin <- file("stdin")
open(stdin, blocking=TRUE)
z<-read.table(stdin,header=T, sep="\t")

### I take the first argument, the outcomes_col_name (the column name of the outcomes of the classification)
outcomes <- shift_fn(args)

### The remaining arguments are the predictors
predictors <- args

# round values
z[,predictors] <- sapply(z[,predictors],function(x){x <- round(as.numeric(x), digits = 4)})
# convert infinite values if any
z[,predictors][sapply(z[,predictors], is.infinite)] <- 1000


table_out<-data.frame(predictor=c(), AUC=c(), WilcoxPval=c())

### For each predictior comute AUC and wilcox.test
for(predictor in predictors){
#	z.auc<-suppressMessages(auc(z[,outcomes], z[,predictor,]))
	roc_obj <- suppressMessages(roc(as.formula(paste(outcomes, "~", predictor)), data = z, direction=opt$direction))
	z.auc <- suppressMessages(auc(roc_obj))
	
	alternative_dir <- "less"
	if(opt$direction == ">") alternative_dir <- "greater"

	z.test<-wilcox.test(data=z, as.formula(paste(predictor,outcomes,sep="~")), alternative=alternative_dir)

	if(z.test$p.value==0){#to avoid printing pvalue==0
		z.test$p.value=.Machine$double.xmin
	}

	write(paste(predictor, z.auc, z.test$p.value, sep="\t"), file="")
}


w=warnings()
sink(stderr())
if(!is.null(w)){
        print(w)
}
sink()

