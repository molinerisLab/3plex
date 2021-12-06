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
z<-read.table(stdin,header=T)

### I take the first argument, the outcomes_col_name (the column name of the outcomes of the classification)
outcomes <- shift_fn(args)

### The remaining arguments are the predictors
predictors <- args

### Definition of ROC curve formula including all predictors
### ROC object
roc_list <- suppressMessages(roc(as.formula(paste(collapse = "~", c(outcomes, paste(collapse = "+", predictors)))), data = z, direction=opt$direction))

### Create label strings for the legend
labels <- predictors
for(i in seq_along(labels)){
  l<-labels[[i]]
  auc<-round(roc_list[[l]]$auc,3)
  labels[[i]]<-paste0(c(l,auc),collapse=": ")
}

### ROC curve plot
ROC_plot <- ggroc(roc_list) + 
      theme_classic(base_size = 10) +
	    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="darkgrey", linetype="dashed") + 
	    theme(plot.margin = margin(1,2,1,2, "cm")
	          , plot.title = element_text(hjust = 0.5)
	          , legend.title = element_blank()
	          , legend.position = c(0.75, 0.22)) +
	    scale_color_discrete(labels=labels) +
	    scale_linetype_manual(breaks="custom_t_pot",values="dashed") +
      coord_fixed() 
# ggtitle("ROC Comparison") +


### Get pdf file name from --pdf option
### Save plot
ggsave(paste0(opt$pdf_name,".pdf"), plot=ROC_plot, device="pdf",width = 6,height=6)

### Create dataframe for stdout print
table_out <- data.frame(pred_1=character(),
                 pred_2=character(), 
                 AUC_1=double(), 
                 AUC_2=double(),
		 p_value_cfr=double(),	
		 stringsAsFactors=FALSE) 

### Fill each row of the dataframe
count=1
for(i in predictors){
  for(j in predictors){
     if(j>i){
	     roc.pred_1 <- suppressMessages(roc(as.formula(paste(outcomes,"~",i)), data= z, direction=opt$direction))
	     roc.pred_2 <- suppressMessages(roc(as.formula(paste(outcomes,"~",j)), data= z, direction=opt$direction))
	     pred_1 <- i
	     pred_2 <- j
	     AUC_1 <- roc.pred_1$auc
	     AUC_2 <- roc.pred_2$auc
	     p_value_cfr <- roc.test(roc.pred_1, roc.pred_2)$p.value
	     table_out[count,] <- c(pred_1,pred_2,AUC_1,AUC_2,p_value_cfr)
	     count = count + 1
    } 
  }
}

### Save dataframe to stdout
write.table(table_out, row.names =FALSE, quote = FALSE, file="", sep="\t") # col.names=NA => frist cell empty
      


w=warnings()
sink(stderr())
if(!is.null(w)){
        print(w)
}
sink()

