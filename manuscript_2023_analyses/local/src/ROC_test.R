#!/usr/bin/env Rscript 

suppressMessages(library("ggplot2"))
suppressMessages(library("pROC"))
suppressMessages(library("optparse"))
suppressMessages(library("ggpubr"))

options(error = function() traceback(2))

option_list <- list(
        make_option(c("-d", "--direction"), action="store", default="auto", dest="direction", help="in which direction to make the comparison? “auto” (default): automatically define in which group the median is higher and take the direction accordingly. “>”: if the predictor values for the control group are higher than the values of the case group (controls > t >= cases). “<”: if the predictor values for the control group are lower or equal than the values of the case group (controls < t <= cases). You should set this explicity to “>” or “<” whenever you are resampling or randomizing the data, otherwise the curves will be biased towards higher AUC values. [default \"%default\"]"),
	make_option(c("-l", "--levelname"), action="store", default=NULL, dest="level_name", help="Name of the level column used to split the table [deafult:\"None\"]"),
	make_option(c("-O", "--pdfname"), action="store", default="ROC_comparison", dest="pdf_name", help="Name of the pdf image [default: \"ROC_comparison\"]")
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

### If specified, split dataframe according to "level_name" column
level_name<-opt$level_name


if(!is.null(level_name)) {
 
 z_split<-split(z,z[level_name])
 
 roc_obj<-lapply(z_split, function(x){
  		suppressMessages(roc(as.formula(paste(collapse = "~", c(outcomes, paste(collapse = "+", predictors)))), data = x, direction=opt$direction))
	})
	
 table_out <- data.frame(level1=character(), pred1=character(), AUC1=double(), level2=character(), pred2=character(), AUC2=double(), p_value_cfr=double(), stringsAsFactors=FALSE)
	
	### Fill dataframe with AUC comparisons
	count=1
	for (l1 in level_name) {
	  for (p1 in predictors) {
	    for (l2 in level_name) {
	      for (p2 in predictors) {
		if(l1>l2 & p1>p2) {
		  roc1<-roc_obj[[l1]][[p1]]
		  roc2<-roc_obj[[l2]][[p2]]
		  table_out[count,]<-c(l1,p1,roc1$auc,l2,p2,roc2$auc,roc.test(roc1,roc2,method="bootstrap")$p.value)
		  count = count +1
		}
	      }
	    }
	  }
	}


} else {
  roc_obj <- suppressMessages(roc(as.formula(paste(collapse = "~", c(outcomes, paste(collapse = "+", predictors)))), data = z, direction=opt$direction))
  table_out <- data.frame(pred_1=character(), pred_2=character(), AUC_1=double(), AUC_2=double(), p_value_cfr=double(), stringsAsFactors=FALSE)

	### Fill each row of the dataframe
	count=1
	for(i in predictors){
	  for(j in predictors){
	    if(i<j){
	      roc.pred_1 <- roc(as.formula(paste(outcomes,"~",i)), data= z)
	      roc.pred_2 <- roc(as.formula(paste(outcomes,"~",j)), data= z)
	      AUC_1 <- roc.pred_1$auc
	      AUC_2 <- roc.pred_2$auc
	      p_value_cfr <- roc.test(roc.pred_1, roc.pred_2)$p.value
	      table_out[count,] <- c(i,j,AUC_1,AUC_2,p_value_cfr)
	      count = count + 1
	    }
	  }
	}

}

### Create label strings for the legend
#labels <- predictors
#for(i in seq_along(labels)){
#  l<-labels[[i]]
#  auc<-round(roc_obj[[l]]$auc,3)
# labels[[i]]<-paste0(c(l,auc),collapse=": ")
#}

### ROC curve plot
if(length(predictors)>1 && !is.null(level_name) && length(level_name)>1){

	ggroc_list<-lapply(1:length(predictors), function(x){
		  ggroc(lapply(roc_obj, `[[`, x)) +
		    theme_classic(base_size = 10) +
		    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="darkgrey", linetype="dashed") +
		    theme(plot.margin = margin(1,2,1,2, "cm")
			  , plot.title = element_text(hjust = 0.5)
			  , legend.title = element_blank()
			  , legend.position = c(0.75, 0.22)) +
		    scale_linetype_manual(breaks="custom_t_pot",values="dashed") +
		    coord_fixed()
		})

	ROC_plot<-ggarrange(plotlist=ggroc_list)

} else {

	ROC_plot <- ggroc(roc_obj) + 
	      theme_classic(base_size = 10) +
		    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="darkgrey", linetype="dashed") + 
		    theme(plot.margin = margin(1,2,1,2, "cm")
			  , plot.title = element_text(hjust = 0.5)
			  , legend.title = element_blank()
			  , legend.position = c(0.75, 0.22)) +
		    scale_linetype_manual(breaks="custom_t_pot",values="dashed") +
	      coord_fixed() 
	# ggtitle("ROC Comparison") +
	# scale_color_discrete(labels=labels) +
}

### Get pdf file name from --pdf option
### Save plot
ggsave(paste0(opt$pdf_name,".pdf"), plot=ROC_plot, device="pdf",width = 6,height=6)


### Save dataframe to stdout
write.table(table_out, row.names =FALSE, quote = FALSE, file="", sep="\t") # col.names=NA => frist cell empty
      


w=warnings()
sink(stderr())
if(!is.null(w)){
        print(w)
}
sink()

