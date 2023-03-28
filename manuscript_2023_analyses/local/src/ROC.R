#!/usr/bin/env Rscript 

suppressMessages(library("ggplot2"))
suppressMessages(library("pROC"))
suppressMessages(library("optparse"))

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


# 0. Utils ----
# Take the first element of a list and remove it
shift_fn <- function(x) {
  if(length(x) == 0) {return(NA)}
  shiftret <- x[1]
  assign(as.character(substitute(x)), x[2:(length(x))], parent.frame())
  return(shiftret)
}


# 1. Assign arguments to variables ----
arguments <- parse_args(parser, positional_arguments = c(1,Inf))
args <- arguments$args
opt<-arguments$options
# Read stdin 
stdin <- file("stdin")
open(stdin, blocking=TRUE)
z<-read.table(stdin, header=T, sep="\t")
# Take the first argument, the outcomes_col_name (the column name of the outcomes of the classification)
outcomes <- shift_fn(args)
z[,outcomes] <- as.factor(z[,outcomes])
outcomes_levels_ordered <- levels(z[,outcomes])[order(levels(z[,outcomes]))]
controls_level <- outcomes_levels_ordered[1]
cases_level <- outcomes_levels_ordered[2]

message("The first two levels of the \"outcomes\" variable are taken as \"control\" and \"cases\" respectively. The remaining levels are ignored.\nThe specified oredered outcomes are: ")
message(paste0("controls: ",controls_level))
message(paste0("cases: ",cases_level))
# The remaining arguments are the predictors
predictors <- args
# order outcomes column
z <- z[order(z[,outcomes]),]
# Coerce predictor values to 4 significant digits numbers
z[,predictors] <- sapply(z[,predictors],function(x){x <- round(as.numeric(x), digits = 4)})
# convert infinite values if any
z[,predictors][sapply(z[,predictors], is.infinite)] <- 1000


# 2. Define ROC obj ----
roc_formula <- as.formula(paste(collapse = "~", c(outcomes, paste(collapse = "+", predictors))))

if(opt$negative_subsampling > 0){
  z_cases <- z[z[,outcomes]==cases_level,]
  z_controls <- z[z[,outcomes]==controls_level,]
  roc_list <- list()
  for (idx in seq(1:opt$negative_subsampling)) {
    z_sub <- rbind(z_cases, z_controls[sample(1:nrow(z_controls),nrow(z_cases)),])
    roc_list[[idx]] <- suppressMessages(roc(roc_formula, data = z_sub, direction = opt$direction, na.rm = T))
  }
} else {
  roc_list[[1]] <- suppressMessages(roc(roc_formula, data = z, direction = opt$direction, na.rm = T))
}
# 3. ROC curve plot ----
message(paste0("Saving ROC curves pdf plot as ", opt$pdf_name))
# check smoothing opt
if(!is.null(opt$smoothing)){
  if(!opt$smoothing %in%c("binormal", "density", "fitdistr", "logcondens", "logcondens.smooth")){
    stop("Smoothing method not available!")
  }
  roc_list <- lapply(roc_list, function(x){smooth(x, method=opt$smoothing)})
}
# Create label strings for the legend
labels <- predictors
for(i in seq_along(labels)){
  l<-labels[[i]]
  auc<-round(roc_list[[l]]$auc,3) # WARNING: this line do not consider the case in which only one preditor is specified
  labels[[i]]<-paste0(c(l,auc),collapse=": ")
}
# plot
ROC_plot <- ggroc(roc_list) + 
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="darkgrey", linetype="dashed") +
  theme_classic(base_size = 10) +
  theme(plot.margin = margin(1,2,1,2, "cm")
        , plot.title = element_text(hjust = 0.5)
        , legend.title = element_blank()
        , legend.position = c(0.75, 0.22)) +
  scale_color_discrete(labels=labels) +
  coord_fixed()
  #scale_linetype_manual(breaks="custom_t_pot", values="dashed")
# Save plot 
ggsave(paste0(opt$pdf_name,".pdf"), plot=ROC_plot, device="pdf",width = 6,height=6)


# 4. Create stdout dataframe ----
message("Saving ROC AUC comparison table as stdout")
table_out <- data.frame(pred_1=character(),
                        pred_2=character(),
                        AUC_1=double(), 
                        AUC_2=double(),
                        p_value_cfr=double(),
                        stringsAsFactors=FALSE) 
# Fill each row of the dataframe
out_row <- 1
for (roc_objs in roc_list) {
  for(pred_i in names(roc_objs)){
    for(pred_j in names(roc_objs)){
       if(pred_i>pred_j){
  	     roc.pred_1 <- roc_objs[[pred_i]]
  	     roc.pred_2 <- roc_objs[[pred_j]]
  	     AUC_1 <- roc.pred_1$auc
  	     AUC_2 <- roc.pred_2$auc
  	     if(opt$negative_subsampling == 0) {
  	       p_value_cfr <- roc.test(roc.pred_1, roc.pred_2)$p.value
  	     } else {
  	       p_value_cfr <- NA
  	     }
  	     table_out[out_row,] <- c(pred_i,pred_j,AUC_1,AUC_2,p_value_cfr)
  	     out_row <- out_row + 1
      } 
    }
  }
}
# Save dataframe to stdout
write.table(table_out, row.names =FALSE, quote = FALSE, file="", sep="\t") # col.names=NA => first cell empty
      

# sink stderr
#w=warnings()
#sink(stderr())
#if(!is.null(w)){
#  print(w)
#}
#sink()

