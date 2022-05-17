#!/usr/bin/env Rscript
# --default-packages=methods,utils,stats
#options(error = function() traceback(2))
library("optparse")
suppressMessages(library(pROC))

option_list <- list(
	make_option(c("-c", "--coef"), action="store_true", default=FALSE, help="Print coefficients"),
	make_option(c("-u", "--univariate"), action="store_true", default=FALSE, help="For each predictor usePrint coefficients"),
	make_option(c("-U", "--univariate_only"), action="store_true", default=FALSE, help="For each predictor usePrint coefficients"),
	make_option(c("-r", "--residuals"), action="store_true", default=FALSE, help="Print residuals"),
	make_option(c("-s", "--summary"), action="store_true", default=FALSE, help="Print the summary"),
	make_option(c("-l", "--lasso"), action="store_true", default=FALSE, help="Use lasso regularization [default \"%default\"]"),
	make_option(c("-p", "--predict"), action="store", default="", dest="new_data", help="Apply the fitted model to the data in PREDICT filename"),
	make_option(c("-y", "--response_file"), action="store", default="", dest="response_file", help="Response are taken from this file and not in INPUT_FILE, in this case all columns of the INPUT_FILE are used as predictors [default \"%default\"]"),
	make_option(c("-M", "--save_model"), action="store", default="", dest="out_model_file", help="save the model in a .R file [default \"%default\"]"),
	make_option(c("-n", "--col_names"), action="store_true", default=FALSE, help="Print col names in predict section of the output"),
	make_option(c("-b", "--bootstrap"), type="integer", default=0, help="Perform a bootstrap validation with this number of repetitions [default \"%default\"]"),
	make_option(c("-B", "--balancing"), action="store_true", default=FALSE, help="Perform a balanced bootstrap, meaning that in each run the numer of chose positive and negative are the same"),
	make_option(c("-H", "--header"), action="store_false", default=TRUE, help="The input files do not have an header in the first row"),
	make_option(c("-i", "--linear"), action="store_true", default=FALSE, help="fit linear models instead of logistic ones"),
	make_option(c("-S", "--stability"), action="store_true", default=FALSE, help="Performs lasso and stability selection, but returns only the maximum of selection probabilities for each entry. If the lasso coefficients are needed use option -l. Reads .gz files.")
)
parser<-OptionParser(usage = "%prog [options] INPUT_FILE FORMULA
INPUT_FILE must be tab delimited with a single header row to name the columns
FORMULA is a classical R formula

IUNPT_FILE can be gzipped.

OUTPUT:
.META: -s
	1	DF_null		1999
	2	DF_residual	1996
	3	AIC		1844.01419152863
	4	AUC		0.862790977998057
	5	WILCOX_PVALUE	correspond to the AUC statistics, do not take into account the DF
	5	ANOVA_PVALUE	compare the model whit a model with only the intercept anova(model, update(model, ~1)), take into account the number of predictors and DF

.META: -c
	1	label
	2	coeff
	3	std_err
	4	t_value
	5	p_value


.META: -u (logistic model only)
	1	label
	2	coeff
	3	std_err
	4	t_value
	5	p_value
	6	AUC for the univariate model


.META: -r
	1	row_num
	2	residual

.META: -p
	1	predicted_probability	[last_row]

.META: -b 100
	1 param_label	coefficent or AUC or Pvalue
	2 number_of_repetitions
	3 estimated_value
	4 lower_bound 95% confidence interval
	5 upper_bound 95% confidence interval

if multiple option are given, separate each portion of the output with fasta header

HOW TO INTERPRET THE COEFFICIENTS?

The dependent variable factors are sorted in numeric or lexicografic order, depending on wether the values of the factors are integer or strings.
The first factor is taken as negative (if the factors are [0,1] the tool behave as expected).
The coefficient of a independent variable is positive (negative) if values of the variable are bigger in the case of positive (negative) factor.


The logistic regression coefficients give the change in the log odds of the outcome for a one unit increase in the predictor variable.
When a binary outcome variable is modeled using logistic regression, it is assumed that the logit transformation of the outcome variable has a linear relationship with the predictor variables.
Let's say that the probability of success of some event is .8.  Then the probability of failure is 1- .8 = .2.  The odds of success are defined as the ratio of the probability of success over the probability of failure.  In our example, the odds of success are .8/.2 = 4.  That is to say that the odds of success are  4 to 1.  If the probability of success is .5, i.e., 50-50 percent chance, then the odds of success is 1 to 1. 

The transformation from probability to odds is a monotonic transformation, meaning the odds increase as the probability increases or vice versa.  
The transformation from odds to log of odds is the log transformation.
http://www.ats.ucla.edu/stat/r/dae/logit.htm
http://www.ats.ucla.edu/stat/mult_pkg/faq/general/odds_ratio.htm

WARING: if a variable is exactly colinear the relative coefficient row is not printed.

ROC/AUC: we use an explicit direction when calculating AUCs: we want predictor values for the control group lower or equal than 
the values of the case group (controls < t <= cases)

BOOTSTRAP

The reported 95% confidence intervals are computed with \"perc\" method, while \"norm\" assumes an asymptotic normal distribution
and use the bootstrap estimates of bias and variance for the parameters of the distribution, the \"perc\" has less restrictive
assumptions.
[http://cran.r-project.org/doc/Rnews/Rnews_2002-3.pdf also availabe at https://www.dropbox.com/s/9422h1itk27yid4/Rnews_2002-3_bootstrap.pdf?dl=0]

OUTPUT EXAMPLE: (with otion -s -c -b)
>summary
      7.1e+03 7.1e+03 6.1e+03    0.55 4.9e-08 5e-08
>coefficients
  (Intercept)    -1.8   0.087     -21   1e-96
        coexp    0.69     0.1     6.6 3.2e-11
     tba_rank    0.17    0.12     1.4    0.16
>bootstrap
  (Intercept)      50    -1.8    -0.3   0.016
        coexp      50    0.69    0.48    0.87
     tba_rank      50    0.17  -0.084     0.4
          AUC      50    0.55    0.53    0.57
       Pvalue      50 4.9e-08 7.4e-09  0.0088


", option_list=option_list)

balancing_mle <- function(data,mle){
  t<-table(data[response_name])
  df<-as.data.frame(t)
	min<-min(df$Freq)
	out<-data.frame()
	for(i in levels(data[[response_name]])){
		data1=data[data[response_name]==i,]
		data1=data1[sample(nrow(data1), min), ]
		out<-rbind(out,data1)
	}
	return(out)
}

balancing <- function(data){
	t<-table(data[response_name])
	df<-as.data.frame(t)
	min<-min(df$Freq)
	out<-data.frame()
	for(i in levels(as.factor(data[[response_name]]))){
		data1=data[data[response_name]==i,]
		data1=data1[sample(nrow(data1), min), ]
		out<-rbind(out,data1)
	}
	return(out)
}

bs <- function(data, indices, formula) {
	d <- data[indices,] # allows boot to select sample
	model <- glm(formula, data=d, family = "binomial")
	data$PROB__ <- predict(model,type=c("response"))
	ROC <- roc(as.formula(paste(response_name,"~PROB__")), data = data, direction="<")
	AUC <- as.numeric(ROC$auc)
	Pvalue = wilcox.test(as.formula(paste("PROB__~",response_name)), data=data, alt="less")$p.value
	v=c(AUC,Pvalue)
	names(v)=c("AUC","Pvalue")
	retval=c(coef(model),v)
	return(retval)
}

#bs_linear <- function(formula, data, indices) {
#  d <- data[indices,] # allows boot to select sample
#  fit <- lm(formula, data=d)
#  return(coef(fit))
#}

arguments <- parse_args(parser, positional_arguments = 2)
opt<-arguments$options

#if(!opt$summary && !opt$coef && !nchar(opt$predict)>0){
#	stop("-s or -c or -p required")
#}
if (opt$linear && opt$bootstrap) {
	stop("Bootstrap not available for linear models")
}

infile <- arguments$args[1]
tryCatch( {
	formula <- as.formula(arguments$args[2])
}, error = function(err) {
	# todo differentiate between errors. For a single placeholder string is Error in eval(expr, envir, enclos) : object 'opt$formula' not found
	# would be nice to translate R cryptic errors about wrong formulas.
	if (nchar(opt$response_file) <= 0) {
		stop("You have not given a suitable formula and did not use -y!")
	}
})

response_name <- gsub("~.*","",arguments$args[2])

data <- read.table(infile, head=opt$header, sep="\t", row.names = NULL)
if(!opt$bootstrap && opt$balancing){
	data<-balancing(data)
}
if(nchar(opt$response_file)>0){
	response_data <- read.table(opt$response_file, head=opt$header, sep="\t", row.names = NULL)
	response_name <- colnames(response_data)
	if (!opt$lasso) {
		newdata <- cbind(data, response_data)
		if (opt$header) {
			formula <- as.formula(paste0(response_name, "~ ."))
		} else {
			colnames(newdata)[length(colnames(newdata))] <- "y"
			formula <- as.formula("y ~ .")
			response_name <- "y"
		}
		data <- newdata
	}
} else if (!opt$header) {
	stop("You can use files without headers only with -y!")
}

if (opt$linear) {
	family <- "gaussian"
} else {
	family <- "binomial"
}
if(!opt$lasso){
	if(!opt$univariate_only){
		if (opt$linear) {
			model <- lm(formula, data=data) # or we want to use glmet...many things need to be traced.
		} else {
			model <- glm(formula, data=data, family = family)
		}
		if (nchar(opt$out_model_file) > 0) {
			save(model,file=opt$out_model_file)
		}

		if(!opt$linear && (opt$summary || opt$bootstrap)) {
			# add a column with the predicted probabilities
			data$PROB__ <- predict(model,type=c("response"))
		}


		if(opt$summary){
			write(">summary",file="")
			if (opt$linear) {
				s <- summary(model)
				f <- s$fstatistic[1]
				df1 <- s$fstatistic[2]
				df2 <- s$fstatistic[3]
				Pvalue <- pf(f, df1, df2, low=FALSE)
				Pvalue_anova  <- anova(model, update(model, ~1))$`Pr(>F)`[2] 	# compare the model with a 
														# model that use only intercept, 
														# anova take correclty into account the number of predictor and DF
				out <- c(df1, df2, s$r.squared, s$adj.r.squared, Pvalue, Pvalue_anova)
			} else {
				ROC <- roc(as.formula(paste(response_name,"~PROB__")), data = data, direction="<")
				AUC <- as.numeric(ROC$auc)
				Pvalue_wilcox <- wilcox.test(as.formula(paste("PROB__~",response_name)), data=data, alt="less")$p.value
				Pvalue_anova  <- anova(model, update(model, ~1), test = "Chisq")$`Pr(>Chi)`[2] 	# compare the model with a 
														# model that use only intercept, 
														# anova take correclty into account the number of predictor and DF
				out <- c(model$df.null, model$df.residual, model$aic, AUC, Pvalue_wilcox, Pvalue_anova)
			}
			write(out, ncolumns=6, file="", sep="\t")
		}

		if(opt$coef){
			write(">coefficients",file="")
			#coeff <- data.frame(column,model$coeff)
			coeff <- summary(model)$coefficients
			write.table(coeff, file="", sep="\t", quote=FALSE, col=F, append=T)
		}

		# obtain the AUC also when the model is applied to new_data	
		# the results of the prediction on new_data should be the last thing in the output file
		if(opt$new_data!=""){
			new_data <- read.table(opt$new_data, head=T, sep="\t")
			# add a column with the predicted probabilities
			new_data$logistic_model_predicted_value <- predict(model, newdata=new_data, type = "response")
			write(">predict",file="")
			write.table(new_data, col.names = opt$col_names, row.names = FALSE, quote = FALSE, file="", sep="\t")
			# suppose that the response_name is the same in new_data
			ROC_predict<-roc(as.formula(paste(response_name,"~logistic_model_predicted_value")), data=new_data, direction="<")
			AUC_predict<-as.numeric(ROC_predict$auc)
			write(">predict_auc", file="")
			write(AUC_predict, file="")
		}

	}
	if(opt$univariate || opt$univariate_only){
	    write(">coefficients_univariate",file="")
	    for(x in attr(model$terms,"term.labels")){
		formula=sprintf("%s~%s",response_name,x)
		formula = as.formula(formula)
		if (opt$linear) {
			model <- lm(formula, data=data) # or we want to use glmet...many things need to be traced.
		} else {
			model <- glm(formula, data=data, family = family)
		}
		c <- summary(model)$coefficients
		c <- as.data.frame(c)
		c <- c[2,]#do not print the intercept
		if(!opt$linear){#if the model is logistic compute the AUC for each predictors
			prob_label <- sprintf('PROB__%s__only', x)
			data[prob_label] <- predict(model,type=c("response"))
			ROC <- roc(as.formula(sprintf("%s~%s",response_name,prob_label)), data = data, direction="<")
			AUC <- as.numeric(ROC$auc)
			c <- cbind(c,c(AUC))
		}
		write.table(c, file="", sep="\t", quote=FALSE, col=F, append=T)
	    }
	}


	if(opt$residuals){
		if(opt$summary || opt$coeff || opt$bootstrap > 0){
			write(">residuals",file="")
		}
		#residuals <- data.frame(column,model$residuals)
		residuals <- model$residuals
		write.table(residuals, file="", sep="\t", quote=FALSE, col=F, append=T)
	}

	if(opt$bootstrap>0 && !opt$linear){
	  if(opt$coef || opt$residuals || opt$summary){
	    write(">bootstrap",file="")
	  }
	  library(boot)
	  if(opt$balancing){
		  b <- boot(data=data, statistic=bs, formula=formula, R=opt$bootstrap, ran.gen=balancing, sim="parametric")
	  }else{
		  b <- boot(data=data, statistic=bs, formula=formula, R=opt$bootstrap)
	  }
	  labels <- names(b$t0)
	  for(i in 1:length(labels)){
	    ci <- boot.ci(b, type="perc", index=i)$percent
	    write.table(as.list(t(c(labels[i],opt$bootstrap,b$t0[i],ci[4],ci[5]))), file="", sep="\t", quote=FALSE, col=F, append=T, row.names = F)
	  }
	}
}else{
	suppressMessages(library(glmnet))
	if (nchar(opt$response_file) > 0) {
		response_data <- as.matrix(response_data)	
	} else {
		response_data <- as.matrix(data[, colnames(data) == response_name])
		data <- data[, colnames(data) != response_name]
	}
	data <- as.matrix(data)
	tryCatch  (
	   {
	      	model <- cv.glmnet(data, response_data, alpha=1,family=family)
	      	#model <- glmnet(data, response_data, alpha=1, lambda=cvfit$lambda.min, family=family)
		if(opt$coef){
			write(">coefficients",file="")
		      	coefm <- as.matrix(coef(model, s = "lambda.min"))[,1]   # s="lambda.min" recover the global model fitted with lambda.min, 
																		# it is not necessary to recoumpute the mode, it is inside model procucted by cv.glmnet
		      	NAs <- rep(NA, length(coefm))
	      		res <- data.frame(coef=coefm, stderr=NAs, tvalue=NAs, pvalue=NAs)
		      	write.table(res, sep="\t", quote=FALSE, col=FALSE, row=TRUE)
		}
		if(opt$summary){
			write(">summary",file="")
			if (opt$linear) {
				pvalue_ridge_glmnet <- function(ridge_fit, y, xm) {
				   SSE <- sum((y-predict(ridge_fit, xm))^2)
				   SSM <- sum((predict(ridge_fit, xm)-mean(y))^2)
				   non_zero <- predict(ridge_fit, type="nonzero")
				   DFM <- length(non_zero[[1]])
				   #or directly n. of x?
				   DFE <- length(y) - DFM - 1
				   MSM <- SSM/DFM
				   MSE <- SSE/DFE
				   Fs <- MSM/MSE
				   pv <- pf(MSM/MSE,DFM,DFE,low=F)
				   return(c(Fs, DFM, DFE, pv))
				}
				# todo CHECK PVALUE OF LASSO MODELS!
				out <- c(NA, model$glmnet.fit$df[model$lambda==model$lambda.min], model$glmnet.fit$dev.ratio[model$lambda==model$lambda.min], model$glmnet.fit$aic, NA)
			} else {
				predicted_probabilities <- as.vector(predict(model,newx=data,type=c("response"),s = "lambda.min"))
				response_data <- as.factor(response_data)
				#R2 <- model$dev.ratio  #consider this value instead of AUC to improove performance
				ROC <- roc(as.vector(response_data), predicted_probabilities, direction="<")
				AUC <- as.numeric(ROC$auc)
				x <- predicted_probabilities[response_data==levels(response_data)[1]]
				y <- predicted_probabilities[response_data==levels(response_data)[2]]
				Pvalue <- wilcox.test(x, y, alt="less")$p.value
				out <- c(NA, model$glmnet.fit$df[model$lambda==model$lambda.min], NA, AUC, Pvalue)
			}
			write(out, ncolumns=6, file="", sep="\t")
		}
	   },
	   error=function(cond) {
	      message(cond)
	      write.table(as.data.frame(matrix(nrow=ncol(data)-1,ncol=3)), sep="\t", quote=FALSE, col=FALSE, row=FALSE)
	   }
	)
}

w=warnings()
sink(stderr())
if(!is.null(w)){
	print(w)
}
sink()
