#!/usr/bin/env Rscript

#options(error = function() traceback(2))
#suppressMessages(library("optparse"))

library(ggplot2)

#parser<-OptionParser(usage = "%prog < INFILE1
#.META: INFILE1
#        tpx_analysis/ssRNA_NAME/cCRE.tpx.best.complete.CONDITION_neg_pos.fisher_select_cutoff
#
#", option_list=option_list)

stdin <- file("stdin")
open(stdin, blocking=TRUE) # http://stackoverflow.com/questions/9370609/piping-stdin-to-r

z<-read.table(stdin,header=F,col.names=c("lncRNA","cCRE","tpx_score","greater_pos","lower_pos","greater_neg","lower_neg","oddsratio","pvalue","TPXcCRE_score"))

p<-ggplot(data=z, aes(x=tpx_score,y=TPXcCRE_score, fill=cCRE))+
  geom_col()+
  facet_grid(.~cCRE, scale="free_y")+
  theme_linedraw()+
  xlab("TPX score") + ylab("TPX-cCRE score")+
  theme(text=element_text(size=11),axis.text.x=element_text(angle=90,vjust=0.5,hjust=1),legend.position='none')+
  xlim(1,NA)+
  geom_hline(yintercept=-log10(0.05),linetype="dashed")

ggsave("/dev/stdout",plot=p,device="pdf")
                
w=warnings()
sink(stderr())
if(!is.null(w)){
        print(w)
}
sink()

