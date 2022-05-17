#!

library("pROC")
library("ggplot2")
library("ggpubr")
library("ggsci")

#setwd("/sto1/epigen/TPXcCRE/dataset/v8_ChIRP_neg_rand")

z<-read.table("raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.logistic.bestParams_bestPredictor.matrix_reduce.header_added"
              , header = T)

filter_table<-read.table("bestParams_bestPredictor.tsv"
                   , header = F
                   , col.names = c("params","ssRNA","predictor"))

### Split dataframe according to levels
z_split<-split(z,z$ssRNA)

### Select the best predictor for each ssRNA through filter file 
z_split_predictors<-lapply(z_split,function(x){
  ssRNA_name<-levels(as.factor(x$ssRNA))
  predictor_name<-filter_table[filter_table$ssRNA==ssRNA_name,c("predictor")]
  z_split_filt<-x[,c("ssRNA","neg_pos",predictor_name)]
})

### Create roc objects with defined direction
roc_matrix<-lapply(z_split_predictors, function(x){
  ssRNA_name<-levels(as.factor(x$ssRNA))
  predictor_name<-filter_table[filter_table$ssRNA==ssRNA_name,c("predictor")]
  suppressMessages(roc(as.formula(paste(collapse = "~", c("neg_pos", predictor_name))), data= x, direction="<"))
})

### AUC values to print
labels <- levels(as.factor(z$ssRNA))
for(i in seq_along(labels)){
  l<-labels[[i]]
  auc<-round(roc_matrix[[l]]$auc,3)
  labels[[i]]<-paste0(c(l,auc),collapse=" - ")
}

### ROC plot
p<-ggroc(roc_matrix) +
  theme_classic(base_size = 10) +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="darkgray", linetype="dashed",size=0.3) +
  theme(plot.margin = margin(1,2,1,2, "cm")
        , plot.title = element_text(hjust = 0.5)
        , legend.title = element_blank()
        , legend.position = c(0.85, 0.3)
        , legend.key.size = unit(0.4, "cm"))+
  coord_fixed() +
  scale_color_igv(labels = labels)

### Save plot
ggsave("ROC_bestParams.pdf", device = "pdf", plot = p, width = 6, height = 6)


