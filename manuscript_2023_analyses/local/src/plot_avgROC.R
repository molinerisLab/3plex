#!/usr/bin/env Rscript 

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(pROC))

# Read stdin 
stdin <- file("stdin")
open(stdin, blocking=TRUE)
roc<-read.table(stdin, header=T, sep="\t")
roc<- read.delim("HOTAIR.triplexAligner.summary.clean.sens_spec.tsv",head=T)

#bawk 'NR>1 {print $specificities}' HOTAIR.triplexAligner.summary.clean.sens_spec.tsv | 
#  sort -k1,1g | 
#  uniq | 
#  enumerate_rows | 
#  perl -lane 'BEGIN{$,="\t"; $i=0; print ">$i"; $i++} $F[0] = $F[0] % int((2216/20)-1); 
#if($F[0]==0){print ">$i"; $i++} print $F[1]' | 
#  fasta2tab

roc$spec_quant<-quantile(roc$specificities)

write.table(roc,file="pippo",sep="\t")
roc$spec_rank <- rank(roc$specificities, ties.method = "first") 
N<-length(unique(roc$spec_rank))
B<-100
n<-N/B
roc$spec_bin<-as.integer(roc$spec_rank/n)

#roc_bin <- data.frame(specificities=unique(roc[,"specificities",drop=F] %>% arrange(specificities))$specificities)


avg_roc2 <- avg_roc %>% 
  group_by(iteration, spec_bin) %>% 
  summarise(mean_spec = mean(spec_bin_avg), 
            mean_sens = mean(sensitivities)) 

avg_roc2 <- avg_roc2 %>%
  group_by(spec_bin) %>%
  summarise(mean_mean_sens = mean(mean_sens)
            , mean_spec = mean(mean_spec)
            , ymin = mean(mean_sens) - sd(mean_sens)/2
            , ymax = mean(mean_sens) + sd(mean_sens)/2
            , ymin2 = min(mean_sens)
            , ymax2 = max(mean_sens)
            , median_sens = median(mean_sens)
            , ymin3 = quantile(mean_sens, probs = c(0.25))[1]
            , ymax3 = quantile(mean_sens, probs = c(0.75))[1])

# aggiungere 0 e 1
avg_roc2 <- rbind(avg_roc2, c(0,1,0,1,1,1,1,1,1,1))
avg_roc2 <- rbind(avg_roc2, c(1,0,1,0,0,0,0,0,0,0))



# plot
p <- ggplot(avg_roc2) +
  scale_x_reverse() +
  geom_ribbon(aes(x = mean_spec, y = mean_mean_sens, ymin = ymin2, ymax = ymax2), fill = "grey70", alpha = 0.4) +
  geom_ribbon(aes(x = mean_spec, y = mean_mean_sens, ymin = ymin, ymax = ymax), fill = "grey70", alpha = 0.4) +
  geom_line(aes(x = mean_spec, y = mean_mean_sens)) +
  theme_classic() +
  coord_fixed() +
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="darkgray", linetype="dashed", size = 0.3) +
  theme(plot.margin = margin(1,2,1,2, "cm"))

print(p)

ggsave(filename = "average_TFs_ROC.smoothed.pdf", plot = p, device = "pdf")
