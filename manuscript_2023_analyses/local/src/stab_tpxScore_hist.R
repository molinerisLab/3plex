##############################
#
#   Plot best TPXcCRE score calcolati con score di triplexator e con stability
#
library("ggplot2")
library("ggpubr")
library("ggsci")

setwd("/sto1/epigen/TPXcCRE/dataset")

z1<-read.table("v1/cCRE.tpx.best.complete.fisher_select_cutoff.ALLcond.best.selected_ssRNA_id.header_added",header=T)
z2<-read.table("v5_meanStability/cCRE.tpx.best.complete.fisher_select_cutoff.ALLcond.best.selected_ssRNA_id.header_added",header=T)

z1$type<-'tpx_score'
z2$type<-'longTarget_stab'

z_merge<-rbind(z2,z1)




################
# Plots TPXcCRE score

# Dodged hist
ggplot(z_merge,aes(x=TPXcCRE_score, fill=type)) +
  geom_histogram(binwidth = 4,alpha=0.7,position = "dodge") +
  theme_pubr(base_size = 10) +
  labs(fill="",x="TPXcCRE score") +
  theme(legend.key.size = unit(3.5,"mm"))


#************************
# DRAFT
# Density
ggplot(z_merge,aes(x=TPXcCRE_score, fill=type,color=type)) +
  geom_density(lwd=0.5,adjust = 0.1,alpha=0.5) +
  theme_pubr(base_size = 10) +
  labs(fill="",x="TPXcCRE score",color="") +
  theme(legend.key.size = unit(3.5,"mm"))
#*************************


# Hist pos vs neg values
z_merge %>%
  filter(TPXcCRE_score<0) %>%
  ggplot(aes(x=TPXcCRE_score, fill=type,position = "dodge")) +
    geom_histogram(binwidth = 3,alpha=0.7) +
    theme_pubr(base_size = 10) +
    labs(fill="",x="TPXcCRE score") +
    theme(legend.key.size = unit(3.5,"mm"))

# Hist divided by conditions and cCRE type
ggplot(z_merge,aes(x=TPXcCRE_score, fill=type)) +
  geom_histogram(binwidth = 20,alpha=0.7,position = "dodge") +
  theme_classic(base_size = 10) +
  labs(fill="",x="TPXcCRE score") +
  theme(legend.key.size = unit(3.5,"mm")) +
  facet_grid(condition~cCRE,scales = "free_y") +
  scale_fill_igv()
#theme_pubclean(base_size = 10)

ggplot(z_merge,aes(x=TPXcCRE_score, fill=type)) +
  geom_density(lwd=0.2,adjust = 0.3,alpha=0.5) +
  theme_classic(base_size = 10) +
  labs(fill="",x="TPXcCRE score") +
  theme(legend.key.size = unit(3.5,"mm")) +
  facet_grid(condition~cCRE,scales = "free_y") +
  scale_fill_igv()
  



################
#  Scatterplot

# x --> TPXcCRE score da tpx_score
# y --> TPXcCRE score da stab

z_scatter<-as.data.frame(cbind(z1$TPXcCRE_score,z2$TPXcCRE_score))
ggplot(z_scatter,aes(x=V1, y=V2)) +
  geom_point(size=0.8)+
  labs(x="TPXcCRE score from tpx score", y="TPXcCRE score from LongTarget stability") +
  theme_pubr(base_size = 10) +
  geom_smooth(method="lm", se=TRUE, fullrange=FALSE, level=0.95,size=0.5,col="#d62828")



################
# Plots odds ratio

# Dodged hist
ggplot(z_merge%>%filter(odds_ratio>1),aes(x=odds_ratio, fill=type)) +
  geom_histogram(binwidth=0.05,alpha=0.7,position = "dodge") +
  theme_pubr(base_size = 10) +
  labs(fill="") +
  theme(legend.key.size = unit(3.5,"mm"))+
  facet_wrap(~type)

# Density
ggplot(z_merge,aes(x=odds_ratio, fill=type,color=type)) +
  geom_density(lwd=0.5,adjust = 0.5,alpha=0.6) +
  theme_pubr(base_size = 10) +
  labs(fill="",color="",x="odds ratio") +
  theme(legend.key.size = unit(3.5,"mm"))

#  facet_wrap(~type)
