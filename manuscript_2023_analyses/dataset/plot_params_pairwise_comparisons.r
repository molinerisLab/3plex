#######################################
#
# params_pairwise_comparisons.pdf 
#

# 0. Resources ----

library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggsci)
library(tidyr)
library(reshape2)



# 1. Read and manipulate dataframe ----

df <- read.table("tpx_paramspace_AUC.method.Npeaks_version.technique.header_added.gz", header=T)

# Set params as factors and relevel 
groups_factors <- c("peak_method_all","ssRNA","singleStrandedness","min_len","error_rate","guanine_rate","repeat_filter","consecutive_error","predictor","technique")
df[,groups_factors] <- lapply(df[,groups_factors] , factor)
df$predictor <- relevel(df$predictor, ref = "t_pot_norm")
df$technique <- relevel(df$technique, ref = "ChIRP-seq")
df$min_len <- relevel(df$min_len, ref = "8")
df$singleStrandedness <- relevel(df$singleStrandedness, ref = "ss50")

# Filter df for TPX ssRNA
triplex_evidence_ssRNA <- c("TERC", "NEAT1", "HOTAIR", "MEG3", "Meg3", "CDKN2B-AS1", "lncSmad7")
df <- df %>% filter(ssRNA %in% triplex_evidence_ssRNA)

# Create df_min_len with default_param dataframe with min len 8, 10, 12, 16
df_min_len <- df %>% filter(singleStrandedness == "ss0"
                            , guanine_rate == 40
                            , error_rate == 20
                            , repeat_filter == "on"
                            , consecutive_error == 1
                            , predictor == "t_pot_norm")

# Remove min len 12, 16 from df
df <- df %>% filter(min_len != 12 & min_len != 16 & predictor != "TTS_covered_frac")
df$min_len <- droplevels(df$min_len)
df$predictor <- droplevels(df$predictor)





##################################################
#
# 2. Tukey's test version (not used in the paper) ----
#

# anova model
groups <- colnames(df)[!colnames(df) %in% c("AUC","technique")]
formula_ssRNA <- as.formula(paste("AUC ~ ", paste(groups, collapse="+")))
aov_ssRNA <- aov(formula_ssRNA, data = df)

# Compute tukey's diff
# https://stackoverflow.com/questions/33644034/how-to-visualize-pairwise-comparisons-with-ggplot2
pairwise_groups <- c("peak_method_all", "singleStrandedness", "predictor", "error_rate", "guanine_rate", "repeat_filter", "consecutive_error")
tky_df <- data.frame(diff=double(), lwr=double(), upr=double(), p_adj=double(), pair=factor(), param=factor())

for(group in pairwise_groups){
  print(group)
  tky <- as.data.frame(TukeyHSD(aov_ssRNA, group, ordered = F)[[group]])
  tky$pair <- rownames(tky)
  tky <- cbind(tky, rep(group, length(tky$diff)))
  print(tky)
  tky_df <- rbind(tky_df, tky)
}
colnames(tky_df)[colnames(tky_df) == 'rep(group, length(tky$diff))'] <- 'param'

# Add min_len values (8, 10, 12, 16)
groups_min_len <- c("peak_method_all", "ssRNA", "n_peaks", "min_len")
formula_min_len <- as.formula(paste("AUC ~ ", paste(groups_min_len, collapse="+")))
model_min_len <- lm(formula_min_len, data = df_min_len)
aov_min_len <- aov(formula_min_len, data = df_min_len)

tky <- as.data.frame(TukeyHSD(aov_min_len, "min_len", ordered = F)[["min_len"]])
tky$pair <- rownames(tky)
tky <- cbind(tky, param=rep("min_len", length(tky$diff)))
tky_df <- rbind(tky_df, tky)

# reorder tky_df
tky_df$pair <- reorder(tky_df$pair, -tky_df$diff)

# dotplot
p <- ggplot(tky_df) +
  geom_hline(yintercept=0, lty="11", colour="grey30") +
  geom_errorbar(aes(pair, ymin=lwr, ymax=upr), width=0.2, color = "grey10") +
  geom_point(aes(pair, diff), color = "#6F99ADFF") +
  labs(colour="") +
  facet_grid(.~param, scales = "free_x", space = "free_x") +
  theme_classic() +
  xlab("") + ylab("") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust=1))

# not saved
# ggsave("params_pairwise_comparisons.pdf", plot = p, device = "pdf", width = 7, height = 5)






###############################
#
# 3. Compute AUC mean and sd on all TPX paramspace for fixed level of each param ----
#

df_param <- df_min_len %>% 
  group_by_("min_len") %>% 
  summarise(mean_AUC = mean(AUC), sd_AUC = sd(AUC)) %>%
  gather(key="param", value="min_len", 1)
colnames(df_param)[colnames(df_param) == "min_len"] <- 'values'

pairwise_groups <- c("peak_method_all", "singleStrandedness", "predictor", "error_rate", "guanine_rate", "repeat_filter", "consecutive_error")

for (group in pairwise_groups) {
  nr <- df %>% 
    group_by_(group) %>% 
    summarise(mean_AUC = mean(AUC), sd_AUC = sd(AUC)) %>%
    gather(key="param", value=group, 1)
  colnames(nr)[colnames(nr) == "group"] <- 'values'
  df_param <- rbind(df_param, nr)
}

# rename params and levels
df_param$param <- gsub("min_len", "minimal length", df_param$param)
df_param$param <- gsub("peak_method_all", "peak filtering method", df_param$param)
df_param$param <- gsub("predictor", "scoring method", df_param$param)
df_param$param <- gsub("singleStrandedness", "single strandedness", df_param$param)
df_param$param <- gsub("_", " ", df_param$param)

df_param$values <- gsub("idr_overlap_top1000", "IDR + overlap (top 1000)", df_param$values)
df_param$values <- gsub("idr", "IDR", df_param$values)
df_param$values <- gsub("t_pot_norm", "Triplexator potential", df_param$values)
df_param$values <- gsub("Score_best", "Triplexator best score", df_param$values)
df_param$values <- gsub("Stability_best", "best stability", df_param$values)
df_param$values <- gsub("Stability_norm_undercount", "normalised stability", df_param$values)
df_param$values <- gsub("_", " ", df_param$values)

df_param$values <- as.factor(df_param$values)
df_param$values <- reorder(df_param$values, -df_param$mean_AUC)

# dotplot
p <- ggplot(df_param) +
  geom_errorbar(aes(y = mean_AUC, x =values, ymin=mean_AUC-sd_AUC/2, ymax=mean_AUC+sd_AUC/2), width=0.2, color = "grey20") +
  geom_point(aes(y=mean_AUC, x = values),color = "#6F99ADFF", size = 3) +
  facet_grid(.~param, scales = "free_x", space = "free_x") +
  theme_pubclean() +
  xlab("") + ylab("mean AUC") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

# save dotplot
ggsave("params_pairwise_comparisons.pdf", plot = p, device = "pdf", width = 7, height = 5)



####################
#
# Boxplot paired points
#

# Dataframe
df <- read.table("tpx_paramspace_AUC.method.Npeaks_version.technique.header_added.gz", header=T)

# Filter df for TPX ssRNA
triplex_evidence_ssRNA <- c("TERC", "NEAT1", "HOTAIR", "MEG3", "Meg3", "CDKN2B-AS1", "lncSmad7")
df <- df %>% filter(ssRNA %in% triplex_evidence_ssRNA)

# dafult params: 
# single strandedness --> ss0 
# error rate --> 20
# guanine rate --> 40
# repeat filter --> on
# consecutive errors --> 1 
# min len --> 8 or 10
df_min_len <- df %>% filter(singleStrandedness == "ss0"
                            , guanine_rate == 40
                            , error_rate == 20
                            , repeat_filter == "on"
                            , consecutive_error == 1
                            , predictor == "t_pot_norm") %>% group_by(ssRNA,min_len) %>% summarise(AUC = max(AUC))
df_min_len <- reshape2::melt(df_min_len, measure.var="min_len")
df_min_len$value <- factor(df_min_len$value, levels = c("8","10","12","16"))

  
df_error_rate <- df %>% filter(singleStrandedness == "ss0"
                            , guanine_rate == 40
                            , min_len == 8 | 10
                            , repeat_filter == "on"
                            , consecutive_error == 1
                            , predictor == "t_pot_norm") %>% group_by(ssRNA,error_rate) %>% summarise(AUC = max(AUC))
df_error_rate <- reshape2::melt(df_error_rate, measure.var="error_rate")

df_cons_error <- df %>% filter(singleStrandedness == "ss0"
                               , guanine_rate == 40
                               , min_len == 8 | 10
                               , repeat_filter == "on"
                               , error_rate ==20
                               , predictor == "t_pot_norm") %>% group_by(ssRNA, consecutive_error) %>% summarise(AUC = max(AUC))
df_cons_error <- reshape2::melt(df_cons_error, measure.var="consecutive_error")

df_ss <- df %>% filter(consecutive_error == 1
                       , guanine_rate == 40
                       , min_len == 8 | 10
                       , repeat_filter == "on"
                       , error_rate ==20
                       , predictor == "t_pot_norm") %>% group_by(ssRNA,singleStrandedness) %>% summarise(AUC = max(AUC))
df_ss <- reshape2::melt(df_ss, measure.var="singleStrandedness")
df_ss$value <- factor(df_ss$value, levels = c("ss0","ss10","ss20","ss50"))

df_g_rate <- df %>% filter(consecutive_error == 1
                       , singleStrandedness == "ss0"
                       , min_len == 8 | 10
                       , repeat_filter == "on"
                       , error_rate ==20
                       , predictor == "t_pot_norm") %>% group_by(ssRNA,guanine_rate) %>% summarise(AUC = max(AUC))
df_g_rate <- reshape2::melt(df_g_rate, measure.var="guanine_rate")

df_rep_filt <- df %>% filter(consecutive_error == 1
                           , singleStrandedness == "ss0"
                           , min_len == 8 | 10
                           , guanine_rate == 40
                           , error_rate ==20
                           , predictor == "t_pot_norm") %>% group_by(ssRNA,repeat_filter) %>% summarise(AUC = max(AUC))
df_rep_filt <- reshape2::melt(df_rep_filt, measure.var="repeat_filter")

df_score <- df %>% filter(consecutive_error == 1
                             , singleStrandedness == "ss0"
                             , min_len == 8 | 10
                             , guanine_rate == 40
                             , error_rate ==20
                             , repeat_filter == "on"
                             , predictor != "TTS_covered_frac") %>% group_by(ssRNA,predictor) %>% summarise(AUC = max(AUC))
df_score <- reshape2::melt(df_score, measure.var="predictor")
df_score$value <- factor(df_score$value)

df_peak_met <- df %>% filter(consecutive_error == 1
                          , singleStrandedness == "ss0"
                          , min_len == 8 | 10
                          , guanine_rate == 40
                          , error_rate ==20
                          , repeat_filter == "on"
                          , predictor == "t_pot_norm") %>% group_by(ssRNA, peak_method_all) %>% summarise(AUC = max(AUC))
df_peak_met <- reshape2::melt(df_peak_met, measure.var="peak_method_all")

df_boxplot <- rbind(df_min_len, df_error_rate, df_cons_error, df_ss, df_g_rate, df_rep_filt, df_score, df_peak_met)

p <- ggplot(df_boxplot, aes(x=value,y=AUC)) +
  facet_grid(.~variable, scales = "free") +
  geom_boxplot(alpha=0.6) +
  geom_line(aes(group = ssRNA), size = 0.2, alpha = 0.7) +
  theme_pubclean() +
  xlab("") + ylab("") +
  theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust=1))

# not saved









