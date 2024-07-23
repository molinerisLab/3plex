library(rjson)
library(tidyverse)
library(reshape2)
library(ggridges)

transcr_length <- 2107

# Read data ----
dbd <- read.delim("SRA1_ssmasked-test.dbd",header = F)
profile_json <- fromJSON(file="SRA1-test-profile_flat.json")
profile_random_json <- fromJSON(file="SRA1-test__profile_range.random.json")


# Observed profile ----
profile <- reshape2::melt(profile_json$profiles)%>%
  mutate_at("L1", as.numeric) %>%
  complete(L1 = seq(1, transcr_length, by=1), 
           fill=list(value=0))

dbd$mean_TTS_count <- sapply(1:nrow(dbd), function(i){
  mean(profile[dbd[i,"V2"]:dbd[i,"V3"],"value",drop=T])
})


# Random profile ----
profile_random <- dcast(melt(profile_random_json), L1 ~ L2)
profile_random <- profile_random %>%
  mutate_at("L1", as.numeric) %>%
  complete(L1 = seq(1, transcr_length, by=1))
  

# boxplot for each nt ---
# pl <- ggplot(profile_random_2[1:50,], aes(x=L1,group=L1)) +
#   geom_boxplot(aes(ymin = min, lower = lower_quartile, middle = median, upper = upper_quartile, ymax = max),
#                stat = "identity")
# LENTO !!!

# Mean, stdev, p-value
dbd$mean_TTS_count_random <- sapply(1:nrow(dbd), function(i){
  row_idx <- which(profile_random$L1==dbd[i,"V2"]):which(profile_random$L1==dbd[i,"V3"])
  mean(profile_random[row_idx, "avg",drop=T])
})

dbd$var_TTS_count_random <- sapply(1:nrow(dbd), function(i){
  row_idx <- which(profile_random$L1==dbd[i,"V2"]):which(profile_random$L1==dbd[i,"V3"])
  sum_scarti <- (sum(profile_random[row_idx,"avg",drop=T]-dbd[i,"mean_TTS_count_random"]))**2
  var_sum <- sum(profile_random[row_idx,"var",drop=T])
  sqrt((var_sum + sum_scarti)/transcr_length)
})

dbd$pval_TTS_count <- pnorm(dbd$mean_TTS_count, 
                            mean=dbd$mean_TTS_count_random, 
                            sd=dbd$var_TTS_count_random, 
                            lower.tail=FALSE)

dbd$V1 <- factor(dbd$V1, levels = dbd$V1[order(dbd$V2,decreasing = T)]) 


# plot norm dist and observed val ----
max_lim <- max(dbd$mean_TTS_count)
toplot <- dbd %>% 
  uncount(500) %>%
  mutate(value = rnorm(n(), mean_TTS_count_random, var_TTS_count_random))
ggplot(toplot, aes(x = value, y = V1, fill=-log10(pval_TTS_count+.0000000001))) + 
  geom_density_ridges(linewidth=.5) +
  geom_point(data=dbd, mapping=aes(x=mean_TTS_count, y=V1), color = "red") +
  theme_bw() +
  geom_hline(yintercept = dbd$V1)
  


# Stability ----
profile <- profile_json$best_stability
dbd$mean_stability <- sapply(1:nrow(dbd), function(i){
  mean(profile[(dbd[i,"V2"]+1):(dbd[i,"V3"]+1)])
})

profile_random_json <- fromJSON(file="SRA1-test__stability_profile.random.json")

profile_random <- reshape2::melt(profile_random_json$all)
profile_random_avg <- profile_random %>% 
  filter(profile_random$L2=="avg") %>%
  mutate_at("L1", as.numeric) %>%
  complete(L1 = seq(1, transcr_length, by=1)) %>%
  select(-L2)

profile_random_avg$L1 <- profile_random_avg$L1-1





