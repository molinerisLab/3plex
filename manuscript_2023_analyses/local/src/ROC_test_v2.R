library(pROC)

z <- read.table("v8_ChIRP_neg_rand/tpx_paramspace/TERC_ss10_singleNt/TERC.neg_pos_rand.bed/min_length~10/max_length~30/error_rate~20/guanine_rate~70/filter_repeat~on/consecutive_errors~3/raw.tpx.custom_summary.neg_pos.covered_by_tts.stability.logistic"
                , header = T
                )

# With this correction, the AUC will be 0.5 if non discriminant and 1.0
# if maximal, whatever the region defined
# Note that this correction is undefined for curves below the diagonal (auc < min). 
# Attempting to correct such an AUC will return NA with a warning.

roc_partial_auc <- roc(z$neg_pos
            , z$Stability_best
            , percent=FALSE
            # arguments for auc
            , partial.auc=c(1, .9)
            , partial.auc.correct=TRUE
            , partial.auc.focus="spec"
            # arguments for ci
            , ci=TRUE
            , boot.n=100
            , ci.alpha=0.9
            , stratified=FALSE
            # arguments for plot
            , plot=TRUE
            , auc.polygon=TRUE
            , max.auc.polygon=TRUE
            , grid=TRUE
            , print.auc=TRUE 
            , show.thres=TRUE
)

roc_partial_auc_formula <- roc(z$neg_pos ~ z$Stability_best + z$Stability_tot_undercount
                       , percent = F
                       # arguments for auc
                       , partial.auc=c(1, .9)
                       , partial.auc.correct=TRUE
                       , partial.auc.focus="spec"
                       # arguments for ci
                       , ci=TRUE
                       , boot.n=100
                       , ci.alpha=0.9
                       , stratified=FALSE
                       # arguments for plot
                       , plot=TRUE
                       , auc.polygon=TRUE
                       , max.auc.polygon=TRUE
                       , grid=TRUE
                       , print.auc=TRUE, 
                       show.thres=TRUE
)


roc_partial_auc_formula2 <- roc(z$neg_pos ~ z$Stability_best + z$Stability_tot_undercount
                               , percent = F
                               # arguments for auc
                               , partial.auc=c(1, .9)
                               , partial.auc.correct=TRUE
                               , partial.auc.focus="spec"
                               # arguments for ci
                               , ci=TRUE
                               , boot.n=100
                               , ci.alpha=0.9
                               , stratified=FALSE)

#ggroc(roc_partial_auc_formula2)

plot.roc(df$neg_pos
         , df$PROB_stability_comb_rank
         , partial.auc = c(1,.9)
         , partial.auc.correct = T
         , print.auc=T
         , auc.polygon=TRUE
         , print.auc.y=.2
         , print.auc.col="#1c61b6"
         , col="#1c61b6"
         , print.auc.pattern="comb pAUC (100-90%% SP):\n%.3f%%"
)







