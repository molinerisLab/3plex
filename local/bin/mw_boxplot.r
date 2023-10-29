#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(rstatix))

option_list <- list( 
  make_option(c("-t", "--tpxlist"), action="store", default=NULL,
              help="Matrix of ssRNA tpx calculated with 3plex on target and background regions."),
  make_option(c("-l", "--tblist"), action="store", default=NULL,
              help="Target/background list path [ssRNA, gene_name, class]."),
  make_option(c("-s", "--score"), action="store", default="Stability_best",
              help="TPX score used to rank target/background regions. Available: Stability_best, Stability_tot, Score_best, Stability_norm."),
)
opt <- parse_args(OptionParser(option_list = option_list))

# read 3plex output 
message(paste0("--- Reading: ",opt$tpxlist))
tpx <- read.delim(opt$tpxlist, header=T)

# read target and background list
message(paste0("--- Reading: ", opt$tblist))
tblist <- read.delim(opt$tblist, header = F)

# retrieve ssRNA name and extract duplex_id name
ssRNA <- unique(tpx$Sequence_ID)
tpx$Duplex_ID_gene <- gsub(pattern = "(.*)::(.*)","\\1",tpx$Duplex_ID)

# subset target-background list (tpx is not to be subset because "*.add_seros.gz" files ar ssRNA specific)
tblist_ssRNA <- tblist[tblist$V1==ssRNA,]
# add target/background labels
# this merging operation removes target/background genes without available tss
m <- merge(tpx,tblist_ssRNA, by.x="Duplex_ID_gene", by.y="V2")
# stat.test
stat.test <- m %>% 
  wilcox_test(as.formula(paste0(opt$score," ~ V3")), 
                               ref.group = "background",
                               alternative = "less") %>%
  add_significance() %>%
  add_xy_position(x = "V3")
# boxplot
bxp <- ggboxplot(m, x = "V3", y = opt$score, 
                 fill = "V3", palette = "d3", outlier.shape = NA) +
  stat_pvalue_manual(stat.test, label = "p.signif") +
  theme_bw() + theme(legend.position = "none") +
  #geom_jitter(width = .1, height = NULL, size =.5, alpha=.2) +
  xlab("") + ggtitle(ssRNA)

# save image
outfile <- sub(opt$tpxlist,pattern = ".gz",replacement = ".boxplot.pdf")
message(paste0("--- Saving to: ",outfile))
pdf(file = outfile, paper="a4", width = 3, height = 3)
bxp
dev.off()