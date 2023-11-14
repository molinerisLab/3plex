#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(rstatix))

option_list <- list( 
  make_option(c("-t", "--tpxlist"), action="store", default=NULL,
              help="Matrix of ssRNA tpx calculated with 3plex on target and background regions."),
  make_option(c("-l", "--tblist"), action="store", default=NULL,
              help="Target/background list path [ gene_name, class ]."),
  make_option(c("-s", "--score"), action="store", default="Stability_best",
              help="TPX score used to rank target/background regions. Available: Stability_best, Stability_tot, Score_best, Stability_norm. [ Default: Stability_norm ]."),
  make_option(c("-d", "--directory"), action ="store", default=".",
              help="Output directory [ Default: . ].")
)
opt <- parse_args(OptionParser(option_list = option_list))

# create output directory
if(!dir.exists(opt$directory)) dir.create(opt$directory, recursive = T, showWarnings = FALSE)

# read 3plex output 
message(paste0("--- Reading: ",opt$tpxlist))
tpx <- read.delim(opt$tpxlist, header=T)

# read target and background list
message(paste0("--- Reading: ", opt$tblist))
tblist <- read.delim(opt$tblist, header=F, col.names=c("GeneID","class"))

# retrieve ssRNA name and extract duplex_id name
ssRNA <- unique(tpx$Sequence_ID)
tpx$GeneID <- gsub(pattern = "(.*)::(.*)","\\1",tpx$Duplex_ID)

# add target/background labels
# this merging operation removes target/background genes without available tss
m <- merge(tpx, tblist, by="GeneID")

# stat.test ----
stat.test <- m %>% 
  wilcox_test(as.formula(paste0(opt$score," ~ class")), 
                               ref.group = "background",
                               alternative = "less") %>%
  add_significance() %>%
  add_xy_position(x = "class")
# save table
outfile <- paste0(opt$directory, "/stability_comp.tsv")
message(paste0("--- Saving to: ",outfile))
write.table(stat.test[,2:8], outfile, sep = "\t", row.names = F, quote = F)

# boxplot ----
bxp <- ggboxplot(m, x = "class", y = opt$score, 
                 fill = "class", palette = "d3", outlier.shape = NA) +
  stat_pvalue_manual(stat.test, label = "p.signif") +
  theme_bw() + theme(legend.position = "none") +
  #geom_jitter(width = .1, height = NULL, size =.5, alpha=.2) +
  xlab("") + ggtitle(ssRNA)
# save image
outfile <- paste0(opt$directory, "/stability_comp_boxplot.pdf")
message(paste0("--- Saving to: ",outfile))
pdf(file = outfile, paper="a4", width = 3, height = 3.5)
bxp
dev.off()