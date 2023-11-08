#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggpubr))

option_list <- list( 
  make_option(c("-t", "--tpxlist"), action="store", default=NULL,
              help="Matrix of all ssRNA tpx calculated with 3plex on target and background regions."),
  make_option(c("-l", "--tblist"), action="store", default=NULL,
              help="Target/background list path [ssRNA, gene_name, class]."),
  make_option(c("-s", "--score"), action="store", default="Stability_best",
              help="TPX score used to rank target/background regions. Available: Stability_best, Stability_tot, Score_best, Stability_norm."),
  make_option(c("-d", "--outdir"), action="store", default=".",
              help="Output directiory.")
)
opt <- parse_args(OptionParser(option_list = option_list))

# read 3plex output 
message(paste0("--- Reading: ",opt$tpxlist))
z <- read.delim(opt$tpxlist, header=T)
z$Duplex_ID_gene <- gsub(pattern = "(.*)::(.*)","\\1",z$Duplex_ID)

# read target and background list
message(paste0("--- Reading: ", opt$tblist))
tblist <- read.delim(opt$tblist, header = F)

# create outdf and boxplot list
outdf <- data.frame(ssRNA=character(), pval=double(), median_diff=double())
outdf_cols <- colnames(outdf)

for (ssRNA in unique(z$Sequence_ID)) {
  #filter for ssRNA
  tblist_ssRNA <- tblist[tblist$V1==ssRNA,]
  z_ssRNA <- z[z$Sequence_ID==ssRNA,]
  # remove backgrounds that are also targets
  tblist_ssRNA_ordered <- tblist_ssRNA[order(tblist_ssRNA$V2, decreasing = TRUE), ]
  tblist_ssRNA <- tblist_ssRNA_ordered[!duplicated(tblist_ssRNA_ordered$V2), ]
  # add target/background labels
  m <- merge(z_ssRNA,tblist_ssRNA, by.x="Duplex_ID_gene", by.y="V2")
  if(any(m$V3=="target")){
    if(any(m$V3=="background")){
      # run mann whitney test
      target_stab <- m[m$V3=="target",opt$score]
      background_stab <- m[m$V3=="background",opt$score]
      mw <-wilcox.test(target_stab, background_stab, paired = F, alternative = "greater")
      # rbind new row
      outdf <- rbind(c(ssRNA, mw$p.value, median(target_stab)-median(background_stab)),outdf)
    }
  }
}

colnames(outdf)<-outdf_cols

# adjust pvalues
outdf$padj <- p.adjust(outdf$pval, method = "fdr")

# save tabular results as stdout
write.table(x = outdf[order(outdf$pval), ], file = "", quote = F, sep = "\t", row.names = F)

