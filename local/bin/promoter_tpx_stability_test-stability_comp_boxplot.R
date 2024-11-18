#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(ggpubr)))
suppressMessages(suppressWarnings(library(rstatix)))
suppressMessages(suppressWarnings(library(dplyr)))


# check infile ----
check_infile <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(paste(file_path, "does not exist."))
  }
  file_info <- file.info(file_path)
  if (file_info$size == 0) {
    stop(paste(file_path, "is empty."))
  } else {
    message(paste("--- Reading:", file_path))
  }
}


# options ---
option_list <- list( 
  make_option(c("-t", "--tpxlist"), action="store", default=NULL,
              help="Matrix of ssRNA tpx calculated with 3plex on target and background regions."),
  make_option(c("-l", "--tblist"), action="store", default=NULL,
              help="Target/background list path [ gene_name, class ]."),
  make_option(c("-s", "--score"), action="store", default="Stability_best",
              help="TPX score used to rank target/background regions. Available: Stability_best, Stability_tot, Score_best, Stability_norm. [ Default: Stability_norm ]."),
  make_option(c("-d", "--directory"), action ="store", default=".",
              help="Output directory [ Default: . ]."),
  make_option(c("-f","--format"), action="store", default="pdf",
              help="Plots format: pdf or png. [Default: pdf]")
)
opt <- parse_args(OptionParser(option_list = option_list))



# create output directory
if(!dir.exists(opt$directory)) dir.create(opt$directory, recursive = T, showWarnings = FALSE)

# read 3plex output 
check_infile(opt$tpxlist)
tpx <- read.delim(opt$tpxlist, header=T)

# read target and background list
check_infile(opt$tblist)
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
boxplot_stats <- m %>%
  group_by(class) %>%
  summarise(
    min = min(!!sym(opt$score)),
    Q1 = quantile(!!sym(opt$score), 0.25),
    median = median(!!sym(opt$score)),
    Q3 = quantile(!!sym(opt$score), 0.75),
    max = max(!!sym(opt$score)),
    IQR = IQR(!!sym(opt$score))
  ) %>%
  mutate(
    lower_whisker = pmax(Q1 - 1.5 * IQR, min),
    upper_whisker = pmin(Q3 + 1.5 * IQR, max)
  )

# Find the maximum upper whisker value
max_whisker <- max(boxplot_stats$upper_whisker)

# Add a small offset to the maximum whisker value for the comparison bar position
y_position <- max_whisker + 0.1 * (max(m[[opt$score]]) - min(m[[opt$score]]))

bxp <- ggboxplot(m, x = "class", y = opt$score, 
                 fill = "class", palette = "d3", outlier.shape = NA) +
  stat_pvalue_manual(stat.test, label = "p.signif", y.position = y_position) +
  theme_bw() + theme(legend.position = "none") +
  #geom_jitter(width = .1, height = NULL, size =.5, alpha=.2) +
  xlab("") + ggtitle(ssRNA) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(min(m[[opt$score]]), y_position + 0.1 * (max(m[[opt$score]]) - min(m[[opt$score]])))) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate(geom = "text", label = paste0("Wilcoxon's p-value: ", stat.test$p), x = -Inf, y = Inf, hjust = -.05, vjust = 2, size=3)

# save image
outfile <- paste0(opt$directory, "/stability_comp_boxplot.pdf")
message(paste0("--- Saving to: ",outfile))
pdf(file = outfile, paper="a4", width = 3, height = 3.5)
bxp
rm <- dev.off()
