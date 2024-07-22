#!/usr/bin/env Rscript

suppressMessages(suppressWarnings((library(fgsea))))
suppressMessages(suppressWarnings((library(ggplot2))))
suppressMessages(suppressWarnings((library(optparse))))

# Functions ----
read.gmt <- function (file) 
{
  if (!grepl("\\.gmt$", file)[1]) {
    stop("Pathway information must be a .gmt file! Exit")
  }
  geneSetDB = readLines(file)
  geneSetDB = strsplit(geneSetDB, "\t")
  names(geneSetDB) = sapply(geneSetDB, "[", 1)
  geneSetDB = lapply(geneSetDB, "[", -1:-2)
  geneSetDB = lapply(geneSetDB, function(x) {
    x[which(x != "")]
  })
  return(geneSetDB)
}

# Options ----
option_list <- list( 
  make_option(c("-g", "--gmtfile"), action="store", default=NULL,
              help="Download the .gmt file"),
  make_option(c("-r", "--rnkfile"), action="store", default=NULL, 
              help="Download the .rnk file"),
  make_option(c("-n", "--nperm"), type="integer", default=1000, 
              help="Number of permutations to compute p-value",
              metavar="number"),
  make_option(c("-b", "--biggest"), type="integer", default=500, 
              help="Maximum genes in gene set",
              metavar="number"),
  make_option(c("-s", "--smallest"), type="integer", default=15, 
              help="Minimum genes in gene set",
              metavar="number"),
  make_option(c("-d", "--directory"), action ="store", default=".",
              help="Directory containing fgsea results, enrichment plots and fgsea Leading Edges Details"),
  make_option(c("-l", "--leadingedge"), action ="store", default=".",
              help="Subdirectory containing fgsea Leading Edges Details"),
  make_option(c("-f", "--filetype"), action ="store", default="csv",
              help="Determines file type for fgsea results (csv by default) and fgsea Leading Edges Details (tsv by default)"),
  make_option("--gseaParam", type="integer", default=0,
              help="GSEA parameter value, all gene-level stats are raised to the power of gseaParam before calculation of ES."),
  make_option(c("-j","--cores"), type="integer", default=1,
              help="Parallelization.")
)
opt <- parse_args(OptionParser(option_list = option_list))


# 0. Check options and read input files ----
# Check gseaParam
if(!opt$gseaParam %in% c(0,1))
  stop("gseaParam must be 0 or 1! Exit")
# Read pathway file
gmt <- read.gmt(opt$gmtfile)
# Read rank file
rnk <- read.table(opt$rnkfile, header = F, col.names = c("GeneID", "value"))
rnk_vector <- unlist(rnk$value)
names(rnk_vector) <- rnk$GeneID
rnk_vector <- sort(rnk_vector, decreasing = TRUE)
# Create fGSEA an Leading_Edge directory
if(!dir.exists(opt$directory)) dir.create(opt$directory, recursive = T, showWarnings = FALSE)
if(!dir.exists(opt$leadingedge)) dir.create(opt$leadingedge, recursive = T, showWarnings = FALSE)


# 1. Run fGSEA ----
# Add condition: if(simple or multilevel)
fgseaRes <- suppressWarnings(fgsea::fgseaSimple(
  pathways = gmt,
	stats    = rnk_vector,
	nperm    = opt$nperm,
	minSize  = opt$smallest,
	maxSize  = opt$biggest,
	gseaParam = opt$gseaParam, #https://www.biostars.org/p/263074/
	nproc = opt$cores
  ))
# Check fGSEA output
if(!is.null(fgseaRes)){
  message(" -- fgseaSimple done")
} else {
  stop(" -- Error in fgseaSimple! Exit")
}
# save fgseaRes
outfile_tab <- paste0(opt$directory,"/fgseaRes.tsv")
message(" -- saving to: ", outfile_tab)
write.table(subset(fgseaRes, select = -c(leadingEdge)), outfile_tab, col.names = T, quote = F, row.names = F, sep="\t")


# 2. Plot fgsea ----
# fGSEA enrich plot
# # https://www.genekitr.fun/plot-gsea-1
fgseaPlot <- fgsea::plotEnrichment(
  pathway = gmt$genes_of_interest,
  stats = rnk_vector,
  gseaParam = opt$gseaParam
  )
fgseaPlot <- fgseaPlot + 
  xlab(paste0("Genes in descending order of ", basename(opt$directory))) +
  labs(title = "GeneSet Enrichment Analysis",
       subtitle = paste0("NES: ", fgseaRes[which(fgseaRes$pathway=="genes_of_interest"),"NES"], "\nFDR: ", fgseaRes[which(fgseaRes$pathway=="genes_of_interest"),"padj"]))

# save plot fgsea
outfile <- paste0(opt$directory,"/enrichment_plot.pdf")
message(" -- saving to: ", outfile)
pdf(file = outfile, paper = "a4", w = 4, h = 4)
print(fgseaPlot)
dev.off()

# 3. Create LeadingEdge table ----
gene_set <- as.data.frame(gmt$genes_of_interest)
colnames(gene_set) <- "GeneID"
leadingEdge_tab <- as.data.frame(fgseaRes$leadingEdge[[1]])
colnames(leadingEdge_tab) <- "GeneID"
leadingEdge_tab$hit <- "Yes"
leadingEdge_tab <- merge(gene_set, leadingEdge_tab, all.x=T, by="GeneID")
leadingEdge_tab[is.na(leadingEdge_tab)] <- "No"
leadingEdge_tab <- merge(leadingEdge_tab, rnk, by="GeneID")
colnames(leadingEdge_tab)[3]<-"ranking_score"
leadingEdge_tab <- leadingEdge_tab[order(leadingEdge_tab$ranking_score),]
# save LeadingEdge
outfile <- paste0(opt$directory,"/leading_edge.tsv")
message(" -- saving to: ", outfile)
write.table(leadingEdge_tab, outfile, col.names = T, quote = F, row.names = F, sep="\t")

