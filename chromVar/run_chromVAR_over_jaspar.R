## input: 1. summarizedExperiement(SE) obj for chromVAR 2. Jaspar matrix 
## output: 1. motif x cell (z score) 2. plot: ranked           


require(chromVAR)
require(BiocParallel)
register(MulticoreParam(10, progressbar = TRUE))
require(SummarizedExperiment)
require(GenomicRanges)
require(Matrix)
require(BSgenome.Hsapiens.UCSC.hg19)
require(motifmatchr)
require(tidyverse)
#
set.seed(1098)
##------------------------------------------------------------
## inputs
##------------------------------------------------------------
iutput.SE.for.chromVar.filtered <- readRDS("output.SE.for.chromVar.filtered.Rdata")
input.motifs <- getJasparMotifs()

system.time(motif_ix <- matchMotifs(input.motifs, iutput.SE.for.chromVar.filtered, genome = BSgenome.Hsapiens.UCSC.hg19))
dev <- computeDeviations(object = iutput.SE.for.chromVar.filtered, 
                         annotations = motif_ix)
tmp <- colData(dev)
tmp <- tmp%>% 
  as.data.frame()%>%
  separate(cell_type,into = c("cell_type_overall","subtype"),remove=F)
tmp[is.na(tmp)] <- 0
colData(dev) <- DataFrame(tmp)
saveRDS(list(dev=dev,motif_ix=motif_ix),"output.jaspar.dev.res.Rdata")

## variability 
variability <- computeVariability(dev)
plotVariability(variability, use_plotly = FALSE)

variability<- variability%>%arrange(desc(variability))

write.csv(variability,file = "output.jaspar.var.res.csv")


# clustering --------------------------------------------------------------

sample_cor <- getSampleCorrelation(dev)
require(pheatmap)
png(filename = "res.png")
pheatmap(as.dist(sample_cor), 
         annotation_row = colData(dev), 
         clustering_distance_rows = as.dist(1-sample_cor), 
         clustering_distance_cols = as.dist(1-sample_cor))
dev.off()


# tsne --------------------------------------------------------------------
tsne_results <- deviationsTsne(dev, threshold = 1.2, perplexity = 10, 
                               shiny = FALSE)

tsne_plots <- plotDeviationsTsne(dev, tsne_results, annotation = "PDX1", 
                                 sample_column = "cell_type_overall", shiny = FALSE)

require(ggpubr)
ggarrange(tsne_plots[[1]],tsne_plots[[2]])

ggsave(filename = "output.res.jaspar.tsne2.pdf",tsne_plots[[2]])
