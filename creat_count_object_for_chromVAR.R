## input: 1. peaks x cells - sparse counts 2. cell,read_depth 3. cell_type assignement
## output: 1. summarizedExperiement(SE) obj for chromVAR
##              


require(data.table)
require(chromVAR)
require(BiocParallel)
register(MulticoreParam(4, progressbar = TRUE))
require(SummarizedExperiment)
require(GenomicRanges)
require(Matrix)
require(BSgenome.Hsapiens.UCSC.hg19)
require(tidyverse)
#

##------------------------------------------------------------
## inputs
##------------------------------------------------------------

input.peaks.by.cells.count.sparse<-readRDS("output.peaks.by.cells.count.sparse.Rdata")
input.reads.depth.per.cell<- fread("./output.reads.depth.per.cell.csv",col.names = c("cell","depth"))
input.cell.celltypes<- read.table('./output.cell.celltypes.txt',header = T,row.names = 1)

# construct SE object  ----------------------------------------------------

output.SE.for.chromVar <- SummarizedExperiment(assays=list(counts=input.peaks.by.cells.count.sparse),
                    colData=input.reads.depth.per.cell)

rowRanges.gr <- do.call(rbind,sapply(rownames(output.SE.for.chromVar),strsplit,split="_"))
rowRanges.gr <- as.data.frame(rowRanges.gr)
colnames(rowRanges.gr)<- c("seqname","start","end")
rowRanges.gr <- makeGRangesFromDataFrame(rowRanges.gr)
rowRanges(output.SE.for.chromVar) <- rowRanges.gr


## filter sample 
filtering_plot <- filterSamplesPlot(output.SE.for.chromVar, min_depth = 1500, 
                                    min_in_peaks = 0.15, use_plotly = FALSE)
filtering_plot

output.SE.for.chromVar.filtered <- filterSamples(output.SE.for.chromVar, min_depth = 1500, 
                                        min_in_peaks = 0.15, shiny = FALSE)

## filter peaks 
output.SE.for.chromVar.filtered <- sort(output.SE.for.chromVar.filtered)
output.SE.for.chromVar.filtered <- filterPeaks(output.SE.for.chromVar.filtered)


## add GC 
output.SE.for.chromVar.filtered<- addGCBias(output.SE.for.chromVar.filtered, genome = BSgenome.Hsapiens.UCSC.hg19)

## add cell_type assignment 
cd <- colData(output.SE.for.chromVar.filtered)
cd$cell_type <- input.cell.celltypes[rownames(cd),1]
colData(output.SE.for.chromVar.filtered) <- cd

saveRDS(output.SE.for.chromVar.filtered,"output.SE.for.chromVar.filtered.Rdata")
