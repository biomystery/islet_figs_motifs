## input: 1. peaks (255k) 2. filtered reads (128M reads) 3. cell-celltype dict (15k cells)
## output: 1. peaks x cells - sparse counts


## https://stackoverflow.com/questions/27574775/is-it-possible-to-use-the-r-data-table-function-foverlaps-to-find-the-intersecti
## https://stackoverflow.com/questions/48141991/foverlaps-and-within-in-data-table
## algorithm: use data.table::foverlap + dcast function


##------------------------------------------------------------
## inputs
##------------------------------------------------------------


require(chromVAR)
require(BiocParallel)
require(GenomicRanges)
require(Matrix)
register(MulticoreParam(4, progressbar = TRUE))
require(Matrix.utils)
require(tidyverse)
##------------------------------------------------------------
## inputs
##------------------------------------------------------------
input.peaks <- getPeaks('./Islet_123.combined.merged.bed',sort_peaks = TRUE)
input.peaks <- resize(input.peaks,width=500,fix="center")
input.peaks <- as.data.table(input.peaks)[,.(seqnames,start,end)]
setkey(input.peaks,seqnames,start,end)

input.filtered.reads <-  fread('./output.filtered.reads.bed')
input.filtered.reads <- input.filtered.reads[,.(V1,V2,V3,V4)]
input.filtered.reads <- input.filtered.reads[,.(seqnames=V1,start=V2,end=V3,cell=V4)]


input.cells.final <- fread('./output.cell.celltypes.txt')

##------------------------------------------------------------
## main
##------------------------------------------------------------
## preprocess to save memory 
system.time(output.peaks.by.cells.count <- foverlaps(input.filtered.reads,input.peaks,nomatch=0L))
output.peaks.by.cells.count <- output.peaks.by.cells.count[,.(.N),by=.(seqnames,start,end,cell)]
output.peaks.by.cells.count <- output.peaks.by.cells.count[,cell:=as.factor(cell)]
output.peaks.by.cells.count <- output.peaks.by.cells.count %>% unite(col = peak,1:3)
output.peaks.by.cells.count <- output.peaks.by.cells.count[,peak:=as.factor(peak)]
fwrite(output.peaks.by.cells.count,"output.peaks.by.cells.count.csv")



#system.time(tmp <- output.peaks.by.cells.count%>%
#              spread(key=cell,value = N))
## Error: vector memory exhausted (limit reached?)
## Timing stopped at: 148.3 7.347 156.2

#system.time(tmp <- dcast(output.peaks.by.cells.count,seqnames + start + end ~ rowid(cell),
#                         fun.aggregate = sum,value.var = "N"))
## CJ() would result in 5534801450 rows which exceeds .Machine$integer.max == 2147483647

system.time(output.peaks.by.cells.count.sparse <- dMcast(output.peaks.by.cells.count,peak ~ cell,
                          fun.aggregate="I",value.var = "N"))

colnames(output.peaks.by.cells.count.sparse) <- sub("cell","",colnames(output.peaks.by.cells.count.sparse))

#writeMM(output.peaks.by.cells.count.sparse,"output.peaks.by.cells.count.sparse.txt")
saveRDS(output.peaks.by.cells.count.sparse,"output.peaks.by.cells.count.sparse.Rdata")
## create sparse matrix 



