## Input: 1. the raw reads file
##        2. the final cell and cell_type file.


## Output: 1. a sparse matrix for peaks counts.
##         2. a sumarizeReads

##------------------------------------------------------------
## libraries
##------------------------------------------------------------
require(data.table)

##------------------------------------------------------------
## prepare input files
##------------------------------------------------------------
input.raw.reads <- fread('./Islet_123.tagAlign', header=FALSE)
## grep -v unknow Islet_123.MNN_corrected.UMAP.txt  | awk '{print $1,$4}' > Islet_123.cell.celltypes.txt
input.cells.final <-  read.table('./Islet_123.cell.celltypes.txt',header=T,stringsAsFactor=F)

## filter reads
output.filtered.reads <- input.raw.reads[V4 %in% input.cells.final$barcodes]

fwrite(output.filtered.reads,"output.filtered.reads.bed",sep='\t',col.names=F)






