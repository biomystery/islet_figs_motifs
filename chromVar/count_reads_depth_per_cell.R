## input:1.  filtered reads (o128M reads) 3. cell-celltype dict (15k cells)
## output: 1.# reads per cells

##------------------------------------------------------------
## inputs
##------------------------------------------------------------

require(data.table)

##------------------------------------------------------------
## inputs
##------------------------------------------------------------
output.filtered.reads <- fread("./output.filtered.reads.bed")

setkey(output.filtered.reads,V4)

system.time(
    output.reads.depth.per.cell <- output.filtered.reads[,.(.N),by=.(V4)]
)

fwrite(output.reads.depth.per.cell,"output.reads.depth.per.cell.csv",col.names=F)



