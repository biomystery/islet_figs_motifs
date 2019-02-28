## input: 1. summarizedExperiement(SE) obj for chromVAR 2. Jaspar matrix 
## output: 1. motif x cell (z score) 2. plot: ranked           
source("./libs.R")
##------------------------------------------------------------
## inputs
##------------------------------------------------------------

input.chromVar.res.list <- readRDS(file = "../dat/output.jaspar.dev.res.Rdata")
input.chromVar.jaspar.z <- data.table(assays(input.chromVar.res.list$dev)$z,keep.rownames = T)
input.umap.res <- fread('../dat/Islet_123.MNN_corrected.UMAP.txt',header = T)
input.chromVar.jaspar.var <- fread("../dat/output.jaspar.var.res.csv")

# filter unkonwn
input.umap.res <- input.umap.res %>% 
  separate(cluster,into = c("cell_type_overall","subtype"),remove=F)
input.umap.res[is.na(input.umap.res)]<-0



# aggregate data  --------------------------------------------------------------
# melt
input.chromVar.jaspar.z.agg <- melt(input.chromVar.jaspar.z,
                                    id="rn",variable.name = "barcodes",value.name = "zval")

# add celltype
input.chromVar.jaspar.z.agg <- merge(input.chromVar.jaspar.z.agg,input.umap.res)

# average over cell type
output.chromvar.jaspar.z.avg_by_ct <- input.chromVar.jaspar.z.agg[,.(zval_avg=mean(zval)),
                                                                  by=.(rn,cell_type_overall)]
output.chromvar.jaspar.z.avg_by_ct <- output.chromvar.jaspar.z.avg_by_ct[, zval_avg:=(zval_avg-min(zval_avg))/(max(zval_avg)-min(zval_avg)),
                                                                         by=.(rn)]
output.chromvar.jaspar.z.avg_by_ct <- dcast(output.chromvar.jaspar.z.avg_by_ct,rn~cell_type_overall,value.var = "zval_avg")
setDF(output.chromvar.jaspar.z.avg_by_ct)
output.chromvar.jaspar.z.avg_by_ct<- output.chromvar.jaspar.z.avg_by_ct %>%
  separate(rn,into = c("id","name"),sep = "_")%>%
  select(-one_of("id"))%>%
  column_to_rownames("name")



# heatmap -----------------------------------------------------------------

select.motifs <- (input.chromVar.jaspar.var%>%filter(variability>1.2))$name

if(F){
  pheatmap(output.chromvar.jaspar.z.avg_by_ct[select.motifs,c(1:3,6,4:5,7:9)],scale = "none",
           cluster_cols = F, border_color = NA,
           #clustering_method = "ward.D2",
           show_rownames = T,fontsize_row = 6,
           color = cols.hm.avg.tf(30))
  
}

if(T){
  callback = function(hc, mat){
    #sv <- apply(mat, 1, which.max)
    dend =rev(as.dendrogram(hc))
    as.hclust(dend)
  }
  
png(filename = "overall_avg_motif_all_cell_type.png",height = 7.5,width = 3,units = 'in',res = 300)
pdf(file = "overall_avg_motif_all_cell_type.pdf",height = 7.5,width = 3)
  p <- pheatmap(output.chromvar.jaspar.z.avg_by_ct[select.motifs,c(1:3,6,5,8,7,9,4)],scale = "none",
           cluster_cols = T,clustering_callback = callback,
           clustering_method = "ward.D2",treeheight_row = 0,
           show_rownames = T,fontsize_row = 5,border_color = NA,
           color = cols.hm.avg.tf(30))
  print(p)
dev.off()




pd <- output.chromvar.jaspar.z.avg_by_ct[select.motifs,c(1:3,6,5,8,7,9,4)]
pd <- pd[p$tree_row$order,]
require(ggdendro)
hc <- hclust(dist(t(pd)))
dd <- rev(as.dendrogram(hc))
dd.reorder <- reorder(dd,1:9)
hc <- as.hclust(rev(dd.reorder))

hc$order
if(T){
  #dd.reorder <- reorder(dd,c(1,3,2,4,5,6,7,8,9))
  plot(dd.reorder)
  plot(dd)
}

pd <- pd[,hc$order]
p <- pheatmap(pd,scale = "none",cluster_rows = F,
              cluster_cols = F,clustering_callback = callback,
              clustering_method = "ward.D2",treeheight_row = 0,
              show_rownames = T,fontsize_row = 5,border_color = NA,
              color = cols.hm.avg.tf(30))
ggdendrogram(hc,  size = 2)+coord_flip()
saveRDS(list(pd=pd,hc=hc),'../dat/figdata/Fig1E.Rdata')

motif.list <- c("PDX1","RFX2","FOXO4","HNF1A","IRF2","ETS1","MEF2B","MAFK")
row_labs <- sapply(select.motifs, function(x) ifelse(x %in% motif.list,x,""))
png(filename = "overall_avg_motif_all_cell_type_hightlight.png",height = 7.5,width = 3,units = 'in',res = 300)
pdf(file = "overall_avg_motif_all_cell_type_hightlight.pdf",height = 7.5,width = 3)
pheatmap(output.chromvar.jaspar.z.avg_by_ct[select.motifs,c(1:3,6,5,8,7,9,4)],scale = "none",
         cluster_cols = F,clustering_callback = callback,
         clustering_method = "ward.D2",
         labels_row = row_labs,
         show_rownames = T,fontsize_row = 6,border_color = NA,
         color = cols.hm.avg.tf(30))
dev.off()
}




