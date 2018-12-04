## input: 1. summarizedExperiement(SE) obj for chromVAR 2. Jaspar matrix 
## output: 1. motif x cell (z score) 2. plot: ranked           
source("./libs.R")
##------------------------------------------------------------
## inputs
##------------------------------------------------------------

input.chromVar.res.list <- readRDS(file = "./dat/output.jaspar.dev.res.Rdata")
input.chromVar.jaspar.z <- data.table(assays(input.chromVar.res.list$dev)$z,keep.rownames = T)
input.umap.res <- fread('./dat/Islet_123.MNN_corrected.UMAP.txt',header = T)
input.chromVar.jaspar.var <- fread("./dat/output.jaspar.var.res.csv")

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


# Recalculate var on selected celltypes  ----------------------------------
select.celltypes <- c("alpha_1" , "alpha_2",  "alpha_3" , "beta_1",   "beta_2" ,  "delta_1",  "delta_2" , "gamma" ,   "exocrine")
dev <- input.chromVar.res.list$dev[,as.character(colData(input.chromVar.res.list$dev)$cell_type) %in% select.celltypes]

## variability 
variability <- computeVariability(dev)
png(filename = "output.jaspar.var.res.subtype.png",width = 4,height = 3,res = 300,units = "in")
plotVariability(variability, use_plotly = FALSE)
#abline(h=1.2,col=2,lty=2,plot.new=T)
dev.off()

variability<- variability%>%arrange(desc(variability))

write.csv(variability,file = "output.jaspar.var.res.subtype.csv")
select.motifs.sub <- (variability%>%filter(variability>1.2))$name


# average over subcell type -----------------------------------------------



output.chromvar.jaspar.z.avg_by_subct<-  
  input.chromVar.jaspar.z.agg[cluster %in% select.celltypes,.(zval_avg=mean(zval)),by=.(rn,cluster)]%>%
  group_by(rn)%>%
  mutate(zval_avg = (zval_avg-min(zval_avg))/(max(zval_avg)-min(zval_avg)))%>%
  separate(rn,into = c("id","name"),sep = "_")%>%
  select(-one_of("id"))%>%
  spread(key = cluster,value = zval_avg)%>%
  as.data.frame()%>%
  column_to_rownames("name")
  
  
callback = function(hc, mat){
  sv <- apply(mat, 1, which.max)
  dend = rev(hc)
  as.hclust(dend)
}

png(filename = "overall_avg_motif_cell_subtype.png",height = 7.5,width = 3,units = 'in',res = 300)
pheatmap(output.chromvar.jaspar.z.avg_by_subct[select.motifs,select.celltypes],scale = "none",
         cluster_cols = F,
         show_rownames = T,fontsize_row = 5,border_color = NA,
         color = cols.hm.zval.fun(30))
dev.off()

motif.list <- c("SNAI2","NFYA","SPIC","POU6F2","ATF4","MSC","TCF4","CTCF","TEAD3")
row_labs <- sapply(select.motifs, function(x) ifelse(x %in% motif.list,x,""))
png(filename = "overall_avg_motif_cell_subtype_highlight.png",height = 7.5,width = 3,units = 'in',res = 300)
pheatmap(output.chromvar.jaspar.z.avg_by_subct[select.motifs,select.celltypes],scale = "none",
         cluster_cols = F,
         show_rownames = T,fontsize_row = 5,border_color = NA,
         color = cols.hm.avg.tf(30),
         labels_row = row_labs)
dev.off()


