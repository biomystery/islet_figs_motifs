## input: 1. summarizedExperiement(SE) obj for chromVAR 2. Jaspar matrix 
## output: 1. motif x cell (z score) 2. plot: ranked           
source("./libs.R")

# combine -----------------------------------------------------------------

ttest.res <- fread("ttest.res.csv")

dmotifs.list <- sapply(c("alpha","beta","delta"), function(x) subset(ttest.res,celltype==x & selected)$motif)
tmp <- gplots::venn(dmotifs.list)
dmotifs.list.inter <- attr(tmp,"intersections")


# Plot_hm -----------------------------------------------------------------\
input.chromVar.res.list <- readRDS(file = "./dat/output.jaspar.dev.res.Rdata")
input.chromVar.jaspar.z <- data.table(assays(input.chromVar.res.list$dev)$z,keep.rownames = T)
rm(input.chromVar.res.list)
input.umap.res <- fread('./dat/Islet_123.MNN_corrected.UMAP.txt',header = T)


# filter unkonwn
input.umap.res <- input.umap.res %>% 
  separate(cluster,into = c("cell_type_overall","subtype"),remove=F)
input.umap.res[is.na(input.umap.res)]<-0

input.chromVar.jaspar.z.agg <- melt(input.chromVar.jaspar.z,
                                    id="rn",variable.name = "barcodes",value.name = "zval")

# add celltype
input.chromVar.jaspar.z.agg <- merge(input.chromVar.jaspar.z.agg,input.umap.res)

output.chromvar.jaspar.z.avg_by_subct<-  
  input.chromVar.jaspar.z.agg[cell_type_overall %in% c("alpha","beta","delta") & subtype!=3,
                              .(zval_avg=mean(zval)),by=.(rn,cluster)]%>%
  group_by(rn)%>%
  mutate(zval_avg = (zval_avg-min(zval_avg))/(max(zval_avg)-min(zval_avg)))%>%
  separate(rn,into = c("id","name"),sep = "_")%>%
  select(-one_of("id"))%>%
  spread(key = cluster,value = zval_avg)%>%
  as.data.frame()%>%
  column_to_rownames("name")

output.chromvar.jaspar.z.avg_by_subct.noscale<-  
  input.chromVar.jaspar.z.agg[cell_type_overall %in% c("alpha","beta","delta") & subtype!=3,
                              .(zval_avg=mean(zval)),by=.(rn,cluster)]%>%
  group_by(rn)%>%
  separate(rn,into = c("id","name"),sep = "_")%>%
  select(-one_of("id"))%>%
  spread(key = cluster,value = zval_avg)%>%
  as.data.frame()%>%
  column_to_rownames("name")

splt <- unlist(sapply(names(dmotifs.list.inter), function(x){
  tmp <- rep(x,length(dmotifs.list.inter[[x]]))
  names(tmp) <-dmotifs.list.inter[[x]]
  tmp
},USE.NAMES = F))


ht <-Heatmap(output.chromvar.jaspar.z.avg_by_subct[unique(unlist(dmotifs.list.inter)),],
        col=cols.hm.avg.tf(30),cluster_columns = F,
        split = factor(splt,levels = c("alpha","beta","delta","alpha:beta","beta:delta","alpha:delta","alpha:beta:delta")),
        row_names_gp=gpar(fontsize=5),name = "ht"
)
if(T){
  pdf("hmp.subcelltype.pdf",height = 10,width = 3)
  print(ht)
  decorate_heatmap_body("ht", {
    grid.lines(c(1/3,1/3), c(-4,1), gp = gpar(lty = 2, lwd = 2,col="white"))
    grid.lines(c(2/3,2/3), c(-4,1), gp = gpar(lty = 2, lwd = 2,col="white"))
  })
  dev.off()
  system("open hmp.subcelltype.pdf")
}


png(filename = "overall_avg_motif_cell_subtype.png",height = 7.5,width = 3,units = 'in',res = 300)
pheatmap(output.chromvar.jaspar.z.avg_by_subct[unique(unlist(dmotifs.list.inter)),],scale = "none",
         cluster_cols = F,
         show_rownames = T,fontsize_row = 5,border_color = NA,
         color = cols.hm.avg.tf(30))
dev.off()
