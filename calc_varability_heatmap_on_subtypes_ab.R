## input: 1. summarizedExperiement(SE) obj for chromVAR 2. Jaspar matrix 
## output: 1. motif x cell (z score) 2. plot: ranked           
source("./libs.R")

# combine -----------------------------------------------------------------

ttest.res <- fread("./dat/ttest.res.csv")

dmotifs.list <- sapply(c("alpha","beta"), function(x) subset(ttest.res,celltype==x & selected)$motif)
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
  input.chromVar.jaspar.z.agg[cell_type_overall %in% c("alpha","beta") & subtype!=3,
                              .(zval_avg=mean(zval)),by=.(rn,cluster)]%>%
  group_by(rn)%>%
  mutate(zval_avg = (zval_avg-min(zval_avg))/(max(zval_avg)-min(zval_avg)))%>%
  separate(rn,into = c("id","name"),sep = "_")%>%
  select(-one_of("id"))%>%
  spread(key = cluster,value = zval_avg)%>%
  as.data.frame()%>%
  column_to_rownames("name")

output.chromvar.jaspar.z.avg_by_subct.noscale<-  
  input.chromVar.jaspar.z.agg[cell_type_overall %in% c("alpha","beta") & subtype!=3,
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




if(T){
  ht <-Heatmap(output.chromvar.jaspar.z.avg_by_subct[unique(unlist(dmotifs.list.inter)),],
               col=cols.hm.avg.tf(30),cluster_columns = F,
               split = factor(splt,levels = c("alpha","beta","alpha:beta")),
               row_names_gp=gpar(fontsize=5),name = "ht"
  )
  
  
  pdf("hmp.subcelltype.pdf",height = 10,width = 3)
  print(ht)
  decorate_heatmap_body("ht", {
    grid.lines(c(1/2,1/2), c(-4,1), gp = gpar(lty = 2, lwd = 2,col="white"))
  })
  dev.off()
  system("open hmp.subcelltype.pdf")
}

if(T){
  ht.nolab <-Heatmap(output.chromvar.jaspar.z.avg_by_subct[unique(unlist(dmotifs.list.inter)),],
                     col=cols.hm.avg.tf(30),cluster_columns = F,show_column_names = F,show_row_names = F,
                     split = factor(splt,levels = c("alpha","beta","alpha:beta")),
                     heatmap_legend_param = list(title = NULL,legend_width = unit(5, "cm"),
                                                 labels_gp = gpar(fontsize = 6)),
                     name = "ht.nolab",show_heatmap_legend = F,show_row_dend = F,combined_name_fun = NULL
  )
  #pdf("hmp.subcelltype.nolab.pdf",height = 4,width = 1.78)
  png(filename = "./figs/fig2/hmp.subcelltype.nolab.png",height = 4,width = 1.78,units = 'in',res = 300)
  print(ht.nolab)
  decorate_heatmap_body("ht.nolab", {
    grid.lines(c(1/2,1/2), c(-4,1), gp = gpar(lty = 2, lwd = 2,col="white"))
  })
  dev.off()
  system("open ./figs/fig2/hmp.subcelltype.nolab.png")
}

if(T){
  ht.nosplit <-Heatmap(output.chromvar.jaspar.z.avg_by_subct[unique(unlist(dmotifs.list.inter)),],
                       col=cols.hm.avg.tf(30),cluster_columns = F,
                       row_names_gp=gpar(fontsize=5),name = "ht.nosplit",
                       heatmap_legend_param = list(title = NULL,legend_height= unit(4, "cm"),
                                                   labels_gp = gpar(fontsize = 5))
  )
  
  
  pdf("./figs/fig2/hmp.subcelltype.nosplit.pdf",height = 10,width = 3)
  print(ht.nosplit)
  decorate_heatmap_body("ht.nosplit", {
    grid.lines(c(1/2,1/2), c(-4,1), gp = gpar(lty = 2, lwd = 2,col="white"))
  })
  dev.off()
  system("open ./figs/fig2/hmp.subcelltype.nosplit.pdf")
}

if(T){
  
  ht.nosplit.nolab <-Heatmap(output.chromvar.jaspar.z.avg_by_subct[unique(unlist(dmotifs.list.inter)),],
                             col=cols.hm.avg.tf(30),cluster_columns = F,
                             show_column_names = F,show_row_names = F,
                             name = "ht.nosplit.nolab",show_heatmap_legend = F,
                             show_row_dend = F,combined_name_fun = NULL)
 
  #pdf("hmp.subcelltype.nosplit.nolab.pdf",height = 4,width = 1.78)
  png(filename = "hmp.subcelltype.nosplit.nolab.png",height = 4,width = 1.78,units = 'in',res = 300)
  print(ht.nosplit.nolab)
  decorate_heatmap_body("ht.nosplit.nolab", {
    grid.lines(c(1/2,1/2), c(-4,1), gp = gpar(lty = 2, lwd = 2,col="white"))
  })
  dev.off()
  system("open hmp.subcelltype.nosplit.nolab.png")
}



if(T){
  png(filename = "./figs/overall_avg_motif_cell_subtype.png",height = 7.5,width = 3,units = 'in',res = 300)
  pheatmap(output.chromvar.jaspar.z.avg_by_subct[unique(unlist(dmotifs.list.inter)),],scale = "none",
           cluster_cols = F,
           show_rownames = T,fontsize_row = 5,border_color = NA,
           color = cols.hm.avg.tf(30))
  dev.off()
  system("open ./figs/overall_avg_motif_cell_subtype.png")
}


# correlations ------------------------------------------------------------
pd.corr<- cor(output.chromvar.jaspar.z.avg_by_subct[unique(unlist(dmotifs.list.inter)),],
              method = "spearman"
              )

pd.corr<- cor(output.chromvar.jaspar.z.avg_by_subct)

Heatmap(pd.corr, name = "foo",
        cluster_rows = F,cluster_columns = F,
        cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%.1f", pd.corr[i, j]), x, y, gp = gpar(fontsize = 10))
})



# plot violin -------------------------------------------------------------
celltypes <- c("alpha","beta")
output.jaspar.z<- input.chromVar.jaspar.z.agg%>%
  filter(cell_type_overall%in% celltypes)%>%
  filter(subtype!="3")%>%
  select(one_of(c("rn","zval","cluster","cell_type_overall","subtype")))%>%
  separate(rn,into = c("Jaspar.id","Motif.name"),sep = "_")
  
  
for(m in c("RFX3","FOSL1")){
  p<-ggviolin(output.jaspar.z%>%
                filter(Motif.name==m),
              x = "subtype",y = "zval",add = "boxplot",fill="cluster"
  )
  
  p<- p + facet_wrap(~ cell_type_overall)+
    scale_fill_manual(values = cols.celltype)+
    theme_light()+
    coord_cartesian(expand = T)
  fn<- paste0("figs/fig2/",m,"_violion.pdf")
  ggsave(p,filename = fn,width = 1.5,height = 1,scale = 5)
  system(paste0("open ",fn))
  fn <-  paste0("figs/fig2/",m,"_violion_nolab.pdf")
  ggsave(p+theme(text = element_blank(),legend.position = "none"),
         filename = fn,width = 1.5,height = 1,scale = 5)
  system(paste0("open ",fn))
  
}
