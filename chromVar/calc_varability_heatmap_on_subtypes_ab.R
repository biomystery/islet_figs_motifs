## input: 1. summarizedExperiement(SE) obj for chromVAR 2. Jaspar matrix 
## output: 1. motif x cell (z score) 2. plot: ranked           
setwd("./chromVar/")
source("./libs.R")

# combine -----------------------------------------------------------------

ttest.res <- fread("../dat/ttest.res.csv")

dmotifs.list <- sapply(c("alpha","beta","delta"), function(x) subset(ttest.res,celltype==x & selected)$motif)
tmp <- gplots::venn(dmotifs.list)
dmotifs.list.inter <- attr(tmp,"intersections")


# Plot_hm -----------------------------------------------------------------\
input.chromVar.res.list <- readRDS(file = "../dat/output.jaspar.dev.res.Rdata")
input.chromVar.jaspar.z <- data.table(assays(input.chromVar.res.list$dev)$z,keep.rownames = T)
rm(input.chromVar.res.list)
input.umap.res <- fread('../dat/Islet_123.MNN_corrected.UMAP.txt',header = T)


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



# fig2A:hm_split ----------------------------------------------------------------
if(T){
  ht <-Heatmap(output.chromvar.jaspar.z.avg_by_subct[unique(unlist(dmotifs.list.inter)),],
               col=cols.hm.avg.tf(30),cluster_columns = F,
               split = factor(splt,levels = c("alpha","beta","delta","alpha:beta","alpha:delta","beta:delta","alpha:beta:delta")),
               row_names_gp=gpar(fontsize=5),name = "ht"
  )
  
  fn <- ".subfig2A_hm_split.pdf"
  pdf(fn,height = 10,width = 3)
  print(ht)
  decorate_heatmap_body("ht", {
    grid.lines(c(1/2,1/2), c(-4,1), gp = gpar(lty = 2, lwd = 2,col="white"))
  })
  dev.off()
  system(paste0("open ",fn))

  ht.nolab <-Heatmap(output.chromvar.jaspar.z.avg_by_subct[unique(unlist(dmotifs.list.inter)),],
                     col=cols.hm.avg.tf(30),cluster_columns = F,show_column_names = F,show_row_names = F,
                     split = factor(splt,levels = c("alpha","beta","alpha:beta")),
                     heatmap_legend_param = list(title = NULL,legend_width = unit(5, "cm"),
                                                 labels_gp = gpar(fontsize = 6)),
                     name = "ht.nolab",show_heatmap_legend = F,show_row_dend = F,combined_name_fun = NULL
  )
  
  fn <- "./subfig2A_hm_split_nolab.pdf"
  pdf(file = fn,height = 4,width = 1.78)
  #png(filename = "./figs/fig2/hmp.subcelltype.nolab.png",height = 4,width = 1.78,units = 'in',res = 300)
  print(ht.nolab)
  decorate_heatmap_body("ht.nolab", {
    grid.lines(c(1/2,1/2), c(-4,1), gp = gpar(lty = 2, lwd = 2,col="white"))
  })
  dev.off()
  system(paste0("open ",fn))
}


# fig2A:hm_nosplit --------------------------------------------------------


if(T){
  ht.nosplit <-Heatmap(output.chromvar.jaspar.z.avg_by_subct[unique(unlist(dmotifs.list.inter)),],
                       col=cols.hm.avg.tf(30),cluster_columns = F,
                       row_names_gp=gpar(fontsize=4),name = "ht.nosplit",
                       heatmap_legend_param = list(title = NULL,legend_height= unit(4, "cm"),
                                                   labels_gp = gpar(fontsize = 5))
  )
  
  fn="/Users/frank/Dropbox (UCSD_Epigenomics)/Islet_snATAC/panel_pdfs/sfigs/fig_s2.motif.hm.with.detla.pdf"
  pdf(fn,height = 10,width = 3)
  print(ht.nosplit)
  #decorate_heatmap_body("ht.nosplit", {
  #  grid.lines(c(1/2,1/2), c(-4,1), gp = gpar(lty = 2, lwd = 2,col="white"))
  #})
  dev.off()
  system(paste0("open ",fn))
  fwrite(data.frame(output.chromvar.jaspar.z.avg_by_subct[unique(unlist(dmotifs.list.inter)),],cate=splt)[row_order(ht.nosplit)[[1]],],
         "/Users/frank/Dropbox (UCSD_Epigenomics)/Islet_snATAC/panel_pdfs/sfigs/fig_s2.motif.hm.with.detla.csv",row.names = T)

  ht.nosplit.nolab <-Heatmap(output.chromvar.jaspar.z.avg_by_subct[unique(unlist(dmotifs.list.inter)),],
                             col=cols.hm.avg.tf(30),cluster_columns = F,
                             show_column_names = F,show_row_names = F,
                             name = "ht.nosplit.nolab",show_heatmap_legend = F,
                             show_row_dend = F,combined_name_fun = NULL)
  fn="./figs/fig2/subfig2A_hm_nolab.pdf"
  pdf(fn,height = 4,width = 1.78)
  #png(filename = "hmp.subcelltype.nosplit.nolab.png",height = 4,width = 1.78,units = 'in',res = 300)
  print(ht.nosplit.nolab)
  decorate_heatmap_body("ht.nosplit.nolab", {
    grid.lines(c(1/2,1/2), c(-4,1), gp = gpar(lty = 2, lwd = 2,col="white"))
  })
  dev.off()
  system(paste0("open ",fn))
  
  saveRDS(output.chromvar.jaspar.z.avg_by_subct[unique(unlist(dmotifs.list.inter)),],file = "../figcode/fig_2B.motif_heatmap_ab.Rds")
  fwrite(output.chromvar.jaspar.z.avg_by_subct[unique(unlist(dmotifs.list.inter)),][row_order(ht.nosplit.nolab)[[1]],],file ="../figcode/fig_2B.motif_heatmap_ab.csv",row.names = T)
  
}


# correlations ------------------------------------------------------------
pd.corr<- cor(output.chromvar.jaspar.z.avg_by_subct[unique(unlist(dmotifs.list.inter)),],
              method = "spearman"
              )

pd.corr<- cor(output.chromvar.jaspar.z.avg_by_subct)

Heatmap(pd.corr, name = "cor",
        cluster_rows = F,cluster_columns = F,
        cell_fun = function(j, i, x, y, width, height, fill) {
  grid.text(sprintf("%.1f", pd.corr[i, j]), x, y, gp = gpar(fontsize = 10))
})



# fig2A:eg violin -------------------------------------------------------------
celltypes <- c("alpha","beta")
output.jaspar.z<- input.chromVar.jaspar.z.agg%>%
  filter(cell_type_overall%in% celltypes)%>%
  filter(subtype!="3")%>%
  select(one_of(c("rn","zval","cluster","cell_type_overall","subtype")))%>%
  separate(rn,into = c("Jaspar.id","Motif.name"),sep = "_")

output.jaspar.z$cluster <- (Hmisc::capitalize(sub("_"," ", output.jaspar.z$cluster)))
fwrite(output.jaspar.z,file = '../figcode/fig_2B.motif_violin_ab.csv') 

select.gene <- c("FOSL1","FOS::JUN","RFX3","TAL1::TCF3","STAT3","NKX6-1","TEAD1","TEAD3")  
for(m in select.gene){
  p<-ggviolin(output.jaspar.z%>%
                filter(Motif.name==m),
              x = "subtype",y = "zval",add = "boxplot",fill="cluster"
  )
  
  p<- p + facet_wrap(~ cell_type_overall)+
    scale_fill_manual(values = cols.celltype)+
    theme_light()+
    coord_cartesian(expand = T)
  fn<- paste0("figs/fig2/subfig2A_",m,"_violion.pdf")
  ggsave(p,filename = fn,width = 1.5,height = 1,scale = 5,useDingbats = F)
  system(paste0("open ",fn))
  fn <-  paste0("figs/fig2/subfig2A_",m,"_violion_nolab.pdf")
  ggsave(p+theme(text = element_blank(),legend.position = "none"),
         filename = fn,width = 1.5,height = 1,scale = 5,useDingbats = F)
  system(paste0("open ",fn))
  
}
