## input: 1. summarizedExperiement(SE) obj for chromVAR 2. Jaspar matrix 
## output: 1. motif x cell (z score) 2. plot: ranked           
source("./libs.R")
##------------------------------------------------------------
## inputs
##------------------------------------------------------------

ttest.res <- fread("ttest.res.csv")

dmotifs.list <- sapply(c("alpha","beta","delta"), function(x) subset(ttest.res,celltype==x & selected)$motif)
tmp <- gplots::venn(dmotifs.list)
dmotifs.list.inter <- attr(tmp,"intersections")

input.chromVar.res.list <- readRDS(file = "./dat/output.jaspar.dev.res.Rdata")
input.chromVar.jaspar.z <- assays(input.chromVar.res.list$dev)$z

select.celltypes <- c("alpha_1" , "alpha_2" , "beta_1",   "beta_2" ,  "delta_1",  "delta_2" )


if(F){
  input.umap.res <- read.table('./Islet_123.MNN_corrected.UMAP.txt',header = T)
  input.umap.res <- input.umap.res %>% 
    filter(cluster!="islet_unknown")%>%
    separate(cluster,into = c("cell_type_overall","subtype"),remove=F)
  input.umap.res[is.na(input.umap.res)]<-0
  saveRDS(input.umap.res,'res.umap.Rdata')
  
}
input.umap.res <- readRDS("./dat/res.umap.Rdata")


# add Z score of motifs  -------------------------------------------------
input.umap.res.sub <- subset(input.umap.res,cluster %in% select.celltypes)

p.default.cluster.sub <-  ggplot(input.umap.res.sub,aes(UMAP1,UMAP2)) + 
  geom_point(aes(colour=cluster),size=1,alpha=0.6) + 
  scale_color_manual(values = cols.subcluster) +
  theme_light()

png(filename ="umap_cell_subtype.png",width = 4,height = 3,res = 300,units = 'in')
print(p.default.cluster.sub)
dev.off()



motif.list <- c("SNAI2","NFYA","SPIC","POU6F2","ATF4","MSC","TCF4","CTCF","TEAD3")

if(T){
  p.default.cluster.motifs <- lapply(motif.list, function(x) fun.plot.project.motif(motif = x,umap.res = input.umap.res.sub))
  png(filename ="project_overall_avg_motif_selected_cell_subtype.png",width = 10,height = 8,res = 300,units = 'in')
  ggarrange(plotlist=p.default.cluster.motifs,ncol = 3,nrow = 3)
  dev.off()
}


# project_overall_avg_motif_all_cell_type_scaled.png
if(T){
  input.chromVar.jaspar.z.scale <- input.chromVar.jaspar.z
  input.chromVar.jaspar.z.scale[input.chromVar.jaspar.z.scale>5] <- 5
  input.chromVar.jaspar.z.scale[input.chromVar.jaspar.z.scale< -5] <- -5
  p.default.cluster.motifs.2 <- lapply(motif.list, 
                                     function(x) fun.plot.project.motif(motif = x,
                                                                        umap.res = input.umap.res.sub,
                                                                        input.chromVar.z = input.chromVar.jaspar.z.scale))
  
  png(filename ="project_overall_avg_motif_selected_cell_subtype_scaled.png",width = 10,height = 8,res = 300,units = 'in')
  ggarrange(plotlist=p.default.cluster.motifs.2,ncol = 3,nrow = 3)
  dev.off()

}





