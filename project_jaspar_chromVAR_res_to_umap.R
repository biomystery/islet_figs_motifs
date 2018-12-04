## input: 1. summarizedExperiement(SE) obj for chromVAR 2. Jaspar matrix 
## output: 1. motif x cell (z score) 2. plot: ranked           
source("./libs.R")
##------------------------------------------------------------
## inputs
##------------------------------------------------------------

input.chromVar.res.list <- readRDS(file = "./dat/output.jaspar.dev.res.Rdata")
input.chromVar.jaspar.z <- assays(input.chromVar.res.list$dev)$z

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


p.default.cluster <-  ggplot(input.umap.res,aes(UMAP1,UMAP2)) + 
  geom_point(aes(colour=cell_type_overall),size=1,alpha=0.6) + 
  scale_color_brewer(palette = "Set1")+
  theme_light()

p.default.cluster.sub <-  ggplot(input.umap.res,aes(UMAP1,UMAP2)) + 
  geom_point(aes(colour=cluster),size=1,alpha=0.6) + 
  scale_color_manual(values = cols.subcluster) +
  theme_light()


if(F){
  # PDX1/NKX6-1 are high in beta/delta cells
  # RFX is high in islet cells
  # FOX is enriched across the board
  # HNF1 is high in exocrine cells
  # IRF is high in immune cells
  # ETS-like factors are high in endothelial cells
  # MEF2 is high in endothelial/stellate
  # MAF is high in alpha/beta
}
motif.list <- c("PDX1","RFX","FOXO4","HNF1","IRF","ETS","MEF2B","MAFK")
if(T){
  p.default.cluster.motifs <- lapply(motif.list, fun.plot.project.motif)
  png(filename ="project_overall_avg_motif_all_cell_type.png",width = 10,height = 8,res = 300,units = 'in')
  ggarrange(plotlist  =c(list(p.default.cluster),
                         p.default.cluster.motifs),ncol = 3,nrow = 3)
  dev.off()
  png(filename ="project_overall_avg_motif_cell_subtype.png",width = 10,height = 8,res = 300,units = 'in')
  ggarrange(plotlist  =c(list(p.default.cluster.sub+theme(legend.box.spacing = unit(0,'points'),legend.key.size = unit(1,'points'))),
                         p.default.cluster.motifs),ncol = 3,nrow = 3)
  dev.off()
  
}


# project_overall_avg_motif_all_cell_type_scaled.png
if(T){
  input.chromVar.jaspar.z.scale <- input.chromVar.jaspar.z
  input.chromVar.jaspar.z.scale[input.chromVar.jaspar.z.scale>5] <- 5
  input.chromVar.jaspar.z.scale[input.chromVar.jaspar.z.scale< -5] <- -5
  p.default.cluster.motifs.2 <- lapply(motif.list, 
                                     function(x) fun.plot.project.motif(motif = x,
                                                                        input.chromVar.z = input.chromVar.jaspar.z.scale))
  
  png(filename ="project_overall_avg_motif_all_cell_type_scaled.png",width = 10,height = 8,res = 300,units = 'in')
  #pdf(file ="project_overall_avg_motif_all_cell_type_scaled.pdf",width = 9.6,height = 7.2)
  ggarrange(plotlist  =c(list(p.default.cluster),
                         p.default.cluster.motifs.2),ncol = 3,nrow = 3)
  dev.off()
  
  png(filename ="project_overall_avg_motif_cell_subtype_scaled.png",width = 10,height = 8,res = 300,units = 'in')
  ggarrange(plotlist  =c(list(p.default.cluster.sub+theme(legend.box.spacing = unit(0,'points'),legend.key.size = unit(1,'points'))),
                         p.default.cluster.motifs.2),ncol = 3,nrow = 3)
  dev.off()
  
}





