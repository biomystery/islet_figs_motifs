## input: 1. summarizedExperiement(SE) obj for chromVAR 2. Jaspar matrix 
## output: 1. motif x cell (z score) 2. plot: ranked           
source("./libs.R")
##------------------------------------------------------------
## inputs
##------------------------------------------------------------

ttest.res <- fread("./dat/ttest.res.csv")

dmotifs.list <- sapply(c("alpha","beta"), function(x) subset(ttest.res,celltype==x & selected)$motif)
tmp <- gplots::venn(dmotifs.list)
dmotifs.list.inter <- attr(tmp,"intersections")


input.chromVar.res.list <- readRDS(file = "./dat/output.jaspar.dev.res.Rdata")
input.chromVar.jaspar.z <- assays(input.chromVar.res.list$dev)$z

select.celltypes <- c("alpha_1" , "alpha_2" , "beta_1",   "beta_2" )



input.umap.res <- readRDS("./dat/res.umap.Rdata")


# add Z score of motifs  -------------------------------------------------
input.umap.res.sub <- subset(input.umap.res,cluster %in% select.celltypes)

p.default.cluster.sub <-  ggplot(input.umap.res.sub,aes(UMAP1,UMAP2)) + 
  geom_point(aes(colour=cluster),size=1,alpha=0.6) + 
  scale_color_manual(values = cols.celltype) +
  theme_light()

png(filename ="umap_cell_subtype_ab.png",width = 4,height = 3,res = 300,units = 'in')
print(p.default.cluster.sub)
dev.off()
system("open umap_cell_subtype_ab.png")

ord <- c("alpha","beta","alpha:beta")
motif.list<- sapply(dmotifs.list.inter,function(x) x[1])[ord]
select.gene <- c("FOS::JUN","STAT3","NKX6-1","TEAD1","TEAD2","TEAD3")
motif.list <- unlist(sapply(select.gene, function(g){
  
  names(dmotifs.list.inter)[grep(g,dmotifs.list.inter)]
  
}))


if(T){
  p.default.cluster.motifs <- lapply(names(motif.list), function(x) 
    fun.plot.project.motif(motif = x,umap.res = input.umap.res.sub,rescale = F)+
      ggtitle(paste(motif.list[x],x,sep = "_")))
  png(filename ="project_overall_avg_motif_selected_cell_subtype.png",width =10,height = 8,res = 300,units = 'in')
  ggarrange(plotlist=list.prepend(
    p.default.cluster.motifs,p.default.cluster.sub),ncol = 3,nrow = 3)
  dev.off()
  system("open project_overall_avg_motif_selected_cell_subtype.png")
}


# project_overall_avg_motif_all_cell_type_scaled.png
if(T){
  input.chromVar.jaspar.z.scale <- input.chromVar.jaspar.z
  input.chromVar.jaspar.z.scale[input.chromVar.jaspar.z.scale>5] <- 5
  input.chromVar.jaspar.z.scale[input.chromVar.jaspar.z.scale< -5] <- -5
  
  p.default.cluster.motifs.2 <- lapply(names(motif.list), 
                                     function(x) fun.plot.project.motif(motif = x,
                                                                        umap.res = input.umap.res.sub,
                                                                        input.chromVar.z = input.chromVar.jaspar.z.scale)+
    ggtitle(paste(motif.list[x],x,sep = "_")))
  
  png(filename ="project_overall_avg_motif_selected_cell_subtype_scaled2.png",width = 10,height = 8,res = 300,units = 'in')
  ggarrange(plotlist=
              list.prepend(
                p.default.cluster.motifs.2,p.default.cluster.sub)
              ,ncol = 3,nrow = 3)
  dev.off()
  system("open project_overall_avg_motif_selected_cell_subtype_scaled2.png")
  
}





