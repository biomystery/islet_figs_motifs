## input: 1. summarizedExperiement(SE) obj for chromVAR 2. Jaspar matrix 
## output: 1. motif x cell (z score) 2. plot: ranked           
source("./libs.R")
##------------------------------------------------------------
## inputs
## 1. motif x cell (zval) 2. cell, celltype 3.  cell, pseudotime 
## output : 1. motif x (pseudotime,type)
##------------------------------------------------------------
# 1. motif x cell - z
input.chromVar.res.list <- readRDS(file = "./dat/output.jaspar.dev.res.Rdata")
input.chromVar.jaspar.z <- data.table(assays(input.chromVar.res.list$dev)$z,keep.rownames = T)%>%
  separate(rn,into=c("id","motif.name"),sep = "_")%>%
  select(-one_of("id"))
rm(input.chromVar.res.list)
input.chromVar.jaspar.z <- melt(input.chromVar.jaspar.z,id="motif.name",variable.name = "barcodes",value.name = "zval")

# cell,celltype
input.umap.res <- fread('./dat/Islet_123.MNN_corrected.UMAP.txt',header = T)%>% 
  select(one_of("barcodes","cluster"))%>%
  separate(cluster,into = c("cell_type_overall","subtype"),remove=F)
input.umap.res[is.na(input.umap.res)]<-0

input.alpha.pseduotime <- fread("./dat/alpha.pseudotime.txt",skip = 1,col.names = c("barcodes","pt"))
input.beta.pseduotime <- fread("./dat/beta.pseudotime2.txt",skip = 1,col.names = c("barcodes","pt"))
input.pseudotime <- rbindlist(list(input.alpha.pseduotime[,type:="alpha"],
                                   input.beta.pseduotime[,type:="beta"]))
rm(input.alpha.pseduotime);rm(input.beta.pseduotime)


# aggregate data  --------------------------------------------------------------
# melt
output.motif.pt <- merge(input.chromVar.jaspar.z,input.pseudotime)
output.motif.pt <- output.motif.pt[,pt:=round(pt,2)]
output.motif.pt <- merge(output.motif.pt,input.umap.res[,c(1,4)])


png(filename ="alpha_pseudo_time_all_motif.png",width = 8,height = 6,res = 300,units = 'in')
fun.plt.motif.pt()
dev.off()

png(filename ="alpha_pseudo_time_all_motif_labeled.png",width = 8,height = 16,res = 300,units = 'in')
fun.plt.motif.pt(showname = T)
dev.off()

png(filename ="beta_pseudo_time2_all_motif.png",width = 8,height = 6,res = 300,units = 'in')
fun.plt.motif.pt(x = "beta")
dev.off()

png(filename ="beta_pseudo_time2_all_motif_labeled.png",width = 8,height = 16,res = 300,units = 'in')
pdf(file="beta_pseudo_time2_all_motif_labeled.pdf",width = 8,height = 16)
fun.plt.motif.pt(x = "beta",showname = T)
dev.off()


# select_motifs_for_alpha&beta --------------------------------------------

select.celltypes <- c("alpha_1" , "alpha_2")
dev <- input.chromVar.res.list$dev[,as.character(colData(input.chromVar.res.list$dev)$cell_type) %in% select.celltypes]

## variability 
variability <- computeVariability(dev)
png(filename = "output.jaspar.var.res.alpha.png",width = 4,height = 3,res = 300,units = "in")
plotVariability(variability, use_plotly = FALSE)
dev.off()
variability<- variability%>%arrange(desc(variability))
write.csv(variability,file = "output.jaspar.var.res.alpha.csv")
select.motifs.alpha <- as.character((variability%>%filter(variability>1.2))$name)

png(filename ="alpha_pseudo_time_selected_motif.png",width = 8,height = 6,res = 300,units = 'in')
fun.plt.motif.pt(selected.motif = select.motifs.alpha)
dev.off()


# selectedd_motifs for beta  ----------------------------------------------
if(T){
  select.celltypes <- c("beta_1" , "beta_2")
  dev <- input.chromVar.res.list$dev[,as.character(colData(input.chromVar.res.list$dev)$cell_type) %in% select.celltypes]
  
  ## variability 
  variability <- computeVariability(dev)
  png(filename = "output.jaspar.var.res.beta.png",width = 4,height = 3,res = 300,units = "in")
  plotVariability(variability, use_plotly = FALSE)
  dev.off()
  variability<- variability%>%arrange(desc(variability))
  write.csv(variability,file = "output.jaspar.var.res.beta.csv")
  select.motifs.beta <- as.character((variability%>%filter(variability>1.2))$name)
  
  png(filename ="beta_pseudo_time_selected_motif.png",width = 8,height = 6,res = 300,units = 'in')
  fun.plt.motif.pt(selected.motif = select.motifs.beta)
  dev.off()
  
  png(filename ="beta_pseudo_time_selected_motif_labed.png",width = 8,height = 7,res = 300,units = 'in')
  fun.plt.motif.pt(selected.motif = select.motifs.beta,showname = T)
  dev.off()
  
  png(filename ="beta_pseudo_time_selected_motif_labed_noModes.png",width = 8,height = 7,res = 300,units = 'in')
  fun.plt.motif.pt(selected.motif = select.motifs.beta,showname = T,useModes = F)
  dev.off()
}

fun.plot.var <- function(variability){
  ggplot(variability,aes(-log10(p_value_adj),variability))+
    geom_point(col='grey')+
    theme_light()

}

# smoothed curves  --------------------------------------------------------
select.gene <- c("FOSL1","RFX3","NRL")

p.egs <-ggplot(output.motif.pt%>% 
         filter(type=="beta",motif.name%in% select.gene),
       aes(pt,zval,color=motif.name)) +
  geom_smooth()+theme_light()+ facet_wrap(~motif.name,scales = "free_y",ncol = 1)
ggsave(filename = "beta.pt2.egs.png",width = 6,height = 4)       


require(parallel)

smoothed.motif.pt.beta <- lapply(unique(output.motif.pt$motif.name),function(x){
  tc <- output.motif.pt%>% 
    filter(type=="beta",motif.name== x)%>%
    arrange(pt)
  print(x)
  fun.gamSmooth(tc)
})

plot(1:nrow(smoothed.motif.pt.beta.mat),log2(sort(apply(smoothed.motif.pt.beta.mat,1,var))))
var.motifs <- which(apply(smoothed.motif.pt.beta.mat,1,var)>0.25)
smoothed.motif.pt.beta.mat.scaled <- t(apply(smoothed.motif.pt.beta.mat, 1, function(x) (x-min(x))/(max(x)-min(x))))
rownames(smoothed.motif.pt.beta.mat.scaled)<- unique(output.motif.pt$motif.name)
require(cluster)
pam(smoothed.motif.pt.beta.mat.scaled[var.motifs,],k=2)
require(ComplexHeatmap)
mat <- smoothed.motif.pt.beta.mat.scaled[var.motifs,]
pa <- pam(mat,k=4)

Heatmap(mat,split =factor(pa$clustering,levels=c(2,3,4,1)),cluster_columns = F,
        col = cols.hm.avg.tf(50),show_column_names = F,
        raster_device = "CairoPNG")

pheatmap(mat,cluster_cols = F,
         show_rownames = T,show_colnames = F,
         color = cols.hm.avg.tf(50),clust)

saveRDS(list(mat=smoothed.motif.pt.beta.mat.scaled,t=seq(0,20,by = .05)),
            file = 'dat/smoothed.motif.pt.beta.mat.scaled.Rds')





     