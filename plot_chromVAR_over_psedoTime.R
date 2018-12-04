## input: 1. summarizedExperiement(SE) obj for chromVAR 2. Jaspar matrix 
## output: 1. motif x cell (z score) 2. plot: ranked           
source("./libs.R")
##------------------------------------------------------------
## inputs
## 1. motif x cell (zval) 2. cell, celltype 3.  cell, pseudotime 
## output : 1. motif x (pseudotime,type)
##------------------------------------------------------------
# motif x cell - z
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

# pseudotime
input.alpha.pseduotime <- fread("./dat/alpha.pseudotime.txt",skip = 1,col.names = c("barcodes","pt"))
input.beta.pseduotime <- fread("./dat/beta.pseudotime.txt",skip = 1,col.names = c("barcodes","pt"))
input.pseudotime <- rbindlist(list(input.alpha.pseduotime[,type:="alpha"],
                                   input.beta.pseduotime[,type:="beta"]))
rm(input.alpha.pseduotime);rm(input.beta.pseduotime)

# melt
output.motif.pt <- merge(input.chromVar.jaspar.z,input.pseudotime)
output.motif.pt <- output.motif.pt[,pt:=round(pt,2)]
output.motif.pt <- merge(output.motif.pt,input.umap.res[,c(1,4)])

# cell type anno
output.motif.pt.anno <-output.motif.pt%>%
  group_by(pt,type)%>%
  summarise(subtype=Modes(as.numeric(subtype))[1])

# aggregate data  --------------------------------------------------------------



png(filename ="alpha_pseudo_time_all_motif.png",width = 8,height = 6,res = 300,units = 'in')
fun.plt.motif.pt()
dev.off()

png(filename ="alpha_pseudo_time_all_motif_labeled.png",width = 8,height = 16,res = 300,units = 'in')
fun.plt.motif.pt(showname = T)
dev.off()

png(filename ="beta_pseudo_time_all_motif.png",width = 8,height = 6,res = 300,units = 'in')
fun.plt.motif.pt(x = "beta")
dev.off()

png(filename ="beta_pseudo_time_all_motif_labeled.png",width = 8,height = 16,res = 300,units = 'in')
fun.plt.motif.pt(x = "beta",showname = T)
dev.off()


# select_motifs_for_alpha&beta --------------------------------------------

if(F){
  select.celltypes <- c("alpha_1" , "alpha_2",  "alpha_3")
  dev <- input.chromVar.res.list$dev[,as.character(colData(input.chromVar.res.list$dev)$cell_type) %in% select.celltypes]
  
  ## variability 
  variability <- computeVariability(dev)
  png(filename = "output.jaspar.var.res.alpha.png",width = 4,height = 3,res = 300,units = "in")
  plotVariability(variability, use_plotly = FALSE)
  dev.off()
  variability<- variability%>%arrange(desc(variability))
  write.csv(variability,file = "output.jaspar.var.res.alpha.csv")
}

alpha.variability <- read.csv("./dat/output.jaspar.var.res.alpha.csv",row.names = 1)
select.motifs.alpha <- as.character((alpha.variability%>%filter(variability>1.2))$name)
beta.variability <- read.csv("./dat/output.jaspar.var.res.beta.csv",row.names = 1)
all.variability <- read.csv("./dat/output.jaspar.var.res.csv",row.names = 1)
islet.variability <- read.csv("./dat/output.jaspar.var.res.subtype.csv",row.names = 1)
pd.variablility <- rbind(data.frame(alpha.variability,type="alpha"),
                         data.frame(beta.variability,type="beta"),
                         data.frame(all.variability,type="all"),
                         data.frame(islet.variability,type="islet"))%>%
  group_by(type)%>%
  arrange(variability)%>%
  mutate(ord=order(variability,decreasing = T))

ggplot(pd.variablility,
       aes(ord,variability)) + geom_point(aes(colour=type))+
  geom_hline(yintercept = 1.2)+ theme_light()+ scale_x_log10()

write.csv(file = "top10.var.motifs.csv",pd.variablility%>%filter(ord<=10)%>% arrange(ord))

if(F){
  png(filename ="alpha_pseudo_time_selected_motif.png",width = 8,height = 6,res = 300,units = 'in')
  fun.plt.motif.pt(selected.motif = select.motifs.alpha)
  dev.off()
  
}


# smooth 
i <- 1 

smoothed.motif.pt.alpha<- lapply(select.motifs.alpha,function(x){
  tc <- output.motif.pt%>% 
    filter(type=="alpha",motif.name== x)%>%
    arrange(pt)
  print(paste(x,i))
  tmp <- fun.gamSmooth(tc,new.tps = seq(0,max(tc$pt),by=0.01))
  i <<- i+ 1
  (tmp-min(tmp))/(max(tmp)-min(tmp))
})

smoothed.motif.pt.alpha <- do.call(rbind, smoothed.motif.pt.alpha)
rownames(smoothed.motif.pt.alpha)<- select.motifs.alpha

ph <- pheatmap(smoothed.motif.pt.alpha,cluster_cols = F,
         show_rownames = T,show_colnames = F,
         cutree_rows = 4,
         color = cols.hm.avg.tf(50))

pa <- pam(smoothed.motif.pt.alpha,k=3,cluster.only = T)
Heatmap(smoothed.motif.pt.alpha,split =pa,cluster_columns = F,
        col = cols.hm.avg.tf(50),show_column_names = F,row_names_gp = gpar(fontsize = 8))

pd <- smoothed.motif.pt.alpha[,-401]
image(x = seq(0,20,by=0.05),
      y= 0:nrow(pd),
      z=t(pd[ph$tree_row$order,]),
      col = cols.hm.avg.tf(50))
## eg
select.gene <- c("FOSL1","ATF4","MEF2A","CEBPB","IRF8","MAFK")

p.egs <-ggplot(output.motif.pt%>% 
                 filter(type=="alpha",motif.name%in% select.gene),
               aes(pt,zval,color=motif.name)) +
  geom_smooth()+theme_light()+ facet_wrap(~motif.name,scales = "free_y",ncol = 1)


ggplot(pd.anno%>%filter(type=="alpha"),aes(pt,subtype)) + 
  geom_point(alpha=.5,shape=21,size=3,color="grey10",aes(fill=as.factor(subtype)))+
  theme_light()

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
select.gene <- c("FOS","FOXO6")

ggplot(output.motif.pt%>% 
         filter(type=="beta",motif.name%in% select.gene),
       aes(pt,zval))+ geom_smooth(aes(color=motif.name))+
  theme_light()
ggsave(filename = "beta.pt.egs.png",width = 6,height = 4)       


