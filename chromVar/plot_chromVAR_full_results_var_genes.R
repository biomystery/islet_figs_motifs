## input: 1. summarizedExperiement(SE) obj for chromVAR 2. Jaspar matrix 
## output: 1. motif x cell (z score) 
##         2. cell annotation (type,subtype)
source("./libs.R")
##------------------------------------------------------------
## inputs
##------------------------------------------------------------

input.chromVar.res.list <- readRDS(file = "output.jaspar.dev.res.Rdata")
input.chromVar.jaspar.z <- data.table(assays(input.chromVar.res.list$dev)$z,keep.rownames = T)
input.umap.res <- fread('./Islet_123.MNN_corrected.UMAP.txt',header = T)
input.chromVar.jaspar.var <- fread("./output.jaspar.var.res.csv")

# filter unkonwn
input.umap.res <- input.umap.res %>% 
  separate(cluster,into = c("cell_type_overall","subtype"),remove=F)
input.umap.res[is.na(input.umap.res)]<-0



# aggregate data  --------------------------------------------------------------
# melt
input.chromVar.jaspar.z.agg <- melt(input.chromVar.jaspar.z,id="rn",variable.name = "barcodes",value.name = "zval")

# add celltype
input.chromVar.jaspar.z.agg <- merge(input.chromVar.jaspar.z.agg,input.umap.res)

select.motifs <- (input.chromVar.jaspar.var%>%filter(variability>1.2))$name

# all heatmap -------------------------------------------------------------
all.motif.barcodes <- input.chromVar.jaspar.z.agg%>%
  separate(rn,into = c("id","name"),sep = "_")

all.motif.barcodes <- dcast(all.motif.barcodes[,c("barcodes","name","zval")],name~barcodes,value.var="zval")
pd <- all.motif.barcodes[,-1]
pd[pd>5] <- 5
pd[pd < -5] <- -5
setDF(pd); rownames(pd) <- all.motif.barcodes$name
pd <- pd[select.motifs,]
#all.motif.barcodes[,-1] <- pd 

all.motif.barcodes.anno <- input.umap.res[,c(1,4,5)]
all.motif.barcodes.anno<- all.motif.barcodes.anno%>%
  as.data.frame()%>%
  column_to_rownames("barcodes")
all.motif.barcodes.anno<- all.motif.barcodes.anno[colnames(all.motif.barcodes[,-1]),]
#rm(list = grep("input",ls(),value = T))


## cluster on each direction 
system.time({
  d <- parDist(as.matrix(pd))
  hc <- hclust(d, "ward.D2")
  hc$labels <- rownames(pd)
  #rm(d)
  hc1.dend <- as.dendrogram(hc)
  # http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning
  plot(hc1.dend,nodePar=list(lab.cex=0.5,pch=c(NA,NA)),ylab="",yaxt="n",horiz=T,leaflab = "none")
})
 
#d2 <- readRDS("d2.Rdata")
system.time({
  require(parallelDist)
  d2 <- parDist(t(pd),threads = 10)
  #saveRDS(d2,"d2.Rdata")
  system.time(hc2 <- hclust(d2,"ward.D2"))
  hc2.dend <- as.dendrogram(hc2)
})

if(T){
  #par(mfrow=c(1,2))
  #layout(matrix(c(1,2),1,2,byrow = T),widths = c(1,4))
  par(fig=c(0,0.3,0,1), new=F)
  plot(hc1.dend,nodePar=list(lab.cex=0.5,pch=c(NA,NA)),ylab="",yaxt="n",horiz=T,leaflab = "none")
  par(fig=c(0.2,1,0,1), new=T)
  imagesc(C = pd.ordered.anno[,1:1000],col = cols.cluster,ylab="",xlab="",yaxt="n",xaxt="n")
}



pd.ordered <- as.matrix(pd)[hc$order,hc2$order]



if(T){
  #require(ggdendro)
  #ggdendrogram(hc2,leaf_labels = F) 
  ## plot 1: 
  png('./select.motif.barcodes.clust.png',width = 960)
  plot(hc2.dend,leaflab = "none",ylim=c(25,55))   
  dev.off()

  png('./select.motif.barcodes.png',width = 960)
  imagesc(C = pd.ordered,col = cols.hm.zval.fun(50))
  dev.off()
}

# column anno - cell type
tmp <-as.numeric(as.factor(all.motif.barcodes.anno$cell_type_overall))[hc2$order]
pd.ordered.anno <- matrix(rep(tmp,each=nrow(pd.ordered)),nrow=nrow(pd.ordered))

# column anno - subcell type 
tmp <-as.numeric(as.factor(all.motif.barcodes.anno$cluster))[hc2$order]
pd.ordered.anno.sub <- matrix(rep(tmp,each=20),nrow=20)

if(T){
  png('./select.motif.barcodes.anno.png',width = 960,height = 280)
  imagesc(C = pd.ordered.anno,col = cols.cluster)
  dev.off()
}

if(T){
  png('./all.motif.barcodes.anno.sub.png',width = 960)
  imagesc(C = pd.ordered.anno.sub,col = cols.subcluster)
  dev.off()
}




