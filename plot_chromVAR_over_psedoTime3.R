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
#fwrite(input.chromVar.jaspar.z,file = "output.chromVar.jaspar.z.csv")
input.chromVar.jaspar.z <- melt(input.chromVar.jaspar.z,id="motif.name",variable.name = "barcodes",value.name = "zval")

# cell,celltype
input.umap.res <- fread('./dat/Islet_123.MNN_corrected.UMAP.txt',header = T)%>% 
  select(one_of("barcodes","cluster"))%>%
  separate(cluster,into = c("cell_type_overall","subtype"),remove=F)
input.umap.res[is.na(input.umap.res)]<-0

input.alpha.pseduotime <- fread("./dat/alpha.pseudotime2.txt",skip = 1,col.names = c("barcodes","pt"))
input.beta.pseduotime <- fread("./dat/beta.pseudotime2.txt",skip = 1,col.names = c("barcodes","pt"))
input.pseudotime <- rbindlist(list(input.alpha.pseduotime[,type:="alpha"],
                                   input.beta.pseduotime[,type:="beta"]))
rm(input.alpha.pseduotime);rm(input.beta.pseduotime)


# aggregate data  --------------------------------------------------------------
# melt
output.motif.pt <- merge(input.chromVar.jaspar.z,input.pseudotime)
output.motif.pt <- output.motif.pt[,pt:=round(pt,2)]
output.motif.pt <- merge(output.motif.pt,input.umap.res[,c(1,4)])
output.motif.pt <- output.motif.pt%>% 
  filter(subtype!='3')

ttest.res <- fread("./dat/ttest.res.csv")
dmotifs.list <- sapply(c("alpha","beta","delta"), function(x) subset(ttest.res,celltype==x & selected)$motif)
tmp <- gplots::venn(dmotifs.list)
dmotifs.list.inter <- attr(tmp,"intersections")

# fig2B:beta hm -------------------------------------------------------------

selected.motif <- list(`alpha`=c(dmotifs.list.inter[["alpha"]],
                                 dmotifs.list.inter[["alpha:delta"]]),
                       `beta`=c(dmotifs.list.inter[["beta"]],
                                dmotifs.list.inter[["beta:delta"]]),
                       `alpha:beta`=c(dmotifs.list.inter[["alpha:beta"]],
                                      dmotifs.list.inter[["alpha:beta:delta"]]))
pt.max.b <- max((output.motif.pt%>%filter(type=="beta"))$pt)

system.time(smoothed.motif.pt.beta <- mclapply(unlist(selected.motif),function(x){
  tc <- output.motif.pt%>% 
    filter(type=="beta",motif.name== x)%>%
    arrange(pt)
  print(x)
  fun.gamSmooth(tc,new.tps = seq(0,pt.max.b,by = 0.05))
},mc.cores = 10))

smoothed.motif.pt.beta.mat <- do.call(rbind,smoothed.motif.pt.beta)
rownames(smoothed.motif.pt.beta.mat) <- unlist(selected.motif)
smoothed.motif.pt.beta.mat.scaled <- t(apply(smoothed.motif.pt.beta.mat, 1, function(x) (x-min(x))/(max(x)-min(x))))

pd.anno.b <-output.motif.pt%>%
  filter(type=="beta")%>%
  group_by(pt)%>%
  summarise(subtype=Modes(as.numeric(subtype))[1])%>%
  mutate(subtype=as.factor(subtype))%>%
  as.data.frame()

mat.b <- smoothed.motif.pt.beta.mat.scaled[c(selected.motif$beta,selected.motif$`alpha:beta`),]
h <- pheatmap(mat.b,cluster_cols = F,show_colnames = F,cellheight = 5)
ord.b <- (h$tree_row$order)
names(ord.b) <- rownames(mat.b)[ord.b]
pheatmap(mat.b,cluster_cols = F,show_colnames = F,cellheight = 5)

n.clusters <- nrow(mat.b)

for(p.lab in c("label","no_label")){
  for(p.format in c("png","pdf")){
    wd = ifelse(p.lab=="no_label",2,6)
    ht = ifelse(p.lab=="no_label",1.5,4)
    fn<- paste0("./figs/fig2/subfig2B_beta_hm",p.lab,".",p.format)
    if(p.format=="png"){
      png(filename = fn,height =ht,width = wd,units = 'in',res = 600)
    }else{
      pdf(file = fn,width = wd,height = ht)
    }
    
    
    colr <- adjustcolor( brewer.pal(9,"Set1")[1:2],alpha.f = 0.2)
    par(mgp=c(1,.3,0),mar=c(0,.4,0,0)+0.1)
    layout(mat = matrix(c(1,2,3),3,1,byrow = T),
           heights = c(2,n.clusters))
    
    if(p.lab=="label")  par(mar=c(0,4,0,0)+0.1)
    
    plot(pd.anno.b$pt,rep(1,nrow(pd.anno.b)),
         col=colr[pd.anno.b$subtype],xaxs="i",
         yaxt='n',ylab="",xaxt='n',pch=20,
         cex=ifelse(p.lab=="label",2,.5),
         bty="n")
    
    par(mar=c(0,.4,0,0)+0.1)
    if(p.lab=="label") par(mar=c(2.5,4,0,0)+0.1)
    image(x=seq(0,pt.max.b,by = 0.05),y=1:length(ord.b),z= t(mat.b[ord.b,]),col = cols.hm.avg.tf(50),
          xlab='', ylab=NA,yaxt="n",tcl=-0.2,xaxt="n",yaxs='i',bty="o",useRaster=F)
    r_labs <- NA
    t_labs <- NA
    if(p.lab=="label"){
      r_labs=names(ord.b)
      title(xlab = "Pseudo-time");t_labs=seq(0,pt.max.b,by = 5)}
  
  axis(side = 2,at = 1:length(ord.b),tcl=-0.2,
       labels = r_labs,
       las=1,cex.axis=0.5,hadj = 1)
  axis(side = 1,at = seq(0,pt.max.b,by = 5),tcl=-0.2,labels = t_labs)
  
  dev.off()
  system(paste0("open ",fn))
  }
}

pdf(file = './figs/fig2/subfig2B_beta_hm_p.pdf')
pheatmap(mat.b[ord.b,],cluster_rows = F,cluster_cols = F,show_colnames = F)
dev.off()

### individual plots 
if(T){
  fn="./figs/fig2/subfig2B_beta_hm.pdf"
  #png(filename = fn,width = 2,height = 1.5,units = "in",res = 300)
  pdf(file = fn,width = 2,height = 1.5)
  par(mgp=c(2,.3,0),mar=c(0,0,0,0))
  image(x=seq(0,pt.max.b,by = 0.05),y=1:length(ord.b),z= t(mat.b[ord.b,]),col = cols.hm.avg.tf(50),
        xlab='', yaxt="n",tcl=-0.2,ylab=NA,xaxt="n",yaxs='i',bty="o")
  axis(side = 2,at = 1:length(ord.b),tcl=0.05,labels = NA,las=1,cex.axis=0.75,hadj = 1)
  axis(side = 1,at = seq(0,pt.max.b,by = 5),tcl=0.05,labels = NA)
  dev.off()
  system(paste0("open ",fn))
}

if(T){
  
  fn="./figs/fig2/subfig2B_beta_hm_anno.pdf"
  #png(filename = fn,width = 2,height = 1.5,units = "in",res = 300)
  pdf(file = fn,width = 2*4,height = 0.5*2)
  par(mgp=c(2,.3,0),mar=c(0,0,0,0))
  colr <- adjustcolor( brewer.pal(9,"Set1")[1:2],alpha.f = 0.2)# 
  #colr <-  brewer.pal(9,"Set1")[1:2]
  plot(pd.anno.b$pt,pd.anno.b$subtype,col=colr[pd.anno.b$subtype],xaxs="i",
       yaxt='n',ylab="",xaxt='n',pch=20,cex=.5,bty="n",ylim=c(-2,5))
  #ylim(c(-2,5))
  dev.off()
  system(paste0("open ",fn))
}

# fig2B:alpha hm-------------------------------------------------------------------
pt.max.a <- max((output.motif.pt%>%filter(type=="alpha"))$pt)

system.time(smoothed.motif.pt.alpha <- mclapply(unlist(selected.motif),function(x){
  tc <- output.motif.pt%>% 
    filter(type=="alpha",motif.name== x)%>%
    arrange(pt)
  print(x)
  fun.gamSmooth(tc,new.tps = seq(0,pt.max.a,by = 0.05))
},mc.cores = 10))

smoothed.motif.pt.alpha.mat <- do.call(rbind,smoothed.motif.pt.alpha)
rownames(smoothed.motif.pt.alpha.mat) <- unlist(selected.motif)
smoothed.motif.pt.alpha.mat.scaled <- t(apply(smoothed.motif.pt.alpha.mat, 1, function(x) (x-min(x))/(max(x)-min(x))))


pd.anno.a <-output.motif.pt%>%
  filter(type=="alpha")%>%
  group_by(pt)%>%
  summarise(subtype=Modes(as.numeric(subtype))[1])%>%
  mutate(subtype=as.factor(subtype))%>%
  as.data.frame()


mat.a <- smoothed.motif.pt.alpha.mat.scaled[c(selected.motif$alpha,selected.motif$`alpha:beta`),]
pt.max.a <- max((output.motif.pt%>%filter(type=="alpha"))$pt)


h <- pheatmap(mat.a,cluster_cols = F,show_colnames = F)
ord.a <- (h$tree_row$order)
names(ord.a) <- rownames(mat.a)[ord.a]
n.clusters <- nrow(mat.a)

for(p.lab in c("label","no_label")){
  for(p.format in c("png","pdf")){
    wd = ifelse(p.lab=="no_label",2,6*1.5)
    ht = ifelse(p.lab=="no_label",1.5,4*1.5)
    fn<- paste0("./figs/fig2/subfig2B_alpha_hm",p.lab,".",p.format)
    if(p.format=="png"){
      png(filename = fn,height =ht,width = wd,units = 'in',res = 600)
    }else{
      pdf(file = fn,width = wd,height = ht)
    }
    
    
    colr <- adjustcolor( brewer.pal(9,"Set1")[1:2],alpha.f = 0.2)
    par(mgp=c(1,.3,0),mar=c(0,.4,0,0)+0.1)
    layout(mat = matrix(c(1,2,3),3,1,byrow = T),
           heights = c(2,46))
    
    if(p.lab=="label")  par(mar=c(0,4,0,0)+0.1)
    
    plot(pd.anno.a$pt,rep(1,nrow(pd.anno.a)),
         col=colr[pd.anno.a$subtype],xaxs="i",
         yaxt='n',ylab="",xaxt='n',pch=20,
         cex=ifelse(p.lab=="label",2,.5),
         bty="n")
    
    par(mar=c(0,.4,0,0)+0.1)
    if(p.lab=="label") par(mar=c(2.5,4,0,0)+0.1)
    image(x=seq(0,pt.max.a,by = 0.05),y=1:length(ord.a),z= t(mat.a[ord.a,]),col = cols.hm.avg.tf(50),
          xlab='', ylab=NA,yaxt="n",tcl=-0.2,xaxt="n",yaxs='i',bty="o",useRaster=F)
    r_labs <- NA
    t_labs <- NA
    if(p.lab=="label"){
      r_labs=names(ord.a)
      title(xlab = "Pseudo-time");t_labs=seq(0,pt.max.b,by = 5)
      axis(side = 2,at = 1:length(ord.a),tcl=-0.2,
           labels = r_labs,
           las=1,cex.axis=0.5,hadj = 1)
      }
    
    axis(side = 1,at = seq(0,pt.max.b,by = 5),tcl=-0.2,labels = t_labs)
    
    dev.off()
    system(paste0("open ",fn))
  }
}

pdf(file = './figs/fig2/subfig2B_alpha_hm_p.pdf')
pheatmap(mat.a[ord.a,],cluster_rows = F,cluster_cols = F,show_colnames = F)
dev.off()





# fig2B:alpha_egs ----------------------------------------------------------------

if(T){
  select.gene <- c("FOS::JUN","STAT3","NKX6-1","TEAD1")
  select.gene <- c("RFX3","FOSL1")
  p.egs <-ggplot(output.motif.pt%>% 
                   filter(type=="alpha",motif.name%in% select.gene)%>%
                   mutate(motif.name =factor(motif.name,levels = select.gene)),
                 aes(pt,zval)) +
    geom_smooth(color="black")+
    facet_wrap(~motif.name,scales = "free_y",ncol = 1)+
    coord_cartesian(expand = F)+
    scale_x_continuous(breaks = c(0,5))+
  theme_bw()
  
  fn <- "./figs/fig2/subfig2B_alpha_egs_lab.pdf"
  ggsave(filename = fn,width = 2,height = 2,
         plot = p.egs+theme(legend.position = "none"),
        scale = 2)   
  system(paste0("open ",fn))
  fn <- "./figs/fig2/subfig2B_alpha_egs.pdf"
  ggsave(filename = fn,width = 2.5,height = 2.5,
         plot = p.egs+theme(legend.position = "none",
                            text = element_blank(),
                            plot.margin =margin(b=.1,l=.5,t=.1,r = .1,unit = "points"))
         ,
         scale = 2)  
  #par(mar=c(0,.4,0,0)+0.1)
  system(paste0("open ",fn))
}


# fig2B:beta_egs ----------------------------------------------------------------

if(T){
  select.gene <- c("FOS::JUN","STAT3","NKX6-1","TEAD1")
  select.gene <- c("RFX3","FOSL1")
  p.egs <-ggplot(output.motif.pt%>% 
                   filter(type=="beta",motif.name%in% select.gene)%>%
                   mutate(motif.name =factor(motif.name,levels = select.gene)),
                 aes(pt,zval)) +
    geom_smooth(color="black")+
    facet_wrap(~motif.name,scales = "free_y",ncol = 1)+
    coord_cartesian(expand = F)+
    scale_x_continuous(breaks = c(0,5,10,15,20))+
    theme_bw()
  
  fn <- "./figs/fig2/subfig2B_beta_egs_lab.pdf"
  ggsave(filename = fn,width = 2,height = 2,
         plot = p.egs+theme(legend.position = "none"),
         scale = 2)   
  system(paste0("open ",fn))
  fn <- "./figs/fig2/subfig2B_beta_egs.pdf"
  ggsave(filename = fn,width = 2.5,height = 2.5,
         plot = p.egs+theme(legend.position = "none",
                            text = element_blank(),
                            plot.margin =margin(b=.1,l=.5,t=.1,r = .1,
                                                unit = "points"))
         ,
         scale = 2)  
  #par(mar=c(0,.4,0,0)+0.1)
  system(paste0("open ",fn))
}


promoter.cpm <- fread('../../atacMotif/test/promoter.cpm.csv',header = T)

