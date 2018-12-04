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

# beta -------------------------------------------------------------

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

require(cluster)

mat.b <- smoothed.motif.pt.beta.mat.scaled[c(selected.motif$beta,selected.motif$`alpha:beta`),]
pa.b <- pam(mat.b,k=4,cluster.only = T)
pa.b[pa.b==4] <- 1
pa.b["NFYB"] <- 3
#pa.b <- pam(mat.b,k=2,cluster.only = T)
#pa.b[pa.b==3]<- 2

ht.beta <- Heatmap(mat.b,split = pa.b,
                   cluster_columns = F,
                   cluster_rows = T,show_row_dend = F,
                   col = cols.hm.avg.tf(50),
                   show_column_names = F,
                   #row_order = ord,
                   #raster_device = "CairoPNG",
                   row_names_gp=gpar(fontsize=5),name = "ht.beta")

if(T){
  #ord.b <- row_order(ht.beta)
  #n.clusters <- sapply(c(1,2,3),function(x) sum(pa.b==x))
  ord.b <- sapply(c(2,3,1),function(x){
    tmp <- apply(mat.b[pa.b==x,],1,mean)
    names(tmp)[order(tmp)]})
  n.clusters <- sapply(c(2,3,1),function(x) sum(pa.b==x))
  
  
  if(F){
    png(filename = "hmp.pt.beta.smooth.png",height =4,
        width = 6,units = 'in',res = 300)
    layout(mat = matrix(c(1,2,3,4),4,1,byrow = T),
           heights = c(1,n.clusters))
    #n.clusters[2]<- 6
    par(mgp=c(2,.3,0),mar=c(0,4,0,0)+0.1)
    colr <- adjustcolor( brewer.pal(9,"Set1")[1:2],alpha.f = 0.2)
    #colr <- cols.celltype[c('beta_1','beta_2')]
    plot(pd.anno.b$pt,rep(1,nrow(pd.anno.b)),col=colr[pd.anno.b$subtype],xaxs="i",
         yaxt='n',ylab="",xaxt='n',pch=20,cex=2,bty="n")
    #axis(side = 1,at = seq(0,20,by = 5),tcl=-0.2,labels = NA)
    
    par(mgp=c(2,.3,0),mar=c(0.5,4,0,0)+0.1)
    # par(mgp=c(2,.3,0),oma=c(0,0,0,0),mar=c(0.5,4,0,0)+0.1)
    for(i in 1:2){
      image(x=seq(0,pt.max.b,by = 0.05),y=1:length(ord.b[[i]]),z= t(mat.b[ord.b[[i]],]),col = cols.hm.avg.tf(50),
            xlab='', yaxt="n",tcl=-0.2,ylab=NA,xaxt="n",yaxs='i',bty="o")
      axis(side = 2,at = 1:length(ord.b[[i]]),tcl=-0.2,labels = ord.b[[i]],las=1,cex.axis=0.75,hadj = 1)
      axis(side = 1,at = seq(0,pt.max.b,by = 5),tcl=-0.2,labels = NA)
    }
    dev.off()
    system("open hmp.pt.beta.smooth.png")
  }
  
  for(p.lab in c("label","no_label")){
    for(p.format in c("png","pdf")){
      wd = ifelse(p.lab=="no_label",2,6)
      ht = ifelse(p.lab=="no_label",1.5,4)
      fn<- paste0("./figs/fig2/hmp.pt.beta.smooth.",p.lab,".",p.format)
      if(p.format=="png"){
        png(filename = fn,height =ht,width = wd,units = 'in',res = 300)
      }else{
        pdf(file = fn,width = wd,height = ht)
      }
      
      layout(mat = matrix(c(1,2,3,4),4,1,byrow = T),
             heights = c(2,n.clusters))
      n.clusters[2]<- 6
      par(mgp=c(2,.3,0),mar=c(0,4,0,0)+0.1)
      if(p.lab=="no_label") par(mgp=c(2,.3,0),mar=c(0,0,0,0)+0.1)
      colr <- adjustcolor( brewer.pal(9,"Set1")[1:2],alpha.f = 0.2)
      #colr <- cols.celltype[c('beta_1','beta_2')]
      plot(pd.anno.b$pt,rep(1,nrow(pd.anno.b)),col=colr[pd.anno.b$subtype],xaxs="i",
           yaxt='n',ylab="",xaxt='n',pch=20,cex=0.5,bty="n")
      #ylim(c(0.9,1.1))
      #box()
      #axis(side = 1,at = seq(0,20,by = 5),tcl=-0.2,labels = NA)
      
      par(mgp=c(2,.3,0),mar=c(0.5,4,0,0)+0.1)
      if(p.lab=="no_label") par(mgp=c(2,.3,0),mar=c(0,0,0,0)+0.1)
      # par(mgp=c(2,.3,0),oma=c(0,0,0,0),mar=c(0.5,4,0,0)+0.1)
      for(i in 1:3){
        image(x=seq(0,pt.max.b,by = 0.05),y=1:length(ord.b[[i]]),z= t(mat.b[ord.b[[i]],]),col = cols.hm.avg.tf(50),
              xlab='', yaxt="n",tcl=-0.2,ylab=NA,xaxt="n",yaxs='i',bty="o")
        axis(side = 2,at = 1:length(ord.b[[i]]),tcl=-0.2,labels = ord.b[[i]],las=1,cex.axis=0.75,hadj = 1)
        axis(side = 1,at = seq(0,pt.max.b,by = 5),tcl=-0.2,labels = NA)
      }
      dev.off()
      system(paste0("open ",fn))
    }
  }
}

if(F){
  ht.beta <- Heatmap(mat.b,split =factor(pa.b,levels=c(2,3,1)),cluster_columns = F,
                     col = cols.hm.avg.tf(50),show_column_names = F,
                     cluster_rows = F,row_order = ord.b,
                     #raster_device = "CairoPNG",
                     row_names_gp=gpar(fontsize=5),name = "ht")
  #pdf("hmp.pt.beta.smooth.pdf",height = 4,width = 6)
  png(filename = "hmp.pt.beta.smooth.png",height =4,width = 6,units = 'in',res = 300)
  print(ht.beta)
  dev.off()
  system("open hmp.pt.beta.smooth.png")
}



if(T){
  select.gene <- c("FOS::JUN","STAT3","NKX6-1","TEAD1")
  
  p.egs <-ggplot(output.motif.pt%>% 
                   filter(type=="beta",motif.name%in% select.gene)%>%
                   mutate(motif.name =factor(motif.name,levels = select.gene)),
                 aes(pt,zval,color=motif.name)) +
    geom_smooth()+theme_light()+ facet_wrap(~motif.name,scales = "free_y",ncol = 2)+
    coord_cartesian(expand = F)
  ggsave(filename = "beta.pt2.egs.png",width = 6,height = 6,plot = p.egs+theme(legend.position = "none"))   
  system("open beta.pt2.egs.png")
}



# alpha -------------------------------------------------------------------
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

require(cluster)


mat.a <- smoothed.motif.pt.alpha.mat.scaled[c(selected.motif$alpha,selected.motif$`alpha:beta`),]
pa.a <- pam(mat.a,k=2,cluster.only = T)
ht.alpha <- Heatmap(mat.a,split = pa.a,
                    cluster_columns = F,
                    cluster_rows = T,show_row_dend = F,
                    col = cols.hm.avg.tf(50),
                    show_column_names = F,
                    #row_order = ord,
                    #raster_device = "CairoPNG",
                    row_names_gp=gpar(fontsize=5),name = "ht.alpha")

ord <- row_order(ht.alpha)

n.clusters <- sapply(ord,length)
with_lab=F

ord[[2]]<-rev(ord[[2]])
for(p.lab in c("label","no_label")){
  for(p.format in c("png","pdf")){
    wd = ifelse(p.lab=="no_label",2,6)
    ht = ifelse(p.lab=="no_label",1.5*109/66,4*109/66)
    fn<- paste0("./figs/fig2/hmp.pt.alpha.smooth.",p.lab,".",p.format)
    if(p.format=="png"){
      png(filename = fn,height =ht,width = wd,units = 'in',res = 600)
    }else{
      pdf(file = fn,width = wd,height = ht)
    }
    
    
    colr <- adjustcolor( brewer.pal(9,"Set1")[1:2],alpha.f = 0.2)
    
    par(mgp=c(1,.3,0),mar=c(0,.4,0,0)+0.1)
    layout(mat = matrix(c(1,2,3),3,1,byrow = T),
           heights = c(2*109/66,n.clusters))
    
    if(p.lab=="label")  par(mar=c(0,4,0,0)+0.1)
    
    plot(pd.anno.a$pt,rep(1,nrow(pd.anno.a)),
         col=colr[pd.anno.a$subtype],xaxs="i",
         yaxt='n',ylab="",xaxt='n',pch=20,
         cex=ifelse(p.lab=="label",2,.5),
         bty="n")
    
    par(mar=c(0,.4,0,0)+0.1)
    if(p.lab=="label")    par(mar=c(0.5,4,0,0)+0.1)
    
    
    for(i in 1:2){
      if(i==2 & p.lab=="label") par(mar=c(2.5,4,0,0)+0.1)
      image(x=seq(0,pt.max.a,by = 0.05),y=1:length(ord[[i]]),z= t(mat.a[rev(ord[[i]]),]),col = cols.hm.avg.tf(50),
            xlab='', yaxt="n",tcl=-0.2,ylab=NA,xaxt="n",yaxs='i',bty="o",useRaster=F)
      r_labs <- NA
      t_labs <- NA
      if(p.lab=="label"){
        r_labs=rownames(mat.a)[rev(ord[[i]])]
        if(i==2) {title(xlab = "Pseudo-time");t_labs=seq(0,pt.max.a,by = 5)}
      }
      axis(side = 2,at = 1:length(ord[[i]]),tcl=-0.2,
           labels = r_labs,
           las=1,cex.axis=0.5,hadj = 1)
      axis(side = 1,at = seq(0,pt.max.a,by = 5),tcl=-0.2,labels = t_labs)
    }
    dev.off()
    system(paste0("open ",fn))
  }
}



if(T){
  select.gene <- c("FOS::JUN","STAT3","NKX6-1","TEAD1")
  select.gene <- c("MAFK","MAFG","MAFF")
  p.egs <-ggplot(output.motif.pt%>% 
                   filter(type=="alpha",motif.name%in% select.gene)%>%
                   mutate(motif.name =factor(motif.name,levels = select.gene)),
                 aes(pt,zval,color=motif.name)) +
    geom_smooth()+theme_light()+ facet_wrap(~motif.name,scales = "free_y",ncol = 2)+
    coord_cartesian(expand = F)
  
  
  ggsave(filename = "alpha.pt2.egs.png",width = 6,height = 6,plot = p.egs+theme(legend.position = "none"))   
  system("open alpha.pt2.egs.png")
}

promoter.cpm <- fread('../../atacMotif/test/promoter.cpm.csv',header = T)

