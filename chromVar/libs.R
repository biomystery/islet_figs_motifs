rm(list=ls())
tryCatch(dev.off(),error=function(e){})
gc()


# libs --------------------------------------------------------------------
require(SummarizedExperiment)
require(data.table)
require(tidyverse)
require(ggplot2)
require(scales)
require(ggpubr)
require(ComplexHeatmap)
require(pheatmap)
require(RColorBrewer)
require(parallelDist)
require(matlab)
require(cluster)
require(rlist)
require(chromVAR)
require(ggdendro)
require(extrafont)  # fonts 
# https://blog.revolutionanalytics.com/2012/09/how-to-use-your-favorite-fonts-in-r-charts.html
suppressMessages(loadfonts())
require(ggthemes)
# #https://rpubs.com/Koundy/71792 theme_Publication
require(ComplexHeatmap)
cols.hm.avg.tf <- colorRampPalette(c(rgb(249, 249, 212, maxColorValue = 255), rgb(60, 
    181, 195, maxColorValue = 255), rgb(30, 35, 86, maxColorValue = 255)), space = "Lab")



theme_pubr<- function (base_size = 11, base_family = "Arial") 
{
    theme_foundation() + theme(line = element_line(colour = "black", 
        lineend = "round", linetype = "solid"), rect = element_rect(fill = "white", 
        colour = "black", linetype = "solid"), text = element_text(colour = "black", 
        face = "plain", family = base_family, size = base_size, 
        vjust = 0.5, hjust = 0.5, lineheight = 0.5), panel.background = element_blank(), 
        plot.background = element_blank(), panel.border = element_rect(colour = "black", 
            fill = NA), panel.grid = element_blank(), strip.background = element_rect(colour = NA), 
        legend.key = element_rect(colour = NA), title = element_text(size = rel(1)), 
        plot.title = element_text(size = rel(1.2), face = "bold"), 
        strip.text = element_text(), axis.ticks.length = unit(1, 
    "mm"))
}


# colors ------------------------------------------------------------------
cols.subcluster <- c('red4','red3','red1','steelblue4','steelblue1','green4','green1',
                     'purple4','purple1',brewer.pal(9,'Set1')[5:9])
cols.celltype <- c(acinar='grey',alpha_1="darkred",alpha_2="lightpink",beta_1="darkgreen",beta_2="lightgreen",
                   delta_1="orange",delta_2="gold",ductal="cyan",
                   gamma="purple",exocrine="black",endothelial_1="royalblue",endothelial_2="navy",immune="blue",
                   stellate="skyblue",glial="cyan")
cols.celltype["beta"] = rgb(t((col2rgb(cols.celltype["beta_1"]) + col2rgb(cols.celltype["beta_2"]))/2), 
    maxColorValue = 255)
cols.celltype["alpha"] = rgb(t((col2rgb(cols.celltype["alpha_1"]) + col2rgb(cols.celltype["alpha_2"]))/2), 
    maxColorValue = 255)
cols.celltype["delta"] = rgb(t((col2rgb(cols.celltype["delta_1"]) + col2rgb(cols.celltype["delta_2"]))/2), 
    maxColorValue = 255)
cols.celltype["endothelial"] = rgb(t((col2rgb(cols.celltype["endothelial_1"]) + col2rgb(cols.celltype["endothelial_2"]))/2), 
    maxColorValue = 255)


cols.cluster <- brewer.pal(9,'Set1')
cols.hm.avg.tf <- colorRampPalette(c(rgb(249,249,212,maxColorValue = 255),
                                     rgb(60,181,195,maxColorValue = 255),
                                     rgb(30,35,86,maxColorValue = 255)),
                                   space="Lab")
#cols.hm.motif <-fread('BuDRd_18.txt',skip = 2)[,1:3]
#cols.hm.motif <- colorRampPalette(apply(cols.hm.motif,1,function(x) rgb(x[1],x[2],x[3])))
cols.hm.zval.fun <- colorRampPalette(c("dodgerblue4","deepskyblue","darkorange","firebrick4"),space="Lab")
cols.hm.zval.fun.2 <- colorRampPalette(c("dodgerblue4","deepskyblue","darkorange","firebrick4"))
plotColScale <- function(x) imagesc(matrix(rep(1:length(x),each=20),nrow=20),col = x)
if(F){
  par(mfrow=c(2,1))
  plotColScale(cols.hm.zval.fun.2(60))
  plotColScale(cols.hm.zval.fun(60))
par(mfrow=c(1,1))
}

#https://eos.org/features/the-end-of-the-rainbow-color-schemes-for-improved-data-graphics
fun.importRGB<- function(r,g,b) rgb(r,g,b,maxColorValue = 255)
cols.BlGy<- c(fun.importRGB(19,146,218),
                                      fun.importRGB(84,181,233),
                                      fun.importRGB(140,204,231),
                                      fun.importRGB(196,227,242),
                                      fun.importRGB(220,221,227),
                                      fun.importRGB(142,147,153),
                                      fun.importRGB(96,103,108),
                                      fun.importRGB(51,59,60))

cols.Spectrum<- colorRampPalette(c(fun.importRGB(153,25,51), #n=11
                               fun.importRGB(199,53,66),
                               fun.importRGB(229,115,99),
                               fun.importRGB(237,171,121),
                               fun.importRGB(248,217,157),
                               fun.importRGB(242,239,192),
                               fun.importRGB(219,233,239),
                               fun.importRGB(160,210,237),
                               fun.importRGB(96,183,231),
                               fun.importRGB(20,147,217),
                               fun.importRGB(15,84,164)))

cols.sky <- c(fun.importRGB(195,86,26),
              fun.importRGB(201,107,24),
              fun.importRGB(221,162,19),
              fun.importRGB(242,213,9),
              fun.importRGB(251,239,212),
              fun.importRGB(252,252,252),
              fun.importRGB(152,200,236),
              fun.importRGB(96,158,214),
              fun.importRGB(77,133,193),
              fun.importRGB(68,104,171),
              fun.importRGB(62,82,150))
# functions ---------------------------------------------------------------

fun.plot.project.motif <- function(motif, input.chromVar.z = input.chromVar.jaspar.z, 
    umap.res = input.umap.res, rescale = F, cls, bks = c(-4, 0, 4), legend.pos = c(0.2, 
        0.9), ...) {
    require(scales)
    
    motif.idx <- grep(motif, rownames(input.chromVar.z))
    if (length(motif.idx) == 0) {
        message(motif, " is not found!")
        return()
    }
    
    motif.z <- input.chromVar.z[motif.idx[1], ]
    if (rescale) {
        sc <- max(abs(quantile(motif.z, probs = c(0.05, 0.95))))
        motif.z[motif.z > sc] <- sc
        motif.z[motif.z < -sc] <- -sc
    }
    
    motif.z <- motif.z %>% as.data.frame() %>% rownames_to_column("barcodes")
    colnames(motif.z)[2] <- "zval"
    
    tmp <- umap.res %>% right_join(motif.z)
    tmp <- tmp %>% mutate(zval = ifelse(zval > bks[3], bks[3], ifelse(zval < bks[1], 
        bks[1], zval)))
    p.default.cluster.motif <- ggplot(tmp, aes(UMAP1, UMAP2)) + geom_point(aes(colour = zval), 
        shape = 16, ...) + ggtitle(rownames(input.chromVar.jaspar.z)[motif.idx]) + 
        scale_color_gradientn(colours = cls, breaks = bks) + theme_pubr() + theme(text = element_blank(), 
        axis.ticks = element_blank(), legend.position = legend.pos, legend.direction = "horizontal", 
        legend.key.width = unit(2, "mm"), legend.key.height = unit(0.1, "inches"), 
        legend.text = element_text(size = 10, family = "Arial"), legend.background = element_rect(fill = "transparent", 
            colour = "transparent"), plot.margin = unit(c(1, 1, -1, -1), "mm"))
    
    p.default.cluster.motif
}

Modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}

fun.plt.motif.pt <- function(x="alpha",selected.motif=NULL,
                             showname=F,useModes=T,
                            motif.pt=output.motif.pt){

  if(useModes){
    pd.anno <-motif.pt[type==x]%>%
      group_by(pt)%>%
      summarise(subtype=Modes(as.numeric(subtype))[1])%>%
      mutate(subtype=as.factor(subtype))%>%
      as.data.frame()%>%
      column_to_rownames("pt")
  } else{
    pd.anno <-motif.pt[type==x]%>%
      mutate(subtype=as.factor(subtype))%>%
      as.data.frame()%>%
      column_to_rownames("pt")
  } 
  
  motif.pt.wd <- dcast(motif.pt[type==x],motif.name~pt,value.var = "zval",fun.aggregate = mean)
  #motif.pt <- merge(motif.pt,input.umap.res)
  pd <- motif.pt.wd %>% 
    as.data.frame()%>%
    column_to_rownames("motif.name")
  if(length(selected.motif)>1) pd <- pd[selected.motif,]
  
  pd.scaled <-  t(apply(pd,1,function(x) (x-min(x))/(max(x)-min(x))))
  if(showname){
    pheatmap(pd.scaled,cluster_cols = F,show_rownames = T,fontsize_row = 5,
             show_colnames = F,main=x,
             annotation_col = pd.anno,color = cols.hm.zval.fun(50))
  }else{
    pheatmap(pd.scaled,cluster_cols = F,show_rownames = F,show_colnames = F,main=x,
             annotation_col = pd.anno,color = cols.hm.zval.fun(50))
    
  }
  
}

fun.gamSmooth <- function(tc,new.tps = seq(0,20,by=0.05)){
  require(mgcv)
  tc.fit <- gam(zval~s(pt,bs="cs"),data = tc)
  predict(tc.fit,newdata = data.frame(pt=new.tps))
}

plotLegend <- function(cols,bks,fnames){
  require(RColorBrewer)
  #bks <- seq(round(min(rsums)),round(max(rsums))+1,length.out = 6)
  #cols<- colorRampPalette(c( "white", "blueviolet"))(5)
  #fnames<-'tmp.eps'
  setEPS()
  postscript(fnames,onefile = F,width = 0.1,height = .1*2*length(bks))
  par(mar=c(0.02, 0.04, 0.04,0.1))
  barplot(rep(1,length(bks)-1),width=diff(bks),space = 0,border = NA,
          col = cols,axes = F,horiz = T,
          xaxs = "i", yaxs = "i",xlim = c(0,1),
          ylim=c(0,sum(diff(bks)))
          #xlim=c(0,1)
  )
  axis(4,labels = F,tck = 0.2,at=cumsum(diff(bks)))
  box()
  dev.off()
}

#load data ------------------------------------------------------------

input.chromVar.res.list <- readRDS(file = "../dat/1910_v2/output.jaspar.dev.res.Rdata")
input.chromVar.jaspar.z <- assays(input.chromVar.res.list$dev)$z
input.umap.res <- fread('../dat/1910_v2/islet.cluster_labels.filt.txt',header = T)%>%separate(cluster_name, into = c("cell_type_overall", "subtype"),remove = F)%>% rename(barcodes = "index", cluster = "cluster_name")

