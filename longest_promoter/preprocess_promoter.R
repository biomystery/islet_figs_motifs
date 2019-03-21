source("./chromVar/libs.R")

# get tfclss  -------------------------------------------------------------
tfclass.db <- readRDS("../../atacMotif/db/tfclass.rds")
tfclass.dic <- readRDS("../../atacMotif/db/dic_jaspar_tfclass.rds")

get_subfamily <- function(g="RFX3",db=tfclass.db){
  idx <- which(toupper(tfclass.db$merge$tf.symbol)==g)
  sf.name <- tfclass.db$merge$subfamily.name[idx]
  sf.id <- tfclass.db$merge$subfamily.id[idx]
  tfclass.db$merge%>%
    filter(subfamily.id == sf.id)
}
get_promoter_dat <- function(tf_res,prom=dat){
  g <- unique(unlist(sapply(toupper(tf_res$tf.symbol), function(x) colnames(prom)[colnames(prom)%in% x])))
  prom[,g,with=F]
}


plot_tf_prom <- function(tf,add_point=F,free_y=F){
  if(length(tf)>1){# use name directly 
    sub_f<- list()
    sub_f$tf.symbol <- tf
  }else{# get sub_family 
    sub_f <- get_subfamily(tf)
  }

  g.dat <- get_promoter_dat(sub_f)%>%
    add_column(x=1:100,.before = 1)%>%
    gather(key = "gene",value = "percent",-1)
  
  r.tf <- setdiff(sub_f$tf.symbol,unique(g.dat$gene))
  if(length(r.tf)>0)
    g.dat <- rbind(g.dat,
                   do.call(rbind,lapply(r.tf, function(x) data.frame(x=1:100,    gene=x,     percent=0))))
  require(ggplot2)
  p<- ggplot(g.dat,aes(x,percent))+
    geom_smooth(aes(color=gene))+
    coord_cartesian(expand = T)+
    theme_classic()+
    facet_wrap(~gene)+
    theme(legend.position = "none")+
    ggtitle(tf)
  if(add_point) p<- p+geom_point()
  if(free_y) p<- p+facet_wrap(~gene,scales = "free_y") else p <-p+facet_wrap(~gene)
  return(list(plt=p,sf=sub_f))
}


# generate motif pseudo time bin  -----------------------------------------
output.motif.pt.bin<-rbind(output.motif.pt%>%
                             filter(type=='alpha')%>%
                             mutate(pt_bin=as.numeric(cut(pt,breaks = seq(0-0.00001,pt.max.a+0.00001,length.out = 100)))),
                           output.motif.pt%>%
                             filter(type=='beta')%>%
                             mutate(pt_bin=as.numeric(cut(pt,breaks = seq(0-0.00001,pt.max.b+0.00001,length.out = 100)))))%>%
  group_by(motif.name,type,pt_bin)%>%
  summarise(zval_bin= mean(zval))
output.motif.pt.bin.wc <- output.motif.pt.bin%>%
  spread(key = pt_bin,
         value = zval_bin)
pd.a <- t(apply((output.motif.pt.bin.wc%>% filter(type=="alpha"))[,-c(1,2)],1,
              function(x) (x-min(x))/(max(x)-min(x))))
rownames(pd.a) <-(output.motif.pt.bin.wc%>% filter(type=="alpha"))$motif.name
pheatmap(pd.a,
         cluster_cols =  F)
# alpha data  -------------------------------------------------------------

fn <- "./dat/alpha.100_bin_pseudotime_raw.promoter.txt"

dat <- fread(fn)
for(x in c("nop","wp")){
  pdf(file = paste0("alpha.promoter.",x,".pdf"))
  a<- ifelse(x=="nop",F,T)
  #res <- plot_tf_prom(c("GCG","CCKAR","MAPK10","EYS","AQPEP","TNFSF9",
  #                      "ACHE","MARCKSL1","ALOXE3","FOSL1"),free_y = F,add_point = a);print(res$plt)
  res <- plot_tf_prom(c("JUNB","JUND","FOSL1","FOSL2","JUN","FOS"),free_y = F);print(res$plt)
  res <- plot_tf_prom("RFX3",free_y = F,add_point = a);print(res$plt)
  res <- plot_tf_prom("NEUROD2",free_y = F,add_point = a);print(res$plt)
  res <- plot_tf_prom("STAT1",free_y = F,add_point = a);print(res$plt)
  res <- plot_tf_prom("SNAI2",free_y = F,add_point = a);print(res$plt)
  res <- plot_tf_prom("IRF2",free_y = F,add_point = a);print(res$plt)
  res <- plot_tf_prom("FOXP2",free_y = F,add_point = a);print(res$plt)
  res <- plot_tf_prom("FOXO4",free_y = F,add_point = a);print(res$plt)
  res <- plot_tf_prom(c("GATA1","GATA2","GATA3","GATA4","GATA5","GATA6"),free_y = F,add_point = a);print(res$plt)
  dev.off()
}

# subfig2E_alpha_prom_1 ---------------------------------------------------
gs<- get_subfamily(g = "NEUROD2",db = tfclass.db)
gs$genus.name <- sub("NGN[1-3]\\/","",gs$genus.name)
res <- plot_tf_prom(tf = gs$genus.name,free_y = F,add_point = a);print(res$plt)
new_res <- lapply(gs$genus.name ,function(g)
  data.frame(x=1:100,g=g,percent=predict(loess(percent~x, data=res$plt$data%>%filter(gene==g)))))
new_res <- do.call(rbind,new_res)%>%
  spread(key = x,value = percent)%>%
  column_to_rownames("g")



pdf(file = "./figs/subfig2e_alpha_prom_1.pdf")
print(res$plt)
#bks <- c(0,5,seq(10,,length.out = 47),36)
bks <- c(0,seq(3,30,length.out = 48),36)
p<- pheatmap(new_res,scale = "none",
         breaks = bks,
         cluster_cols = F,
         cluster_rows = T,
         #annotation_row = data.frame(apply(new_res, 1, max)),
         show_colnames = F,
         border_color = NA,
         cellheight = 10,cellwidth = 2,
         color = colorRampPalette(c("white",brewer.pal(9,"Reds")))(51))
dev.off()
system("open ./figs/subfig2e_alpha_prom_1.pdf")
if(T){
  setEPS()
  postscript(file = "./figs/subfig2e_alpha_prom_1.eps",width = 2,height =1.5)
  #bks <- c(0,5,seq(10,,length.out = 47),36)
  bks <- c(0,seq(3,30,length.out = 48),36)
  pheatmap(new_res[p$tree_row$order,],scale = "none",
           breaks = bks,
           legend = F,
           cluster_cols = F,
           cluster_rows = F,
           #annotation_row = data.frame(apply(new_res, 1, max)),
           show_colnames = F,
           show_rownames = F,
           border_color = NA,
           cellheight = 5,cellwidth =1,
           color = colorRampPalette(c("white",brewer.pal(9,"Reds")))(51))
  dev.off()
}


sink("./figs/subfig2e_alpha_prom_1.txt")
print(bks)
do.call(rbind,
        apply(new_res, 1, function(x)data.frame(min=min(x),max=max(x),rg=max(x)-min(x))))
sink()


# subfig2E_alpha_prom_2 ---------------------------------------------------
res <- plot_tf_prom(c("JUN","JUNB","JUND","FOS","FOSB","FOSL1","FOSL2"),free_y = F,add_point = a);print(res$plt)
new_res <- lapply(c("FOSL2","FOSL1","FOSB","FOS","JUN","JUND","JUNB"),function(g)
  data.frame(x=1:100,g=g,percent=predict(loess(percent~x, data=res$plt$data%>%filter(gene==g)))))
new_res <- do.call(rbind,new_res)%>%
  spread(key = x,value = percent)%>%
  column_to_rownames("g")

fwrite(new_res,file = "/Users/frank/Dropbox (UCSD_Epigenomics)/projects/islet/slides/2019-03-18_dat_figs/fig2E_fos_jun_alpha_smoothed.ps.csv",row.names  = T)
new_res_scale <- t(apply(new_res,1,function(x) (x-min(x))/(max(x)-min(x))))

if(T){
  pdf(file = "./figs/subfig2e_alpha_prom_2.pdf")
  print(res$plt)
  bks <- c(0,9.9,seq(10,45,length.out = 47),50)
  pheatmap(new_res,scale = "none",
           breaks = bks,
           cluster_cols = F,
           cluster_rows = F,
           #annotation_row = data.frame(apply(new_res, 1, max)),
           show_colnames = F,
           border_color = NA,
           cellheight = 10,cellwidth = 2,
           color = colorRampPalette(c("white",brewer.pal(9,"Reds")))(51))
  pheatmap(new_res_scale,scale = "none",
           cluster_cols = F,
           cluster_rows = T,
           annotation_row = data.frame(apply(new_res, 1, max)),
           show_colnames = F,
           border_color = NA,
           cellheight = 10,cellwidth = 1.5,
           color = colorRampPalette(c("white",brewer.pal(9,"Reds")))(51)) 
  dev.off()
}
sink("./figs/subfig2e_alpha_prom_2.txt")
print(bks)
do.call(rbind,
        apply(new_res, 1, function(x)data.frame(min=min(x),max=max(x),rg=max(x)-min(x))))
sink()


# beta explore--------------------------------------------------------------------

fn <- "./dat/beta.100_bin_pseudotime_raw.promoter.txt"

dat <- fread(fn)
for(x in c("nop","wp")){
  pdf(file = paste0("beta.promoter.",x,".pdf"))
  a<- ifelse(x=="nop",F,T)
  #res <- plot_tf_prom(c("GCG","CCKAR","MAPK10","EYS","AQPEP","TNFSF9",
  #                      "ACHE","MARCKSL1","ALOXE3","FOSL1"),free_y = F,add_point = a);print(res$plt)
  res <- plot_tf_prom(c("JUN","JUNB","JUND","FOS","FOSB","FOSL1","FOSL2"),free_y = F,add_point = a);print(res$plt)
  res <- plot_tf_prom("RFX3",free_y = F,add_point = a);print(res$plt)
  res <- plot_tf_prom("NRL",free_y = F,add_point = a);print(res$plt)
  res <- plot_tf_prom("CTCF",free_y = F,add_point = a);print(res$plt)
  res <- plot_tf_prom("FOXA1",free_y = F,add_point = a);print(res$plt)
  res <- plot_tf_prom("STAT1",free_y = F,add_point = a);print(res$plt)
  res <- plot_tf_prom("NEUROD2",free_y = F,add_point = a);print(res$plt)
  res <- plot_tf_prom("NFATC2",free_y = F,add_point = a);print(res$plt)
  dev.off()
}


# subfig2E_beta_motif_prom_2 ------------------------------------------------
get_subfamily(g = "FRA1",db = tfclass.db)
get_subfamily(g = "JUNB",db = tfclass.db)

res <- plot_tf_prom(c("JUN","JUNB","JUND","FOS","FOSB","FOSL1","FOSL2"),free_y = F,add_point = a);print(res$plt)
new_res <- lapply(c("FOSL2","FOSL1","FOSB","FOS","JUN","JUND","JUNB"),function(g)
  data.frame(x=1:100,g=g,percent=predict(loess(percent~x, data=res$plt$data%>%filter(gene==g)))))
new_res <- do.call(rbind,new_res)%>%
  spread(key = x,value = percent)%>%
  column_to_rownames("g")

fwrite(new_res,file = "/Users/frank/Dropbox (UCSD_Epigenomics)/projects/islet/slides/2019-03-18_dat_figs/fig2E_fos_jun_beta_smoothed.ps.csv",row.names = T)

new_res_scale <- t(apply(new_res,1,function(x) (x-min(x))/(max(x)-min(x))))

pdf(file = "./figs/subfig2e_beta_prom_2.pdf")
print(res$plt)
bks <- c(0,9.9,seq(10,30,length.out = 47),40)
pheatmap(new_res,scale = "none",
         breaks = bks,
         cluster_cols = F,
         cluster_rows = F,
         #annotation_row = data.frame(apply(new_res, 1, max)),
         show_colnames = F,
         border_color = NA,
         cellheight = 10,cellwidth = 2,
         color = colorRampPalette(c("white",brewer.pal(9,"Reds")))(51))
pheatmap(new_res_scale,scale = "none",
         cluster_cols = F,
         cluster_rows = T,
         annotation_row = data.frame(apply(new_res, 1, max)),
         show_colnames = F,
         border_color = NA,
         cellheight = 10,cellwidth = 1.5,
         color = colorRampPalette(c("white",brewer.pal(9,"Reds")))(51)) 
dev.off()
sink("./figs/subfig2e_beta_prom_2.txt")
print(bks)
do.call(rbind,
        apply(new_res, 1, function(x)data.frame(min=min(x),max=max(x),rg=max(x)-min(x))))
sink()

## raw 
if(F){
  raw_res <- res$plt$data%>%
    spread(x,percent)%>%
    column_to_rownames("gene")
  
  pheatmap(raw_res,scale = "none",
           cluster_cols = F,
           cluster_rows = T,
           show_colnames = F,
           border_color = NA,
           cellheight = 10,
           color = colorRampPalette(c("white",brewer.pal(9,"Reds")))(21))                       
  
  raw_res_scale <- t(apply(raw_res,1,function(x) x/max(x)))
  
  pheatmap(raw_res_scale,scale = "none",
           cluster_cols = F,
           cluster_rows = T,
           annotation_row = data.frame(apply(raw_res, 1, max)),
           show_colnames = F,
           border_color = NA,
           cellheight = 10,
           color = colorRampPalette(c("white",brewer.pal(9,"Reds")))(21)) 
}




# subfig2E_beta_motif_prom_1 ------------------------------------------------

gs <- get_subfamily(g = "RFX3",db = tfclass.db)
res <- plot_tf_prom("RFX3",free_y = F,add_point = a);print(res$plt)

new_res <- lapply(paste0("RFX",seq(8,1)),function(g)
  data.frame(x=1:100,g=g,percent=predict(loess(percent~x, data=res$plt$data%>%filter(gene==g)))))
new_res <- do.call(rbind,new_res)%>%
  spread(key = x,value = percent)%>%
  column_to_rownames("g")

new_res_scale <- t(apply(new_res,1,function(x) (x-min(x))/(max(x)-min(x))))
pdf(file = "./figs/subfig2e_beta_prom_1.pdf")
print(res$plt)
bks <- c(0,seq(20,30,length.out = 48),40)
pheatmap(new_res,scale = "none",
         breaks = bks,
         cluster_cols = F,
         cluster_rows = F,
         #annotation_row = data.frame(apply(new_res, 1, max)),
         show_colnames = F,
         border_color = NA,
         cellheight = 10,cellwidth = 2,
         color = colorRampPalette(c("white",brewer.pal(9,"Reds")))(51))

pheatmap(new_res_scale[c(3,7,8,6,2,5,1,4),],scale = "none",
         cluster_cols = F,
         cluster_rows = F,
         annotation_row = data.frame(apply(new_res, 1, max)),
         show_colnames = F,
         border_color = NA,
         cellheight = 10,
         color = colorRampPalette(c("white",brewer.pal(9,"Reds")))(51)) 
dev.off()
sink("./figs/subfig2e_beta_prom_1.txt")
print(bks)
do.call(rbind,
        apply(new_res, 1, function(x)data.frame(min=min(x),max=max(x),rg=max(x)-min(x))))
sink()
## raw 
if(F){
  raw_res <- res$plt$data%>%
    spread(x,percent)%>%
    column_to_rownames("gene")
  
  pheatmap(raw_res,scale = "none",
           cluster_cols = F,
           cluster_rows = T,
           show_colnames = F,
           border_color = NA,
           cellheight = 10,
           color = colorRampPalette(c("white",brewer.pal(9,"Reds")))(21))                       
  
  raw_res_scale <- t(apply(raw_res,1,function(x) x/max(x)))
  
  pheatmap(raw_res_scale,scale = "none",
           cluster_cols = F,
           cluster_rows = T,
           annotation_row = data.frame(apply(raw_res, 1, max)),
           show_colnames = F,
           border_color = NA,
           cellheight = 10,
           color = colorRampPalette(c("white",brewer.pal(9,"Reds")))(21)) 
}



