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


# variablity for alpha --------------------------------------------

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




# variablity  for beta  ----------------------------------------------
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
  
  
}



# aggregate results  ------------------------------------------------------


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





# mannual calc ------------------------------------------------------------

output.jaspar.z  <- merge(input.chromVar.jaspar.z,input.umap.res)
pd.variablility.mannual <- output.jaspar.z%>%
  group_by(cell_type_overall,motif.name)%>%
  summarise(var= var(zval),
            sd = sd(zval))%>%
  group_by(cell_type_overall)%>%
  arrange(var)%>%
  mutate(ord.var = order(var,decreasing = T),
         ord.sd = order(sd,decreasing = T))
  
ggplot(pd.variablility.mannual,
       aes(ord.sd,sd)) + geom_point(aes(colour=cell_type_overall))+
  geom_hline(yintercept = 1.2)+ theme_light()+ scale_x_log10()


# compare subtypes --------------------------------------------------------

fun.ttest.motif.cellsubtypes <- function(celltype="alpha",motif="NFIL3",...){
  test.dat <- output.jaspar.z%>%
    filter(subtype!="3")%>%
    filter(motif.name==motif,cell_type_overall==celltype)
  test.res <-t.test(test.dat%>%filter(subtype==1)%>%
                      select(zval),test.dat%>%filter(subtype==2)%>%
                      select(zval),...)
  return(c(mean_x=test.res$estimate[1],
           mean_y=test.res$estimate[2],
           pval=test.res$p.value/2))
  }
i<- 1
require(parallel)
parallel::mclapply()


alpha.ttest.res <- mclapply(unique(output.jaspar.z$motif.name), function(x){
  print(paste(i,x))
  i<<- i+1
  fun.ttest.motif.cellsubtypes(motif = x)
  },mc.cores = 10)

alpha.ttest.res <- do.call(rbind, alpha.ttest.res)
colnames(alpha.ttest.res) <- c("mean.x","mean.y","pval")
rownames(alpha.ttest.res) <- unique(output.jaspar.z$motif.name)
alpha.ttest.res <-alpha.ttest.res %>%
  as.data.frame()%>%
  mutate(FDR=p.adjust(tmp[,"pval"], method ="BH", n = nrow(tmp)),
         delta=mean.x - mean.y)

alpha.ttest.res <- alpha.ttest.res%>% 
  rownames_to_column("motif")%>%
  arrange(desc(delta))

fwrite(alpha.ttest.res,file = "alpha.ttest.res.csv")
alpha.ttest.res <- fread("alpha.ttest.res.csv")

alpha.ttest.res.pd <- alpha.ttest.res
min.Fdr <- min(alpha.ttest.res$FDR[alpha.ttest.res$FDR!=0])
alpha.ttest.res.pd$FDR[alpha.ttest.res.pd$FDR<min.Fdr]<- min.Fdr


ggplot(alpha.ttest.res.pd,aes(delta,-log10(FDR)))+
  geom_point(aes(color=ifelse(abs(delta)<0.5,"red","grey")
                 ),shape=ifelse(alpha.ttest.res$FDR< min.Fdr,17,16))+
  theme_light()+
  geom_hline(yintercept = 2,linetype=2,color='grey')+
  geom_vline(xintercept = c(-.5,.5),linetype=2,color='grey')+
  scale_y_log10()+
  scale_color_manual(values = c("red","grey"))+
  theme(legend.position = "none")+
  geom_text(data = alpha.ttest.res.pd[alpha.ttest.res$FDR< min.Fdr,],
           aes(label=motif),nudge_x = 0,nudge_y = -0.05,angle=-90,size=3,hjust=0)


table(tmp$FDR<0.01)
sum(tmp$mean.x-tmp$mean.y>0.5)
sum(tmp$mean.x-tmp$mean.y< -0.5)



# beta --------------------------------------------------------------------


system.time(beta.ttest.res <- mclapply(unique(output.jaspar.z$motif.name), function(x){
  fun.ttest.motif.cellsubtypes(motif = x,celltype = "beta")
},mc.cores = 10))

beta.ttest.res <- do.call(rbind, beta.ttest.res)
colnames(beta.ttest.res) <- c("mean.x","mean.y","pval")

beta.ttest.res <-beta.ttest.res %>%
  as.data.frame()%>%
  mutate(FDR=p.adjust(beta.ttest.res[,"pval"], method ="BH", n = nrow(beta.ttest.res)),
         delta=mean.x - mean.y)

rownames(beta.ttest.res) <- unique(output.jaspar.z$motif.name)
beta.ttest.res <- beta.ttest.res%>% 
  rownames_to_column("motif")%>%
  arrange(delta)

fwrite(beta.ttest.res,file = "beta.ttest.res.csv")


beta.ttest.res.pd <- beta.ttest.res
min.Fdr <- min(beta.ttest.res.pd$FDR[beta.ttest.res.pd$FDR!=0])
beta.ttest.res.pd$FDR[beta.ttest.res.pd$FDR<min.Fdr]<- min.Fdr


ggplot(beta.ttest.res.pd,aes(delta,-log10(FDR)))+
  geom_point(aes(color=ifelse(abs(delta)<0.5,"red","grey")
  ),shape=ifelse(beta.ttest.res$FDR< min.Fdr,17,16))+
  theme_light()+
  geom_hline(yintercept = 2,linetype=2,color='grey')+
  geom_vline(xintercept = c(-.5,.5),linetype=2,color='grey')+
  scale_y_log10()+
  scale_color_manual(values = c("red","grey"))+
  theme(legend.position = "none")+
  geom_text(data = beta.ttest.res.pd[beta.ttest.res$FDR< min.Fdr,],
            aes(label=motif),nudge_x = 0,nudge_y = -0.05,angle=-90,size=3,hjust=0)

motif.list <- beta.ttest.res.pd[beta.ttest.res$FDR< min.Fdr,]$motif
plot.ggviolin(motif.list,
              celltype = "beta")+ 
  facet_wrap(~motif.name)+
  theme(legend.position = "none")

p.default.cluster.motifs <- lapply(motif.list, function(x) fun.plot.project.motif(motif = x,umap.res = input.umap.res.sub))
ggarrange(plotlist=p.default.cluster.motifs,ncol = 3,nrow = 4)


# delta -------------------------------------------------------------------
system.time(delta.ttest.res <- mclapply(unique(output.jaspar.z$motif.name), function(x){
  fun.ttest.motif.cellsubtypes(motif = x,celltype = "delta")
},mc.cores = 10))

delta.ttest.res <- do.call(rbind, delta.ttest.res)
colnames(delta.ttest.res) <- c("mean.x","mean.y","pval")

delta.ttest.res <-delta.ttest.res %>%
  as.data.frame()%>%
  mutate(FDR=p.adjust(delta.ttest.res[,"pval"], method ="BH", n = nrow(delta.ttest.res)),
         delta=mean.x - mean.y)


rownames(delta.ttest.res) <- unique(output.jaspar.z$motif.name)
delta.ttest.res <- delta.ttest.res%>% 
  rownames_to_column("motif")%>%
  arrange(delta)

fwrite(delta.ttest.res,file = "delta.ttest.res.csv")


delta.ttest.res.pd <- delta.ttest.res
min.Fdr <- min(delta.ttest.res.pd$FDR[delta.ttest.res.pd$FDR!=0])
delta.ttest.res.pd$FDR[delta.ttest.res.pd$FDR<min.Fdr]<- min.Fdr
delta.ttest.res.pd<- delta.ttest.res.pd%>%
  mutate(selected=ifelse(abs(delta)>0.5 & FDR<0.01,T,F))
ggplot(delta.ttest.res.pd,
       aes(delta,-log10(FDR)))+
  geom_point(aes(color=selected),
             shape=ifelse(delta.ttest.res$FDR< min.Fdr,17,16))+
  theme_light()+
  geom_hline(yintercept = 2,linetype=2,color='grey')+
  geom_vline(xintercept = c(-.5,.5),linetype=2,color='grey')+
  scale_y_log10()+
  scale_color_manual(values = c("grey","red"))+
  theme(legend.position = "none")+
  geom_text(data = delta.ttest.res.pd[delta.ttest.res$FDR< min.Fdr,],
            aes(label=motif),nudge_x = 0,nudge_y = -0.05,angle=-90,size=3,hjust=0)

motif.list <- beta.ttest.res.pd[delta.ttest.res.pd$selected,]$motif
plot.ggviolin(motif.list[1:12],
              celltype = "beta")+ 
  facet_wrap(~motif.name)+
  theme(legend.position = "none")

p.default.cluster.motifs <- lapply(motif.list[1:12], function(x) fun.plot.project.motif(motif = x,umap.res = input.umap.res.sub))
ggarrange(plotlist=p.default.cluster.motifs,ncol = 3,nrow = 4)


# combine -----------------------------------------------------------------

ttest.res <- rbind(data.frame(alpha.ttest.res,celltype="alpha"),
                   data.frame(beta.ttest.res,celltype="beta"),
                   data.frame(delta.ttest.res,celltype="delta"))
ttest.res <- ttest.res%>% 
  group_by(celltype)%>%
  mutate(selected=ifelse(abs(delta)>0.5 & FDR<0.01,T,F))

fwrite(ttest.res,file = "ttest.res.csv")


# plot volcano  -------------------------------------------------------------------
ttest.res.pd <- ttest.res%>%
  group_by(celltype)%>%
  mutate(min.FDR = min(FDR[FDR!=0]),
         FDR.bak =FDR)%>%
  mutate(FDR=ifelse(FDR==0,min.FDR,FDR))

plt.volcano.ttest <- ggplot(ttest.res.pd,
       aes(delta,-log10(FDR)))+
  geom_point(aes(color=selected),
             shape=ifelse(ttest.res.pd$FDR< ttest.res.pd$min.FDR,17,16))+
  theme_light()+
  geom_hline(yintercept = 2,linetype=2,color='grey')+
  geom_vline(xintercept = c(-.5,.5),linetype=2,color='grey')+
  scale_y_log10()+
  scale_color_manual(values = c("grey","red"))+
  theme(legend.position = "none")+
  facet_wrap(~celltype)+
  xlab("\\delta Z (type 1 - type 2)")

if(T){
  png(filename ="ttest.volcano.png",width = 8,height = 3,res = 300,units = 'in')
  print(plt.volcano.ttest)
  dev.off()
  system("open ttest.volcano.png")
}



plt.volcano.ttest+ geom_text(data = ttest.res.pd%>%
            filter(FDR.bak< min.FDR),
          aes(label=motif),nudge_x = 0,nudge_y = -0.05,angle=-90,size=3,hjust=0)+
  


# plot examples violin ----------------------------------------------------

dmotifs.list <- sapply(c("alpha","beta","delta"), function(x) subset(ttest.res,celltype==x & selected)$motif)
tmp <- gplots::venn(dmotifs.list)
if(T){
  pdf(file = "venn.ttest.res.pdf",width = 4,height = 4)
  gplots::venn(dmotifs.list)
  dev.off()
  system("open venn.ttest.res.pdf")  
}
dmotifs.list.inter <- attr(tmp,"intersections")
sapply(dmotifs.list.inter,function(x) x[1])

require(ggpubr)
my_comparisons <- list(c("1","2"))
plot.ggviolin <- function(m,celltypes="alpha"){
  ggviolin(output.jaspar.z%>%
             filter(subtype!="3")%>%
             filter(motif.name%in%m,cell_type_overall%in%celltypes),
           x = "subtype",y = "zval",add = "boxplot",fill="cluster"
  )
}

plot.ggviolin(sapply(dmotifs.list.inter,function(x) x[1]),celltypes = c("alpha","beta","delta"))+ 
  facet_grid(cell_type_overall~motif.name,scales = "free_y")+
  scale_fill_manual(values = cols.celltype)+
  theme_light()
#+  stat_compare_means(comparisons = my_comparisons,method = "t.test",method.args=list(mu=0.5))

plot.ggbox<- function(m,celltypes="alpha"){
  dic<- names(m)
  names(dic)<- m
  pd <- output.jaspar.z%>%
    filter(subtype!="3")%>%
    filter(motif.name%in%m,cell_type_overall%in%celltypes)%>%
    mutate(ctype=factor(dic[motif.name],
                        levels =c("alpha","beta","delta","alpha:beta","beta:delta","alpha:delta","alpha:beta:delta")))
  
  ggboxplot(pd,
           x = "subtype",y = "zval",fill="cluster",
           outlier.shape = NA)
}

# seperate
for(s in names(dmotifs.list.inter)){
  print(s)
  x <- dmotifs.list.inter[[s]][1]
  p<- plot.ggbox(x,celltypes = c("alpha","beta","delta"))+
    facet_grid(motif.name~cell_type_overall,scales = "free_y")+
    scale_fill_manual(values = cols.celltype)+
    theme_light()+
    theme(legend.position = "none")+
    ggtitle(paste(s,x,sep = ":"))
  
  png(filename = paste0(s,"_",x,".png"),width = 5,height = 3,res = 300,units = 'in')
  print(p)
  dev.off()
  system(paste0("open ",s,"_",x,".png"))
}

# merged
ord <- c("alpha","beta","delta","alpha:beta","beta:delta","alpha:delta","alpha:beta:delta")
p.ggboxes <- lapply(ord,function(s){
  x <- dmotifs.list.inter[[s]][1]
  plot.ggbox(x,celltypes = c("alpha","beta","delta"))+
    facet_grid(motif.name~cell_type_overall,scales = "free_y")+
    scale_fill_manual(values = cols.celltype)+
    theme_light()+
    theme(legend.position = "none")+
    ggtitle(paste(s,x,sep = ":"))
})


if(T){
  png(filename = "ttest.egs.boxplot.png",height = 9,width = 4,units = 'in',res = 300)
  print(plot.ggbox(sapply(dmotifs.list.inter[ord],function(x) x[1]),celltypes = c("alpha","beta","delta"))+ 
    facet_grid(ctype~cell_type_overall,scales = "free_y")+
    scale_fill_manual(values = cols.celltype)+
    theme_light()+
      theme(legend.position = "none"))
  dev.off()
  system("open ttest.egs.boxplot.png")
}
  
#+  stat_compare_means(comparisons = my_comparisons,method = "t.test",method.args=list(mu=0.5))


# Plot_hm -----------------------------------------------------------------\
input.chromVar.res.list <- readRDS(file = "./dat/output.jaspar.dev.res.Rdata")
input.chromVar.jaspar.z <- data.table(assays(input.chromVar.res.list$dev)$z,keep.rownames = T)
input.umap.res <- fread('./dat/Islet_123.MNN_corrected.UMAP.txt',header = T)
input.chromVar.jaspar.var <- fread("./dat/output.jaspar.var.res.csv")

# filter unkonwn
input.umap.res <- input.umap.res %>% 
  separate(cluster,into = c("cell_type_overall","subtype"),remove=F)
input.umap.res[is.na(input.umap.res)]<-0

input.chromVar.jaspar.z.agg <- melt(input.chromVar.jaspar.z,
                                    id="rn",variable.name = "barcodes",value.name = "zval")

# add celltype
input.chromVar.jaspar.z.agg <- merge(input.chromVar.jaspar.z.agg,input.umap.res)

output.chromvar.jaspar.z.avg_by_subct<-  
  input.chromVar.jaspar.z.agg[cluster %in% select.celltypes,.(zval_avg=mean(zval)),by=.(rn,cluster)]%>%
  group_by(rn)%>%
  mutate(zval_avg = (zval_avg-min(zval_avg))/(max(zval_avg)-min(zval_avg)))%>%
  separate(rn,into = c("id","name"),sep = "_")%>%
  select(-one_of("id"))%>%
  spread(key = cluster,value = zval_avg)%>%
  as.data.frame()%>%
  column_to_rownames("name")


