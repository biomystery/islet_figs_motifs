source("./libs.R")


# LOAD data  -------------------------------------------------------------

dat.prom <- rbind(data.frame(fread("./dat/alpha.100_bin_pseudotime_raw.promoter.txt")%>%
                               add_column(bin=1:100,.before = 1)%>%
  gather(key="gene",value = "percent",-1), celltype='alpha'),
  data.frame(data.frame(fread("./dat/beta.100_bin_pseudotime_raw.promoter.txt")%>%
                          add_column(bin=1:100,.before = 1)%>%
                          gather(key="gene",value = "percent",-1), celltype='beta')))
  
bin.threshold.alpha <- round(1.4265/1.9208*100)
bin.threshold.beta <- round((10.3107-9.39)/(11.27-9.39)*100)

# add subtype label
dat.prom <- (dat.prom%>%
  mutate(subtype=ifelse(celltype=="alpha" & bin <=bin.threshold.alpha, 1,
                        ifelse(celltype=="alpha" & bin >bin.threshold.alpha,2,
                               ifelse(bin <=bin.threshold.beta,1,2)))))


# explore data ------------------------------------------------------------
dat.prom.max <- dat.prom%>%
  group_by(gene,celltype)%>%
  summarise(max.percent= max(percent))

p <- ggplot(dat.prom.max,aes(max.percent))+
  geom_density()+
  theme_bw()
p+geom_vline(xintercept = quantile(dat.prom.max$max.percent),linetype=2,color='red')

p+ facet_wrap(~celltype)  

require(ggpubr)

p<-ggscatterhist(
  dat.prom.max%>%
    spread(key = "celltype",value = "max.percent"), x = "alpha", y = "beta", 
  size=1,alpha = 0.6
  #color = "Species", # comment out this and last line to remove the split by species
  #margin.plot = "histogram", # I'd suggest removing this line to get density plots
  #margin.params = list(fill = "Species", color = "black", size = 0.2)
)

# filter ------------------------------------------------------------------
dat.prom.keep <- dat.prom%>%
  group_by(gene,celltype)%>%
  summarise(max.p = max(percent))%>%
  filter(max.p>10)

dat.prom.filtered <-dat.prom %>% semi_join(dat.prom.keep,by=c("gene","celltype"))

# ttest  ------------------------------------------------------------------

fun.ttest.gene.cellsubtypes <- function(ct="alpha",g="NFIL3",...){
  test.dat <- dat.prom.filtered %>%
    filter(gene==g,celltype==ct)
  test.res <-t.test(test.dat%>%filter(subtype==1)%>%
                      select(percent),test.dat%>%filter(subtype==2)%>%
                      select(percent),...)
  return(c(mean_x=test.res$estimate[1],
           mean_y=test.res$estimate[2],
           pval=test.res$p.value/2))
}

i<- 1
require(parallel)

a.genes <- unique((dat.prom.filtered%>%filter(celltype=="alpha"))$gene)

system.time(alpha.ttest.res <- mclapply(a.genes, function(x){
  fun.ttest.gene.cellsubtypes(g = x)
},mc.cores = 10))

alpha.ttest.res <- lapply(a.genes, function(x){
  if(mod(i,100)==0) print(paste(i,x))
  i<<- i+1
  fun.ttest.gene.cellsubtypes(g = x)
})

alpha.ttest.res <- do.call(rbind, alpha.ttest.res)
colnames(alpha.ttest.res) <- c("mean.x","mean.y","pval")
rownames(alpha.ttest.res) <- a.genes
tmp <-alpha.ttest.res
alpha.ttest.res <-alpha.ttest.res %>%
  as.data.frame()


alpha.ttest.res <- (alpha.ttest.res)%>%as.data.frame()%>% rownames_to_column("gene")%>%
  mutate(FDR=p.adjust(tmp[,"pval"], method ="BH", n = nrow(tmp)),
         delta=mean.x - mean.y)

fwrite(alpha.ttest.res,file = "./dat/alpha.ttest.gene.csv")

min.Fdr <- 0.05;min.delta=5

ggplot(alpha.ttest.res,aes(mean.x,delta))+
  geom_point(color=ifelse(alpha.ttest.res$FDR<= min.Fdr,"red","grey"),shape=16)+
  theme_bw()+
  geom_hline(yintercept = 0,linetype=2,color='black')+
  #geom_vline(xintercept = c(-5,5,0),linetype=2,color='grey')+
  #scale_y_log10()+
  scale_color_manual(values = c("red","grey"))+
  theme(legend.position = "none") +  
  geom_text(data = data.frame(mean.x=30,
                              delta=c(-15,5),
                              z=c(paste0("n=",sum(alpha.ttest.res$FDR<0.05 &  alpha.ttest.res$delta<0)),
                                  paste0("n=",sum(alpha.ttest.res$FDR<0.05 &  alpha.ttest.res$delta>0)))
  ),aes(label=z))

ggsave(filename = "subfig2s_gene_alpha_ttest.pdf",useDingbats=F)


# beta test ---------------------------------------------------------------

b.genes <- unique((dat.prom.filtered%>%filter(celltype=="beta"))$gene)

system.time(beta.ttest.res <- mclapply(b.genes, function(x){
  fun.ttest.gene.cellsubtypes(g = x,ct = "beta")
},mc.cores = 10))


beta.ttest.res <- do.call(rbind, beta.ttest.res)
colnames(beta.ttest.res) <- c("mean.x","mean.y","pval")
rownames(beta.ttest.res) <- b.genes
tmp <-beta.ttest.res

beta.ttest.res <- (beta.ttest.res)%>%as.data.frame()%>% rownames_to_column("gene")%>%
  mutate(FDR=p.adjust(tmp[,"pval"], method ="BH", n = nrow(tmp)),
         delta=mean.x - mean.y)

fwrite(beta.ttest.res,file = "./dat/beta.ttest.gene.csv")

min.Fdr <- 0.05;min.delta=5

ggplot(beta.ttest.res,aes(mean.x,delta))+
  geom_point(color=ifelse(beta.ttest.res$FDR<= min.Fdr,"red","grey"),shape=16)+
  theme_bw()+
  geom_hline(yintercept = 0,linetype=2,color='black')+
  #geom_vline(xintercept = c(-5,5,0),linetype=2,color='grey')+
  #scale_y_log10()+
  scale_color_manual(values = c("red","grey"))+
  theme(legend.position = "none") +  
  geom_text(data = data.frame(mean.x=30,
                              delta=c(-15,5),
                              z=c(paste0("n=",sum(beta.ttest.res$FDR<0.05 &  beta.ttest.res$delta<0)),
                                  paste0("n=",sum(beta.ttest.res$FDR<0.05 &  beta.ttest.res$delta>0)))
  ),aes(label=z))

ggsave(filename = "subfig2s_gene_beta_ttest.pdf",useDingbats=F)


# genes final  ------------------------------------------------------------
prom.DEG.res <- rbind(data.frame(alpha.ttest.res%>% filter(FDR <0.05),celltype="alpha")
  ,data.frame(beta.ttest.res%>% filter(FDR <0.05),celltype="beta"))

fwrite(prom.DEG.res,"prom.DEG.res.csv")
system("open prom.DEG.res.csv")


# Prom_hm -----------------------------------------------------------------
dat.prom.DEGs <- semi_join(dat.prom.filtered,prom.DEG.res,by=c("gene","celltype"))
ct <- "alpha"
system.time(dat.prom.DEGs.ct <- mclapply(unique(subset(dat.prom.DEGs,celltype==ct)$gene),
                           function(g){data.frame(bin=1:100,
                                                  gene=g,
                                                  percent=predict(loess(percent~bin, data=dat.prom.DEGs%>%filter(gene==g & celltype==ct))),
                                                  celltype=ct)},mc.cores=10))
dat.prom.DEGs.smooth <- do.call(rbind,dat.prom.DEGs.ct)
ct <- "beta"
system.time(dat.prom.DEGs.ct <- mclapply(unique(subset(dat.prom.DEGs,celltype==ct)$gene),
                                         function(g){data.frame(bin=1:100,
                                                                gene=g,
                                                                percent=predict(loess(percent~bin, data=dat.prom.DEGs%>%filter(gene==g & celltype==ct))),
                                                                celltype=ct)},mc.cores=10))
dat.prom.DEGs.smooth <- rbind(dat.prom.DEGs.smooth,do.call(rbind,dat.prom.DEGs.ct))
fwrite(dat.prom.DEGs.smooth,"dat.prom.DEGs.smooth.csv")

dat.prom.DEGs.smooth <- fread('dat.prom.DEGs.smooth.csv')

# Subfig2D_beta_prom_hm ---------------------------------------------------
pdf(file = "subfig2D_DEG_prom_hm.pdf")
for(ct in c("alpha","beta")){
  pd <- dat.prom.DEGs.smooth%>%
    filter(celltype==ct)%>%
    select(-celltype)%>%
    spread(key = bin,value = percent)%>%
    column_to_rownames("gene")
  
  pd.scale <- t(apply(pd, 1, function(x) (x-min(x))/(max(x)-min(x))))
  fwrite(x = as.data.frame(pd.scale), file = paste0(ct,".smoothed.prom.scaled.csv"),row.names = T,col.names = T)
  fwrite(as.data.frame(pd),file = paste0(ct,".smoothed.prom.raw.csv"),row.names = T,col.names = T)
  #pheatmap(pd.scale,cluster_cols = F,scale = "none",show_rownames = F,show_colnames = F,
  #         cutree_rows =9,# ifelse(ct=="beta",9,6),
  #         main = paste(ct,"DEGs' promoter pesudo-time-course"))
}
dev.off()
