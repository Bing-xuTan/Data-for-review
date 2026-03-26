if (T) {
  dir.create("scripts")
  dir.create("results")
  dir.create("files")
  dir.create("figures")
  dir.create("origin_datas/GEO",recursive = T)
  dir.create("origin_datas/TCGA")
}
library(stringr)
library(tidydr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(RisktypeProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggstatsplot)
library(ggcor)
library(ggstance)
options(stringsAsFactors = F)
source('base.R')
my_violin=function(dat,group,group_cols=ggsci::pal_aaas()(10),#test_method='kruskal.test',
                   fill= "Group",label=c("p.format",'p.signif')[1],
                   xlab='',ylab='',title='',x.size=10,y.size=10,legend.position='top'){
  
  # dat = tcga.b.cell$GeneSet,
  # group = tcga.subtype$Cluster
  data=data.frame(Group=group,value=dat)
  data=crbind2DataFrame(data)
  data=melt(data)
  data=data[which(!is.na(data[,1])),]
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=ggplot(data,aes(x=Group, y=value,fill=Group)) +
    geom_violin()+  
    geom_boxplot(width=0.2,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
    scale_fill_manual(values = group_cols)+
    theme_classic(base_size = 20)+labs(x=xlab,y=ylab,title = title)+
    ggpubr::stat_compare_means(aes(group=Group), label = label, method =test_method)+
    theme(legend.position = legend.position,axis.text = element_text(color = 'black'),
          title = element_text(size = 12),text = element_text(family = 'Times'),
          axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))
  return(p)
}
my_mutiboxplot=function(dat,group,group_cols=ggsci::pal_aaas()(10),
                        #test_method=c('t.test','wilcox.test','anova','kruskal.test')[4],
                        bw=T,xlab='',ylab='score',title='',size=10,angle = 45, hjust = 1,
                        legend.position='top',fill='group',notch=F){
  # dat=tcga.est[tcga.subtype.cli$Samples,]
  # group=tcga.subtype.cli$Cluster
  dat.bind=cbind(dat,Cluster=group)
  dat.bind=crbind2DataFrame(dat.bind)
  dat.melt=melt(dat.bind)
  #data=data[which(!is.na(data[,1])),]
  colnames(dat.melt)=c('Group','type','value')
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=dat.melt %>%
    ggplot(aes(x=type, y=value,fill=Group)) +
    geom_boxplot(notch = notch) +  
    scale_fill_manual(values =group_cols)+  
    ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method = test_method)+
    labs(x="", y = ylab, fill =fill,title =title) +
    #theme_light()+
    theme_bw()+
    #theme_classic()
    theme(legend.position = legend.position,                
          plot.title = element_text(hjust = 0.5),text = element_text(family = 'Times'),
          axis.text.x = element_text(size = size,angle = angle, hjust = hjust)) 
  return(p)
}


dir.create('results/01.scRNA')
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(dplyr)

dir_name=list.files('origin_datas/GEO/GSE196756/')
tissue=c('tumor','peri-tumor','tumor','peri-tumor','tumor','peri-tumor')
datalist=list()
for (i in 1:length(dir_name)){
  dir.10x = paste0("origin_datas/GEO/GSE196756/",dir_name[i])
  list.files(dir.10x)
  my.data <- Read10X(data.dir = dir.10x)
  datalist[[i]]<- CreateSeuratObject(counts=my.data,project = dir_name[i],min.cells = 3, min.features = 200)
  datalist[[i]] <- AddMetaData(datalist[[i]] , dir_name[i],col.name = "Samples")
  datalist[[i]] <- AddMetaData(datalist[[i]] , tissue[i],col.name = "Tissue")
}
names(datalist)=dir_name
rm(my.data)


for (i in 1:length(datalist)){
  sce <- datalist[[i]]
  sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
  sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")
  datalist[[i]] <- sce
  rm(sce)
}
sce <- merge(datalist[[1]],y=datalist[2:length(datalist)])
rm(datalist)

raw_meta=sce@meta.data
raw_count <- table(raw_meta$Samples)
raw_count
sum(raw_count)#  36515
pearplot_befor<-VlnPlot(sce,group.by ='Samples',
                        features = c("nFeature_RNA", "nCount_RNA","percent.mt"),
                        pt.size = 0,
                        ncol = 3)
pearplot_befor
ggsave('results/01.scRNA/pearplot_befor.pdf',pearplot_befor,height = 5,width = 10)

#过滤
sce=subset(sce, subset=nFeature_RNA>200 & nFeature_RNA<8000 & percent.mt<10)
clean_meta=sce@meta.data
clean_count <- table(clean_meta$Samples)
clean_count
sum(clean_count)#28387
pearplot_after <- VlnPlot(sce,group.by ='Samples',
                          features = c("nFeature_RNA", "nCount_RNA","percent.mt"),
                          pt.size = 0,
                          ncol = 3)
pearplot_after
ggsave('results/01.scRNA/pearplot_after.pdf',pearplot_after,height = 6,width = 15)


sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 2000)
sce <- ScaleData(sce, features = rownames(sce))

sce <- RunPCA(sce, features = VariableFeatures(sce))
colnames(sce@meta.data)

library(harmony)
sce = RunHarmony(sce, group.by.vars="Samples")

pca.plot=ElbowPlot(sce,ndims = 50)+theme(text = element_text(family = 'Times',size = 12))
pca.plot
ggsave('results/01.scRNA/PCA_plot.pdf',pca.plot,height = 6,width = 6)

sce <- RunUMAP(sce, dims=1:20, reduction="harmony")#harmony

after_batch=DimPlot(sce,group.by = 'Samples',reduction="umap",label = F,pt.size = 0.2)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(text = element_text(family = 'Times',size = 12),panel.grid = element_blank())
after_batch
ggsave('results/01.scRNA/after_batch.pdf',after_batch,height = 5,width = 6)

library(clustree)
sce <- FindNeighbors(sce, dims = 1:20, reduction="harmony")

sce <- FindClusters(object = sce,resolution = 0.1)
DefaultAssay(sce) <- "RNA"
colnames(sce@meta.data)
length(table(sce@meta.data$seurat_clusters))



seurat_clusters_umap=DimPlot(sce,group.by='seurat_clusters',reduction="umap",label = F,pt.size = 0.2)+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(text = element_text(family = 'Times',size = 12),panel.grid = element_blank())+
  NoLegend()
seurat_clusters_umap=LabelClusters(seurat_clusters_umap,id = 'seurat_clusters',family='Times')
seurat_clusters_umap 
ggsave('results/01.scRNA/seurat_clusters_umap.pdf',seurat_clusters_umap,height = 5,width = 6)



Logfc = 0.5

Minpct = 0.25
Idents(sce)<-'seurat_clusters'
unique(Idents(sce))
sce.markers <- FindAllMarkers(object = sce,logfc.threshold = Logfc, min.pct = Minpct,only.pos = T)
head(sce.markers)
sce.markers["pct.diff"]=sce.markers$pct.1-sce.markers$pct.2
sce.markers <- sce.markers[sce.markers$p_val_adj<0.05,]
table(sce.markers$cluster)
length(unique(sce.markers$gene))#157
head(sce.markers)
table(sce.markers$cluster)
write.csv(sce.markers,'results/01.scRNA/DEGs_seurat_clusters.csv')

Top10 <- sce.markers %>% group_by(cluster) %>% slice_max(n =5, order_by = avg_logFC)
length(Top10$gene)
length(unique(Top10$gene))

DotPlot(object = sce, features = unique(Top10$gene),
        cols=c("snow", "red"),scale = T)+ggtitle("Marker Genes")+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(family = 'Times')) +
  xlab('')+ylab('')+coord_flip()
ggsave('results/01.scRNA/gene_dotplot.pdf',height = 12,width = 10)


# VlnPlot(sce,features = c('CD3E','CD3D','NKG7'),pt.size = 0,group.by = 'seurat_clusters')+NoLegend()
# ggsave('results/01.scRNA/T_vln.pdf',height = 3,width = 10)
# #B/Plasma  1,7,8
# VlnPlot(sce,features = c('CD79A','MZB1','IGKC'),pt.size = 0,group.by = 'seurat_clusters')+NoLegend()
# ggsave('results/01.scRNA/B_vln.pdf',height = 3,width = 10)
# #Fibroblast 2
# VlnPlot(sce,features = c('LUM','COL6A3','FBLN1'),pt.size = 0,group.by = 'seurat_clusters')+NoLegend()
# ggsave('results/01.scRNA/Fibroblast_vln.pdf',height = 3,width = 10)
# #Macrophage 3  
# VlnPlot(sce,features = c('C1QA','C1QB'),pt.size = 0,group.by = 'seurat_clusters')+NoLegend()
# ggsave('results/01.scRNA/Macrophage_vln.pdf',height = 3,width = 7)
# #Neutrophil 4,13  
# VlnPlot(sce,features = c('G0S2','CXCL8'),pt.size = 0,group.by = 'seurat_clusters')+NoLegend()
# ggsave('results/01.scRNA/Neutrophil_vln.pdf',height = 3,width = 7)
# ##Mast cell  5
# VlnPlot(sce,features = c('TPSAB1','TPSB2','CPA3','CTSG'),pt.size = 0,group.by = 'seurat_clusters',ncol = 4)+NoLegend()
# ggsave('results/01.scRNA/Mast_vln.pdf',height = 3,width = 12)
# #Epithelial cell	6
# VlnPlot(sce,features = c('KRT5','KRT14','KRT13'),pt.size = 0,group.by = 'seurat_clusters')+NoLegend()
# ggsave('results/01.scRNA/Epithelial_vln.pdf',height = 3,width = 10)
# #Endothelial cell  9,12
# VlnPlot(sce,features = c('VWF','PLVAP','AQP1'),pt.size = 0,group.by = 'seurat_clusters')+NoLegend()
# ggsave('results/01.scRNA/Endothelial_vln.pdf',height = 3,width = 10)
# #Smooth muscle cell	 10
# VlnPlot(sce,features = c('ACTA2','TAGLN'),pt.size = 0,group.by = 'seurat_clusters')+NoLegend()
# ggsave('results/01.scRNA/SMC_vln.pdf',height = 3,width = 10)


marker <- data.frame(Risktype = 0:13,cell = 0:13)
marker[marker$Risktype %in% c(0,11),2] <- 'T cells'
marker[marker$Risktype %in% c(1,7,8),2] <- 'B/Plasma cells'
marker[marker$Risktype %in% c(2),2] <- 'Fibroblast'
marker[marker$Risktype %in% c(3),2] <- 'Macrophage'
marker[marker$Risktype %in% c(4,13),2] <- 'Neutrophil'
marker[marker$Risktype %in% c(5),2] <- 'Mast cells'
marker[marker$Risktype %in% c(6),2] <- 'Epithelial cells'
marker[marker$Risktype %in% c(9,12),2] <- 'Endothelial cells'
marker[marker$Risktype %in% c(10),2] <- 'Smooth muscle cells'
marker
sce@meta.data$cell_type <- sapply(sce@meta.data$seurat_clusters,function(x){marker[x,2]})
cell_type_umap=DimPlot(sce,group.by='cell_type',reduction="umap",label = F,pt.size = 0.2,cols =pal_futurama()(9))+
  theme_dr(xlength = 0.3, ylength = 0.3,arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")) +
  theme(panel.grid = element_blank(),text = element_text(family = 'Times',size = 12))
cell_type_umap
ggsave('results/01.scRNA/cell_type_umap.pdf',cell_type_umap,height = 5,width = 6)

table(sce@meta.data$cell_type,sce@meta.data$Tissue)

feat_gene<-c('C1QA','C1QB',
             'CD79A','MZB1','IGKC',
             'KRT5','KRT14','KRT13',
             'G0S2','CXCL8',
             'CD3E','CD3D','NKG7',
             'LUM','DCN','COL6A3',
             'ACTA2','TAGLN',
             'VWF','PLVAP','AQP1',
             'TPSAB1','TPSB2','CPA3','CTSG')


Idents(sce)='cell_type'
marker.dotplot=DotPlot(object = sce, features = feat_gene,
                            cols=c("snow", "red"),scale = T)+
  ggtitle("Marker Genes")+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(family = 'Times'),axis.text.x = element_text(angle = 30,hjust = 1)) +
  xlab('')+ylab('')+coord_flip()
marker.dotplot
ggsave('results/01.scRNA/marker.dotplot.pdf',marker.dotplot,height = 6,width = 6)


chisq.test(table(sce$cell_type,sce$Tissue))
cell_freq1=data.frame(t(prop.table(table(sce$cell_type,sce$Samples),margin=2)),tissue)
cell_freq1
colnames(cell_freq1)<-c('Samples','cell_type','Freq','Tissue')
cell_prop_fig1=ggplot(cell_freq1,aes(x=Samples,y=Freq,fill=cell_type))+
  scale_fill_manual(values = pal_futurama()(9))+
  facet_grid(~Tissue,scales = 'free',space='free')+
  geom_bar(position = "fill",stat="identity")+
  xlab('')+ylab('Proportion')+
  theme(text = element_text(family = 'Times',size=12),
        axis.text.x = element_text(angle = 30,hjust = 1),
        legend.text =element_text(family = 'Times'),
        legend.title = element_text(family = 'Times'))
cell_prop_fig1

cell_freq2=melt(prop.table(table(sce$cell_type,sce$Tissue),margin=2))
cell_freq2
colnames(cell_freq2)<-c('cell_type','Tissue','Freq')
cell_prop_fig2=ggplot(cell_freq2,aes(x=reorder(cell_type,Freq),y=Freq,fill=Tissue))+
  scale_fill_manual(values = pal_d3()(10))+
  geom_bar(stat="identity", position="dodge")+
  xlab('')+ylab('Proportion')+coord_flip()+
  theme(text = element_text(family = 'Times',size=12),
        legend.position = 'top',legend.text =element_text(family = 'Times'),
        legend.title = element_text(family = 'Times'))
cell_prop_fig2

cell_prop_fig=mg_merge_plot(cell_prop_fig1,cell_prop_fig2)
cell_prop_fig
ggsave('results/01.scRNA/cell_prop_fig.pdf',cell_prop_fig,height = 5,width = 12)

saveRDS(sce,file = 'results/01.scRNA/sce.rds')


Logfc = 0.5

Minpct = 0.25
Idents(sce)<-'cell_type'
unique(Idents(sce))
sce.markers <- FindAllMarkers(object = sce,logfc.threshold = Logfc, min.pct = Minpct,only.pos = T)
head(sce.markers)
sce.markers["pct.diff"]=sce.markers$pct.1-sce.markers$pct.2
sce.markers <- sce.markers[sce.markers$p_val_adj<0.05,]
table(sce.markers$cluster)
length(unique(sce.markers$gene))#157
head(sce.markers)
table(sce.markers$cluster)
write.csv(sce.markers,'results/01.scRNA/DEGs_cell_type.csv')

Top10 <- sce.markers %>% group_by(cluster) %>% slice_max(n =5, order_by = avg_logFC)
length(Top10$gene)
length(unique(Top10$gene))

DotPlot(object = sce, features = unique(Top10$gene),
                            cols=c("snow", "red"),scale = T)+
  ggtitle("Marker Genes")+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(family = 'Times')) +
  xlab('')+ylab('')+coord_flip()+RotatedAxis()
ggsave('results/01.scRNA/diff_marker_dotplot.pdf',height = 10,width = 8)

fig1=mg_merge_plot(cell_type_umap,marker.dotplot,cell_prop_fig1,cell_prop_fig2,nrow=2,ncol=2,labels = LETTERS[1:4])
ggsave('results/01.scRNA/Fig1.pdf',fig1,height = 10,width = 10)



dir.create('results/01.scRNA/Epithelial')
table(sce$cell_type)
Epithelial=subset(sce,cell_type %in% 'Epithelial cells')
Epithelial <- NormalizeData(Epithelial)
Epithelial <- FindVariableFeatures(Epithelial, selection.method = "vst", nfeatures = 2000)
Epithelial <- ScaleData(Epithelial, features = rownames(Epithelial))

Epithelial <- RunPCA(Epithelial, features = VariableFeatures(Epithelial))
colnames(Epithelial@meta.data)

Epithelial = RunHarmony(Epithelial, group.by.vars="Samples")

ElbowPlot(Epithelial,ndims = 50)+theme(text = element_text(family = 'Times',size = 12))

Epithelial <- RunUMAP(Epithelial, dims=1:10, reduction="harmony")


mallignant <- read.delim("results/01.scRNA/Epithelial/Epithelial_copykat_prediction.txt", row.names = 1)
head(mallignant)
mallignant$copykat.pred=ifelse(mallignant$copykat.pred=='aneuploid','aneuploid','diploid')

Epithelial <- AddMetaData(Epithelial, metadata = mallignant)

p1=DimPlot(Epithelial, group.by = "Tissue")+ 
  scale_color_manual(values = c("grey", "blue"))+theme(text = element_text(family = 'Times',size = 12))

p2=DimPlot(Epithelial, group.by = "copykat.pred") + scale_color_manual(values = c("red", "grey"))+
  theme(text = element_text(family = 'Times',size = 12))

p1|p2
ggsave('results/01.scRNA/Epithelial/copykat_umap.pdf',height = 5,width = 12)

library(AUCell)
cells_rankings <- AUCell_buildRankings(Epithelial@assays$RNA@data)
library(clusterProfiler)
DDR.gmt=read.gmt('origin_datas/WP_DNA_DAMAGE_RESPONSE.v2023.1.Hs.gmt')
geneSets <- list(DDR=DDR.gmt$gene)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
AUCs <- as.numeric(getAUC(cells_AUC)['DDR', ])
Epithelial$AUCs  <- AUCs

sc_df=data.frame(Epithelial@meta.data, Epithelial@reductions$umap@cell.embeddings)
head(sc_df)


library(stats)
library(ggpubr)
score_stat <- sc_df  %>%
  wilcox_test(AUCs ~ copykat.pred, paired = FALSE, p.adjust.method = "fdr")%>%
  add_significance() %>% add_xy_position(x = "copykat.pred")
score_stat

p3=ggplot(sc_df,aes(x=copykat.pred,y = AUCs))+
  geom_boxplot(aes(fill = copykat.pred), varwidth = T, alpha = 0.9, outlier.shape = NA)+
  stat_pvalue_manual(score_stat, label = "p.signif", tip.length = 0.02)+
  ylab('DDR score')+ theme(text = element_text(family = 'Times',size = 12))
p3
ggsave('results/01.scRNA/Epithelial/DDR_compare.pdf',height = 5,width = 5)
saveRDS(Epithelial,file = 'results/01.scRNA/Epithelial/Epithelial.rds')


Idents(Epithelial)<-'copykat.pred'
unique(Idents(Epithelial))
Epithelial.markers <- FindAllMarkers(object = Epithelial,logfc.threshold = 0.25, min.pct = 0.25,only.pos = T)
head(Epithelial.markers)
Epithelial.markers["pct.diff"]=Epithelial.markers$pct.1-Epithelial.markers$pct.2
Epithelial.markers <- Epithelial.markers[Epithelial.markers$p_val_adj<0.05,]
table(Epithelial.markers$cluster)
length(unique(Epithelial.markers$gene))#1239
head(Epithelial.markers)
table(Epithelial.markers$cluster)
write.csv(Epithelial.markers,'results/01.scRNA/Epithelial/Epithelial_DEGs.csv')

ep.top.marker=Epithelial.markers %>% group_by(cluster) %>% slice_max(n =10, order_by = avg_logFC)
Idents(Epithelial)='copykat.pred'
p4=DotPlot(object = Epithelial, features = unique(ep.top.marker$gene),
        cols=c("snow", "red"),scale = T)+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(family = 'Times')) +
  xlab('')+ylab('')+coord_flip()
p4
ggsave('results/01.scRNA/Epithelial/Epithelial_DEGs_dotplot.pdf',height = 5,width = 5)


fig2=mg_merge_plot(p1,p2,p3,p4,nrow=2,ncol=2,labels = LETTERS[1:4])
ggsave('results/01.scRNA/Fig2.pdf',fig2,height = 10,width = 10)


rm(list=ls())

Epithelial.markers=read.csv('results/01.scRNA/Epithelial/Epithelial_DEGs.csv')

tcga.pancancer.cli=read.xlsx('origin_datas/TCGA/TCGA_pancancer_cli_PMID_29625055.xlsx')
head(tcga.pancancer.cli)
table(tcga.pancancer.cli$type)
tcga.cli=tcga.pancancer.cli[which(tcga.pancancer.cli$type=='ESCA'),]
head(tcga.cli)
tcga.cli=data.frame(Samples=paste0(tcga.cli$bcr_patient_barcode,'-01'),
                    Age=tcga.cli$age_at_initial_pathologic_diagnosis,
                    Gender=tcga.cli$gender,
                    AJCC_stage=tcga.cli$ajcc_pathologic_tumor_stage,
                    Grade=tcga.cli$histological_grade,
                    tcga.cli[,c('OS','OS.time','DSS','DSS.time','DFI','DFI.time','PFI','PFI.time')])
rownames(tcga.cli)=tcga.cli$Samples

head(tcga.cli)
tcga.cli$OS.time
tcga.cli=tcga.cli %>% drop_na(OS.time)
tcga.cli=tcga.cli[tcga.cli$OS.time>0,]
dim(tcga.cli)
head(tcga.cli)

table(tcga.cli$Gender)
table(tcga.cli$AJCC_stage)
tcga.cli$AJCC_stage[tcga.cli$AJCC_stage=='[Discrepancy]'|
                      tcga.cli$AJCC_stage=='[Not Available]'|
                      tcga.cli$AJCC_stage=='[Unknown]']=NA
tcga.cli$AJCC_stage=gsub('[ABC]','',tcga.cli$AJCC_stage)
tcga.cli$AJCC_stage=gsub('Stage ','',tcga.cli$AJCC_stage)

table(tcga.cli$Grade)
tcga.cli$Grade[tcga.cli$Grade=='GX']=NA



tcga.data=read.delim('origin_datas/TCGA/TCGA_ESCA_TPM.txt',check.names = F,row.names = 1)
tcga.data[1:4,1:4]
table(substr(colnames(tcga.data),14,15))
sample.T=colnames(tcga.data)[substr(colnames(tcga.data),14,15)=='01']
sample.N=colnames(tcga.data)[substr(colnames(tcga.data),14,15)=='11']
tcga.type=data.frame(Samples=c(sample.T,sample.N),type=rep(c('Tumor','Normal'),c(length(sample.T),length(sample.N))))
rownames(tcga.type)=tcga.type$Samples

range(tcga.data)
tcga.data=log2(tcga.data[,tcga.type$Samples]+1)
range(tcga.data)

tcga.sample=intersect(colnames(tcga.data),rownames(tcga.cli))
length(tcga.sample)

tcga.exp=tcga.data[,tcga.sample]
tcga.cli=tcga.cli[tcga.sample,]
tcga.cli=crbind2DataFrame(tcga.cli)
dim(tcga.cli)
fivenum(tcga.cli$Age)
tcga.cli$Age1=ifelse(tcga.cli$Age>60,'>60','<=60')
range(tcga.exp)
dim(tcga.exp)
tcga.exp[1:5,1:5]





dir.create('results/02.DDRscore')
DDR.gmt=read.gmt('origin_datas/WP_DNA_DAMAGE_RESPONSE.v2023.1.Hs.gmt')
DDR.score=t(ssGSEAScore_by_genes(gene.exp = tcga.data,genes = DDR.gmt$gene))
write.csv(DDR.score,'results/02.DDRscore/TCGA_DDRscore.csv')


load('results/tcga.hall.ssGSEA.RData')
tcga.hall.ssGSEA[1:5,1:5]
rownames(tcga.hall.ssGSEA)=gsub('HALLMARK_','',rownames(tcga.hall.ssGSEA))
rownames(tcga.hall.ssGSEA)=gsub('_',' ',rownames(tcga.hall.ssGSEA))
library(ggcorrplot)
library(psych)
hallmark_RS_cor <- corr.test(x =as.numeric(DDR.score[tcga.cli$Samples,1]),
                             y = t(tcga.hall.ssGSEA[,tcga.cli$Samples]),
                             method = "spearman",adjust = "BH",ci = F)

hallmark_RS_cor_res=data.frame(pathway=colnames(t(tcga.hall.ssGSEA[,tcga.cli$Samples])))
hallmark_RS_cor_res$cor<-as.numeric(hallmark_RS_cor$r)
hallmark_RS_cor_res$p.adj<-as.numeric(hallmark_RS_cor$p.adj)
head(hallmark_RS_cor_res)
write.csv(hallmark_RS_cor_res,'results/02.DDRscore/TCGA_hallmark_cor_DDR.csv',row.names = F)
hallmark_RS_cor_res=hallmark_RS_cor_res[order(hallmark_RS_cor_res$cor),]


library(rcartocolor)
ggplot(data=hallmark_RS_cor_res,aes(x=cor,y=reorder(pathway,cor), color = -log10(p.adj))) +
  geom_point(aes(size=abs(cor)),show.legend = F) +
  scale_color_continuous(type = "viridis")+
  geom_segment(aes(yend=pathway,xend=0),size=.5) +
  labs(x='spearman Correlation',y='')+theme_bw()+
  theme(text = element_text(family = 'Times'))

ggsave('results/02.DDRscore/DDRscore_cor_hallmar.pdf',height = 10,width = 7)


df=data.frame(DDR=DDR.score[tcga.type$Samples,1],Type=tcga.type$type)
head(df)
library(stats)
library(ggpubr)
score_stat <- df  %>%
  wilcox_test(DDR ~ Type, paired = FALSE, p.adjust.method = "fdr")%>%
  add_significance() %>% add_xy_position(x = "Type")
score_stat

p1=ggplot(df,aes(x=Type,y = DDR))+
  geom_boxplot(aes(fill = Type), varwidth = T, alpha = 0.9, outlier.shape = NA)+
  stat_pvalue_manual(score_stat, label = "p.signif", tip.length = 0.01)+
  ylab('DDR score')+ theme_classic()+theme(text = element_text(family = 'Times',size = 15),legend.position = 'none')
 



colnames(tcga.immu.ssgsea)
tcga.ifnγ=immu_Th1_IFNγScore(tcga.exp[,tcga.cli$Samples])

df=cbind.data.frame(DDRscore=as.numeric(DDR.score[tcga.cli$Samples,1]),
                    IFNγ=as.numeric(tcga.ifnγ),
                    tcga.immu.ssgsea[tcga.cli$Samples,c(3,5,7,8,17,18,25,26)]
)
cor_res <- Hmisc::rcorr(as.matrix(df),type = 'spearman')
cor_res$P[is.na(cor_res$P)] <- 0
cor_res$r


library(ggpubr)
p2=ggscatter(df,x = "DDRscore", y = "Natural killer cell",
          add = "reg.line", conf.int = TRUE,
          color = "black", shape = 19, size = 2, # Points color, shape and size
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          cor.coef = TRUE, cor.method = "spearman",
          xlab = 'DDR score',ylab = 'Natural killer cell')

p3=ggscatter(df,x = "DDRscore", y = "Activated CD8 T cell",
          add = "reg.line", conf.int = TRUE,
          color = "black", shape = 19, size = 2, # Points color, shape and size
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          cor.coef = TRUE, cor.method = "spearman",
          xlab = 'DDR score',ylab = 'Activated CD8 T cell')

p4=ggscatter(df,x = "DDRscore", y = "IFNγ",
          add = "reg.line", conf.int = TRUE,
          color = "black", shape = 19, size = 2, # Points color, shape and size
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          cor.coef = TRUE, cor.method = "spearman",
          xlab = 'DDR score',ylab = 'IFN-γ')

mg_merge_plot(p1,p2,p3,p4,ncol=2,nrow=2,labels = LETTERS[1:4])
ggsave('results/02.DDRscore/DDRscore_TvsN.pdf',height = 10,width =10)

###########03.WGCNA##########################
dir.create('results/03.WGCNA')

# Risktype gene
library(WGCNA)
allowWGCNAThreads(nThreads = 36)
enableWGCNAThreads(nThreads = 36)


genecode=read.delim('origin_datas/GeneTag.genecode.v32.txt',sep='\t',header = T)
table(genecode$TYPE)
mrna_genecode=genecode[which(genecode$TYPE=='protein_coding'),]
head(mrna_genecode)

tcga.exp=tcga.exp[intersect(mrna_genecode$SYMBOL,rownames(tcga.exp)),]
tpm_T2=tcga.exp[which(apply(tcga.exp,1,sd)>0.5),]
dim(tpm_T2)
tpm_T2=t(tpm_T2)
dim(tpm_T2)
range(tpm_T2)

pdf('results/03.WGCNA/1.pdf',width = 8,height = 8)
tpm_T2.power=mg_wgcna_get_power(tpm_T2)
dev.off()


tpm_T2.power$cutPower
tpm_T2.module=mg_WGCNA_getModule(tpm_T2,
                                 power = tpm_T2.power$cutPower,
                                 deepSplit=2,
                                 mergeCutHeight=0.25,
                                 minModuleSize=50)

table(tpm_T2.module$Modules[,2])
length(table(tpm_T2.module$Modules[,2]))

pdf('results/03.WGCNA/2.pdf',height = 6,width = 10)
plotDendroAndColors(tpm_T2.module$Tree, tpm_T2.module$Modules,
                    c("Dynamic Module",'Merged Module'),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
# 
#writeMatrix(tpm_T2.module$Modules,outpath = 'results/03.WGCNA/tcga.wgcna.module.genes.txt')
pdf('results/03.WGCNA/3.pdf',height = 6,width = 6)
mg_barplot_point(labels = names(table(tpm_T2.module$Modules[,2]))
                 ,values = as.numeric(table(tpm_T2.module$Modules[,2]))
                 ,point_sizes = 2
                 ,point_cols = names(table(tpm_T2.module$Modules[,2]))
                 ,xlab = 'Number of Genes',legend.pos = NULL)
dev.off()


# Calculate eigengenes
MEs = tpm_T2.module$MEs
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Risktype module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf('results/03.WGCNA/4.pdf',height = 6,width = 12,onefile = T)
plot(METree, main = "Risktypeing of module eigengenes",xlab = "", sub = "")
dev.off()




tcga_cli_use=cbind(DDR.score[tcga.cli$Samples,],tcga.cli)
rownames(tcga_cli_use)=tcga_cli_use$Samples
colnames(tcga_cli_use)[1]='DDR'
colnames(tcga_cli_use)
colnames(tcga_cli_use)[c(1,5,6)]


tcga_cli_use.part=tcga_cli_use[tcga.cli$Samples,c(1,5,6)]
str(tcga_cli_use.part)

tcga_cli_use.part=sapply(tcga_cli_use.part, function(x)as.numeric(as.factor(x)))
spms=tcga_cli_use.part

MEs_col<-tpm_T2.module$MEs
dim(MEs_col)

modTraitCor = cor(MEs_col[,rownames(MEDiss)[METree$order]]
                  , spms
                  ,use = 'pairwise.complete.obs')
modTraitP = corPvalueStudent(modTraitCor, dim(spms)[1])

textMatrix = paste(signif(modTraitCor, 2), " (", format(modTraitP,scientific =TRUE,digits = 3), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
dim(textMatrix)

rownames(modTraitCor)=gsub("ME","",rownames(modTraitCor))
rownames(textMatrix)=gsub("ME","",rownames(textMatrix))
colnames(modTraitCor)

pdf('results/03.WGCNA/5.pdf',width =10,height = 10)
labeledHeatmap(Matrix = data.frame(modTraitCor),
               xLabels = colnames(modTraitCor),
               yLabels = rownames(modTraitCor),
               cex.lab = 1,
               ySymbols = colnames(t(modTraitCor)), colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = data.frame(textMatrix),
               setStdMargins = FALSE,xLabelsAngle = 0,
               cex.text = 0.8, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


geneModuleMembership <- signedKME(tpm_T2
                                  , data.frame(tpm_T2.module$MEs)
                                  , outputColumnName = "")
head(geneModuleMembership)
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership)
                                           , nrow(tpm_T2.module$MEs)))

geneTraitSignificance <- as.data.frame(cor(tpm_T2
                                           , spms
                                           , use = 'pairwise.complete.obs'))

GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance)
                                           , nrow(spms)))

modNames<-colnames(geneModuleMembership)
modNames

module = "brown"
column = match(module, modNames)
column
moduleGenes <- (tpm_T2.module$Modules[,'mergedColors']==module)
brown.gene=names(which(moduleGenes))
length(brown.gene)
#1092

pdf('results/03.WGCNA/6.pdf',height = 5,width = 6)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 'DDR']),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for DDR",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2
                   , col = module,lwd=2)
dev.off()

wgcna.enrich=mg_clusterProfiler(brown.gene)
enrichplot::dotplot(wgcna.enrich$KEGG)+ggtitle('KEGG')+theme(text = element_text(family = 'Times'))
ggsave('results/03.WGCNA/wgcna.enrich.pdf',height = 5,width = 6)
write.csv(wgcna.enrich$KEGG@result,'results/03.WGCNA/wgcna_KEGG_res.csv',row.names = F)


dir.create('results/04.model')
pre.gene=intersect(Epithelial.markers$gene,brown.gene)
length(pre.gene)
#125
library(eulerr)
v=list(Epithelial.markers$gene,brown.gene)
names(v)=c('DEGs','DDR related-genes')
pdf('results/04.model/venn.pdf',height = 4,width = 7)
plot(venn(v),
     labels = list(col = "gray20", font = 2), 
     edges = list(col="gray60", lex=1),
     fills = list(fill = c("#297CA0", "#E9EA77"), alpha = 0.6),
     quantities = list(cex=.8, col='gray20'))
dev.off()


gene.cox=cox_batch(dat = tcga.exp[pre.gene,tcga.cli$Samples],time = tcga.cli$OS.time,event = tcga.cli$OS)
gene.cox=na.omit(gene.cox)
table(gene.cox$p.value<0.05)
gene.cox=gene.cox[gene.cox$p.value<0.05,]
gene.cox

tcga_model_data <- cbind(tcga.cli[, c("OS.time", "OS")],
                         t(tcga.exp[rownames(gene.cox), tcga.cli$Samples]))
colnames(tcga_model_data) <- gsub('-', '_', colnames(tcga_model_data))



fmla <- as.formula(paste0("Surv(OS.time, OS) ~"
                          ,paste0(rownames(gene.cox),collapse = '+')))
cox <- coxph(fmla, data =as.data.frame(tcga_model_data))
cox=step(cox)
lan <- coef(cox)
lan
paste0(round(lan, 3), '*', names(lan),collapse = '+')
ggforest(cox, data = tcga_model_data, 
                          main = "Hazardratio", fontsize =1.0, 
                          noDigits = 2)+theme(text = element_text(family = 'Times'))
ggsave('results/04.model/module.coxforest.pdf',height = 4,width = 7)

risktype.col=pal_jco()(10)[9:10]
risk.tcga=as.numeric(lan%*%as.matrix(t(tcga_model_data[tcga.cli$Samples,names(lan)])))
tcga.risktype.cli=data.frame(tcga.cli,Riskscore=risk.tcga)
tcga.point <- surv_cutpoint(tcga.risktype.cli, time = "OS.time", event = "OS",
                            variables = 'Riskscore')
tcga.point.cutoff <- as.numeric(summary(tcga.point)[1])
tcga.point.cutoff

tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>tcga.point.cutoff,'High','Low')
tcga.roc.OS=ggplotTimeROC(tcga.risktype.cli$OS.time,
                          tcga.risktype.cli$OS,
                          tcga.risktype.cli$Riskscore,mks = c(1:4))
tcga.roc.OS
ggsave('results/04.model/tcga.roc.OS.pdf',height = 5,width = 5)

tcga.km.OS=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ Risktype,
                                  data = tcga.risktype.cli),
                      data=tcga.risktype.cli,
                      conf.int = T,pval = T,risk.table = T, 
                      fun = "pct",size = 1,surv.median.line = 'hv',
                      title='Overall Survival(OS)',legend.title='Risktype',
                      legend.labs = c('High','Low'),
                      linetype = c("solid", "dashed","strata")[1],
                      palette = risktype.col,#ylab='Overall Survival(OS)',
                      legend.position='top',
                      ggtheme = theme_bw(base_size = 12))
tcga.km.OS=mg_merge_plot(tcga.km.OS$plot,tcga.km.OS$table,nrow=2,heights = c(3,1),align = 'v')
tcga.km.OS


tcga.km.PFS=ggsurvplot(fit=survfit(Surv(PFI.time/365, PFI) ~ Risktype,
                                   data = tcga.risktype.cli),
                       data=tcga.risktype.cli,
                       conf.int = T,pval = T,risk.table = T, 
                       fun = "pct",size = 1,surv.median.line = 'hv',
                       title='Progression Free Interval(PFI)',legend.title='Risktype',
                       legend.labs = c('High','Low'),
                       linetype = c("solid", "dashed","strata")[1],
                       palette = risktype.col,#ylab='Progression Free Interval(PFI)',
                       legend.position='top',
                       ggtheme = theme_bw(base_size = 12))
tcga.km.PFS=mg_merge_plot(tcga.km.PFS$plot,tcga.km.PFS$table,nrow=2,heights = c(3,1),align = 'v')
tcga.km.PFS



tcga.km.DSS=ggsurvplot(fit=survfit(Surv(DSS.time/365, DSS) ~ Risktype,
                                   data = tcga.risktype.cli),
                       data=tcga.risktype.cli,
                       conf.int = T,pval = T,risk.table = T, 
                       fun = "pct",size = 1,surv.median.line = 'hv',
                       title='Disease Specific Survival(DSS)',legend.title='Risktype',
                       legend.labs = c('High','Low'),
                       linetype = c("solid", "dashed","strata")[1],
                       palette = risktype.col,#ylab='Disease Specific Survival(DSS)',
                       legend.position='top',
                       ggtheme = theme_bw(base_size = 12))
tcga.km.DSS=mg_merge_plot(tcga.km.DSS$plot,tcga.km.DSS$table,nrow=2,heights = c(3,1),align = 'v')
tcga.km.DSS

tcga.km=mg_merge_plot(tcga.km.OS,tcga.km.PFS,tcga.km.DSS,ncol=3,labels = LETTERS[3:5])
ggsave('results/04.model/tcga.km.pdf',tcga.km,height = 5,width = 15)

library(ggbiplot)
risktype.pca <- prcomp(t(tcga.exp[names(lan),tcga.risktype.cli$Samples]), scale=T)
pca_plot <- ggbiplot(risktype.pca, scale=1, groups = tcga.risktype.cli$Risktype,
                       ellipse = TRUE,ellipse.prob=0.5, circle = F,var.axes=F) +
  scale_color_manual(values =risktype.col) + 
  # xlim(-5, 5) + ylim(-5,5) +
  theme_light() +
  theme(legend.direction = 'horizontal', legend.position = 'top') +
  xlab('PCA1') + ylab('PCA2')+theme(text = element_text(family = 'Times'))
pca_plot
ggsave('results/04.model/PCA.pdf',pca_plot,height = 5,width = 5)



dir.create('results/05.clinical')
tcga_cox_datas=tcga.risktype.cli
colnames(tcga_cox_datas)
table(tcga_cox_datas$Age1)


table(tcga_cox_datas$AJCC_stage)
tcga_cox_datas$AJCC_stage[tcga_cox_datas$AJCC_stage=='I'|tcga_cox_datas$AJCC_stage=='II']<-'I+II'
tcga_cox_datas$AJCC_stage[tcga_cox_datas$AJCC_stage=='III'|tcga_cox_datas$AJCC_stage=='IV']<-'III+IV'


table(tcga_cox_datas$Grade)
tcga_cox_datas$Grade[tcga_cox_datas$Grade=='G1'|tcga_cox_datas$Grade=='G2']<-'G1+G2'
tcga_cox_datas$Grade[tcga_cox_datas$Grade=='G3']<-'G3'


#Age
tcga_cox_datas=crbind2DataFrame(tcga_cox_datas)
Age_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age,
                             data=tcga_cox_datas))
Age_sig_cox_dat <- data.frame(Names=rownames(Age_sig_cox[[8]]),
                              HR = round(Age_sig_cox[[7]][,2],3),
                              lower.95 = round(Age_sig_cox[[8]][,3],3),
                              upper.95 = round(Age_sig_cox[[8]][,4],3),
                              p.value=round(Age_sig_cox[[7]][,5],3))
Age_sig_cox_dat

#Gender
Gender_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Gender,
                                data=tcga_cox_datas))
Gender_sig_cox_dat <- data.frame(Names=rownames(Gender_sig_cox[[8]]),
                                 HR = round(Gender_sig_cox[[7]][,2],3),
                                 lower.95 = round(Gender_sig_cox[[8]][,3],3),
                                 upper.95 = round(Gender_sig_cox[[8]][,4],3),
                                 p.value=round(Gender_sig_cox[[7]][,5],3))
Gender_sig_cox_dat
#AJCC_stage
Stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~AJCC_stage,
                               data=tcga_cox_datas))
Stage_sig_cox_dat <- data.frame(Names=rownames(Stage_sig_cox[[8]]),
                                HR = round(Stage_sig_cox[[7]][,2],3),
                                lower.95 = round(Stage_sig_cox[[8]][,3],3),
                                upper.95 = round(Stage_sig_cox[[8]][,4],3),
                                p.value=round(Stage_sig_cox[[7]][,5],3))
Stage_sig_cox_dat

#Grade
Grade_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Grade,
                               data=tcga_cox_datas))
Grade_sig_cox_dat <- data.frame(Names=rownames(Grade_sig_cox[[8]]),
                                HR = round(Grade_sig_cox[[7]][,2],3),
                                lower.95 = round(Grade_sig_cox[[8]][,3],3),
                                upper.95 = round(Grade_sig_cox[[8]][,4],3),
                                p.value=round(Grade_sig_cox[[7]][,5],3))
Grade_sig_cox_dat

#riskscore
riskscore_sig_cox<-summary(coxph(formula=Surv(OS.time, OS)~Riskscore,
                                 data=tcga_cox_datas))
riskscore_sig_cox_dat <- data.frame(Names=rownames(riskscore_sig_cox[[8]]),
                                    HR = round(riskscore_sig_cox[[7]][,2],3),
                                    lower.95 = round(riskscore_sig_cox[[8]][,3],3),
                                    upper.95 = round(riskscore_sig_cox[[8]][,4],3),
                                    p.value=round(riskscore_sig_cox[[7]][,5],3))
riskscore_sig_cox_dat

sig_cox_dat <- rbind(Age_sig_cox_dat,
                     Gender_sig_cox_dat,
                     Stage_sig_cox_dat,
                     Grade_sig_cox_dat,
                     riskscore_sig_cox_dat)
data.sig <- data.frame(Features=sig_cox_dat$Names,
                       p.value=sig_cox_dat$p.value,
                       sig_cox_dat$HR,
                       sig_cox_dat$lower.95,
                       sig_cox_dat$upper.95)
data.sig <- crbind2DataFrame(data.sig)
rownames(data.sig) <- c("Age","Gender","AJCC stage","Grade","RiskScore")
data.sig$Features=rownames(data.sig) 
data.sig
data.sig$p.value=ifelse(data.sig$p.value<0.001,'<0.001',data.sig$p.value)
pdf('results/04.model/Univariate.pdf',height = 3,width = 6,onefile = F)
mg_forestplot_v2(data.sig,xlog = T,colgap =5,lineheight = 10,zero = 1,
                 boxsize = 0.3,lwd.zero=1,lwd.ci=1,lwd.xaxis=1
                 ,box_col='blue',summary_col="black",lines_col='black',zero_col='grey'
                 ,xlab='Hazard Ratio',lty.ci = 6,graph.pos = 4)
dev.off()

#T.stage +  N.stage +M.stage
muti_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~AJCC_stage+Riskscore, 
                              data=tcga_cox_datas))
muti_cox_dat <- data.frame(Names=rownames(muti_sig_cox[[8]]),
                           HR = round(muti_sig_cox[[7]][,2],3),
                           lower.95 = round(muti_sig_cox[[8]][,3],3),
                           upper.95 = round(muti_sig_cox[[8]][,4],3),
                           p.value=round(muti_sig_cox[[7]][,5],3))
data.muti <- data.frame(Features=muti_cox_dat$Names,
                        p.value=muti_cox_dat$p.value,
                        muti_cox_dat$HR,
                        muti_cox_dat$lower.95,
                        muti_cox_dat$upper.95)
data.muti <- crbind2DataFrame(data.muti)
data.muti
rownames(data.muti) <- c("AAJCC stage","RiskScore")
data.muti$Features=rownames(data.muti)
data.muti
data.muti$p.value=ifelse(data.muti$p.value<0.001,'<0.001',data.muti$p.value)
pdf('results/04.model/Multivariate.pdf',height = 2.5,width = 6,onefile = F)
mg_forestplot_v2(data.muti,xlog = T,colgap =5,lineheight = 10,zero = 1,
                 boxsize = 0.3,lwd.zero=1,lwd.ci=1,lwd.xaxis=1
                 ,box_col='blue',summary_col="black",lines_col='black',zero_col='grey'
                 ,xlab='Hazard Ratio',lty.ci = 6,graph.pos = 4)
dev.off()

head(tcga.risktype.cli)





head(tcga_cox_datas)
cli.km1=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                data = tcga_cox_datas[which(tcga_cox_datas$AJCC_stage=='I+II'),]),
                   data=tcga_cox_datas[which(tcga_cox_datas$AJCC_stage=='I+II'),],
                   conf.int = F,pval = T,fun = "pct",risk.table = F, size = 0.7,
                   title='AJCC stage I+II',#ggtheme=custom_theme(),
                   linetype = c("solid", "dashed","strata")[1],
                   legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                   legend.title = "Risktype",palette = risktype.col,
                   legend.labs = c("High","Low"))
cli.km1

cli.km2=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                data = tcga_cox_datas[which(tcga_cox_datas$AJCC_stage=='III+IV'),]),
                   data=tcga_cox_datas[which(tcga_cox_datas$AJCC_stage=='III+IV'),],
                   conf.int = F,pval = T,fun = "pct",risk.table = F, size = 0.7,
                   title='AJCC stage III+IV',#ggtheme=custom_theme(),
                   linetype = c("solid", "dashed","strata")[1],
                   legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                   legend.title = "Risktype",palette = risktype.col,
                   legend.labs = c("High","Low"))
cli.km2



table(tcga_cox_datas$Grade)
cli.km3=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                data = tcga_cox_datas[which(tcga_cox_datas$Grade=='G1+G2'),]),
                   data=tcga_cox_datas[which(tcga_cox_datas$Grade=='G1+G2'),],
                   conf.int = F,pval = T,fun = "pct",risk.table = F, size = 0.7,
                   title='Grade G1+G2',#ggtheme=custom_theme(),
                   linetype = c("solid", "dashed","strata")[1],
                   legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                   legend.title = "Risktype",palette = risktype.col,
                   legend.labs = c("High","Low"))
cli.km3

cli.km4=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                data = tcga_cox_datas[which(tcga_cox_datas$Grade=='G3'),]),
                   data=tcga_cox_datas[which(tcga_cox_datas$Grade=='G3'),],
                   conf.int = F,pval = T,fun = "pct",risk.table = F, size = 0.7,
                   title='Grade G3',#ggtheme=custom_theme(),
                   linetype = c("solid", "dashed","strata")[1],
                   legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                   legend.title = "Risktype",palette = risktype.col,
                   legend.labs = c("High","Low"))
cli.km4



table(tcga_cox_datas$Gender)
cli.km5=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                data = tcga_cox_datas[which(tcga_cox_datas$Gender=='MALE'),]),
                   data=tcga_cox_datas[which(tcga_cox_datas$Gender=='MALE'),],
                   conf.int = F,pval = T,fun = "pct",risk.table = F, size = 0.7,
                   title='MALE',#ggtheme=custom_theme(),
                   linetype = c("solid", "dashed","strata")[1],
                   legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                   legend.title = "Risktype",palette = risktype.col,
                   legend.labs = c("High","Low"))
cli.km5

cli.km6=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                data = tcga_cox_datas[which(tcga_cox_datas$Gender=='FEMALE'),]),
                   data=tcga_cox_datas[which(tcga_cox_datas$Gender=='FEMALE'),],
                   conf.int = F,pval = T,fun = "pct",risk.table = F, size = 0.7,
                   title='FEMALE',#ggtheme=custom_theme(),
                   linetype = c("solid", "dashed","strata")[1],
                   legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                   legend.title = "Risktype",palette = risktype.col,
                   legend.labs = c("High","Low"))
cli.km6


table(tcga_cox_datas$Age1)
cli.km7=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                data = tcga_cox_datas[which(tcga_cox_datas$Age1=='<=60'),]),
                   data=tcga_cox_datas[which(tcga_cox_datas$Age1=='<=60'),],
                   conf.int = F,pval = T,fun = "pct",risk.table = F, size = 0.7,
                   title='Age<=60',#ggtheme=custom_theme(),
                   linetype = c("solid", "dashed","strata")[1],
                   legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                   legend.title = "Risktype",palette = risktype.col,
                   legend.labs = c("High","Low"))
cli.km7

cli.km8=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
                                data = tcga_cox_datas[which(tcga_cox_datas$Age1=='>60'),]),
                   data=tcga_cox_datas[which(tcga_cox_datas$Age1=='>60'),],
                   conf.int = F,pval = T,fun = "pct",risk.table = F, size = 0.7,
                   title='Age>60',#ggtheme=custom_theme(),
                   linetype = c("solid", "dashed","strata")[1],
                   legend = c('top', 'bottom', 'left', 'right', 'none')[1],
                   legend.title = "Risktype",palette = risktype.col,
                   legend.labs = c("High","Low"))
cli.km8

clinical.km=mg_merge_plot(cli.km8$plot,cli.km7$plot,cli.km6$plot,cli.km5$plot,cli.km3$plot,cli.km4$plot,cli.km1$plot,cli.km2$plot,
                          nrow=2,ncol=4,labels = LETTERS[1:8])
ggsave('results/05.clinical/clinical.km.pdf',clinical.km,height = 9,width = 18)


dir.create('results/06.TME')
load('results/06.TME/tcga.cellscore.RData')
tcga.cell.ssgsea[1:5,1:5]
p1=my_mutiboxplot(t(tcga.cell.ssgsea[,tcga.risktype.cli$Samples]),
                   group = tcga.risktype.cli$Risktype,legend.pos = 'top',group_cols = risktype.col,
                   ylab = 'Score',fill = 'Risktype',angle = 0,hjust = .5)
p1

tcga.estimate=read.delim('results/06.TME/ESTIMATE_score.txt')
head(tcga.estimate)
p2=my_mutiboxplot( tcga.estimate[tcga.risktype.cli$Samples,1:3],
                   group = tcga.risktype.cli$Risktype,legend.pos = 'top',group_cols = risktype.col,
                   ylab = 'Score',fill = 'Risktype',angle = 0,hjust = .5)

p2

p3=my_violin(dat = as.numeric(tcga.ifnγ[tcga.risktype.cli$Samples]),group = tcga.risktype.cli$Risktype,
          group_cols = risktype.col,ylab = 'IFN-γ',legend.position = 'none',xlab = 'Risktype')
p3

load('results/06.TME/tcga.immu.ssgsea.RData')
p4=my_mutiboxplot(tcga.immu.ssgsea[tcga.risktype.cli$Samples,],
                   group = tcga.risktype.cli$Risktype,#test_method = 'wilcox.test',
                   legend.pos = 'none',group_cols = risktype.col,
                   ylab = 'Score')
p4




tcga.mcp=immu_MCPcounter(tcga.exp)
p5=my_mutiboxplot( tcga.mcp[tcga.risktype.cli$Samples,],
                   group = tcga.risktype.cli$Risktype,legend.pos = 'none',group_cols = risktype.col,
                   ylab = 'Score')
p5

icgs.gene=read.xlsx('origin_datas/79ICGs_PMID32814346.xlsx',check.names = F)
head(icgs.gene)
icgs.gene=icgs.gene[,c(1,3)]
table(icgs.gene$Role.with.Immunity)
icgs.gene=icgs.gene[icgs.gene$Symbol%in%rownames(tcga.exp),]
icgs.gene=icgs.gene[icgs.gene$Role.with.Immunity!='TwoSide',]
p6=my_mutiboxplot(t(tcga.exp[icgs.gene$Symbol[icgs.gene$Role.with.Immunity=='Inhibit'],tcga.risktype.cli$Samples]),
                   group = tcga.risktype.cli$Risktype,legend.pos = 'none',group_cols = risktype.col,
                   ylab = 'Gene Expression')
p6

immu.fig=mg_merge_plot(mg_merge_plot(p1,p2,p3,ncol=3,widths = c(1.8,1,0.7),align = 'h',common.legend = T,labels = c('A','B','C')),
                       p4,mg_merge_plot(p5,p6,labels = c('E','F'),align = 'h',widths = c(1.2,2)),
                       nrow = 3,ncol=1,labels = c('','D',''))
ggsave('results/06.TME/immu.fig.pdf',immu.fig,height = 14,width = 20)



dir.create('results/07.pathway')
getGeneFC=function(gene.exp,group,ulab=ulab,dlab=dlab){
  degs_C1_C3=mg_limma_DEG(gene.exp, 
                          group,
                          ulab=ulab,
                          dlab=dlab)
  
 
  degs_C1_C3_sig<-degs_C1_C3$DEG[which(degs_C1_C3$DEG$adj.P.Val <= 1),]
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  degs_C1_C3_sig_gene<-rownames(degs_C1_C3_sig)
  
  degs_gene_entz=bitr(degs_C1_C3_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)
  
  gene_df <- data.frame(logFC=degs_C1_C3_sig$logFC,
                        SYMBOL = rownames(degs_C1_C3_sig))
  gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
  head(gene_df)
  
  geneList<-gene_df$logFC
  names(geneList)=gene_df$ENTREZID 
  head(geneList)
  
  geneList=sort(geneList,decreasing = T)
  head(geneList) 
  return(geneList)
}
tcga.geneList=getGeneFC(gene.exp=tcga.exp,group=tcga.risktype.cli$Risktype
                        ,ulab='High',dlab = 'Low')

h.all.gmt<-read.gmt("origin_datas/h.all.v2023.1.Hs.entrez.gmt")
tcga.hallmark.gsea<-GSEA(tcga.geneList,TERM2GENE = h.all.gmt,seed=T)
library(enrichplot)
library(ggplot2)
library(GseaVis)

tcga.hallmark.gsea.res=tcga.hallmark.gsea@result
write.csv(tcga.hallmark.gsea.res,'results/07.pathway/tcga.hallmark.gsea.res.csv',row.names = F)
gsea.dotplot=dotplotGsea(data = tcga.hallmark.gsea,topn = 10,
                         order.by = 'NES',pajust = 0.05)

p1=gsea.dotplot$plot+theme(text = element_text(family = 'Times'))+
  scale_y_discrete(labels=function(x) str_remove(x,"HALLMARK_"))
pdf('results/07.pathway/GSEA.pdf',height = 6,width = 8)
p1
dev.off()




dir.create('results/08.mutation')
tcga.sub.all=readMatrix('PMC5982584_supplement_2.txt')
table(tcga.sub.all$`TCGA Study`)
tcga.immu.sub<-tcga.sub.all[tcga.sub.all$`TCGA Study`=='ESCA',]
head(tcga.immu.sub)
tcga.immu.sub$Samples=paste0(rownames(tcga.immu.sub),'-01')
rownames(tcga.immu.sub)=tcga.immu.sub$Samples
tcga.immu.sub=merge(tcga.immu.sub,tcga.risktype.cli,by='Samples')
dim(tcga.immu.sub)
# 403  63
colnames(tcga.immu.sub)
col.selected=c('Aneuploidy Score','Homologous Recombination Defects','Fraction Altered','Number of Segments')
p2=my_violin(dat = tcga.immu.sub[,19],group =  tcga.immu.sub$Risktype,group_cols = risktype.col,
          legend.position = 'none',xlab = '',ylab = colnames(tcga.immu.sub)[19])
p3=my_violin(dat = tcga.immu.sub[,20],group =  tcga.immu.sub$Risktype,group_cols = risktype.col,
          legend.position = 'none',xlab = 'Risktype',ylab = colnames(tcga.immu.sub)[20])

mg_merge_plot(p1,mg_merge_plot(p2,p3,nrow = 2,labels = c('B','C')),
              widths = c(2,1),labels = c('A',''))
ggsave('results/07.pathway/GSEA&mutation.pdf',height = 8,width = 14)

tcga.maf=getTCGAMAFByCode('ESCA')
tcga.risktype.cli.use=tcga.risktype.cli
table(tcga.risktype.cli.use$Risktype)
colnames(tcga.risktype.cli.use)[1]='Tumor_Sample_Barcode'
tcga.risktype.cli.use$Tumor_Sample_Barcode=substr(tcga.risktype.cli.use$Tumor_Sample_Barcode,1,12)
tcga.risktype.cli.use.high=tcga.risktype.cli.use[which(tcga.risktype.cli.use$Risktype=='High'),]
tcga.risktype.cli.use.low=tcga.risktype.cli.use[which(tcga.risktype.cli.use$Risktype=='Low'),]

write.table(tcga.risktype.cli.use.high,file='results/08.mutation/tcga.risktype.cli.high.txt')
write.table(tcga.risktype.cli.use.low,file='results/08.mutation/tcga.risktype.cli.low.txt')

tcga.maf1=subsetMaf(tcga.maf,tsb=intersect(tcga.maf@data$Tumor_Sample_Barcode,tcga.risktype.cli.use.high$Tumor_Sample_Barcode))
tcga.maf1<-read.maf(tcga.maf1@data,isTCGA=T,clinicalData = 'results/08.mutation/tcga.risktype.cli.high.txt')
tcga.maf1@clinical.data

tcga.maf2=subsetMaf(tcga.maf,tsb=intersect(tcga.maf@data$Tumor_Sample_Barcode,tcga.risktype.cli.use.low$Tumor_Sample_Barcode))
tcga.maf2<-read.maf(tcga.maf2@data,isTCGA=T,clinicalData = 'results/08.mutation/tcga.risktype.cli.low.txt')
tcga.maf2@clinical.data


mdf_cp=mafCompare(tcga.maf1,tcga.maf2,m1Name = 'High',m2Name = 'Low')
mdf_cp$results
write.csv(mdf_cp$results,'results/08.mutation/maf_compare.csv')

pdf('results/08.mutation/maf_compare.pdf',height = 5,width = 7,onefile = F)
forestPlot(mafCompareRes = mdf_cp,color = c('royalblue', 'maroon'))
dev.off()

maf.res=mdf_cp$results
maf.res=maf.res[maf.res$pval<0.05,]
maf.res=maf.res[order(maf.res$High,decreasing = T),]
pdf('results/08.mutation/maf_oncoplot.pdf',height = 5,width = 7)
coOncoplot(m1=tcga.maf1, m2=tcga.maf2, m1Name = 'High',m2Name = 'Low',
           genes =  maf.res$Hugo_Symbol)
dev.off()

save.image(file='bulk_project.RData')
