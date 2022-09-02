library(Seurat)
library(ggplot2)
library(ggrepel)
library(scales)
library(ggpubr)

### Function-----
plotfigure<-function(object=NULL,feature=NULL,x_adj=F){
  if(grepl("TotalSeqC",feature)){
    D<-data.frame(row.names  = colnames(object),
                  IGHD=object@assays$Protein@data[feature,],
                  celltype=Idents(object))
    feature_name<-gsub("-","_",feature)
    colnames(D)<-c(feature_name,"cluster")
    # D$cluster<-factor(D$cluster,levels = sort(as.numeric(levels(D$cluster))))
  } 
  else if (feature %in% rownames(object@assays$RNA@counts)){
    D<-data.frame(row.names  = colnames(object),
                  IGHD=object@assays$RNA@data[feature,],
                  celltype=Idents(object))
    feature_name<-gsub("-","_",feature)
    colnames(D)<-c(feature_name,"cluster")
    # D$cluster<-factor(D$cluster,levels = sort(as.numeric(levels(D$cluster))))
    
  }
  
  else {
    D<-data.frame(row.names  = colnames(object),
                  IGHD=object[[feature]],
                  celltype=Idents(object))
    feature_name<-gsub("-","_",feature)
    colnames(D)<-c(feature_name,"cluster")
    # D$cluster<-factor(D$cluster,levels = sort(as.numeric(levels(D$cluster))))
  }
  # return(D)
  if(x_adj){
    ggplot(D, aes_string(x="cluster", y=feature_name,color="cluster",fill="cluster")) +
      geom_jitter(shape=16, position=position_jitter(0.2),size=0.8,color="black")+
      geom_violin(alpha=0.5,scale = "width",size=1)+
      stat_boxplot(geom="errorbar",width=0.1,color="red",size=0.8)+
      geom_boxplot(width=0.2,color="red",outlier.alpha = 0)+
      geom_hline(yintercept=median(D[,feature_name]),color="grey",size=1)+
      stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.")+
      theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(),axis.line = element_line(colour = "black"),
            text = element_text(size = 15,family = "ArialMT"),axis.text.x = element_text(vjust=0.5,angle = 45),plot.title = element_text(hjust = 0.5),
            legend.position = "none")+
      ylab("Expression level")+ggtitle(feature_name)+xlab("")
  }
  else{
    ggplot(D, aes_string(x="cluster", y=feature_name,color="cluster",fill="cluster")) +
      geom_jitter(shape=16, position=position_jitter(0.2),size=0.8,color="black")+
      geom_violin(alpha=0.5,scale = "width",size=1)+
      stat_boxplot(geom="errorbar",width=0.1,color="red",size=0.8)+
      geom_boxplot(width=0.2,color="red",outlier.alpha = 0)+
      geom_hline(yintercept=median(D[,feature_name]),color="grey",size=1)+
      stat_compare_means(label = "p.signif", method = "t.test",ref.group = ".all.")+   
      theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(),axis.line = element_line(colour = "black"),
            text = element_text(size = 15,family = "ArialMT"),plot.title = element_text(hjust = 0.5),
            legend.position = "none")+
      ylab("Expression level")+ggtitle(feature_name)+xlab("")
  }
}

### Load the PBMC samples----
PBMC1<-Read10X("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/Data/supperseq/PBMC1/GEX+ADT/filtered")
PBMC2<-Read10X("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/Data/supperseq/PBMC2/GEX+ADT/filtered")
PBMC3<-Read10X("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/Data/supperseq/PBMC3/GEX+ADT/filtered")
colnames(PBMC1$`Gene Expression`)<-sub("-1","_1",colnames(PBMC1$`Gene Expression`))
colnames(PBMC2$`Gene Expression`)<-sub("-1","_2",colnames(PBMC2$`Gene Expression`))
colnames(PBMC3$`Gene Expression`)<-sub("-1","_3",colnames(PBMC3$`Gene Expression`))

PBMC1_object<-CreateSeuratObject(counts = PBMC1$`Gene Expression`, project = "PBMC")
PBMC2_object<-CreateSeuratObject(counts = PBMC2$`Gene Expression`, project = "PBMC")
PBMC3_object<-CreateSeuratObject(counts = PBMC3$`Gene Expression`, project = "PBMC")

Total_PBMC<-merge(x=PBMC1_object,y = c(PBMC2_object,PBMC3_object))
rm(PBMC1,PBMC2,PBMC3,PBMC1_object,PBMC2_object,PBMC3_object)
Total_PBMC[['percent.RIB']]<-PercentageFeatureSet(Total_PBMC, pattern = "^RP[SL]")




### Load the bone marrow samples----
BM2<-Read10X("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/Data/supperseq/BM2/GEX+ADT/filtered")
BM3<-Read10X("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/Data/supperseq/BM3/GEX+ADT/filtered")
colnames(BM2$`Gene Expression`)<-sub("-1","_1",colnames(BM2$`Gene Expression`))
colnames(BM3$`Gene Expression`)<-sub("-1","_2",colnames(BM3$`Gene Expression`))
BM2_object<-CreateSeuratObject(counts = BM2$`Gene Expression`, project = "BM")
BM3_object<-CreateSeuratObject(counts = BM3$`Gene Expression`, project = "BM")
Total_BM<-merge(x=BM2_object,y = BM3_object)
rm(BM2,BM3,BM2_object,BM3_object)
Total_BM[['percent.RIB']]<-PercentageFeatureSet(Total_BM, pattern = "^RP[SL]")

### Load gating strategy----
PBMCs_major_lineages <- readRDS("~/GitHub/SuPERR-Seq/data/Major lineages/PBMCs_major_lineages.rds")
BMs_major_lineages <- readRDS("~/GitHub/SuPERR-Seq/data/Major lineages/BMs_major_lineages.rds")
names(PBMCs_major_lineages)<-sub("PBMC","",names(PBMCs_major_lineages))
names(BMs_major_lineages)<-sub("BM","",names(BMs_major_lineages))
Total_PBMC<-subset(Total_PBMC,cells = names(PBMCs_major_lineages))
Total_BM<-subset(Total_BM,cells = names(BMs_major_lineages))
# cell_classification<-factor(PBMCs_major_lineages,levels = c(""))
Idents(Total_PBMC)<-PBMCs_major_lineages
Idents(Total_BM)<-BMs_major_lineages
# A<-subset(Total_PBMC,idents = a)

### Figure 4 ABC----
### PBMC
setwd("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/Junkai/Junkai_update/PBMC_variation")
RNA_PBMC<-plotfigure(Total_PBMC,feature = "nCount_RNA")+
  ylab("UMI count")+ggtitle("")#+stat_compare_means(method = "anova",label.x=7,label.y = 37000)
ggsave(filename = "PBMC_nUMI.pdf",width = 118.5,height = 75.1,units = "mm",device = "pdf",plot = RNA_PBMC,scale = 1.4)
RP_PBMC<-plotfigure(Total_PBMC,feature = "percent.RIB")+
  ylab("Percent of Ribosomal")+ggtitle("")#+stat_compare_means(method = "anova",label.y = 60,label.x=7 )
ggsave(filename = "PBMC_RP.pdf",width = 118.5,height = 75.1,units = "mm",device = "pdf",plot = RP_PBMC,scale = 1.4)
nGene_PBMC<-plotfigure(Total_PBMC,feature = "nFeature_RNA")+
  ylab("Gene count")+ggtitle("")#+stat_compare_means(method = "anova",label.x=7,label.y = 4400)
ggsave(filename = "PBMC_nGene.pdf",width = 118.5,height = 75.1,units = "mm",device = "pdf",plot = nGene_PBMC,scale = 1.4)

### BM----
setwd("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/Junkai/Junkai_update/BM_varitaion")
RNA_BM<-plotfigure(Total_BM,feature = "nCount_RNA")+
  ylab("UMI count")+ggtitle("")#+stat_compare_means(method = "anova",label.x=4,label.y = 67000)
ggsave(filename = "BM_nUMI.pdf",width = 118.5,height = 75.1,units = "mm",device = "pdf",plot = RNA_BM,scale = 1.4)
RP_BM<-plotfigure(Total_BM,feature = "percent.RIB")+
  ylab("Percent of Ribosomal")+ggtitle("")#+stat_compare_means(method = "anova",label.y = 62,label.x=4 )
ggsave(filename = "BM_RP.pdf",width = 118.5,height = 75.1,units = "mm",device = "pdf",plot = RP_BM,scale = 1.4)
nGene_BM<-plotfigure(Total_BM,feature = "nFeature_RNA")+
  ylab("Gene count")+ggtitle("")#+stat_compare_means(method = "anova",label.x=4,label.y = 6400)
ggsave(filename = "BM_nGene.pdf",width = 118.5,height = 75.1,units = "mm",device = "pdf",plot = nGene_BM,scale = 1.4)


### Figure 4 D----
setwd("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/Junkai/Junkai_update/PBMC_singlet/HVGs")
# This is the HVGs of the integrated total PBMC from Seurat
Total_gene<-read.csv("PBMC_HVGs_meta_all_cells.csv")
Total_gene<-Total_gene[order(Total_gene$vst.variance.standardized,decreasing = T),]
rownames(Total_gene)<-NULL

# This is the HVGs of the integrated PBMC B cells from Seurat
B_gene<-read.csv("PBMC_HVGS_meta_total_B_cells.csv")
B_gene<-B_gene[order(B_gene$vst.variance.standardized,decreasing = T),]
rownames(B_gene)<-NULL

big_gene<-union(Total_gene$X[1:300],B_gene$X[1:30])
refined_gene<-Total_gene$X[1:300]
Total_gene<-Total_gene[Total_gene$X %in% big_gene,]

### B on Total
options(ggrepel.max.overlaps = Inf)
setwd("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/Junkai/Junkai_update/PBMC_B/HVG")
a<-ggplot(Total_gene, aes(x= vst.mean, y = vst.variance.standardized))+
  geom_point(color = "grey", size = 3)+
  geom_point(data=Total_gene[Total_gene$X %in% intersect(Total_gene$X,B_gene$X[1:30]),], aes(x=vst.mean, y=vst.variance.standardized), color="green", size=3)+
  geom_point(data=Total_gene[Total_gene$X %in% setdiff(big_gene,refined_gene),], aes(x=vst.mean, y=vst.variance.standardized), color="green", size=3)+
  geom_label_repel(data=Total_gene[Total_gene$X %in% B_gene$X[1:30],],aes(label = X),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  scale_y_continuous(trans='log10',
                     breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x)),limits = c(10^-1,10^3)) +
  scale_x_continuous(trans='log10',
                     breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x)),limits = c(10^-3,10^1)) +
  xlab("Average expression")+ylab("Standardized variance")+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        text = element_text(size = 15,family = "ArialMT"),plot.title = element_text(hjust = 0.5))+annotation_logticks()

ggsave(filename = "B_HVG.pdf",width = 118.5,height = 80,units = "mm",device = "pdf",scale = 1.4,plot = a)

### Top 10 PBMC
b<-ggplot(Total_gene[1:300,], aes(x= vst.mean, y = vst.variance.standardized))+
  geom_point(color = "grey", size = 3)+
  geom_point(data=Total_gene[1:30,], aes(x=vst.mean, y=vst.variance.standardized), color="red", size=3)+
  geom_label_repel(data=Total_gene[1:30,],aes(label = X),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')+
  scale_y_continuous(trans='log10',
                     breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x)),limits = c(10^-1,10^3)) +
  scale_x_continuous(trans='log10',
                     breaks=trans_breaks('log10', function(x) 10^x),
                     labels=trans_format('log10', math_format(10^.x)),limits = c(10^-3,10^1)) +
  xlab("Average expression")+ylab("Standardized variance")+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"),
        text = element_text(size = 15,family = "ArialMT"),plot.title = element_text(hjust = 0.5))+annotation_logticks()
ggsave(filename = "Total_HVG.pdf",width = 118.5,height = 80,units = "mm",device = "pdf",scale = 1.4,plot = b)
