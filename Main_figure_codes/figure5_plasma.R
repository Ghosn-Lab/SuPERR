library(Seurat)
library(ggpubr)

plotfigure<-function(object=NULL,feature=NULL,x_adj=F){
  if(grepl("TotalSeqC",feature)){
    D<-data.frame(row.names  = colnames(object),
                  IGHD=object@assays$Protein@data[feature,],
                  celltype=Idents(object))
    feature_name<-gsub("-","_",feature)
    colnames(D)<-c(feature_name,"cluster")
    D$cluster<-factor(D$cluster,levels = sort(as.numeric(levels(D$cluster))))
  } 
  else if (feature %in% rownames(object@assays$RNA@counts)){
    D<-data.frame(row.names  = colnames(object),
                  IGHD=object@assays$RNA@data[feature,],
                  celltype=Idents(object))
    feature_name<-gsub("-","_",feature)
    colnames(D)<-c(feature_name,"cluster")
    D$cluster<-factor(D$cluster,levels = sort(as.numeric(levels(D$cluster))))
    
  }
  
  else {
    D<-data.frame(row.names  = colnames(object),
                  IGHD=object[[feature]],
                  celltype=Idents(object))
    feature_name<-gsub("-","_",feature)
    colnames(D)<-c(feature_name,"cluster")
    D$cluster<-factor(D$cluster,levels = sort(as.numeric(levels(D$cluster))))
  }
  # return(D)
  if(x_adj){
    ggplot(D, aes_string(x="cluster", y=feature_name,color="cluster",fill="cluster")) +
      geom_jitter(shape=16, position=position_jitter(0.2),size=0.8,color="black")+
      geom_violin(alpha=0.5,scale = "width",size=1)+
      stat_boxplot(geom="errorbar",width=0.1,color="red",size=0.8)+
      geom_boxplot(width=0.2,color="red",outlier.alpha = 0)+
      geom_hline(yintercept=median(D[,feature_name]),color="grey",size=1)+
      theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),legend.position = "None",
            panel.background = element_blank(),axis.line = element_line(colour = "black"),
            text = element_text(size = 15,family = "ArialMT"),axis.text.x = element_text(vjust=0.5,angle = 90,hjust = 1),plot.title = element_text(hjust = 0.5))+
      ylab("Expression level")+ggtitle(feature_name)+xlab("")
  }
  else{
    ggplot(D, aes_string(x="cluster", y=feature_name,color="cluster",fill="cluster")) +
      geom_jitter(shape=16, position=position_jitter(0.2),size=0.8,color="black")+
      geom_violin(alpha=0.5,scale = "width",size=1)+
      stat_boxplot(geom="errorbar",width=0.1,color="red",size=0.8)+
      geom_boxplot(width=0.2,color="red",outlier.alpha = 0)+
      geom_hline(yintercept=median(D[,feature_name]),color="grey",size=1)+
      theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),legend.position = "None",
            panel.background = element_blank(),axis.line = element_line(colour = "black"),
            text = element_text(size = 15,family = "ArialMT"),plot.title = element_text(hjust = 0.5))+ylab("Expression level")+ggtitle(feature_name)+xlab("")
  }
}


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

### Load the clustering and UMAP results from the previous steps----

###Load S3 UMAP
setwd("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/UMAP_Embedding/singlet_2lowmt")

UMAP_cord<-read.csv("Seurat_cell_umap_loading_total_plasma_cells.csv",row.names = "X")
rownames(UMAP_cord)<-sub("BM2","1",rownames(UMAP_cord))
rownames(UMAP_cord)<-sub("BM3","2",rownames(UMAP_cord))
Total_BM<-subset(Total_BM,cells = rownames(UMAP_cord))
Total_BM <- NormalizeData(Total_BM, normalization.method = "LogNormalize", scale.factor = 10000)
Total_BM <- FindVariableFeatures(Total_BM, selection.method = "vst", nfeatures = 2000)
Total_BM <- ScaleData(Total_BM, features = VariableFeatures(object = Total_BM),verbose = F)
Total_BM <- RunPCA(Total_BM, features = VariableFeatures(object = Total_BM),verbose = F)
Total_BM <- RunUMAP(Total_BM, dims = 1:10,verbose = F)

colnames(UMAP_cord)<-colnames(Total_BM@reductions$umap@cell.embeddings)
Total_BM@reductions$umap@cell.embeddings<-as.matrix(UMAP_cord)

###Load Seurat3 clusters
setwd("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/superrseq_clusters/BM/singlet_2lowmt")
S3_cluster_138hi<-read.csv("Seurat_clustering_result_total_CD138high_plasma_cells.csv",row.names = "X")
S3_cluster_138lo<-read.csv("Seurat_clustering_result_total_CD138low_plasma_cells.csv",row.names = "X")
S3_cluster_138lo$x<-S3_cluster_138lo$x+1
S3_cluster_138hi$x<-S3_cluster_138hi$x+3
S3_cluster<-rbind(S3_cluster_138hi,S3_cluster_138lo)
rownames(S3_cluster)<-sub("BM2","1",rownames(S3_cluster))
rownames(S3_cluster)<-sub("BM3","2",rownames(S3_cluster))

cluster_label<-factor(S3_cluster$x,levels = 1:max(S3_cluster$x))
names(cluster_label)<-rownames(S3_cluster)
Total_BM$S3_cluster<-cluster_label
Idents(Total_BM)<-Total_BM$S3_cluster

# Calculate percent of Ig genes
Total_BM[['percent.Ig']]<-PercentageFeatureSet(Total_BM, pattern = "^IG[HLK]V")
#Calculat PC-scores after removing Ig[HKL]V genes
GEX_matrix<-Total_BM@assays$RNA@counts
GEX_matrix<-GEX_matrix[-grep("^IG[HLK]V",rownames(GEX_matrix)),]
Plasma_cell_genes<-c("ITGB7", "IRF4", "CD9", "PRDM1", "XBP1", "SDC1", "VCAM1", "CD38")
PC_expreesion<-colSums(GEX_matrix[Plasma_cell_genes,])
Total_BM[["PC_scores"]]<-log1p(10000*PC_expreesion/colSums(GEX_matrix))
#Renormalize data after removing Ig[HKL]V genes
Total_BM<-subset(Total_BM,features = rownames(GEX_matrix))
Total_BM <- NormalizeData(Total_BM, normalization.method = "LogNormalize", scale.factor = 10000)

### pannel a----
setwd("C:\\Users\\jyan399\\OneDrive - Emory University\\Box\\Junkai Yang Ghosn Lab\\updates\\SuPERR-seq\\Junkai\\Junkai_update\\BM_plasma\\figure6")
DimPlot(Total_BM,label = T,pt.size = 3,label.size = 6)
ggsave(filename = "UMAP_plasma.pdf",width = 118.5,height = 100,units = "mm",device = "pdf",scale = 1.4)

### pannel b----
plotfigure(Total_BM,feature = "percent.Ig")+ylab("Percent\nIg-specific genes (%)")+
  ggtitle("Before removing Ig-specific genes")
ggsave(filename = "violin_Ig.pdf",width =150,height = 70,units = "mm",device = "pdf",scale = 1.4)

plotfigure(Total_BM,feature = "PC_scores")+ylab("PC-scores")+
  ggtitle("After removing Ig-specific genes")
ggsave(filename = "violin_plasma.pdf",width = 150,height = 70,units = "mm",device = "pdf",scale = 1.4)

### pannel c----
plotfigure(Total_BM,feature = "ITGB7")
ggsave(filename = "violin_ITGB7.pdf",width = 55,height = 65,units = "mm",device = "pdf",scale = 1.4)

plotfigure(Total_BM,feature = "IRF4")
ggsave(filename = "violin_IRF4.pdf",width = 55,height = 65,units = "mm",device = "pdf",scale = 1.4)

plotfigure(Total_BM,feature = "CD9")
ggsave(filename = "violin_CD9.pdf",width = 55,height = 65,units = "mm",device = "pdf",scale = 1.4)

plotfigure(Total_BM,feature = "PRDM1")
ggsave(filename = "violin_PRDM1.pdf",width = 55,height = 65,units = "mm",device = "pdf",scale = 1.4)

plotfigure(Total_BM,feature = "XBP1")
ggsave(filename = "violin_XBP1.pdf",width = 55,height = 65,units = "mm",device = "pdf",scale = 1.4)

plotfigure(Total_BM,feature = "SDC1")
ggsave(filename = "violin_SDC1.pdf",width = 55,height = 65,units = "mm",device = "pdf",scale = 1.4)

plotfigure(Total_BM,feature = "HLA-DR-TotalSeqC")
ggsave(filename = "violin_HLA_DR.pdf",width = 55,height = 65,units = "mm",device = "pdf",scale = 1.4)

plotfigure(Total_BM,feature = "log_g2m")
ggsave(filename = "violin_g2m.pdf",width = 55,height = 65,units = "mm",device = "pdf",scale = 1.4)

plotfigure(Total_BM,feature = "PDL1-TotalSeqC")
ggsave(filename = "violin_PDL1ADT.pdf",width = 55,height = 65,units = "mm",device = "pdf",scale = 1.4)

plotfigure(Total_BM,feature = "CD138-TotalSeqC")
ggsave(filename = "violin_CD138ADT.pdf",width = 55,height = 65,units = "mm",device = "pdf",scale = 1.4)

setwd("H:/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/Junkai/Junkai_update/BM_plasma/figureS9")
plotfigure(Total_BM,feature = "HLA-DR-TotalSeqC")
ggsave(filename = "violin_HLA_DR.pdf",width = 100,height = 80,units = "mm",device = "pdf",scale = 1.4)

plotfigure(Total_BM,feature = "log_g2m")
ggsave(filename = "violin_g2m.pdf",width = 100,height = 80,units = "mm",device = "pdf",scale = 1.4)


