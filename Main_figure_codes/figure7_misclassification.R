library(Seurat)
library(ggpubr)
library(VennDiagram)

### Function----
source("Ben_CFS.R")
misclassification<-function(object=NULL,cell_type=NULL,clusters=NULL){
  classification<-data.frame(row.names = colnames(object),
                             ADT=object@meta.data[,cell_type],
                             cluster=object@meta.data[,clusters],
                             misclassification=FALSE)
  
  for (i in 1:max(as.numeric(levels(classification$cluster)))){
    sum<-table(classification$ADT[classification$cluster==i])
    if (sum(sum)==0) next
    cell_name<-names(sum)[sum==max(sum)][1]
    classification$misclassification[classification$cluster==i&classification$ADT!=cell_name]<-TRUE
  }
  return(classification)
}

###Load gating strategy and samples----
# PBMCs
PBMCs_major_lineages <- readRDS("~/GitHub/SuPERR-Seq/data/Major lineages/PBMCs_major_lineages.rds")
BMs_major_lineages <- readRDS("~/GitHub/SuPERR-Seq/data/Major lineages/BMs_major_lineages.rds")


PBMC1<-Read10X("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/Data/supperseq/PBMC1/GEX+ADT/filtered")
PBMC2<-Read10X("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/Data/supperseq/PBMC2/GEX+ADT/filtered")
PBMC3<-Read10X("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/Data/supperseq/PBMC3/GEX+ADT/filtered")
colnames(PBMC1$`Gene Expression`)<-sub("-1","_1",colnames(PBMC1$`Gene Expression`))
colnames(PBMC2$`Gene Expression`)<-sub("-1","_2",colnames(PBMC2$`Gene Expression`))
colnames(PBMC3$`Gene Expression`)<-sub("-1","_3",colnames(PBMC3$`Gene Expression`))
DSB_PBMC1<-read.csv("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/Data/supperseq/DSB_PBMC1.csv",row.names = "X")
DSB_PBMC2<-read.csv("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/Data/supperseq/DSB_PBMC2.csv",row.names = "X")
DSB_PBMC3<-read.csv("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/Data/supperseq/DSB_PBMC3.csv",row.names = "X")

PBMC1_object<-CreateSeuratObject(counts = PBMC1$`Gene Expression`, project = "PBMC")
PBMC2_object<-CreateSeuratObject(counts = PBMC2$`Gene Expression`, project = "PBMC")
PBMC3_object<-CreateSeuratObject(counts = PBMC3$`Gene Expression`, project = "PBMC")
PBMC1_object[["batch"]]<-"PBMC1"
PBMC2_object[["batch"]]<-"PBMC2"
PBMC3_object[["batch"]]<-"PBMC3"
PBMC1_object[["Protein"]]<-CreateAssayObject(counts = DSB_PBMC1)
PBMC2_object[["Protein"]]<-CreateAssayObject(counts = DSB_PBMC2)
PBMC3_object[["Protein"]]<-CreateAssayObject(counts = DSB_PBMC3)

Total_PBMC<-merge(x=PBMC1_object,y = c(PBMC2_object,PBMC3_object))
rm(PBMC1,PBMC2,PBMC3,DSB_PBMC1,DSB_PBMC2,DSB_PBMC3,PBMC1_object,PBMC2_object,PBMC3_object)

Total_PBMC<-subset(Total_PBMC,cells = rownames(UMAP_cord))
Total_PBMC <- NormalizeData(Total_PBMC, normalization.method = "LogNormalize", scale.factor = 10000)
Total_PBMC <- FindVariableFeatures(Total_PBMC, selection.method = "vst", nfeatures = 2000)
Total_PBMC <- ScaleData(Total_PBMC, features = VariableFeatures(object = Total_PBMC),verbose = F)
Total_PBMC <- RunPCA(Total_PBMC, features = VariableFeatures(object = Total_PBMC),verbose = F)
Total_PBMC <- RunUMAP(Total_PBMC, dims = 1:10,verbose = F)

Total_PBMC<-subset(Total_PBMC,cells = names(PBMCs_major_lineages))
# Merge Treg with CD4T becaue Treg is a subset of CD4 T cell. 
# Separating it generates more misclassification for the tested tools
PBMCs_major_lineages[PBMCs_major_lineages=="Treg"]<-"CD4T"
PBMCs_major_lineages<-droplevels(PBMCs_major_lineages,"Treg")
Total_PBMC$SuPERR_cells<-PBMCs_major_lineages

#BMs
BM2<-Read10X("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/Data/supperseq/BM2/GEX+ADT/filtered")
BM3<-Read10X("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/Data/supperseq/BM3/GEX+ADT/filtered")
colnames(BM2$`Gene Expression`)<-sub("-1","_1",colnames(BM2$`Gene Expression`))
colnames(BM3$`Gene Expression`)<-sub("-1","_2",colnames(BM3$`Gene Expression`))
BM2_object<-CreateSeuratObject(counts = BM2$`Gene Expression`, project = "BM")
BM3_object<-CreateSeuratObject(counts = BM3$`Gene Expression`, project = "BM")
BM2_object[["batch"]]<-"BM1"
BM3_object[["batch"]]<-"BM2"
Total_BM<-merge(x=BM2_object,y = BM3_object)
rm(BM2,BM3,BM2_object,BM3_object)

Total_BM <- NormalizeData(Total_BM, normalization.method = "LogNormalize", scale.factor = 10000)
Total_BM <- FindVariableFeatures(Total_BM, selection.method = "vst", nfeatures = 2000)
Total_BM <- ScaleData(Total_BM, features = VariableFeatures(object = Total_BM),verbose = F)
Total_BM <- RunPCA(Total_BM, features = VariableFeatures(object = Total_BM),verbose = F)
Total_BM <- RunUMAP(Total_BM, dims = 1:10,verbose = F)

Total_BM<-subset(Total_BM,cells = names(BMs_major_lineages))

Total_BM$SuPERR_cells<-BMs_major_lineages

### Pannel A----

# Load Seurat3 clusters
S3_cluster<-read.csv("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/superrseq_clusters/PBMC/singlet_2lowmt/Seurat_clustring_result_Anchor_PBMC_v1.csv",row.names = "X")
S3_cluster$x<-S3_cluster$x+1
rownames(S3_cluster)<-sub("PBMC","",rownames(S3_cluster))
cluster_label<-factor(S3_cluster$x,levels = 1:max(S3_cluster$x))
names(cluster_label)<-rownames(S3_cluster)
Total_PBMC$S3_cluster<-cluster_label

# Load Seurat4 clusters
S4_cluster<-read.csv("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/Junkai/Junkai_update/Seurat4/Seurat_clustring_result_Anchor_wnn_PBMC.csv",row.names = "X")
S4_cluster$x<-S4_cluster$x+1
rownames(S4_cluster)<-sub("PBMC","",rownames(S4_cluster))
cluster_label<-factor(S4_cluster$x,levels = 1:max(S4_cluster$x))
names(cluster_label)<-rownames(S4_cluster)
Total_PBMC$S4_cluster<-cluster_label

# Load CiteFuse clustes
Cite_cluster<-read.csv("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/Junkai/Junkai_update/Seurat4/CiteFuse_hires_clustering_result_PBMC_v1.csv",row.names = "X")
rownames(Cite_cluster)<-sub("PBMC","",rownames(Cite_cluster))
Cite_cluster<-Cite_cluster[rownames(Cite_cluster)%in%colnames(Total_PBMC),,drop=F]
cluster_label<-factor(Cite_cluster$SNF_W_louvain,levels = 1:max(Cite_cluster$SNF_W_louvain))
names(cluster_label)<-rownames(Cite_cluster)
Total_PBMC$Cite_cluster<-cluster_label

# Load Seurat3 UMAP and calculate CMS----
UMAP_cord<-read.csv("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/UMAP_Embedding/singlet_2lowmt/umap_cell_embeddings_Anchor_PBMC_v1.csv",row.names = "X")
rownames(UMAP_cord)<-sub("PBMC","",rownames(UMAP_cord))
colnames(UMAP_cord)<-colnames(Total_PBMC@reductions$umap@cell.embeddings)
Total_PBMC@reductions$umap@cell.embeddings<-as.matrix(UMAP_cord)

CMS<-1-getCFS(Total_PBMC,ident.1 = "SuPERR_cells",ident.2 = "S3_cluster")
misclass_cells<-misclassification(Total_PBMC,cell_type = "SuPERR_cells",clusters = "S3_cluster")

DimPlot(Total_PBMC,label = T,
        cells.highlight = rownames(misclass_cells)[misclass_cells$misclassification==TRUE],repel = T,
        pt.size = 1,label.size = 5)+theme(legend.position = "")+ggtitle(CMS)
ggsave(filename = "PBMC_misclassify_S3.pdf",width = 118.5,height = 100.1,units = "mm",device = "pdf",scale = 1.4)

# Load Seurat4 UMAP and calculate CMS----
UMAP_cord<-read.csv("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/Junkai/Junkai_update/Seurat4/umap_cell_embeddings_Anchor_wnn_PBMC.csv",row.names = "X")
rownames(UMAP_cord)<-sub("PBMC","",rownames(UMAP_cord))
colnames(UMAP_cord)<-colnames(Total_PBMC@reductions$umap@cell.embeddings)
Total_PBMC@reductions$umap@cell.embeddings<-as.matrix(UMAP_cord)

CMS<-1-getCFS(Total_PBMC,ident.1 = "SuPERR_cells",ident.2 = "S4_cluster")
misclass_cells<-misclassification(Total_PBMC,cell_type = "SuPERR_cells",clusters = "S4_cluster")

DimPlot(Total_PBMC,label = T,
        cells.highlight = rownames(misclass_cells)[misclass_cells$misclassification==TRUE],repel = T,
        pt.size = 1,label.size = 5)+theme(legend.position = "")+ggtitle(CMS)
ggsave(filename = "PBMC_misclassify_S4.pdf",width = 118.5,height = 100.1,units = "mm",device = "pdf",scale = 1.4)

# Load CiteFuse UMAP and calculate CMS----
UMAP_cord<-read.csv("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/Junkai/Junkai_update/Seurat4/CiteFuse_umap_PBMC.csv",row.names = "X")
rownames(UMAP_cord)<-sub("PBMC","",rownames(UMAP_cord))
UMAP_cord<-UMAP_cord[rownames(UMAP_cord)%in%colnames(Total_PBMC),]
rownames(UMAP_cord)<-sub("PBMC","",rownames(UMAP_cord))
colnames(UMAP_cord)<-colnames(Total_PBMC@reductions$umap@cell.embeddings)
Total_PBMC@reductions$umap@cell.embeddings<-as.matrix(UMAP_cord)

CMS<-1-getCFS(Total_PBMC,ident.1 = "SuPERR_cells",ident.2 = "Cite_cluster")
misclass_cells<-misclassification(Total_PBMC,cell_type = "SuPERR_cells",clusters = "Cite_cluster")

DimPlot(Total_PBMC,label = T,
        cells.highlight = rownames(misclass_cells)[misclass_cells$misclassification==TRUE],repel = T,
        pt.size = 1,label.size = 5)+theme(legend.position = "")+ggtitle(CMS)

ggsave(filename = "PBMC_misclassify_CiteFuse.pdf",width = 118.5,height = 100.1,units = "mm",device = "pdf",scale = 1.4)



### Pannel B----
# Load Seurat3 clusters
S3_cluster<-read.csv("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/superrseq_clusters/BM/singlet_2lowmt/Seurat_clustring_result_Anchor_BM.csv",row.names = "X")
S3_cluster$x<-S3_cluster$x+1
rownames(S3_cluster)<-sub("BM2","1",rownames(S3_cluster))
rownames(S3_cluster)<-sub("BM3","2",rownames(S3_cluster))
cluster_label<-factor(S3_cluster$x,levels = 1:max(S3_cluster$x))
names(cluster_label)<-rownames(S3_cluster)
Total_BM$S3_cluster<-cluster_label

# Load Seurat4 clusters
S4_cluster<-read.csv("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/Junkai/Junkai_update/Seurat4/Seurat_clustring_result_Anchor_wnn_BM.csv",row.names = "X")
S4_cluster$x<-S4_cluster$x+1
rownames(S4_cluster)<-sub("BM2","1",rownames(S4_cluster))
rownames(S4_cluster)<-sub("BM3","2",rownames(S4_cluster))
cluster_label<-factor(S4_cluster$x,levels = 1:max(S4_cluster$x))
names(cluster_label)<-rownames(S4_cluster)
Total_BM$S4_cluster<-cluster_label

# Load CiteFuse clustes
Cite_cluster<-read.csv("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/Junkai/Junkai_update/Seurat4/CiteFuse_clustering_result_BM.csv",row.names = "X")
rownames(Cite_cluster)<-sub("BM2","1",rownames(Cite_cluster))
rownames(Cite_cluster)<-sub("BM3","2",rownames(Cite_cluster))
Cite_cluster<-Cite_cluster[rownames(Cite_cluster)%in%colnames(Total_BM),,drop=F]
cluster_label<-factor(Cite_cluster$SNF_W_louvain,levels = 1:max(Cite_cluster$SNF_W_louvain))
names(cluster_label)<-rownames(Cite_cluster)
Total_BM$Cite_cluster<-cluster_label




# Load Seurat3 UMAP and calculate CMS----
UMAP_cord<-read.csv("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/UMAP_Embedding/singlet_2lowmt/umap_cell_embeddings_Anchor_BM.csv",row.names = "X")
rownames(UMAP_cord)<-sub("BM2","1",rownames(UMAP_cord))
rownames(UMAP_cord)<-sub("BM3","2",rownames(UMAP_cord))
colnames(UMAP_cord)<-colnames(Total_BM@reductions$umap@cell.embeddings)
Total_BM@reductions$umap@cell.embeddings<-as.matrix(UMAP_cord)

CMS<-1-getCFS(Total_BM,ident.1 = "SuPERR_cells",ident.2 = "S3_cluster")
misclass_cells<-misclassification(Total_BM,cell_type = "SuPERR_cells",clusters = "S3_cluster")

DimPlot(Total_BM,label = T,
        cells.highlight = rownames(misclass_cells)[misclass_cells$misclassification==TRUE],repel = T,
        pt.size = 1,label.size = 5)+theme(legend.position = "")+ggtitle(CMS)
ggsave(filename = "PBMC_misclassify_S3.pdf",width = 118.5,height = 100.1,units = "mm",device = "pdf",scale = 1.4)

# Load Seurat4 UMAP and calculate CMS----
UMAP_cord<-read.csv("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/Junkai/Junkai_update/Seurat4/umap_cell_embeddings_Anchor_wnn_BM.csv",row.names = "X")
rownames(UMAP_cord)<-sub("BM2","1",rownames(UMAP_cord))
rownames(UMAP_cord)<-sub("BM3","2",rownames(UMAP_cord))
colnames(UMAP_cord)<-colnames(Total_BM@reductions$umap@cell.embeddings)
Total_BM@reductions$umap@cell.embeddings<-as.matrix(UMAP_cord)

CMS<-1-getCFS(Total_BM,ident.1 = "SuPERR_cells",ident.2 = "S4_cluster")
misclass_cells<-misclassification(Total_BM,cell_type = "SuPERR_cells",clusters = "S4_cluster")

DimPlot(Total_BM,label = T,
        cells.highlight = rownames(misclass_cells)[misclass_cells$misclassification==TRUE],repel = T,
        pt.size = 1,label.size = 5)+theme(legend.position = "")+ggtitle(CMS)
ggsave(filename = "PBMC_misclassify_S4.pdf",width = 118.5,height = 100.1,units = "mm",device = "pdf",scale = 1.4)

# Load CiteFuse UMAP and calculate CMS----
UMAP_cord<-read.csv("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/Junkai/Junkai_update/Seurat4/CiteFuse_umap_BM.csv",row.names = "X")
rownames(UMAP_cord)<-sub("BM2","1",rownames(UMAP_cord))
rownames(UMAP_cord)<-sub("BM3","2",rownames(UMAP_cord))
UMAP_cord<-UMAP_cord[rownames(UMAP_cord)%in%colnames(Total_BM),]
colnames(UMAP_cord)<-colnames(Total_BM@reductions$umap@cell.embeddings)
Total_BM@reductions$umap@cell.embeddings<-as.matrix(UMAP_cord)

CMS<-1-getCFS(Total_BM,ident.1 = "SuPERR_cells",ident.2 = "Cite_cluster")
misclass_cells<-misclassification(Total_BM,cell_type = "SuPERR_cells",clusters = "Cite_cluster")

DimPlot(Total_BM,label = T,
        cells.highlight = rownames(misclass_cells)[misclass_cells$misclassification==TRUE],repel = T,
        pt.size = 1,label.size = 5)+theme(legend.position = "")+ggtitle(round(CMS,4))

ggsave(filename = "PBMC_misclassify_CiteFuse.pdf",width = 118.5,height = 100.1,units = "mm",device = "pdf",scale = 1.4)








### panel c------
# Left
Idents(Total_PBMC)<-Total_PBMC$S4_cluster
cluster4<-subset(Total_PBMC,subset = S4_cluster==4)
FeatureScatter(cluster4,feature1 = "CD3-TotalSeqC",pt.size = 2,feature2 = "CD56-TotalSeqC")+
  ggtitle("")+xlab("CD3-surface")+ylab("CD56-surface")+theme(legend.position = 0)
ggsave(filename = "cluster3_contamination.pdf",width = 118.5,height = 100,units = "mm",device = "pdf",scale = 1.4)
# Middle

a<-FindMarkers(Total_PBMC,ident.1 = 4,only.pos = T)
a<-a[a$p_val_adj<0.05,]
a<-a[order(a$avg_log2FC,decreasing = T),]

cluster4<-colnames(Total_PBMC)[Idents(Total_PBMC)==4]
Idents(Total_PBMC)<-Total_PBMC$SuPERR_cells
cluster4<-subset(Total_PBMC,cells = cluster4,idents = "NK")
test<-SetIdent(Total_PBMC,cells = colnames(cluster4),value = "test")
b<-FindMarkers(test,ident.1 = "test",only.pos = T)
b<-b[b$p_val_adj<0.05,]
b<-b[order(b$avg_log2FC,decreasing = T),]


input1<- list(Seurat_cluster4=rownames(a),clean_cluster4=rownames(b))
grid.newpage()
temp<-venn.diagram(input1, fill = c("#FBB4AE", "#B3CDE3"), 
                   alpha = c(0.5, 0.5), lwd =3,cex = 2,cat.fontface = 2,
                   lty =1, filename =NULL)
grid.draw(temp)

pdf(file="venn.pdf")
grid.draw(temp)
dev.off()

ggsave(filename = "cluster4_DEG.pdf",width = 118.5,height = 100,units = "mm",device = "pdf",scale = 1.4)



# Right

cluster4<-subset(Total_PBMC,subset = S4_cluster==4)
visual<-data.frame(row.names = colnames(cluster4),CD3_surface=cluster4@assays$Protein@counts["CD3-TotalSeqC",],
                   CD56_surface=cluster4@assays$Protein@counts["CD56-TotalSeqC",],
                   TRGV10=cluster4@assays$RNA@data["TRDC",])
ggplot() +  geom_point(data=visual, aes(x=CD3_surface, y=CD56_surface, colour=TRGV10),size=2)+
  scale_color_gradient(low="grey", high="red")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        panel.background = element_blank(),axis.text = element_text(size = 13),axis.title = element_text(size = 15))
ggsave(filename = "cluster3_TRGV10.pdf",width = 118.5,height = 100,units = "mm",device = "pdf",scale = 1.4)
