library(Seurat)
library(ggplot2)
library(ggvenn)
library(scDblFinder)
library(RColorBrewer)
library(scales)
library(VennDiagram)
### Functions ----
#Violin
plotfigure<-function(object=NULL,feature=NULL,x_adj=F){
  if(grepl("TotalSeqC",feature)){
    D<-data.frame(row.names  = colnames(object),
                  IGHD=object@assays$Protein@data[feature,],
                  celltype=Idents(object))
    feature_name<-gsub("-","_",feature)
    colnames(D)<-c(feature_name,"cluster")
    D$cluster<-factor(D$cluster,levels = levels(D$cluster))
  } 
  else if (feature %in% rownames(object@assays$RNA@counts)){
    D<-data.frame(row.names  = colnames(object),
                  IGHD=object@assays$RNA@data[feature,],
                  celltype=Idents(object))
    feature_name<-gsub("-","_",feature)
    colnames(D)<-c(feature_name,"cluster")
    D$cluster<-factor(D$cluster,levels = levels(D$cluster))
    
  }
  
  else {
    D<-data.frame(row.names  = colnames(object),
                  IGHD=object[[feature]],
                  celltype=Idents(object))
    feature_name<-gsub("-","_",feature)
    colnames(D)<-c(feature_name,"cluster")
    D$cluster<-factor(D$cluster,levels = levels(D$cluster))
  }
  # return(D)
  if(x_adj){
    ggplot(D, aes_string(x="cluster", y=feature_name,color="cluster",fill="cluster")) +
      geom_jitter(shape=16, position=position_jitter(0.2),size=0.5,color="black")+
      geom_violin(alpha=0.5,scale = "width",size=1)+
      stat_boxplot(geom="errorbar",width=0.1,color="red",size=0.8)+
      geom_boxplot(width=0.2,color="red",outlier.alpha = 0)+
      geom_hline(yintercept=median(D[,feature_name]),color="grey",size=1)+
      theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(),axis.line = element_line(colour = "black"),
            text = element_text(size = 15,family = "ArialMT"),axis.text.x = element_text(vjust=0.5,angle = 90,hjust = 1),plot.title = element_text(hjust = 0.5))+
      ylab("Expression level")+ggtitle(feature_name)+xlab("")
  }
  else{
    ggplot(D, aes_string(x="cluster", y=feature_name,color="cluster",fill="cluster")) +
      geom_jitter(shape=16, position=position_jitter(0.2),size=0.5,color="black")+
      geom_violin(alpha=0.5,scale = "width",size=1)+
      stat_boxplot(geom="errorbar",width=0.1,color="red",size=0.8)+
      geom_boxplot(width=0.2,color="red",outlier.alpha = 0)+
      geom_hline(yintercept=median(D[,feature_name]),color="grey",size=1)+
      theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(),axis.line = element_line(colour = "black"),
            text = element_text(size = 15,family = "ArialMT"),plot.title = element_text(hjust = 0.5))+ylab("Expression level")+ggtitle(feature_name)+xlab("")
  }
}
#Scatter 1
scatter_DB<-function(object=NULL,Dbfinder=NULL,doublet=NULL,feature1=NULL,feature2=NULL){
  SuPERR<-colnames(object)[Idents(object)==doublet]
  doublet_object<-subset(object,cells=union(setdiff(Dbfinder,PBMCs_doublets),SuPERR))
  doublet_object<-Seurat::SetIdent(doublet_object,cells=setdiff(Dbfinder,PBMCs_doublets),value="scDblFinder")
  doublet_object<-Seurat::SetIdent(doublet_object,cells=intersect(Dbfinder,SuPERR),value="shared")
  doublet_object<-Seurat::SetIdent(doublet_object,cells=setdiff(SuPERR,Dbfinder),value="SuPERR")
  finalplot<-FeatureScatter(doublet_object,feature1 = feature1,feature2 = feature2,
                            cols = list("red","green2","skyblue2"),pt.size = 1.5)+
    # ggtitle(doublet)+xlim(min(PBMC_reference_level@assays$Protein@counts[feature1,]),max(PBMC_reference_level@assays$Protein@counts[feature2,]))+ylim(min(PBMC_reference_level@assays$Protein@counts[feature2,]),max(PBMC_reference_level@assays$Protein@counts[feature2,]))
    ggtitle(doublet)+
    xlim(0,max(PBMC_reference_level[[feature1]]))+
    ylim(0,max(PBMC_reference_level[[feature2]]))
  
  return(finalplot)
}


#Scatter 2
scatter_DB_1<-function(object=NULL,Dbfinder=NULL,doublet=NULL,feature1=NULL,feature2=NULL){
  SuPERR<-colnames(object)[Idents(object)==doublet]
  doublet_object<-subset(object,cells=union(setdiff(Dbfinder,PBMCs_doublets),SuPERR))
  doublet_object<-Seurat::SetIdent(doublet_object,cells=setdiff(Dbfinder,PBMCs_doublets),value="scDblFinder")
  doublet_object<-Seurat::SetIdent(doublet_object,cells=intersect(Dbfinder,SuPERR),value="shared")
  doublet_object<-Seurat::SetIdent(doublet_object,cells=setdiff(SuPERR,Dbfinder),value="SuPERR")
  finalplot<-FeatureScatter(doublet_object,feature1 = feature1,feature2 = feature2,
                            cols = list("red","green2","skyblue2"),pt.size = 1.5)+
    ggtitle(doublet)+
    xlim(min(PBMC_reference_level@assays$Protein@counts[feature1,]),max(PBMC_reference_level@assays$Protein@counts[feature1,]))+
    # xlim(0,max(PBMC_reference_level[[feature1]]))+
    ylim(min(PBMC_reference_level@assays$Protein@counts[feature2,]),max(PBMC_reference_level@assays$Protein@counts[feature2,]))
  # ylim(0,max(PBMC_reference_level[[feature2]]))
  return(finalplot)
}
###Load gating strategy----
PBMCs_major_lineages <- readRDS("~/GitHub/SuPERR-Seq/data/Major lineages/PBMCs_major_lineages.rds")
PBMCs_doublets <- readRDS("~/GitHub/SuPERR-Seq/data/Major lineages/PBMCs_doublets.rds")
BMs_major_lineages <- readRDS("~/GitHub/SuPERR-Seq/data/Major lineages/BMs_major_lineages.rds")
BMs_doublets <- readRDS("~/GitHub/SuPERR-Seq/data/Major lineages/BMs_doublets.rds")

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
PBMC1_object[["batch"]]<-"PBMC1"
PBMC2_object[["batch"]]<-"PBMC2"
PBMC3_object[["batch"]]<-"PBMC3"

Total_PBMC<-merge(x=PBMC1_object,y = c(PBMC2_object,PBMC3_object))
rm(PBMC1,PBMC2,PBMC3,PBMC1_object,PBMC2_object,PBMC3_object)
Total_PBMC <- NormalizeData(Total_PBMC, normalization.method = "LogNormalize", scale.factor = 10000)
Total_PBMC <- FindVariableFeatures(Total_PBMC, selection.method = "vst", nfeatures = 2000)
Total_PBMC <- ScaleData(Total_PBMC, features = VariableFeatures(object = Total_PBMC),verbose = F)
Total_PBMC <- RunPCA(Total_PBMC, features = VariableFeatures(object = Total_PBMC),verbose = F)
Total_PBMC <- RunUMAP(Total_PBMC, dims = 1:10,verbose = F)
Total_PBMC[['percent.RIB']]<-PercentageFeatureSet(Total_PBMC, pattern = "^RP[SL]")

# Seurat4 UMAP PBMC
UMAP_cord<-read.csv("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/Junkai/Junkai_update/Seurat4/umap_cell_embeddings_Anchor_wnn_PBMC.csv",row.names = "X")
rownames(UMAP_cord)<-sub("PBMC","",rownames(UMAP_cord))
Total_PBMC<-subset(Total_PBMC,cells = rownames(UMAP_cord))
colnames(UMAP_cord)<-colnames(Total_PBMC@reductions$umap@cell.embeddings)
Total_PBMC@reductions$umap@cell.embeddings<-as.matrix(UMAP_cord)




### Pannel A PBMC------
doublets<-names(PBMCs_doublets)
Total_PBMC<-SetIdent(Total_PBMC,cells = doublets,value = "doublets")

Idents(Total_PBMC)<-factor(Idents(Total_PBMC),levels = c("PBMC","doublets"))

setwd("H:/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/Junkai/Junkai_update/PBMC_doublets/NEW")
# Doublets PBMC
doublets_qc<-FeatureScatter(Total_PBMC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",cells = doublets,group.by = "batch" )+
  ggtitle("")+ylim(0,6500)+xlim(0,55000)+
  geom_vline(xintercept=c(4*sd(Total_PBMC$nCount_RNA)+mean(Total_PBMC$nCount_RNA)),color=c("salmon2"),size=0.5,linetype=c("dashed"))+
  geom_hline(yintercept=c(4*sd(Total_PBMC$nFeature_RNA)+mean(Total_PBMC$nFeature_RNA)),color=c("salmon2"),size=0.5,linetype=c("dashed"))
ggsave(filename = "doublets_qc.pdf",width = 118.5,height = 75.1,units = "mm",device = "pdf",plot = doublets_qc,scale = 1.4)

#Singlets PBMC 
Real_PBMC<-subset(Total_PBMC,cells = names(PBMCs_major_lineages))
Idents(Real_PBMC)<-PBMCs_major_lineages
real_qc<-FeatureScatter(Real_PBMC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",cells = colnames(Real_PBMC))+
  ylim(0,6500)+xlim(0,55000)+ggtitle("")+
  geom_vline(xintercept=c(4*sd(Real_PBMC$nCount_RNA)+mean(Real_PBMC$nCount_RNA)),color=c("salmon2"),size=0.5,linetype=c("dashed"))+
  geom_hline(yintercept=c(4*sd(Real_PBMC$nFeature_RNA)+mean(Real_PBMC$nFeature_RNA)),color=c("salmon2"),size=0.5,linetype=c("dashed"))
ggsave(filename = "real_qc.pdf",width = 118.5,height = 75.1,units = "mm",device = "pdf",plot = real_qc,scale = 1.4)


preserved<-subset(Total_PBMC,cells = doublets,subset = nCount_RNA<c(4*sd(Total_PBMC$nCount_RNA)+mean(Total_PBMC$nCount_RNA)))
preserved<-subset(preserved,subset = nFeature_RNA<c(4*sd(Total_PBMC$nFeature_RNA)+mean(Total_PBMC$nFeature_RNA)))


### Load the BM samples----
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
Total_BM[['percent.RIB']]<-PercentageFeatureSet(Total_BM, pattern = "^RP[SL]")
Total_BM <- NormalizeData(Total_BM, normalization.method = "LogNormalize", scale.factor = 10000)
Total_BM <- FindVariableFeatures(Total_BM, selection.method = "vst", nfeatures = 2000)
Total_BM <- ScaleData(Total_BM, features = VariableFeatures(object = Total_BM),verbose = F)
Total_BM <- RunPCA(Total_BM, features = VariableFeatures(object = Total_BM),verbose = F)
Total_BM <- RunUMAP(Total_BM, dims = 1:10,verbose = F)

# Seurat4 UMAP BM
UMAP_cord<-read.csv("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/Junkai/Junkai_update/Seurat4/umap_cell_embeddings_Anchor_wnn_BM.csv",row.names = "X")
rownames(UMAP_cord)<-sub("BM","",rownames(UMAP_cord))
rownames(UMAP_cord)<-sub("_2","_1",rownames(UMAP_cord))
rownames(UMAP_cord)<-sub("_3","_2",rownames(UMAP_cord))
Total_BM<-subset(Total_BM,cells = rownames(UMAP_cord))
colnames(UMAP_cord)<-colnames(Total_BM@reductions$umap@cell.embeddings)
Total_BM@reductions$umap@cell.embeddings<-as.matrix(UMAP_cord)

### Pannel A BM------
doublets<-names(BMs_doublets)
Total_BM<-SetIdent(Total_BM,cells = doublets,value = "doublets")
Idents(Total_BM)<-factor(Idents(Total_BM),levels = c("BM","doublets"))

setwd("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/Junkai/Junkai_update/BM_doublets/NEW")
# Doublets BM
doublets_qc<-FeatureScatter(Total_BM, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",cells = doublets,group.by = "batch" )+
  ggtitle("")+ylim(0,8000)+xlim(0,70000)+
  geom_vline(xintercept=c(4*sd(Total_BM$nCount_RNA)+mean(Total_BM$nCount_RNA)),color=c("salmon2"),size=0.5,linetype=c("dashed"))+
  geom_hline(yintercept=c(4*sd(Total_BM$nFeature_RNA)+mean(Total_BM$nFeature_RNA)),color=c("salmon2"),size=0.5,linetype=c("dashed"))
ggsave(filename = "doublets_qc.pdf",width = 118.5,height = 75.1,units = "mm",device = "pdf",plot = doublets_qc,scale = 1.4)

# Singlets BM
Real_BM<-subset(Total_BM,cells = names(BMs_major_lineages))
Idents(Real_BM)<-BMs_major_lineages
real_qc<-FeatureScatter(Real_BM, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",cells = colnames(Real_BM))+
  ylim(0,8000)+xlim(0,70000)+ggtitle("")+
  geom_vline(xintercept=c(4*sd(Real_BM$nCount_RNA)+mean(Real_BM$nCount_RNA)),color=c("salmon2"),size=0.5,linetype=c("dashed"))+
  geom_hline(yintercept=c(4*sd(Real_BM$nFeature_RNA)+mean(Real_BM$nFeature_RNA)),color=c("salmon2"),size=0.5,linetype=c("dashed"))
ggsave(filename = "real_qc.pdf",width = 118.5,height = 75.1,units = "mm",device = "pdf",plot = real_qc,scale = 1.4)


### pannel B BM -----
setwd("H:/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/Junkai/Junkai_update/BM_doublets/NEW")
singlet_doublet<-c(BMs_major_lineages,BMs_doublets)
BM_doublets_real<-subset(Total_BM,cells = intersect(colnames(Total_BM),names(singlet_doublet)))

Idents(BM_doublets_real)<-factor(singlet_doublet,levels = c("B","Myeloid","HSPCs","138-PC","138+PC","B-HSPC","B-Myeloid","B-T",
                                                                     "Myeloid-NK","PC-HSPC","PC-NK","Others"))
plotfigure(BM_doublets_real,feature = "nCount_RNA",x_adj = T)
ggsave(filename = "violin_nUMI.pdf",width = 118.5,height = 75.1,units = "mm",device = "pdf",scale = 1.4)
plotfigure(BM_doublets_real,feature = "nFeature_RNA",x_adj = T)
ggsave(filename = "violin_nGene.pdf",width = 118.5,height = 75.1,units = "mm",device = "pdf",scale = 1.4)

### pannel B PBMC -----
setwd("H:/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/Junkai/Junkai_update/PBMC_doublets/NEW")
singlet_doublet<-c(PBMCs_major_lineages,PBMCs_doublets)
PBMC_doublets_real<-subset(Total_PBMC,cells = names(singlet_doublet))

Idents(PBMC_doublets_real)<-factor(singlet_doublet,levels = c("PCs","B","NK","NKT","CD4T",'CD8T',"Monocytes","Treg","B-Myeloid","B-NK",
                                                                         "B-T","CD4T-CD8T","CD4T-Myeloid","CD8T-Myeloid","NK-Myeloid","NKT-Myeloid",                                                                        "Others"))
plotfigure(PBMC_doublets_real,feature = "nCount_RNA",x_adj = T)+ylab("")+theme(legend.position = "")
ggsave(filename = "violin_nUMI.pdf",width = 118.5,height = 75.1,units = "mm",device = "pdf",scale = 1.4)
plotfigure(PBMC_doublets_real,feature = "nFeature_RNA",x_adj = T)+ylab("")+theme(legend.position = "")
ggsave(filename = "violin_nGene.pdf",width = 118.5,height = 75.1,units = "mm",device = "pdf",scale = 1.4)

### pannel c -----
setwd("H:/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/Junkai/Junkai_update/BM_doublets/NEW")


DimPlot(Total_BM,cells.highlight = names(BMs_doublets),pt.size = 0.5)+theme(legend.position = "None")
ggsave(filename = "UMAP_doublets.pdf",width = 118.5,height = 100,units = "mm",device = "pdf",scale = 1.4)


setwd("H:/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/Junkai/Junkai_update/PBMC_doublets/NEW")
DimPlot(Total_PBMC,cells.highlight = names(PBMCs_doublets),pt.size = 0.5)+theme(legend.position = "None")
ggsave(filename = "UMAP_doublets.pdf",width = 118.5,height = 100,units = "mm",device = "pdf",scale = 1.4)


### pannel d-------------------------------------------------------------------

# PBMC scdblfinder
PBMC1<-subset(Total_PBMC,subset = batch=="PBMC1")
PBMC2<-subset(Total_PBMC,subset = batch=="PBMC2")
PBMC3<-subset(Total_PBMC,subset = batch=="PBMC3")

PBMC<-list(PBMC1,PBMC2,PBMC3)
scd_doubletPBMC<-c()

for (i in 1:length(PBMC)){
  PBMC[[i]]<-NormalizeData(PBMC[[i]],scale.factor = 10000)
  PBMC[[i]] <- FindVariableFeatures(PBMC[[i]], selection.method = "vst", nfeatures = 2000,verbose = F)
  PBMC[[i]] <- ScaleData(PBMC[[i]],verbose = F)
  PBMC[[i]] <- RunPCA(PBMC[[i]],verbose = F)
  sce<-as.SingleCellExperiment(PBMC[[i]])
  sce <- scDblFinder(sce)
  # doublet2<-head(colnames(sce2)[order(sce2$scDblFinder.score,decreasing = T)],n=length(grep("_2",names(PBMCs_doublets))))
  scd_doubletPBMC<-c(scd_doubletPBMC,colnames(sce)[sce$scDblFinder.class=="doublet"])
}

Real_PBMC<-subset(Total_PBMC,cells = names(PBMCs_major_lineages))
scd_doubletPBMC<-intersect(scd_doubletPBMC,union(names(PBMCs_doublets),colnames(Real_PBMC)))


### Venn diagram
input1<- list(scDblfinder=scd_doubletPBMC,SuPERR_seq=names(PBMCs_doublets))
grid.newpage()
temp<-venn.diagram(input1, fill = c("coral2", "skyblue1"), 
                   alpha = c(0.5, 0.5), lwd =3,cex = 2,cat.fontface = 2,
                   lty =1, filename =NULL)
grid.draw(temp)

pdf(file="C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/Junkai/Junkai_update/PBMC_doublets/NEW/PBMCscdblvenn.pdf")
grid.draw(temp)
dev.off()

# BM scdblFinder
BM1<-subset(Total_BM,subset = batch=="BM1")
BM2<-subset(Total_BM,subset = batch=="BM2")


BM<-list(BM1,BM2)
scd_doubletBM<-c()

for (i in 1:length(BM)){
  BM[[i]]<-NormalizeData(BM[[i]],scale.factor = 10000)
  BM[[i]] <- FindVariableFeatures(BM[[i]], selection.method = "vst", nfeatures = 2000,verbose = F)
  BM[[i]] <- ScaleData(BM[[i]],verbose = F)
  BM[[i]] <- RunPCA(BM[[i]],verbose = F)
  sce<-as.SingleCellExperiment(BM[[i]])
  sce <- scDblFinder(sce)
  # doublet2<-head(colnames(sce2)[order(sce2$scDblFinder.score,decreasing = T)],n=length(grep("_2",names(BMs_doublets))))
  scd_doubletBM<-c(scd_doubletBM,colnames(sce)[sce$scDblFinder.class=="doublet"])
}

Real_BM<-subset(Total_BM,cells = names(BMs_major_lineages))
scd_doubletBM<-intersect(scd_doubletBM,union(names(BMs_doublets),colnames(Real_BM)))

input1<- list(scDblfinder=scd_doubletBM,SuPERR_seq=names(BMs_doublets))
grid.newpage()
temp<-venn.diagram(input1, fill = c("coral2", "skyblue1"), 
                   alpha = c(0.5, 0.5), lwd =3,cex = 2,cat.fontface = 2,
                   lty =1, filename =NULL)
grid.draw(temp)

pdf(file="C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/Junkai/Junkai_update/BM_doublets/NEW/BMscdblvenn.pdf")
grid.draw(temp)
dev.off()


### pannel E ----
# SuPERR-seq specific doubelts

b.data<-data.frame(table(PBMCs_doublets),group="SuPERR")
SuPER_fra<-ggplot(b.data, aes(fill=PBMCs_doublets,y=Freq, x=group)) + 
  geom_bar(position="fill", stat="identity",lwd = 1,color="white")+
  xlab("")+ylab("Fraction")+labs(fill="Doublets")
ggsave(filename = "PBMC_frac.pdf",width = 45,height = 110,units = "mm",device = "pdf",scale = 1.4,plot = SuPER_fra)

### pannel F-----
# Scatter plot for SuPERR-seq doublets
BnoIg_sum<-colSums(Total_PBMC@assays$RNA@counts[grep("^CD79A$|^MS4A1$|^CD19$|^VPREB3$",rownames(Total_PBMC)),])
CD4T_sum<-colSums(Total_PBMC@assays$RNA@counts[grep("^CD4$|^CD3D$|^CD3E$|^CD5$|^IL7R$",rownames(Total_PBMC)),]) #Please remove CD5 and IL7R
CD8T_sum<-colSums(Total_PBMC@assays$RNA@counts[grep("^CD8A$|^CD8B$|^GZMK$|^GZMH$|^CD3D$|^CD3E$",rownames(Total_PBMC)),])
Myeloid_sum<-colSums(Total_PBMC@assays$RNA@counts[grep("^LYZ$|^S100A8$|^S100A9$|^S100A12$|^CD68$|^CD14$|^CYBB$",rownames(Total_PBMC)),])
NK_sum<-colSums(Total_PBMC@assays$RNA@counts[grep("^GNLY$|^NKG7$|^GZMB$|^KLRD1$|^GZMA$",rownames(Total_PBMC)),])
NKT_sum<-colSums(Total_PBMC@assays$RNA@counts[grep("^CCL5$|^GNLY$|^TRDV2$|^NKG7$|^GZMH$|^KLRB1$|^CST7$",rownames(Total_PBMC)),])
Treg_sum<-colSums(Total_PBMC@assays$RNA@counts[grep("^FOXP3$|^CTLA4$|^IL2RA$|^IL32$",rownames(Total_PBMC)),])
T_sum<-colSums(Total_PBMC@assays$RNA@counts[grep("^CD3D$|^CD3E$|^CD4$|^CD8A$|^CD8B$|^CD5$|^IL7R$|^GZMK$|^GZMH$",rownames(Total_PBMC)),])
CD4TT_sum<-colSums(Total_PBMC@assays$RNA@counts[grep("^CD4$|^CD5$|^IL7R$",rownames(Total_PBMC)),])   #Please remove CD5 and IL7R
CD8TT_sum<-colSums(Total_PBMC@assays$RNA@counts[grep("^CD8A$|^CD8B$|^GZMK$|^GZMH$",rownames(Total_PBMC)),])

# Total_PBMC$B_score<-log1p(10000*B_sum/Total_PBMC$nCount_RNA)
Total_PBMC$B_score<-log1p(10000*BnoIg_sum/Total_PBMC$nCount_RNA)
Total_PBMC$CD4T_score<-log1p(10000*CD4T_sum/Total_PBMC$nCount_RNA)
Total_PBMC$CD8T_score<-log1p(10000*CD8T_sum/Total_PBMC$nCount_RNA)
Total_PBMC$Myeloid_score<-log1p(10000*Myeloid_sum/Total_PBMC$nCount_RNA)
Total_PBMC$NK_score<-log1p(10000*NK_sum/Total_PBMC$nCount_RNA)
Total_PBMC$NKT_score<-log1p(10000*NKT_sum/Total_PBMC$nCount_RNA)
Total_PBMC$Treg_score<-log1p(10000*Treg_sum/Total_PBMC$nCount_RNA)
Total_PBMC$CD4TT_score<-log1p(10000*CD4TT_sum/Total_PBMC$nCount_RNA)
Total_PBMC$CD8TT_score<-log1p(10000*CD8TT_sum/Total_PBMC$nCount_RNA)
Total_PBMC$T_score<-log1p(10000*T_sum/Total_PBMC$nCount_RNA)

doubletPBMC<-subset(Total_PBMC,cells = names(PBMCs_doublets))
Idents(doubletPBMC)<-PBMCs_doublets
Real_PBMC<-subset(Total_PBMC,cells = names(PBMCs_major_lineages))
Idents(Real_PBMC)<-PBMCs_major_lineages

p1<-scatter_DB(object=doubletPBMC,Dbfinder=scd_doubletPBMC,doublet="B-T",
               feature1='B_score',feature2="T_score")+
  geom_rect(xmin=1.8,xmax=4.1,ymin=1.5,ymax=4,alpha=0,color="grey",linetype=8)
p2<-FeatureScatter(doubletPBMC,feature1 = "B_score",feature2 = "T_score",pt.size = 1.5)+ggtitle("")+
  geom_rect(xmin=1.8,xmax=4.1,ymin=1.5,ymax=4,alpha=0,color="grey",linetype=8)
p3<-scatter_DB(object=doubletPBMC,Dbfinder=scd_doubletPBMC,doublet="CD4T-Myeloid",
               feature1='CD4T_score',feature2="Myeloid_score")+
  geom_rect(xmin=1.6,xmax=3.6,ymin=2.3,ymax=6,alpha=0,color="grey",linetype=8)
p4<-FeatureScatter(Real_PBMC,feature1 = "CD4T_score",feature2 = "Myeloid_score",pt.size = 1.5)+ggtitle("")+
  geom_rect(xmin=1.6,xmax=3.6,ymin=2.3,ymax=6,alpha=0,color="grey",linetype=8)


# p4<-FeatureScatter(PBMC_reference_level,feature1 = "CD4T_score",feature2 = "Myeloid_score",pt.size = 1.5)+ggtitle("")+
#   geom_rect(xmin=1.6,xmax=3.6,ymin=2.3,ymax=6,alpha=0,color="grey",linetype=8)

a<-ggarrange(p1,p2,p3,p4,nrow=2,ncol=2)
b<-ggarrange(p1,p2,p3,p4,nrow=1,ncol=4)
c<-ggarrange(p1,p2,p3,p4,nrow=1,ncol=4)
setwd("C:/Users/jyan399/OneDrive - Emory University/Box/Junkai Yang Ghosn Lab/updates/SuPERR-seq/Junkai/Junkai_update/PBMC_doublets/NEW")
ggsave(filename = "BnoIg_T_doublets.pdf",width = 150,height = 60,units = "mm",device = "pdf",scale = 1.4,plot = a)
ggsave(filename = "BT_CD4TMyeloid_doublets.pdf",width = 300,height = 60,units = "mm",device = "pdf",scale = 1.4,plot = b)
ggsave(filename = "BT_adt_doublets.pdf",width = 300,height = 60,units = "mm",device = "pdf",scale = 1.4,plot = c)
