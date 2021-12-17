
library(plyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(stringr)
library(Seurat)
library(DelayedArray)
library(DelayedMatrixStats)
require(gridExtra)
library(tidyr)
library("dendextend")
library(DESeq2)
#can be downloaded from https://github.com/Berlin-Hamster-Single-Cell-Consortium/Single-cell-sequencing-of-COVID-19-pathogenesis-in-golden-Hamsters-
source("./smooth_DimPlot.R")
library(akima)
library(pheatmap)
#see http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)
source("./summarySE.R")
library(cowplot)
library(SeuratObject)
library(lme4)
library(ComplexHeatmap)
library(patchwork)

##################
#read integrated object (output of merge_integrate.R)
DexAb.int <- readRDS("./DexAb_combined_integrated.rds")

#Object containing the Nr3c1 information
DexAb.int.b <- readRDS("./DexAb_combined_integrated_b.rds")

#read the annotated object
DexAb.int <- readRDS("./DexAb.int.rds")

#read in Neutrophil object
DexAb.neu <- readRDS("/Applications/stuff/DexAb/DexAb.neu.rds")


################
#supplementary code to get the Nr3c1 expression
celltype_info <- DexAb.int@meta.data %>% rownames_to_column("cell_id") %>% select("cell_id", "celltype")
DexAb.int.b@meta.data$cell_id <- rownames(DexAb.int.b@meta.data)
DexAb.int.b@meta.data <- left_join(DexAb.int.b@meta.data , celltype_info , by = c('cell_id')) %>% replace(., is.na(.), "unknown")
rownames(DexAb.int.b@meta.data) <- DexAb.int.b@meta.data$cell_id
DexAb.int.b@meta.data$orig.ident <- gsub("Untr","aaUntr",DexAb.int.b@meta.data$orig.ident)
DexAb.int.b@meta.data$treatment <- gsub("_[1-4].*","",DexAb.int.b@meta.data$orig.ident)
DexAb.int.b@meta.data$hamster <- gsub(".*_([1-4].*)","Ha\\1",DexAb.int.b@meta.data$orig.ident)
DefaultAssay(DexAb.int.b) <- "SCT"
DexAb.int.b@meta.data = cbind(DexAb.int.b@meta.data, DexAb.int.b@reductions$umap@cell.embeddings)


##############
#annotation of the object, this as already been done with the object DexAb.int.rds.gz on https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE191080
DefaultAssay(DexAb.int) <- 'RNA'
SCoV2_rawcounts <- FetchData(DexAb.int, grep("SCoV2", DexAb.int@assays$RNA@counts@Dimnames[[1]], value="TRUE"), slot="counts")
DexAb.int@meta.data$SCoV2_sum <- rowSums(SCoV2_rawcounts)
SCoV2_rawcounts$SCoV2sum <- rowSums(SCoV2_rawcounts)
DexAb.int@meta.data$SCoV2_load <- SCoV2_rawcounts$SCoV2sum/DexAb.int@meta.data$nCount_RNA*100
DexAb.int@meta.data = cbind(DexAb.int@meta.data, DexAb.int@reductions$umap@cell.embeddings)
DefaultAssay(DexAb.int) <- 'SCT'

DexAb.int@meta.data$orig.ident <- gsub("Untr","aaUntr",DexAb.int@meta.data$orig.ident)
DexAb.int@meta.data$treatment <- gsub("_[1-4].*","",DexAb.int@meta.data$orig.ident)
DexAb.int@meta.data$hamster <- gsub(".*_([1-4].*)","Ha\\1",DexAb.int@meta.data$orig.ident)

ggplot()+geom_point(data=subset(DexAb.int@meta.data, SCoV2_load<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", size=0.5, shape = 16)+
  geom_point(data=subset(DexAb.int@meta.data, SCoV2_load>0), aes(x=UMAP_1, y=UMAP_2, colour=log10(SCoV2_load)), size=0.5, shape = 16)+
  scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle("log10 SCoV2 load")
ggsave(paste("log10SCoV2_load", "pdf", sep="."), useDingbats=FALSE)



#Plot with clusters
UMAPPlot(DexAb.int, label=TRUE)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle("Clusters")
ggsave(paste("clusters", "pdf", sep="."), useDingbats=FALSE)

#To adjust cluster granularity
DefaultAssay(DexAb.int) <- 'integrated'
DexAb.int <- FindClusters(DexAb.int, resolution = 1)
DefaultAssay(DexAb.int) <- 'SCT'


#Plots with treatment and hamsters
ggplot()+
  geom_point(data=DexAb.int@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=treatment), size=0.5, shape=16, alpha=0.5)+
  ggtitle("treatment")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())
ggsave("treatment.pdf", useDingbats=FALSE)
ggplot()+
  geom_point(data=DexAb.int@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=hamster), size=0.2, shape=16, alpha=0.5)+
  ggtitle("hamster")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())
ggsave("hamsters.pdf", useDingbats=FALSE)




###############
#cell type annotation

markers <- c("Ace2", "Ackr1", "Acta2", "Adgre1", "Arg1", "Aspn", "C1qb", "Camp", "Ccr2", "Ccr5", "Cd3e", "Cd4", "Cd79b", "Cd8a", "Cldn5", "Col1a2", "Cx3cr1", "Cxcr2", "Dcn", "Ednrb", "Fcgr4", "Flt3", "Foxj1", "Gzma", "Il7r", "Irf8", "Lamp3", "Marco", "Mrc1", "Ms4a1", "Nkg7", "Plvap", "Retn", "Rtkn2", "S100a8", "Siglecf", "Tagln", "Tcf4", "Treml4", "Pdpn")
#for (gene in markers) {
#for (gene in subset(top10, cluster %in% c(19, 22, 29, 37, 40) & gene !="4930438A08Rik")$gene) {
for (gene in c("Csf1r", "Ccr1", "Ccr4", "Ccr5")) {
  df <- FetchData(DexAb.int, gene)
  gene <- str_replace_all(gene, "-", "_")
  colnames(df) <- gene
  df = cbind(df, DexAb.int@reductions$umap@cell.embeddings)
  ggplot()+geom_point(data=subset(df, eval(parse(text = gene))<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", shape=16, size=0.5)+
    geom_point(data=subset(df, eval(parse(text = gene)) > 0), aes(x=UMAP_1, y=UMAP_2, colour=eval(parse(text = gene))), shape=16, size=0.5)+
    scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    coord_fixed(ratio=1)+
    ggtitle(gene)
  ggsave(paste(gene, "pdf", sep="."), useDingbats=FALSE)
}

avg <- AverageExpression(DexAb.int, assays=c("SCT"), features = markers, return.seurat = T, slot="data") 

DoHeatmap(avg, size=5, features=markers)+scale_fill_gradient2(low="midnightblue", mid= "white", high="red", na.value = "white")+NoLegend()
ggsave("cluster_heatmap.pdf")

Idents(DexAb.int) <- DexAb.int@meta.data$seurat_clusters
cluster.markers <- FindAllMarkers(DexAb.int, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cluster.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(DexAb.int, features = top10$gene) + NoLegend()
saveRDS(cluster.markers, "/Applications/stuff/DexAb/clustermarkers.rds")
ggplot()+geom_violin(data=DexAb.int@meta.data, aes(x=seurat_clusters, y=nFeature_RNA))
ggplot()+geom_violin(data=DexAb.int@meta.data, aes(x=seurat_clusters, y=nCount_RNA))+ylim(0,20000)
View(subset(top10, cluster %in% c(19, 22, 29, 37, 40)))

#assign clusters to cell types
Idents(DexAb.int) <- DexAb.int@meta.data$seurat_clusters
DexAb.int <- RenameIdents(DexAb.int, 
                          '0'='Endothelial',
                          '1'='Neutrophils',
                          '2'='Neutrophils',
                          '3'='MonocyticMacrophages',
                          '4'='Bcells',
                          '5'='AT2',
                          '6'='AlveolarMacrophages',
                          '7'='MonocyticMacrophages',
                          '8'='Treml4+Macrophages',
                          '9'='TNKcells',
                          '10'='Neutrophils',
                          '11'='TNKcells',
                          '12'='InterstitialMacrophages',
                          '13'='InterstitialMacrophages',
                          '14'='Endothelial',
                          '15'='TNKcells',
                          '16'='Endothelial',
                          '17'='AlveolarMacrophages',
                          '18'='Myofibroblast',
                          '19'='unknown',
                          '20'='Neutrophils',
                          '21'='TNKcells',
                          '22'='unknown',
                          '23'='Endothelial',
                          '24'='AlveolarMacrophages',
                          '25'='SmoothMuscle',
                          '26'='unknown',
                          '27'='Endothelial',
                          '28'='unknown',
                          '29'='unknown',
                          '30'='TNKcells',
                          '31'='unknown',
                          '32'='pDC',
                          '33'='MyeloidDendritic',
                          '34'='AT1',
                          '35'='unknown',
                          '36'='Ciliated',
                          '37'='unknown',
                          '38'='Bcells',
                          '39'='Fibroblasts',
                          '40'='unknown',
                          '41'='unknown')
DexAb.int@meta.data$celltype <- Idents(DexAb.int)

#Order cell types
the_celltypes = c("AlveolarMacrophages",
                  "InterstitialMacrophages",
                  "MonocyticMacrophages",
                  "Treml4+Macrophages",
                  "Neutrophils",
                  "MyeloidDendritic",
                  "pDC",
                  "TNKcells",
                  "Bcells",
                  "AT1",
                  "AT2",
                  "Fibroblasts",
                  "Ciliated",
                  "Endothelial",
                  "Myofibroblast",
                  "SmoothMuscle",
                  "unknown")

#Generate color schemes
celltypecolors = c(
  "AlveolarMacrophages" = "#DFACC4",
  "InterstitialMacrophages" = "#B97C9D",
  "MonocyticMacrophages" = "#B7245C",
  "Treml4+Macrophages" = "#3E2F5B",
  "Neutrophils" = "#0081AF",
  "MyeloidDendritic" = "#4F6D7A",
  "pDC" = "#7C6A0A",
  "TNKcells" = "#368F8B",
  "Bcells" = "#62C370",
  "AT1" = "#F7C548",
  "AT2" = "#F97E44",
  "Fibroblasts" = "#B2675E",
  "Ciliated" = "#FB3640",
  "Endothelial" = "#0D3B66",
  "Myofibroblast" = "#C4A381",
  "SmoothMuscle" = "#644536",
  "unknown" = "#CAD2C5"
)
#create named vector for colors, which is later used to color the plots etc.
#this makes colors darker by treatment, to make them lighter use  r=(col2rgb(cbright)+40)[[1]] and r=ifelse(r>255,255,r) etc.
expr <- list()
for (ct in names(celltypecolors)) {
  cbright = celltypecolors[[ct]]
  r=(col2rgb(cbright)-40)[[1]]
  r=ifelse(r<0,0,r)
  g=(col2rgb(cbright)-40)[[2]]
  g=ifelse(g<0,0,g)
  b=(col2rgb(cbright)-40)[[3]]
  b=ifelse(b<0,0,b)
  cdark = rgb(r, g, b, maxColorValue = 255)
  the_scale <- scales::seq_gradient_pal(cbright, cdark, "Lab")(seq(0,1,length.out=4))
  expr[[paste(ct, "aaUntr", sep="_")]] <- the_scale[[1]]
  expr[[paste(ct, "ab", sep="_")]] <- the_scale[[2]]
  expr[[paste(ct, "dex", sep="_")]] <- the_scale[[3]]
  expr[[paste(ct, "dexab",  sep="_")]] <- the_scale[[4]]
}
the_colors = unlist(expr)

############
#UMAP with cell types
means <- DexAb.int@meta.data %>%
  group_by(celltype) %>%
  dplyr::summarise(mean_U1 = mean(UMAP_1), mean_U2 = mean(UMAP_2)) %>%
  filter(celltype != "unknown")

ggplot()+
  geom_point(data=DexAb.int@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=celltype), size=0.5, shape=16)+
  geom_text(data=subset(means, celltype!="Unclear"), aes(x=mean_U1, y=mean_U2, label=celltype), size=3)+
  coord_fixed(ratio=1)+
  ggtitle("celltype")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  scale_colour_manual(values = celltypecolors)
ggsave("celltypes.pdf", useDingbats=FALSE)


###################
#Density plots
#smooth dimplots
cond1="aaUntr"
cond2="ab"
smooth_DimPlot(subset(DexAb.int, treatment %in% c(cond1, cond2)),
               group.by='treatment',
               reduction='umap', colors.use=c('blue3','red3'), label=TRUE, repel=TRUE)
ggsave(paste("density", cond1, "vs", cond2, "pdf", sep="."), useDingbats=FALSE)
cond1="aaUntr"
cond2="dex"
smooth_DimPlot(subset(DexAb.int, treatment %in% c(cond1, cond2)),
               group.by='treatment',
               reduction='umap', colors.use=c('blue3','red3'), label=TRUE, repel=TRUE)
ggsave(paste("density", cond1, "vs", cond2, "pdf", sep="."), useDingbats=FALSE)
cond1="aaUntr"
cond2="dexab"
smooth_DimPlot(subset(DexAb.int, treatment %in% c(cond1, cond2)),
               group.by='treatment',
               reduction='umap', colors.use=c('blue3','red3'), label=TRUE, repel=TRUE)
ggsave(paste("density", cond1, "vs", cond2, "pdf", sep="."), useDingbats=FALSE)
cond1="ab"
cond2="dex"
smooth_DimPlot(subset(DexAb.int, treatment %in% c(cond1, cond2)),
               group.by='treatment',
               reduction='umap', colors.use=c('blue3','red3'), label=TRUE, repel=TRUE)
ggsave(paste("density", cond1, "vs", cond2, "pdf", sep="."), useDingbats=FALSE)
cond1="ab"
cond2="dexab"
smooth_DimPlot(subset(DexAb.int, treatment %in% c(cond1, cond2)),
               group.by='treatment',
               reduction='umap', colors.use=c('blue3','red3'), label=TRUE, repel=TRUE)
ggsave(paste("density", cond1, "vs", cond2, "pdf", sep="."), useDingbats=FALSE)
cond1="dex"
cond2="dexab"
smooth_DimPlot(subset(DexAb.int, treatment %in% c(cond1, cond2)),
               group.by='treatment',
               reduction='umap', colors.use=c('blue3','red3'), label=TRUE, repel=TRUE)
ggsave(paste("density", cond1, "vs", cond2, "pdf", sep="."), useDingbats=FALSE)


##################
#some code for auxiliary plots

#percent cell type per treatment
df1 <- cbind.data.frame(DexAb.int@meta.data$SCoV2_load,
                        DexAb.int@meta.data$celltype,
                        DexAb.int@meta.data$treatment,
                        DexAb.int@meta.data$hamster)
colnames(df1) <- c("SCoV2_load", "celltype", "treatment", "hamster")

a = df1 %>% group_by(treatment, hamster) %>% tally(name="tot") 
b = df1 %>% group_by(celltype, treatment, hamster) %>% tally(name="pos") 
c =
  left_join(a , b , by = c('treatment', 'hamster')) %>% 
  replace(., is.na(.), 0) %>%
  mutate(fraction = pos / tot)
write.table(c, "celltypefractions.tsv", sep="\t", quote=FALSE, row.names = FALSE)

tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("celltype", "treatment"))
tgc = tgc %>% mutate(celltype = forcats::fct_relevel(celltype, the_celltypes))
legendcolors <- c("aaUntr" ="gray80", "ab"="gray65", "dex"="gray50", "dexab"="gray35")
ggplot(tgc, aes(x=celltype, y=fraction, fill=paste(celltype, treatment, sep="_")))+
  geom_bar(aes(colour=treatment), position=position_dodge(), stat="identity", size=0.2)+
  geom_point(data=c, aes(x=celltype, y=fraction, fill=paste(celltype, treatment, sep="_")), position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size=0.2)+
  geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.8), size=0.25)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle(paste("celltypes per treatment", sep=" "))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
  scale_fill_manual(values=the_colors)+
  scale_color_manual(values = legendcolors)+
  guides(fill=FALSE, color = guide_legend(override.aes = list(size=5, fill="white")))

ggsave("barplot_celltypepercentage_pertimepoint.pdf", useDingbats=FALSE)


#percent virus positive by treatment
a = df1 %>% group_by(celltype, treatment, hamster) %>% tally(name="tot") 
b = df1 %>% filter(SCoV2_load>0) %>% group_by(celltype, treatment, hamster) %>% tally(name="pos") 
c =
  left_join(a , b , by = c('celltype', 'treatment', 'hamster')) %>% 
  replace(., is.na(.), 0) %>%
  mutate(fraction = pos / tot)
#write.table(c, paste("days", gene, "tsv", sep="."), sep="\t", quote=FALSE, row.names = FALSE)
tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("celltype", "treatment"))
tgc = tgc %>% mutate(celltype = forcats::fct_relevel(celltype, the_celltypes))
legendcolors <- c("aaUntr" ="gray80", "ab"="gray65", "dex"="gray50", "dexab"="gray35")
ggplot(tgc, aes(x=celltype, y=fraction, fill=paste(celltype, treatment, sep="_")))+
  geom_bar(aes(colour=treatment), position=position_dodge(), stat="identity", size=0.2)+
  geom_point(data=c, aes(x=celltype, y=fraction, fill=paste(celltype, treatment, sep="_")), position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size=0.2)+
  geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.8), size=0.25)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle(paste("virus positive per treatment", sep=" "))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
  scale_fill_manual(values=the_colors)+
  scale_color_manual(values = legendcolors)+
  guides(fill=FALSE, color = guide_legend(override.aes = list(size=5, fill="white")))
ggsave("barplot_viruspositivepercelltype_pertimepoint.pdf", useDingbats=FALSE)

#several genes in one cell type
for (the_celltype in c("Endothelial")) {
  genes_to_display=c('Mki67', 'Top2a', 'Mx2')
  
  expr <- list()
  for (the_gene in genes_to_display) {
    df1 <- cbind.data.frame(df <- FetchData(DexAb.int, the_gene),
                            DexAb.int@meta.data$celltype,
                            DexAb.int@meta.data$treatment,
                            DexAb.int@meta.data$hamster)
    colnames(df1) <- c("gene", "celltype", "treatment", "hamster")
    df1 <- subset(df1, celltype == the_celltype)
    df1$gene_name <- the_gene
    expr[[the_gene]] <- df1
  }
  df1 <- as.data.frame(do.call(rbind, expr))
  
  a = df1 %>% group_by(gene_name, treatment, hamster) %>% tally(name="tot") 
  b = df1 %>% filter(gene>0) %>% group_by(gene_name, treatment, hamster) %>% tally(name="pos") 
  c =
    left_join(a , b , by = c('gene_name', 'treatment', 'hamster')) %>% 
    replace(., is.na(.), 0) %>%
    mutate(fraction = pos / tot)
  tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("gene_name", "treatment"))
  legendcolors <- c("aaUntr" ="gray80", "ab"="gray65", "dex"="gray50", "dexab"="gray35")
  tgc = tgc %>% mutate(gene_name = forcats::fct_relevel(gene_name, genes_to_display))
  df1 = df1 %>% mutate(gene_name = forcats::fct_relevel(gene_name, genes_to_display))
  ggplot(tgc, aes(x=gene_name, y=fraction, fill=paste(the_celltype, treatment, sep="_")))+
    geom_bar(aes(colour=treatment), position=position_dodge(0.9), stat="identity", size=0.2)+
    geom_point(data=c, aes(x=gene_name, y=fraction, fill=paste(the_celltype, treatment, sep="_")), position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), size=0.2)+
    geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(0.9), size=0.2)+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    ggtitle(paste("positive for some genes per celltype per treatment", sep=" "))+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
    scale_fill_manual(values=the_colors)+
    scale_color_manual(values = legendcolors)+
    guides(fill=FALSE, color = guide_legend(override.aes = list(size=5, fill="white")))
  ggsave(paste("somegenes", the_celltype, "pdf", sep="."), useDingbats=FALSE, width=9, height = 6)
}

#percent positive for gene and boxplot for gt0
for (gene in c("Csf1", "Ccl3")) {
  celltypes_to_display=c('TNKcells', 'Bcells', 'AlveolarMacrophages', 'Endothelial', 'Myofibroblasts', 'MonocyticMacrophages', 'AT2', 'InterstitialMacrophages', 'Treml4+Macrophages', 'Neutrophils', 'SmoothMuscle', 'MyeloidDendritic')
  #celltypes_to_display=c('Endothelial')
  
  df1 <- cbind.data.frame(df <- FetchData(DexAb.int, gene),
                          DexAb.int@meta.data$celltype,
                          DexAb.int@meta.data$treatment,
                          DexAb.int@meta.data$hamster)
  colnames(df1) <- c("gene", "celltype", "treatment", "hamster")
  df1 <- subset(df1, celltype %in% celltypes_to_display)
  
  a = df1 %>% group_by(celltype, treatment, hamster) %>% tally(name="tot") 
  b = df1 %>% filter(gene>0) %>% group_by(celltype, treatment, hamster) %>% tally(name="pos") 
  c =
    left_join(a , b , by = c('celltype', 'treatment', 'hamster')) %>% 
    replace(., is.na(.), 0) %>%
    mutate(fraction = pos / tot)
  
  expr <- list()
  #p-value calculation
  for (ct in celltypes_to_display) {
    for (d in c("aaUntr", "ab", "dex", "dexab")) {
      for (dd in c("aaUntr", "ab", "dex", "dexab")) {
        df2 <- subset(df1, treatment %in% c(d, dd) & celltype == ct)
        df2$expr <- ifelse(df2$gene>0, TRUE, FALSE)
        df2$celltype <- NULL
        df2$load <- NULL
        df2$expression <- NULL
        df2$cluster <- NULL
        df2$gene <- NULL
        res <- try(p <- summary(glmer(df2$expr ~ df2$treatment + (1 | df2$hamster), family=binomial))$coefficients[2,4])
        if(inherits(res, "try-error")) {
          p <- "NA"
        }
        expr[[paste(ct, d, dd, sep="_")]] <- p
      }
    }
  }
  pvalues <- as.data.frame(do.call(rbind, expr))
  write.table(pvalues, paste("pvalues_treatments_barplots", gene, "txt", sep="."), quote=FALSE, sep="\t", col.names=FALSE)
  
  expr <- list()
  #p-value calculation
  for (ct in celltypes_to_display) {
    for (d in c("aaUntr", "ab", "dex", "dexab")) {
      for (dd in c("aaUntr", "ab", "dex", "dexab")) {
        df2 <- subset(df1, celltype == ct)
        res <- try(p <- wilcox.test(x=subset(df2, treatment == d & gene > 0)$gene, y=subset(df2, treatment == dd & gene > 0)$gene)$p.value)
        if(inherits(res, "try-error")) {
          p <- "NA"
        }
        expr[[paste(ct, d, dd, sep="_")]] <- p
      }
    }
  }
  pvalues <- as.data.frame(do.call(rbind, expr))
  write.table(pvalues, paste("pvalues_treatments_boxplots", gene, "txt", sep="."), quote=FALSE, sep="\t", col.names=FALSE)
  write.table(c, paste("days", gene, "tsv", sep="."), sep="\t", quote=FALSE, row.names = FALSE)
  tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("celltype", "treatment"))
  tgc = tgc %>% mutate(celltype = forcats::fct_relevel(celltype, the_celltypes))
  df1 = df1 %>% mutate(celltype = forcats::fct_relevel(celltype, the_celltypes))
  legendcolors <- c("aaUntr" ="gray80", "ab"="gray65", "dex"="gray50", "dexab"="gray35")
  p1 <- ggplot(tgc, aes(x=celltype, y=fraction, fill=paste(celltype, treatment, sep="_")))+
    geom_bar(aes(colour=treatment), position=position_dodge(0.9), stat="identity", size=0.2)+
    geom_point(data=c, aes(x=celltype, y=fraction, fill=paste(celltype, treatment, sep="_")), position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), size=0.2)+
    geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.9), size=0.25)+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    ggtitle(paste("positive for", gene, "per celltype per treatment", sep=" "))+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
    scale_fill_manual(values=the_colors)+
    scale_color_manual(values = legendcolors)+
    guides(fill=FALSE, color = guide_legend(override.aes = list(size=5, fill="white")))
  p2 <- ggplot()+
    geom_boxplot(data=subset(df1, gene>0), aes(x=celltype, y=gene, fill=paste(celltype, treatment, sep="_"), colour=treatment), outlier.size=0.2, position=position_dodge(.8), width = 0.7, lwd=0.25, fatten=4)+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    ggtitle(paste("expression for", gene, "in cells positive for it",  sep=" "))+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
    scale_fill_manual(values=the_colors)+
    scale_color_manual(values = legendcolors)+
    guides(fill=FALSE, color = guide_legend(override.aes = list(size=2, fill="white")))+
    ylim(0, NA)
  pdf(file=paste(gene, "percentpositive_boxplotgt0_in_celltypes", "pdf", sep="."), width=length(celltypes_to_display)+2, height=8, useDingbats=FALSE)
  grid.arrange(p1, p2, ncol=1)
  dev.off()
}

#To show variability by animal, do:
for (gene in c("Cxcl10")) {
  celltype_to_display="AlveolarMacrophages"
  
  df1 <- cbind.data.frame(df <- FetchData(DexAb.int, gene),
                          DexAb.int@meta.data$celltype,
                          DexAb.int@meta.data$treatment,
                          DexAb.int@meta.data$hamster)
  colnames(df1) <- c("gene", "celltype", "treatment", "hamster")
  df1 <- subset(df1, celltype %in% celltype_to_display)
  
  a = df1 %>% group_by(celltype, treatment, hamster) %>% tally(name="tot") 
  b = df1 %>% filter(gene>0) %>% group_by(celltype, treatment, hamster) %>% tally(name="pos") 
  c =
    left_join(a , b , by = c('celltype', 'treatment', 'hamster')) %>% 
    replace(., is.na(.), 0) %>%
    mutate(fraction = pos / tot)
  legendcolors <- c("aaUntr" ="gray80", "ab"="gray65", "dex"="gray50", "dexab"="gray35")
  p1 <- ggplot(c, aes(x=treatment, y=fraction, fill=hamster))+
    geom_bar(aes(colour=treatment), position=position_dodge(), stat="identity", size=0.4)+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    ggtitle(paste("positive for", gene, "in", celltype_to_display, "per hamster", sep=" "))+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
    scale_color_manual(values = legendcolors)
  p2 <- ggplot(subset(df1, gene>0), aes(x=treatment, y=gene, fill=hamster))+
    geom_boxplot(aes(colour=treatment), outlier.size=0.2, position=position_dodge(.8), width = 0.7, lwd=0.25, fatten=4)+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    ggtitle(paste("expression for", gene, "in", celltype_to_display, "positive for it by hamster", sep=" "))+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
    scale_color_manual(values = legendcolors)+
    ylim(0,NA)
  pdf(file=paste(gene, "by_hamster.percentpositive", "boxplotgt0", celltype_to_display, "pdf", sep="."), width=8, height=8, useDingbats=FALSE)
  grid.arrange(p1, p2, ncol=1)
  dev.off()
}


###########################
#Pseudobulk
#Threshold of 5 cells is very low, consider increasing

expr <- list()
for (cluster in unique(Idents(DexAb.int))) {
  for (sample in unique(DexAb.int@meta.data$orig.ident)) {
    cells <- Cells(DexAb.int)[(DexAb.int@meta.data$orig.ident==sample) & (Idents(DexAb.int)==cluster)]
    if (length(cells) > 5) { 
      expr[[paste0(cluster,'_',sample)]] <- rowSums(DexAb.int@assays$RNA@counts[,cells])
    }
  }
}
for (sample in unique(DexAb.int@meta.data$orig.ident)) {
  cells <- Cells(DexAb.int)[(DexAb.int@meta.data$orig.ident==sample)]
  expr[[paste0('all_',sample)]] <- rowSums(DexAb.int@assays$RNA@counts[,cells])
}

colData <- data.frame(cell.type=factor(gsub('_.*','',names(expr))),
                      donor=gsub('[^_]*_([^_]*)_([0-9])','hamster_\\2',names(expr)),
                      treatment=gsub('[^_]*_([^_]*)_([0-9])','\\1',names(expr)),
                      group=names(expr),
                      row.names=names(expr))

counts <- do.call(cbind,expr)

clusters_to_check=c('TNKcells', 'Bcells', 'AlveolarMacrophages', 'Endothelial', 'Myofibroblasts', 'MonocyticMacrophages', 'AT2', 'InterstitialMacrophages', 'Treml4+Macrophages', 'Neutrophils', 'SmoothMuscle', 'MyeloidDendritic')

res <- list()
for (cluster in c('all',clusters_to_check)) {
  take_row <- rowSums(counts) > 0
  take_col <- (colSums(counts) > 0) & (colData[,'cell.type']==cluster)
  try({
    dds <- DESeqDataSetFromMatrix(countData=counts[take_row,take_col],
                                  colData=colData[take_col,,drop=FALSE],
                                  design=~treatment)
    if (cluster!='all')
      dds <- estimateSizeFactors(dds, type='poscounts')
    dds <- DESeq(dds)
    res[[paste0(cluster,'_abvsuntr')]] <- lfcShrink(dds,
                                                    contrast=c('treatment','ab','aaUntr'),
                                                    type='normal',
                                                    format="DataFrame") %>%
      as.data.frame() %>%
      tibble::rownames_to_column('gene_name') %>%
      dplyr::mutate(cluster=cluster,contrast='abvsuntr')
    res[[paste0(cluster,'_dexvsuntr')]] <- lfcShrink(dds,
                                                    contrast=c('treatment','dex','aaUntr'),
                                                    type='normal',
                                                    format="DataFrame") %>%
      as.data.frame() %>%
      tibble::rownames_to_column('gene_name') %>%
      dplyr::mutate(cluster=cluster,contrast='dexvsuntr')
    res[[paste0(cluster,'_dexabvsuntr')]] <- lfcShrink(dds,
                                                    contrast=c('treatment','dexab','aaUntr'),
                                                    type='normal',
                                                    format="DataFrame") %>%
      as.data.frame() %>%
      tibble::rownames_to_column('gene_name') %>%
      dplyr::mutate(cluster=cluster,contrast='dexabvsuntr')
    res[[paste0(cluster,'_dexabvsdex')]] <- lfcShrink(dds,
                                                       contrast=c('treatment','dexab','dex'),
                                                       type='normal',
                                                       format="DataFrame") %>%
      as.data.frame() %>%
      tibble::rownames_to_column('gene_name') %>%
      dplyr::mutate(cluster=cluster,contrast='dexabvsdex')
    res[[paste0(cluster,'_dexabvsab')]] <- lfcShrink(dds,
                                                      contrast=c('treatment','dexab','ab'),
                                                      type='normal',
                                                      format="DataFrame") %>%
      as.data.frame() %>%
      tibble::rownames_to_column('gene_name') %>%
      dplyr::mutate(cluster=cluster,contrast='dexabvsab')
    res[[paste0(cluster,'_dexvsab')]] <- lfcShrink(dds,
                                                     contrast=c('treatment','dex','ab'),
                                                     type='normal',
                                                     format="DataFrame") %>%
      as.data.frame() %>%
      tibble::rownames_to_column('gene_name') %>%
      dplyr::mutate(cluster=cluster,contrast='dexvsab')
  })
}

pseudobulk_dexab <- do.call(rbind,res)

write.table(pseudobulk_dexab, "./pseudobulk_dexab.txt", row.names = TRUE, sep="\t", quote=FALSE)
pseudobulk_dexab <- read.table("./pseudobulk_dexab.txt", sep="\t", header=TRUE)

#Pseudobulk with separating cells with and without virus
#Threshold of 5 cells is very low




#clusters=c('TNKcells', 'Bcells', 'AlveolarMacrophages', 'Endothelial', 'Myofibroblasts', 'MonocyticMacrophages', 'AT2', 'InterstitialMacrophages', 'Treml4+Macrophages', 'Neutrophils', 'SmoothMuscle', 'MyeloidDendritic')
clusters=c( 'AlveolarMacrophages', 'MonocyticMacrophages', 'InterstitialMacrophages', 'Treml4+Macrophages', 'Neutrophils')


for (contrasts_to_use in c("abvsuntr", "dexvsuntr", "dexabvsuntr", "dexabvsdex", "dexabvsab", "dexvsab")) {
  
  #Select top genes
  #For the small version  use n=4 in slice_min and save as pseudobulk_small
  genes <- pseudobulk_dexab %>%
    dplyr::filter(cluster %in% clusters) %>%
    dplyr::filter(contrast %in% contrasts_to_use) %>%
    dplyr::filter(!grepl('SCoV2',gene_name)) %>%
    dplyr::filter(!(is.na(padj)) & (padj < .01)) %>%
    dplyr::group_by(contrast, cluster) %>%
    dplyr::slice_min(order_by=padj,n=30,with_ties=FALSE) %>%
    dplyr::pull(gene_name)
  
  #cluster values to get proper order for genes (rows) and cell types (columns)
  df <- pseudobulk_dexab %>%
    dplyr::filter(cluster %in% clusters) %>%
    dplyr::filter(contrast %in% contrasts_to_use) %>%
    dplyr::filter(gene_name %in% genes) %>%
    dplyr::mutate(group=paste0(cluster,'_',contrast))  %>%
    dplyr::select(group,gene_name,log2FoldChange) %>%
    spread(group,log2FoldChange) %>%
    tibble::column_to_rownames('gene_name')
  
  df[is.na(df)] <- 0
  hc <- hclust(dist(df))
  gene.order <- row.names(df)[order.hclust(hc)]
  hc <- hclust(dist(t(df)))
  group.order <- colnames(df)[order.hclust(hc)]
  group.order = sub("(.+)_.*", "\\1", group.order)
  
  ggplot(pseudobulk_dexab %>%
           dplyr::filter(cluster %in% clusters) %>%
           dplyr::filter(contrast %in% contrasts_to_use) %>%
           dplyr::filter(gene_name %in% genes) %>%
           dplyr::mutate(padj=ifelse(padj < 1E-20, 1E-20, padj)) %>%
           dplyr::mutate(padj=ifelse(is.na(padj), 1, padj)) %>%
           dplyr::mutate(gene=factor(gene_name,levels=gene.order)) %>%
           dplyr::mutate(cluster=factor(cluster,levels=group.order)),
         aes(y=gene,x=cluster,size=-log10(padj),color=log2FoldChange)) +
    geom_point(shape=16) +
    scale_size_continuous(name='adj. p-value',
                          breaks=c(5,10,15),
                          labels=c("1e-5","1e-10", "<1e-15"),
                          limits=c(0,20)) +
    scale_color_gradient2(low='blue3',mid='gray',high='red3',
                          limits=c(-6,6),oob=scales::squish) +
    facet_wrap(~contrast) +
    coord_cartesian(clip = 'off') +
    theme_minimal(base_size=10) +
    guides(color=guide_colorbar(barheight=4)) +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10)) +
    labs(x='',y='')
  
  #save 8x12 inches
  h=length(gene.order) / 5
  ggsave(paste("pseudobulk", contrasts_to_use, "pdf", sep="."), width=6, height = h+2, units="in", useDingbats=FALSE)
  
}


#Option 1: select genes
#When doing specific subsets, such as everything containing "Tgf", do these two lines and comment out the filter step dplyr::filter(!(is.na(padj)) & (padj < .0001) & (abs(log2FoldChange) >= 1)) %>% below
#cytokines <- as.data.frame(grep("Tgf", row.names(hamster.integrated), value = TRUE))
colnames(cytokines) <- "gene"
#contrasts_to_use <- c("abvsuntr", "dexvsuntr", "dexabvsuntr", "dexabvsdex", "dexabvsab", "dexvsab")
contrasts_to_use <- c("abvsuntr", "dexvsuntr")
contrasts_to_use <- c("abvsuntr", "dexvsuntr", "dexabvsuntr")
clusters=c( 'AlveolarMacrophages', 'MonocyticMacrophages', 'InterstitialMacrophages', 'Treml4+Macrophages', 'Neutrophils')
clusters=c('TNKcells', 'Bcells', 'AlveolarMacrophages', 'Endothelial', 'MonocyticMacrophages', 'AT2', 'InterstitialMacrophages', 'Treml4+Macrophages', 'Neutrophils', 'SmoothMuscle', 'MyeloidDendritic')

#For Fig. 3A, use:
#contrasts_to_use <- c("abvsuntr", "dexvsuntr", "dexabvsuntr")
#clusters=c('TNKcells', 'Bcells', 'AlveolarMacrophages', 'Endothelial', 'MonocyticMacrophages', 'AT2', 'InterstitialMacrophages', 'Treml4+Macrophages', 'Neutrophils', 'SmoothMuscle', 'MyeloidDendritic')
#Option 2 with
#dplyr::filter(!(is.na(padj)) & (padj < .01)) %>%
#dplyr::filter(log2FoldChange > 2) %>%
#dplyr::group_by(contrast, cluster) %>%
#dplyr::slice_min(order_by=padj,n=3,with_ties=FALSE) %>%

#For Fig. S3D
#contrasts_to_use <- c("abvsuntr", "dexvsuntr", "dexabvsuntr", "dexabvsdex", "dexabvsab", "dexvsab")
#clusters=c('TNKcells', 'Bcells', 'AlveolarMacrophages', 'Endothelial', 'MonocyticMacrophages', 'AT2', 'InterstitialMacrophages', 'Treml4+Macrophages', 'Neutrophils', 'SmoothMuscle', 'MyeloidDendritic')
#
#
#

#For Fig. 3B
clusters=c( 'AlveolarMacrophages', 'MonocyticMacrophages', 'InterstitialMacrophages', 'Treml4+Macrophages', 'Neutrophils')
contrasts_to_use <- c("abvsuntr", "dexvsuntr", "dexabvsuntr")


#These are the gene sets defined in supplementary figure 4F
genes_for_dex_infection <- c("Ifit3", "Ifit2", "Irf7", "Ube2l6", "Apobec1", "Usp18", "Ido1", "Aire", "Fcgr4", "Rtp4", "Ifi47", "ENSMAUG00000016439""Mx2", "Ddx60", "Herc6", "Slfn11", "ENSMAUG00000009900", "Nlrc5", "Oas2", "Oas1", "Isg15", "Oasl1", "ENSMAUG00000002119", "Gbp2", "Gbp7", "Gbp5", "Cmpk2", "Cfb", "Tctex1d1", "Aicda")
genes_for_dex_virus <- ("Cd274", "Ms4a7", "Ms4a15", "Plaat3", "Ctsd", "Serping1", "Plk2", "Arg1", "Slc22a2", "ENSMAUG00000020829", "Gem", "Cd83", "Tst", "Galnt6", "Gpr84", "Cxcl10", "CUNH15orf48", "Mertk", "Il1a", "Flrt3", "Sec14l4", "Renbp", "Acod1", "Acsl3", "Cd86", "Tm4sf19", "Srxn1", "Cpq", "Rgs1", "Il18bp", "Ibsp", "Spp1", "Ccl4", "Ccl3", "Ccl5", "Ccl8", "Ccl12", "Ccl2", "Fam25a", "Isg20", "Cfap161", "Dab", "Slamf9", "Slamf7", "Ctsv", "Tnfrsf9", "Zp3r", "Vgf", "Optn", "Ptprr", "ENSMAUG00000020943", "Tgm2", "Pdss1", "Nebl", "Gng11", "Asns","Upp1", "ENSMAUG00000018670", "Plau", "Il1rn", "Socs1", "Akap2", "Susd1", "Ttc28", "Ralgds", "B3galnt1", "Metazoa-SRP.33", "Btbd16", "Tmprss6", "Afp", "Gm15056")

#recode ma gene names
genes_for_dex_infection <- genes_for_dex_infection %>% recode("ENSMAUG00000004294" = "Ifit3", "ENSMAUG00000010948"="Usp18", "ENSMAUG00000008425"="Rtp4", "ENSMAUG00000013430"="Slfn11", "ENSMAUG00000015073"="Oas1", "ENSMAUG00000007782"="Oasl1", "ENSMAUG00000005083"="Cfb")
genes_for_dex_virus <- genes_for_dex_virus %>% recode("ENSMAUG00000015397" = "Dab")
the_filename="pseudobulk_infection_virus_genes"

#define file name without suffix
the_filename=("the_filename")


#To select top x genes per contrast, do:
#if only upregulated genes do log2FoldChange > 0 filtering
genes <- pseudobulk_dexab %>%
  dplyr::filter(cluster %in% clusters) %>%
  dplyr::filter(contrast %in% contrasts_to_use) %>%
  dplyr::filter(!grepl('SCoV2',gene_name)) %>%
  dplyr::filter(!(is.na(padj)) & (padj < .01)) %>%
  #dplyr::filter(log2FoldChange > 2) %>%
  dplyr::group_by(contrast, cluster) %>%
  dplyr::slice_min(order_by=padj,n=10,with_ties=FALSE) %>%
  dplyr::pull(gene_name)
the_filename="pseudobulk_top3_upregulated_facet"

#Otherwise define some genes, e.g.
genes <- c("Cd274", "Ccl5","Spp1", "Cxcl11","Tnfsf10", "Ccl2", "Ccl8", "Ccl12", "Cxcl10", "Ccl3", "Ccl4", "Mx2", "Ifit2", "Isg15", "Siglec1", "Mx1", "Sod2", "Irf7", "Ddx58")

#clustering
df <- pseudobulk_dexab %>%
  dplyr::filter(cluster %in% clusters) %>%
  dplyr::filter(contrast %in% contrasts_to_use) %>%
  dplyr::filter(!(is.na(padj)) & (padj < .01)) %>%
  dplyr::filter(gene_name %in% genes) %>%
  dplyr::mutate(group=paste0(cluster,'_',contrast))  %>%
  dplyr::select(group,gene_name,log2FoldChange) %>%
  spread(group,log2FoldChange) %>%
  tibble::column_to_rownames('gene_name')
df[is.na(df)] <- 0
hc <- hclust(dist(df))
gene.order <- row.names(df)[order.hclust(hc)]
hc <- hclust(dist(t(df)))
group.order <- colnames(df)[order.hclust(hc)]
group.order = sub("(.+)_.*", "\\1", group.order)
genes <- row.names(df)

#when clustering two categories separately (as in Figure 4B), do
df <- pseudobulk_dexab %>%
  dplyr::filter(cluster %in% clusters) %>%
  dplyr::filter(contrast %in% contrasts_to_use) %>%
  dplyr::filter(!(is.na(padj)) & (padj < .1)) %>%
  dplyr::filter(gene_name %in% genes_for_dex_virus) %>%
  dplyr::mutate(group=paste0(cluster,'_',contrast))  %>%
  dplyr::select(group,gene_name,log2FoldChange) %>%
  spread(group,log2FoldChange) %>%
  tibble::column_to_rownames('gene_name')
df[is.na(df)] <- 0
hc <- hclust(dist(df))
gene.order1 <- row.names(df)[order.hclust(hc)]
hc <- hclust(dist(t(df)))
group.order <- colnames(df)[order.hclust(hc)]
group.order = sub("(.+)_.*", "\\1", group.order)
genes1 <- row.names(df)
df <- pseudobulk_dexab %>%
  dplyr::filter(cluster %in% clusters) %>%
  dplyr::filter(contrast %in% contrasts_to_use) %>%
  dplyr::filter(!(is.na(padj)) & (padj < .1)) %>%
  dplyr::filter(gene_name %in% genes_for_dex_infection) %>%
  dplyr::mutate(group=paste0(cluster,'_',contrast))  %>%
  dplyr::select(group,gene_name,log2FoldChange) %>%
  spread(group,log2FoldChange) %>%
  tibble::column_to_rownames('gene_name')
df[is.na(df)] <- 0
hc <- hclust(dist(df))
gene.order2 <- row.names(df)[order.hclust(hc)]
hc <- hclust(dist(t(df)))
group.order <- colnames(df)[order.hclust(hc)]
group.order = sub("(.+)_.*", "\\1", group.order)
genes2 <- row.names(df)
gene.order <- c(gene.order1, gene.order2)
genes <- c(genes1, genes2)

#now plot
ggplot(pseudobulk_dexab %>%
         dplyr::filter(cluster %in% clusters) %>%
         dplyr::filter(contrast %in% contrasts_to_use) %>%
         dplyr::filter(gene_name %in% genes) %>%
         dplyr::mutate(padj=ifelse(padj < 1E-20, 1E-20, padj)) %>%
         dplyr::mutate(padj=ifelse(is.na(padj), 1, padj)) %>%
         dplyr::mutate(gene=factor(gene_name,levels=gene.order)) %>%
         dplyr::mutate(cluster=factor(cluster,levels=the_celltypes)),
       aes(y=gene,x=cluster,size=-log10(padj),color=log2FoldChange)) +
  geom_point(shape=16) +
  scale_size_continuous(name='adj. p-value',
                        breaks=c(2,4,6),
                        labels=c("1e-2","1e-4", "<1e-6"),
                        limits=c(0,NA)) +
  scale_color_gradient2(low='blue3',mid='gray',high='red3',
                        limits=c(-3,3),oob=scales::squish) +
  facet_wrap(~contrast, nrow=1) +
  coord_cartesian(clip = 'off') +
  theme_minimal(base_size=10) +
  guides(color=guide_colorbar(barheight=4)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10)) +
  labs(x='',y='')

#save file
ggsave(paste(the_filename, "pdf", sep="."), width=length(clusters)+1, height = length(genes)/6+2, units="in", useDingbats=FALSE)




#######################
#Neutrophils cells subclustering

#read in DexAb.neu from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE191080 or define as following

DexAb.neu <- subset(DexAb.int, subset = (celltype == "Neutrophils"))
#Remove these unclear outlier cells
DexAb.neu <- subset(DexAb.neu, UMAP_2 < -4)


DefaultAssay(DexAb.neu) <- "integrated"                  
DexAb.neu <- RunPCA(DexAb.neu, verbose = FALSE)
DexAb.neu <- RunUMAP(DexAb.neu, dims = 1:30)
DexAb.neu <- FindNeighbors(DexAb.neu, dims = 1:30)
DexAb.neu <- FindClusters(DexAb.neu, resolution = 0.9)
DexAb.neu@meta.data$UMAP_1 <- NULL
DexAb.neu@meta.data$UMAP_2 <- NULL
UMAPPlot(DexAb.neu)
DexAb.neu@meta.data = cbind(DexAb.neu@meta.data, DexAb.neu@reductions$umap@cell.embeddings)
DefaultAssay(DexAb.neu) <- "SCT"   



#Plots for Figure 5 (clusters, viral load, relative density)
#Adjust the original color scale to make the viral load a bit brighter
p1 <- ggplot()+geom_point(data=DexAb.neu@meta.data, aes(x=UMAP_1, y=UMAP_2, colour=seurat_clusters), size=1, shape = 16)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle("clusters")
p2 <- ggplot()+geom_point(data=subset(DexAb.neu@meta.data, SCoV2_load<=0), aes(x=UMAP_1, y=UMAP_2), colour="grey90", size=1, shape = 16)+
  geom_point(data=subset(DexAb.neu@meta.data, SCoV2_load>0), aes(x=UMAP_1, y=UMAP_2, colour=log10(SCoV2_load)), size=1, shape = 16)+
  #scale_colour_gradientn(colours=c("#00004C", "#0C005C", "#1D006B", "#31007A", "#49008A", "#660099", "#8700A8", "#AB00B8", "#C700BA", "#D600AB", "#E60099", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
  scale_colour_gradientn(colours=c("#87a5ff", "#7075ff", "#7a69ff", "#7f54ff", "#9c38ff", "#a742ff", "#bf00ed", "#df07f0", "#fc03ec", "#fa00c8", "#ff00aa", "#F50083", "#FF0569", "#FD1754", "#FF2441", "#FF3333"), guide="colourbar")+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle("log10 SCoV2 load")
cond1="aaUntr"
cond2="dex"
p3 <- smooth_DimPlot(subset(DexAb.neu, treatment %in% c(cond1, cond2)),
                     group.by='treatment',
                     reduction='umap', colors.use=c('blue3','red3'), label=FALSE, repel=TRUE, pt.size=1.25)
p1 + p2 + p3

ggsave(paste("neu_clusters_log10SCoV2load_density", "pdf", sep="."), useDingbats=FALSE, width = 24, height=8)

#For Figure S5, do the other density plots individually
cond1="aaUntr"
cond2="ab"
p1 <- smooth_DimPlot(subset(DexAb.neu, treatment %in% c(cond1, cond2)),
               group.by='treatment',
               reduction='umap', colors.use=c('blue3','red3'),  repel=TRUE, pt.size=1.25)
cond1="aaUntr"
cond2="dexab"
p2 <- smooth_DimPlot(subset(DexAb.neu, treatment %in% c(cond1, cond2)),
               group.by='treatment',
               reduction='umap', colors.use=c('blue3','red3'),  repel=TRUE, pt.size=1.25)
p1 + p2

ggsave(paste("neu_supp_densities", "pdf", sep="."), useDingbats=FALSE, width = 16, height=8)



#cluster composition by treatment as ratio barplot
expr <- list()
for (the_treatment in c("aaUntr", "ab", "dex", "dexab")) {
  
  
  df1 <- cbind.data.frame(DexAb.neu@meta.data$seurat_clusters,
                          DexAb.neu@meta.data$treatment,
                          DexAb.neu@meta.data$hamster)
  colnames(df1) <- c("clusters", "treatment", "hamster")
  for (the_cluster in unique(DexAb.neu@meta.data$seurat_clusters)) {
    a <- dim(subset(DexAb.neu@meta.data, seurat_clusters==the_cluster))[[1]]
    b <- dim(subset(DexAb.neu@meta.data, seurat_clusters==the_cluster & treatment == the_treatment))[[1]]
    expr[[paste(the_treatment, the_cluster, sep="_")]] <- b/a
  }
}
positives <- as.data.frame(do.call(rbind,expr))
positives$name <- rownames(positives)
positives <- positives %>% separate("name", into=c("treatment", "cluster"),  sep="_")
positives$V1 <- as.numeric(positives$V1)
positives$cluster <- as.integer(positives$cluster)

df <- cbind.data.frame(subset(positives, treatment=="aaUntr")$cluster, subset(positives, treatment=="aaUntr")$V1, subset(positives, treatment=="ab")$V1, subset(positives, treatment=="dex")$V1, subset(positives, treatment=="dexab")$V1)
colnames(df) <- c("cluster", "aaUntr", "ab", "dex", "dexab")
df <- df[with(df, order(cluster)),]
df$ab <- log2(df$ab/df$aaUntr)
df$dex <- log2(df$dex/df$aaUntr)
df$dexab <- log2(df$dexab/df$aaUntr)
df$aaUntr <- NULL
df <- reshape2::melt(df, id.vars = c("cluster"))

ggplot()+
  geom_bar(data=df, aes(x=as.factor(cluster), y=value, fill=variable), position=position_dodge(), stat="identity", size=0)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle(paste("celltypes per time point", sep=" "))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))
ggsave("neu_clusters_ratios.pdf", width=5, height = 2, units="in", useDingbats=FALSE)


#virus positive per cluster and treatment
a = df1 %>% group_by(clusters, treatment, hamster) %>% tally(name="tot") 
b = df1 %>% filter(SCoV2_load>0) %>% group_by(clusters, treatment, hamster) %>% tally(name="pos") 
c =
  left_join(a , b , by = c('clusters', 'treatment', 'hamster')) %>% 
  replace(., is.na(.), 0) %>%
  mutate(fraction = pos / tot)
#write.table(c, paste("days", gene, "tsv", sep="."), sep="\t", quote=FALSE, row.names = FALSE)
tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("clusters", "treatment"))
legendcolors <- c("aaUntr" ="gray80", "ab"="gray65", "dex"="gray50", "dexab"="gray35")
ggplot(tgc, aes(x=clusters, y=fraction, fill=paste(clusters, treatment, sep="_")))+
  geom_bar(aes(colour=treatment), position=position_dodge(), stat="identity", size=0.2)+
  geom_point(data=c, aes(x=clusters, y=fraction, fill=paste(clusters, treatment, sep="_")), position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size=0.2)+
  geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.8), size=0.25)+
  theme_bw()+
  theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
  ggtitle(paste("virus positive per treatment", sep=" "))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
  scale_color_manual(values = legendcolors)+
  guides(fill=FALSE, color = guide_legend(override.aes = list(size=5, fill="white")))
ggsave("neu_barplot_viruspositivepercluster_pertreatment.pdf", useDingbats=FALSE)

#genes in clusters
for (gene in c("Cxcl10", "Mx2", "Tnfsf10")) {

  df1 <- cbind.data.frame(df <- FetchData(DexAb.neu, gene),
                          DexAb.neu@meta.data$seurat_clusters,
                          DexAb.neu@meta.data$treatment,
                          DexAb.neu@meta.data$hamster)
  colnames(df1) <- c("gene", "clusters", "treatment", "hamster")

  a = df1 %>% group_by(clusters, treatment, hamster) %>% tally(name="tot") 
  b = df1 %>% filter(gene>0) %>% group_by(clusters, treatment, hamster) %>% tally(name="pos") 
  c =
    left_join(a , b , by = c('clusters', 'treatment', 'hamster')) %>% 
    replace(., is.na(.), 0) %>%
    mutate(fraction = pos / tot)
  tgc <- summarySE(as.data.frame(c), measurevar="fraction", groupvars=c("clusters", "treatment"))
  legendcolors <- c("aaUntr" ="gray80", "ab"="gray65", "dex"="gray50", "dexab"="gray35")
  p1 <- ggplot(tgc, aes(x=clusters, y=fraction, fill=paste(clusters, treatment, sep="_")))+
    geom_bar(aes(colour=treatment), position=position_dodge(), stat="identity", size=0.2)+
    geom_point(data=c, aes(x=clusters, y=fraction, fill=paste(clusters, treatment, sep="_")), position=position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8), size=0.2)+
    geom_errorbar(aes(ymin=fraction-sd, ymax=fraction+sd), width=.2, position=position_dodge(.8), size=0.25)+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    ggtitle(paste("positive for", gene, "per clusters per treatment", sep=" "))+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
    scale_color_manual(values = legendcolors)+
    guides(fill=FALSE, color = guide_legend(override.aes = list(size=5, fill="white")))
  p2 <- ggplot()+
    geom_boxplot(data=subset(df1, gene>0), aes(x=clusters, y=gene, fill=paste(clusters, treatment, sep="_"), colour=treatment), outlier.size=0.2, position=position_dodge(.8), width = 0.7, lwd=0.25, fatten=4)+
    theme_bw()+
    theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
    ggtitle(paste("expression for", gene, "in cells positive for it",  sep=" "))+
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))+
    scale_color_manual(values = legendcolors)+
    guides(fill=FALSE, color = guide_legend(override.aes = list(size=2, fill="white")))+
    ylim(0, NA)
  pdf(file=paste(gene, "neu_percentpositive_boxplotgt0_in_clusters", "pdf", sep="."), width=16, height=8, useDingbats=FALSE)
  grid.arrange(p1, p2, ncol=1)
  dev.off()
}

#For Fig. 5 and S5D (width and midpoint need to be adjusted dependent on the genes)
#
expr <- list()
for (gene in c("Il1r2", "Isg20")) {
  
  
  df1 <- cbind.data.frame(df <- FetchData(DexAb.neu, gene),
                          DexAb.neu@meta.data$seurat_clusters,
                          DexAb.neu@meta.data$treatment,
                          DexAb.neu@meta.data$hamster)
  colnames(df1) <- c("gene", "clusters", "treatment", "hamster")
  for (the_cluster in unique(DexAb.neu@meta.data$seurat_clusters)) {
    a <- dim(subset(df1, clusters==the_cluster & treatment == "aaUntr"))[[1]]
    b <- dim(subset(df1, clusters==the_cluster & treatment == "aaUntr" & gene >0))[[1]]
    c <- mean(subset(df1, clusters==the_cluster & treatment == "aaUntr" & gene >0)$gene)
    expr[[paste(gene, the_cluster, sep="_")]] <- paste(b/a, c, sep="_")
  }
}
positives <- as.data.frame(do.call(rbind,expr))
positives$name <- rownames(positives)
positives <- positives %>% separate("name", into=c("gene", "cluster"),  sep="_") %>% separate("V1", into=c("pos", "expr"),  sep="_")
positives$pos <- as.numeric(positives$pos)
positives$expr <- as.numeric(positives$expr)
positives$cluster <- as.integer(positives$cluster)

ggplot(positives,
       aes(x=gene,y=as.factor(cluster),size=pos,color=expr)) +
  geom_point(shape=16) +
  scale_size_continuous(name='fraction positive',
                        breaks=c(0.25,0.5,0.75),
                        labels=c("0.25","0.5", "0.75"),
                        limits=c(0,1)) +
  scale_color_gradient2(low='blue3',mid='gray',high='red3', midpoint=1.2) +
  coord_cartesian(clip = 'off') +
  theme_minimal(base_size=10) +
  guides(color=guide_colorbar(barheight=4)) +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5, size=10), axis.text.y=element_text(size=10)) +
  labs(x='',y='')

ggsave("neu_expr_in_clusters.pdf", width=2.3, height = 5, units="in", useDingbats=FALSE)


