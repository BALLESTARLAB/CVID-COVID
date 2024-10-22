# NK RECLUSTERING COHORT 1 ---------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)

library(tidytext)
library(Cairo)
library(SummarizedExperiment)
library(monocle3)
library(SeuratWrappers)
library(harmony)
library(ggalluvial)
library(tidyr)

load("../subset_nk_coh1.RData")

DefaultAssay(nk_1) <- "RNA"
nk.list <- SplitObject(nk_1, split.by = "sample_names")
nk.feature <- SelectIntegrationFeatures(object.list = nk.list)
nk.anchors <- FindIntegrationAnchors(object.list = nk.list,
                                        anchor.features = nk.feature, reference = 1)
nk_1 <- IntegrateData(anchorset = nk.anchors)

save(nk_1, file = "nk_cohort1_befscale.RData")

nk_1 <- ScaleData(nk_1)
nk_1 <- RunPCA(nk_1)
nk_1 <- RunUMAP(nk_1, reduction = "pca", dims = 1:30)
nk_1 <- FindNeighbors(nk_1, reduction = "pca", dims = 1:30)
nk_1 <- FindClusters(nk_1, resolution = 2, graph.name = "integrated_snn")

save(nk_1, file = "nk_cohort1_samplenames.RData")

# Annotation

Idents(nk_1) <- "seurat_clusters"
current.cluster <- levels(nk_1)
current.cluster


new.cluster <- c("HLAhigh_CD56dim", #cl 0
                 "Inflamed_CD56dim", #cl 1
                 "Cytokine_CD56dim", #cl 2
                 "SAP_CD56dim", #cl 3
                 "CD56dim", #cl 4 
                 "HLAhigh_CD56dim", #cl 5
                 "CD56bright", #cl 6
                 "CD56dim", #cl 7 
                 "CD56bright" #cl 8
)

names(x = new.cluster) <- levels(x = nk_1)
nk_1 <- RenameIdents(object = nk_1, new.cluster)
nk_1$general_cluster <- Idents(object = nk_1)
table(nk_1$general_cluster)
Idents(nk_1) <- "general_cluster"
nk_1@active.ident <- factor(nk_1@active.ident, 
                            levels = rev(c("CD56bright",
                                           "Inflamed_CD56dim",
                                           "Cytokine_CD56dim", 
                                           "HLAhigh_CD56dim", 
                                           "CD56dim", "SAP_CD56dim"
                            )))

#Supp Figure 5b
DotPlot(nk_1, features = c("NCAM1", "SELL", "GZMK", "TCF7", "XCL1", #CD56 bright
                           "IFI6", "ISG15", "XAF1", "MX1", #Inflamed_CD56dim
                           "CCL4", "CCL4L2", "CCL3", "IFNG", #Cytokine_CD56dim
                           "HLA-DPA1", "HLA-DRA", "HLA-DRB1", "HLA-DPB1",  #HLAhigh_CD56dim
                           "PRF1", "GZMB", "FCGR3A", #CD56 dim
                           "EIF3G", "MATK", "SH2D1A" #SAP CD56dim
), scale.max = 100, scale.min = 0, dot.scale = 6, scale.by = "size")  +
  scale_color_gradientn(colors = c("skyblue3", "white", "#770000")) +
  theme_pubr(legend = "right") + border()+ rotate_x_text(45)

#Figure 5a
COL <- c("CD56dim" = "#A6D854",
         "Inflamed_CD56dim" = "#FC8D62",
         "CD56bright" = "#FFD92F",
         "Cytokine_CD56dim" = "#8DA0CB",
         "SAP_CD56dim" = "#E78AC3",
         "HLAhigh_CD56dim" = "#66C2A5")

DimPlot(nk_1, label = F, group.by = "general_cluster", pt.size = 0.5, cols = COL) +
  labs(title = "NK (n = 19,243)")

save(nk_1, file = "newdata/integration_samples/cohort1/nuevonk/nk_annot2.RData")

# Figure 5b
nk_markers <- FindAllMarkers(nk_1, assay = "RNA", logfc.threshold = 0.25, min.pct = 0.3)

top.genes.all <- data.frame()
nk_1@active.ident <- factor(nk_1@active.ident, 
                            levels = rev(c("CD56bright",
                                           "Inflamed_CD56dim",
                                           "Cytokine_CD56dim", 
                                           "HLAhigh_CD56dim", 
                                           "CD56dim", "SAP_CD56dim")))

nk.clusters <- as.vector(as.data.frame(table(nk_1$general_cluster))$Var1)

for (i in nk.clusters) {
  top.genes <- nk_markers %>% 
    mutate(sign = ifelse((p_val_adj < 0.05 & cluster == i & pct.1 > 0.3 & avg_log2FC > 0), paste("sign", i, sep = "_"), "others")) %>% 
    dplyr::filter(sign != "others")
  top.genes.all <- rbind(top.genes.all, top.genes)
}

gene.markers.nk <- (table(top.genes.all$gene)) %>% as.data.frame() %>% 
  dplyr::filter(Freq == 1) %>% 
  dplyr::filter(!grepl("^MT-", Var1))

cluster <- c("CD56bright",
             "Inflamed_CD56dim",
             "Cytokine_CD56dim", 
             "HLAhigh_CD56dim", 
             "CD56dim", "SAP_CD56dim")

top.genes.common <- top.genes.all %>% 
  dplyr::filter(gene %in% gene.markers.nk$Var1) %>% 
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) 

DefaultAssay(nk_1) <- "RNA"
DotPlot(nk_1, features = top.genes.common$gene, scale.max = 100, scale.min = 0, dot.scale = 6, scale.by = "size")  +
  scale_color_gradientn(colors = c("skyblue3", "white", "#770000")) +
  theme_pubr(legend = "right") + border()+ rotate_x_text(45) 

# Supp Figure 5a

nk_1@meta.data <- nk_1@meta.data %>% 
  mutate(PID = ifelse(StudyName == "Javi", "CVID", "non-CVID") ) %>% 
  select(-StudyName)
table(nk_1$PID)

set.seed(123)

#CVID
sample1_rna <- subset(nk_1, subset = PID == "CVID")
metadata1 <- as.data.frame(sample1_rna@meta.data)
nuevo_m1 <- row.names(sample_n(metadata1, 2677))
cvid_nk_subsample <- subset(nk_1, cells = nuevo_m1)

noncvid_nk <- subset(nk_1, subset = PID == "non-CVID") 
#non-CVID mild
sample2_rna <- subset(noncvid_nk, subset = CovidSeverity == c("Control", "Mild"))
metadata2 <- as.data.frame(sample2_rna@meta.data)
nuevo_m2 <- row.names(sample_n(metadata2, 2677))
noncvid_mild_nk_subsample <- subset(nk_1, cells = nuevo_m2)

#non-CVID severe
sample3_rna <- subset(noncvid_nk, subset = CovidSeverity == c("Control", "Severe"))
metadata3 <- as.data.frame(sample3_rna@meta.data)
nuevo_m3 <- row.names(sample_n(metadata3, 2677))
noncvid_severe_nk_subsample <- subset(nk_1, cells = nuevo_m3)

p1 <- DimPlot(cvid_nk_subsample, group.by = "general_cluster", pt.size = 0.75, cols = pal(6), raster = F) + labs(title = "CVID") + NoLegend()
p2 <- DimPlot(noncvid_mild_nk_subsample, group.by = "general_cluster", pt.size = 0.75, cols = pal(6), raster = F) + labs(title = "non-CVID mild") + NoLegend()
p3 <- DimPlot(noncvid_severe_nk_subsample, group.by = "general_cluster", pt.size = 0.75, cols = pal(6), raster = F) + labs(title = "non-CVID severe")

p1+p2+p3

# NK RECLUSTERING COHORT 2 -------------
load("~/nk_coh2_annot.RData")

nk_2 <- FindClusters(nk_2, resolution = 0.5, graph.name = "integrated_snn") 

Idents(nk_2) <- "seurat_clusters"
current.cluster <- levels(nk_2)
current.cluster

new.cluster <- c("Cytokine_CD56dim", #cl 0
                 "Inflamed_CD56dim", #cl 1
                 "CD56dim", #cl 2
                 "SAP_CD56dim", #cl 3
                 "CD56bright", #cl 4 
                 "Cytokine_CD56dim", #cl 5
                 "HLAhigh_CD56dim", #cl 6
                 "HLAhigh_CD56dim",  #cl 7
                 "CD56bright", #cl 8
                 "NA" #cl9
)

names(x = new.cluster) <- levels(x = nk_2)
nk_2 <- RenameIdents(object = nk_2, new.cluster)
nk_2$general_cluster <- Idents(object = nk_2)
table(nk_2$general_cluster)
Idents(nk_2) <- "general_cluster"

nk_2@active.ident <- factor(nk_2@active.ident, 
                            levels = rev(c("CD56bright",
                                           "Inflamed_CD56dim",
                                           "Cytokine_CD56dim", 
                                           "HLAhigh_CD56dim", 
                                           "CD56dim", "SAP_CD56dim", "NA")))

# Supp Fig 5c 
new.colors <- rev(c( "#A6D854",  #"CD56dim"
                     "#66C2A5",  #"HLAhigh_CD56dim", 
                     "#FFD92F", #"CD56bright", 
                     "#E78AC3",  #"SAP_CD56dim",
                     "#8DA0CB", #"Cytokine_CD56dim", 
                     "#FC8D62", #"Inflamed_CD56dim",
                     "gray" #NA
))

DimPlot(nk_2, label = F, group.by = "general_cluster", pt.size = 0.5, order = c(
  "CD56dim","HLAhigh_CD56dim", "CD56bright", "SAP_CD56dim",                                                               
  "Cytokine_CD56dim", "Inflamed_CD56dim"), cols = new.colors) +
  labs(title = "NK (n = 25,999)")

# Supp Fig 5d 
DotPlot(subset(nk_2, subset = general_cluster == "NA", invert = T), features = c("GZMK", "IL7R", "COTL1", "CAPG", "SPTSSB",
                                                                                 "IGFBP7", "AKR1C3", "CX3CR1", "LAIR2", "LINC00299", 
                                                                                 "KLRB1", "JUN", "AC007952.4", "CCL3", "ICAM2",
                                                                                 "CD3E", "IL32", "TRG-AS1", "HLA-DRB5", "CLEC2D",
                                                                                 "MALAT1", "RNF213", "SYNE1", "GK5", "DIP2A",
                                                                                 "EIF3G", "MATK", "SH2D1A"),
        scale.max = 100, scale.min = 0, dot.scale = 6, scale.by = "size")  +
  scale_color_gradientn(colors = c("skyblue3", "white", "#770000")) +
  theme_pubr(legend = "right") + border()+ rotate_x_text(45)

# Pseudotime analysis --------
set.seed(578)
cds <- as.cell_data_set(nk_1);
fData(cds)$gene_short_name <- rownames(fData(cds))

reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)
cds@clusters$UMAP$partitions <- reacreate.partition

# Assign the cluster info 
list_cluster <- nk_1@active.ident
cds@clusters$UMAP$clusters <- list_cluster
cds@int_colData@listData$reducedDims$UMAP <- nk_1@reductions$umap@cell.embeddings

cds <- learn_graph(cds, use_partition = FALSE)

plot_cells(cds,
           label_groups_by_cluster = T,
           label_branch_points = T,
           label_roots = T,
           label_leaves = T,
           group_label_size = 5, label_principal_points = TRUE)

cds <- order_cells(cds, reduction_method = 'UMAP')

#Fig 5e
plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE, 
           trajectory_graph_segment_size = 0.75, cell_size = 0.5,
           trajectory_graph_color = "black", alpha = 5) +
  scale_color_gradientn(colors = colorRampPalette(c("#d9ed92", "#b5e48c", "#99d98c",
                                                    "#76c893", "#52b69a", "#34a0a4","#168aad",
                                                    "#1a759f", "#1e6091", "#184e77"))(250))

pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo, 
       aes(monocle3_pseudotime, reorder(general_cluster, monocle3_pseudotime, median), fill = general_cluster), color = general_cluster) +
  geom_boxplot(size = 0.1, outlier.size = 0.1) +
  theme_pubr(border = T, legend = "right") +
  scale_fill_manual(values = c("CD56dim" = "#A6D854",
                               "Inflamed_CD56dim" = "#FC8D62",
                               "CD56bright" = "#FFD92F",
                               "Cytokine_CD56dim" = "#8DA0CB",
                               "SAP_CD56dim" = "#E78AC3",
                               "HLAhigh_CD56dim" = "#66C2A5") )

#DEGs pseudotime in seurat
deg <- readxl::read_xlsx(path = "psuedotime/nk_pseudotime_all.xlsx")

deg %>% 
  arrange(q_value) %>% 
  filter(status == 'OK') -> deg_ok

deg_ok %>% 
  arrange("morans_I") -> ok_arrange

ok_gene <- ok_arrange %>% 
  dplyr::filter(morans_I > 0) %>% 
  pull(gene_short_name)

ok_gene_neg <- ok_arrange %>% 
  dplyr::filter(morans_I < 0) %>% 
  pull(gene_short_name)

nk_1$pseudotime <- pseudotime(cds)
Idents(nk_1) <- nk_1$general_cluster
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

# Fig 5f
data("dorothea_hs", package = "dorothea")
regulons = dorothea_hs %>%
  filter(confidence %in% c("A", "B")) %>% 
  mutate(direction = ifelse(mor == "-1", "DOWN", "UP"))

regulon.n <- regulons %>% 
  filter(tf %in% c("IRF9")) %>% 
  dplyr::filter(direction == "UP")
regulon.n <- regulons %>% 
  filter(tf %in% c("STAT1")) %>% 
  dplyr::filter(direction == "UP")
regulon.n <- regulons %>% 
  filter(tf %in% c("STAT2")) %>% 
  dplyr::filter(direction == "UP")
nk_1 <- AddModuleScore(nk_1, features = list(regulon.n$target), name = "n")

plot.pseudo <- nk_1@meta.data %>%
  mutate(PID = ifelse(StudyName == "Javi", "CVID", "non-CVID")) %>% 
  mutate(Stage_CVID = paste0(PID, "_", Stage, "_", CovidSeverity)) %>% 
 # rename(FT = paste0("n", "1")) %>% 
  dplyr::select(Stage_CVID, pseudotime, n1, general_cluster, PID, Stage, CovidSeverity) %>%
  rownames_to_column("cell") 

plot.pseudo$Stage <- factor(plot.pseudo$Stage, levels = c("baseline", "progression", "convalescence"))

colors <-  c("CD56dim" = "#A6D854",
             "Inflamed_CD56dim" = "#FC8D62",
             "CD56bright" = "#FFD92F",
             "Cytokine_CD56dim" = "#8DA0CB",
             "SAP_CD56dim" = "#E78AC3",
             "HLAhigh_CD56dim" = "#66C2A5",
             "CVID_progression_Mild" = "#D9472A", 
             "non-CVID_progression_Mild" = "#D9472A", 
             "non-CVID_progression_Severe" = "#D9472A",
             "CVID_convalescence_Mild" = "#656ECD", 
             "non-CVID_convalescence_Mild" = "#656ECD", 
             "non-CVID_convalescence_Severe" = "#656ECD",
             "CVID_baseline_Control" = "#237A00", 
             "non-CVID_baseline_Control" = "#237A00", 
             "non-CVID_baseline_Mild" = "#237A00" #588B8B
)

plot.pseudo.baseline.noncvid <- plot.pseudo %>% 
  dplyr::filter(PID == "non-CVID", Stage == "baseline") %>% 
  mutate(Stage_CVID = "non-CVID_baseline_Mild") 


plot.pseudo.nuevo <- rbind(plot.pseudo, plot.pseudo.baseline.noncvid) %>% 
  mutate(PID_severity = "CVID") %>% 
  mutate(PID_severity = ifelse(CovidSeverity == "Mild" & PID == "non-CVID", "non-CVID_Mild",
                               ifelse(CovidSeverity == "Severe" & PID == "non-CVID", "non-CVID_Severe", PID_severity))) %>% 
  mutate(PID_severity = ifelse(Stage_CVID == "non-CVID_baseline_Mild", "non-CVID_Mild", PID_severity)) %>% 
  mutate(PID_severity = ifelse(Stage_CVID == "non-CVID_baseline_Control", "non-CVID_Severe", PID_severity))


ggplot(plot.pseudo.nuevo , aes(x = pseudotime, n1, color = general_cluster, fill = general_cluster)) +
  geom_point(alpha = 0.1, size = 0.75, shape = 16) +
  geom_smooth(aes(fill = Stage_CVID, color = Stage_CVID), linewidth = 0.75) +
  theme_light() +
  facet_grid(~PID_severity) +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  theme_pubr(border = T, legend = "right") +
  NoLegend() +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) +
  ylim(-0.5, 0.8) + 
  xlab("") + ylab("")


#Boxplots scores --------------

nk_1@meta.data <- nk_1@meta.data %>% 
  mutate(StageSeverity = paste0(Stage, "_", CovidSeverity),
         PID = ifelse(StudyName == "Javi", "CVID", "non-CVID"),
         PID_StageSeverity = ifelse(PID == "CVID", paste("CVID", StageSeverity, sep = "_"), paste("non-CVID", StageSeverity, sep = "_")), 
         Patient_Stage = paste0(patientID, "_", PID_StageSeverity)
  )

clusters <- c("CD56bright", "CD56dim", "HLAhigh_CD56dim")

## Cytotoxic response -----
markers <- read.table("markers/GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY.v2023.2.Hs.grp", header = T, sep = ",") 
#Cohort1
#NON-CVID
nk_1 <- subset(nk_1, subset = general_cluster %in% clusters)
non_cvid <- subset(nk_1, subset = StudyName == "China")
non_cvid <- AddModuleScore(non_cvid, assay = "RNA",
                           features =list(markers$GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY), #list(unique(markers$genes)), 
                           name = c("MARKER"), search = T)

non_cvid@meta.data <- non_cvid@meta.data %>% 
  mutate(Stage_n = ifelse(Stage == "baseline", "1baseline", ifelse(Stage == "convalescence", "3convalescence", "2progression")),
         PID = ifelse(StudyName == "Javi", "CVID", "non-CVID"),
         Stage_Severity = paste0(PID, "_", Stage, "_", CovidSeverity),
         Patient_Stage = paste0(Stage_Severity, ":", patientID, ":", Stage_n))

dot_noncvid <- data.frame()
for (n in clusters) {
  Idents(non_cvid) <- "Patient_Stage"
  a <- DotPlot(subset(non_cvid,  subset = general_cluster %in% n), features = "MARKER1", group.by = "Patient_Stage") +
    scale_color_gradientn(colors = c("skyblue3", "white", "#770000")) +
    theme_pubr(legend = "right") + border()+ rotate_x_text(45) +
    coord_flip() 
  
  a <- a$data %>% 
    mutate(cluster = n)
  dot_noncvid <- rbind(dot_noncvid, a) 
}  

dot_noncvid <- dot_noncvid %>% 
  separate(col = "id", sep = ":", into = c("STAGE", "resto"), remove = F, extra = "merge")  
#########################################################
###############################################
#CVID
cvid <- subset(nk_1, subset = StudyName == "Javi")

cvid <- AddModuleScore(cvid, assay = "RNA",
                       features = list(markers$GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY), #list(unique(markers$genes)), 
                       name = c("MARKER"), search = T)
cvid@meta.data <- cvid@meta.data %>% 
  mutate(Stage_n = ifelse(Stage == "baseline", "1baseline", ifelse(Stage == "convalescence", "3convalescence", "2progression")),
         PID = ifelse(StudyName == "Javi", "CVID", "non-CVID"),
         Stage_Severity = paste0(PID, "_", Stage, "_", CovidSeverity),
         Patient_Stage = paste0(Stage_Severity, ":", patientID, ":", Stage_n))

dot_cvid <- data.frame()
for (n in clusters) {
  Idents(cvid) <- "Patient_Stage"
  a <- DotPlot(subset(cvid,  subset = general_cluster %in% n), features = "MARKER1", group.by = "Patient_Stage") +
    scale_color_gradientn(colors = c("skyblue3", "white", "#770000")) +
    theme_pubr(legend = "right") + border()+ rotate_x_text(45) +
    coord_flip() 
  
  a <- a$data %>% 
    mutate(cluster = n)
  dot_cvid <- rbind(dot_cvid, a) 
}  

dot_cvid <- dot_cvid %>% 
  separate(col = "id", sep = ":", into = c("STAGE", "resto"), remove = F, extra = "merge") 
dot_coh1 <- rbind(dot_cvid, dot_noncvid)

writexl::write_xlsx(dot_coh1, path = "markers/nk_clusters/CYTOTOX_PATIENT_COH1_boxplot.xlsx")


subset_coh2 <- subset(nk_2, subset = general_cluster %in% clusters)
subset_coh2 <- subset(subset_coh2, subset = PID == "non-CVID") 

subset_coh2@meta.data <- subset_coh2@meta.data %>% 
  mutate(patient_stage = paste0(Stage, ":", patient_id))

subset_coh2 <- AddModuleScore(subset_coh2,
                              features = list(markers$GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY), #list(unique(markers$genes)), list(unique(markers$genes)), 
                              name = c("MARKERS"), search = T)

dot <- data.frame()
for (n in clusters) {
  Idents(subset_coh2) <- "patient_stage"
  a <- DotPlot(subset(subset_coh2,  subset = general_cluster %in% n),
               features = "MARKERS1") 
  
  a <- a$data %>% 
    mutate(cluster = n)
  dot <- rbind(dot, a) 
}  

dot_coh2 <- dot  %>% 
  separate(col = "id", sep = ":", into = c("STAGE", "resto"), remove = F, extra = "merge")

# BOXPLOT
dot_all <- rbind(dot_coh1, dot_coh2)

writexl::write_xlsx(dot_all, path = "markers/nk_clusters/CYTOX_PATIENT_ALLT_boxplot.xlsx")
dot_all <- readxl::read_xlsx(path = "markers/nk_clusters/CYTOX_PATIENT_ALLT_boxplot.xlsx")

dot_all$cluster <- factor(dot_all$cluster, levels = clusters)

dot_all$STAGE <- factor(dot_all$STAGE, levels = (c("CVID_baseline_Control", "CVID_progression_Mild", "CVID_convalescence_Mild",
                                                   "non-CVID_baseline_Control", "non-CVID_progression_Mild", "non-CVID_convalescence_Mild",
                                                   "non-CVID_progression_Severe", "non-CVID_convalescence_Severe",
                                                   "baseline", "progression_mild", "convalescence_mild",
                                                   "progression_severe", "convalescence_severe"
)))


# Fig 5g
colors <- c(RColorBrewer::brewer.pal(12, "Set3")[1:8], RColorBrewer::brewer.pal(12, "Set3")[4:8])
names(colors) <- levels(dot_all$STAGE)

ggplot(dot_all, aes(x=STAGE, y=avg.exp), fill = avg.exp, color = STAGE) +
  geom_boxplot(aes(fill = STAGE), outlier.size = 0.5) +
  geom_jitter(width = 0.1, size = 0.5) +
  theme_pubr() + border() + 
  facet_wrap(~cluster, scales = "free_y", ncol = 3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 6)) +
  scale_fill_manual(values = colors) +
  ylab(label = "Avg expression") + xlab(label = "Stage") + labs(title=" ") +
  geom_vline(xintercept = c(8.5), linetype="solid") +
  geom_vline(xintercept = c(3.5), linetype="dashed") +
  NoLegend() 

## Profibrotic response -----
markers <- c("AREG", "DUSP2", "ZFP36L2", "CXCR4", "TSC22D3")
#Cohort1
#NON-CVID
nk_1 <- subset(nk_1, subset = general_cluster %in% clusters)
non_cvid <- subset(nk_1, subset = StudyName == "China")
non_cvid <- AddModuleScore(non_cvid, assay = "RNA",
                           features =list(markers), #list(unique(markers$genes)), 
                           name = c("MARKER"), search = T)

non_cvid@meta.data <- non_cvid@meta.data %>% 
  mutate(Stage_n = ifelse(Stage == "baseline", "1baseline", ifelse(Stage == "convalescence", "3convalescence", "2progression")),
         PID = ifelse(StudyName == "Javi", "CVID", "non-CVID"),
         Stage_Severity = paste0(PID, "_", Stage, "_", CovidSeverity),
         Patient_Stage = paste0(Stage_Severity, ":", patientID, ":", Stage_n))

dot_noncvid <- data.frame()
for (n in clusters) {
  Idents(non_cvid) <- "Patient_Stage"
  a <- DotPlot(subset(non_cvid,  subset = general_cluster %in% n), features = "MARKER1", group.by = "Patient_Stage") +
    scale_color_gradientn(colors = c("skyblue3", "white", "#770000")) +
    theme_pubr(legend = "right") + border()+ rotate_x_text(45) +
    coord_flip() 
  
  a <- a$data %>% 
    mutate(cluster = n)
  dot_noncvid <- rbind(dot_noncvid, a) 
}  

dot_noncvid <- dot_noncvid %>% 
  separate(col = "id", sep = ":", into = c("STAGE", "resto"), remove = F, extra = "merge")  
#########################################################
###############################################
#CVID
cvid <- subset(nk_1, subset = StudyName == "Javi")

cvid <- AddModuleScore(cvid, assay = "RNA",
                       features = list(markers), #list(unique(markers$genes)), 
                       name = c("MARKER"), search = T)
cvid@meta.data <- cvid@meta.data %>% 
  mutate(Stage_n = ifelse(Stage == "baseline", "1baseline", ifelse(Stage == "convalescence", "3convalescence", "2progression")),
         PID = ifelse(StudyName == "Javi", "CVID", "non-CVID"),
         Stage_Severity = paste0(PID, "_", Stage, "_", CovidSeverity),
         Patient_Stage = paste0(Stage_Severity, ":", patientID, ":", Stage_n))

dot_cvid <- data.frame()
for (n in clusters) {
  Idents(cvid) <- "Patient_Stage"
  a <- DotPlot(subset(cvid,  subset = general_cluster %in% n), features = "MARKER1", group.by = "Patient_Stage") +
    scale_color_gradientn(colors = c("skyblue3", "white", "#770000")) +
    theme_pubr(legend = "right") + border()+ rotate_x_text(45) +
    coord_flip() 
  
  a <- a$data %>% 
    mutate(cluster = n)
  dot_cvid <- rbind(dot_cvid, a) 
}  

dot_cvid <- dot_cvid %>% 
  separate(col = "id", sep = ":", into = c("STAGE", "resto"), remove = F, extra = "merge") 
dot_coh1 <- rbind(dot_cvid, dot_noncvid)

writexl::write_xlsx(dot_coh1, path = "markers/nk_clusters/PROFIBROTIC_PATIENT_COH1_boxplot.xlsx")


subset_coh2 <- subset(nk_2, subset = general_cluster %in% clusters)
subset_coh2 <- subset(subset_coh2, subset = PID == "non-CVID") 

subset_coh2@meta.data <- subset_coh2@meta.data %>% 
  mutate(patient_stage = paste0(Stage, ":", patient_id))

subset_coh2 <- AddModuleScore(subset_coh2,
                              features = list(markers), #list(unique(markers$genes)), list(unique(markers$genes)), 
                              name = c("MARKERS"), search = T)

dot <- data.frame()
for (n in clusters) {
  Idents(subset_coh2) <- "patient_stage"
  a <- DotPlot(subset(subset_coh2,  subset = general_cluster %in% n),
               features = "MARKERS1") 
  
  a <- a$data %>% 
    mutate(cluster = n)
  dot <- rbind(dot, a) 
}  

dot_coh2 <- dot  %>% 
  separate(col = "id", sep = ":", into = c("STAGE", "resto"), remove = F, extra = "merge")

# BOXPLOT
dot_all <- rbind(dot_coh1, dot_coh2)

writexl::write_xlsx(dot_all, path = "markers/nk_clusters/FIBROTIC_PATIENT_ALLT_boxplot.xlsx")
dot_all <- readxl::read_xlsx(path = "markers/nk_clusters/FIBROTIC_PATIENT_ALLT_boxplot.xlsx")

dot_all$cluster <- factor(dot_all$cluster, levels = clusters)

dot_all$STAGE <- factor(dot_all$STAGE, levels = (c("CVID_baseline_Control", "CVID_progression_Mild", "CVID_convalescence_Mild",
                                                   "non-CVID_baseline_Control", "non-CVID_progression_Mild", "non-CVID_convalescence_Mild",
                                                   "non-CVID_progression_Severe", "non-CVID_convalescence_Severe",
                                                   "baseline", "progression_mild", "convalescence_mild",
                                                   "progression_severe", "convalescence_severe"
)))


# Supp Fig 5f
colors <- c(RColorBrewer::brewer.pal(12, "Set3")[1:8], RColorBrewer::brewer.pal(12, "Set3")[4:8])
names(colors) <- levels(dot_all$STAGE)

ggplot(dot_all, aes(x=STAGE, y=avg.exp), fill = avg.exp, color = STAGE) +
  geom_boxplot(aes(fill = STAGE), outlier.size = 0.5) +
  geom_jitter(width = 0.1, size = 0.5) +
  theme_pubr() + border() + 
  facet_wrap(~cluster, scales = "free_y", ncol = 3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 6)) +
  scale_fill_manual(values = colors) +
  ylab(label = "Avg expression") + xlab(label = "Stage") + labs(title=" ") +
  geom_vline(xintercept = c(8.5), linetype="solid") +
  geom_vline(xintercept = c(3.5), linetype="dashed") +
  NoLegend() 

#----------------------------

# Supp Fig 5e

## SCORES HEATMAP

# NK_CD56bright -----------
bright <- subset(nk_1, subset = general_cluster == "CD56bright", invert = F)
#GO:0009615_response to virus CLUSTER 1
IFN <- read.table("markers/nk_clusters/GOBP_RESPONSE_TO_TYPE_I_INTERFERON.v2023.1.Hs.grp", header = T)
bright <- AddModuleScore(bright, features = list(unique(IFN$GOBP_RESPONSE_TO_TYPE_I_INTERFERON)),
                         name = c("cl1_IFN_"), search = T)
#GO:0045954_positive regulation of natural killer cell mediated cytotoxicity CLUSTER 7
positive <- read.table("markers/nk_clusters/GOBP_POSITIVE_REGULATION_OF_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXCITY.Hs.grp", header = T)
bright <- AddModuleScore(bright, features = list(unique(positive$GENES)),
                         name = c("cl7_positive_"), search = T)
#GO:0032606_type I interferon production #CLUSTER 9 
production <- read.table("markers/nk_clusters/GOBP_TYPE_I_INTERFERON_PRODUCTION.v2023.2.Hs.grp", header = T)
bright <- AddModuleScore(bright, features = list(unique(production$GOBP_TYPE_I_INTERFERON_PRODUCTION)),
                        name = c("cl9_production_"), search = T)
#GO:0042267_natural killer cell mediated cytotoxicity CLUSTER 13
nk <- read.table("markers/nk_clusters/GO_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY.v2023.2.Hs.grp", header = T)
bright <- AddModuleScore(bright, features = list(unique(nk$GENE)),
                         name = c("cl13_nk_"), search = T)
#GO:0009615_response to virus
VIRUS <- read.table("markers/nk_clusters/GOBP_RESPONSE_TO_VIRUS.v2023.2.Hs.grp", header = T)
bright <- AddModuleScore(bright, features = list(unique(VIRUS$GOBP_RESPONSE_TO_VIRUS)),
                         name = c("VIRUS_"), search = T)

plot <- bright@meta.data %>%
  mutate(StudyName = ifelse(StudyName == "China", "nonCVID", "CVID")) %>% 
  mutate(sample = paste0(Stage, "_", CovidSeverity, "-", StudyName, ":", general_cluster)) %>%
  mutate(patient = paste0(sample, ":", StageCVID)) %>%
  select(general_cluster, sample, StageCVID, patient, cl1_IFN_1, cl7_positive_1, cl9_production_1, cl13_nk_1, VIRUS_1)


#Plot por los clusters de los scores
## Summary plot intento
#Mean per stage  
cl1_IFN_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl1_IFN_1), sd=sd(cl1_IFN_1)) %>% 
  mutate(cluster = "cl1_IFN_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl7_positive_1 <- plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl10_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl7_positive_1), sd=sd(cl7_positive_1)) %>% 
  mutate(cluster = "cl7_positive_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl9_production_1 <- plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl10_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl9_production_1), sd=sd(cl9_production_1)) %>% 
  mutate(cluster = "cl9_production_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl13_nk_1 <- plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl13_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl13_nk_1), sd=sd(cl13_nk_1)) %>% 
  mutate(cluster = "cl13_nk_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

VIRUS_1 <- plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl13_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(VIRUS_1), sd=sd(VIRUS_1)) %>% 
  mutate(cluster = "VIRUS_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))
# 
all_pop_plot_median <- rbind ( cl1_IFN_1 , cl7_positive_1, cl9_production_1, cl13_nk_1, VIRUS_1 )

all_pop_plot_median$stage <- factor(all_pop_plot_median$stage, levels = c("baseline", "progression", "convalescence"))
writexl::write_xlsx(all_pop_plot_median, path = "new_images_def/NEW/nkbright_GO_nocorrect.xlsx")

#PLOTS cl1_IFN_1 , cl7_positive_1, cl9_production_1, cl13_nk_1 VIRUS_1
VIRUS_1$PID_stage <- factor(VIRUS_1$PID_stage, levels = c("CVID:baseline_Control", "CVID:progression_Mild", "CVID:convalescence_Mild",
                                                          "nonCVID:baseline_Control", "nonCVID:progression_Mild", "nonCVID:convalescence_Mild",
                                                          "nonCVID:progression_Severe", "nonCVID:convalescence_Severe"))

ggplot(VIRUS_1, aes(x=PID_stage, y=cluster, fill = mean) ) +
  geom_tile(aes(fill=mean), color = "black")  +
  facet_grid(~cluster, scales = "free") +
  scale_fill_gradientn(colours =  c("#04508b", "#087ca7", "#E3B264","#FAE7C8")) +
  ggpubr::theme_pubr(legend = "right", border = T, x.text.angle = 45) 


# NK_CD56dim -----------------------
dim <- subset(nk_1, subset = general_cluster == "CD56dim", invert = F)

#GO:0034340_response to type I interferon CLUSTER 1, 5
IFN <- read.table("markers/nk_clusters/GOBP_RESPONSE_TO_TYPE_I_INTERFERON.v2023.1.Hs.grp", header = T)
dim <- AddModuleScore(dim, features = list(unique(IFN$GOBP_RESPONSE_TO_TYPE_I_INTERFERON)),
                      name = c("cl1_cl5_IFN_"), search = T)

#GO:0035455_response to interferon-alpha CLUSTER 9, 14
ALPHA <- read.table("markers/nk_clusters/GOBP_RESPONSE_TO_INTERFERON_ALPHA.v2023.2.Hs.grp", header = T)
dim <- AddModuleScore(dim, features = list(unique(ALPHA$GOBP_RESPONSE_TO_INTERFERON_ALPHA)),
                      name = c("cl9_cl14_ALPHA_"), search = T)

#GO:00422671_natural killer cell mediated cytotoxicity CLUSTER 12, 13
nk <- read.table("markers/nk_clusters/GO_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY.v2023.2.Hs.grp", header = T)
dim <- AddModuleScore(dim, features = list(unique(nk$GENE)),
                      name = c("cl12_cl13_nk_"), search = T)

#GO:0043320_natural killer cell degranulation #CLUSTER 12 
degranulation <- read.table("markers/nk_clusters/GOBP_NATURAL_KILLER_CELL_DEGRANULATION.grp", header = T)
dim <- AddModuleScore(dim, features = list(unique(degranulation$GO.0043320.naturalkillercelldegranulation)),
                      name = c("cl12_degranulation_"), search = T)

#GO:0009615_response to virus
virus <- read.table("markers/nk_clusters/GOBP_RESPONSE_TO_VIRUS.v2023.2.Hs.grp", header = T)
dim <- AddModuleScore(dim, features = list(unique(virus$GOBP_RESPONSE_TO_VIRUS)),
                      name = c("virus_"), search = T)

plot <- dim@meta.data %>%
  mutate(StudyName = ifelse(StudyName == "China", "nonCVID", "CVID")) %>% 
  mutate(sample = paste0(Stage, "_", CovidSeverity, "-", StudyName, ":", general_cluster)) %>%
  mutate(patient = paste0(sample, ":", StageCVID)) %>%
  select(general_cluster, sample, StageCVID, patient, cl1_cl5_IFN_1, cl9_cl14_ALPHA_1, cl12_cl13_nk_1, cl12_degranulation_1, virus_1)


#Plot por los clusters de los scores
## Summary plot intento
#Mean per stage  
cl1_cl5_IFN_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl1_cl5_IFN_1), sd=sd(cl1_cl5_IFN_1)) %>% 
  mutate(cluster = "cl1_cl5_IFN_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl9_cl14_ALPHA_1 <- plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl8_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl9_cl14_ALPHA_1), sd=sd(cl9_cl14_ALPHA_1)) %>% 
  mutate(cluster = "cl9_cl14_ALPHA_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl12_cl13_nk_1 <- plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl10_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl12_cl13_nk_1), sd=sd(cl12_cl13_nk_1)) %>% 
  mutate(cluster = "cl12_cl13_nk_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl12_degranulation_1 <- plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl10_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl12_degranulation_1), sd=sd(cl12_degranulation_1)) %>% 
  mutate(cluster = "cl12_degranulation_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

virus_1 <- plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl10_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(virus_1), sd=sd(virus_1)) %>% 
  mutate(cluster = "virus_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

all_pop_plot_median <- rbind ( cl1_cl5_IFN_1 , cl9_cl14_ALPHA_1, cl12_cl13_nk_1, cl12_degranulation_1, virus_1 )

all_pop_plot_median$stage <- factor(all_pop_plot_median$stage, levels = c("baseline", "progression", "convalescence"))
writexl::write_xlsx(all_pop_plot_median, path = "new_images_def/NEW/nkdim_GO_nocorrect.xlsx")

#PLOTS cl1_cl5_IFN_1 , cl9_cl14_ALPHA_1, cl12_cl13_nk_1, cl12_degranulation_1 virus_1
virus_1$PID_stage <- factor(virus_1$PID_stage, levels = c("CVID:baseline_Control", "CVID:progression_Mild", "CVID:convalescence_Mild",
                                                          "nonCVID:baseline_Control", "nonCVID:progression_Mild", "nonCVID:convalescence_Mild",
                                                          "nonCVID:progression_Severe", "nonCVID:convalescence_Severe"))

ggplot(virus_1, aes(x=PID_stage, y=cluster, fill = mean) ) +
  geom_tile(aes(fill=mean), color = "black")  +
  facet_grid(~cluster, scales = "free") +
  scale_fill_gradientn(colours = c("#04508b", "#087ca7", "#E3B264","#FAE7C8")) +
  ggpubr::theme_pubr(legend = "right", border = T, x.text.angle = 45) 


# HLAhigh_CD56dim --------------
hla <- subset(nk_1, subset = general_cluster == "HLAhigh_CD56dim", invert = F)

#GO:0034340_response to type I interferon CLUSTER 4
IFN <- read.table("markers/nk_clusters/GOBP_RESPONSE_TO_TYPE_I_INTERFERON.v2023.1.Hs.grp", header = T)
hla <- AddModuleScore(hla, features = list(unique(IFN$GOBP_RESPONSE_TO_TYPE_I_INTERFERON)),
                      name = c("cl4_IFN_"), search = T)

#GO:0042267_natural killer cell mediated cytotoxicity CLUSTER 3, 8, 13
nk <- read.table("markers/nk_clusters/GO_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY.v2023.2.Hs.grp", header = T)
hla <- AddModuleScore(hla, features = list(unique(nk$GENE)),
                      name = c("cl3_cl8_cl13_nk_"), search = T)

#GO:0035456_response to interferon-beta CLUSTER 6
beta <- read.table("markers/nk_clusters/GOBP_RESPONSE_TO_INTERFERON_BETA.v2023.2.Hs.grp", header = T)
hla <- AddModuleScore(hla, features = list(unique(beta$GOBP_RESPONSE_TO_INTERFERON_BETA)),
                      name = c("cl6_beta_"), search = T)


#GO:0043320_natural killer cell degranulation CLUSTER 13
degranulation <- read.table("markers/nk_clusters/GOBP_NATURAL_KILLER_CELL_DEGRANULATION.grp", header = T)
hla <- AddModuleScore(hla, features = list(unique(degranulation$GO.0043320.naturalkillercelldegranulation)),
                      name = c("cl13_degranulation_"), search = T)

#response to virus
virus <- read.table("markers/nk_clusters/GOBP_RESPONSE_TO_VIRUS.v2023.2.Hs.grp", header = T)
hla <- AddModuleScore(hla, features = list(unique(virus$GOBP_RESPONSE_TO_VIRUS)),
                      name = c("virus_"), search = T)

plot <- hla@meta.data %>%
  mutate(StudyName = ifelse(StudyName == "China", "nonCVID", "CVID")) %>% 
  mutate(sample = paste0(Stage, "_", CovidSeverity, "-", StudyName, ":", general_cluster)) %>%
  mutate(patient = paste0(sample, ":", StageCVID)) %>%
  select(general_cluster, sample, StageCVID, patient, cl4_IFN_1,cl3_cl8_cl13_nk_1, cl6_beta_1, cl13_degranulation_1, virus_1)

#Plot por los clusters de los scores
## Summary plot intento
#Mean per stage  
cl4_IFN_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl4_IFN_1), sd=sd(cl4_IFN_1)) %>% 
  mutate(cluster = "cl4_IFN_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl3_cl8_cl13_nk_1 <- plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl8_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl3_cl8_cl13_nk_1), sd=sd(cl3_cl8_cl13_nk_1)) %>% 
  mutate(cluster = "cl3_cl8_cl13_nk_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl6_beta_1 <- plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl10_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl6_beta_1), sd=sd(cl6_beta_1)) %>% 
  mutate(cluster = "cl6_beta_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl13_degranulation_1 <- plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl10_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl13_degranulation_1), sd=sd(cl13_degranulation_1)) %>% 
  mutate(cluster = "cl6_beta_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

virus_1 <- plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl10_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(virus_1), sd=sd(virus_1)) %>% 
  mutate(cluster = "virus_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))


all_pop_plot_median <- rbind ( cl4_IFN_1 , cl3_cl8_cl13_nk_1, cl6_beta_1, cl13_degranulation_1, virus_1 )

all_pop_plot_median$stage <- factor(all_pop_plot_median$stage, levels = c("baseline", "progression", "convalescence"))
writexl::write_xlsx(all_pop_plot_median, path = "new_images_def/NEW/hlaNK_GO_nocorrect.xlsx")

#PLOTS cl4_IFN_1, cl3_cl8_cl13_nk_1, cl6_beta_1, cl13_degranulation_1
cl4_IFN_1$PID_stage <- factor(cl4_IFN_1$PID_stage, levels = c("CVID:baseline_Control", "CVID:progression_Mild", "CVID:convalescence_Mild",
                                                          "nonCVID:baseline_Control", "nonCVID:progression_Mild", "nonCVID:convalescence_Mild",
                                                          "nonCVID:progression_Severe", "nonCVID:convalescence_Severe"))

ggplot(cl4_IFN_1, aes(x=PID_stage, y=cluster, fill = mean) ) +
  geom_tile(aes(fill=mean), color = "black")  +
  facet_grid(~cluster, scales = "free") +
  scale_fill_gradientn(colours = c("#04508b", "#087ca7", "#E3B264","#FAE7C8")) +
  ggpubr::theme_pubr(legend = "right", border = T, x.text.angle = 45) 
