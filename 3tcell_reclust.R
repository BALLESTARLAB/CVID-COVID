# T CELL RECLUSTERING COHORT 1 ---------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)

load("../COVID_REANNOT.RData")

Idents(covid_cvid) <- "integrated_snn_res.2"

tcell_1 <- subset(covid_cvid, idents = c(0, 1, 2, 15, 17, 18, 20:27, 31, 34, 35, 36, 38, 39, 40,
                                         6, 14, 29, 33, 37))

save(tcell_1, file = "subset_tcell_coh1.RData")

DefaultAssay(tcell_1) <- "RNA"
tcell.list <- SplitObject(tcell_1, split.by = "sample_names")
tcell.feature <- SelectIntegrationFeatures(object.list = tcell.list)
tcell.anchors <- FindIntegrationAnchors(object.list = tcell.list,
                                        anchor.features = tcell.feature, reference = 1)
tcell_1 <- IntegrateData(anchorset = tcell.anchors)

save(tcell_1, file = "tcell_cohort1_befscale.RData")

tcell_1 <- ScaleData(tcell_1)
tcell_1 <- RunPCA(tcell_1)
tcell_1 <- RunUMAP(tcell_1, reduction = "pca", dims = 1:30)
tcell_1 <- FindNeighbors(tcell_1, reduction = "pca", dims = 1:30)
tcell_1 <- FindClusters(tcell_1, resolution = 0.2, graph.name = "integrated_snn")

save(tcell_1, file = "tcell_cohort1_samplenames.RData")

# Annotation
tcell_1 <- FindClusters(tcell_1, resolution = 0.9, graph.name = "integrated_snn")

Idents(tcell_1) <- "seurat_clusters"
current.cluster <- levels(tcell_1)
current.cluster

new.cluster <- c("CD8_TEM", #cl 0 -
                 "CD4_naive", #cl 1 -
                 "CD4_TCM", #cl 2 -
                 "CD8_Teff", #cl 3
                 "CD4_TEM", #cl 4 -
                 "CD8_TCM", #cl 5 -
                 "Doublets", #cl 6 -
                 "CD4_CTL_1", #cl 7 -
                 "Proliferating_MKI67", #cl 8 -
                 "CD8_exhausted", #cl 9 -
                 "CD8_naive", #cl 10 - 
                 "CD8_NKT_like", #cl 11
                 "CD4_TCM", #cl 12 -
                 "TFH", #cl 13 -
                 "CD8_NKT_like", #cl 14 
                 "gdTC", #cl 15 -
                 "Treg", #cl 16 -
                 "Doublets", #cl 17 -
                 "CD4_CTL_2", #cl 18 -
                 "CD8_TEM", #cl 19 -
                 "MAIT", #cl 20 -
                 "Doublets", #cl 21 -
                 "CD4_CTL_1" #cl 22 -
)

names(x = new.cluster) <- levels(x = tcell_1)
tcell_1 <- RenameIdents(object = tcell_1, new.cluster)
tcell_1$general_cluster <- Idents(object = tcell_1)
table(tcell_1$general_cluster)
DimPlot(tcell_1, raster = F, label = T)

save(tcell_1, file = "newdata/integration_samples/cohort1/tcell/reannot/tcell_cohort1_annot.RData")

#Figure 3a
unique(tcell_1$general_cluster) %>% 
  as.vector() -> clusters

new.colors <- c("Doublets" = "#E7D68C", #"Doublets
                "CD4_CTL_2" = "#D38C22", #"CD4+ CTL 2"
                "CD4_CTL_1" = "#E36056", #"CD4+ CTL 1
                "Treg"= "#F19D7E", #"Treg",
                "CD4_TCM" = "#FF7033", #"CD4 TCM
                "CD4_TEM" = "#78759C", #"CD4_TEM
                "MAIT" = "#40A082", #MAIT
                "CD8_TCM" = "#4DD57F", #CD8_TEM_1
                "CD8_Teff" = "#DDB092", #cl3
                "gdTC" = "#C9B5DB", #gdTC
                "CD8_NKT_like" = "#8FC0E0", #cl11
                "CD8_naive" = "#AAD852", #"CD8+ Naive"
                "CD8_TEM" = "#8780DB", #CD8_TEM_2
                "Proliferating_MKI67" = "#BABF77", #Proliferating_MKI67
                "CD4_naive" = "#B496C7", #"CD4+ Naive"
                "CD8_exhausted" = "#FAD53E", #CD8_exhausted
                "TFH" = "#E28BC3" #"CD4_exhausted"
)

DimPlot(tcell_1, label = F, group.by = "general_cluster", pt.size = 0.2, raster = F, 
        order = c(  "MAIT", "gdTC", "CD4_naive", "CD8_naive", "CD8_exhausted", "CD8_TEM", "CD4_CTL_1",
                    "CD4_TEM", "CD4_TCM", "TFH", "Treg",
                    "CD8_TCM", "Proliferating_MKI67", "CD4_CTL_2", "Doublets"),
        cols = rev(new.colors)       ) + labs(title = "")


#Figure 3b
patchwork::wrap_plots( FeaturePlot(tcell_1,
                                   features = c("CD4", "CD8A", "CD38", 
                                                "FOXP3", "FCGR3A", "MKI67"), ncol = 3, label = F, repel = TRUE, keep.scale = "all", pt.size = 0.01, order =  T, raster = F)) &
  theme_minimal() & scale_color_viridis_c(option = "C", limits = c(0, 4)) & theme_pubr(border = F, legend = "none") &
  NoAxes()

# Supp figure 3a
tcell_1@meta.data <- tcell_1@meta.data %>% 
  mutate(PID = ifelse(StudyName == "Javi", "CVID", "non-CVID") ) %>% 
  select(-StudyName)
table(tcell_1$PID)

set.seed(123)

#CVID
sample1_rna <- subset(tcell_1, subset = PID == "CVID")
metadata1 <- as.data.frame(sample1_rna@meta.data)
nuevo_m1 <- row.names(sample_n(metadata1, 44769))
cvid_tcell_subsample <- subset(tcell_1, cells = nuevo_m1)

noncvid_tcell <- subset(tcell_1, subset = PID == "non-CVID") 
#non-CVID mild
sample2_rna <- subset(noncvid_tcell, subset = CovidSeverity == c("Severe"), invert = T)
metadata2 <- as.data.frame(sample2_rna@meta.data)
nuevo_m2 <- row.names(sample_n(metadata2, 44769))
noncvid_mild_tcell_subsample <- subset(tcell_1, cells = nuevo_m2)

#non-CVID severe
sample3_rna <- subset(noncvid_tcell, subset = CovidSeverity == c("Mild"), invert = T)
metadata3 <- as.data.frame(sample3_rna@meta.data)
nuevo_m3 <- row.names(sample_n(metadata3, 44769))
noncvid_severe_tcell_subsample <- subset(tcell_1, cells = nuevo_m3)

new.colors <- c("Doublets" = "#E7D68C", #"Doublets
                "CD4_CTL_2" = "#D38C22", #"CD4+ CTL 2"
                "CD4_CTL_1" = "#E36056", #"CD4+ CTL 1
                "Treg"= "#F19D7E", #"Treg",
                "CD4_TCM" = "#FF7033", #"CD4 TCM
                "CD4_TEM" = "#78759C", #"CD4_TEM
                "MAIT" = "#40A082", #MAIT
                "CD8_TCM" = "#4DD57F", #CD8_TEM_1
                "CD8_Teff" = "#DDB092", #cl3
                "gdTC" = "#C9B5DB", #gdTC
                "CD8_NKT_like" = "#8FC0E0", #cl11
                "CD8_naive" = "#AAD852", #"CD8+ Naive"
                "CD8_TEM" = "#8780DB", #CD8_TEM_2
                "Proliferating_MKI67" = "#BABF77", #Proliferating_MKI67
                "CD4_naive" = "#B496C7", #"CD4+ Naive"
                "CD8_exhausted" = "#FAD53E", #CD8_exhausted
                "TFH" = "#E28BC3" #"CD4_exhausted"
)
clusters <- c(  "MAIT", "gdTC", "CD4_naive", "CD8_naive", "CD8_exhausted", "CD8_TEM", "CD4_CTL_1",
                "CD4_TEM", "CD4_TCM", "TFH",  "Treg",
                "CD8_TCM", "Proliferating_MKI67", "CD4_CTL_2", "Doublets")

p1 <- DimPlot(cvid_tcell_subsample, group.by = "general_cluster", pt.size = 0.25, cols = new.colors, raster = F, order = clusters) + labs(title = "CVID") + NoLegend()
p2 <- DimPlot(noncvid_mild_tcell_subsample, group.by = "general_cluster", pt.size = 0.25, cols = new.colors, raster = F, order = clusters) + labs(title = "non-CVID mild") + NoLegend()
p3 <- DimPlot(noncvid_severe_tcell_subsample, group.by = "general_cluster", pt.size = 0.25, cols = new.colors, raster = F, order = clusters) + labs(title = "non-CVID severe")
p1+p2+p3

DimPlot(subset_tcell_1_subsample, pt.size = 0.25, group.by = "general_cluster", cols = new.colors, split.by = "PID", order = clusters, raster = F) + labs(title = "T cells (n = 112,288)")

# PHENOGRAPH
markers.tcell <- FindAllMarkers(tcell_1, assay = "RNA", logfc.threshold = 0.25, min.pct = 0.3)
writexl::write_xlsx(markers.tcell,  "newdata/integration_samples/cohort1/tcell/reannot/markers_pheno_reannotation.xlsx")

top.genes.all <- data.frame()
tcell_1@active.ident <- factor(tcell_1@active.ident, 
                               levels = rev(c("Treg", #cl 16
                                              "CD4_naive", #cl 1
                                              "CD4_TCM", #cl 13
                                              "CD4_TEM", #cl 2
                                              "CD4_CTL_1", #cl 12
                                              "CD4_CTL_2", #cl 4 
                                              "TFH",
                                              "Proliferating_MKI67", #cl 7 y 18 y 22
                                              "CD8_naive", #cl 18
                                              "CD8_TCM",
                                              "CD8_TEM",
                                              "CD8_Teff",
                                              "CD8_exhausted",
                                              "CD8_NKT_like",
                                              "gdTC",
                                              "MAIT",
                                              "Doublets" #cl 17
                               )))


tcell.clusters <- as.data.frame(table(tcell_1$general_cluster))$Var1 %>% as.vector()

for (i in tcell.clusters) {
  top.genes <- markers.tcell %>% 
    mutate(sign = ifelse(p_val_adj < 0.05 & cluster == i & pct.1 > 0.5 & avg_log2FC > 0, paste("sign", i, sep = "_"), "others")) %>% 
    dplyr::filter(sign != "others")
  top.genes.all <- rbind(top.genes.all, top.genes)
}

gene.markers.tcell <- (table(top.genes.all$gene)) %>% as.data.frame() %>% 
  dplyr::filter(Freq == 1) %>% 
  dplyr::filter(!grepl("^MT-", Var1))

cluster <- c("Treg", #cl 16
             "CD4_naive", #cl 1
             "CD4_TCM", #cl 13
             "CD4_TEM", #cl 2
             "CD4_CTL_1", #cl 12
             "CD4_CTL_2", #cl 4 
             "TFH",
             "Proliferating_MKI67", #cl 7 y 18 y 22
             "CD8_naive", #cl 18
             "CD8_TCM",
             "CD8_TEM",
             "CD8_Teff",
             "CD8_exhausted",
             "CD8_NKT_like",
             "gdTC",
             "MAIT",
             "Doublets" #cl 17
)

top.genes.common <- top.genes.all %>% 
  dplyr::filter(gene %in% gene.markers.tcell$Var1) %>% 
  group_by(cluster) %>%
  top_n(n = 15, wt = avg_log2FC) 

writexl::write_xlsx(top.genes.common.order, path = "new_images_def/cohort1/tcell/reannot/pheno_commongenes_reannot.xlsx")

top.genes.common.order <-  top.genes.common.order %>% arrange(factor( cluster, 
                                                                      levels = (c("Treg", #cl 16
                                                                                  "CD4_naive", #cl 1
                                                                                  "CD4_TCM", #cl 13
                                                                                  "CD4_TEM", #cl 2
                                                                                  "CD4_CTL_1", #cl 12
                                                                                  "CD4_CTL_2", #cl 4 
                                                                                  "TFH",
                                                                                  "Proliferating_MKI67", #cl 7 y 18 y 22
                                                                                  "CD8_naive", #cl 18
                                                                                  "CD8_TCM",
                                                                                  "CD8_TEM",
                                                                                  "CD8_Teff",
                                                                                  "CD8_exhausted",
                                                                                  "CD8_NKT_like",
                                                                                  "gdTC",
                                                                                  "MAIT",
                                                                                  "Doublets" #cl 17
                                                                      ))) )

top.genes.common.order.new <- c(top.genes.common.order$gene[1:5], "FAU", "UBA52", "CD37", "EEF1A1", "PIK3IP1",
                                top.genes.common.order$gene[11:75]) %>% as.vector()


top.genes.common.order <- readxl::read_xlsx(path = "new_images_def/cohort1/tcell/reannot/pheno_commongenes_reannot.xlsx")

#Supp Fig 3c
DefaultAssay(tcell_1) <- "RNA"
DotPlot(tcell_1, features = top.genes.common.order.new, scale.max = 100, scale.min = 0, dot.scale = 5)  +
  scale_color_gradientn(colors = c("skyblue3", "white", "#770000")) +
  theme_pubr(legend = "right") + border()+ rotate_x_text(45) +
  theme(axis.text.x=element_text(size=7))

# T CELL RECLUSTERING COHORT 2 -------------
load("../coh2_covid_samples_annot.RData")

Idents(covid_2) <- "general_cluster"
tcell <- subset(covid_2, idents = "T cells")
tcell@meta.data$seurat_clusters <- droplevels(tcell@meta.data$seurat_clusters)

save(tcell, file = "tcell_clean_new.RData")
Idents(tcell) <- "sample_names"
DefaultAssay(tcell) <- "RNA"
#rerun clustering
tcell.list <- SplitObject(tcell, split.by = "sample_names")
tcell.feature <- SelectIntegrationFeatures(object.list = tcell.list)
tcell.anchors <- FindIntegrationAnchors(object.list = tcell.list, anchor.features = tcell.feature)
tcell <- IntegrateData(anchorset = tcell.anchors)

save(tcell, file = "tcell_cohort2_befscale_clean.RData")

tcell <- ScaleData(tcell)
tcell <- RunPCA(tcell)
tcell <- harmony::RunHarmony(tcell, "sample_names")
tcell <- RunUMAP(tcell, reduction = "pca", dims = 1:30)
tcell <- FindNeighbors(tcell, reduction = "pca", dims = 1:30)
tcell_2 <- FindClusters(tcell, resolution = 0.7, graph.name = "integrated_snn")

save(tcell_2, file = "tcell_cohort2_samplenames_clean.RData")

Idents(tcell_2) <- "seurat_clusters"
current.cluster <- levels(tcell_2)
current.cluster
new.cluster <- c("CD8_TEM", #cl 0
                 "CD4_naive", #cl 1
                 "CD4_TCM", #cl 2
                 "CD8_NKT_like", #cl 3
                 "CD8_naive", #cl 4 
                 "CD4_TEM", #cl 5
                 "CD8_exhausted", #cl 6
                 "Doublets", #cl 7
                 "CD8_TCM", #cl 8
                 "CD8_Teff", #cl 9
                 "CD4_TCM", #cl 10
                 "MAIT", #cl 11
                 "CD4_TEM", #cl 12
                 "Sample_015_CV", #cl 13
                 "CD4_CTL_1", #cl 14 
                 "TFH", #cl 15 ESTAS SON TFH
                 "Treg", #cl 16
                 "Proliferating_MKI67", #cl 17
                 "Sample_015_CV", #cl 18
                 "CD4_CTL_2", #cl 19
                 "CD8_TCM", #cl 20
                 "Doublets", #cl 21
                 "CD4_naive" #cl 22
)

names(x = new.cluster) <- levels(x = tcell_2)
tcell_2 <- RenameIdents(object = tcell_2, new.cluster)
tcell_2$general_cluster <- Idents(object = tcell_2)

#Supp Fig 3b

tcell_2@meta.data <- tcell_2@meta.data %>% 
  mutate(seurat_clusters_07 = seurat_clusters, 
         old_general_cluster = general_cluster)
table(tcell_2$seurat_clusters_07)

tcell_2 <- FindClusters(tcell_2, resolution = 1, graph.name = "integrated_snn")

tcell_2@meta.data <- tcell_2@meta.data %>% 
  mutate(general_subcluster = ifelse(seurat_clusters_1 == "8", "gdTC", paste(general_cluster)))
table(tcell_2$general_subcluster)

new.colors <- c("Doublets" = "#E7D68C", #"Doublets
                "CD4_CTL_2" = "#D38C22", #"CD4+ CTL 2"
                "CD4_CTL_1" = "#E36056", #"CD4+ CTL 1
                "Treg"= "#F19D7E", #"Treg",
                "CD4_TCM" = "#FF7033", #"CD4 TCM
                "CD4_TEM" = "#78759C", #"CD4_TEM
                "MAIT" = "#40A082", #MAIT
                "CD8_TCM" = "#4DD57F", #CD8_TEM_1
                "CD8_Teff" = "#DDB092", #cl3
                "gdTC" = "#C9B5DB", #gdTC
                "CD8_NKT_like" = "#8FC0E0", #cl11
                "CD8_naive" = "#AAD852", #"CD8+ Naive"
                "CD8_TEM" = "#8780DB", #CD8_TEM_2
                "Proliferating_MKI67" = "#BABF77", #Proliferating_MKI67
                "CD4_naive" = "#B496C7", #"CD4+ Naive"
                "CD8_exhausted" = "#FAD53E", #CD8_exhausted
                "TFH" = "#E28BC3", #"CD4_exhausted"
                "Sample_015_CV" = "gray"
)

DimPlot(tcell_2, label = F, group.by = "general_subcluster", pt.size = 0.2, raster = F, 
        order = c( "TFH", "Treg", 
                   "MAIT", "CD4_naive", "CD8_naive",  "CD8_TEM", "CD4_CTL_1",
                   "CD4_TCM","CD8_NKT_like", "gdTC", "CD8_Teff","CD8_TCM", "CD4_TEM","CD8_exhausted",
                   "Proliferating_MKI67","CD4_CTL_2", "Doublets"),
        cols = rev(new.colors)       ) + labs(title = "")

# Supp Fig 3d
Idents(tcell_2) <- "general_subcluster"
tcell_2@active.ident <- factor(tcell_2@active.ident, 
                               levels = rev(c("Treg", #cl 16
                                              "CD4_naive", #cl 1
                                              "CD4_TCM", #cl 13
                                              "CD4_TEM", #cl 2
                                              "CD4_CTL_1", #cl 12
                                              "CD4_CTL_2", #cl 4 
                                              "TFH",
                                              "Proliferating_MKI67", #cl 7 y 18 y 22
                                              "CD8_naive", #cl 18
                                              "CD8_TCM",
                                              "CD8_TEM",
                                              "CD8_Teff",
                                              "CD8_exhausted",
                                              "CD8_NKT_like",
                                              "gdTC",
                                              "MAIT",
                                              "Doublets", #cl 17
                                              "Sample_015_CV"
                               )))


top.genes.common.order <- readxl::read_xlsx("new_images_def/cohort1/tcell/reannot/pheno_commongenes_reannot.xlsx")
top.genes.common.order.new <- c(top.genes.common.order$gene[1:5], "FAU", "UBA52", "CD37", "EEF1A1", "PIK3IP1",
                                top.genes.common.order$gene[11:70]) %>% as.vector()

DefaultAssay(tcell_2) <- "RNA"
DotPlot(subset(tcell_2, idents = c("Sample_015_CV", "Doublets"), invert = T), features = top.genes.common.order.new, scale.max = 100, scale.min = 0, dot.scale = 5)  +
  scale_color_gradientn(colors = c("skyblue3", "white", "#770000")) +
  theme_pubr(legend = "right") + border()+ rotate_x_text(45) +
  theme(axis.text.x=element_text(size=7))

#----------------------------
# T helper 
Idents(tcell_1) <- "general_cluster"

cd4 <- subset(tcell_1, idents = c("CD4_naive", "CD4_TCM", "CD4_TEM") )
save(cd4, file = "cd4_tcell.RData")

#Annotation of T helper
cd4 <- AddModuleScore(cd4, features = list(c("TBX21", "TNF", "IFNG")), name = "Th1_")
cd4 <- AddModuleScore(cd4, features = list(c("GATA3", "IL5", "IL4")), name = "Th2_")
cd4 <- AddModuleScore(cd4, features = list(c("RORC", "IL17F")), name = "Th17_")

score <- cd4@meta.data %>% 
  select(Th1_1, Th2_1, Th17_1) 

score_max <- data.frame()
for (row in 1:nrow(score)) {
  cell <- score[row, ] 
  
  a <- cell %>% mutate(column_max = max.col(cell), 
                       max_value = cell[,max.col(cell)]) 
  
  score_max <- rbind(score_max, a)
}

scores_def <- score_max %>% 
  mutate(Th = ifelse(column_max == 1, "Th1", 
                     ifelse(column_max == 2, "Th2",
                            ifelse(column_max == 3, "Th17", 
                                   "NA"))) ) %>% 
  mutate(Th_neg = ifelse(max_value <= 0.4, "Th0", Th))
table(scores_def$Th_neg)

cd4 <- AddMetaData(cd4, metadata = scores_def$Th_neg, col.name = "Th")

#Plot 
plot <- cd4@meta.data %>% 
  mutate(sample = paste0(Stage, "_", CovidSeverity, "_", StudyName, "_", general_cluster)) %>% 
  mutate(sample_th = paste0(Th, "_", sample)) %>% 
  mutate(patient = paste0(sample, ":", StageCVID)) %>% 
  select(general_cluster, Th, sample, sample_th, StageCVID, patient)

pop_tb <- table(plot$patient, plot$Th)
pop_tb <- as.data.frame.matrix(pop_tb) 
pop_tb$sample <- rownames(pop_tb)

pop_plot <- tidyr::gather(pop_tb, key = general_cluster, value = counts, -sample)

pop_plot <- pop_plot %>%
  group_by(sample,general_cluster) %>%
  summarise(counts = sum(counts)) %>%
  mutate(percent= counts*100/sum(counts)) %>%
  mutate(percent.text = paste0(sprintf("%.2f",percent),"%")) %>%
  mutate(position = 100-(cumsum(percent) - (0.5*percent))) 

pop_plot_sum <- pop_plot %>% 
  filter(general_cluster != "Th0") %>% 
  separate(sample, into = c("stage", "severity", "PID", "cluster_pat"), sep = "_", remove = F, extra = "merge") %>% 
  separate(cluster_pat, into = c("cluster", "patient"), sep = ":", remove = F) %>% 
  mutate(stage_sev = paste0(stage, "_", severity),
         PID_sev = ifelse(PID == "Javi", "CVID", "nonCVID_Mild") ) %>% 
  mutate(PID_sev= ifelse(severity == "Severe", "nonCVID_Severe", PID_sev),
         Th_cluster = paste0(cluster, ":", general_cluster, ":", stage_sev, ":", PID_sev))  %>% 
  filter(general_cluster %in% c("Th1", "Th2", "Th17"))

pop_plot_median <- plyr::ddply(pop_plot_sum, "Th_cluster", summarise, mean=mean(percent), sd=sd(percent)) %>% 
  separate(Th_cluster, into = c("cluster", "general_cluster", "stage_sev", "PID_Sev"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c("PID", "severity"), sep = "_", remove = F) %>% 
  separate(stage_sev, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")

rep_baseline <- pop_plot_median %>% 
  filter(PID == "nonCVID", stage_sev == "baseline_Control") %>% 
  mutate(stage_sev = paste0(stage_sev, "_Severe"),
         PID_Sev = "nonCVID_Severe") 

pop_plot_median <- rbind(pop_plot_median, rep_baseline) %>% 
  mutate(Th_stage = paste0(cluster, ":", general_cluster))

pop_plot_median <- pop_plot_median %>% 
  mutate ( stage = factor(stage, levels = c("baseline", "progression", "convalescence") )  )

# Figure 3e

# Dot plot
th <- subset(cd4, subset = Th %in% c("Th1", "Th2", "Th17"))

features <- c("TBX21", "IFNG", "TNF", "GATA3", "IL4", "IL5", "RORC", "IL17F")
data <- data.frame()
for (f in features) {
  
  Idents(th) <- "Th"
  a <- DotPlot(th, features = f) 
  a <- a$data %>% 
    mutate(markers = f)
  
  data <- rbind(data, a)
  
}

data$markers <- factor(data$markers, levels = rev(features) )

ggplot(data, aes(markers, id, fill = avg.exp.scaled, color = avg.exp.scaled)) +
  geom_tile() +
  scale_color_gradientn(colors = c("skyblue3", "white", "#770000")) +
  scale_fill_gradientn(colors = c("skyblue3", "white", "#770000")) +
  ggpubr::theme_pubr(border = T, legend = "right")

# Line graph
std.error <- function(x) sd(x)/sqrt(length(x)) #mean

pop_plot_median <- plyr::ddply(pop_plot_sum, "Th_cluster", summarise, mean=mean(percent), std=std.error(percent)) %>% 
  separate(Th_cluster, into = c("cluster", "general_cluster", "stage_sev", "PID_Sev"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c("PID", "severity"), sep = "_", remove = F) %>% 
  separate(stage_sev, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")

rep_baseline <- pop_plot_median %>% 
  filter(PID == "nonCVID", stage_sev == "baseline_Control") %>% 
  mutate(stage_sev = paste0(stage_sev, "_Severe"),
         PID_Sev = "nonCVID_Severe") 

pop_plot_median <- rbind(pop_plot_median, rep_baseline) %>% 
  mutate(Th_stage = paste0(cluster, ":", general_cluster))

pop_plot_median <- pop_plot_median %>% 
  mutate ( stage = factor(stage, levels = c("baseline", "progression", "convalescence") )  )

ggplot(pop_plot_median, aes( x = stage, y = mean, fill = general_cluster, color = general_cluster, group = Th_stage)) +
  geom_point(alpha = 1) +
  geom_line(alpha = 0.9)  +
  geom_errorbar(aes(ymin=mean-std, ymax=mean+std), width=.2) +
  facet_grid(cluster~PID_Sev) +
  scale_colour_manual(values = c("#2BA84A", "#F77F00","#815CAD")) +
  ylim(0,32) + 
  theme_pubr(border = T) 

#----------------------------
# CD4 CELLS  -------------------------

# Supp Fig 3e

## SCORES HEATMAP
# CD4 NAIVE --------------
cd4naive <- subset(tcell_1, subset = general_cluster == "CD4_naive", invert = F)

#GO:0060337_type I interferon signaling pathway CLUSTER 3
signaling <- read.table("markers/tcell_cluster/GOBP_TYPEI_IFN_SIGNALING_PATHWAY.grp", header = T)
cd4naive <- AddModuleScore(cd4naive, features = list(unique(signaling$GENE)),
                           name = c("cl4_signaling_"), search = T)

#GO:0032606_type I interferon production CLUSTER 7
production <- read.table("markers/tcell_cluster/GOBP_TYPE_I_INTERFERON_PRODUCTION.v2023.2.Hs.grp", header = T)
cd4naive <- AddModuleScore(cd4naive, features = list(unique(production$GOBP_TYPE_I_INTERFERON_PRODUCTION)),
                           name = c("cl7_production_"), search = T)

#GO:0050852_T cell receptor signaling pathway CLUSTER 15 Y 16 
TCR <- read.table("markers/tcell_cluster/GO:0050852.grp", header = T)
cd4naive <- AddModuleScore(cd4naive, features = list(unique(TCR$GENE)),
                           name = c("cl15_cl16_TCR_"), search = T)

#GO:0009615GOBP_RESPONSE_TO_VIRUS CLUSTER 8,3 y 2
virus <- read.table("markers/tcell_cluster/GOBP_RESPONSE_TO_VIRUS.v2023.2.Hs.grp", header = T)
cd4naive <- AddModuleScore(cd4naive, features = list(unique(virus$GOBP_RESPONSE_TO_VIRUS)),
                           name = c("virus_"), search = T)

plot <- cd4naive@meta.data %>%
  mutate(StudyName = ifelse(StudyName == "China", "nonCVID", "CVID")) %>% 
  mutate(sample = paste0(Stage, "_", CovidSeverity, "-", StudyName, ":", general_cluster)) %>%
  mutate(patient = paste0(sample, ":", StageCVID)) %>%
  select(general_cluster, sample, StageCVID, patient, cl1_differentiation_1, cl4_signaling_1, cl7_production_1, cl8_activation_1, cl15_cl16_TCR_1, virus_1)

#Plot por los clusters de los scores
## Summary plot  
# #Mean per stage  

#TYPE I IFN SIGNALING PATHWAY
cl4_signaling_1 <- plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl10_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl4_signaling_1), sd=sd(cl4_signaling_1)) %>% 
  mutate(cluster = "cl4_signaling_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl7_production_1 <-plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl13_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl7_production_1), sd=sd(cl7_production_1)) %>% 
  mutate(cluster = "cl7_production_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl15_cl16_TCR_1 <-plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl13_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl15_cl16_TCR_1), sd=sd(cl15_cl16_TCR_1)) %>% 
  mutate(cluster = "cl15_cl16_TCR_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

virus_1 <-plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl13_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(virus_1), sd=sd(virus_1)) %>% 
  mutate(cluster = "virus_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

all_pop_plot_median <- rbind (cl1_differentiation_1, cl4_signaling_1, cl7_production_1, cl8_activation_1, cl15_cl16_TCR_1, virus_1 )

all_pop_plot_median$stage <- factor(all_pop_plot_median$stage, levels = c("baseline", "progression", "convalescence"))
writexl::write_xlsx(all_pop_plot_median, path = "new_images_def/NEW/cd4naive_GO_nocorrect.xlsx")

#PLOTS  cl4_signaling_1, cl7_production_1, cl15_cl16_TCR_1 virus_1
cl15_cl16_TCR_1$PID_stage <- factor(cl15_cl16_TCR_1$PID_stage, levels = c("CVID:baseline_Control", "CVID:progression_Mild", "CVID:convalescence_Mild",
                                                          "nonCVID:baseline_Control", "nonCVID:progression_Mild", "nonCVID:convalescence_Mild",
                                                          "nonCVID:progression_Severe", "nonCVID:convalescence_Severe"))

ggplot(cl15_cl16_TCR_1, aes(x=PID_stage, y=cluster, fill = mean) ) +
  geom_tile(aes(fill=mean), color = "black")  +
  facet_grid(~cluster, scales = "free") +
  scale_fill_gradientn(colours = c("#04508b", "#087ca7", "#E3B264","#FAE7C8")) +
  ggpubr::theme_pubr(legend = "right", border = T, x.text.angle = 45) 

# CD4 TCM --------------
cd4tcm <- subset(tcell_1, subset = general_cluster == "CD4_TCM", invert = F)

#GO:0034340_response to type I interferon CLUSTER 2
IFN <- read.table("markers/tcell_cluster/GOBP_RESPONSE_TO_TYPE_I_INTERFERON.v2023.1.Hs.grp", header = T)
cd4tcm <- AddModuleScore(cd4tcm, features = list(unique(IFN$GOBP_RESPONSE_TO_TYPE_I_INTERFERON)),
                         name = c("cl2_IFN_"), search = T)

#GO:0035456_response to interferon-beta CLUSTER 10
BETA <- read.table("markers/tcell_cluster/GOBP_RESPONSE_TO_INTERFERON_BETA.v2023.2.Hs.grp", header = T)
cd4tcm <- AddModuleScore(cd4tcm, features = list(unique(BETA$GOBP_RESPONSE_TO_INTERFERON_BETA)),
                         name = c("cl10_beta_"), search = T)

#GO:0050852_T cell receptor signaling pathway cluster 15
TCR <- read.table("markers/tcell_cluster/GO:0050852.grp", header = T)
cd4tcm <- AddModuleScore(cd4tcm, features = list(unique(TCR$GENE)),
                         name = c("cl15_TCR_"), search = T)

#GO:0009615GOBP_RESPONSE_TO_VIRUS CLUSTER 
VIRUS <- read.table("markers/tcell_cluster/GOBP_RESPONSE_TO_VIRUS.v2023.2.Hs.grp", header = T)
cd4tcm <- AddModuleScore(cd4tcm, features = list(unique(VIRUS$GOBP_RESPONSE_TO_VIRUS)),
                         name = c("VIRUS_"), search = T)

#GO:0072678_T cell migration
migration <- read.table("markers/tcell_cluster/GO:0072678.grp", header = T)
cd4tcm <- AddModuleScore(cd4tcm, features = list(unique(migration$gene)),
                         name = c("migration_"), search = T)

plot <- cd4tcm@meta.data %>%
  mutate(StudyName = ifelse(StudyName == "China", "nonCVID", "CVID")) %>% 
  mutate(sample = paste0(Stage, "_", CovidSeverity, "-", StudyName, ":", general_cluster)) %>%
  mutate(patient = paste0(sample, ":", StageCVID)) %>%
  select(general_cluster, sample, StageCVID, patient, cl2_IFN_1, cl10_beta_1, cl7_differentiation_1, cl15_TCR_1, VIRUS_1, migration_1)

#Plot por los clusters de los scores
## Summary plot  
#Mean per stage  
cl2_IFN_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl2_IFN_1), sd=sd(cl2_IFN_1)) %>% 
  mutate(cluster = "cl2_IFN_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl10_beta_1 <- plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl10_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl10_beta_1), sd=sd(cl10_beta_1)) %>% 
  mutate(cluster = "cl10_beta_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl15_TCR_1 <-plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl13_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl15_TCR_1), sd=sd(cl15_TCR_1)) %>% 
  mutate(cluster = "cl15_TCR_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

VIRUS_1 <-plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl13_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(VIRUS_1), sd=sd(VIRUS_1)) %>% 
  mutate(cluster = "VIRUS_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

migration_1 <-plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl13_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(migration_1), sd=sd(migration_1)) %>% 
  mutate(cluster = "migration_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

all_pop_plot_median <- rbind (cl2_IFN_1, cl10_beta_1, cl7_differentiation_1, cl15_TCR_1, VIRUS_1, migration_1 )

all_pop_plot_median$stage <- factor(all_pop_plot_median$stage, levels = c("baseline", "progression", "convalescence"))
writexl::write_xlsx(all_pop_plot_median, path = "new_images_def/NEW/cd4TCM_GO_nocorrect.xlsx")

#PLOTS cl2_IFN_1, cl10_beta_1, cl7_differentiation_1, cl15_TCR_1 VIRUS_1 migration_1
cl10_beta_1$PID_stage <- factor(cl10_beta_1$PID_stage, levels = c("CVID:baseline_Control", "CVID:progression_Mild", "CVID:convalescence_Mild",
                                                          "nonCVID:baseline_Control", "nonCVID:progression_Mild", "nonCVID:convalescence_Mild",
                                                          "nonCVID:progression_Severe", "nonCVID:convalescence_Severe"))

ggplot(cl10_beta_1, aes(x=PID_stage, y=cluster, fill = mean) ) +
  geom_tile(aes(fill=mean), color = "black")  +
  facet_grid(~cluster, scales = "free") +
  scale_fill_gradientn(colours = c("#04508b", "#087ca7", "#E3B264","#FAE7C8")) +
  ggpubr::theme_pubr(legend = "right", border = T, x.text.angle = 45) 

# CD4 TEM --------------
cd4tem <- subset(tcell_1, subset = general_cluster == "CD4_TEM", invert = F)

#GO:0034340_response to type I interferon CLUSTER 3
IFN <- read.table("markers/tcell_cluster/GOBP_RESPONSE_TO_TYPE_I_INTERFERON.v2023.1.Hs.grp", header = T)
cd4tem <- AddModuleScore(cd4tem, features = list(unique(IFN$GOBP_RESPONSE_TO_TYPE_I_INTERFERON)),
                         name = c("cl3_IFN_"), search = T)

#GO:0042093_T-helper cell differentiation CLUSTER 8
th <- read.table("markers/tcell_cluster/GO:0042093.grp", header = T)
cd4tem <- AddModuleScore(cd4tem, features = list(unique(th$GENE)),
                         name = c("cl8_TH_"), search = T)

#GO:0050852_T cell receptor signaling pathway CLUSTER 16
signaling <- read.table("markers/tcell_cluster/GO:0050852.grp", header = T)
cd4tem <- AddModuleScore(cd4tem, features = list(unique(signaling$GENE)),
                         name = c("cl16_signaling_"), search = T)

#GO:0060337_type I interferon signaling pathway CLUSTER 17
typei <- read.table("markers/tcell_cluster/GOBP_TYPEI_IFN_SIGNALING_PATHWAY.grp", header = T)
cd4tem <- AddModuleScore(cd4tem, features = list(unique(typei$GENE)),
                         name = c("cl17_typei_"), search = T)

#GO:0009615GOBP_RESPONSE_TO_VIRUS CLUSTER 8,3 y 2
virus <- read.table("markers/tcell_cluster/GOBP_RESPONSE_TO_VIRUS.v2023.2.Hs.grp", header = T)
cd4tem <- AddModuleScore(cd4tem, features = list(unique(virus$GOBP_RESPONSE_TO_VIRUS)),
                         name = c("virus_"), search = T)

plot <- cd4tem@meta.data %>%
  mutate(StudyName = ifelse(StudyName == "China", "nonCVID", "CVID")) %>% 
  mutate(sample = paste0(Stage, "_", CovidSeverity, "-", StudyName, ":", general_cluster)) %>%
  mutate(patient = paste0(sample, ":", StageCVID)) %>%
  select(general_cluster, sample, StageCVID, patient, cl3_IFN_1, cl8_TH_1, cl16_signaling_1, cl13_differentiation_1, cl17_typei_1, virus_1)

#Plot por los clusters de los scores
## Summary plot  
#Mean per stage  
cl3_IFN_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl3_IFN_1), sd=sd(cl3_IFN_1)) %>% 
  mutate(cluster = "cl3_IFN_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl8_TH_1 <- plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl10_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl8_TH_1), sd=sd(cl8_TH_1)) %>% 
  mutate(cluster = "cl8_TH_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl16_signaling_1 <-plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl13_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl16_signaling_1), sd=sd(cl16_signaling_1)) %>% 
  mutate(cluster = "cl16_signaling_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl17_typei_1 <-plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl13_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl17_typei_1), sd=sd(cl17_typei_1)) %>% 
  mutate(cluster = "cl17_typei_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))


virus_1 <-plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl13_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(virus_1), sd=sd(virus_1)) %>% 
  mutate(cluster = "virus_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

all_pop_plot_median <- rbind (cl3_IFN_1,cl8_TH_1,cl16_signaling_1,  cl17_typei_1, virus_1)

all_pop_plot_median$stage <- factor(all_pop_plot_median$stage, levels = c("baseline", "progression", "convalescence"))
writexl::write_xlsx(all_pop_plot_median, path = "new_images_def/NEW/cd4TEM_GO_nocorrect.xlsx")

#PLOTS cl3_IFN_1, cl8_TH_1, ,cl16_signaling_1,  cl17_typei_1
cl13_differentiation_1$PID_stage <- factor(cl13_differentiation_1$PID_stage, levels = c("CVID:baseline_Control", "CVID:progression_Mild", "CVID:convalescence_Mild",
                                                          "nonCVID:baseline_Control", "nonCVID:progression_Mild", "nonCVID:convalescence_Mild",
                                                          "nonCVID:progression_Severe", "nonCVID:convalescence_Severe"))

ggplot(cl13_differentiation_1, aes(x=PID_stage, y=cluster, fill = mean) ) +
  geom_tile(aes(fill=mean), color = "black")  +
  facet_grid(~cluster, scales = "free") +
  scale_fill_gradientn(colours =  c("#04508b", "#087ca7", "#E3B264","#FAE7C8")) +
  ggpubr::theme_pubr(legend = "right", border = T, x.text.angle = 45) 

## -----------------
#IFNG SCORE
# Supp Fig 3g 

#COH 1 
cd4 <- subset(tcell_1, idents = c("CD4_naive", "CD4_TCM", "CD4_TEM"))

cd4@meta.data <- cd4@meta.data %>% 
  mutate(Stage_PID = paste0(StudyName, "_", Stage, "_", CovidSeverity))

a <- DotPlot(cd4, features = "IFNG", group.by = "Stage_PID", assay = "RNA", split.by = "general_cluster", cols = c("red", "blue", "green")) 
CD4_ifng <- a$data %>% 
  separate(id, into = c("PID", "stage", "severity", "cluster"), sep = "_", remove = F, extra = "merge") %>% 
  mutate(stage_order = paste0(PID, "_", stage, "_", severity))

CD4_ifng$stage_order <- factor(CD4_ifng$stage_order, levels = c("Javi_baseline_Control", "Javi_progression_Mild", "Javi_convalescence_Mild",
                                                                "China_baseline_Control", "China_progression_Mild", "China_convalescence_Mild",
                                                                "China_progression_Severe", "China_convalescence_Severe"))

CD4_ifng$cluster <- factor(CD4_ifng$cluster, levels = rev(c("CD4_naive", "CD4_TCM", "CD4_TEM")) )

writexl::write_xlsx(a$data, path = "markers/CD4_IFNG_COH1.xlsx")

#COH 2
Idents(tcell_2) <- "general_subcluster"
tcell_2 <- subset(tcell_2, idents = c("CD4_naive", "CD4_TCM", "CD4_TEM"))

tcell_2 <- subset(tcell_2, subset = PID == "non-CVID")
tcell_2@meta.data <- tcell_2@meta.data %>% 
  mutate(Stage = ifelse(Stage == "baseline", "baseline_control", paste(Stage)) ) %>% 
  mutate(Stage_PID = paste0("coh2_", Stage)) 

b <- DotPlot(tcell_2, features = "IFNG", group.by = "Stage_PID", assay = "RNA", split.by = "general_cluster", cols = c("red", "blue", "green")) 
CD4_ifng_coh2 <- b$data %>% 
  separate(id, into = c("PID", "stage", "severity", "cluster"), sep = "_", remove = F, extra = "merge") %>% 
  mutate(stage_order = paste0(PID, "_", stage, "_", severity)) 

CD4_ifng_coh2$stage_order <- factor(CD4_ifng_coh2$stage_order, levels = c("coh2_baseline_control", "coh2_progression_mild", "coh2_convalescence_mild",
                                                                          "coh2_progression_severe", "coh2_convalescence_severe"))

CD4_ifng_coh2$cluster <- factor(CD4_ifng_coh2$cluster, levels = rev(c("CD4_naive", "CD4_TCM", "CD4_TEM")) )

writexl::write_xlsx(b$data, path = "markers/CD4_IFNG_COH2.xlsx")
# CD4_ifng_coh2 <- readxl::read_xlsx("../bsc/nk_reclustering/scores/scores/markers/CD4_IFNG_COH2.xlsx")
# CD4_ifng <- readxl::read_xlsx("../bsc/nk_reclustering/scores/scores/markers/CD4_IFNG_COH1.xlsx")

all <- rbind(CD4_ifng, CD4_ifng_coh2) %>% 
  separate(id, into = c("PID", "stage", "severity", "cluster"), sep = "_", remove = F, extra = "merge") %>% 
  mutate(stage_order = paste0(PID, "_", stage, "_", severity)) 
  
ggplot(all, aes(x=stage_order, y=cluster, size = pct.exp, fill = avg.exp, color = avg.exp)) +
  geom_point(shape = 21) +
  theme_pubr(legend = "right") + border() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  scale_fill_gradientn(colors = c( "white", "darkred")) +
  scale_color_gradientn(colors = c("white", "darkred")) +
  ylab(label = "Clusters") + xlab(label = "Stage") + labs(title="IFNG") +
  geom_vline(xintercept = c(3.5), linetype="dashed") +
  geom_vline(xintercept = c(8.5), linetype="solid") 


#######################
#######################

# CD8 CELLS ----------------

#Cytototoxic response 
DefaultAssay(tcell_1) <- "RNA"
tcell_1@meta.data <- tcell_1@meta.data %>% 
  mutate(StageSeverity = paste0(Stage, "_", CovidSeverity),
         PID = ifelse(StudyName == "Javi", "CVID", "non-CVID"),
         PID_StageSeverity = ifelse(PID == "CVID", paste("CVID", StageSeverity, sep = "_"), paste("non-CVID", StageSeverity, sep = "_")), 
         Patient_Stage = paste0(patientID, "_", PID_StageSeverity)
  )

clusters_cd8 <- c("CD8_naive", "CD8_TCM", "CD8_TEM")
#MARKERS
markers <- read.table("/markers/GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY.v2023.2.Hs.grp", header = T, sep = ",") 

#BOXPLOT PER PATIENT COH1
tcell_1 <- subset(tcell_1, subset = general_cluster %in% clusters_cd8)
#NON-CVID
non_cvid <- subset(tcell_1, subset = StudyName == "China")
non_cvid <- AddModuleScore(non_cvid, assay = "RNA",
                           features = list(markers$GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY), name = c("MARKER"), search = T)

non_cvid@meta.data <- non_cvid@meta.data %>% 
  mutate(Stage_n = ifelse(Stage == "baseline", "1baseline", ifelse(Stage == "convalescence", "3convalescence", "2progression")),
         PID = ifelse(StudyName == "Javi", "CVID", "non-CVID"),
         Stage_Severity = paste0(PID, "_", Stage, "_", CovidSeverity),
         Patient_Stage = paste0(Stage_Severity, ":", patientID, ":", Stage_n))

dot_noncvid <- data.frame()
for (n in clusters_cd8) {
  Idents(non_cvid) <- "Patient_Stage"
  a <- DotPlot(subset(non_cvid,  subset = general_cluster %in% n), features = "MARKER1", group.by = "Patient_Stage") 
  a <- a$data %>% 
    mutate(cluster = n)
  dot_noncvid <- rbind(dot_noncvid, a) 
}  

dot_noncvid <- dot_noncvid %>% 
  separate(col = "id", sep = ":", into = c("STAGE", "resto"), remove = F, extra = "merge")  

#CVID
cvid <- subset(tcell_1, subset = StudyName == "Javi")
cvid <- AddModuleScore(cvid, assay = "RNA",
                           features = list(markers$GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY), name = c("MARKER"), search = T)

cvid@meta.data <- cvid@meta.data %>% 
  mutate(Stage_n = ifelse(Stage == "baseline", "1baseline", ifelse(Stage == "convalescence", "3convalescence", "2progression")),
         PID = ifelse(StudyName == "Javi", "CVID", "non-CVID"),
         Stage_Severity = paste0(PID, "_", Stage, "_", CovidSeverity),
         Patient_Stage = paste0(Stage_Severity, ":", patientID, ":", Stage_n))

dot_cvid <- data.frame()
for (n in clusters_cd8) {
  Idents(cvid) <- "Patient_Stage"
  a <- DotPlot(subset(cvid,  subset = general_cluster %in% n), features = "MARKER1", group.by = "Patient_Stage") 
  a <- a$data %>% 
    mutate(cluster = n)
  dot_cvid <- rbind(dot_cvid, a) 
}  

dot_cvid <- dot_cvid %>% 
  separate(col = "id", sep = ":", into = c("STAGE", "resto"), remove = F, extra = "merge") 
dot_coh1 <- rbind(dot_cvid, dot_noncvid)

writexl::write_xlsx(dot_coh1, path = "markers/outs/cd8_NEW_COH1_boxplot.xlsx")

# BOXPLOT PER PATIENT COH1
subset_coh2 <- subset(tcell_2, subset = general_subcluster %in% clusters_cd8)
subset_coh2 <- subset(subset_coh2, subset = PID == "non-CVID") 

subset_coh2@meta.data <- subset_coh2@meta.data %>% 
  mutate(patient_stage = paste0(Stage, ":", patient_id))

subset_coh2 <- AddModuleScore(subset_coh2,
                              features = list(markers$GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY
                              ), 
                              name = c("MARKERS"), search = T)

dot <- data.frame()
for (n in clusters_cd8) {
  Idents(subset_coh2) <- "patient_stage"
  a <- DotPlot(subset(subset_coh2,  subset = general_subcluster %in% n),
               features = "MARKERS1") 
  
  a <- a$data %>% 
    mutate(cluster = n)
  dot <- rbind(dot, a)
}  

dot_coh2 <- dot  %>% 
  separate(col = "id", sep = ":", into = c("STAGE", "resto"), remove = F, extra = "merge")

# BOXPLOT
dot_all <- rbind(dot_coh1, dot_coh2)

writexl::write_xlsx(dot_all, path = "markers/outs/tcell/ifn_PATIENT_ALL_PATIENT_boxplot_cd4.xlsx")
dot_all <- readxl::read_xlsx("markers/outs/tcell/cd8_ifn_PATIENT_ALL_PATIENT_boxplot.xlsx")
dot_all$cluster <- factor(dot_all$cluster, levels = clusters_cd8)

dot_all$STAGE <- factor(dot_all$STAGE, levels = (c("CVID_baseline_Control", "CVID_progression_Mild", "CVID_convalescence_Mild",
                                                   "non-CVID_baseline_Control", "non-CVID_progression_Mild", "non-CVID_convalescence_Mild",
                                                   "non-CVID_progression_Severe", "non-CVID_convalescence_Severe",
                                                   "baseline", "progression_mild", "convalescence_mild",
                                                   "progression_severe", "convalescence_severe"
)))

colors <- c(RColorBrewer::brewer.pal(12, "Set3")[1:8], RColorBrewer::brewer.pal(12, "Set3")[4:8])
names(colors) <- levels(dot_all$STAGE)

#Figure 4c 
ggplot(dot_all, aes(x=STAGE, y=avg.exp.scaled), fill = avg.exp.scaled, color = STAGE) +
  geom_boxplot(aes(fill = STAGE), outlier.size = 0.5) +
  geom_jitter(width = 0.1, size = 0.5) +
  theme_pubr() + border() + 
  facet_wrap(~cluster, scales = "free_y",
             ncol = 3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 6)) +
  scale_fill_manual(values = colors) +
  ylab(label = "Score") + xlab(label = "Stage") + labs(title=" ") +
  geom_vline(xintercept = c(8.5), linetype="solid") +
  geom_vline(xintercept = c(3.5), linetype="dashed") +
  NoLegend() 

############ -------
#TCR Clonotypes
load("newdata/integration_samples/cohort1/tcell/reannot/tcell_cohort1_annot.RData")

tcell_1@meta.data <- tcell_1@meta.data %>% 
  mutate(PID = ifelse(StudyName == "Javi", "CVID", "non-CVID"),
         Patient_Stage = paste0(patientID, "_", Stage)) 

## CVID patients ------------------
Idents(tcell_1) <- "StudyName"

cvid <- subset(tcell_1, subset = StudyName == "Javi")
S1 <- read.csv("newdata/VDJ/vdj_t_filtered_contig_annotations.csv") %>% #78096 cells total
  filter(barcode %in% colnames(cvid)) #58817 only in T cell

annot <- cvid@meta.data %>% 
  tibble::rownames_to_column("barcode") %>% 
  filter(barcode %in% S1$barcode) %>% select(barcode, general_cluster, Stage, CovidSeverity, PID, patientID, orig.ident, Patient_Stage) 

S1 <- left_join(S1, annot, "barcode")

#LIST OF CONTIG FOR EACH SAMPLE
contig.list <- createHTOContigList(S1, cvid, group.by = "Patient_Stage")

#WE NEED TO SIMPLIFY THE BARCODES NAMES
samples <- names(contig.list)
contig_list <- list()
for (i in samples) {
  contig_list[[i]] <- separate(contig.list[[i]], col = "barcode", into = c("sample", "barcode"), sep = "-0_") %>% 
    select(-sample)
  # contig_list[[i]] <- stripBarcode(contig.list[[i]], column = 1, connector = "0_", num_connects = 3)
}

saveRDS(contig_list, file = "vdj/tcr_noncvid_contiglist.RDS")

# COMBINING THE CONTIGS
combined <- combineTCR(contig_list, 
                       samples = samples, 
                       removeNA = FALSE,
                       removeMulti = TRUE,
                       filterMulti = TRUE) 

combined_order <- combined[order(match(names(combined), c("CVID1_baseline", "CVID1_progression", "CVID1_convalescence",
                                                          "CVID2_baseline", "CVID2_progression", "CVID2_convalescence",
                                                          "CVID3_baseline", "CVID3_progression", "CVID3_convalescence",
                                                          "CVID4_baseline", "CVID4_progression", "CVID4_convalescence",
                                                          "CVID5_baseline", "CVID5_progression", "CVID5_convalescence")) )]

combined_order <- addVariable(combined_order, 
                              name = "Stage", 
                              variables = rep(c("baseline", "progression", "convalescence"), 5))

df <- plyr::ldply(contig_list, data.frame) %>% 
  mutate(barcode = paste0(Patient_Stage, "_", barcode)) %>% 
  select(barcode, general_cluster)

combined_order_cluster <- plyr::ldply(combined_order, data.frame) %>% 
  left_join(df, by = "barcode") 

################
clusters <- combined_order_cluster$general_cluster %>% unique()
plot <- data.frame()
for (n in clusters) {
  combined_cluster <- combined_order_cluster %>% 
    filter(general_cluster == n) %>% 
    split(.$.id)
  
  a <- clonalDiversity(combined_cluster, 
                       cloneCall = "strict", 
                       group.by = "sample", 
                       x.axis = "Stage", 
                       n.boots = 100, return.boots = F)

  a <- a$data %>% 
    mutate(cluster = n)
  
  plot <- rbind(plot, a)
  
}

clusters.cd8 <- c("CD8_naive", "CD8_TCM", "CD8_TEM", "CD8_Teff", "Proliferating_MKI67", "CD8_exhausted")

plot_filter <- plot %>% 
  filter(variable == "Shannon",
         cluster %in% clusters.cd8)

plot_filter$Stage <- factor(plot_filter$Stage, levels = c("baseline", "progression", "convalescence"))

colors <- c(RColorBrewer::brewer.pal(12, "Set3")[1:8], RColorBrewer::brewer.pal(12, "Set3")[4:8])
names(colors) <- levels(plot_filter$Stage)

# Supp Fig 4e ------------
a <- ggplot(plot_filter, aes(x=Stage, y=value), color = Stage) +
  #ylim(0.5, 1.8) +
  geom_boxplot(aes(fill = Stage), # outlier.size = 0.5, 
               outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.5) +
  theme_pubr() + border() +
  facet_wrap(~cluster, scales = "free_y", 
             ncol = 14
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 6),
        axis.text.y = element_text(hjust = 1, vjust = 1, size = 9)) +
  scale_fill_manual(values = colors) +
  
  ylab(label = "Avg expression") + xlab(label = "Stage") + labs(title=" ") +
  NoLegend() 

# ADD TO SEURAT OBJECT
meta.data <- cvid@meta.data %>%
  tibble::rownames_to_column("barcode.new") %>% separate(col = "barcode.new", into = c("sample", "barcode"), sep = "-0_", extra = "drop") %>% 
  mutate(barcodes = paste0(patientID, "_", Stage, "_", barcode)) %>% 
  tibble::column_to_rownames("barcodes")

cvid <- RenameCells(cvid, new.names = rownames(meta.data))

seurat <- combineExpression(combined, cvid, 
                            cloneCall="strict", 
                            group.by = "sample", 
                            proportion = FALSE, 
                            cloneTypes = c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500) )

table(seurat$cloneType)



##Single Cell Experiment Format
seurat@meta.data <- seurat@meta.data %>% 
  mutate(cloneType = ifelse(is.na(cloneType) == T, "NA", cloneType)) %>% 
  mutate(cloneType = ifelse(cloneType == "1", "Hyperexpanded (100 < X <= 500)", 
                            ifelse(cloneType == "2", "Large (20 < X <= 100)", 
                                   ifelse(cloneType == "3", "Medium (5 < X <= 20)", 
                                          ifelse(cloneType == "4", "Small (1 < X <= 5)", "Single (0 < X <= 1)")
                                          ))) )


cols <-  c("Hyperexpanded (100 < X <= 500)" = "#ffbe0b",
           "Large (20 < X <= 100)" = "#fb5607",
           "Medium (5 < X <= 20)" = "#ff006e",
           "Single (0 < X <= 1)" = "#3a86ff",
           "Small (1 < X <= 5)" = "#8338ec",
           "NA" = "#d2d2d2")

## SUBSAMPLING
min <- table(seurat$Stage) %>%  min()

all_metadata <- vector()
stages <- table(seurat$Stage) %>% names()
for (c in stages) {
  sample_n <- subset(seurat, subset = Stage == c)
  metadata <- as.data.frame(sample_n@meta.data)
  nuevo <- row.names(sample_n(metadata, min))
  all_metadata <- c(all_metadata, nuevo) 
}

subset <- subset(seurat, cells = all_metadata)

# Figure 4e  -------------
DimPlot(subset, group.by = "cloneType", cols = cols, order = c("Hyperexpanded (100 < X <= 500)",
                                                               "Large (20 < X <= 100)",
                                                               "Medium (5 < X <= 20)",
                                                               "Small (1 < X <= 5)",
                                                               "Single (0 < X <= 1)",
                                                               "NA" ),
        split.by = "Stage")



# BARPLOT 
#Select populations with expansions (CD8+)
metadata <- seurat@meta.data %>% 
  mutate(nStage = ifelse(Stage == "baseline", "1baseline",
                         ifelse(Stage == "progression", "2progression",
                                "3convalescence"))) %>% 
  mutate(cluster_stage = paste0(general_cluster, ":", nStage)) 
metadata <- metadata %>% 
  filter(cloneType != "NA")

pop_tb <- table(metadata$cluster_stage, metadata$cloneType)
pop_tb <- as.data.frame.matrix(pop_tb) 
pop_tb$sample <- rownames(pop_tb)


pop_plot <- tidyr::gather(pop_tb, key = general_cluster, value = counts, -sample)
plot_order <- c("Hyperexpanded (100 < X <= 500)",
                "Large (20 < X <= 100)", 
                "Medium (5 < X <= 20)", 
                "Small (1 < X <= 5)",
                "Single (0 < X <= 1)",
                "NA" ) 
pop_plot$general_cluster <- factor(pop_plot$general_cluster, levels = plot_order)


pop_plot <- pop_plot %>%
  group_by(sample,general_cluster) %>%
  summarise(counts = sum(counts)) %>%
  mutate(percent= counts*100/sum(counts)) %>%
  mutate(percent.text = paste0(sprintf("%.2f",percent),"%")) %>%
  mutate(position = 100-(cumsum(percent) - (0.5*percent))) 

pop_plot_filter <- pop_plot %>% 
  separate(sample, into = c("cluster", "stage"), sep = ":", remove = F) %>% 
  filter(!cluster %in% c("Doublets", "CD4_CTL_2", "gdTC"))


plot_order <- c('CD4_naive:1baseline', "CD4_naive:2progression", "CD4_naive:3convalescence",
                'CD4_TCM:1baseline',"CD4_TCM:2progression", "CD4_TCM:3convalescence",
                'CD4_TEM:1baseline',"CD4_TEM:2progression", "CD4_TEM:3convalescence",
                "Treg:1baseline",   "Treg:2progression","Treg:3convalescence",
                "TFH:1baseline",    "TFH:2progression", "TFH:3convalescence",
                "MAIT:1baseline",   "MAIT:2progression","MAIT:3convalescence",
                "CD8_naive:1baseline", "CD8_naive:2progression", "CD8_naive:3convalescence",
                "CD8_TCM:1baseline","CD8_TCM:2progression", "CD8_TCM:3convalescence",          
                "CD8_TEM:1baseline","CD8_TEM:2progression", "CD8_TEM:3convalescence",
                "CD8_Teff:1baseline","CD8_Teff:2progression",              "CD8_Teff:3convalescence",
                "Proliferating_MKI67:1baseline",      "Proliferating_MKI67:2progression",   "Proliferating_MKI67:3convalescence",
                "CD4_CTL_1:1baseline","CD4_CTL_1:2progression",             "CD4_CTL_1:3convalescence",
                "CD8_NKT_like:1baseline",            "CD8_NKT_like:2progression",          "CD8_NKT_like:3convalescence",
                'CD8_exhausted:1baseline',            "CD8_exhausted:2progression",         "CD8_exhausted:3convalescence"   
)

pop_plot_filter <- pop_plot_filter %>% 
  mutate ( sample = factor(sample, levels = plot_order )  )

#Fig 4e barplot
library(RColorBrewer)
library(ggpubr)
ggplot(pop_plot_filter, aes(x= sample, y = percent)) +
  geom_bar(aes(fill = general_cluster), stat = "identity", width = 0.5, color = "black") +
  scale_fill_manual(values = cols) +
  geom_vline(xintercept = seq(3.5, 40, 3), linetype="dashed") +
  theme_pubr(legend = "right") + border() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab(label = "Percentage") + xlab(label = "")

# Supp Fig 4a -------------
## SCORES HEATMAP
#CD8 NAIVE -----------------
cd8naive <- subset(tcell_1, subset = general_cluster == "CD8_naive", invert = F)

#GO:0060337_type I interferon signaling pathway CLUSTER 6 Y 7
signaling <- read.table("markers/tcell_cluster/GOBP_TYPEI_IFN_SIGNALING_PATHWAY.grp", header = T)
cd8naive <- AddModuleScore(cd8naive, features = list(unique(signaling$GENE)),
                           name = c("cl6_cl7_signaling_"), search = T)

#GO:0035455_response to interferon-alpha CLUSTER 14
alpha <- read.table("markers/tcell_cluster/GOBP_RESPONSE_TO_INTERFERON_ALPHA.v2023.2.Hs.grp", header = T)
cd8naive <- AddModuleScore(cd8naive, features = list(unique(alpha$GOBP_RESPONSE_TO_INTERFERON_ALPHA)),
                           name = c("cl14_alpha_"), search = T)

#GO:0050863_regulation of T cell activation CLUSTER 16 
activation <- read.table("markers/tcell_cluster/GO:0050863.grp", header = T)
cd8naive <- AddModuleScore(cd8naive, features = list(unique(activation$gene)),
                           name = c("cl16_activation_"), search = T)

#GO:0050852_T cell receptor signaling pathway CLUSTER 17 
TCR <- read.table("markers/tcell_cluster/GO:0050852.grp", header = T)
cd8naive <- AddModuleScore(cd8naive, features = list(unique(TCR$GENE)),
                           name = c("cl17_TCR_"), search = T)

#GO:0009615GOBP_RESPONSE_TO_VIRUS CLUSTER 8,3 y 2
virus <- read.table("markers/tcell_cluster/GOBP_RESPONSE_TO_VIRUS.v2023.2.Hs.grp", header = T)
cd8naive <- AddModuleScore(cd8naive, features = list(unique(virus$GOBP_RESPONSE_TO_VIRUS)),
                           name = c("virus_"), search = T)

plot <- cd8naive@meta.data %>%
  mutate(StudyName = ifelse(StudyName == "China", "nonCVID", "CVID")) %>% 
  mutate(sample = paste0(Stage, "_", CovidSeverity, "-", StudyName, ":", general_cluster)) %>%
  mutate(patient = paste0(sample, ":", StageCVID)) %>%
  select(general_cluster, sample, StageCVID, patient, cl6_cl7_signaling_1, cl14_alpha_1, cl16_activation_1, cl17_TCR_1, virus_1 )

#Plot por los clusters de los scores
## Summary plot  
#Mean per stage  
cl6_cl7_signaling_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl6_cl7_signaling_1), sd=sd(cl6_cl7_signaling_1)) %>% 
  mutate(cluster = "cl6_cl7_signaling_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl14_alpha_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl14_alpha_1), sd=sd(cl14_alpha_1)) %>% 
  mutate(cluster = "cl14_alpha_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl16_activation_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl16_activation_1), sd=sd(cl16_activation_1)) %>% 
  mutate(cluster = "cl16_activation_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl17_TCR_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl17_TCR_1), sd=sd(cl17_TCR_1)) %>% 
  mutate(cluster = "cl17_TCR_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

virus_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(virus_1), sd=sd(virus_1)) %>% 
  mutate(cluster = "virus_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))
# 
all_pop_plot_median <- rbind (cl6_cl7_signaling_1, cl14_alpha_1, cl16_activation_1, cl17_TCR_1, virus_1 )

all_pop_plot_median$stage <- factor(all_pop_plot_median$stage, levels = c("baseline", "progression", "convalescence"))
writexl::write_xlsx(all_pop_plot_median, path = "new_images_def/NEW/cd8naive_GO_nocorrect.xlsx")

#PLOTS  cl6_cl7_signaling_1, cl14_alpha_1, cl16_activation_1, cl17_TCR_1 virus_1
cl6_cl7_signaling_1$PID_stage <- factor(cl6_cl7_signaling_1$PID_stage, levels = c("CVID:baseline_Control", "CVID:progression_Mild", "CVID:convalescence_Mild",
                                                          "nonCVID:baseline_Control", "nonCVID:progression_Mild", "nonCVID:convalescence_Mild",
                                                          "nonCVID:progression_Severe", "nonCVID:convalescence_Severe"))

ggplot(cl6_cl7_signaling_1, aes(x=PID_stage, y=cluster, fill = mean) ) +
  geom_tile(aes(fill=mean), color = "black")  +
  facet_grid(~cluster, scales = "free") +
  scale_fill_gradientn(colours =  c("#04508b", "#087ca7", "#E3B264","#FAE7C8")) +
  ggpubr::theme_pubr(legend = "right", border = T, x.text.angle = 45) 

#CD8 TCM ---------------------
cd8tcm <- subset(tcell_1, subset = general_cluster == "CD8_TCM", invert = F)

#GO:0034340_response to type I interferon CLUSTER 6 
IFN <- read.table("markers/tcell_cluster/GOBP_RESPONSE_TO_TYPE_I_INTERFERON.v2023.1.Hs.grp", header = T)
cd8tcm <- AddModuleScore(cd8tcm, features = list(unique(IFN$GOBP_RESPONSE_TO_TYPE_I_INTERFERON)),
                         name = c("cl6_IFN_"), search = T)

#GO:0044194_cytolytic granule CLUSTER 7 Y 8
granule <- read.table("markers/tcell_cluster/GO:0044194.grp", header = T)
cd8tcm <- AddModuleScore(cd8tcm, features = list(unique(granule$gene)),
                         name = c("cl7_8_granule_"), search = T)

#GGO:0042608_T cell receptor binding CLUSTER 9
TCR <- read.table("markers/tcell_cluster/GO:0042608.grp", header = T)
cd8tcm <- AddModuleScore(cd8tcm, features = list(unique(TCR$gene)),
                         name = c("cl9_tcr_"), search = T)

#GO:0051607_defense response to virus CLUSTER 8
defense <- read.table("markers/tcell_cluster/GO:0051607.grp", header = T)
cd8tcm <- AddModuleScore(cd8tcm, features = list(unique(defense$GENE)),
                         name = c("defense_"), search = T)

plot <- cd8tcm@meta.data %>%
  mutate(StudyName = ifelse(StudyName == "China", "nonCVID", "CVID")) %>% 
  mutate(sample = paste0(Stage, "_", CovidSeverity, "-", StudyName, ":", general_cluster)) %>%
  mutate(patient = paste0(sample, ":", StageCVID)) %>%
  select(general_cluster, sample, StageCVID, patient, cl6_IFN_1, cl7_8_granule_1, cl9_tcr_1, defense_1 )

#Plot por los clusters de los scores
## Summary plot 
#Mean per stage  cl6_IFN_1, cl7_8_granule_1, cl9_tcr_1
cl6_IFN_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl6_IFN_1), sd=sd(cl6_IFN_1)) %>% 
  mutate(cluster = "cl6_IFN_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl7_8_granule_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl7_8_granule_1), sd=sd(cl7_8_granule_1)) %>% 
  mutate(cluster = "cl7_8_granule_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl9_tcr_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl9_tcr_1), sd=sd(cl9_tcr_1)) %>% 
  mutate(cluster = "cl9_tcr_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

defense_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(defense_1), sd=sd(defense_1)) %>% 
  mutate(cluster = "defense_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))
# 
# 
all_pop_plot_median <- rbind (cl6_IFN_1, cl7_8_granule_1, cl9_tcr_1, defense_1)

all_pop_plot_median$stage <- factor(all_pop_plot_median$stage, levels = c("baseline", "progression", "convalescence"))
writexl::write_xlsx(all_pop_plot_median, path = "new_images_def/NEW/cd8tcm_GO_nocorrect.xlsx")

#PLOTS  cl6_IFN_1, cl7_8_granule_1, cl9_tcr_1 defense_1
cl6_IFN_1$PID_stage <- factor(cl6_IFN_1$PID_stage, levels = c("CVID:baseline_Control", "CVID:progression_Mild", "CVID:convalescence_Mild",
                                                              "nonCVID:baseline_Control", "nonCVID:progression_Mild", "nonCVID:convalescence_Mild",
                                                              "nonCVID:progression_Severe", "nonCVID:convalescence_Severe"))

ggplot(cl6_IFN_1, aes(x=PID_stage, y=cluster, fill = mean) ) +
  geom_tile(aes(fill=mean), color = "black")  +
  facet_grid(~cluster, scales = "free") +
  scale_fill_gradientn(colours = c("#04508b", "#087ca7", "#E3B264","#FAE7C8")) +
  ggpubr::theme_pubr(legend = "right", border = T, x.text.angle = 45) 


#CD8 TEM ---------------------
cd8tem <- subset(tcell_1, subset = general_cluster == "CD8_TEM", invert = F)

#GO:0034340_response to type I interferon CLUSTER 3
IFN <- read.table("markers/tcell_cluster/GOBP_RESPONSE_TO_TYPE_I_INTERFERON.v2023.1.Hs.grp", header = T)
cd8tem <- AddModuleScore(cd8tem, features = list(unique(IFN$GOBP_RESPONSE_TO_TYPE_I_INTERFERON)),
                         name = c("cl3_IFN_"), search = T)

#GO:0071357_cellular response to type I interferon CLUSTER 6
cellular <- read.table("markers/tcell_cluster/GO:0071357.grp", header = T)
cd8tem <- AddModuleScore(cd8tem, features = list(unique(cellular$gene)),
                         name = c("cl6_cellular_"), search = T)

#GO:0001906_cell killing cluster 13 y 16
killing <- read.table("markers/tcell_cluster/GOBP_CELL_KILLING.v2023.2.Hs.grp", header = T)
cd8tem <- AddModuleScore(cd8tem, features = list(unique(killing$GOBP_CELL_KILLING)),
                         name = c("cl13_16_killing_"), search = T)

#GO:0044194_cytolytic granule 13 
granule <- read.table("markers/tcell_cluster/GO:0044194.grp", header = T)
cd8tem <- AddModuleScore(cd8tem, features = list(unique(granule$gene)),
                         name = c("cl13_granule_"), search = T)

#GO:0035456_response to interferon-beta CLUSTER 17
BETA <- read.table("markers/tcell_cluster/GOBP_RESPONSE_TO_INTERFERON_BETA.v2023.2.Hs.grp", header = T)
cd8tem <- AddModuleScore(cd8tem, features = list(unique(BETA$GOBP_RESPONSE_TO_INTERFERON_BETA)),
                         name = c("cl17_BETA_"), search = T)

#GO:0001909_leukocyte mediated cytotoxicity CLUSTER 16
leukocyte <- read.table("markers/tcell_cluster/GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY.v2023.2.Hs.grp", header = T)
cd8tem <- AddModuleScore(cd8tem, features = list(unique(leukocyte$GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY)),
                         name = c("cl16_leukocyte_"), search = T)

#GO:0050863_regulation of T cell activation CLUSTER 18
activation <- read.table("markers/tcell_cluster/GO:0050863.grp", header = T)
cd8tem <- AddModuleScore(cd8tem, features = list(unique(activation$gene)),
                         name = c("cl18_activation_"), search = T)

#GO:0009615GOBP_RESPONSE_TO_VIRUS CLUSTER 8,3 y 2
virus <- read.table("markers/tcell_cluster/GOBP_RESPONSE_TO_VIRUS.v2023.2.Hs.grp", header = T)
cd8tem <- AddModuleScore(cd8tem, features = list(unique(virus$GOBP_RESPONSE_TO_VIRUS)),
                         name = c("virus_"), search = T)


plot <- cd8tem@meta.data %>%
  mutate(StudyName = ifelse(StudyName == "China", "nonCVID", "CVID")) %>% 
  mutate(sample = paste0(Stage, "_", CovidSeverity, "-", StudyName, ":", general_cluster)) %>%
  mutate(patient = paste0(sample, ":", StageCVID)) %>%
  select(general_cluster, sample, StageCVID, patient, cl3_IFN_1,  cl6_cellular_1, cl13_16_killing_1, cl13_granule_1, cl17_BETA_1, cl16_leukocyte_1, cl18_activation_1, virus_1)

#Plot por los clusters de los scores
## Summary plot  
#Mean per stage  cl3_IFN_1,  cl6_cellular_1, cl13_16_killing_1, cl13_granule_1, cl17_BETA_1, cl16_leukocyte_1, cl18_activation_1
cl3_IFN_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl3_IFN_1), sd=sd(cl3_IFN_1)) %>% 
  mutate(cluster = "cl3_IFN_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl6_cellular_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl6_cellular_1), sd=sd(cl6_cellular_1)) %>% 
  mutate(cluster = "cl6_cellular_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl13_16_killing_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl13_16_killing_1), sd=sd(cl13_16_killing_1)) %>% 
  mutate(cluster = "cl13_16_killing_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl13_granule_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl13_granule_1), sd=sd(cl13_granule_1)) %>% 
  mutate(cluster = "cl13_granule_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl16_leukocyte_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl16_leukocyte_1), sd=sd(cl16_leukocyte_1)) %>% 
  mutate(cluster = "cl16_leukocyte_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl17_BETA_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl17_BETA_1), sd=sd(cl17_BETA_1)) %>% 
  mutate(cluster = "cl17_BETA_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl18_activation_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl18_activation_1), sd=sd(cl18_activation_1)) %>% 
  mutate(cluster = "cl18_activation_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

virus_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(virus_1), sd=sd(virus_1)) %>% 
  mutate(cluster = "virus_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))
# 
all_pop_plot_median <- rbind (cl3_IFN_1,  cl6_cellular_1, cl13_16_killing_1, cl13_granule_1, cl17_BETA_1, cl16_leukocyte_1, cl18_activation_1, virus_1)

all_pop_plot_median$stage <- factor(all_pop_plot_median$stage, levels = c("baseline", "progression", "convalescence"))
writexl::write_xlsx(all_pop_plot_median, path = "new_images_def/NEW/cd8tem_GO_nocorrect.xlsx")

#PLOTS  cl3_IFN_1,  cl6_cellular_1, cl13_16_killing_1, cl13_granule_1, cl17_BETA_1, cl16_leukocyte_1, cl18_activation_1
cl13_16_killing_1$PID_stage <- factor(cl13_16_killing_1$PID_stage, levels = c("CVID:baseline_Control", "CVID:progression_Mild", "CVID:convalescence_Mild",
                                                          "nonCVID:baseline_Control", "nonCVID:progression_Mild", "nonCVID:convalescence_Mild",
                                                          "nonCVID:progression_Severe", "nonCVID:convalescence_Severe"))

ggplot(cl13_16_killing_1, aes(x=PID_stage, y=cluster, fill = mean) ) +
  geom_tile(aes(fill=mean), color = "black")  +
  facet_grid(~cluster, scales = "free") +
  scale_fill_gradientn(colours = c("#04508b", "#087ca7", "#E3B264","#FAE7C8")) +
  ggpubr::theme_pubr(legend = "right", border = T, x.text.angle = 45) 


# TC percentages compared to Georg paper --------
load("~/Escritorio/datos_covid/newdata/integration_samples/cohort1/tcell/tcell_cohort1_annot.RData")

georg <- c("CD38+ MKI67++ FCGR3A+ CM")

# COHORT 1 
#Percentages T cells per cluster 
tcell_metadata <- tcell_1@meta.data %>% 
  mutate(StageSeverity = paste0(Stage, "_", CovidSeverity),
         PID = ifelse(StudyName == "Javi", "CVID", "non-CVID"),
         PID_StageSeverity = ifelse(PID == "CVID", paste("CVID", StageSeverity, sep = "_"), paste("non-CVID", StageSeverity, sep = "_")), 
         Patient_Stage = paste0(patientID, "_", PID_StageSeverity)
  )

#General general_cluster
pop_tb <- table(tcell_metadata$Patient_Stage, tcell_metadata$general_cluster)
pop_tb <- as.data.frame.matrix(pop_tb)
pop_tb$sample <- rownames(pop_tb)
pop_plot <- tidyr::gather(pop_tb, key = general_cluster, value = counts, -sample)

samples <- unique(pop_plot$sample)
pop_plot_new <- data.frame()
for (patient in samples) {
  
  pops <- pop_plot %>%
    filter(sample == patient) %>% 
    group_by(sample,general_cluster) %>%
    summarise(counts = sum(counts)) %>%
    mutate(percent= counts*100/sum(counts)) %>%
    mutate(percent.text = paste0(sprintf("%.2f",percent),"%")) %>%
    mutate(position = 100-(cumsum(percent) - (0.5*percent))) %>% 
    mutate(PID = ifelse(str_detect(sample, "non-CVID"), "non-CVID", "CVID") ) %>% 
    mutate(Severity = ifelse(str_detect(sample, "baseline"), "baseline", 
                             (ifelse(str_detect(sample, "Mild"), "Mild", "Severe")) ) ) %>% 
    mutate(Stage = ifelse(str_detect(sample, "progression"), "progression", 
                          (ifelse(str_detect(sample, "baseline"), "", "convalescence")) ) )  %>% 
    mutate(Stage_Severity = paste0(PID, "_", Stage, Severity)) 
  
  pop_plot_new <- rbind(pops, pop_plot_new)
  
}

pop_plot_new_coh1 <- pop_plot_new %>%
  filter(general_cluster %in% georg) %>% 
  mutate(Stage_Severity = factor(Stage_Severity,
                                 levels = c("CVID_baseline", "CVID_progressionMild", "CVID_convalescenceMild",
                                            "non-CVID_baseline", "non-CVID_progressionMild", "non-CVID_convalescenceMild",
                                            "non-CVID_progressionSevere", "non-CVID_convalescenceSevere" )) ) %>% 
  mutate(COH = "COH1")


pop_plot_new_coh1$general_cluster <- factor(pop_plot_new_coh1$general_cluster, levels = rev(georg))

#COH2 PERCENTAGES
load("~/Escritorio/datos_covid/newdata/integration_samples/cohort2/filter_10000/tcell/tcell_coh2_annot.RData")

#COH1 AND COH2
#Percentages T cells per cluster 
tcell_metadata_2 <- tcell_2@meta.data 

tcell_metadata_2_new <- tcell_metadata_2 %>% 
  separate(Stage, into = c("Stage", "nada"), sep = "_") %>% 
  mutate(StageSeverity = paste0(Stage, "_", CovidSeverity),
         PID_StageSeverity = ifelse(PID == "CVID", paste("CVID", StageSeverity, sep = "_"), paste("non-CVID", StageSeverity, sep = "_")), 
         Patient_Stage = paste0(patient_id, "_", PID_StageSeverity)
  )

#General general_cluster
pop_tb <- table(tcell_metadata_2_new$Patient_Stage, tcell_metadata_2_new$general_cluster)
pop_tb <- as.data.frame.matrix(pop_tb)
pop_tb$sample <- rownames(pop_tb)
pop_plot <- tidyr::gather(pop_tb, key = general_cluster, value = counts, -sample)

samples <- unique(pop_plot$sample)
pop_plot_new_coh2 <- data.frame()
for (patient in samples) {
  
  pops <- pop_plot %>%
    filter(sample == patient) %>% 
    group_by(sample,general_cluster) %>%
    summarise(counts = sum(counts)) %>%
    mutate(percent= counts*100/sum(counts)) %>%
    mutate(percent.text = paste0(sprintf("%.2f",percent),"%")) %>%
    mutate(position = 100-(cumsum(percent) - (0.5*percent))) %>% 
    mutate(PID = ifelse(str_detect(sample, "non-CVID"), "non-CVID", "CVID") ) %>% 
    mutate(Severity = ifelse(str_detect(sample, "baseline"), "baseline", 
                             (ifelse(str_detect(sample, "Mild"), "Mild", "Severe")) ) ) %>% 
    mutate(Stage = ifelse(str_detect(sample, "progression"), "progression", 
                          (ifelse(str_detect(sample, "baseline"), "", "convalescence")) ) )  %>% 
    mutate(Stage_Severity = paste0(PID, "_", Stage, Severity)) 
  
  pop_plot_new_coh2 <- rbind(pops, pop_plot_new_coh2)
  
}


pop_plot_new_coh2 <- pop_plot_new_coh2 %>%
  filter(general_cluster %in% georg,
         PID %in% "non-CVID") %>% 
  mutate(Stage_Severity = factor(Stage_Severity,
                                 levels = c("CVID_baseline", "CVID_progressionMild", "CVID_convalescenceMild",
                                            "non-CVID_baseline", "non-CVID_progressionMild", "non-CVID_convalescenceMild",
                                            "non-CVID_progressionSevere", "non-CVID_convalescenceSevere"
                                 ))) %>% 
  mutate(COH = "COH2")

pop_plot_new <- rbind(pop_plot_new_coh1, pop_plot_new_coh2) %>% 
  mutate(coh_stage = paste0(COH, ":", Stage_Severity))

pop_plot_new <- pop_plot_new %>%
  mutate(coh_stage = factor(coh_stage,
                            levels = c("COH1:CVID_baseline", "COH1:CVID_progressionMild", "COH1:CVID_convalescenceMild",
                                       "COH1:non-CVID_baseline", "COH1:non-CVID_progressionMild", "COH1:non-CVID_convalescenceMild",
                                       "COH1:non-CVID_progressionSevere", "COH1:non-CVID_convalescenceSevere",
                                       
                                       "COH2:non-CVID_baseline", "COH2:non-CVID_progressionMild", "COH2:non-CVID_convalescenceMild",
                                       "COH2:non-CVID_progressionSevere", "COH2:non-CVID_convalescenceSevere"
                                       
                            ))
  )


# Supp Fig 4c ------------
pal <- c(RColorBrewer::brewer.pal(12, "Set3")[1:8], RColorBrewer::brewer.pal(12, "Set3")[4:8])
ggplot(pop_plot_new, aes(x=coh_stage, y = percent), fill = Stage_Severity, color = Stage_Severity) +
  geom_boxplot(aes(fill = Stage_Severity), outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 0.5) +
  theme_pubr() + border() + 
  # ylim(0,15) +
  facet_wrap(~general_cluster, scales = "free_y", ncol = 3) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6)) +
  scale_fill_manual(values = pal) +
  ylab(label = "Percentage") + xlab(label = "Stage") + labs(title=" ") +
  geom_vline(xintercept = c(8.5), linetype="solid") +
  geom_vline(xintercept = c(3.5), linetype="dashed") +
  NoLegend() +
  geom_signif(test = "wilcox.test", comparisons = list(c("COH1:CVID_baseline", "COH1:CVID_progressionMild"),
                                                       c("COH1:CVID_baseline", "COH1:CVID_convalescenceMild"),
                                                       c("COH1:non-CVID_baseline", "COH1:non-CVID_progressionMild"),
                                                       c("COH1:non-CVID_baseline", "COH1:non-CVID_convalescenceMild"),
                                                       c("COH1:non-CVID_baseline", "COH1:non-CVID_progressionSevere"),
                                                       c("COH1:non-CVID_baseline", "COH1:non-CVID_convalescenceSevere"),
                                                       c("COH2:non-CVID_baseline", "COH2:non-CVID_progressionMild"),
                                                       c("COH2:non-CVID_baseline", "COH2:non-CVID_convalescenceMild"),
                                                       c("COH2:non-CVID_baseline", "COH2:non-CVID_progressionSevere"),
                                                       c("COH2:non-CVID_baseline", "COH2:non-CVID_convalescenceSevere")),
              map_signif_level = T,
              textsize = 4, 
              y_position = c(15, 18, 20, 15, 18, 20, 22, 15, 18, 20, 22),
  )


