# B CELL RECLUSTERING COHORT 1 --------
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)

load("../COVID_REANNOT.RData")

Idents(covid_cvid) <- "integrated_snn_res.2"

bcell_1 <- subset(covid_cvid, idents = c(5, 19, 47, 48, 44))

DefaultAssay(bcell_1) <- "RNA"
bcell.list <- SplitObject(bcell_1, split.by = "sample_names")
bcell.feature <- SelectIntegrationFeatures(object.list = bcell.list)
bcell.anchors <- FindIntegrationAnchors(object.list = bcell.list, anchor.features = bcell.feature)
bcell_1 <- IntegrateData(anchorset = bcell.anchors)

save(bcell_1, file = "bcell_cohort1_befscale.RData")

bcell_1 <- ScaleData(bcell_1)
bcell_1 <- RunPCA(bcell_1)
bcell_1 <- RunUMAP(bcell_1, reduction = "pca", dims = 1:30)
bcell_1 <- FindNeighbors(bcell_1, reduction = "pca", dims = 1:30)
bcell_1 <- FindClusters(bcell_1, resolution = 0.6, graph.name = "integrated_snn")

# Annotation
Idents(bcell_1) <- "integrated_snn_res.0.6"
current.cluster <- levels(bcell_1)
current.cluster

new.cluster <- c("Naive_BC", #cl 0
                 "US_Memory_BC", #cl 1
                 "S_Memory_BC", #cl 2
                 "Doublets B:Myeloids", #cl 3
                 "Transitional_BC", #cl 4 
                 "CD21low_BC", #cl 5
                 "Plasma cells", #cl 6
                 "Doublets B:Myeloids", #cl 7
                 "Doublets B:T cells", #cl 8
                 "Plasma cells" #cl9
)

names(x = new.cluster) <- levels(x = bcell_1)
bcell_1 <- RenameIdents(object = bcell_1, new.cluster)
bcell_1$general_cluster <- Idents(object = bcell_1)

## Figure 2a
DimPlot(bcell_1, label = F, group.by = "general_cluster", pt.size = 0.5, cols = c("Naive_BC" ="#DC9F2D",
                                                                                  "US_Memory_BC" = "#D58EC4", 
                                                                                  "S_Memory_BC" = "#E3A39B", 
                                                                                  "Doublets B:Myeloids" = "#5CB32E", 
                                                                                  "Transitional_BC" = "#F6D832", 
                                                                                  "CD21low_BC" = "#969DCA", 
                                                                                  "Plasma cells" = "#66C2A5", 
                                                                                  "Doublets B:T cells" = "#B7D84C") ) +
  labs(title = "B cells (n = 14,967)")

## Supp figure 2a
bcell_1@meta.data <- bcell_1@meta.data %>% 
  mutate(PID = ifelse(StudyName == "Javi", "CVID", "non-CVID") ) %>% 
  select(-StudyName)
min <- table(bcell_1$PID) %>%  min()

set.seed(123)
DefaultAssay(bcell_1) <- "RNA"
#CVID
sample1_rna <- subset(bcell_1, subset = PID == "CVID")
metadata1 <- as.data.frame(sample1_rna@meta.data)
nuevo_m1 <- row.names(sample_n(metadata1, min))
cvid_bcell_subsample <- subset(bcell_1, cells = nuevo_m1)

noncvid_bcell <- subset(bcell_1, subset = PID == "non-CVID") 
#non-CVID mild
sample2_rna <- subset(noncvid_bcell, subset = CovidSeverity == c("Severe"), invert = T)
metadata2 <- as.data.frame(sample2_rna@meta.data)
nuevo_m2 <- row.names(sample_n(metadata2, min))
noncvid_mild_bcell_subsample <- subset(bcell_1, cells = nuevo_m2)
#non-CVID severe
sample3_rna <- subset(noncvid_bcell, subset = CovidSeverity == c("Mild"), invert = T)
metadata3 <- as.data.frame(sample3_rna@meta.data)
nuevo_m3 <- row.names(sample_n(metadata3, min))
noncvid_severe_bcell_subsample <- subset(bcell_1, cells = nuevo_m3)

new.col.bcell <- c("#DC9F2D",#Naive_BC
                   "#D58EC4", # US_Memory_BC
                   "#E3A39B", # S_Memory_BC
                   "#5CB32E", # Doublets B:bcelloids
                   "#F6D832", # Transitional_BC
                   "#969DCA", # CD21low_BC
                   "#66C2A5", # Plasma cells
                   "#B7D84C" # Doublets B:T cells
)

p1 <- DimPlot(cvid_bcell_subsample, group.by = "general_cluster", pt.size = 0.25, cols = new.col.bcell, raster = F) + labs(title = "CVID") + NoLegend()
p2 <- DimPlot(noncvid_mild_bcell_subsample, group.by = "general_cluster", pt.size = 0.25, cols = new.col.bcell, raster = F) + labs(title = "non-CVID mild") + NoLegend()
p3 <- DimPlot(noncvid_severe_bcell_subsample, group.by = "general_cluster", pt.size = 0.25, cols = new.col.bcell, raster = F) + labs(title = "non-CVID severe")
p1+p2+p3

# PHENOGRAPH
Idents(bcell_1) <- "general_cluster"
markers.bcell <- FindAllMarkers(bcell_1, assay = "RNA", logfc.threshold = 0.25, min.pct = 0.3)

top.genes.all <- data.frame()
bcell_1@active.ident <- factor(bcell_1@active.ident, 
                               levels = rev(c("Naive_BC", #cl 0
                                              "Transitional_BC", #cl 4 
                                              "US_Memory_BC", #cl 1
                                              "S_Memory_BC", #cl 3
                                              "CD21low_BC", #cl 5
                                              "Plasma cells", #cl 6
                                              "Doublets B:T cells", #cl 7
                                              "Doublets B:Myeloids" #cl 2
                               )))

clusters <- c("Naive_BC", #cl 0
                   "Transitional_BC", #cl 4 
                   "US_Memory_BC", #cl 1
                   "S_Memory_BC", #cl 3
                   "CD21low_BC", #cl 5
                   "Plasma cells", #cl 6
                   "Doublets B:T cells", #cl 7
                   "Doublets B:Myeloids" #cl 2
)

for (i in clusters) {
  top.genes <- markers.bcell %>% 
    mutate(sign = ifelse(p_val_adj < 0.05 & cluster == i & pct.1 > 0.4 & avg_log2FC > 0, paste("sign", i, sep = "_"), "others")) %>% 
    dplyr::filter(sign != "others")
  top.genes.all <- rbind(top.genes.all, top.genes)
}

gene.markers <- (table(top.genes.all$gene)) %>% as.data.frame() %>% 
  dplyr::filter(Freq == 1) %>% 
  dplyr::filter(!grepl("^MT-", Var1))

top.genes.common <- top.genes.all %>% 
  dplyr::filter(gene %in% gene.markers$Var1) %>% 
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) 

## Figure 2b 
DefaultAssay(bcell_1) <- "RNA"
DotPlot(bcell_1, features = top.genes.common$gene, scale.max = 100, scale.min = 0, dot.scale = 5)  +
  scale_color_gradientn(colors = c("skyblue3", "white", "#770000")) +
  theme_pubr(legend = "right") + border()+ rotate_x_text(45) 

# B CELL RECLUSTERING COHORT 2 -----------
load("../coh2_covid_samples_annot.RData")

Idents(covid_2) <- "general_cluster"
bcell <- subset(covid_2, idents = "B cells")

#Filter out all samples with <100 cells
bcell <- subset(bcell, idents = c("progression_70_progression_mild",
                                  "progression_39_progression_severe",
                                  "progression_29_progression_mild", "progression_325_progression_mild",
                                  "progression_21_progression_severe", "progression_12_progression_severe",
                                  "progression_108_progression_severe", "convalescence_70_convalescence_mild",
                                  "convalescence_50029_convalescence_mild",
                                  "convalescence_39_convalescence_severe", "convalescence_325_convalescence_mild",
                                  "convalescence_29_convalescence_mild", "convalescence_136_convalescence_severe", "convalescence_108_convalescence_severe",
                                  "baseline_91_baseline", "baseline_442_baseline", "baseline_358_baseline", "baseline_385_baseline", #MENOS DE 50 CELULAS
                                  "baseline_79_baseline" #MENOS DE 100 CELULAS
), invert = T)

DefaultAssay(bcell) <- "RNA"
#rerun clustering
bcell.list <- SplitObject(bcell, split.by = "sample_names")
bcell.feature <- SelectIntegrationFeatures(object.list = bcell.list)
bcell.anchors <- FindIntegrationAnchors(object.list = bcell.list, anchor.features = bcell.feature)
bcell <- IntegrateData(anchorset = bcell.anchors)

bcell <- ScaleData(bcell)
bcell <- RunPCA(bcell)
bcell <- harmony::RunHarmony(bcell, "sample_names")
bcell <- RunUMAP(bcell, reduction = "pca", dims = 1:30)
bcell <- FindNeighbors(bcell, reduction = "pca", dims = 1:30)
bcell_2 <- FindClusters(bcell, resolution = 0.5, graph.name = "integrated_snn")


bcell_2 <- FindClusters(bcell_2, resolution = 0.2, graph.name = "integrated_snn")
Idents(bcell_2) <- "seurat_clusters"
current.cluster <- levels(bcell_2)
current.cluster

new.cluster <- c("Naive_BC", #cl 0
                 "US_Memory_BC", #cl 1
                 "Transitional_BC", #cl 2
                 "Plasma cells", #cl 3
                 "Doublets B:T cells", #cl 4 
                 "S_Memory_BC", #cl 5
                 "CD21low_BC", #cl 6
                 "Doublets B:Myeloids", #cl 7
                 "US_Memory_BC" #cl 8
)

names(x = new.cluster) <- levels(x = bcell_2)
bcell_2 <- RenameIdents(object = bcell_2, new.cluster)
bcell_2$general_cluster <- Idents(object = bcell_2)
table(bcell_2$general_cluster)

#Supp Fig 2 b
new.col.bcell <- c("#DC9F2D",#Naive_BC
                   "#D58EC4", # US_Memory_BC
                   "#F6D832", # Transitional_BC
                   "#66C2A5", # Plasma cells
                   "#B7D84C", # Doublets B:T cells
                   "#E3A39B", # S_Memory_BC
                   "#969DCA", # CD21low_BC
                   "#5CB32E" # Doublets B:Myeloids
)

DimPlot(bcell_2, label = F, group.by = "general_cluster", pt.size = 0.25, cols = new.col.bcell) +
  labs(title = "B cells")

#Supp Fig 2 c
DefaultAssay(bcell_2) <- "RNA"
DotPlot(bcell_2, features = top.genes.common$gene, scale.max = 100, scale.min = 0, dot.scale = 5)  +
  scale_color_gradientn(colors = c("skyblue3", "white", "#770000")) +
  theme_pubr(legend = "right") + border()+ rotate_x_text(45) 

#--------------------
# Supp Fig 2 b 

## SCORES HEATMAP
# NAIVE B CELLLS
naive <- subset(bcell_1, subset = general_cluster == "Naive_BC", invert = F)

#GO:0002377_immunoglobulin production CLUSTER 3, 13
IMMUNOGLOBULIN <- read.table("markers/bcell_clusters/GOBP_IMMUNOGLOBULIN_PRODUCTION.v2023.2.Hs.grp", header = T)
naive <- AddModuleScore(naive, features = list(unique(IMMUNOGLOBULIN$GOBP_IMMUNOGLOBULIN_PRODUCTION)),
                        name = c("cl3_13_IFN_"), search = T)

#GO:0034341_response to interferon-gamma CLUSTER 6
IFN <- read.table("markers/bcell_clusters/GOBP_RESPONSE_TO_TYPE_I_INTERFERON.v2023.1.Hs.grp", header = T)
naive <- AddModuleScore(naive, features = list(unique(IFN$GOBP_RESPONSE_TO_TYPE_I_INTERFERON)),
                        name = c("cl6_IFN_"), search = T)

#GO:0035455 response to interferon-alpha CLUSTER 8
ALPHA <- read.table("markers/bcell_clusters/GOBP_RESPONSE_TO_INTERFERON_ALPHA.v2023.2.Hs.grp", header = T)
naive <- AddModuleScore(naive, features = list(unique(ALPHA$GOBP_RESPONSE_TO_INTERFERON_ALPHA)), name = "cl8_alpha")

#GO:0032606_type I interferon production CLUSTER 10
PRODUCTION <- read.table("markers/bcell_clusters/GOBP_TYPE_I_INTERFERON_PRODUCTION.v2023.2.Hs.grp", header = T)
naive <- AddModuleScore(naive, features = list(unique(PRODUCTION$GOBP_TYPE_I_INTERFERON_PRODUCTION)), name = "cl10_production")

#GO:00301831_B cell differentiation CLUSTER 13
DIFFERENTIATION <- read.table("markers/bcell_clusters/GOBP_B_CELL_DIFFERENTIATION.v2023.2.Hs.grp", header = T)
naive <- AddModuleScore(naive, features = list(unique(DIFFERENTIATION$GOBP_B_CELL_DIFFERENTIATION)), name = "cl13_cl15_differentiation")

#GO:00096151_response to virus
virus <- read.table("markers/bcell_clusters/GOBP_RESPONSE_TO_VIRUS.v2023.2.Hs.grp", header = T)
naive <- AddModuleScore(naive, features = list(unique(virus$GOBP_RESPONSE_TO_VIRUS)), name = "virus_")

#GO:0098586_cellular response to virus
cellular <- read.table("markers/bcell_clusters/GO:0098586.grp", header = T)
naive <- AddModuleScore(naive, features = list(unique(cellular$GENE)), name = "cellular_")

plot <- naive@meta.data %>%
  mutate(StudyName = ifelse(StudyName == "China", "nonCVID", "CVID")) %>% 
  mutate(sample = paste0(Stage, "_", CovidSeverity, "-", StudyName, ":", general_cluster)) %>%
  mutate(patient = paste0(sample, ":", StageCVID)) %>%
  select(general_cluster, sample, StageCVID, patient, cl6_IFN_1, cl8_alpha1, cl10_production1, cl13_cl15_differentiation1, cl3_13_IFN_1, virus_1, cellular_1)


#Plot por los clusters de los scores
## Summary plot intento
#Mean per stage  
cl3_13_IFN <-  plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl3_13_IFN_1), sd=sd(cl3_13_IFN_1)) %>% 
  mutate(cluster = "cl3_13_IFN_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))


pop_plot_cl6 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl6_IFN_1), sd=sd(cl6_IFN_1)) %>% 
  mutate(cluster = "cl6")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

pop_plot_cl8 <- plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl8_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl8_alpha1), sd=sd(cl8_alpha1)) %>% 
  mutate(cluster = "cl8")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

pop_plot_cl10 <- plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl10_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl10_production1), sd=sd(cl10_production1)) %>% 
  mutate(cluster = "cl10")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

pop_plot_cl13_15 <-plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl13_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl13_cl15_differentiation1), sd=sd(cl13_cl15_differentiation1)) %>% 
  mutate(cluster = "cl13")  %>% 
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

cellular_1 <-plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl13_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cellular_1), sd=sd(cellular_1)) %>% 
  mutate(cluster = "cellular_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))
# 
all_pop_plot_median <- rbind ( pop_plot_cl6 , pop_plot_cl8, pop_plot_cl10, pop_plot_cl13_15, virus_1, cellular_1 )

# all_pop_plot_median <- all_pop_plot_median %>% 
#   separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
#   separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
#   mutate(group = paste0(PID_stage, "_", cluster))

all_pop_plot_median$stage <- factor(all_pop_plot_median$stage, levels = c("baseline", "progression", "convalescence"))
writexl::write_xlsx(all_pop_plot_median, path = "new_images_def/NEW/bcell_naive_GO_nocorrect.xlsx")

#PLOTS pop_plot_cl6 pop_plot_cl8 pop_plot_cl10 pop_plot_cl13_15
cellular_1$PID_stage <- factor(cellular_1$PID_stage, levels = c("CVID:baseline_Control", "CVID:progression_Mild", "CVID:convalescence_Mild",
                                                                "nonCVID:baseline_Control", "nonCVID:progression_Mild", "nonCVID:convalescence_Mild",
                                                                "nonCVID:progression_Severe", "nonCVID:convalescence_Severe"))

ggplot(cellular_1, aes(x=PID_stage, y=cluster, fill = mean) ) +
  geom_tile(aes(fill=mean), color = "black")  +
  facet_grid(~cluster, scales = "free") +
  #scale_fill_viridis() +GO:00096151_response to virus
  scale_fill_gradientn(colours = c("#4267AC", "#1982C4", "#8AC926", "#FFCA3A", "#FF924C", "#FF595E")) +
  #                        c("white", "darkred") ) +
  ggpubr::theme_pubr(legend = "right", border = T, x.text.angle = 45) 


# US MEMORY
memory <- subset(bcell_1, subset = general_cluster == "US_Memory_BC", invert = F)

#clusters scores
#GO:00343401_response to type I interferon CLUSTER 1, 8 
IFN <- read.table("markers/bcell_clusters/GOBP_RESPONSE_TO_TYPE_I_INTERFERON.v2023.1.Hs.grp", header = T)
memory <- AddModuleScore(memory, features = list(unique(IFN$GOBP_RESPONSE_TO_TYPE_I_INTERFERON)), name = "cl1_cl8_response")


#GO:0002377_immunoglobulin production CLUSTER 5, 9
IMMUNOGLOBULIN <- read.table("markers/bcell_clusters/GOBP_IMMUNOGLOBULIN_PRODUCTION.v2023.2.Hs.grp", header = T)
memory <- AddModuleScore(memory, features = list(unique(IMMUNOGLOBULIN$GOBP_IMMUNOGLOBULIN_PRODUCTION)),
                         name = c("cl5_9_IFN_"), search = T)


#GO:0034340_response to type I interferon CLUSTER 11
BETA <- read.table("markers/bcell_clusters/GOBP_INTERFERON_BETA_PRODUCTION.v2023.2.Hs.grp", header = T)
memory <- AddModuleScore(memory, features = list(unique(BETA$GOBP_INTERFERON_BETA_PRODUCTION)), name = "cl11_beta")

#GO:00301831_B cell differentiation CLUSTER 15 
DIFFERENTIATION <- read.table("markers/bcell_clusters/GOBP_B_CELL_DIFFERENTIATION.v2023.2.Hs.grp", header = T)
memory <- AddModuleScore(memory, features = list(unique(DIFFERENTIATION$GOBP_B_CELL_DIFFERENTIATION)), name = "cl15_differentiation")


#GO:00096151_response to virus
virus <- read.table("markers/bcell_clusters/GOBP_RESPONSE_TO_VIRUS.v2023.2.Hs.grp", header = T)
memory <- AddModuleScore(memory, features = list(unique(virus$GOBP_RESPONSE_TO_VIRUS)), name = "virus_")

score <- memory@meta.data %>% 
  select(cl1_cl8_response1, cl11_beta1, cl15_differentiation1, cl5_9_IFN_1) 
plot <- memory@meta.data %>%
  mutate(StudyName = ifelse(StudyName == "China", "nonCVID", "CVID")) %>% 
  mutate(sample = paste0(Stage, "_", CovidSeverity, "-", StudyName, ":", general_cluster)) %>%
  mutate(patient = paste0(sample, ":", StageCVID)) %>%
  select(general_cluster, sample, StageCVID, patient, cl1_cl8_response1, cl11_beta1, cl15_differentiation1, cl5_9_IFN_1, virus_1)


#Plot por los clusters de los scores
#Mean per stage and patient
cl1_1 <- plot %>% #plyr::ddply(plot, "patient", summarise, media=mean(cl1_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl1_cl8_response1), sd=sd(cl1_cl8_response1)) %>% 
  mutate(cluster = "cl1") %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))


cl5_9_IFN_1 <- plot %>% #plyr::ddply(plot, "patient", summarise, media=mean(cl1_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl5_9_IFN_1), sd=sd(cl5_9_IFN_1)) %>% 
  mutate(cluster = "cl5_9_IFN_1") %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl11_1 <- plot %>% # plyr::ddply(plot, "patient", summarise, media=mean(cl11_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl11_beta1), sd=sd(cl11_beta1)) %>% 
  mutate(cluster = "cl11") %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl15_1 <- plot %>% #plyr::ddply(plot, "patient", summarise, media=mean(cl15_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl15_differentiation1), sd=sd(cl15_differentiation1)) %>% 
  mutate(cluster = "cl15") %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

virus_1 <- plot %>% #plyr::ddply(plot, "patient", summarise, media=mean(cl15_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(virus_1), sd=sd(virus_1)) %>% 
  mutate(cluster = "cl15") %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

# 
all_pop_plot_median <- rbind ( cl1_1, cl11_1, cl15_1, cl5_9_IFN_1, virus_1 )


all_pop_plot_median$stage <- factor(all_pop_plot_median$stage, levels = c("baseline", "progression", "convalescence"))

writexl::write_xlsx(all_pop_plot_median, path = "new_images_def/NEW/bcell_US_GO_mean_nocorreciton.xlsx")

#PLOTS cl1_1 cl11_1 cl15_1 cl5_9_IFN_1, cl15_1
virus_1$PID_stage <- factor(virus_1$PID_stage, levels = c("CVID:baseline_Control", "CVID:progression_Mild", "CVID:convalescence_Mild",
                                                          "nonCVID:baseline_Control", "nonCVID:progression_Mild", "nonCVID:convalescence_Mild",
                                                          "nonCVID:progression_Severe", "nonCVID:convalescence_Severe"))

ggplot(virus_1, aes(x=PID_stage, y=cluster, fill = mean) ) +
  geom_tile(aes(fill=mean), color = "black")  +
  facet_grid(~cluster, scales = "free") +
  #scale_fill_viridis() +
  scale_fill_gradientn(colours = c("#4267AC", "#1982C4", "#8AC926", "#FFCA3A", "#FF924C", "#FF595E")) +
  #                        c("white", "darkred") ) +
  ggpubr::theme_pubr(legend = "right", border = T, x.text.angle = 45) 

