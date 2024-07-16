library(tidyverse)
library(dplyr)
library(Seurat)
library(ggpubr)

# Integration cohort 1

load("/ijc/LABS/BALLESTAR/DATA/DATA_CELIA/DATOS_COVID/raw/covidcvid_cohort1.RData")

Idents(object = covidcvid) <- "sample_names" 
all.list <-  SplitObject(covidcvid, split.by = "sample_names")

for (i in 1:length(all.list)) {
  all.list[[i]] <- NormalizeData(all.list[[i]], verbose = FALSE)
  all.list[[i]] <- FindVariableFeatures(all.list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
}

features <- SelectIntegrationFeatures(object.list = all.list)
anchors <- FindIntegrationAnchors(object.list = all.list, reference = 1)
covid <- IntegrateData(anchorset = anchors)
gc()

save(covid, file="cohort1_integrated_samplenames.RData")

covid <- ScaleData(covid, verbose = FALSE)
covid <- RunPCA(covid, npcs = 50, #50 PCA by default
                features = features, verbose = FALSE)
covid_cvid <- RunUMAP(covid_cvid, dims = 1:50)
covid_cvid <- FindNeighbors(covid_cvid, reduction = "pca", dims = 1:50)

covid_cvid <- FindClusters(covid_cvid, resolution = c(1,2), graph.name = "integrated_snn")

# Add metadata
covid_cvid@meta.data <- covid_cvid@meta.data %>% 
  mutate(StageSeverity = paste0(Stage, "_", CovidSeverity),
         PID = ifelse(StudyName == "Javi", "CVID", "non-CVID"),
         StageStudy = paste0(Stage, "_", CovidSeverity, "_", PID),
         PID_StageSeverity = ifelse(PID == "CVID", paste("CVID", StageSeverity, sep = "_"), paste("non-CVID", StageSeverity, sep = "_")), 
         Patient_Stage = paste0(patientID, "_", PID_StageSeverity)
  )

# QC 
DefaultAssay(covid_cvid) <- "RNA"
Idents(covid_cvid) <- "StudyName"

covid_cvid[["percent.mt"]] <- PercentageFeatureSet(covid_cvid, pattern = "^MT-") 
plot1 <- FeatureScatter(covid_cvid, feature1 = "nFeature_RNA", feature2 = "percent.mt", group.by = "StudyName")
plot2 <- FeatureScatter(covid_cvid, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "StudyName")
plot1 + plot2

covid_cvid <- subset(covid_cvid, subset = nFeature_RNA > 500 & percent.mt < 10 & nFeature_RNA < 5000 & nCount_RNA < 7500)

## SUBSAMBLING CVID and nonCVID
set.seed(123)
#CVID
sample1_rna <- subset(covid_cvid, subset = StudyName == "China")
metadata1 <- as.data.frame(sample1_rna@meta.data)
nuevo_m1 <- row.names(sample_n(metadata1, 70643))
noncvid_subsample <- subset(covid_cvid, cells = nuevo_m1)

#non-CVID 
sample2_rna <- subset(covid_cvid, subset = StudyName == "Javi")
metadata2 <- as.data.frame(sample2_rna@meta.data)
nuevo_m2 <- row.names(sample_n(metadata2, 70643))
cvid_subsample <- subset(covid_cvid, cells = nuevo_m2)

all_metadata <- c(nuevo_m1, nuevo_m2) 
subset <- subset(covid_cvid, cells = all_metadata)
subset 

subset@meta.data$integrated_snn_res.1 <- droplevels(subset@meta.data$integrated_snn_res.1)
p1 <- DimPlot(cvid_subsample, group.by = "StudyName", raster = F, cols = "#8FB8ED") + labs(title = "nonCVID") + NoLegend()
p2 <- DimPlot(noncvid_subsample, group.by = "StudyName", raster = F, cols =  "#ae2012") + labs(title = "CVID") + NoLegend()
p1+p2

## Supplementary Fig 1 a, b, c, d, e
DimPlot(covid_cvid, group.by = "Sex", pt.size = 0.3)
DimPlot(covid_cvid, group.by = "PMID", pt.size = 0.3, order = c("PMID: 32759967", "PMID: 32788748", "PMID: 32377375", "CVID"), 
        cols = c("CVID" = "darkseagreen4", "PMID: 32759967" = "darksalmon", "PMID: 32788748" = "cornflowerblue", "PMID: 32377375" = "khaki3"))
DimPlot(covid_cvid, group.by = "StageSeverity", pt.size = 0.3)
DimPlot(covid_cvid, group.by = "Age", pt.size = 0.3)

# Annotation of general clusters
DefaultAssay(covid_cvid) <- "RNA"
Idents(covid_cvid) <- "seurat_clusters"
current.cluster <- levels(covid_cvid)
current.cluster

new.cluster <- c("CD8_memory", #cl 0
                 "CD4_memory", #cl 1
                 "Classical Mono", #cl 2
                 "NK CD56dim", #cl 3
                 "CD4_naive", #cl 4
                 "Classical Mono", #cl 5
                 "Classical Mono", #cl 6
                 "Naive BC", #cl 7
                 "Doublets", #cl 8
                 "Non-classical Mono", #cl 9
                 "Classical Mono", #cl 10
                 "CD8_memory", #cl 11
                 "CD8_memory", #cl 12
                 "CD8_NKT_like", #cl 13
                 "Memory BC", #cl 14
                 "Proliferating TC", #cl 15
                 "CD8_naive", #cl 16
                 "CD8_memory", #cl 17
                 "gdTC", #cl 18 
                 "CD8_memory", #cl 19
                 "Platelets", #cl 20
                 "Treg", #cl 21
                 "DC", #cl 22
                 "CD8_memory", #cl 23
                 "Classical Mono", #cl 24
                 "CD8_memory", #cl 25
                 "CD8_memory", #cl 26
                 "CD4_CTL_P3", #cl 27
                 "Plasma cells", #cl 28
                 "MAIT", #cl 29
                 "CD8_memory", #cl 30
                 "NK CD56bright", #cl 31
                 "Doublets", #cl 32
                 "Doublets", #cl 33
                 "Doublets", #cl 34
                 "hSC" #cl 35
)
names(x = new.cluster) <- levels(x = covid_cvid)
covid_cvid <- RenameIdents(object = covid_cvid, new.cluster)
covid_cvid$clusters <- Idents(object = covid_cvid)
table(covid_cvid$clusters)

# Annotate pDC with the metadata of myeloids reclustering
metadata_myel <- myel_1@meta.data %>% 
  select(general_cluster) %>% 
  mutate(myel_annotation = general_cluster) %>% 
  tibble::rownames_to_column("barcodes")

metadata <- covid_cvid@meta.data %>% 
  mutate(clusters = as.character(clusters)) %>% 
  tibble::rownames_to_column("barcodes")

metadata_new <- left_join(metadata, metadata_myel, "barcodes")
metadata_new <- metadata_new %>% 
  mutate(myel_annot = ifelse(is.na(myel_annotation) == T, "NA", paste(myel_annotation))) %>% 
  mutate(general_cluster = ifelse(myel_annot == "pDC", "pDC", clusters)) %>% 
  tibble::column_to_rownames("barcodes")
table(metadata_new$general_cluster)

covid_cvid <- AddMetaData(covid_cvid, metadata_new$general_cluster, col.name = "CLUSTERS")

save(covid_cvid, file = "../../newdata/integration_samples/cohort1/COVID_REANNOT.RData")

## Figure 1 b
DimPlot(covid_cvid, raster = F, group.by = "CLUSTERS",
        cols = c("pDC" = "#03045e",
                 "Naive BC" = "#faa275",
                 "Classical Mono" = "#61a1d7",
                 "NK CD56dim" = "#f4d06f",
                 "hSC" = "#ffd6ff",
                 "NK CD56bright" = "#90be6d",
                 "Plasma cells" = "#43aa8b",
                 "DC" = "#bbdef0",
                 "Platelets" = "#219ebc",
                 "Proliferating TC" = "#e76f51", 
                 "Memory BC" = "#ffe94e",
                 "Non-classical Mono" = "#2660a4",
                 "CD4_naive" = "#E2CCF5",
                 "CD4_memory" = "#B884D7",
                 "CD8_naive" = "#8D9933", 
                 "CD8_memory" = "#83c5be",
                 "CD8_NKT_like" = "#BE8F0E",
                 "gdTC" = "#3e5c76",
                 "MAIT" = "#8F0000",
                 "Treg" = "#d00000",
                 "CD4_CTL_P3" = "#A0B6BB", 
                 "Doublets" = "gray89"
        ))

## Figure 1 c
Idents(covid_cvid) <- "CLUSTERS"
DefaultAssay(covid_cvid) <- "RNA"

covid_cvid@active.ident <- factor(covid_cvid@active.ident, 
                                  levels = (c("Treg", "CD4_naive", "CD4_memory", "CD4_CTL_P3", 
                                              "CD8_naive", "CD8_memory", "CD8_NKT_like", "gdTC", "MAIT", "Proliferating TC", 
                                              "NK CD56bright", "NK CD56dim",
                                              "Naive BC", "Memory BC", "Plasma cells", 
                                              "Classical Mono", "Non-classical Mono", "DC", "pDC",
                                              "Platelets", "hSC", "Doublets"
                                  )))

covid_cvid <- subset(covid_cvid, idents = "Doublets", invert = T)

pal<- colorRampPalette(brewer.pal(12, "Spectral"))(19)
VlnPlot(covid_cvid, features = c("CD3D", "FOXP3", "CD4", "CD8A", "FOXP3", "TRDV2", "SLC4A10", "MKI67", 
                                 "NCAM1", "FCGR3A", "CD19", "CCR7", "IGHM", "JCHAIN", 
                                 "CD33", "FCER1A", "CD14", "PPBP", "CD34"), 
        assay = "RNA",
        raster = F, pt.size = 0, stack = T, flip = T, adjust = T, cols = pal) + NoLegend()

## Figure 1 d
covid_cvid <- subset(covid_cvid, idents = c("CD4_CTL_P3", "hSC"), invert = T)

covid_cvid@meta.data <- covid_cvid@meta.data %>% 
  mutate(StageSeverity = paste0(Stage, "_", CovidSeverity),
         PID = ifelse(StudyName == "Javi", "CVID", "non-CVID"))  %>% 
  mutate(PID_StageSeverity = ifelse(PID == "CVID", paste("CVID", StageSeverity, sep = "_"), paste("non-CVID", StageSeverity, sep = "_") ) )

#General clusters
pop_tb <- table(covid_cvid$PID_StageSeverity, covid_cvid$CLUSTERS)
pop_tb <- as.data.frame.matrix(pop_tb)
pop_tb$sample <- rownames(pop_tb)

pop_plot <- tidyr::gather(pop_tb, key = CLUSTERS, value = counts, -sample)
plot_order <- c("CVID_baseline_Control", "CVID_progression_Mild", "CVID_convalescence_Mild",
                "non-CVID_baseline_Control", "non-CVID_progression_Mild", "non-CVID_convalescence_Mild", "non-CVID_progression_Severe", "non-CVID_convalescence_Severe")
pop_plot$sample <- factor(pop_plot$sample, levels = rev(plot_order))

pop_plot <- pop_plot %>%
  group_by(sample,CLUSTERS) %>%
  summarise(counts = sum(counts)) %>%
  mutate(percent= counts*100/sum(counts)) %>%
  mutate(percent.text = paste0(sprintf("%.2f",percent),"%")) %>%
  mutate(position = 100-(cumsum(percent) - (0.5*percent))) %>% 
  mutate(general_cluster = factor(CLUSTERS, levels = rev(c("Treg", "CD4_naive", "CD4_memory",
                                                        "CD8_naive", "CD8_memory", "CD8_NKT_like", "gdTC", "MAIT", "Proliferating TC", 
                                                        "NK CD56bright", "NK CD56dim",
                                                        "Naive BC", "Memory BC", "Plasma cells", 
                                                        "Classical Mono", "Non-classical Mono", "DC", "pDC",
                                                        "Platelets" ))) ) %>% 
  mutate(PID = ifelse(sample %in% c("CVID_baseline_Control", "CVID_progression_Mild", "CVID_convalescence_Mild"), "CVID", "non-CVID")) 

ggplot(pop_plot, aes(x= sample, y = percent)) +
  geom_bar(aes(fill = general_cluster), stat = "identity", width = 0.8, color = "gray38") +
  scale_fill_manual(values = c("pDC" = "#03045e",
                               "Naive BC" = "#faa275",
                               "Classical Mono" = "#61a1d7",
                               "NK CD56dim" = "#f4d06f",
                               "hSC" = "#ffd6ff",
                               "NK CD56bright" = "#90be6d",
                               "Plasma cells" = "#43aa8b",
                               "DC" = "#bbdef0",
                               "Platelets" = "#219ebc",
                               "Proliferating TC" = "#e76f51", 
                               "Memory BC" = "#ffe94e",
                               "Non-classical Mono" = "#2660a4",
                               "CD4_naive" = "#E2CCF5",
                               "CD4_memory" = "#B884D7",
                               "CD8_naive" = "#8D9933", 
                               "CD8_memory" = "#83c5be",
                               "CD8_NKT_like" = "#BE8F0E",
                               "gdTC" = "#3e5c76",
                               "MAIT" = "#8F0000",
                               "Treg" = "#d00000",
                               "CD4_CTL_P3" = "#A0B6BB", 
                               "Doublets" = "gray89")) + coord_flip() +
  # facet_grid(~PID, scales = "free_x") + 
  theme_pubr(legend = "right") + border() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab(label = "Percentage") + xlab(label = "") + labs(fill="NK cells") 

