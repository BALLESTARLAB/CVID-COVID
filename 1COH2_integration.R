library(tidyverse)
library(dplyr)
library(Seurat)
library(ggpubr)

# Integration cohort 2

# Clean and prepare data
load("newdata/integration_samples/cohort2/filter_10000/coh2_covid_samples.RData")


covid_2@meta.data <- covid_2@meta.data %>% 
  mutate(StageCVID = batch) %>% 
  mutate(Days_sample_collection = "Control",
         Days_sample_collection = ifelse(batch == "P1 Prog", "27", Days_sample_collection),
         Days_sample_collection = ifelse(batch == "P1 Conv", "69", Days_sample_collection),
         Days_sample_collection = ifelse(batch == "P2 Prog", "64", Days_sample_collection),
         Days_sample_collection = ifelse(batch == "P2 Conv", "77", Days_sample_collection),
         Days_sample_collection = ifelse(batch == "P3 Prog", "34", Days_sample_collection),
         Days_sample_collection = ifelse(batch == "P3 Conv", "136", Days_sample_collection),
         Days_sample_collection = ifelse(batch == "P4 Prog", "5", Days_sample_collection),
         Days_sample_collection = ifelse(batch == "P1 Conv", "103", Days_sample_collection),
         Days_sample_collection = ifelse(batch == "P5 Prog", "15", Days_sample_collection),
         Days_sample_collection = ifelse(batch == "P1 Conv", "112", Days_sample_collection) ) %>%
  select(-c("souporcell", "batch", "full_clustering", "initial_clustering", "Fraction", "scrublet_pred", "scrublet_local_pred", "scrublet_score", "scrublet_cluster_score", "S_score", "G2M_score", "leiden_sampl_Hani","leidenres2_sampl_Hani",
            'leiden_sampl_China', 'leidenres2_sampl_China', 'leiden_sampl',
            'leidenres2_sampl', 'sampl_Hani', 'sampl', 'leiden_bcell',
            'leidenres2_bcell', 'leiden_mono', 'leidenres2_mono','sampl_China', 'leidenres2_sampl_totvi',
            'leiden_sampl_totvi', 'totvi', 'general', 'general_China',
            'general_Hani', 'prog_FDR','prog_logFC', 'prog_spFDR', 'conv_FDR', 'conv_logFC', 'conv_spFDR',
            'leiden_totvi_tcells', 'leidenres2_totvi_tcells',
            'leidenres3_totvi_tcells', 'leiden_totvi_bcells',
            'leidenres2_totvi_bcells', 'leiden_totvi_mono', 'leidenres2_totvi_mono', 'totvi_bcells',
            'totvi_mono', 'leiden_sampl_everything', 'leidenres2_sampl_everything',
            'leiden_sampl_baseline', 'neigh', 'neigh_javi', 'general_integrated','leiden_every_tcell', 'leidenres2_every_tcell',
            'leiden_every_bcell', 'leidenres2_every_bcell', 'leiden_every_mono',
            'leidenres2_every_mono', 'leiden_every64', 'leidenres2_every64',
            'leiden_every_tcell20','leidenres2_every_tcell20', 'leiden_every_tcell20unpool',
            'leidenres2_every_tcell20unpool', 'leiden_every_bcell20unpool',
            'leidenres2_every_bcell20unpool', 'leiden_every_mono20unpool',
            'leidenres2_every_mono20unpool', 'leiden_everything_unpool',
            'leidenres2_everything_unpool', 'leiden_totvi20unpool','leidenres2_totvi20unpool', 'bcell20unpool', 'tcell20unpool',
            'mono20unpool','neigh_Stage', 'neigh_Enrich',
            'neigh_Enrich_logPval','tmp','leiden_totvi20unpool_2', 'SampleID',
            'leidenres2_totvi20unpool_2', 'leiden_totvi_monocytes', 'demultiplexed', 
            'leidenres2_totvi_monocytes', 'filtered_cells', 'clonotype_class', 'TRA_V_Gene', 'TRB_V_Gene', 'TRA_J_Gene', 'TRB_J_Gene', 'TR_C_Genes', 'TR_D_Genes'
  ))

DefaultAssay(covid_2) <- "RNA"
covid.list <- SplitObject(covid_2, split.by = "sample_names")
covid.feature <- SelectIntegrationFeatures(object.list = covid.list)
covid.anchors <- FindIntegrationAnchors(object.list = covid.list, anchor.features = covid.feature)
covid <- IntegrateData(anchorset = covid.anchors)

covid <- ScaleData(covid, npcs = 50, #50 PCA by default
                   features = covid.feature, verbose = FALSE)
covid <- RunPCA(covid)
covid <- RunUMAP(covid, reduction = "pca", dims = 1:50)
covid <- FindNeighbors(covid, reduction = "pca", dims = 1:50)
covid_2 <- FindClusters(covid, resolution = 1, graph.name = "integrated_snn")

# QC
DefaultAssay(covid_2) <- "RNA"

covid_2[["percent.mt"]] <- PercentageFeatureSet(covid_2, pattern = "^MT-") 

covid_2@meta.data <- covid_2@meta.data %>% 
  mutate(percent.mt = ifelse(is.na(percent.mt) == "TRUE", 0, percent.mt))

plot1 <- FeatureScatter(covid_2, feature1 = "nFeature_RNA", feature2 = "percent.mt", group.by = "Author")
plot2 <- FeatureScatter(covid_2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Author")
plot1 + plot2

covid_2 <- subset(covid_2, subset = nFeature_RNA > 500 & percent.mt < 10 & nFeature_RNA < 5000 & nCount_RNA < 7500)

load("newdata/integration_samples/cohort2/coh2_covid_samples.RData")

# Annotation of clusters
DefaultAssay(covid_2) <- "RNA"
Idents(covid_2) <- "seurat_clusters"
current.cluster <- levels(covid_2)
current.cluster

new.cluster <- c("CD4_naive", #cl 0
                 "CD8_memory", #cl 1
                 "CD4_memory", #cl 2
                 "CD4_memory", #cl 3 
                 "NK CD56dim", #cl 4
                 "CD8_naive", #cl 5
                 "Naive BC", #cl 6
                 "CD8_NKT_like", #cl 7
                 "Classical Mono", #cl 8
                 "Classical Mono", #cl 9
                 "CD8_memory", #cl 10
                 "Doublets", #cl 11
                 "Doublets", #cl 12
                 "MAIT", #cl 13
                 "Memory BC", #cl 14
                 "CD8_memory", #cl 15
                 "Non-classical Mono", #cl 16
                 "Doublets", #cl 17 Naive BC
                 "Treg", #cl 18 
                 "NK CD56bright", #cl 19
                 "Proliferating TC", #cl 20
                 "gdTC", #cl 21
                 "CD4_memory", #cl 22
                 "DC", #cl 23
                 "Platelets", #cl 24
                 "pDC", #cl 25
                 "CD8_memory", #cl 26
                 "Naive BC", #cl 27
                 "CD4_CTL_P3", #cl 28
                 "CD4_CTL_P3", #cl 29
                 "hSC", #cl 30
                 "Plasma cells", #cl 31
                 "Naive BC", #cl 32
                 "Classical Mono" #cl 33
)

names(x = new.cluster) <- levels(x = covid_2)
covid_2 <- RenameIdents(object = covid_2, new.cluster)
covid_2$clusters <- Idents(object = covid_2)
table(covid_2$clusters)

## Supplementary Figure 1 h
DimPlot(covid_2, raster = F, group.by = "clusters",
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
        ), 
        order = c("Plasma cells" , 
                  "pDC" ,
                  "Memory BC" ,
                  "Naive BC",
                  "Classical Mono",
                  "NK CD56dim" ,
                  "hSC" ,
                  "NK CD56bright" ,
                  "CD4_naive" ,
                  "CD8_naive",
                  "CD8_memory" ,
                  "CD8_NKT_like",
                  "gdTC" ,
                  "MAIT" ,
                  "CD4_memory" ,
                  "Proliferating TC",
                  "Non-classical Mono",
                  "Treg",
                  "DC",
                  "Platelets" ,
                  "CD4_CTL_P3" ,
                  "Doublets"),
        pt.size = 0.01)


## Supplementary Figure 1 i
Idents(covid_2) <- "clusters"
DefaultAssay(covid_2) <- "RNA"

covid_2@active.ident <- factor(covid_2@active.ident, 
                               levels = (c("Treg", "CD4_naive", "CD4_memory", "CD4_CTL_P3", 
                                           "CD8_naive", "CD8_memory", "CD8_NKT_like", "gdTC", "MAIT", "Proliferating TC", 
                                           "NK CD56bright", "NK CD56dim",
                                           "Naive BC", "Memory BC", "Plasma cells", 
                                           "Classical Mono", "Non-classical Mono", "DC", "pDC",
                                           "Platelets", "hSC", "Doublets"
                               )))

pal<- colorRampPalette(brewer.pal(12, "Spectral"))(19)
VlnPlot(covid_2, features = c("CD3D", "FOXP3", "CD4", "CD8A", "FOXP3", "TRDV2", "SLC4A10", "MKI67", 
                              "NCAM1", "FCGR3A", "CD19", "CCR7", "IGHM", "JCHAIN", 
                              "CD33", "FCER1A", "CD14", "PPBP", "CD34"), 
        assay = "RNA",
        raster = F, pt.size = 0, stack = T, flip = T, adjust = T, cols = pal) + NoLegend()

# Annotation of general clusters
covid_2 <- FindNeighbors(covid_2, reduction = "pca", dims = 1:50)
covid_2 <- FindClusters(covid_2, resolution = 1.5)

DefaultAssay(covid_2) <- "RNA"
Idents(covid_2) <- "seurat_clusters"
current.cluster <- levels(covid_2)
current.cluster

new.cluster <- c("T cells", #cl 0
                 "T cells", #cl 1
                 "NK", #cl 2
                 "T cells", #cl 3
                 "T cells", #cl 4
                 "T cells", #cl 5
                 "B cells", #cl 6
                 "T cells", #cl 7
                 "Myeloids", #cl 8
                 "T cells", #cl 9
                 "T cells", #cl 10
                 "T cells", #cl 11
                 "T cells", #cl 12
                 "Myeloids", #cl 13
                 "T cells", #cl 14
                 "T cells", #cl 15
                 "Myeloids", #cl 16
                 "T cells", #cl 17
                 "B cells", #cl 18 
                 "T cells", #cl 19
                 "Myeloids", #cl 20
                 "B cells", #cl 21
                 "T cells", #cl 22
                 "NK", #cl 23
                 "T cells", #cl 24
                 "T cells", #cl 25
                 "T cells", #cl 26
                 "T cells", #cl 27
                 "Myeloids", #cl 28
                 "T cells", #cl 29
                 "Myeloids", #cl 30
                 "B cells", #cl 31
                 "T cells", #cl 32
                 "Myeloids", #cl 33
                 "B cells", #cl 34
                 "T cells", #cl 35
                 "NK", #cl 36
                 "B cells", #cl 37
                 "B cells", #cl 38 
                 "hSC", #cl 39
                 "Myeloids" #cl 40
)

names(x = new.cluster) <- levels(x = covid_2)
covid_2 <- RenameIdents(object = covid_2, new.cluster)
covid_2$general_cluster <- Idents(object = covid_2)
table(covid_2$general_cluster)


metadata <- covid_2@meta.data %>% 
  mutate(StageStudy = paste0(PID, "_", Stage)) %>% 
  mutate(Patient_Stage = ifelse(is.na(patientID) == T, paste0(patient_id, ":", StageStudy), paste0(patientID, ":", StageStudy) )) 

covid_2 <- AddMetaData(covid_2, metadata = metadata$Patient_Stage, col.name = "Patient_Stage") 
covid_2 <- AddMetaData(covid_2, metadata = metadata$StageStudy, col.name = "StageStudy" )

