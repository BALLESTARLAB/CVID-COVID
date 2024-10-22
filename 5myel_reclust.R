# MYELOID RECLUSTERING COHORT 1 ---------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)

load("~/myel_cohort1_samplenames_annot.RData")

#Annotation
myel_1 <- FindClusters(myel_1, resolution = 0.5, graph.name = "integrated_snn")

Idents(myel_1) <- "seurat_clusters"
current.cluster <- levels(myel_1)
current.cluster

new.cluster <- c("Classical_mono", #cl 0
                 "HLADRlow_S100Ahigh_mono", #cl 1
                 "HLADRlow_S100Ahigh_mono", #cl 2
                 "HLADRlow_S100Ahigh_mono", #cl 3
                 "Intermediate_mono", #cl 4 
                 "Non_classical_mono", #cl 5
                 "HLADRlow_CD163high_mono", #cl 6
                 "DC3", #cl 7 
                 "DC2", #cl 8
                 "Doublets myeloids/platelets", #cl 9
                 "HLADRhigh_CD83high_mono", #cl 10
                 "Non_classical_mono", #cl 11
                 "pDC", #cl 12
                 "DC1", #cl 13
                 "HLADRlow_S100Ahigh_mono" #cl 14
)

DefaultAssay(myel) <- "RNA"
names(x = new.cluster) <- levels(x = myel_1)
myel_1 <- RenameIdents(object = myel_1, new.cluster)
myel_1$general_cluster <- Idents(object = myel_1)
table(myel_1$general_cluster)
Idents(myel_1) <- "general_cluster"
# Figure 6a
new.col.myel <- c("#5CB32E", # "Classical_mono", 
                  "#969DCA", #  "HLADRlow_S100Ahigh_mono",
                  "#F6D832",#Intermediate_mono
                  "#DC9F2D", #  "Non_classical_mono",
                  "#9E8546", #"HLADRlow_CD163high_mono",
                  "#E3A39B", # DC3
                  "#66C2A5", # DC2
                  "#F3DAAC", # Doublets
                  "#B7D84C", # "HLADRhigh_CD83high_mono",
                  "#D58EC4", # pDC
                  "#BCBF2E" # DC1
                  
)

DimPlot(myel_1, label = F, group.by = "general_cluster", pt.size = 0.5, cols = new.col.myel) +
  labs(title = "Myeloids (n = 65,375)")

#Fig 6b
DefaultAssay(myel_1) <- "RNA"
p1 <- FeaturePlot(myel_1, features = c("FCGR3A"))  +
  scale_color_viridis_c(option = "C")
p2 <- FeaturePlot(myel_1, features = c("CD14"))  +
  scale_color_viridis_c(option = "C")
p1 + p2

#Fig 6c
markers.myel <- FindAllMarkers(myel_1, assay = "RNA", logfc.threshold = 0.25, min.pct = 0.3)

top.genes.all <- data.frame()
myel_1@active.ident <- factor(myel_1@active.ident, 
                              levels = rev(c("HLADRlow_CD163high_mono",
                                             "HLADRlow_S100Ahigh_mono", 
                                             "HLADRhigh_CD83high_mono",
                                             "Classical_mono", 
                                             "Intermediate_mono",
                                             "Non_classical_mono",
                                             "DC1", "DC2", "DC3", "pDC", 
                                             "Doublets myeloids/platelets")))

myel.clusters <- (levels(myel_1$general_cluster))

for (i in myel.clusters) {
  top.genes <- markers.myel %>% 
    mutate(sign = ifelse(p_val_adj < 0.05 & cluster == i & pct.1 > 0.4 & avg_log2FC > 0, paste("sign", i, sep = "_"), "others")) %>% 
    dplyr::filter(sign != "others")
  top.genes.all <- rbind(top.genes.all, top.genes)
}

gene.markers.myel <- (table(top.genes.all$gene)) %>% as.data.frame() %>% 
  dplyr::filter(Freq == 1) %>% 
  dplyr::filter(!grepl("^MT-", Var1))

top.genes.common <- top.genes.all %>% 
  dplyr::filter(gene %in% gene.markers.myel$Var1) %>% 
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) 


markers.myel.new <- c(top.genes.common$gene[1:20],"CD14",
                      top.genes.common$gene[21:25], "FCGR3A",
                      top.genes.common$gene[26:length(top.genes.common$gene)])

DefaultAssay(myel_1) <- "RNA"
a <- DotPlot(myel_1, features = markers.myel.new, scale.max = 100, scale.min = 0, dot.scale = 5)  +
  scale_color_gradientn(colors = c("skyblue3", "white", "#770000")) +
  theme_pubr(legend = "right") + border()+ rotate_x_text(45) 

# Subsampling  ------------
myel_1@meta.data <- myel_1@meta.data %>% 
  mutate(PID = ifelse(StudyName == "Javi", "CVID", "non-CVID") ) %>% 
  select(-StudyName)
table(myel_1$PID)

set.seed(123)

#CVID
sample1_rna <- subset(myel_1, subset = PID == "CVID")
metadata1 <- as.data.frame(sample1_rna@meta.data)
nuevo_m1 <- row.names(sample_n(metadata1, 19451))
cvid_myel_subsample <- subset(myel_1, cells = nuevo_m1)

noncvid_myel <- subset(myel_1, subset = PID == "non-CVID") 
#non-CVID mild
sample2_rna <- subset(noncvid_myel, subset = CovidSeverity == c("Severe"), invert = T)
metadata2 <- as.data.frame(sample2_rna@meta.data)
nuevo_m2 <- row.names(sample_n(metadata2, 19451))
noncvid_mild_myel_subsample <- subset(myel_1, cells = nuevo_m2)

#non-CVID severe
sample3_rna <- subset(noncvid_myel, subset = CovidSeverity == c("Mild"), invert = T)
metadata3 <- as.data.frame(sample3_rna@meta.data)
nuevo_m3 <- row.names(sample_n(metadata3, 19451))
noncvid_severe_myel_subsample <- subset(myel_1, cells = nuevo_m3)

# Supp Fig 6a
new.col.myel <- c("#5CB32E", # "Classical_mono", 
                  "#969DCA", #  "HLADRlow_S100Ahigh_mono",
                  "#F6D832",#Intermediate_mono
                  "#DC9F2D", #  "Non_classical_mono",
                  "#9E8546", #"HLADRlow_CD163high_mono",
                  "#E3A39B", # DC3
                  "#66C2A5", # DC2
                  "#F3DAAC", # Doublets
                  "#B7D84C", # "HLADRhigh_CD83high_mono",
                  "#D58EC4", # pDC
                  "#BCBF2E" # DC1
                  
)
p1 <- DimPlot(cvid_myel_subsample, group.by = "general_cluster", pt.size = 0.25, cols = new.col.myel, raster = F) + labs(title = "CVID") + NoLegend()
p2 <- DimPlot(noncvid_mild_myel_subsample, group.by = "general_cluster", pt.size = 0.25, cols = new.col.myel, raster = F) + labs(title = "non-CVID mild") + NoLegend()
p3 <- DimPlot(noncvid_severe_myel_subsample, group.by = "general_cluster", pt.size = 0.25, cols = new.col.myel, raster = F) + labs(title = "non-CVID severe")
p1+p2+p3

# Supp Fig6b
a<- DotPlot(myel_1, features = c("IFIT3", "MX1", "IFI6", "ISG15", "TNFSF10", #CD163 high mono
                             "S100A8","S100A9", "S100A12", "MAFB", "ATP2B1-AS1", #HLA-DRlow S100Ahigh mono
                             "IL1B", "SGK1", "CLEC7A", "CD83", #HLA DRhigh_CD83high_mono
                             "AP1S2", "CD14", "CSTA", "TYROBP", "LYZ", "CST3", "FCN1", "MS4A6A", #Classical_mono
                             "FCGR3A", "MS4A7", "CDKN1C", "COTL1", "RHOC", #Non_clasical mono
                             "CLEC9A","THBD", #DC1
                             "HLA-DRA", "HLA-DRB1", "CD1C", "CLEC10A", #dc2
                             "CCL5", "IL32", "GNLY", #DC3
                             "JCHAIN", "GZMB", "ITM2C", "LILRA4", #pDC
                             "PPBP", "TUBB1", "PF4", "CAVIN2"
), scale.max = 100, scale.min = 0, dot.scale = 4, scale.by = "size")  +
  scale_color_gradientn(colors = c("skyblue3", "white", "#770000")) +
  theme_pubr(legend = "right") + border()+ rotate_x_text(45)

#Heatmap scores -----------------------
#Supp Fig 6e  

# Clasical mono -----------
classical <- subset(myel_1, subset = general_cluster == "Classical_mono", invert = F)

#GO:0009615_response TYPE_I_INTERFERON CLUSTER 10
IFN <- read.table("markers/myel_cluster//GOBP_RESPONSE_TO_TYPE_I_INTERFERON.v2023.1.Hs.grp", header = T)
classical <- AddModuleScore(classical, features = list(unique(IFN$GOBP_RESPONSE_TO_TYPE_I_INTERFERON)),
                            name = c("cl10_IFN_"), search = T)

#GO:0035455_response to interferon-alpha CLUSTER 11
alpha <- read.table("markers/myel_cluster/GOBP_RESPONSE_TO_INTERFERON_ALPHA.v2023.2.Hs.grp", header = T)
classical <- AddModuleScore(classical, features = list(unique(alpha$GOBP_RESPONSE_TO_INTERFERON_ALPHA)),
                            name = c("cl11_alpha_"), search = T)

#GO:0002675_positive regulation of acute inflammatory response CLUSTER 6
inflammatory <- read.table("markers/myel_cluster/GOBP_POSITIVE_REGULATION_OF_ACUTE_INFLAMMATORY_RESPONSE_TO_ANTIGENIC_STIMULUS.v2023.2.Hs.grp", header = T)
classical <- AddModuleScore(classical, features = list(unique(inflammatory$GOBP_POSITIVE_REGULATION_OF_ACUTE_INFLAMMATORY_RESPONSE_TO_ANTIGENIC_STIMULUS)),
                            name = c("cl6_inflammatory_"), search = T)

#GO:0140632_inflammasome complex assembly CLUSTER 13
assembly <- read.table("markers/myel_cluster/GOBP_INFLAMMASOME_COMPLEX_ASSEMBLY.grp", header = T)
classical <- AddModuleScore(classical, features = list(unique(assembly$GENE)),
                            name = c("cl13_assembly_"), search = T)

#GO:0032611_interleukin-1 beta production CLUSTER 13
interleukin <- read.table("markers/myel_cluster/GOBP_INTERLEUKIN_1_PRODUCTION.v2023.2.Hs.grp", header = T)
classical <- AddModuleScore(classical, features = list(unique(interleukin$GOBP_INTERLEUKIN_1_PRODUCTION)),
                            name = c("cl13_interleukin_"), search = T)

score <- classical@meta.data %>% 
  select(cl10_IFN_1, cl11_alpha_1, cl6_inflammatory_1, cl13_assembly_1, cl13_interleukin_1) 
plot <- classical@meta.data %>%
  mutate(StudyName = ifelse(StudyName == "China", "nonCVID", "CVID")) %>% 
  mutate(sample = paste0(Stage, "_", CovidSeverity, "-", StudyName, ":", general_cluster)) %>%
  mutate(patient = paste0(sample, ":", StageCVID)) %>%
  select(general_cluster, sample, StageCVID, patient, cl10_IFN_1, cl11_alpha_1, cl6_inflammatory_1, cl13_assembly_1, cl13_interleukin_1)

#Plot por los clusters de los scores
#Mean per stage  
cl10_IFN_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl10_IFN_1), sd=sd(cl10_IFN_1)) %>% 
  mutate(cluster = "cl10_IFN_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl11_alpha_1 <- plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl10_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl11_alpha_1), sd=sd(cl11_alpha_1)) %>% 
  mutate(cluster = "cl11_alpha_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl6_inflammatory_1 <-plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl13_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl6_inflammatory_1), sd=sd(cl6_inflammatory_1)) %>% 
  mutate(cluster = "cl6_inflammatory_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl13_assembly_1 <-plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl13_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl13_assembly_1), sd=sd(cl13_assembly_1)) %>% 
  mutate(cluster = "cl13_assembly_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl13_interleukin_1 <-plot  %>%  # plyr::ddply(plot, "patient", summarise, media=mean(cl13_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl13_interleukin_1), sd=sd(cl13_interleukin_1)) %>% 
  mutate(cluster = "cl13_interleukin_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

# 
# 
all_pop_plot_median <- rbind (cl10_IFN_1, cl11_alpha_1 , cl6_inflammatory_1, cl13_assembly_1, cl13_interleukin_1 )

all_pop_plot_median$stage <- factor(all_pop_plot_median$stage, levels = c("baseline", "progression", "convalescence"))
writexl::write_xlsx(all_pop_plot_median, path = "new_images_def/NEW/classical_GO_nocorrect.xlsx")

#PLOTS cl10_IFN_1, cl11_alpha_1 , cl6_inflammatory_1, cl13_assembly_1, cl13_interleukin_1
cl10_IFN_1$PID_stage <- factor(cl10_IFN_1$PID_stage, levels = c("CVID:baseline_Control", "CVID:progression_Mild", "CVID:convalescence_Mild",
                                                                                "nonCVID:baseline_Control", "nonCVID:progression_Mild", "nonCVID:convalescence_Mild",
                                                                                "nonCVID:progression_Severe", "nonCVID:convalescence_Severe"))

ggplot(cl10_IFN_1, aes(x=PID_stage, y=cluster, fill = mean) ) +
  geom_tile(aes(fill=mean), color = "black")  +
  facet_grid(~cluster, scales = "free") +
  scale_fill_gradientn(colours = c("#04508b", "#087ca7", "#E3B264","#FAE7C8")) +
  ggpubr::theme_pubr(legend = "right", border = T, x.text.angle = 45) 


# Intermediate mono -----------
intermediate <- subset(myel_1, subset = general_cluster == "Intermediate_mono", invert = F)

#GO:0035455_response to interferon-alpha CLUSTER 9
alpha <- read.table("markers/myel_cluster/GOBP_RESPONSE_TO_INTERFERON_ALPHA.v2023.2.Hs.grp", header = T)
intermediate <- AddModuleScore(intermediate, features = list(unique(alpha$GOBP_RESPONSE_TO_INTERFERON_ALPHA)),
                               name = c("cl9_alpha_"), search = T)

#GO:0061702_inflammasome complex CLUSTER 13 
inflammasome <- read.table("markers/myel_cluster/GOCC_CANONICAL_INFLAMMASOME_COMPLEX.v2023.2.Hs.grp", header = T)
intermediate <- AddModuleScore(intermediate, features = list(unique(inflammasome$GOCC_CANONICAL_INFLAMMASOME_COMPLEX)),
                               name = c("cl13_inflammasome_"), search = T)

#GO:0032611_interleukin-1 beta production CLUSTER 16 
interleukin <- read.table("markers/myel_cluster/GOBP_INTERLEUKIN_1_PRODUCTION.v2023.2.Hs.grp", header = T)
intermediate <- AddModuleScore(intermediate, features = list(unique(interleukin$GOBP_INTERLEUKIN_1_PRODUCTION)),
                               name = c("cl16_interleukin_"), search = T)

#GO:00096151_response to virus
virus <- read.table("markers/myel_cluster/GOBP_RESPONSE_TO_VIRUS.v2023.2.Hs.grp", header = T)
intermediate <- AddModuleScore(intermediate, features = list(unique(virus$GOBP_RESPONSE_TO_VIRUS)),
                               name = c("virus_"), search = T)

#GO:00516072_defense response to virus
defense <- read.table("markers/myel_cluster/GO:0051607.grp", header = T)
intermediate <- AddModuleScore(intermediate, features = list(unique(defense$GENE)),
                               name = c("defense_"), search = T)

plot <- intermediate@meta.data %>%
  mutate(StudyName = ifelse(StudyName == "China", "nonCVID", "CVID")) %>% 
  mutate(sample = paste0(Stage, "_", CovidSeverity, "-", StudyName, ":", general_cluster)) %>%
  mutate(patient = paste0(sample, ":", StageCVID)) %>%
  select(general_cluster, sample, StageCVID, patient, cl9_alpha_1, cl13_inflammasome_1, cl16_interleukin_1, virus_1, defense_1)

## Summary plot intento
#Mean per stage  
cl9_alpha_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl9_alpha_1), sd=sd(cl9_alpha_1)) %>% 
  mutate(cluster = "cl9_alpha_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl13_inflammasome_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl13_inflammasome_1), sd=sd(cl13_inflammasome_1)) %>% 
  mutate(cluster = "cl13_inflammasome_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl16_interleukin_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl16_interleukin_1), sd=sd(cl16_interleukin_1)) %>% 
  mutate(cluster = "cl16_interleukin_1")  %>% 
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

all_pop_plot_median <- rbind (cl9_alpha_1, cl13_inflammasome_1 , cl16_interleukin_1, virus_1, defense_1)

all_pop_plot_median$stage <- factor(all_pop_plot_median$stage, levels = c("baseline", "progression", "convalescence"))
writexl::write_xlsx(all_pop_plot_median, path = "new_images_def/NEW/intermediate_GO_nocorrect.xlsx")

#PLOTS (cl9_alpha_1, cl13_inflammasome_1 , cl16_interleukin_1, virus_1, defense_1)
cl9_alpha_1$PID_stage <- factor(cl9_alpha_1$PID_stage, levels = c("CVID:baseline_Control", "CVID:progression_Mild", "CVID:convalescence_Mild",
                                                              "nonCVID:baseline_Control", "nonCVID:progression_Mild", "nonCVID:convalescence_Mild",
                                                              "nonCVID:progression_Severe", "nonCVID:convalescence_Severe"))

ggplot(cl9_alpha_1, aes(x=PID_stage, y=cluster, fill = mean) ) +
  geom_tile(aes(fill=mean), color = "black")  +
  facet_grid(~cluster, scales = "free") +
  scale_fill_gradientn(colours = c("#04508b", "#087ca7", "#E3B264","#FAE7C8")) +
  ggpubr::theme_pubr(legend = "right", border = T, x.text.angle = 45) 

# Non-classical mono -----------
nonclassical <- subset(myel_1, subset = general_cluster == "Non_classical_mono", invert = F)

#GO:0032606_type I interferon production CLUSTER 9
IFN <- read.table("markers/myel_cluster/GOBP_TYPE_I_INTERFERON_PRODUCTION.v2023.2.Hs.grp", header = T)
nonclassical <- AddModuleScore(nonclassical, features = list(unique(IFN$GOBP_TYPE_I_INTERFERON_PRODUCTION)),
                               name = c("cl9_IFN_"), search = T)

#GO:0061702_inflammasome complex CLUSTER 11
inflammasome <- read.table("markers/myel_cluster/GOCC_CANONICAL_INFLAMMASOME_COMPLEX.v2023.2.Hs.grp", header = T)
nonclassical <- AddModuleScore(nonclassical, features = list(unique(inflammasome$GOCC_CANONICAL_INFLAMMASOME_COMPLEX)),
                               name = c("cl11_inflammasome_"), search = T)

#GO:0032611_interleukin-1 beta production CLUSTER 11
interleukin <- read.table("markers/myel_cluster/GOBP_INTERLEUKIN_1_PRODUCTION.v2023.2.Hs.grp", header = T)
nonclassical <- AddModuleScore(nonclassical, features = list(unique(interleukin$GOBP_INTERLEUKIN_1_PRODUCTION)),
                               name = c("cl11_interleukin_"), search = T)

#GO:0009615_response to virus 
virus <- read.table("markers/myel_cluster/GOBP_RESPONSE_TO_VIRUS.v2023.2.Hs.grp", header = T)
nonclassical <- AddModuleScore(nonclassical, features = list(unique(virus$GOBP_RESPONSE_TO_VIRUS)),
                               name = c("VIRUS_"), search = T)

plot <- nonclassical@meta.data %>%
  mutate(StudyName = ifelse(StudyName == "China", "nonCVID", "CVID")) %>% 
  mutate(sample = paste0(Stage, "_", CovidSeverity, "-", StudyName, ":", general_cluster)) %>%
  mutate(patient = paste0(sample, ":", StageCVID)) %>%
  select(general_cluster, sample, StageCVID, patient, cl9_IFN_1, cl11_inflammasome_1, cl11_interleukin_1, VIRUS_1)

## Summary plot intento
#Mean per stage  
cl9_IFN_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl9_IFN_1), sd=sd(cl9_IFN_1)) %>% 
  mutate(cluster = "cl9_IFN_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl11_inflammasome_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl11_inflammasome_1), sd=sd(cl11_inflammasome_1)) %>% 
  mutate(cluster = "cl11_inflammasome_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))

cl11_interleukin_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
  separate(patient, into = c("PID_Sev", "cluster", "patient"), sep = ":", remove = F, extra = "merge") %>% 
  separate(PID_Sev, into = c( "severity", "PID"), sep = "-", remove = F) %>% 
  separate(severity, into = c("stage", "severity"), sep = "_", remove = F, extra = "merge")  %>% 
  mutate(PID_stage = paste0(PID, ":", stage, "_", severity)) %>% 
  plyr::ddply("PID_stage", summarise, mean=mean(cl11_interleukin_1), sd=sd(cl11_interleukin_1)) %>% 
  mutate(cluster = "cl11_interleukin_1")  %>% 
  separate(PID_stage, into = c("PID", "Stage"), sep = ":",remove = F) %>% 
  separate(Stage, into = c("stage", "severity"), sep = "_") %>% 
  mutate(group = paste0(PID_stage, "_", cluster))


VIRUS_1 <- plot  %>%  #plyr::ddply(plot, "patient", summarise, media =mean(cl6_1)) %>% 
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
# 
all_pop_plot_median <- rbind (cl9_IFN_1, cl11_inflammasome_1 , cl11_interleukin_1, VIRUS_1)

all_pop_plot_median$stage <- factor(all_pop_plot_median$stage, levels = c("baseline", "progression", "convalescence"))
writexl::write_xlsx(all_pop_plot_median, path = "new_images_def/NEW/nonclassical_GO_nocorrect.xlsx")


#PLOTS cl9_IFN_1, cl11_inflammasome_1 , cl11_interleukin_1, VIRUS_1)
cl9_IFN_1$PID_stage <- factor(cl9_IFN_1$PID_stage, levels = c("CVID:baseline_Control", "CVID:progression_Mild", "CVID:convalescence_Mild",
                                                          "nonCVID:baseline_Control", "nonCVID:progression_Mild", "nonCVID:convalescence_Mild",
                                                          "nonCVID:progression_Severe", "nonCVID:convalescence_Severe"))

ggplot(cl9_IFN_1, aes(x=PID_stage, y=cluster, fill = mean) ) +
  geom_tile(aes(fill=mean), color = "black")  +
  facet_grid(~cluster, scales = "free") +
  scale_fill_gradientn(colours = c("#04508b", "#087ca7", "#E3B264","#FAE7C8")) +
  ggpubr::theme_pubr(legend = "right", border = T, x.text.angle = 45) 

# Pseudotime analysis ---------
Idents(myel_1) <- "general_cluster"
myel_1 <- subset(myel_1, idents = c("Doublets myeloids/platelets", "DC3", "pDC", "DC1", "DC2"), invert = T)
myel_1@meta.data$general_cluster <- droplevels(myel_1@meta.data$general_cluster)

cds <- as.cell_data_set(myel_1); gc() 
fData(cds)$gene_short_name <- rownames(fData(cds))

reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)
cds@clusters$UMAP$partitions <- reacreate.partition

list_cluster <- myel_1@active.ident
cds@clusters$UMAP$clusters <- list_cluster
cds@int_colData@listData$reducedDims$UMAP <- myel_1@reductions$umap@cell.embeddings
cds <- learn_graph(cds, use_partition = FALSE)

cds <- order_cells(cds, reduction_method = 'UMAP')

save(cds, file = "newdata/integration_samples/cohort1/myel/myel_coh1_cds.RData")

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

cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

data.pseudo <- data.pseudo %>% 
  mutate(Trajectory = ifelse(general_cluster == "Non_classical_mono" | general_cluster ==  "Intermediate_mono" | general_cluster == "Classical_mono", "Traj1", "Traj2"))

data.pseudo$general_cluster <- factor(data.pseudo$general_cluster, levels = c("HLADRhigh_CD83high_mono", 
                                                                              "HLADRlow_S100Ahigh_mono",
                                                                              "HLADRlow_CD163high_mono",
                                                                              "Classical_mono",
                                                                              "Intermediate_mono",
                                                                              "Non_classical_mono") )


# Supp Figure 6g
ggplot(data.pseudo, 
       aes(monocle3_pseudotime, 
           general_cluster, fill = general_cluster), color = general_cluster) +
  geom_boxplot(size = 0.1, outlier.size = 0.1) +
  theme_pubr(border = T, legend = "right") +
  scale_fill_manual(values = c("Classical_mono" = "#5CB32E", 
                               "HLADRlow_S100Ahigh_mono" = "#969DCA",
                               "Intermediate_mono" = "#F6D832",
                               "Non_classical_mono" = "#DC9F2D", 
                               "HLADRlow_CD163high_mono" = "#9E8546",
                               "DC3" = "#E3A39B", 
                               "DC2" = "#66C2A5",
                               "Doublets" = "#F3DAAC", 
                               "HLADRhigh_CD83high_mono" = "#B7D84C",
                               "pDC" = "#D58EC4",
                               "DC1" = "#BCBF2E") ) +
  NoLegend() +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) +
  xlab("") + ylab("")

#DEGs pseudotime in seurat
deg <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 10)

deg_wr <- deg %>% 
  tibble::rownames_to_column("genes")
writexl::write_xlsx(deg_wr, path = "psuedotime/myel_cvid_pseudotime.xlsx")
deg <- readxl::read_xlsx("psuedotime/myel_cvid_pseudotime.xlsx")

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

myel_1$pseudotime <- pseudotime(cds)
Idents(myel_1) <- myel_1$general_cluster

# PSEUDOTIME REGULONS
data("dorothea_hs", package = "dorothea")
regulons = dorothea_hs %>%
  filter(confidence %in% c("A", "B")) %>% 
  mutate(direction = ifelse(mor == "-1", "DOWN", "UP"))

regulon.n <- regulons %>% 
  filter(tf %in% c("SPI1")) %>% 
  dplyr::filter(direction == "UP")

myel_1 <- AddModuleScore(myel_1, features = list(regulon.n$target), name = "n")

plot.pseudo <- myel_1@meta.data %>%
  mutate(PID = ifelse(StudyName == "Javi", "CVID", "non-CVID")) %>% 
  mutate(Stage_CVID = paste0(PID, "_", Stage, "_", CovidSeverity)) %>% 
  dplyr::select(Stage_CVID, pseudotime, n1, general_cluster, PID, Stage, CovidSeverity) %>%
  rownames_to_column("cell") 

plot.pseudo$Stage <- factor(plot.pseudo$Stage, levels = c("baseline", "progression", "convalescence"))

colors <- c("Classical_mono" = "#5CB32E", # "Classical_mono", 
            "HLADRlow_S100Ahigh_mono" = "#969DCA", #  "HLADRlow_S100Ahigh_mono",
            "Intermediate_mono" = "#F6D832", #Intermediate_mono
            "Non_classical_mono" = "#DC9F2D", #"Non_classical_mono",
            "HLADRlow_CD163high_mono" = "#9E8546", #"HLADRlow_CD163high_mono",
            "DC3" = "#E3A39B", # DC3
            "DC2" = "#66C2A5", # DC2
            "Doublets" = "#F3DAAC", # Doublets
            "HLADRhigh_CD83high_mono" = "#B7D84C", # "HLADRhigh_CD83high_mono",
            "pDC" = "#D58EC4", # pDC
            "DC1" = "#BCBF2E", # DC1
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
  mutate(PID_severity = ifelse(Stage_CVID == "non-CVID_baseline_Control", "non-CVID_Severe", PID_severity)) %>% 
  mutate(Trajectory = ifelse(general_cluster == "Non_classical_mono" | general_cluster ==  "Intermediate_mono" | general_cluster == "Classical_mono", "Traj1", "Traj2"))

#TRAYECTORIA 1
plot.pseudo.nuevo.traj1 <- plot.pseudo.nuevo %>% 
  dplyr::filter(Trajectory == "Traj1") %>% 
  dplyr::filter(pseudotime > 0 & pseudotime < 15) 

#TRAYECTORIA 2
plot.pseudo.classical <- plot.pseudo.nuevo %>% 
  dplyr::filter(general_cluster == "Classical_mono") %>% 
  mutate(Trajectory = "Traj2") 

plot.pseudo.nuevo.traj2 <- plot.pseudo.nuevo %>% 
  dplyr::filter(Trajectory == "Traj2") %>% 
  dplyr::filter(pseudotime > 0 & pseudotime < 20) 

plot.pseudo.nuevo.traj2 <- rbind(plot.pseudo.nuevo.traj2, plot.pseudo.classical)

#Fig 6g
ggplot(plot.pseudo.nuevo.traj1 , aes(x = pseudotime, n1, color = general_cluster, fill = general_cluster)) +
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
  ylim(0, 0.4) + 
  xlab("") + ylab("")

#Supp Fig 7c
ggplot(plot.pseudo.nuevo.traj2 , aes(x = pseudotime, n1, color = general_cluster, fill = general_cluster)) +
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
  ylim(0, 0.4) + 
  xlab("") + ylab("")

# Cohort 2 ------------
load("~/myel_cohort2_samplenames.RData")

myel_2 <- FindClusters(myel, resolution = 0.6, graph.name = "integrated_snn") 

Idents(myel_2) <- "seurat_clusters"
current.cluster <- levels(myel_2)
current.cluster

new.cluster <- c("HLADRlow_S100Ahigh_mono", #cl 0
                 "HLADRlow_S100Ahigh_mono", #cl 1
                 "Classical_mono", #cl 2
                 "Intermediate_mono", #cl 3
                 "HLADRhigh_CD83high_mono", #cl 4 
                 "DC3", #cl 5
                 "Non_classical_mono", #cl 6
                 "HLADRlow_CD163high_mono", #cl 7 
                 "DC2", #cl 8
                 "pDC", #cl 9
                 "Doublets myeloids/platelets", #cl 10
                 "Non_classical_mono", #cl 11
                 "HLADRlow_CD163high_mono", #cl 12
                 "Classical_mono", #cl 13
                 "HLADRlow_S100Ahigh_mono" #cl 14
)

DefaultAssay(myel_2) <- "RNA"
names(x = new.cluster) <- levels(x = myel_2)
myel_2 <- RenameIdents(object = myel_2, new.cluster)
myel_2$general_cluster <- Idents(object = myel_2)
table(myel_2$general_cluster)
Idents(myel_2) <- "general_cluster"
myel_2@active.ident <- factor(myel_2@active.ident, 
                              levels = rev(c("HLADRlow_CD163high_mono",
                                             "HLADRlow_S100Ahigh_mono", 
                                             "HLADRhigh_CD83high_mono",
                                             "Classical_mono", 
                                             "Intermediate_mono",
                                             "Non_classical_mono",
                                             "DC2", "DC3", "pDC", 
                                             "Doublets myeloids/platelets")))

g <- readxl::read_xlsx(path = "source_data/fig6c.xlsx")
#Supp fig 6d
DotPlot(myel_2, features = c(top.genes.common$gene[1:50], "PPBP"), scale.max = 100, scale.min = 0, dot.scale = 4, scale.by = "size")  +
  scale_color_gradientn(colors = c("skyblue3", "white", "#770000")) +
  theme_pubr(legend = "right") + border()+ rotate_x_text(45)
