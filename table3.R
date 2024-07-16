library(tidyverse)
library(rstatix)
load("~/Escritorio/datos_covid/newdata/integration_samples/cohort1/COVID_REANNOT.RData")
#Number of cells in total
covid_cvid@meta.data <- covid_cvid@meta.data %>% 
mutate(StageSeverity = paste0(Stage, "_", CovidSeverity),
       PID = ifelse(StudyName == "Javi", "CVID", "non-CVID"), 
       Patient_Stage = paste0(patientID, "_", PID, "_", StageSeverity))
  
cvid <- subset(covid_cvid, subset = PID == "CVID")

cvid_cells <- table(cvid$CLUSTERS, cvid$Stage) %>%  as.data.frame() %>% 
  tidyr::spread(value = Freq, key = Var2)

noncvid <- subset(covid_cvid, subset = PID == "non-CVID")

NONcvid_cells <- table(noncvid$CLUSTERS, noncvid$StageSeverity) %>%  as.data.frame() %>% 
  tidyr::spread(value = Freq, key = Var2)

number_cells <- left_join(cvid_cells, NONcvid_cells, "Var1") %>% 
  rename(cluster = Var1)

#Number of cells per patient
table(cvid$Patient_Stage)

cvid_cells_patient <- table(cvid$CLUSTERS, cvid$Patient_Stage) %>%  as.data.frame() %>% 
  separate(Var2, sep = "_", into = c("patient", "PID", "Stage"), remove = F, extra = "merge") %>% 
  filter(!Var1 %in% c("CD4_CTL_P3", "Doublets", "hSC"))
  
#CVID
CLUSTERS <- unique(cvid_cells_patient$Var1) %>% as.vector()
stage <- unique(cvid_cells_patient$Stage)
all_freq_patient_cvid <- data.frame()
all_freq_stages_cvid <- data.frame()
patients <- unique(cvid_cells_patient$patient)

for (s in stage) {
#Por cada stage
  filter_stage <- cvid_cells_patient %>% 
    filter(Stage == s) 
  
  for (p in patients) {
    #POR CADA PACIENTE
    filter_patient <-  filter_stage %>% 
      filter(patient == p)
  total <- sum(as.numeric(filter_patient$Freq))
  
    for (n in CLUSTERS) {
  #Mean por patients and CLUSTERS 
  x <- filter_patient %>% 
    filter(Var1 == n) %>% 
    mutate(percentage = Freq/total)
  
  all_freq_patient_cvid <- rbind(all_freq_patient_cvid, x)     
  }
  
  }

  #ira por stage
  f <- tapply(all_freq_patient_cvid$percentage, all_freq_patient_cvid$Var1, mean) %>% as.data.frame() %>% 
    mutate(stage = s) %>% 
  tibble::rownames_to_column("cluster")
  
  all_freq_stages_cvid <- rbind(all_freq_stages_cvid, f)
  
}


all_freq_stages_cvid <- all_freq_stages_cvid %>% na.omit() %>% 
  mutate(pat = paste0("patientwise_cvid", stage)) %>% 
  select(pat, ".", cluster) %>% 
  tidyr::spread(value = ".", key = pat)

#Number of cells per patient NONcvid_cells
table(noncvid$Patient_Stage)

noncvid_cells_patient <- table(noncvid$CLUSTERS, noncvid$Patient_Stage) %>%  as.data.frame() %>% 
  separate(Var2, sep = "_", into = c("patient", "PID", "Stage"), remove = F, extra = "merge") %>% 
  filter(!Var1 %in% c("CD4_CTL_P3", "Doublets", "hSC"))

#NONcvid_cells
CLUSTERS <- unique(noncvid_cells_patient$Var1)  %>% as.vector()
stage <- unique(noncvid_cells_patient$Stage)
all_freq_patient_noncvid <- data.frame()
all_freq_stages_non_cvid <- data.frame()
for (s in stage) {
  #Por cada stage
  filter_stage <- noncvid_cells_patient %>% 
    filter(Stage == s) 
  
patients <- unique(filter_stage$patient)

  for (p in patients) {
    #POR CADA PACIENTE
    filter_patient <-  filter_stage %>% 
      filter(patient == p)
    total <- sum(as.numeric(filter_patient$Freq))
    
    for (n in CLUSTERS) {
      #Mean por patients and CLUSTERS 
      x <- filter_patient %>% 
        filter(Var1 == n) %>% 
        mutate(percentage = Freq/total)
      
      all_freq_patient_noncvid <- rbind(all_freq_patient_noncvid, x)     
    }
    
  }
  
  #ira por stage
  f <- tapply(all_freq_patient_noncvid$percentage, all_freq_patient_noncvid$Var1, mean) %>% as.data.frame() %>% 
    mutate(stage = s) %>% 
  tibble::rownames_to_column("cluster")
  
  all_freq_stages_non_cvid <- rbind(all_freq_stages_non_cvid, f)
  
}

all_freq_stages_non_cvid <- all_freq_stages_non_cvid %>% na.omit() %>% 
  mutate(pat = paste0("patientwise_non_cvid", stage)) %>% 
  select(pat, ".", cluster) %>% 
  tidyr::spread(value = ".", key = pat)


#WILCOXON CALCULATE
all <- rbind(all_freq_patient_cvid, all_freq_patient_noncvid) %>% 
  mutate(PID_STAGE = paste0(PID, "_", Stage))
CLUSTERS <- unique(all$Var1)
#Calculate significance
sign_ttest <- data.frame()
wilcox <- data.frame()
for (n in CLUSTERS) {
  all_filter <- all %>% filter(Var1 == n) 
   # 
   # a <- rstatix::wilcox_test(formula = percentage ~ PID_STAGE, data = all_filter,
   #                          comparisons = list(c("CVID_baseline_Control", "non-CVID_baseline_Control"),
   #                                             c("CVID_progression_Mild", "non-CVID_progression_Mild"),
   #                                             c("CVID_convalescence_Mild", "non-CVID_convalescence_Mild"),
   #                                             c("CVID_progression_Mild","non-CVID_progression_Severe"),
   #                                             c("CVID_convalescence_Mild","non-CVID_convalescence_Severe") )) %>%
   #                                               mutate(cluster = n,
   #                                                      cluster_stage = paste0(n, "_", group2))
   # wilcox <- rbind(wilcox, a)
   
  
  #WELCH TEST
  stat.test <- all_filter %>% 
    t_test(percentage ~ PID_STAGE, alternative = "two.sided", var.equal = F, p.adjust.method = "none",
           comparisons = list(c("CVID_baseline_Control", "non-CVID_baseline_Control"),
                              c("CVID_progression_Mild", "non-CVID_progression_Mild"), 
                              c("CVID_convalescence_Mild", "non-CVID_convalescence_Mild"), 
                              c("CVID_progression_Mild","non-CVID_progression_Severe"), 
                              c("CVID_convalescence_Mild","non-CVID_convalescence_Severe") ) ) %>%
    mutate(cluster = n)  %>%
    add_significance()
  
  sign_ttest <- rbind(sign_ttest, stat.test)
  
}

test <- sign_ttest %>% 
  mutate(comparativa = paste0(group1, "_", group2) ) %>% 
  select(cluster, comparativa, p) %>% 
  tidyr::spread(value = p, key = comparativa)

bonito <- left_join(test, all_freq_stages_cvid, by = "cluster") %>% 
  left_join(all_freq_stages_non_cvid, by = "cluster") %>%
  left_join(number_cells, by = "cluster")  %>% 
  select(cluster, baseline, baseline_Control, patientwise_cvidbaseline_Control, patientwise_non_cvidbaseline_Control, "CVID_baseline_Control_non-CVID_baseline_Control",
         progression, progression_Mild, patientwise_cvidprogression_Mild, patientwise_non_cvidprogression_Mild, "CVID_progression_Mild_non-CVID_progression_Mild",
         progression_Severe, patientwise_non_cvidprogression_Severe, "CVID_progression_Mild_non-CVID_progression_Severe", 
         convalescence, convalescence_Mild, patientwise_cvidconvalescence_Mild, patientwise_non_cvidconvalescence_Mild, "CVID_convalescence_Mild_non-CVID_convalescence_Mild",
         convalescence_Severe, patientwise_non_cvidconvalescence_Severe, "CVID_convalescence_Mild_non-CVID_convalescence_Severe") 


writexl::write_xlsx(bonito, path = "supp_table3.xlsx")

#BARPLOT PERCENTAGES

plot <-tidyr::gather(bonito, key = "stage", value = "percentage", 
                c(patientwise_cvidbaseline_Control, patientwise_non_cvidbaseline_Control, patientwise_cvidprogression_Mild, patientwise_non_cvidprogression_Mild,
                  patientwise_non_cvidprogression_Severe, patientwise_cvidconvalescence_Mild, patientwise_non_cvidconvalescence_Mild, 
                  patientwise_non_cvidconvalescence_Severe)) %>% 
  select(cluster, stage, percentage)


plot$cluster <- factor(plot$cluster, 
                                  levels = (c("Treg", "CD4_naive", "CD4_memory", "CD4_CTL_P3", 
                                              "CD8_naive", "CD8_memory", "CD8_NKT_like", "gdTC", "MAIT", "Proliferating TC", 
                                              "NK CD56bright", "NK CD56dim",
                                              "Naive BC", "Memory BC", "Plasma cells", 
                                              "Classical Mono", "Non-classical Mono", "DC", "pDC",
                                              "Platelets", "hSC", "Doublets"
                                  )))


plot$stage <- factor(plot$stage, 
                       levels = rev(c("patientwise_cvidbaseline_Control", "patientwise_cvidprogression_Mild", "patientwise_cvidconvalescence_Mild", 
                                   "patientwise_non_cvidbaseline_Control", "patientwise_non_cvidprogression_Mild", "patientwise_non_cvidconvalescence_Mild", 
                                   "patientwise_non_cvidprogression_Severe", "patientwise_non_cvidconvalescence_Severe"
                       )) )

ggplot(plot, aes(x= stage, y = percentage*100)) +
  geom_bar(aes(fill = cluster), stat = "identity", width = 0.8, color = "gray38") +
  scale_fill_manual(values =  c("pDC" = "#03045e",
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
  )) + coord_flip() +
  # facet_grid(~PID, scales = "free_x") + 
  theme_pubr(legend = "right") + border() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ylab(label = "Percentage") + xlab(label = "") + labs(fill=" ") 

