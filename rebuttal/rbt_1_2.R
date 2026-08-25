#!/usr/bin/env Rscript

# Author: Wei (Adelyn) Tsai; tsai.wei@mayo.edu
# If you are using this script, please cite our study
# This script addresses comment 2 about hippocampus and entorhinal cortex specificity from reviewer 1.
# Please note that based on the 2nd rebuttal, regional specificity results have been updated. Please see rebuttal_2nd/rbt_1_3.Rmd for more details.

suppressPackageStartupMessages(library(tidyverse)) #v2.0.0
suppressPackageStartupMessages(library(WGCNA)) #v1.72
suppressPackageStartupMessages(library(ggrepel)) #v0.9.8
suppressPackageStartupMessages(library(openxlsx)) #v4.2.8

source("./Codes/fns.R") #source a list of frequently used functions

# list of input/output directories
indir <- "./indir"
residir <- "./outdir/resids"
blwgcna_dir <- "./outdir/bl_wgcna"
brwgcna_dir <- "./outdir/br_wgcna"

outdir <- "./Rebuttal/outdir"

#######################################
# we test whether blood M3 is most preserved in PHG compared to other brain regions using modulePreservation function in WGCNA
comRes <- readRDS(paste0(residir, "/comRes.rds")) # from 1_CreateResids.R 
comRes_ls <- list("MCSA"=comRes[1:105,],
                  "ADNI"=comRes[106:196,])
comTrait_br <- readRDS(paste0(indir, "/comTrait_br.rds")) # a file with organized AMP-AD metadata
consNet_bl <- readRDS(paste0(blwgcna_dir, "/consNet_object_2023-08-15.rds")) #this is the blood consensus network created from 2_bl_wgcna.R
comRes_br <- readRDS(paste0(residir,"/comRes_br_lmer.rds")) # this is the brain residual files created from 3_br_wgcna.R

mp_bb_byReg <-vector("list",7)
names(mp_bb_byReg) <- unique(comTrait_br$Region)
for (r in names(mp_bb_byReg)){
  print(paste0("start:",r))
  comRes_br_sub <- comRes_br[which(rownames(comRes_br) %in% subset(comTrait_br, Region==r)$SampleID),]
  
  mp = WGCNA::modulePreservation(
    multiData = list(
      Blood = list(data = comRes),
      Brain = list(data = comRes_br_sub)
    ),
    multiColor = list(Blood = consNet_bl$consNet$colors), # Blood is the reference network and its colors are being tested in each brain dataset
    networkType = "signed", corFnc = "bicor", randomSeed = 1, nPermutations=100, parallelCalculation = T, verbose=3
  )
  mp_bb_byReg[[r]] <- mp
  
  rm(comRes_br_sub, mp)
  print(paste0("finish:",r))
}
saveRDS(mp_bb_byReg, paste0(outdir, "/1_2/mp_bb_byReg.rds"))


# plot the results 
# palette for the blood WGCNA modules
palette <- WGCNA::labels2colors(c(1:36)) 
names(palette) <- paste0("M", seq_along(1:36))

# prepare dataframe for plotting
mp_bb_byReg_df <- list(
  `Mayo STG` = mp_bb_byReg[[1]]$preservation$Z$ref.Blood$inColumnsAlsoPresentIn.Brain, 
  `Mayo CER` = mp_bb_byReg$CER$preservation$Z$ref.Blood$inColumnsAlsoPresentIn.Brain,
  `MSSM FP` = mp_bb_byReg$FP$preservation$Z$ref.Blood$inColumnsAlsoPresentIn.Brain,
  `MSSM STG` = mp_bb_byReg[[4]]$preservation$Z$ref.Blood$inColumnsAlsoPresentIn.Brain,
  `MSSM PHG` = mp_bb_byReg$PHG$preservation$Z$ref.Blood$inColumnsAlsoPresentIn.Brain,
  `MSSM IFG` = mp_bb_byReg$IFG$preservation$Z$ref.Blood$inColumnsAlsoPresentIn.Brain,
   `ROSMAP DLPFC` = mp_bb_byReg$DLPFC$preservation$Z$ref.Blood$inColumnsAlsoPresentIn.Brain
) %>%
  purrr::map(tibble::rownames_to_column, var = "Module") %>%
  bind_rows(.id = "Study") %>%
  filter(!Module %in% c("0", "0.1")) %>%
  mutate(Module = as.numeric(Module)) %>%
  arrange(desc(Study), Module)  %>%
  mutate(
    Study = factor(Study, levels=c("Mayo STG", "Mayo CER", "MSSM FP", "MSSM STG", "MSSM PHG", "MSSM IFG", "ROSMAP DLPFC")),
    moduleSize = rep(as.numeric(table(consNet_bl$consNet$colors)[2:37]),7),
    Module = str_c("M", as.character(Module)),
    Module = factor(Module, levels = paste0("M",seq_along(1:36))),
    mod_fill_color = palette[Module],
    mod_txt_color = get_max_contrast_color(mod_fill_color))


pdf(paste0(outdir, "/1_2/mp_bb_byReg.pdf"), width = 15, height = 7)
mp_bb_byReg_df %>% 
  ggplot(aes(x = moduleSize, y = Zsummary.pres)) +
  geom_point(aes(fill = Module), shape = 21, color = "black", stroke = 1, size = 4) +
  coord_cartesian(clip = "off") +
  geom_label_repel(aes(label = Module, fill = Module, color = mod_txt_color), data = subset(mp_bb_byReg_df, Module == "M3"), max.overlaps = Inf, size=6/.pt) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 5, linetype = "dashed", color = "darkgreen") +
  geom_hline(yintercept = 2, linetype = "dashed", color = "blue") +
  scale_color_identity() +
  scale_fill_manual(values = palette) +
  facet_grid(cols =vars(Study)) +
  labs(
    x= "Module Size",
    y ="Preservation Z Summary"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 6.5),
    strip.background = element_rect(color = "black", fill = "white", linewidth = 1),
    panel.background = element_rect(color = "black"),
    axis.line = element_blank(),
    text = element_text(size=7)
  )
dev.off() # this is in first rebuttal Rebuttal Fig.3
#######################################

#######################################
# We plotted the module-trait association in brain networks using forest plot
load(paste0(brwgcna_dir, "/METr.RData"), verbose = T) # created from 3_br_wgcna.R
pdf(paste0(outdir, "/1_2/BrMETr_forestplot.pdf"), width = 9, height = 5, onefile = T)
METr_allNoNA %>% 
  mutate(
    Study = case_when(
      Region == "CER" ~ "Mayo\nCER",
      Region %in% c("FP", "PHG", "IFG") ~ paste0("MSSM\n", Region),
      Region == "DLPFC" ~ "ROSMAP\nDLPFC",
      TRUE ~ Region
    ),
    Study =factor(Study, levels = c("ROSMAP\nDLPFC", "MSSM\nIFG", "MSSM\nPHG", "MSSM\nSTG", "MSSM\nFP", "Mayo\nCER", "Mayo\nSTG")),
    #Study =factor(Region, levels = c("Mayo\nSTG", "CER", "FP", "MSSM\nSTG", "PHG", "IFG", "DLPFC")),
    Shape = case_when(
      p.value > 0.05 ~ "NS",
      p.value < 0.05 & q.value > 0.05 ~ "Nominal",
      q.value < 0.05 ~ "Significant"
    ),
    term = if_else(term=="DX", "CTRL DX", term),
    term = factor(term, levels = c("CTRL DX", "CERAD", "Thal", "Braak"))
  ) %>% 
  ggplot(aes(x=estimate, y=Study,col=Module)) + 
  geom_errorbar(aes(xmax = conf.high, xmin = conf.low, color = Module), size = 0.5, width = 0.5, position=position_dodge(width = 0.5)) + 
  #specify position here too
  geom_point(aes(color = Module, shape = Shape), size = 2, position=position_dodge(width = 0.5)) +
  scale_shape_manual(values = c(NS = 4, Nominal = 17, Significant = 16)) +
  #geom_vline(xintercept = c(0.6, 0.65, 0.7), linetype = "dashed", color="grey40") +
  facet_grid(~term, scales = "free", space = "free") +
  theme_bw() +
  theme(text = element_text(size = 7))
dev.off() # this is in first rebuttal Rebuttal Fig.2
#######################################

#######################################
# to test the blood M3 association with hippocampal subfield

# prepare files
adni_hp_subfield <- read.csv(paste0(indir, "/ADNI_Hippo_CA.csv")) # this file contains information of baseline hippocampal volume subfield (BL_CA1,BL_CA2_3,BL_CA4_DG). Rows are subjects columns are volumes.
adni_mri <- read.csv(paste0(indir, "/ADNI_s91_MRI.csv")) # this file contains magnetic field, ICV and whole hippocampal volume information across multiple time lines
bl_ME <- read.csv(paste0(blwgcna_dir, "/MEs1_SP12.csv")) #generated from 2_bl_wgcna.R
comTrait <- readRDS(paste0(indir,"/comTrait.rds")) # a file having cross-sectional endophenotypes for MCSA (from row 1-105) and ADNI (from row 106-196); row names are donorID

adni_hp_subfield %>% 
  inner_join(bl_ME %>% dplyr::select(SubjectID, ME3), by=c("SubjID" = "SubjectID")) %>% 
  inner_join(adni_mri %>% dplyr::select(SubjectID, BL_Mag:BL_HippVol), by=c("SubjID" = "SubjectID")) %>% # select baseline magnetic field, ICV and hippocampal volume
  inner_join(comTrait %>% rownames_to_column("SubjID") %>% dplyr::select(SubjID, Educ, Age, Sex)) %>% # 42 subjects have available data, all subjects have gene expression data from baseline
  rename_with(~ str_replace(., "BL_", ""), contains("BL_")) %>% 
  pivot_longer(cols = 2:4, names_to = "Subfield", values_to = "Measures") %>% 
  named_group_split(Subfield) %>% 
  map(~ lm(paste0("ME3 ~ Measures + Age + Sex + Educ + ICV"), data = .x)) %>% 
  map_dfr(broom::tidy, conf.int = TRUE, .id = "Subfield") %>%
  filter(term == "Measures") %>% 
  write.csv(paste0(outdir, "/1_2/ADNI_HP_Subfield.csv"), row.names = F) # these results are now in supplementary table 10b
