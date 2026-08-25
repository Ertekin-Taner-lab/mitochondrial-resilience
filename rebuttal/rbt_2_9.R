# Author: Wei (Adelyn) Tsai; tsai.wei@mayo.edu
# If you are using this script, please cite our study
# This script addresses comment 9 about hippocampal volume from reviewer 2

suppressPackageStartupMessages(library(tidyverse)) #v2.0.0
suppressPackageStartupMessages(library(ggfortify)) #v0.4.19
suppressPackageStartupMessages(library(openxlsx)) #v4.2.8
suppressPackageStartupMessages(library(broom)) #v1.0.12
suppressPackageStartupMessages(library(patchwork)) #v1.3.2
suppressPackageStartupMessages(library(meta)) #v6.5
suppressPackageStartupMessages(library(pbapply)) #v1.7


# list of input/output directories
indir <- "./indir"
residir <- "./outdir/resids"
blwgcna_dir <- "./outdir/bl_wgcna"
nduf_assoc_dir <- "./outdir/nduf_assoc_de"

# output directory
outdir <- "./Rebuttal/outdir"

source("./Codes/fns.R") #source a list of frequently used functions

# Prepare dataframe for associations ----
comTrait <- readRDS(paste0(indir,"/comTrait.rds")) # a file having cross-sectional endophenotypes for MCSA (from row 1-105) and ADNI (from row 106-196); row names are donorID
adni_harm_hv <- read_csv(paste0(indir, "/adni_harm_hv.csv")) #this is a file containing harmonized ADNI hippocampal volume data harmonized to MCSA's
# Relevant columns include:
# ptid = SubjectID
# fieldstrength: 1.5 or 3T
# hpvol: harmonized hippocampal volume
# tiv: total intracranial volume
MEs1 <- read.csv(paste0(blwgcna_dir, "/MEs1_SP12.csv")) # from 2_bl_wgcna.R
nduf_bb <- read.csv(paste0(nduf_assoc_dir, "/nduf_bb.csv")) #generated in 4_nduf_assoc_de.R

# redo analysis for ADNI using volume harmonized to MCSA's
# read in microarray file and organize the microarray file to contain only NDUFs and rename the columns from gene id to gene name
adni_mc <- readRDS(paste0(indir,"/ADNI_s91_microarray.rds")) #microarray file for 91 donors from ADNI. This has 10116 protein-coding genes, gene id as column names, donorID as row names
adni_covars <- read_delim(paste0(residir,"/ADNI_s91_covars.txt"), delim = "\t") #covariates file for 91 donors from ADNI. Rows are donors columns are covariates
mc_ndufs <- adni_mc %>% dplyr::select(nduf_bb$gene_id)
mc_ndufs <- mc_ndufs %>% rename_at(vars(1:(ncol(.))), ~ nduf_bb$gene_name[match(names(mc_ndufs)[1:(ncol(mc_ndufs))], nduf_bb$gene_id)]) 

pdf(paste0(outdir, "/2_9/ADNI_HV_distribution.pdf"), width = 4, height = 4)
adni_anadf %>% 
  dplyr::select(SubjectID, HippVol, hpvol) %>% 
  pivot_longer(cols = -SubjectID, names_to = "DataType", values_to = "Volume") %>% 
  mutate(DataType = ifelse(DataType == "HippVol", "Unharmonized", "Harmonized")) %>% 
  mutate(DataType = factor(DataType, levels = c("Unharmonized", "Harmonized"))) %>% 
  ggplot(aes(x = Volume, fill = DataType, color = DataType)) +
  geom_histogram(position="identity", alpha=0.5, linewidth = 0.1) + 
  theme_classic(base_size = 6.5) +
  theme(legend.key.width = unit(3, "mm"),
        legend.key.height = unit(4, "mm"))
dev.off() # this is in the first rebuttal Rebuttal Fig.24a


pdf(paste0(outdir, "/2_9/old_new_volume_compare.pdf"), width = 8, height = 4)
adni_anadf %>% 
  ggplot(aes(x = HippVol, y = hpvol)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Original HippVol", y = "Harmonized HippVol", title = "HippVol comparison") +
  theme_classic(base_size = 6.5) +
  adni_anadf %>% 
    ggplot(aes(x = Hipp_ICV, y = tiv)) +
    geom_point() +
    geom_smooth(method = "lm") +
    labs(x = "Original ICV", y = "Harmonized ICV", title = "ICV comparison") +
    theme_classic(base_size = 6.5)
dev.off()# this is in the first rebuttal Rebuttal Fig.24b


# conduct analysis using standardized volume in ADNI----
## give the change in scale in harmonized data, we standardized the dependent variable (M3/NDUFs) and independent variable (hpvol)
## so we can compare between estimates
standardize_FUN <- function(x) {
  
  x_med <- median(x, na.rm = TRUE)
  x_mad <- mad(x, na.rm = TRUE)
  x <- ((x - x_med) / x_mad)
  x
  
}
adni_anadf_std <- comTrait[106:196, ] %>% 
  rownames_to_column("SubjectID") %>% 
  dplyr::select(SubjectID, Age, Sex, Educ, HippVol, Hipp_ICV, Mag_HippICV) %>% 
  left_join(adni_covars %>% 
              dplyr::select(SubjectID, RIN, AffyPlate.meanREE, SITE.meanREE)) %>% 
  left_join(adni_harm_hv %>% 
              dplyr::select(ptid, fieldstrength, hpvol, tiv) %>% 
              mutate(fieldstrength = factor(fieldstrength)), by=c("SubjectID" = "ptid")) %>% 
  filter(!is.na(HippVol)) %>% # there're 4 people with data in the harmonized dataset but not in the original dataset, 
  # so we filtered to retain those that have the original HV data to allow fair comparison between the old and new results
  left_join(MEs1 %>% dplyr::select(SubjectID, ME3)) %>% 
  left_join(mc_ndufs %>% tibble::rownames_to_column("SubjectID")) %>% 
  mutate_at(.vars = c("hpvol", "ME3", nduf_bb$gene_name), .funs = standardize_FUN) #hpvol = harmonized HV

adni_stdharm_res <- list()
for (i in c("ME3", nduf_bb$gene_name)){
  if (i == "ME3"){
    adni_stdharm_res[[i]] <- tidy(lm(ME3 ~ hpvol + tiv + fieldstrength + Age + Sex + Educ, data=adni_anadf_std), conf.int = TRUE) %>% filter(term == "hpvol")
  }
  else{
    adni_stdharm_res[[i]] <- tidy(lm(paste0(i, " ~ hpvol + tiv + fieldstrength + Age + Sex + Educ + RIN + AffyPlate.meanREE + SITE.meanREE"), data=adni_anadf_std), conf.int = TRUE) %>% filter(term == "hpvol")
  }
}
adni_stdharm_res = bind_rows(adni_stdharm_res, .id = "ME_NDUF")

# meta-analyzed results from standardized harmonized volume ----
mcsa_stdres <- read.xlsx("./Rebuttal/outdir/1_12/std_assoc.xlsx", sheet = "MCSA_Combined") # from rbt_1_12.R
mcsa_stdres_use <- mcsa_stdres %>% 
  filter(Phenotype == "HippVol" & Type %in% c("Bl_M3", nduf_bb$gene_name)) %>% 
  mutate(Type = ifelse(Type == "Bl_M3", "ME3", Type)) %>% 
  rename(ME_NDUF = Type)
# run meta-analysis with new est and std.error
meta_df = inner_join(mcsa_stdres_use %>% dplyr::select(-Phenotype), adni_stdharm_res %>% dplyr::select(-term),
                     by = c("ME_NDUF"), suffix = c("_MCSA", "_ADNI")) %>% 
  column_to_rownames("ME_NDUF")
# this indicates column indices corresponding to estimates and standard error of MCSA and ADNI
eff_idx1 = which(colnames(meta_df) == "estimate_MCSA")
eff_idx2 = which(colnames(meta_df) == "estimate_ADNI")
err_idx1 = which(colnames(meta_df) == "std.error_MCSA")
err_idx2 = which(colnames(meta_df) == "std.error_ADNI")
# run meta-analysis using the function meta_gen_fn defined above
meta_raw = meta_df %>% pbapply(1, meta_gen_fn, eff_idx1 = eff_idx1, eff_idx2 = eff_idx2, err_idx1 = err_idx1, err_idx2 = err_idx2)
# Get all fixed-effect results
fix_res = meta_raw %>%
  pbsapply(function(x) x$fixed) %>%
  t() %>%
  as.data.frame() %>%
  rename_with(~paste0(.x, "_fixed")) %>%
  rownames_to_column(var = "ME_NDUF") %>%
  mutate_all(.funs = ~unlist(.x)) %>% 
  mutate_all(.funs = ~unname(.x))
# Get all random-effect results
rand_res = meta_raw %>%
  pbsapply(function(x) x$random) %>%
  t() %>%
  data.frame() %>%
  rename_with(~paste0(.x, "_random")) %>%
  rownames_to_column(var = "ME_NDUF") %>%
  mutate_all(.funs = ~unlist(.x)) %>%
  mutate_all(.funs = ~unname(.x))
# Get heterogenity meassurments Q, I2, tau2
het_res = meta_raw %>%
  map(~.x[c("Q","df.Q","pval.Q","I2", "lower.I2", "upper.I2", "tau2", "se.tau2", "H", "lower.H", "upper.H")]) %>%
  bind_rows(.id = "ME_NDUF")
# Join the results to the original regression coeffcients
meta_df = meta_df %>%
  rownames_to_column(var = "ME_NDUF") %>% 
  inner_join(fix_res) %>%
  inner_join(rand_res) %>% 
  rename_with(~sub("^TE", "beta", .x), everything()) %>%
  rename_with(~sub("^seTE", "se", .x), everything()) %>%
  inner_join(het_res, by = "ME_NDUF")
# final meta-analyze dataframe, including adjusting for multiple test correction, determining meta-analyzed beta, p and q-values based on heterogeneity measures I2
stdharm_hv_meta = meta_df %>%
  mutate(q_random = p.adjust(p_random, method = "fdr"),
         q_fixed = p.adjust(p_fixed, method = "fdr"),
         model_meta = ifelse(I2 < 0.25, "Fixed", "Random"),
         beta_meta = ifelse(model_meta =="Fixed", beta_fixed, beta_random),
         se_meta = ifelse(model_meta =="Fixed", se_fixed, se_random),
         p_meta = ifelse(model_meta =="Fixed", p_fixed, p_random),
         q_meta = ifelse(model_meta =="Fixed", q_fixed, q_random),
         CI_H_meta = ifelse(model_meta =="Fixed", upper_fixed, upper_random),
         CI_L_meta = ifelse(model_meta == "Fixed", lower_fixed, lower_random)) %>%
  dplyr::select(-c(ends_with("_random"), ends_with("_fixed"))) 

# Plot standardized effect size and 95%CI for results using unharmonized vs harmonized HV data----
meta_m3_stdres_og <- read.xlsx("./Rebuttal/outdir/1_12/std_assoc.xlsx", sheet = "Bl_M3_meta") # from rbt_1_12.R
meta_ndufs_stdres_og <- read.xlsx("./Rebuttal/outdir/1_12/std_assoc.xlsx", sheet = "Bl_NDUFs_meta") # from rbt_1_12.R
## original meta-analyzed results
meta_og_use <- meta_m3_stdres_og %>% 
  filter(Phenotype == "HippVol") %>% 
  mutate(ME_NDUF = "ME3", .before = "Phenotype") %>% 
  rbind(meta_ndufs_stdres_og %>% 
          filter(Phenotype == "HippVol") %>% 
          rename(ME_NDUF = Gene) %>% 
          dplyr::select(-q.value_MCSA, -q.value_ADNI, -q_meta)) %>% 
  dplyr::select(-Phenotype) %>% 
  mutate(Type = "Unharmonized")
## combining original and new meta-analyzed results
meta_combine <- meta_og_use %>% 
  arrange(ME_NDUF) %>% 
  relocate(CI_L_meta, .before = "CI_H_meta") %>% 
  rbind(stdharm_hv_meta %>% 
          arrange(ME_NDUF) %>% 
          dplyr::select(-q_meta) %>% mutate(Type = "Harmonized"))
write.xlsx(list("ADNI" = adni_stdharm_res, "Meta" = meta_combine), paste0(outdir, "/2_9/hv_harmonization_check.xlsx")) 
# these are in the first rebuttal Rebuttal Table 36 and 37.

## forest plot
pdf(paste0(outdir, "/2_9/forestp_est_compare.pdf"), width = 5, height = 6)
meta_combine %>% 
  dplyr::select(Type, ME_NDUF, estimate_ADNI, p.value_ADNI, conf.low_ADNI, conf.high_ADNI) %>% 
  rename_with(~sub("_ADNI", "", .x), everything()) %>% 
  mutate(Cohort = "ADNI") %>%
  rbind(
    meta_combine %>% 
      dplyr::select(Type, ME_NDUF, beta_meta, p_meta, CI_L_meta, CI_H_meta) %>% 
      dplyr::rename(
        estimate = beta_meta,
        p.value = p_meta,
        conf.low = CI_L_meta,
        conf.high = CI_H_meta
      ) %>% 
      mutate(Cohort = "Meta-analyzed")
  ) %>% 
  mutate(Cohort = factor(Cohort, levels = c("ADNI", "Meta-analyzed"))) %>% 
  mutate(Type = factor(Type, levels = c("Unharmonized", "Harmonized"))) %>% 
  mutate(ME_NDUF = ifelse(ME_NDUF == "ME3", "Bl_M3", ME_NDUF)) %>%
  mutate(sig = ifelse(p.value < 0.05, "Nominal", "NS")) %>% 
  ggplot(aes(x=estimate, y=ME_NDUF)) + 
  geom_errorbar(aes(xmax = conf.high, xmin = conf.low, color = Type), size = 0.5, width = 0.3, position=position_dodge(width = 0.5)) + 
  geom_point(aes(color = Type, shape = sig), size = 2, position=position_dodge(width = 0.5)) +
  scale_y_discrete(limits = rev) + 
  scale_shape_manual(values = c(NS = 4, Nominal = 16)) +
  facet_grid(.~Cohort, scales = "free", space = "free") +
  labs(y = "") +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic"))
dev.off() # this is in the first rebuttal Rebuttal Fig.25