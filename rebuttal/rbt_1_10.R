# Author: Wei (Adelyn) Tsai; tsai.wei@mayo.edu
# If you are using this script, please cite our study
# This script addresses comment 10 about genetic risk interaction from reviewer 1

suppressPackageStartupMessages(library(tidyverse)) #v2.0.0
suppressPackageStartupMessages(library(broom)) #v1.0.12
suppressPackageStartupMessages(library(broom.mixed)) #v0.2.9
suppressPackageStartupMessages(library(ggeffects)) #v1.7
suppressPackageStartupMessages(library(patchwork)) #v1.3.2
suppressPackageStartupMessages(library(openxlsx)) #v4.2.8
suppressPackageStartupMessages(library(lmerTest)) #v3.1
suppressPackageStartupMessages(library(meta)) #v6.5
suppressPackageStartupMessages(library(pbapply)) #v1.7

source("./Codes/fns.R") #source a list of frequently used functions

#######################################

# list of input directories where files are
indir <- "./indir"
residir <- "./outdir/resids"
blwgcna_dir <- "./outdir/bl_wgcna"
nduf_assoc_dir <- "./outdir/nduf_assoc_de"

# output directory
outdir <- "./Rebuttal/outdir"

#######################################

# prepare files
comTrait <- readRDS(paste0(indir,"/comTrait.rds")) # a file having cross-sectional endophenotypes for MCSA (from row 1-105) and ADNI (from row 106-196); row names are donorID
MEs1 <- read.csv(paste0(blwgcna_dir, "/MEs1_SP12.csv")) #from 2_blwgcna.R
comRes <- readRDS(paste0(residir, "/comRes.rds")) # from 1_CreateResids.R, 
nduf_bb <- read.csv(paste0(nduf_assoc_dir,"/nduf_bb.csv")) # generated in 4_nduf_assoc_de.R; also in supplementary table 12
mcsa_covars <- read_delim(paste0(indir,"/R01_Resilience_PAXgene_Covars.txt"), delim = "\t") %>% mutate(SubjectID=as.character(ptnum)) #covariates file for MCSA, rows are donor columns are covariates
adni_covars <- read_delim(paste0(residir,"/ADNI_s91_covars.txt"), delim = "\t") #covariates file for 91 donors from ADNI. Rows are donors columns are covariates

# prepare files for MCSA
mcsa_df  <- comTrait[1:105,] %>% 
  rownames_to_column("SubjectID") %>% 
  dplyr::select(SubjectID, Memory, LMDR, HippVol, Hipp_ICV, Mag_HippICV, Sex, Educ, Age) %>% 
  left_join(mcsa_covars %>% 
              dplyr::select(SubjectID, APOE4) 
              ) %>% 
  left_join(MEs1 %>% dplyr::select(SubjectID, ME3)) %>% 
  left_join(
    comRes[1:105, ] %>% as.data.frame() %>% rownames_to_column("SubjectID") %>% dplyr::select(SubjectID, any_of(nduf_bb$gene_id))
  )
# rename geneID to gene name of NDUFs
mcsa_df  <- mcsa_df  %>% rename_at(vars(12:(ncol(.))), ~ nduf_bb$gene_name[match(names(mcsa_df)[12:(ncol(mcsa_df))], nduf_bb$gene_id)])
# transform to long files for analysis
mcsa_df <- mcsa_df %>% 
  pivot_longer(cols = 11:ncol(.), names_to = "Gene", values_to = "resids") %>% 
  named_group_split(Gene)

# prepare files for ADNI
adni_df  <- comTrait[106:196,] %>% 
  rownames_to_column("SubjectID") %>% 
  dplyr::select(SubjectID, Memory, LMDR, HippVol, Hipp_ICV, Mag_HippICV, Sex, Educ, Age) %>% 
  left_join(adni_covars %>% 
              dplyr::select(SubjectID, APOE4)) %>% 
  left_join(MEs1 %>% dplyr::select(SubjectID, ME3)) %>% 
  left_join(
    comRes[106:196, ] %>% as.data.frame() %>% rownames_to_column("SubjectID") %>% dplyr::select(SubjectID, any_of(nduf_bb$gene_id))
  )
# rename geneID to gene name of NDUFs
adni_df  <- adni_df  %>% rename_at(vars(12:(ncol(.))), ~ nduf_bb$gene_name[match(names(adni_df)[12:(ncol(adni_df))], nduf_bb$gene_id)])
# transform to long files for analysis
adni_df <- adni_df %>% 
  pivot_longer(cols = 11:ncol(.), names_to = "Gene", values_to = "resids") %>% 
  named_group_split(Gene)

#######################################

#######################################

# MCSA: run interaction analysis
# We can't run hippocampal volume analysis for MCSA because only 1 person with APOE4/E4 has the data
mcsa_fit <- list()
for (i in c("LMDR", "Memory")){
  df <- mcsa_df %>%
      map(~ lm(paste0(i, "~ resids*APOE4 + Age + Sex + Educ"), data = .x)) %>%
      map_dfr(broom::tidy, conf.int = TRUE, .id = "Gene") %>%
      filter(term %in% c("resids", "APOE4", "resids:APOE4")) %>%
      group_by(term) %>%
      mutate(q.value = p.adjust(p.value, method = "fdr")) %>% 
      ungroup()
    
  mcsa_fit[[i]] <- df
  rm(df, i)
}
# merge results
mcsa_fit_df <- bind_rows(mcsa_fit, .id="Phenotype")
#save models for plotting
mcsa_models <- list()
for (i in c("LMDR", "Memory")){
  m <- mcsa_df %>%
    map(~ lm(paste0(i, "~ resids*APOE4 + Age + Sex + Educ"), data = .x))
  
  mcsa_models[[i]] <- m
  rm(m, i)
}

# ADNI: run interaction analysis
adni_fit <- list()
for (i in c("LMDR", "Memory", "HippVol")){
  if(i != "HippVol"){
    df <- adni_df %>%
      map(~ lm(paste0(i, "~ resids*APOE4 + Age + Sex + Educ"), data = .x)) %>%
      map_dfr(broom::tidy, conf.int = TRUE, .id = "Gene") %>%
      filter(term %in% c("resids", "APOE4", "resids:APOE4")) %>%
      group_by(term) %>% 
      mutate(q.value = p.adjust(p.value, method = "fdr"))  %>% 
      ungroup()
    
    adni_fit[[i]] <- df
  }else{
    df <- adni_df %>%
      map(~ lm(paste0(i, "~ resids*APOE4 + Age + Sex + Educ + Hipp_ICV + Mag_HippICV"), data = .x)) %>%
      map_dfr(broom::tidy, conf.int = TRUE, .id = "Gene") %>%
      filter(term %in% c("resids", "APOE4", "resids:APOE4")) %>%
      group_by(term) %>%
      mutate(q.value = p.adjust(p.value, method = "fdr"))  %>% 
      ungroup()
    
    adni_fit[[i]] <- df
  }
  rm(df, i)
}
# merge_results
adni_fit_df <- bind_rows(adni_fit, .id="Phenotype")
#save models for plotting
adni_models <- list()
for (i in c("LMDR", "Memory", "HippVol")){
  if(i != "HippVol"){
    m <- adni_df %>%
      map(~ lm(paste0(i, "~ resids*APOE4 + Age + Sex + Educ"), data = .x))
    
    adni_models[[i]] <- m
  }else{
    m <- adni_df %>%
      map(~ lm(paste0(i, "~ resids*APOE4 + Age + Sex + Educ + Hipp_ICV + Mag_HippICV"), data = .x))
    
    adni_models[[i]] <- m
  }
  rm(m, i)
}
save(mcsa_models, adni_models, file = paste0(outdir, "/1_10/models.RData"))

#######################################

# meta-analyze interaction analysis results from MCSA and ADNI

# prepare files
options(stringsAsFactors = FALSE)
pboptions(type = "timer")
meta_df = inner_join(mcsa_fit_df %>% filter(term == "resids:APOE4"), adni_fit_df %>% filter(term == "resids:APOE4"), by = c("Phenotype", "Gene"), suffix = c("_MCSA", "_ADNI")) %>% 
  mutate(Pheno_NDUF = paste0(Phenotype, "_", Gene)) %>% 
  column_to_rownames("Pheno_NDUF")

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
  rownames_to_column(var = "Pheno_Gene") %>%
  mutate_all(.funs = ~unlist(.x)) %>% 
  mutate_all(.funs = ~unname(.x))

# Get all random-effect results
rand_res = meta_raw %>%
  pbsapply(function(x) x$random) %>%
  t() %>%
  data.frame() %>%
  rename_with(~paste0(.x, "_random")) %>%
  rownames_to_column(var = "Pheno_Gene") %>%
  mutate_all(.funs = ~unlist(.x)) %>%
  mutate_all(.funs = ~unname(.x))

# Get heterogenity meassurments Q, I2, tau2
het_res = meta_raw %>%
  map(~.x[c("Q","df.Q","pval.Q","I2", "lower.I2", "upper.I2", "tau2", "se.tau2", "H", "lower.H", "upper.H")]) %>%
  bind_rows(.id = "Pheno_Gene")

# Join the results to the original regression coeffcients
meta_df = meta_df %>%
  rownames_to_column(var = "Pheno_Gene") %>% 
  inner_join(fix_res) %>%
  inner_join(rand_res) %>% 
  rename_with(~sub("^TE", "beta", .x), everything()) %>%
  rename_with(~sub("^seTE", "se", .x), everything()) %>%
  inner_join(het_res, by = "Pheno_Gene")


# final meta-analyze dataframe, including adjusting for multiple test correction, determining meta-analyzed beta, p and q-values based on heterogeneity measures I2
meta_final = meta_df %>%
  group_by(Phenotype) %>%
  mutate(q_random = p.adjust(p_random, method = "fdr"),
         q_fixed = p.adjust(p_fixed, method = "fdr"),
         model_meta = ifelse(I2 < 0.25, "Fixed", "Random"),
         beta_meta = ifelse(model_meta =="Fixed", beta_fixed, beta_random),
         se_meta = ifelse(model_meta =="Fixed", se_fixed, se_random),
         p_meta = ifelse(model_meta =="Fixed", p_fixed, p_random),
         q_meta = ifelse(model_meta =="Fixed", q_fixed, q_random),
         CI_H_meta = ifelse(model_meta =="Fixed", upper_fixed, upper_random),
         CI_L_meta = ifelse(model_meta == "Fixed", lower_fixed, lower_random)) %>%
  dplyr::select(-c(ends_with("_random"), ends_with("_fixed"))) %>% 
  dplyr::select(-Pheno_Gene)


rm(fix_res, rand_res, het_res, meta_df, meta_raw, eff_idx1, eff_idx2, err_idx1, err_idx2)

write.xlsx(list("Meta" = meta_final, "MCSA" = mcsa_fit_df, "ADNI" = adni_fit_df), paste0(outdir,"/1_10/APOE4_NDUF_interaction.xlsx"))
# significant results are shown in the first rebuttal Rebuttal Table 21 & 22

#######################################

#######################################
# Plot significant interaction results

p1 <- predict_response(mcsa_models$Memory$NDUFS5, terms = c("resids", "APOE4")) %>% 
  plot()  + labs(title = "NDUFS5*APOE4_dose on memory (MCSA)", x = "residuals") + theme(text = element_text(size = 5), axis.text = element_text(size = 5))
p2 <- predict_response(adni_models$Memory$NDUFS5, terms = c("resids", "APOE4")) %>% 
  plot()  + labs(title = "NDUFS5*APOE4_dose on memory (ADNI)", x = "residuals") + theme(text = element_text(size = 5), axis.text = element_text(size = 5))
p3 <- predict_response(adni_models$HippVol$NDUFAF2, terms = c("resids", "APOE4")) %>% 
  plot() + labs(title = "NDUFAF2*APOE4_dose on hippocampal volume (ADNI)", x = "residuals") + theme(text = element_text(size = 5), axis.text = element_text(size = 5))

pdf(paste0(outdir, "/1_10/interaction_p.pdf"), width = 6.2, height = 2)
p1 + p2 + p3 + plot_layout(guides = "collect") & theme(legend.key.size = unit(0.1, "in"))
dev.off() # this is in the first rebuttal Rebuttal Fig.14
#######################################