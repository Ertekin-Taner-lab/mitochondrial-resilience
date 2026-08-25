# Author: Wei (Adelyn) Tsai; tsai.wei@mayo.edu
# If you are using this script, please cite our study
# This script addresses comment 12 about effect size reporting and forest plots from reviewer 1
# we presented standardized effect size, 95% CI and forestplot for:
#  1) Blood M2, M3 and M30 with memory, LMDR, hippocample volume
#  2) Blood NDUFs with memory, LMDR, HippVol
#  3) NDUFB9 longitudinal results

suppressPackageStartupMessages(library(tidyverse)) #v2.0.0
suppressPackageStartupMessages(library(lmerTest)) #v3.1
suppressPackageStartupMessages(library(openxlsx)) #v4.2.8
suppressPackageStartupMessages(library(broom)) #v1.0.12
suppressPackageStartupMessages(library(broom.mixed)) #v0.2.9
suppressPackageStartupMessages(library(pbapply)) #v1.7
suppressPackageStartupMessages(library(patchwork)) #v1.3.2



source("./Codes/fns.R") #source a list of frequently used functions

# function for standardization
standardize_FUN <- function(x) {
  
  x_med <- median(x, na.rm = TRUE)
  x_mad <- mad(x, na.rm = TRUE)
  x <- ((x - x_med) / x_mad)
  x
  
}

# list of input/output directories
indir <- "./indir"
blwgcna_dir <- "./outdir/bl_wgcna"
brwgcna_dir <- "./outdir/br_wgcna"
residir <- "./outdir/resids"
nduf_assoc_dir <- "./outdir/nduf_assoc_de"

outdir <- "./Rebuttal/outdir"

# 1) Blood M2, M3, M30 with memory, LMDR, hippocample volume ----
# prepare data frame for analysis
comTrait <- readRDS(paste0(indir,"/comTrait.rds")) # a file having cross-sectional endophenotypes for MCSA (from row 1-105) and ADNI (from row 106-196); row names are donorID
MEs1 <- read.csv(paste0(blwgcna_dir, "/MEs1_SP12.csv")) # from 2_bl_wgcna.R
analysis_df <- comTrait %>% 
  rownames_to_column("SubjectID") %>% 
  left_join(MEs1 %>% dplyr::select(SubjectID, ME2, ME3, ME30)) %>% 
  mutate(Mag_HippICV = factor(Mag_HippICV))
# standardize the dependent and independent variables
analysis_df[1:105,] <- analysis_df[1:105, ] %>% 
  mutate_at(.vars = c("Memory", "LMDR", "HippVol", "ME2", "ME3", "ME30"), .funs = standardize_FUN)
analysis_df[106:196,] <- analysis_df[106:196, ] %>% 
  mutate_at(.vars = c("Memory", "LMDR", "HippVol", "ME2", "ME3", "ME30"), .funs = standardize_FUN)

# run regression to get new standardized estimates
mcsa_me <- list()
for (i in c("Memory", "LMDR", "HippVol")){
  for (j in c("ME2", "ME3", "ME30")){
    if (i != "HippVol"){
      mcsa_me[[i]][[j]] <- tidy(lm(paste0(j, " ~", i, " + Age + Sex + Educ"), data=analysis_df[1:105, ]), conf.int = TRUE) %>% filter(term == i)
    }
    else{
      mcsa_me[[i]][[j]] <- tidy(lm(paste0(j, " ~", i, " + Age + Sex + Educ + Hipp_ICV"), data=analysis_df[1:105, ]), conf.int = TRUE) %>% filter(term == i)
    }
  }
}
adni_me <- list()
for (i in c("Memory", "LMDR", "HippVol")){
  for (j in c("ME2", "ME3", "ME30")){
      if (i != "HippVol"){
      adni_me[[i]][[j]] <- tidy(lm(paste0(j, " ~", i, " + Age + Sex + Educ"), data=analysis_df[106:196, ]), conf.int = TRUE) %>% filter(term == i)
    }
    else{
      adni_me[[i]][[j]] <- tidy(lm(paste0(j, " ~", i, " + Age + Sex + Educ + Hipp_ICV + Mag_HippICV"), data=analysis_df[106:196, ]), conf.int = TRUE) %>% filter(term == i)
    }
  }
}
#bind results
mcsa_me_dfs <- mcsa_me %>% map(~ bind_rows(.x, .id = "Module")) %>% bind_rows(.id = "Phenotype")
adni_me_dfs <- adni_me %>% map(~ bind_rows(.x, .id = "Module")) %>% bind_rows(.id = "Phenotype")
rm(i, j)

# do meta-analysis using new est and std.error
meta_df = inner_join(mcsa_me_dfs, adni_me_dfs, by = c("Phenotype", "Module", "term"), suffix = c("_MCSA", "_ADNI")) %>% 
  mutate(Pheno_Module = paste0(Phenotype, "_", Module)) %>% 
  column_to_rownames("Pheno_Module")
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
  rownames_to_column(var = "Pheno_Module") %>%
  mutate_all(.funs = ~unlist(.x)) %>% 
  mutate_all(.funs = ~unname(.x))
# Get all random-effect results
rand_res = meta_raw %>%
  pbsapply(function(x) x$random) %>%
  t() %>%
  data.frame() %>%
  rename_with(~paste0(.x, "_random")) %>%
  rownames_to_column(var = "Pheno_Module") %>%
  mutate_all(.funs = ~unlist(.x)) %>%
  mutate_all(.funs = ~unname(.x))
# Get heterogenity meassurments Q, I2, tau2
het_res = meta_raw %>%
  map(~.x[c("Q","df.Q","pval.Q","I2", "lower.I2", "upper.I2", "tau2", "se.tau2", "H", "lower.H", "upper.H")]) %>%
  bind_rows(.id = "Pheno_Module")
# Join the results to the original regression coeffcients
meta_df = meta_df %>%
  rownames_to_column(var = "Pheno_Module") %>% 
  inner_join(fix_res) %>%
  inner_join(rand_res) %>% 
  rename_with(~sub("^TE", "beta", .x), everything()) %>%
  rename_with(~sub("^seTE", "se", .x), everything()) %>%
  inner_join(het_res, by = "Pheno_Module")
# final meta-analyze dataframe, including adjusting for multiple test correction, determining meta-analyzed beta, p and q-values based on heterogeneity measures I2
me_meta = meta_df %>%
  mutate(model_meta = ifelse(I2 < 0.25, "Fixed", "Random"),
         beta_meta = ifelse(model_meta =="Fixed", beta_fixed, beta_random),
         se_meta = ifelse(model_meta =="Fixed", se_fixed, se_random),
         p_meta = ifelse(model_meta =="Fixed", p_fixed, p_random),
         CI_H_meta = ifelse(model_meta =="Fixed", upper_fixed, upper_random),
         CI_L_meta = ifelse(model_meta == "Fixed", lower_fixed, lower_random)) %>%
  dplyr::select(-c(ends_with("_random"), ends_with("_fixed"))) %>% 
  dplyr::select(-term, -Pheno_Module)

rm(fix_res, rand_res, het_res, meta_df, meta_raw, eff_idx1, eff_idx2, err_idx1, err_idx2)


# 2) Blood NDUFs with memory, LMDR, HippVol ----
# prepare files
# MCSA
nduf_bb <- read.csv(paste0(nduf_assoc_dir, "/nduf_bb.csv")) #generated in 4_nduf_assoc_de.R
cqn_raw <- read_delim(paste0(indir, "/PaX108_R01resilience_gene_CQN_neg3_postQC_s105.txt"), delim = "\t") # MCSA cqn file for 18046 genes. Each row is a gene and each column is a donor
cqn_use <- cqn_raw %>% filter(GeneId %in% nduf_bb$gene_id)
mcsa_covars <- read_delim(paste0(indir,"/R01_Resilience_PAXgene_Covars.txt"), delim = "\t") %>% mutate(SubjectID=as.character(ptnum)) #covariates file for MCSA, rows are donor columns are covariates
# prepare a long-format input file for association
cqn_long <- cqn_use %>%
  dplyr::select(-Chromosome, -Start, -End, -Length, -GeneId, -GeneBiotype) %>%
  tidyr::pivot_longer(-GeneName, names_to = "SubjectID", values_to = "Expr") %>%
  inner_join(dplyr::select(mcsa_covars, SubjectID, RIN, PAXgene_flowcell, PAXgene_Batch), by = "SubjectID") %>%
  mutate(PAXgene_Batch = factor(PAXgene_Batch),
         PAXgene_flowcell = as.factor(PAXgene_flowcell))  %>%
  inner_join(comTrait %>% 
               tibble::rownames_to_column("SubjectID") %>% 
               dplyr::select(SubjectID, Sex, Educ, Age, LMDR, Memory, HippVol, Hipp_ICV), by="SubjectID") %>% 
  named_group_split(GeneName) %>% # split by gene name, each gene name contains a dataframe with 105 subjects
  # standardize dependent and independent variables
  map(~mutate_at(.x, .vars = c("Expr", "Memory", "LMDR", "HippVol"), .funs = standardize_FUN))
# run regression analysis to get new estimates
mcsa_ndufs_ls <- list()
for (i in c("LMDR", "Memory", "HippVol")){
  if (i != "HippVol"){
    mcsa_ndufs_ls[[i]] <- cqn_long %>%
      map(~ lm(paste0("Expr ~ ", i, "+ Age + Sex + Educ + RIN + PAXgene_flowcell + PAXgene_Batch"), data = .x)) %>%
      map_dfr(broom::tidy, conf.int = TRUE, .id = "Gene") %>%
      filter(term == i)
  }else{
    mcsa_ndufs_ls[[i]] <- cqn_long %>%
      map(~ lm(paste0("Expr ~ ", i, "+ Age + Sex + Educ + Hipp_ICV + RIN + PAXgene_flowcell + PAXgene_Batch"), data = .x)) %>%
      map_dfr(broom::tidy, conf.int = TRUE, .id = "Gene") %>%
      filter(term == i)
  }
}
mcsa_ndufs_df <- bind_rows(mcsa_ndufs_ls) %>% 
  group_by(term) %>% 
  mutate(q.value = p.adjust(p.value, method = "fdr"))
rm(i)

# ADNI
# prepare files
adni_mc <- readRDS(paste0(indir,"/ADNI_s91_microarray.rds")) #microarray file for 91 donors from ADNI. This has 10116 protein-coding genes, gene id as column names, donorID as row names
adni_covars <- read_delim(paste0(residir,"/ADNI_s91_covars.txt"), delim = "\t") #covariates file for 91 donors from ADNI. Rows are donors columns are covariates
mc_long <- adni_mc %>% dplyr::select(nduf_bb$gene_id)
mc_long <- mc_long %>% rename_at(vars(1:(ncol(.))), ~ nduf_bb$gene_name[match(names(mc_long)[1:(ncol(mc_long))], nduf_bb$gene_id)]) 
mc_use <- mc_long %>%
  tibble::rownames_to_column("SubjectID") %>% 
  filter(SubjectID %in% adni_covars$SubjectID) %>% 
  pivot_longer(2:21, names_to = "Gene", values_to = "Expr") %>% 
  inner_join(adni_covars %>% dplyr::select(SubjectID, RIN, AffyPlate.meanREE, SITE.meanREE), by="SubjectID") %>% 
  inner_join(comTrait %>% 
               tibble::rownames_to_column("SubjectID") %>% 
               dplyr::select(SubjectID, Sex, Educ, Age, LMDR, Memory, HippVol, Hipp_ICV, Mag_HippICV) %>%
               mutate(Mag_HippICV=factor(Mag_HippICV)), by="SubjectID") %>% 
  named_group_split(Gene) %>% # split by gene name, each gene name contains a dataframe with 91 subjects
  # standardize dependent and independent variables
  map(~mutate_at(.x, .vars = c("Expr", "Memory", "LMDR", "HippVol"), .funs = standardize_FUN))


#run regression analysis to get new estimates
adni_ndufs_ls <- list()
for (i in c("LMDR", "Memory", "HippVol")){
  if(i != "HippVol"){
    adni_ndufs_ls[[i]] <- mc_use %>%
      map(~ lm(paste0("Expr ~ ", i, "+ Age + Sex + Educ + RIN + AffyPlate.meanREE + SITE.meanREE"), data = .x)) %>%
      map_dfr(broom::tidy, conf.int = TRUE, .id = "Gene") %>%
      filter(term == i)
  }else{
    adni_ndufs_ls[[i]] <- mc_use %>%
      map(~ lm(Expr ~ HippVol + Age + Sex + Educ + Hipp_ICV + Mag_HippICV + RIN + AffyPlate.meanREE + SITE.meanREE, data = .x)) %>%
      map_dfr(broom::tidy, conf.int = TRUE, .id = "Gene") %>%
      filter(term == i)
  }
}
adni_ndufs_df <- bind_rows(adni_ndufs_ls)  %>% 
  group_by(term) %>% 
  mutate(q.value = p.adjust(p.value, method = "fdr"))

rm(i)


# run meta-analysis with new est and std.error
meta_df = inner_join(mcsa_ndufs_df, adni_ndufs_df, by = c("Gene", "term"), suffix = c("_MCSA", "_ADNI")) %>% 
  mutate(Pheno_NDUF = paste0(term, "_", Gene)) %>% 
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
  rownames_to_column(var = "Pheno_NDUF") %>%
  mutate_all(.funs = ~unlist(.x)) %>% 
  mutate_all(.funs = ~unname(.x))
# Get all random-effect results
rand_res = meta_raw %>%
  pbsapply(function(x) x$random) %>%
  t() %>%
  data.frame() %>%
  rename_with(~paste0(.x, "_random")) %>%
  rownames_to_column(var = "Pheno_NDUF") %>%
  mutate_all(.funs = ~unlist(.x)) %>%
  mutate_all(.funs = ~unname(.x))
# Get heterogenity meassurments Q, I2, tau2
het_res = meta_raw %>%
  map(~.x[c("Q","df.Q","pval.Q","I2", "lower.I2", "upper.I2", "tau2", "se.tau2", "H", "lower.H", "upper.H")]) %>%
  bind_rows(.id = "Pheno_NDUF")
# Join the results to the original regression coeffcients
meta_df = meta_df %>%
  rownames_to_column(var = "Pheno_NDUF") %>% 
  inner_join(fix_res) %>%
  inner_join(rand_res) %>% 
  rename_with(~sub("^TE", "beta", .x), everything()) %>%
  rename_with(~sub("^seTE", "se", .x), everything()) %>%
  inner_join(het_res, by = "Pheno_NDUF")
# final meta-analyze dataframe, including adjusting for multiple test correction, determining meta-analyzed beta, p and q-values based on heterogeneity measures I2
bl_ndufs_meta = meta_df %>%
  group_by(term) %>%
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
  dplyr::rename(Phenotype = term) %>% 
  dplyr::select(-Pheno_NDUF)


rm(fix_res, rand_res, het_res, meta_df, meta_raw, eff_idx1, eff_idx2, err_idx1, err_idx2)



# 3) NDUFB9 longitudinal results----
# prepare files
comRes <- readRDS(paste0(residir, "/comRes.rds")) # from 1_CreateResids.R

# load in longitudinal files of mcsa and adni
# for both files, rows are donors, and there are at least two or more rows per donor because there are data from different time point.
# columns are cognitive endophenotypes/hippcampal volume at each time point. 
# there is also a column of difference in years of each time point from the baseline. 
# and there are columns of demographic information including sex, education yrs, age at baseline, intracranial volume or MRI magnetic field
mcsa_long <- read.csv(paste0(indir,"/mcsa_long.csv")) %>% mutate(SubjectID = as.character(SubjectID))
adni_long <- read.csv(paste0(indir, "/adni_long.csv"))

# check for each subject, how many missing across visits there are for different endophenotypes 
mcsa_na_check = mcsa_long %>% group_by(SubjectID) %>% summarise_all(list(~sum(!is.na(.))))
adni_na_check = adni_long %>% group_by(SubjectID) %>% summarise_all(list(~sum(!is.na(.))))

# prepare files for association (MCSA)
mcsa_long_df <- mcsa_long %>% 
  dplyr::select(SubjectID, YrDiff, Sex, Educ, Age_bl, APOE4, LMDR, Memory, HippVol, ICV) %>% 
  left_join(comRes[1:105, ] %>% 
              as.data.frame() %>% 
              rownames_to_column("SubjectID") %>% 
              dplyr::select(SubjectID, ENSG00000147684) %>% 
              mutate(ENSG00000147684 = standardize_FUN(ENSG00000147684))) %>% 
  dplyr::rename(NDUFB9  = ENSG00000147684) 

# run regression analysis to get new estimates (MCSA)
mcsa_long_assoc <- list()
for (i in c("Memory", "LMDR", "HippVol")){
  print(i)
  
  if (i != "HippVol"){
    inputdf <- mcsa_long_df %>% 
      filter(! SubjectID %in% mcsa_na_check$SubjectID[which(mcsa_na_check[[i]]<=1)]) %>%  # only retain subjects with more than 1 longitudinal measures
      mutate_at(.vars = i, .funs = standardize_FUN)
    mcsa_long_assoc[[i]] <- tidy(lmer(paste0(i, "~", "NDUFB9*YrDiff + Age_bl + Sex + Educ + (1 + YrDiff|SubjectID)"), data=inputdf), conf.int = T) %>% 
      filter(term == "NDUFB9:YrDiff")
  }
  
  else {
    inputdf <- mcsa_long_df %>% 
      filter(! SubjectID %in% mcsa_na_check$SubjectID[which(mcsa_na_check[[i]]<=1)]) %>% 
      mutate_at(.vars = i, .funs = standardize_FUN)
    mcsa_long_assoc[[i]] <- tidy(lmer(paste0(i, "~", "NDUFB9*YrDiff + Age_bl + Sex + Educ + ICV + (1 + YrDiff|SubjectID)"), data=inputdf), conf.int = T) %>% 
      filter(term == "NDUFB9:YrDiff")
  }
}
# bind results
mcsa_long_res <- mcsa_long_assoc %>% bind_rows(.id="Phenotype") %>% dplyr::select(-group, -effect)

# prepare files for association (ADNI)
adni_long_df <- adni_long %>% 
  dplyr::select(SubjectID, YrDiff, Sex, Educ, Age_bl, APOE4, LMDR, Memory, HippVol, ICV, Mag) %>% 
  left_join(comRes[106:196, ] %>% 
              as.data.frame() %>% 
              rownames_to_column("SubjectID") %>% 
              dplyr::select(SubjectID, ENSG00000147684) %>% 
              mutate(ENSG00000147684 = standardize_FUN(ENSG00000147684))) %>% 
  dplyr::rename(NDUFB9 = ENSG00000147684)
# run regression analysis to get new estimates (ADNI)
adni_long_assoc <- list()
for (i in c("Memory", "LMDR", "HippVol")){
  print(i)
  
  if (i != "HippVol"){
    inputdf <- adni_long_df %>% 
      filter(! SubjectID %in% adni_na_check$SubjectID[which(adni_na_check[[i]]<=1)]) %>%  # only retain subjects with more than 1 longitudinal measures
      mutate_at(.vars = i, .funs = standardize_FUN)
    adni_long_assoc[[i]] <- tidy(lmer(paste0(i, "~", "NDUFB9*YrDiff + Age_bl + Sex + Educ + (1 + YrDiff|SubjectID)"), data=inputdf), conf.int = T) %>% 
      filter(term == "NDUFB9:YrDiff")
  }
  else {
    inputdf <- adni_long_df %>% 
      filter(! SubjectID %in% adni_na_check$SubjectID[which(adni_na_check[[i]]<=1)]) %>% 
      mutate_at(.vars = i, .funs = standardize_FUN)
    adni_long_assoc[[i]] <-  tidy(lmer(paste0(i, "~", "NDUFB9*YrDiff + Age_bl + Sex + Educ + ICV + Mag + (1 + YrDiff|SubjectID)"), data=inputdf), conf.int = T) %>% 
      filter(term == "NDUFB9:YrDiff")
  }
}
# bind results
adni_long_res <- adni_long_assoc %>% bind_rows(.id="Phenotype") %>% dplyr::select(-group, -effect)

# conduct meta-analysis on new estimates and std.error
meta_df = inner_join(mcsa_long_res, adni_long_res, by = c("Phenotype", "term"), suffix = c("_MCSA", "_ADNI")) %>% 
  column_to_rownames("Phenotype")
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
  rownames_to_column(var = "Phenotype") %>%
  mutate_all(.funs = ~unlist(.x)) %>% 
  mutate_all(.funs = ~unname(.x))
# Get all random-effect results
rand_res = meta_raw %>%
  pbsapply(function(x) x$random) %>%
  t() %>%
  data.frame() %>%
  rename_with(~paste0(.x, "_random")) %>%
  rownames_to_column(var = "Phenotype") %>%
  mutate_all(.funs = ~unlist(.x)) %>%
  mutate_all(.funs = ~unname(.x))
# Get heterogenity meassurments Q, I2, tau2
het_res = meta_raw %>%
  map(~.x[c("Q","df.Q","pval.Q","I2", "lower.I2", "upper.I2", "tau2", "se.tau2", "H", "lower.H", "upper.H")]) %>%
  bind_rows(.id = "Phenotype")
# Join the results to the original regression coeffcients
meta_df = meta_df %>%
  rownames_to_column(var = "Phenotype") %>% 
  inner_join(fix_res) %>%
  inner_join(rand_res) %>% 
  rename_with(~sub("^TE", "beta", .x), everything()) %>%
  rename_with(~sub("^seTE", "se", .x), everything()) %>%
  inner_join(het_res, by = "Phenotype")
# final meta-analyze dataframe, including adjusting for multiple test correction, determining meta-analyzed beta, p and q-values based on heterogeneity measures I2
long_meta = meta_df %>%
  mutate(q_random = p.adjust(p_random, method = "fdr"),
         model_meta = ifelse(I2 < 0.25, "Fixed", "Random"),
         beta_meta = ifelse(model_meta =="Fixed", beta_fixed, beta_random),
         se_meta = ifelse(model_meta =="Fixed", se_fixed, se_random),
         p_meta = ifelse(model_meta =="Fixed", p_fixed, p_random),
         CI_H_meta = ifelse(model_meta =="Fixed", upper_fixed, upper_random),
         CI_L_meta = ifelse(model_meta == "Fixed", lower_fixed, lower_random)) %>%
  dplyr::select(-c(ends_with("_random"), ends_with("_fixed"))) 


# save analysis results dataframes----
mcsa_combine_dfs <- mcsa_me_dfs %>% 
  rename(Type = Module) %>% 
  dplyr::select(-term) %>% 
  rbind(mcsa_ndufs_df %>% 
          rename(Type = Gene, Phenotype = term) %>% 
          dplyr::select(-q.value)) %>% 
  rbind(mcsa_long_res %>% 
          mutate(Type = "Longitudinal", .before="Phenotype") %>% 
          dplyr::select(-df, -term))
adni_combine_dfs <- adni_me_dfs %>% 
  rename(Type = Module) %>% 
  dplyr::select(-term) %>% 
  rbind(adni_ndufs_df %>% 
          rename(Type = Gene, Phenotype = term) %>% 
          dplyr::select(-q.value)) %>% 
  rbind(adni_long_res %>% 
          mutate(Type = "Longitudinal", .before="Phenotype") %>% 
          dplyr::select(-df, -term))


# save files to excel and plot
write.xlsx(list("MCSA_Combined" = mcsa_combine_dfs,
                "ADNI_Combined" = adni_combine_dfs,
                "Bl_ME_meta" = me_meta, 
                "Bl_NDUFs_meta" = bl_ndufs_meta, 
                "Longitudinal_meta" = long_meta), 
           paste0(outdir, "/1_12/std_assoc.xlsx")) # the organized results of these tables are in supplementary table 5
#######################################

#######################################
# make forest plot
# these are in extended data fig.4
pdf(paste0(outdir, "/1_12/BlME_std.pdf"), width = 3.8, height = 2)
me_meta %>% 
  dplyr::select(Phenotype, Module, estimate_MCSA, p.value_MCSA, conf.low_MCSA, conf.high_MCSA) %>% 
  rename_with(~sub("_MCSA", "", .x), everything()) %>% 
  mutate(Cohort = "MCSA") %>% 
  bind_rows(
    me_meta %>% 
      dplyr::select(Phenotype, Module, estimate_ADNI, p.value_ADNI, conf.low_ADNI, conf.high_ADNI) %>% 
      rename_with(~sub("_ADNI", "", .x), everything()) %>% 
      mutate(Cohort = "ADNI")
  ) %>% 
  bind_rows(
    me_meta %>% 
      dplyr::select(Phenotype, Module, beta_meta, p_meta, CI_L_meta, CI_H_meta) %>% 
      dplyr::rename(
        estimate = beta_meta,
        p.value = p_meta,
        conf.low = CI_L_meta,
        conf.high = CI_H_meta
      ) %>% 
      mutate(Cohort = "Meta-analyzed")
  ) %>% 
  mutate(Cohort = factor(Cohort, levels = c("ADNI", "MCSA", "Meta-analyzed")),
         Phenotype = factor(Phenotype, levels = c("Memory", "LMDR", "HippVol"))) %>% 
  ggplot(aes(x=estimate, y=Module)) + 
  geom_errorbar(aes(xmax = conf.high, xmin = conf.low, color = Cohort), size = 0.5, width = 0.3, position=position_dodge(width = 0.5)) + 
  geom_point(aes(color = Cohort), size = 2, position=position_dodge(width = 0.5)) +
  scale_color_manual(values = c("ADNI" = "green4", "MCSA" = "blue4", "Meta-analyzed" = "cyan2")) +
  labs(title = "Blood module eigengene-endophenotype association") +
  facet_grid(~Phenotype, scales = "free", space = "free") +
  scale_y_discrete(limits = rev) +
  theme_bw(base_size = 6.5)
dev.off()

pdf(paste0(outdir, "/1_12/BlNDUFs_std.pdf"), width = 4.5, height = 7)
bl_ndufs_meta %>% 
  dplyr::select(Gene, Phenotype, estimate_MCSA, p.value_MCSA, conf.low_MCSA, conf.high_MCSA) %>% 
  rename_with(~sub("_MCSA", "", .x), everything()) %>% 
  mutate(Cohort = "MCSA") %>% 
  bind_rows(
    bl_ndufs_meta %>% 
      dplyr::select(Gene, Phenotype, estimate_ADNI, p.value_ADNI, conf.low_ADNI, conf.high_ADNI) %>% 
      rename_with(~sub("_ADNI", "", .x), everything()) %>% 
      mutate(Cohort = "ADNI")
  ) %>% 
  bind_rows(
    bl_ndufs_meta %>% 
      dplyr::select(Gene, Phenotype, beta_meta, p_meta, CI_L_meta, CI_H_meta) %>% 
      dplyr::rename(
        estimate = beta_meta,
        p.value = p_meta,
        conf.low = CI_L_meta,
        conf.high = CI_H_meta
      ) %>% 
      mutate(Cohort = "Meta-analyzed")
  ) %>% 
  mutate(Cohort = factor(Cohort, levels = c("ADNI", "MCSA", "Meta-analyzed"))) %>% 
  mutate(Phenotype = factor(Phenotype, levels = c("Memory", "LMDR", "HippVol"))) %>% 
  ggplot(aes(x=estimate, y=Gene)) + 
  geom_errorbar(aes(xmax = conf.high, xmin = conf.low, color = Cohort), size = 0.5, width = 0.3, position=position_dodge(width = 0.5)) + 
  geom_point(aes(color = Cohort), size = 1.5, position=position_dodge(width = 0.5)) +
  scale_color_manual(values = c("ADNI" = "green4", "MCSA" = "blue4", "Meta-analyzed" = "cyan2")) +
  labs(title = expression("Blood " * italic("NDUFs") * "-endophenotype association")) +
  facet_grid(~Phenotype, scales = "free", space = "free") +
  theme_bw(base_size = 6.5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
        axis.text.y = element_text(face = "italic"))
dev.off()

pdf(paste0(outdir, "/1_12/NDUFB9_std.pdf"), width = 2.5, height = 2.2)
long_meta %>% 
  dplyr::select(Phenotype, estimate_MCSA, p.value_MCSA, conf.low_MCSA, conf.high_MCSA) %>% 
  rename_with(~sub("_MCSA", "", .x), everything()) %>% 
  mutate(Cohort = "MCSA") %>% 
  bind_rows(
    long_meta %>% 
      dplyr::select(Phenotype, estimate_ADNI, p.value_ADNI, conf.low_ADNI, conf.high_ADNI) %>% 
      rename_with(~sub("_ADNI", "", .x), everything()) %>% 
      mutate(Cohort = "ADNI")
  ) %>% 
  bind_rows(
    long_meta %>% 
      dplyr::select(Phenotype, beta_meta, p_meta, CI_L_meta, CI_H_meta) %>% 
      dplyr::rename(
        estimate = beta_meta,
        p.value = p_meta,
        conf.low = CI_L_meta,
        conf.high = CI_H_meta
      ) %>% 
      mutate(Cohort = "Meta-analyzed")
  ) %>% 
  mutate(Cohort = factor(Cohort, levels = c("ADNI", "MCSA", "Meta-analyzed"))) %>% 
  mutate(Phenotype = factor(Phenotype, levels = c("HippVol", "LMDR", "Memory"))) %>% 
  ggplot(aes(x=estimate, y=Phenotype)) + 
  geom_errorbar(aes(xmax = conf.high, xmin = conf.low, color = Cohort), size = 0.5, width = 0.3, position=position_dodge(width = 0.5)) + 
  geom_point(aes(color = Cohort), size = 2, position=position_dodge(width = 0.5)) +
  scale_color_manual(values = c("ADNI" = "green4", "MCSA" = "blue4", "Meta-analyzed" = "cyan2")) +
  labs(title = expression(atop("Blood " * italic("NDUFB9") * "-endophenotype", "longitudinal association"))) +
  theme_bw(base_size = 6.5) +
  plot_layout(guides = "collect")
dev.off()