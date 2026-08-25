# Author: Wei (Adelyn) Tsai; tsai.wei@mayo.edu
# If you are using this script, please cite our study
# This script addresses comment 13 about outlier diagnostics from reviewer 1

suppressPackageStartupMessages(library(tidyverse)) #v2.0.0
suppressPackageStartupMessages(library(patchwork)) #v1.3.2
suppressPackageStartupMessages(library(broom)) #v1.0.12
suppressPackageStartupMessages(library(broom.mixed)) #v0.2.9
suppressPackageStartupMessages(library(openxlsx)) #v4.2.8
suppressPackageStartupMessages(library(meta)) #v6.5
suppressPackageStartupMessages(library(pbapply)) #v1.7
suppressPackageStartupMessages(library(lmerTest)) #v3.1
suppressPackageStartupMessages(library(ggpubr)) #v0.6.3


source("./Codes/fns.R") #source a list of frequently used functions


# list of input/output directories
indir <- "./indir"
blwgcna_dir <- "./outdir/bl_wgcna"
brwgcna_dir <- "./outdir/br_wgcna"
residir <- "./outdir/resids"
nduf_assoc_dir <- "./outdir/nduf_assoc_de"

outdir <- "./Rebuttal/outdir"


# 1) Blood M3 with memory, LMDR, hippocample volume ----
# prepare data frame for analysis
comTrait <- readRDS(paste0(indir,"/comTrait.rds")) # a file having cross-sectional endophenotypes for MCSA (from row 1-105) and ADNI (from row 106-196); row names are donorID
MEs1 <- read.csv(paste0(blwgcna_dir, "/MEs1_SP12.csv")) # from 2_bl_wgcna.R
analysis_df <- comTrait %>% 
  rownames_to_column("SubjectID") %>% 
  left_join(MEs1 %>% dplyr::select(SubjectID, ME3)) %>% 
  mutate(Mag_HippICV = factor(Mag_HippICV))
# we performed rank-transformation for the dependent and independent variables
analysis_df[1:105,] <- analysis_df[1:105, ] %>% 
  mutate_at(.vars = c("Memory", "LMDR", "HippVol", "ME3"), .funs = rank, na.last = "keep")
analysis_df[106:196,] <- analysis_df[106:196, ] %>% 
  mutate_at(.vars = c("Memory", "LMDR", "HippVol", "ME3"), .funs = rank, na.last = "keep")

# run analysis to get new estimates after rank transformation
mcsa_m3 <- list()
for (i in c("Memory", "LMDR", "HippVol")){
  if (i != "HippVol"){
    mcsa_m3[[i]] <- tidy(lm(paste0("ME3 ~", i, " + Age + Sex + Educ"), data=analysis_df[1:105, ]), conf.int = TRUE) %>% filter(term == i)
  }
  else{
    mcsa_m3[[i]] <- tidy(lm(paste0("ME3 ~", i, " + Age + Sex + Educ + Hipp_ICV"), data=analysis_df[1:105, ]), conf.int = TRUE) %>% filter(term == i)
  }
}
adni_m3 <- list()
for (i in c("Memory", "LMDR", "HippVol")){
  if (i != "HippVol"){
    adni_m3[[i]] <- tidy(lm(paste0("ME3 ~", i, " + Age + Sex + Educ"), data=analysis_df[106:196, ]), conf.int = TRUE) %>% filter(term == i)
  }
  else{
    adni_m3[[i]] <- tidy(lm(paste0("ME3 ~", i, " + Age + Sex + Educ + Hipp_ICV + Mag_HippICV"), data=analysis_df[106:196, ]), conf.int = TRUE) %>% filter(term == i)
  }
}
#bind results
mcsa_m3_dfs <- mcsa_m3 %>% bind_rows(.id="Phenotype")
adni_m3_dfs <- adni_m3 %>% bind_rows(.id="Phenotype")

# do meta-analysis using new est and std.error
meta_df = inner_join(mcsa_m3_dfs, adni_m3_dfs, by = c("Phenotype", "term"), suffix = c("_MCSA", "_ADNI")) %>% 
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
m3_meta = meta_df %>%
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
  dplyr::select(-term)


rm(fix_res, rand_res, het_res, meta_df, meta_raw, eff_idx1, eff_idx2, err_idx1, err_idx2)



# 2) Brain modules with phenotype in AMP-AD cohort ----
# prepare dataframe for association
consMEs_br <- read.csv(paste0(brwgcna_dir, "/MEs1_SP12_Br.csv")) #created from 3_br_wgcna.R
comTrait_br <- readRDS(paste0(indir, "/comTrait_br.rds")) # a file containing organized phenotypes & covariates of brain datasets; rows are donors columns are phenotypes/covariates
consMEs_br_long_sub <- consMEs_br %>%
  dplyr::select(SubjectID, ME1, ME17, ME26) %>%
  pivot_longer(-SubjectID, names_to = "Module", values_to = "ME") %>%
  inner_join(comTrait_br, by=c("SubjectID"="SampleID")) %>%
  named_group_split(Region) %>%
  map(., ~mutate(.x, 
                 DX=if_else(DX=="Control", 1L, 0L), 
                 Batch = factor(Batch, levels=unique(Batch)),
                 Source = factor(Source, levels=unique(Source)))) %>%
  map(., ~named_group_split(.x, Module)) %>% # this is split so in each list is a dataframe for each region and module
  map(~map(.x, ~mutate_at(.x, .vars = c("ME"), .funs = rank, na.last = "keep"))) # rank transform module eigengenes


# run association 
br_mods_ls <- list()
for (region in c("Mayo\nSTG", "CER", "FP", "MSSM\nSTG", "PHG", "IFG", "DLPFC")){
  print(region)
  if (region %in% c("Mayo\nSTG", "CER", "DLPFC")){
    for (pheno in c("DX", "Braak", "Thal")){ 
      print(pheno)
      formula <- paste0("ME ~ ", pheno, "+ AOD + Sex_M1F0")
      df <- consMEs_br_long_sub[[region]] %>% 
        map(~ lm(formula, data = .x)) %>% 
        map_dfr(broom::tidy, conf.int = TRUE, .id = "Module") %>%
        filter(term == pheno)
      br_mods_ls[[region]][[pheno]] <- df
    }
  } else{
    for (pheno in c("DX", "Braak", "CERAD")){
      print(pheno)
      formula <- paste0("ME ~ ", pheno, "+ AOD + Sex_M1F0")
      df <- consMEs_br_long_sub[[region]] %>% 
        map(~ lm(formula, data = .x)) %>% 
        map_dfr(broom::tidy, conf.int = TRUE, .id = "Module") %>%
        filter(term == pheno)
      br_mods_ls[[region]][[pheno]] <- df
    }
  } 
}
# bind results
br_mods_df <- br_mods_ls %>% 
  map(~ bind_rows(.x)) %>% 
  bind_rows(.id = "Study") %>% 
  mutate(Study = case_when(
    Study == "CER" ~ "Mayo\nCER",
    Study %in% c("FP", "PHG", "IFG") ~ paste0("MSSM\n", Study),
    Study == "DLPFC" ~ "ROSMAP\nDLPFC",
    TRUE ~ Study
  )) %>% 
group_by(Study, term) %>% 
mutate(q.value = p.adjust(p.value, method = "fdr"))


# 3) Blood NDUFs with memory, LMDR, HippVol ----
# prepare files
# MCSA
nduf_bb <- read.csv(paste0(nduf_assoc_dir, "/nduf_bb.csv")) #generated in 4_nduf_assoc_de.R
cqn_raw <- read_delim(paste0(indir, "/PaX108_R01resilience_gene_CQN_neg3_postQC_s105.txt"), delim = "\t") # MCSA cqn file for 18046 genes. Each row is a gene and each column is a donor
cqn_use <- cqn_raw %>% filter(GeneId %in% nduf_bb$gene_id)
mcsa_covars <- read_delim(paste0(indir,"/R01_Resilience_PAXgene_Covars.txt"), delim = "\t") %>% mutate(SubjectID=as.character(ptnum)) #covariates file for MCSA, rows are donor columns are covariates
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
  # rank dependent and independent variables
  map(~mutate_at(.x, .vars = c("Expr", "Memory", "LMDR", "HippVol"), .funs = rank, na.last = "keep"))
# run analysis
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
  # rank dependent and independent variables
  map(~mutate_at(.x, .vars = c("Expr", "Memory", "LMDR", "HippVol"), .funs = rank, na.last = "keep"))


#run analysis
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



# 4) NDUFB9 longitudinal association using interaction model ----
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
              mutate(ENSG00000147684 = rank(ENSG00000147684, na.last = "keep"))) %>% 
  dplyr::rename(NDUFB9  = ENSG00000147684) 

# run association (MCSA)
mcsa_long_assoc <- list()
for (i in c("Memory", "LMDR", "HippVol")){
  print(i)
  
  if (i != "HippVol"){
    inputdf <- mcsa_long_df %>% 
      filter(! SubjectID %in% mcsa_na_check$SubjectID[which(mcsa_na_check[[i]]<=1)]) %>%  # only retain subjects with more than 1 longitudinal measures
      mutate_at(.vars = i, .funs = rank, na.last = "keep")
    mcsa_long_assoc[[i]] <- tidy(lmer(paste0(i, "~", "NDUFB9*YrDiff + Age_bl + Sex + Educ + (1 + YrDiff|SubjectID)"), data=inputdf), conf.int = T) %>% 
      filter(term == "NDUFB9:YrDiff")
  }
  
  else {
    inputdf <- mcsa_long_df %>% 
      filter(! SubjectID %in% mcsa_na_check$SubjectID[which(mcsa_na_check[[i]]<=1)]) %>% 
      mutate_at(.vars = i, .funs = rank, na.last = "keep")
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
              mutate(ENSG00000147684 = rank(ENSG00000147684, na.last = "keep"))) %>% 
  dplyr::rename(NDUFB9 = ENSG00000147684)
# run association (ADNI)
adni_long_assoc <- list()
for (i in c("Memory", "LMDR", "HippVol")){
  print(i)
  
  if (i != "HippVol"){
    inputdf <- adni_long_df %>% 
      filter(! SubjectID %in% adni_na_check$SubjectID[which(adni_na_check[[i]]<=1)]) %>%  # only retain subjects with more than 1 longitudinal measures
      mutate_at(.vars = i, .funs = rank, na.last = "keep")
    adni_long_assoc[[i]] <- tidy(lmer(paste0(i, "~", "NDUFB9*YrDiff + Age_bl + Sex + Educ + (1 + YrDiff|SubjectID)"), data=inputdf), conf.int = T) %>% 
      filter(term == "NDUFB9:YrDiff")
  }
  else {
    inputdf <- adni_long_df %>% 
      filter(! SubjectID %in% adni_na_check$SubjectID[which(adni_na_check[[i]]<=1)]) %>% 
      mutate_at(.vars = i, .funs = rank, na.last = "keep")
    adni_long_assoc[[i]] <-  tidy(lmer(paste0(i, "~", "NDUFB9*YrDiff + Age_bl + Sex + Educ + ICV + Mag + (1 + YrDiff|SubjectID)"), data=inputdf), conf.int = T) %>% 
      filter(term == "NDUFB9:YrDiff")
  }
}
# bind results
adni_long_res <- adni_long_assoc %>% bind_rows(.id="Phenotype") %>% dplyr::select(-group, -effect)

# conduct meta-analysis
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

write.xlsx(list("Bl_M3" = m3_meta, "Br_Mods" = br_mods_df, "Bl_NDUFs" = bl_ndufs_meta, "Longitudinal" = long_meta), paste0(outdir, "/1_13/ranked_assoc_df.xlsx")) # these are organized in supplementary table 7a-d


# outlier diagnosis for longitudinal model ----
suppressPackageStartupMessages(library(HLMdiag)) #v0.5.1

# load the models saved in 5_longitudinal.R
long_dir <- "./outdir/longitudinal"
load(paste0(long_dir, "/longitudinal_models.RData"), verbose = T)

# 1. Cook's distance (fixed effects)
infl_l1_mcsa <- hlm_influence(mcsa_models$HippVol_NDUFB9, level = 1) # first level influence; this shows whether a specific observation may be an outlier
summary(infl_l1_mcsa) # Cutoff: 4/n = 4/365 (365 observations) 
infl_l2_mcsa <- hlm_influence(mcsa_models$HippVol_NDUFB9, level = "SubjectID")
summary(infl_l2_mcsa) # Cutoff: 4/n = 4/86 (86 subjects)

infl_l1_adni <- hlm_influence(adni_models$HippVol_NDUFB9, level = 1) # first level influence
summary(infl_l1_adni) # Cutoff: 4/n = 4/596 (596 observations) 
infl_l2_adni <- hlm_influence(adni_models$HippVol_NDUFB9, level = "SubjectID")
summary(infl_l2_adni) # cutoff: 4/n = 4/90 (90 subjects)


# 2. Boxplot to visualize residuals
l2_resids_mcsa <- hlm_resid(mcsa_models$HippVol_NDUFB9, level = "SubjectID", include.ls = F)
l2_resids_adni <- hlm_resid(adni_models$HippVol_NDUFB9, level = "SubjectID", include.ls = F)

pdf(paste0(outdir, "/1_13/LongResidsP_lvl2.pdf"), height = 3, width = 3)
boxplot(l2_resids_mcsa$.ranef.intercept, main = "MCSA: random intercept distribution", ylab = "Random intercept", cex.axis = 0.5, cex.lab = 0.6, cex.main = 0.7)
boxplot(l2_resids_mcsa$.ranef.yr_diff, main = "MCSA: random slope distribution", ylab = "Random slope", cex.axis = 0.5, cex.lab = 0.6, cex.main = 0.7)
boxplot(l2_resids_adni$.ranef.intercept, main = "ADNI: random intercept distribution", ylab = "Random intercept", cex.axis = 0.5, cex.lab = 0.6, cex.main = 0.7)
boxplot(l2_resids_adni$.ranef.yr_diff, main = "ADNI: random slope distribution", ylab = "Random slope", cex.axis = 0.5, cex.lab = 0.6, cex.main = 0.7)
dev.off() # this is in Extended data fig.7e

## repeated analyses without outliers ----
### we found 1 outlier from random slope distribution in MCSA and 3 in ADNI - remove those and repeated the analyses
ol_long_mcsa <- l2_resids_mcsa %>% 
  mutate(OL = ifelse(findoutlier(.ranef.yr_diff), SubjectID, NA)) %>% 
  filter(!is.na(OL)) %>% 
  pull(OL)
ol_long_adni <- l2_resids_adni %>% 
  mutate(OL = ifelse(findoutlier(.ranef.yr_diff), SubjectID, NA)) %>% 
  filter(!is.na(OL)) %>% 
  pull(OL)

inputdf <- mcsa_long %>% 
  dplyr::select(SubjectID, YrDiff, Sex, Educ, Age_bl, APOE4, LMDR, Memory, HippVol, ICV) %>% 
  left_join(comRes[1:105, ] %>% 
              as.data.frame() %>% 
              rownames_to_column("SubjectID") %>% 
              dplyr::select(SubjectID, ENSG00000147684)
            ) %>% 
  dplyr::rename(NDUFB9  = ENSG00000147684) %>% 
  filter(! SubjectID %in% mcsa_na_check$SubjectID[which(mcsa_na_check[["HippVol"]]<=1)]) %>% 
  filter(! SubjectID %in% ol_long_mcsa)
mcsa_long_hv_noOL <- tidy(lmer(HippVol ~ NDUFB9*YrDiff + Age_bl + Sex + Educ + ICV + (1 + YrDiff|SubjectID), data=inputdf), conf.int = T) %>%
  filter(term == "NDUFB9:YrDiff") %>% 
  dplyr::select(-effect, -group)

rm(inputdf)

inputdf <- adni_long %>% 
  dplyr::select(SubjectID, YrDiff, Sex, Educ, Age_bl, APOE4, LMDR, Memory, HippVol, ICV, Mag) %>% 
  left_join(comRes[106:196, ] %>% 
              as.data.frame() %>% 
              rownames_to_column("SubjectID") %>% 
              dplyr::select(SubjectID, ENSG00000147684)
            ) %>% 
  dplyr::rename(NDUFB9 = ENSG00000147684) %>% 
  filter(! SubjectID %in% adni_na_check$SubjectID[which(adni_na_check[["HippVol"]]<=1)]) %>% 
  filter(! SubjectID %in% ol_long_adni)
adni_long_hv_noOL <-  tidy(lmer(HippVol ~ NDUFB9*YrDiff + Age_bl + Sex + Educ + ICV + Mag + (1 + YrDiff|SubjectID), data=inputdf), conf.int = T) %>% 
  filter(term == "NDUFB9:YrDiff") %>% 
  dplyr::select(-effect, -group)
rm(inputdf)

# conduct meta-analysis
meta_df = inner_join(mcsa_long_hv_noOL, adni_long_hv_noOL, by = c("term"), suffix = c("_MCSA", "_ADNI"))
rownames(meta_df) = "HippVol"
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
long_meta_noOL = meta_df %>%
  mutate(q_random = p.adjust(p_random, method = "fdr"),
         model_meta = ifelse(I2 < 0.25, "Fixed", "Random"),
         beta_meta = ifelse(model_meta =="Fixed", beta_fixed, beta_random),
         se_meta = ifelse(model_meta =="Fixed", se_fixed, se_random),
         p_meta = ifelse(model_meta =="Fixed", p_fixed, p_random),
         CI_H_meta = ifelse(model_meta =="Fixed", upper_fixed, upper_random),
         CI_L_meta = ifelse(model_meta == "Fixed", lower_fixed, lower_random)) %>%
  dplyr::select(-c(ends_with("_random"), ends_with("_fixed"))) 
rm(meta_df, eff_idx1, eff_idx2, err_idx1, err_idx2, meta_raw, fix_res, rand_res, het_res)
write.xlsx(list("HV_Long_noOL" = long_meta_noOL), paste0(outdir, "/1_13/HV_noOL_results.xlsx")) # this is in supplementary table 7e