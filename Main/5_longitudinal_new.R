# Author: Wei (Adelyn) Tsai; tsai.wei@mayo.edu
# If you are using this script, please cite our study

suppressPackageStartupMessages(library(tidyverse)) #v2.0.0
suppressPackageStartupMessages(library(broom)) #v1.0.12
suppressPackageStartupMessages(library(broom.mixed)) #v0.2.9
suppressPackageStartupMessages(library(lmerTest)) #v3.1
suppressPackageStartupMessages(library(pbapply)) #v1.7
suppressPackageStartupMessages(library(meta)) #v6.5
suppressPackageStartupMessages(library(openxlsx)) #v4.2.8
suppressPackageStartupMessages(library(patchwork)) #v1.3.2

source("./Codes/fns.R") #source a list of frequently used functions

# list of input directories where files are
indir <- "./indir"
residir <- "./outdir/resids"
nduf_assoc_dir <- "./outdir/nduf_assoc_de"

# output directory
outdir <- "./outdir/longitudinal"


# prepare files ----
comRes <- readRDS(paste0(residir, "/comRes.rds")) # from 1_CreateResids.R, 
nduf_bb <- read.csv(paste0(nduf_assoc_dir,"/nduf_bb.csv")) # generated in 4_nduf_assoc_de.R; also in supplementary table 12

# load in longitudinal files of mcsa and adni
# for both files, rows are donors, and there are at least two or more rows per donor because there are data from different time point.
# columns are cognitive endophenotypes/hippcampal volume at each time point. 
# there is also a column of difference in years of each time point from the baseline (YrDiff). 
# and there are columns of demographic information including sex, education yrs (Educ), age at baseline (Age_bl), intracranial volume (ICV) or MRI magnetic field (Mag)
mcsa_long <- read.csv(paste0(indir,"/mcsa_long.csv")) %>% mutate(SubjectID = as.character(SubjectID))
adni_long <- read.csv(paste0(indir, "/adni_long.csv"))

# check for each subject, how many missing across visits there are for different endophenotypes 
mcsa_na_check = mcsa_long %>% group_by(SubjectID) %>% summarise_all(list(~sum(!is.na(.))))
adni_na_check = adni_long %>% group_by(SubjectID) %>% summarise_all(list(~sum(!is.na(.))))

# mcsa input dataframe for longitudinal association
mcsa_df <- mcsa_long %>% 
  dplyr::select(SubjectID, YrDiff, Sex, Educ, Age_bl, APOE4, LMDR, Memory, HippVol, ICV) %>% 
  left_join(comRes[1:105, ] %>% as.data.frame() %>% rownames_to_column("SubjectID") %>% dplyr::select(SubjectID, any_of(nduf_bb$gene_id)))
# rename the columns from NDUF gene ID to gene name
mcsa_df <- mcsa_df %>% rename_at(vars(11:(ncol(.))), ~ nduf_bb$gene_name[match(names(mcsa_df)[11:(ncol(mcsa_df))], nduf_bb$gene_id)])

# adni input dataframe for longitudinal association
adni_df <- adni_long %>% 
  dplyr::select(SubjectID, YrDiff, Sex, Educ, Age_bl, APOE4, LMDR, Memory, HippVol, ICV, Mag) %>% 
  left_join(comRes[106:196, ] %>% as.data.frame() %>% rownames_to_column("SubjectID") %>% dplyr::select(SubjectID, any_of(nduf_bb$gene_id)))
# rename the columns from NDUF gene ID to gene name
adni_df <- adni_df %>% rename_at(vars(12:(ncol(.))), ~ nduf_bb$gene_name[match(names(adni_df)[12:(ncol(adni_df))], nduf_bb$gene_id)])


# run associations ----
# we save both the lmer models as well as the coefficients from the lmer models
# run association (MCSA)
mcsa_models <- list()
mcsa_coefs <- list()
for (i in c("Memory", "LMDR", "HippVol")){
  print(i)
  
  if (i != "HippVol"){
    inputdf <- mcsa_df %>% 
      filter(! SubjectID %in% mcsa_na_check$SubjectID[which(mcsa_na_check[[i]]<=1)]) # only retain subjects with more than 1 longitudinal measures
    for (j in c(nduf_bb$gene_name)){
      model <- lmer(paste0(i, "~", j , "*YrDiff + Age_bl + Sex + Educ + (1 + YrDiff|SubjectID)"), data=inputdf)
      
      df <- tidy(model, conf.int = T)
      
      mcsa_models[[paste0(i,"_",j)]] <- model
      mcsa_coefs[[paste0(i,"_",j)]] <- df
    }
    
  } else {
    inputdf <- mcsa_df %>% 
      filter(! SubjectID %in% mcsa_na_check$SubjectID[which(mcsa_na_check[[i]]<=1)])
    for (j in c(nduf_bb$gene_name)){
      model <- lmer(paste0(i, "~", j, "*YrDiff + Age_bl + Sex + Educ + ICV + (1 + YrDiff|SubjectID)"), data=inputdf)
      
      df <- tidy(model, conf.int = T)
      
      mcsa_models[[paste0(i,"_",j)]] <- model
      mcsa_coefs[[paste0(i,"_",j)]] <- df
    }
  }
}
# bind results
mcsa_coefs_dfs <- mcsa_coefs %>% 
  bind_rows(.id="Phenotype") %>% 
  separate(col = "Phenotype", into = c("Phenotype", "NDUF"), sep = "_") %>%
  filter(str_detect(term, ":")) %>%
  group_by(Phenotype) %>% 
  mutate(q.value = p.adjust(p.value, "fdr")) %>% 
  dplyr::select(-effect, -group)

# run association (ADNI)
adni_models <- list()
adni_coefs <- list()
for (i in c("Memory", "LMDR", "HippVol")){
  print(i)
  
  if (i != "HippVol"){
    inputdf <- adni_df %>% 
      filter(! SubjectID %in% adni_na_check$SubjectID[which(adni_na_check[[i]]<=1)]) # only retain subjects with more than 1 longitudinal measures
    for (j in c(nduf_bb$gene_name)){
      model <- lmer(paste0(i, "~", j , "*YrDiff + Age_bl + Sex + Educ + (1 + YrDiff|SubjectID)"), data=inputdf)
      
      df <- tidy(model, conf.int = T)
      
      adni_models[[paste0(i,"_",j)]] <- model
      adni_coefs[[paste0(i,"_",j)]] <- df
    }
    
  } else {
    inputdf <- adni_df %>% 
      filter(! SubjectID %in% adni_na_check$SubjectID[which(adni_na_check[[i]]<=1)])
    for (j in c(nduf_bb$gene_name)){
      model <- lmer(paste0(i, "~", j, "*YrDiff + Age_bl + Sex + Educ + ICV + Mag + (1 + YrDiff|SubjectID)"), data=inputdf)
      
      df <- tidy(model, conf.int = T)
      
      adni_models[[paste0(i,"_",j)]] <- model
      adni_coefs[[paste0(i,"_",j)]] <- df
    }
  }
}
# bind results
adni_coefs_dfs <- adni_coefs %>% 
  bind_rows(.id="Phenotype") %>% 
  separate(col = "Phenotype", into = c("Phenotype", "NDUF"), sep = "_") %>%
  filter(str_detect(term, ":")) %>%
  group_by(Phenotype) %>% 
  mutate(q.value = p.adjust(p.value, "fdr")) %>% 
  dplyr::select(-effect, -group)

# Meta-analysis ----
# prepare files
options(stringsAsFactors = FALSE)
pboptions(type = "timer")
meta_df = inner_join(mcsa_coefs_dfs, adni_coefs_dfs, by = c("Phenotype", "NDUF"), suffix = c("_MCSA", "_ADNI")) %>% 
  mutate(Pheno_NDUF = paste0(Phenotype, "_", NDUF)) %>% 
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
  dplyr::select(-Pheno_NDUF)


rm(fix_res, rand_res, het_res, meta_df, meta_raw, eff_idx1, eff_idx2, err_idx1, err_idx2)

# save results
write.xlsx(list("Meta" = meta_final, "MCSA" = mcsa_coefs_dfs, "ADNI" = adni_coefs_dfs), paste0(outdir,"/longitudinal_new.xlsx"))
# also save all the models
save(list =(ls(pattern = "_models")), file = paste0(outdir, "/longitudinal_models.RData"))