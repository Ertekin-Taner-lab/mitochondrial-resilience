# Author: Wei (Adelyn) Tsai; tsai.wei@mayo.edu
# If you are using this script, please cite our study
# This script addresses comment 3 about mediation from reviewer 1

suppressPackageStartupMessages(library(tidyverse)) #v2.0.0
suppressPackageStartupMessages(library(mediation)) #v4.5.0
suppressPackageStartupMessages(library(broom)) #v1.0.12
suppressPackageStartupMessages(library(openxlsx)) #v4.2.8

#######################################
# list of input/output directories
indir <- "./indir"
nduf_assoc_dir <- "./outdir/nduf_assoc_de"
residir <- "./outdir/resids"


outdir <- "./Rebuttal/outdir"
#######################################

#######################################
# prepare data frame for analysis
comTrait <- readRDS(paste0(indir,"/comTrait.rds")) # a file having cross-sectional endophenotypes for MCSA (from row 1-105) and ADNI (from row 106-196); row names are donorID
comRes <- readRDS(paste0(residir, "/comRes.rds")) # from 1_CreateResids.R, 
nduf_bb <- read.csv(paste0(nduf_assoc_dir, "/nduf_bb.csv")) #generated in 4_nduf_assoc_de.R


# prepare data frames combining phenotype and residuals of NDUFs
## MCSA
mcsa_med_df  <- comTrait[1:105,] %>% 
  rownames_to_column("SubjectID") %>% 
  dplyr::select(SubjectID, Sex, Educ, Age, LMDR, Memory, HippVol, Hipp_ICV) %>% 
  left_join(
    comRes[1:105, ] %>% as.data.frame() %>% rownames_to_column("SubjectID") %>% dplyr::select(SubjectID, any_of(nduf_bb$gene_id))
  )
### rename geneID to gene name of NDUFs
mcsa_med_df  <- mcsa_med_df  %>% rename_at(vars(9:(ncol(.))), ~ nduf_bb$gene_name[match(names(mcsa_med_df)[9:(ncol(mcsa_med_df))], nduf_bb$gene_id)])


## ADNI
adni_med_df  <- comTrait[106:196,] %>% 
  rownames_to_column("SubjectID") %>% 
  dplyr::select(SubjectID, Sex, Educ, Age, LMDR, Memory, HippVol, Hipp_ICV, Mag_HippICV) %>% 
  left_join(
    comRes[106:196, ] %>% as.data.frame() %>% rownames_to_column("SubjectID") %>% dplyr::select(SubjectID, any_of(nduf_bb$gene_id))
  )
### rename geneID to gene name of NDUFs
adni_med_df  <- adni_med_df  %>% rename_at(vars(10:(ncol(.))), ~ nduf_bb$gene_name[match(names(adni_med_df)[10:(ncol(adni_med_df))], nduf_bb$gene_id)])


#######################################
# perform mediation to see whether cell type fraction is a mediator between phenotype and ME
# function for standardization
standardize_FUN <- function(x) {
  
  x_med <- median(x, na.rm = TRUE)
  x_mad <- mad(x, na.rm = TRUE)
  x <- ((x - x_med) / x_mad)
  x
  
}

mcsa_med_df <- mcsa_med_df %>% 
  mutate_at(.vars = c("Memory", "LMDR", "HippVol"), .funs = standardize_FUN) %>% 
  mutate_at(vars(NDUFAB1:NDUFA7), .funs = standardize_FUN)
adni_med_df <- adni_med_df %>% 
  mutate_at(.vars = c("Memory", "LMDR", "HippVol"), .funs = standardize_FUN) %>% 
  mutate_at(vars(NDUFAB1:NDUFA7), .funs = standardize_FUN)


# mediation function
run_mediation_all <- function(data, 
                              treatments, 
                              mediators, 
                              outcomes, 
                              covariates,
                              Cohort) {
  # Helper function to fit models and extract summaries
  run_mediation_one <- function(treat, med, out) {
    print(paste0(treat, "; ", med, ";", out))
    data <- data %>% 
      dplyr::select(any_of(c(treat, med, out, covariates))) %>% 
      dplyr::rename(X = all_of(treat), M= all_of(med), Y = all_of(out)) %>% 
      drop_na()
    if (Cohort == "MCSA"){
        # Fit the model
        med.fit <- lm(M ~ X + Age + Sex + Educ + Hipp_ICV, data = data)
        out.fit <- lm(Y ~ M + X + Age + Sex + Educ + Hipp_ICV, data = data)
        int.fit <- lm(Y ~ M * X + Age + Sex + Educ + Hipp_ICV, data = data)
        
        # Mediation analysis
        med.out <- mediate(med.fit, out.fit, treat = "X", mediator = "M", boot = TRUE, boot.ci.type = "bca", sims = 1000)
        
        # Sensitivity analysis
        sens.out <- medsens(med.out, rho.by = 0.1, effect.type = "both", sims = 1000)
      }
    else {
        # Fit the model
        med.fit <- lm(M ~ X + Age + Sex + Educ + Hipp_ICV + Mag_HippICV, data = data)
        out.fit <- lm(Y ~ M + X + Age + Sex + Educ + Hipp_ICV + Mag_HippICV, data = data)
        int.fit <- lm(Y ~ M * X + Age + Sex + Educ + Hipp_ICV + Mag_HippICV, data = data)
        
        # Mediation analysis
        med.out <- mediate(med.fit, out.fit, treat = "X", mediator = "M", boot = TRUE, boot.ci.type = "bca", sims = 1000)
        
        # Sensitivity analysis
        sens.out <- medsens(med.out, rho.by = 0.1, effect.type = "both", sims = 1000)
    }
    # Extract mediation stats
    df <- tibble(
      ACME_est = med.out$d0,
      ACME_p = med.out$d0.p,
      ACME_L95 = med.out$d0.ci[1],
      ACME_H95 = med.out$d0.ci[2],
      ADE_est = med.out$z0,
      ADE_p = med.out$z0.p,
      ADE_L95 = med.out$z0.ci[1],
      ADE_H95 = med.out$z0.ci[2],
      TotalEffect_est = med.out$tau.coef,
      TotalEffect_p = med.out$tau.p,
      TotalEffect_L95 = med.out$tau.ci[1],
      TotalEffect_H95 = med.out$tau.ci[2],
      PropMediated = med.out$n0,
      PropMediated_p = med.out$n0.p,
      PropMediated_L95 = med.out$n0.ci[1],
      PropMediated_H95 = med.out$n0.ci[2],
      Sensitivity_rho_ACME0 = sens.out$err.cr.d,
      Sensitivity_rho_ADE0 = sens.out$err.cr.z,
      Int_est = last(tidy(int.fit)$estimate),
      Int_p = last(tidy(int.fit)$p.value),
      Int_L95 = last(tidy(int.fit, conf.int = T)$conf.low),
      Int_H95 = last(tidy(int.fit, conf.int = T)$conf.high)
    )
    
    return(df)
  }
  
  # Run for all combinations
  results <- expand.grid(
    Treatment = treatments,
    Mediator = mediators,
    Outcome = outcomes,
    stringsAsFactors = F
  ) %>%
    mutate(result = pmap(list(Treatment, Mediator, Outcome),
                         ~ run_mediation_one(..1, ..2, ..3))) %>%
    unnest(result)
  
  return(results)
}

set.seed(123)
mcsa_medRes <- run_mediation_all(
      data = mcsa_med_df,
      treatments = nduf_bb$gene_name,
      mediators = "HippVol",
      outcomes = c("Memory", "LMDR"),
      covariates = c("Age", "Sex", "Educ", "Hipp_ICV"),
      Cohort = "MCSA"
    )
adni_medRes <- run_mediation_all(
  data = adni_med_df,
  treatments = nduf_bb$gene_name,
  mediators = "HippVol",
  outcomes = c("Memory", "LMDR"),
  covariates = c("Age", "Sex", "Educ", "Hipp_ICV", "Mag_HippICV"),
  Cohort = "ADNI"
)
# sensitivity analysis to see the results of swapping
mcsa_medRes_sens <- run_mediation_all(
  data = mcsa_med_df,
  treatments = nduf_bb$gene_name,
  mediators = c("Memory", "LMDR"),
  outcomes = "HippVol",
  covariates = c("Age", "Sex", "Educ", "Hipp_ICV"),
  Cohort = "MCSA"
)
adni_medRes_sens <- run_mediation_all(
  data = adni_med_df,
  treatments = nduf_bb$gene_name,
  mediators = c("Memory", "LMDR"),
  outcomes =  "HippVol",
  covariates = c("Age", "Sex", "Educ", "Hipp_ICV", "Mag_HippICV"),
  Cohort = "ADNI"
)


write.xlsx(list(
  "MCSA" = mcsa_medRes,
  "ADNI" = adni_medRes,
  "MCSA_sens" = mcsa_medRes_sens,
  "ADNI_sens" = adni_medRes_sens),
  paste0(outdir, "/1_3/medRes_ndufs_HV_cog.xlsx")
) # the organized table from these results are in supplementary table 15

#######################################
