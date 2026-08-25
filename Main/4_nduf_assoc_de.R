# Author: Wei (Adelyn) Tsai; tsai.wei@mayo.edu
# If you are using this script, please cite our study

suppressPackageStartupMessages(library(broom)) #v1.0.12
suppressPackageStartupMessages(library(broom.mixed)) #v0.2.9
suppressPackageStartupMessages(library(tidyverse)) #v2.0.0
suppressPackageStartupMessages(library(lmerTest)) #v3.1
suppressPackageStartupMessages(library(openxlsx)) #v4.2.8
suppressPackageStartupMessages(library(meta)) # v6.5
suppressPackageStartupMessages(library(pbapply)) # v1.7

indir <- ""
outdir <- ""

source("./Codes/fns.R") #source a list of frequently used functions

MM_tbl <- read.csv(paste0(outdir,"/bl_wgcna/MM_tbl_2023-08-15.csv")) # from 2_bl_wgcna.R
MM_tbl_br <- read.csv(paste0(outdir,"/br_wgcna/MM_tbl_br.csv")) # from 3_br_wgcna.R
# create a table focusing on NDUFs
nduf_bb <- MM_tbl %>% 
  filter(Module ==3) %>% 
  dplyr::select(gene_id, gene_name, Module, MM) %>% 
  filter(grepl("NDUF", .$gene_name)) %>% 
  left_join(MM_tbl_br %>% 
              dplyr::select(gene_id, Module, MM), by=c("gene_id"), suffix=c("_blood", "_consbrain"))
write.csv(nduf_bb, paste0(outdir,"/nduf_assoc_de/nduf_bb.csv")) # this is supplementary table 12

comTrait <- readRDS(paste0(indir,"/comTrait.rds")) # a file having cross-sectional endophenotypes for MCSA (from row 1-105) and ADNI (from row 106-196); row names are donorID

# Blood ----
## mcsa ----
cqn_raw <- read_delim(paste0(indir, "/PaX108_R01resilience_gene_CQN_neg3_postQC_s105.txt"), delim = "\t") # MCSA cqn file for 18046 genes. Each row is a gene and each column is a donor
cqn_use <- cqn_raw %>% filter(GeneId %in% nduf_bb$gene_id)

mcsa_covars <- read_delim(paste0(indir,"/R01_Resilience_PAXgene_Covars.txt"), delim = "\t") %>% mutate(SubjectID=as.character(ptnum)) #covariates file for MCSA, rows are donor columns are covariates
## transform to long format
cqn_long <- cqn_use %>%
  dplyr::select(-Chromosome, -Start, -End, -Length, -GeneName, -GeneBiotype) %>%
  tidyr::pivot_longer(-GeneId, names_to = "SubjectID", values_to = "Expr") %>%
  inner_join(dplyr::select(mcsa_covars, SubjectID, RIN, PAXgene_flowcell, PAXgene_Batch), by = "SubjectID") %>%
  mutate(PAXgene_Batch = factor(PAXgene_Batch),
         PAXgene_flowcell = as.factor(PAXgene_flowcell))  %>%
  inner_join(comTrait %>% 
               tibble::rownames_to_column("SubjectID") %>% 
               dplyr::select(SubjectID, Sex, Educ, Age, LMDR, Memory, HippVol, Hipp_ICV), by="SubjectID") %>% 
  named_group_split(GeneId) 
## run regression
mcsa_fit_coefs <- list()
for (i in c("LMDR", "Memory", "HippVol")){
  if (i != "HippVol"){
    df <- cqn_long %>%
      map(~ lm(paste0("Expr ~ ", i, "+ Age + Sex + Educ + RIN + PAXgene_flowcell + PAXgene_Batch"), data = .x)) %>%
      map_dfr(broom::tidy, conf.int = TRUE, .id = "GeneID") %>%
      filter(term == i) %>%
      mutate(q.value = p.adjust(p.value, method = "fdr"))
    
    mcsa_fit_coefs[[i]] <- df
  }else{
    df <- cqn_long %>%
      map(~ lm(paste0("Expr ~ ", i, "+ Age + Sex + Educ + Hipp_ICV + RIN + PAXgene_flowcell + PAXgene_Batch"), data = .x)) %>%
      map_dfr(broom::tidy, conf.int = TRUE, .id = "GeneID") %>%
      filter(term == i) %>%
      mutate(q.value = p.adjust(p.value, method = "fdr"))
    
    mcsa_fit_coefs[[i]] <- df
  }
  
}
mcsa_fit_coefs_df <- bind_rows(mcsa_fit_coefs, .id="Phenotype") %>% dplyr::select(-term)

rm(i, df)

## adni ----
adni_mc <- readRDS(paste0(indir,"/ADNI_s91_microarray.rds")) #microarray file for 91 donors from ADNI. This has 10116 protein-coding genes, gene id as column names, donorID as row names
adni_covars <- read_delim(paste0(indir,"/ADNI_s91_covars.txt"), delim = "\t") #covariates file for 91 donors from ADNI. Rows are donors columns are covariates
## transform to long format
mc_long <- adni_mc %>% 
  dplyr::select(nduf_bb$gene_id) %>%
  tibble::rownames_to_column("SubjectID") %>% 
  filter(SubjectID %in% adni_covars$SubjectID) %>% 
  pivot_longer(2:21, names_to = "Gene", values_to = "Expr") %>% 
  inner_join(adni_covars %>% dplyr::select(SubjectID, RIN, AffyPlate.meanREE, SITE.meanREE), by="SubjectID") %>% 
  inner_join(comTrait %>% 
               tibble::rownames_to_column("SubjectID") %>% 
               dplyr::select(SubjectID, Sex, Educ, Age, LMDR, Memory, HippVol, Hipp_ICV, Mag_HippICV) %>%
               mutate(Mag_HippICV=factor(Mag_HippICV)), by="SubjectID") %>% 
  named_group_split(Gene) 
## run regression
adni_fit_coefs <- list()
for (i in c("LMDR", "Memory", "HippVol")){
  if(i != "HippVol"){
    df <- mc_long %>%
      map(~ lm(paste0("Expr ~ ", i, "+ Age + Sex + Educ + RIN + AffyPlate.meanREE + SITE.meanREE"), data = .x)) %>%
      map_dfr(broom::tidy, conf.int = TRUE, .id = "GeneID") %>%
      filter(term == i) %>%
      mutate(q.value = p.adjust(p.value, method = "fdr"))
    
    adni_fit_coefs[[i]] <- df
  }else{
    df <- mc_long %>%
      map(~ lm(Expr ~ HippVol + Age + Sex + Educ + Hipp_ICV + Mag_HippICV + RIN + AffyPlate.meanREE + SITE.meanREE, data = .x)) %>%
      map_dfr(broom::tidy, conf.int = TRUE, .id = "GeneID") %>%
      filter(term == "HippVol") %>%
      mutate(q.value = p.adjust(p.value, method = "fdr"))
    
    adni_fit_coefs[[i]] <- df
  }
}
adni_fit_coefs_df <- bind_rows(adni_fit_coefs, .id="Phenotype") %>% dplyr::select(-term)

rm(i, df)

## meta-analyze(MCSA&ADNI)----
options(stringsAsFactors = FALSE)
pboptions(type = "timer")
## prepare meta-analysis input df
meta_df = inner_join(mcsa_fit_coefs_df, adni_fit_coefs_df, by = c("Phenotype", "GeneID"), suffix = c("_MCSA", "_ADNI")) %>% 
  mutate(Pheno_NDUF = paste0(Phenotype, "_", GeneID)) %>% 
  column_to_rownames("Pheno_NDUF")

eff_idx1 = which(colnames(meta_df) == "estimate_MCSA")
eff_idx2 = which(colnames(meta_df) == "estimate_ADNI")
err_idx1 = which(colnames(meta_df) == "std.error_MCSA")
err_idx2 = which(colnames(meta_df) == "std.error_ADNI")

meta_raw = meta_df %>% pbapply(1, meta_gen_fn, eff_idx1 = eff_idx1, eff_idx2 = eff_idx2, err_idx1 = err_idx1, err_idx2 = err_idx2)
## Get all fixed-effect results
fix_res = meta_raw %>%
  pbsapply(function(x) x$fixed) %>%
  t() %>%
  as.data.frame() %>%
  rename_with(~paste0(.x, "_fixed")) %>%
  rownames_to_column(var = "Pheno_NDUF") %>%
  mutate_all(.funs = ~unlist(.x)) %>% 
  mutate_all(.funs = ~unname(.x))

## Get all random-effect results
rand_res = meta_raw %>%
  pbsapply(function(x) x$random) %>%
  t() %>%
  data.frame() %>%
  rename_with(~paste0(.x, "_random")) %>%
  rownames_to_column(var = "Pheno_NDUF") %>%
  mutate_all(.funs = ~unlist(.x)) %>%
  mutate_all(.funs = ~unname(.x))

## Get heterogenity meassurments Q, I2, tau2
het_res = meta_raw %>%
  map(~.x[c("Q","df.Q","pval.Q","I2", "lower.I2", "upper.I2", "tau2", "se.tau2", "H", "lower.H", "upper.H")]) %>%
  bind_rows(.id = "Pheno_NDUF")

## Join the results to the original regression coeffcients
meta_df = meta_df %>%
  rownames_to_column(var = "Pheno_NDUF") %>% 
  inner_join(fix_res) %>%
  inner_join(rand_res) %>% 
  rename_with(~sub("^TE", "beta", .x), everything()) %>%
  rename_with(~sub("^seTE", "se", .x), everything()) %>%
  inner_join(het_res, by = "Pheno_NDUF")


## meta-analyze
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
  inner_join(MM_tbl %>% dplyr::select(gene_id, gene_name, Module, MM), by=c("GeneID"="gene_id")) %>%
  dplyr::select(Phenotype, GeneID, gene_name, Module, MM, everything()) %>%
  dplyr::select(-Pheno_NDUF)


rm(fix_res, rand_res, het_res, meta_df, meta_raw, eff_idx1, eff_idx2, err_idx1, err_idx2)

write.csv(meta_final, paste0(outdir,"/nduf_assoc_de/Bl_Assoc.csv"),row.names = F)


# Brain ----
cohort= c("Mayo", "MSSM", "ROSMAP")
region= c("Mayo\nSTG", "CER", "FP", "MSSM\nSTG", "PHG", "IFG", "DLPFC")

comTrait_br <- readRDS(paste0(indir, "/comTrait_br.rds")) # a file containing organized phenotypes & covariates of brain datasets; rows are donors columns are phenotypes/covariates
br_cqn2use <- readRDS(paste0(indir, "/br_cqn2use.rds")) # a list containing cqn of brain datasets. In each list is a dataframe where rows are genes and columns are donors
br_cqn_long_ls <- br_cqn2use %>% map(., ~(pivot_longer(.x, cols = c(4:ncol(.x)), names_to="SampleID", values_to="CQN")))

for (reg in names(br_cqn_long_ls)){
  br_cqn_long_ls[[reg]] <- br_cqn_long_ls[[reg]] %>%
    filter(gene_id %in% nduf_bb$gene_id) %>% 
    inner_join(comTrait_br, by=c("Cohort", "Region", "SampleID")) %>%
    mutate(DX=if_else(DX=="Control", 1L, 0L), 
           Batch = factor(Batch, levels=unique(Batch)),
           Source = factor(Source, levels = unique(Source))) %>%
    named_group_split(gene_id) # this creates a list of 20 dataframes for 20 NDUFs per dataset
}
## regressio function
assoc_br <- function(region){
  if (region %in% c("Mayo\nSTG", "CER")){
    fit_ls <- vector("list", 3)
    names(fit_ls) <- c("DX", "Braak", "Thal")
    for (pheno in c("DX", "Braak", "Thal")){
      if (pheno == "DX"){
        print(pheno)
        formula <- paste0("CQN ~ ", pheno, "+ AOD + Sex_M1F0 + RIN + Source + (1|Batch)")
        fit <- br_cqn_long_ls[[region]] %>%
          map(~ lmer(formula, data = .x)) %>% 
          map_dfr(broom.mixed::tidy, conf.int = TRUE, .id = "gene_id") %>%
          filter(term == pheno) %>%
          mutate(q.value = p.adjust(p.value, method = "fdr")) %>%
          dplyr::select(-effect, -group)
        fit_ls[[pheno]] <- fit}
      else {
        print(pheno)
        formula <- paste0("CQN ~ ", pheno, "+ AOD + Sex_M1F0 + RIN + (1|Batch)")
        fit <- br_cqn_long_ls[[region]] %>%
          map(~ lmer(formula, data = .x)) %>% 
          map_dfr(broom.mixed::tidy, conf.int = TRUE, .id = "gene_id") %>%
          filter(term == pheno) %>%
          mutate(q.value = p.adjust(p.value, method = "fdr")) %>%
          dplyr::select(-effect, -group)
        fit_ls[[pheno]] <- fit }
    }
  } else if (region %in% c("FP", "MSSM\nSTG", "PHG", "IFG")){
    fit_ls <- vector("list", 3)
    names(fit_ls) <- c("DX", "Braak", "CERAD")
    for (pheno in c("DX", "Braak", "CERAD")){
      print(pheno)
      formula <- paste0("CQN ~ ", pheno, "+ AOD + Sex_M1F0 + RIN + (1|Batch)")
      fit <- br_cqn_long_ls[[region]] %>%
        map(~ lmer(formula, data = .x)) %>% 
        map_dfr(broom.mixed::tidy, conf.int = TRUE, .id = "gene_id") %>%
        filter(term == pheno) %>%
        mutate(q.value = p.adjust(p.value, method = "fdr")) %>%
        dplyr::select(-effect, -group)
      fit_ls[[pheno]] <- fit
    }
  } else {
    fit_ls <- vector("list", 4)
    names(fit_ls) <- c("DX", "Braak", "Thal")
    for (pheno in c("DX", "Braak", "Thal")){
      print(pheno)
      formula <- paste0("CQN ~ ", pheno, "+ AOD + Sex_M1F0 + RIN + Source + (1|Batch)")
      fit <- br_cqn_long_ls[[region]] %>%
        map(~ lmer(formula, data = .x)) %>% 
        map_dfr(broom.mixed::tidy, conf.int = TRUE, .id = "gene_id") %>%
        filter(term == pheno) %>%
        mutate(q.value = p.adjust(p.value, method = "fdr")) %>%
        dplyr::select(-effect, -group)
      fit_ls[[pheno]] <- fit
    }
  }
  return(fit_ls)
}  

Assoc_MayoSTG <- assoc_br("Mayo\nSTG")
Assoc_CER <- assoc_br("CER")
Assoc_FP <- assoc_br("FP")
Assoc_MSSMSTG <- assoc_br("MSSM\nSTG")
Assoc_PHG <- assoc_br("PHG")
Assoc_IFG <- assoc_br("IFG")
Assoc_DLPFC <- assoc_br("DLPFC")

# Join the results together
## Create empty dfs so all METr lists have the same length for binding later
Assoc_MayoSTG$CERAD <- tibble(gene_id=names(br_cqn_long_ls[[1]]), term="CERAD", estimate=NA, std.error=NA, statistic=NA, p.value=NA, conf.low=NA, conf.high=NA, q.value=NA)
Assoc_CER$CERAD <- tibble(gene_id=names(br_cqn_long_ls[[1]]), term="CERAD", estimate=NA, std.error=NA, statistic=NA, p.value=NA, conf.low=NA, conf.high=NA, q.value=NA)
Assoc_FP$Thal <- tibble(gene_id=names(br_cqn_long_ls[[1]]), term="Thal", estimate=NA, std.error=NA, statistic=NA, p.value=NA, conf.low=NA, conf.high=NA, q.value=NA)
Assoc_MSSMSTG$Thal <- tibble(gene_id=names(br_cqn_long_ls[[1]]), term="Thal", estimate=NA, std.error=NA, statistic=NA, p.value=NA, conf.low=NA, conf.high=NA, q.value=NA)
Assoc_PHG$Thal <- tibble(gene_id=names(br_cqn_long_ls[[1]]), term="Thal", estimate=NA, std.error=NA, statistic=NA, p.value=NA, conf.low=NA, conf.high=NA, q.value=NA)
Assoc_IFG$Thal <- tibble(gene_id=names(br_cqn_long_ls[[1]]), term="Thal", estimate=NA, std.error=NA, statistic=NA, p.value=NA, conf.low=NA, conf.high=NA, q.value=NA)
Assoc_DLPFC$CERAD <- tibble(gene_id=names(br_cqn_long_ls[[1]]), term="CERAD", estimate=NA, std.error=NA, statistic=NA, p.value=NA, conf.low=NA, conf.high=NA, q.value=NA)

Assoc_all_br <- list(Assoc_MayoSTG, Assoc_CER, Assoc_FP, Assoc_MSSMSTG, Assoc_PHG, Assoc_IFG, Assoc_DLPFC)
names(Assoc_all_br) <- region
Assoc_all_br <- map_dfr(Assoc_all_br, ~ bind_rows(.x),.id = 'Region')
Assoc_NoNA_all_br <- drop_na(Assoc_all_br)

write.xlsx(list("withNA"=Assoc_all_br, "noNA"=Assoc_NoNA_all_br), paste0(outdir,"/nduf_assoc_de/NDUF_Br_Assoc.xlsx"))

rm(Assoc_MayoSTG, Assoc_CER, Assoc_FP, Assoc_MSSMSTG, Assoc_PHG, Assoc_IFG, Assoc_DLPFC)


# Mouse ----
## DEG function fixed effect model
run_de <- function(tbl, splitby, formula, filter_term){
  results = tbl %>%
    named_group_split(!!! syms(splitby)) %>% #https://stackoverflow.com/questions/52437463/function-calling-variable-names-for-group-by-in-dplyr-how-do-i-vectorise-this
    map(~ lm(formula, data = .x)) %>% 
    map_dfr(broom::tidy, conf.int = TRUE, .id = "Gene") %>%
    filter(term == filter_term) %>%
    mutate(q.value = p.adjust(p.value, method = "fdr"))
  
  return(results)
} # for nested data
## DEG function for mixed model
run_de_mixed <- function(tbl, splitby, formula, filter_term){
  results = tbl %>%
    named_group_split(!!! syms(splitby)) %>% 
    map(~ lmer(formula, data = .x)) %>% 
    map_dfr(broom.mixed::tidy, effects="fixed", conf.int = TRUE, .id = "Gene") %>%
    filter(term == filter_term) %>%
    mutate(q.value = p.adjust(p.value, method = "fdr"))
  
  return(results)
}

m_bm <- readRDS(paste0(indir, "/mouse_bm.rds")) #mouse database from biomart
m_nduf <- m_bm %>% filter(`MGI symbol` %in% stringr::str_to_title(nduf_bb$gene_name)) #all the 20 NDUFs have homologs in mouse

## ADBXD ----
adbxd <- read.csv(paste0(indir, "/adbxd.csv")) %>% #a table of nduf genes log2-transformed expression from hippocampus of 14months ADBXD; 
  #rows are each mouse columns are strain/genotype/sex/age information and ndufs expression
  #provided by Drs. Catherine Kaczorowski and Amy Dunn
  mutate(Genotype=factor(Genotype, levels=c("Ntg", "ADBXD")),
         Sex = factor(Sex, levels = c("Female", "Male")), 
         Strain=factor(Strain, levels=unique(Strain)))

adbxd_de <- adbxd %>%
  tidyr::pivot_longer(-c(1:5), names_to = "Gene", values_to = "Expr") %>%
  named_group_split(Gene) %>%
  map(~ lmer(Expr ~ Genotype + Sex + (1|Strain), data = .x)) %>% 
  map_dfr(broom.mixed::tidy, conf.int = TRUE, .id = "Gene") %>%
  filter(term == "GenotypeADBXD") %>%
  mutate(q.value = p.adjust(p.value, method = "fdr")) %>%
  dplyr::select(-group, -effect, -df)

## UCI 5xFAD----
UCI5xfad_rnaseq <- read.csv(paste0(indir,"/UCI_5xFAD_RNAseq_metadata_UCI.csv")) %>% #downloaded from syn16798076
  dplyr::select(specimenID, RIN, rnaBatch) %>%
  mutate(ID=gsub("([0-9]+).*$", "\\1", specimenID),
         Region=if_else(grepl("C", specimenID), "cortex", "hippocampus"))

UCI5xfad_tpm <- read_delim(paste0(indir,"/GSE168137_expressionList.txt"), delim = "\t") %>% #downloaded from GSE168137
  mutate(gene_id = gsub("\\..*", "",gene_id))

UCI5xfad_meta <- colnames(UCI5xfad_tpm[,-1]) %>%
  as.data.frame() %>%
  separate(col = 1, into = c("Genotype", "Region", "Age", "Sex", "ID"), sep = "_", remove = F) %>%
  left_join(UCI5xfad_rnaseq %>% dplyr::select(-specimenID), by=c("ID","Region")) %>%
  mutate(Age=as.numeric(stringr::str_remove(Age, "mon")),
         Sex_M1_F0 = if_else(Sex == "Male", 1L, 0L),
         Gt_ad1_ctrl0 = ifelse(Genotype == "5xFAD;BL6", 1L, 0L)) %>%
  dplyr::rename(CombinedID=1)

UCI5xfad_nduf_deg <- UCI5xfad_tpm %>%
  filter(gene_id %in% m_nduf$`Gene stable ID`) %>%
  mutate_at(.vars = c(2:ncol(.)), .funs = ~log(.x + 1, 2)) %>%
  pivot_longer(cols = 2:ncol(.), names_to = "CombinedID", values_to = "Expr") %>%
  inner_join(m_nduf, by=c("gene_id"="Gene stable ID")) %>%
  inner_join(UCI5xfad_meta, by="CombinedID") %>%
  filter(Age==18) %>%
  dplyr::select(gene_id, `MGI symbol`, CombinedID, Genotype:Gt_ad1_ctrl0, Expr) %>%
  tidyr::nest(.by = Region) %>%
  mutate(DEG=map(data, run_de, splitby="MGI symbol", formula="Expr ~ Gt_ad1_ctrl0 + Sex_M1_F0 + RIN", filter_term="Gt_ad1_ctrl0")) %>%
  tidyr::unnest(DEG) %>%
  dplyr::select(-data)

## P301S ----
p301s_meta <- read_delim(paste0(indir,"/GSE90693_metaData_SampleSheet_TPR50_mRNAseq.txt"), delim="\t", col_select = c(2:6)) %>% #downloaded from GSE90693
  separate(col = 5, sep = "\t", into = c("tmp1", "tmp2"), convert = T) 
names(p301s_meta) <- c("ID", "Age", "Strain", "Genotype", "Region","RIN")
p301s_meta <- p301s_meta %>% 
  filter(Strain=="C57BL6") %>% # they used different strains, we only focused on B6 strain
  mutate(ID = as.character(ID),
         Region = if_else(Region=="stem", "brainstem", Region),
         Age=stringr::str_replace(Age, " month", "m"),
         Age=factor(Age, levels=c("3m","6m")),
         Gt_ad1_ctrl0 = ifelse(Genotype == "Tg", 1L, 0L))


expr_ls <- list()
for (i in c("Cortex_6months", "Hippocampus_6months")){
    tmp_fpkm <- read_delim(paste0(indir,"/GSE90693_normalizedFPKM_TPR50_",i,".txt"), delim = "\t") #downloaded from GSE90693
    name <- names(tmp_fpkm)
    tmp_fpkm <- tmp_fpkm %>% separate(col = ncol(.), sep = "\t", into = c("tmp1", "tmp2"), convert = T)
    names(tmp_fpkm) <- c("gene_id", name)
    
    expr_ls[[i]] <- tmp_fpkm
}
rm(i, tmp_fpkm, name)


p301s_nduf_deg <- expr_ls %>%
  map(., ~filter(.x, gene_id %in% m_nduf$`Gene stable ID`)) %>%
  #map(., ~mutate_at(.x, .vars = c(2:ncol(.)), .funs = ~log(.x + 1, 2))) %>% #GEO says it's already log2 transformed!
  map(., ~pivot_longer(.x, cols = 2:ncol(.), names_to = "ID", values_to = "Expr")) %>%
  map(., ~inner_join(.x, m_nduf, by=c("gene_id"="Gene stable ID"))) %>%
  map(., ~inner_join(.x, p301s_meta, by="ID")) %>%
  map(., ~dplyr::select(.x, gene_id, `MGI symbol`, ID, Age:Gt_ad1_ctrl0, Expr)) %>%
  map(., ~tidyr::nest(.x, .by = Region)) %>%
  map(., ~mutate(.x, DEG=map(data, run_de, splitby="MGI symbol", formula="Expr ~ Gt_ad1_ctrl0 + RIN", filter_term="Gt_ad1_ctrl0"))) %>%
  map(., ~tidyr::unnest(.x, DEG)) %>%
  map(., ~dplyr::select(.x, -data))
p301s_nduf_deg_df <- bind_rows(p301s_nduf_deg, .id = "Region") %>% mutate(Region = gsub("_.*","", Region))

## APOE-TR ----
apoe_tr_cqn <- read.csv(paste0(indir, "/APOETR_CQN.csv")) #available on syn21198888
apoe_tr_meta <- read.csv(paste0(indir, "/APOETR_Covars.csv")) %>% #provided by Drs. Na Zhao and Yingxue Ren
  filter(Age==24 & Genotype!="APOE2") %>%
  mutate(Sex_M1_F0 = if_else(Sex == "Male", 1L, 0L), 
         Gt_ad1_ctrl0 = ifelse(Genotype == "APOE4", 1L, 0L),
         BIC_batch = factor(BIC_batch), #there's only 1 level of batch
         Flowcell = factor(Flowcell))

apoe_tr_nduf_deg <- apoe_tr_cqn %>%
  filter(GeneID %in% m_nduf$`Gene stable ID`) %>%
  dplyr::select(GeneID, any_of(apoe_tr_meta$ID)) %>%
  pivot_longer(cols = 2:ncol(.), names_to = "ID", values_to = "Expr") %>%
  inner_join(m_nduf, by=c("GeneID"="Gene stable ID")) %>%
  inner_join(apoe_tr_meta, by="ID") %>%
  dplyr::select(GeneID, `MGI symbol`, ID:Gt_ad1_ctrl0, Expr) %>%
  tidyr::nest(-Age) %>%
  mutate(DEG=map(data, run_de, splitby="MGI symbol", formula="Expr ~ Gt_ad1_ctrl0 + Sex_M1_F0 + RIN + Flowcell", filter_term="Gt_ad1_ctrl0")) %>%
  tidyr::unnest(DEG) %>%
  dplyr::select(-data)

mice_deg <- bind_rows(
  adbxd_de %>% mutate(Model = "ADBXD", Region = "hippocampus",.before = "Gene"),
  UCI5xfad_nduf_deg %>% mutate(Model = "5xFAD", .before = "Region"),
  p301s_nduf_deg_df %>% mutate(Region = str_to_lower(Region)) %>% mutate(Model = "P301S", .before = "Region"),
  apoe_tr_nduf_deg %>% dplyr::select(-Age) %>% mutate(Model = "APOE-TR", Region = "cortex", .before = "Gene")
)
write.csv(mice_deg, paste0(outdir,"/nduf_assoc_de/mice_deg.csv"), row.names = F)
