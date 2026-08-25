# Author: Wei (Adelyn) Tsai; tsai.wei@mayo.edu
# If you are using this script, please cite our study

suppressPackageStartupMessages(library(broom)) #v1.0.12
suppressPackageStartupMessages(library(tidyverse)) #v2.0.0
suppressPackageStartupMessages(library(lmerTest)) #v3.1
suppressPackageStartupMessages(library(WGCNA)) #v1.72-1


indir <- ""
outdir <- ""

named_group_split <- function(.tbl, ...) {
  grouped <- group_by(.tbl, ...)
  names <- rlang::eval_bare(rlang::expr(paste(!!!group_keys(grouped), sep = " / ")))
  
  grouped %>% 
    group_split() %>% 
    rlang::set_names(names)
}


# Build residuals for MCSA----
cqn_raw <- read_delim(paste0(indir, "/PaX108_R01resilience_gene_CQN_neg3_postQC_s105.txt"), delim = "\t") # MCSA cqn file for 18046 genes. Each row is a gene and each column is a donor
mcsa_covars <- read_delim(paste0(indir,"/R01_Resilience_PAXgene_Covars.txt"), delim = "\t") %>% mutate(ptnum=as.character(ptnum)) #covariates file for MCSA, rows are donor columns are covariates
## change to long format 
cqn_long <- cqn_raw %>%
  dplyr::select(-Chromosome, -Start, -End, -Length, -GeneName, -GeneBiotype) %>%
  tidyr::pivot_longer(-GeneId, names_to = "ptnum", values_to = "expression") %>%
  inner_join(dplyr::select(mcsa_covars, ptnum, Sex, RIN, PAXgene_flowcell, PAXgene_Batch, PAXgene_Age, educ), by = "ptnum") %>%
  mutate(PAXgene_Batch = factor(PAXgene_Batch),
         Sex = factor(Sex, levels = c("Male","Female")),
         PAXgene_flowcell = as.factor(PAXgene_flowcell)) 
## create residuals
mcsa_res <- cqn_long %>%
  named_group_split(GeneId) %>%
  purrr::map(~ lm(expression ~ PAXgene_Age + Sex + RIN + PAXgene_flowcell + PAXgene_Batch, data = .x)) %>%
  purrr::map_dfr(broom::augment, .id = "GeneId") %>%
  dplyr::select(GeneId, .resid) %>%
  bind_cols(dplyr::select(cqn_long, ptnum))
saveRDS(mcsa_res,paste0(outdir, "/resids/mcsa_res.rds"))


# Build residuals for ADNI----
adni_mc <- readRDS(paste0(indir,"/ADNI_s91_microarray.rds")) #microarray file for 91 donors from ADNI. This has 10116 protein-coding genes, gene id as column names, donorID as row names
adni_covars <- read_delim(paste0(indir,"/ADNI_s91_covars.txt"), delim = "\t") %>% #covariates file for 91 donors from ADNI. Rows are donors columns are covariates
  mutate(Affy.Plate = factor(Affy.Plate),
         SITE = factor(SITE),
         PTGENDER = factor(PTGENDER, levels = c("Male","Female")))

## calculate random effect estimates for SITES and Affyplate
seq <- c(1:102)
seq_vector <- dput(as.character(seq))
for (i in 1:102){
  itr <- as.numeric(strsplit(seq_vector,",")[[i]])
  print(itr)

  stopCol <- if(itr*100<ncol(adni_mc)) itr*100 else ncol(adni_mc)
  startCol <- itr*100 - 99
  print(startCol)
  print(stopCol)
  
  genEx <- adni_mc[,startCol:stopCol] %>%
    rownames_to_column("SubjectID") %>%
    pivot_longer(!SubjectID, names_to = "Gene", values_to = "GenExpr")
  
  myData <- right_join(adni_covars %>% select(SubjectID, AGE_EXAM, PTGENDER, RIN, Affy.Plate, SITE), genEx)  %>% arrange(Gene)
  ## -1 + Gene: uppresses the global intercept and estimates a separate fixed-effect mean for each gene. 
  ## This is necessary because genes within the same batch have different expression, so a single shared intercept would be inappropriate.
  ##   (1|Affy.Plate) + (1|SITE) : Random intercepts for array plate and collection site. 
  ##      We extract these random effect estimates (ranef) and average them across batches (see affyPlate.REEs, site.REEs below).
  ##      These mean REEs are later added back as fixed effects in the per-gene regression models.
  ##       By doing so, we avoid the computation burden of a full mixed model for each gene individually.
  fit.lmer <- lmer(GenExpr ~ -1 + Gene + AGE_EXAM + PTGENDER + RIN + (1|Affy.Plate) + (1|SITE), data = myData) 
  ## extract random effect
  ranef.fit.lmer <- ranef(fit.lmer)
  
  saveRDS(ranef.fit.lmer, file = paste0(outdir, "/resids_batch_s91/itr_", itr, ".rds"))
}

for(i in 1:102){
  ranef.fit.lmer <- readRDS(paste0(outdir, "/resids_batch_s91/itr_", i, ".rds"))
  if(i==1){
    resid_affy <- cbind.data.frame(ranef.fit.lmer$Affy.Plate)
    resid_site <- cbind.data.frame(ranef.fit.lmer$SITE)
  }else{
    resid_affy <- cbind.data.frame(resid_affy, ranef.fit.lmer$Affy.Plate)
    resid_site <- cbind.data.frame(resid_site, ranef.fit.lmer$SITE)
  }
}

## Calculate mean random effect estimates (REE)
affyPlate.REEs <- rowMeans(resid_affy) %>% data.frame() %>% rownames_to_column("Affy.Plate") %>% mutate(Affy.Plate=as.factor(Affy.Plate))
colnames(affyPlate.REEs)<- c("Affy.Plate", "AffyPlate.meanREE")
site.REEs <- rowMeans(resid_site) %>% data.frame() %>% rownames_to_column("SITE") %>% mutate(SITE=as.factor(SITE))
colnames(site.REEs)<- c("SITE", "SITE.meanREE")

## Merge REE with covars 
covar.v2 <- adni_covars %>% left_join(affyPlate.REEs) %>% left_join(site.REEs)
write_delim(covar.v2, paste0(outdir,"/resids/ADNI_s91_covars.txt"), delim = "\t")

## Transform data into long format
genEx.long <- adni_mc %>% rownames_to_column("SubjectID") %>%  pivot_longer(!SubjectID, names_to = "Gene", values_to = "GenExpr") %>% inner_join(covar.v2)  %>% drop_na()

## calculate residuals
genEx.long_wResids <- genEx.long %>%
  group_by(Gene) %>%
  do(data.frame(., resid = residuals(lm(GenExpr ~ AGE_EXAM + PTGENDER + RIN + AffyPlate.meanREE + SITE.meanREE, data=.))))

adni_res <- genEx.long_wResids %>% dplyr::select(SubjectID, Gene, resid) %>% pivot_wider(names_from = Gene, values_from = resid)

saveRDS(adni_res, file = paste0(outdir, "/resids/adni_res.rds"))

# Combine residual files for later use ----
## we focus on genes commonly present in both MCSA and ADNI datasets
com_expr_genes <- unique(intersect(mcsa_res$GeneId, colnames(adni_res)[2:10117])) # 9744 genes
## filter the residuals to only commonly present genes
mcsa_res_use <- mcsa_res %>% filter(GeneId %in% com_expr_genes) %>% pivot_wider(names_from = "GeneId", values_from = ".resid") %>% column_to_rownames("ptnum") %>% data.matrix() 
adni_res_use <- adni_res %>% column_to_rownames("SubjectID") %>% dplyr::select(matches(com_expr_genes)) %>% data.matrix()

ResList <- list(MCSA = list(data = (mcsa_res_use)),
                ADNI = list(data = (adni_res_use)))
saveRDS(ResList, paste0(outdir, "/resids/ResList.rds"))
checkSets(ResList)

comRes <- lapply(ResList, function(set) set$data) %>%
  do.call(what = rbind, args = .)
saveRDS(comRes, paste0(outdir, "/resids/comRes.rds"))