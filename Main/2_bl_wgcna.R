# Author: Wei (Adelyn) Tsai; tsai.wei@mayo.edu
# If you are using this script, please cite our study


suppressPackageStartupMessages(library(broom)) #v1.0.12
suppressPackageStartupMessages(library(tidyverse)) #v2.0.0
suppressPackageStartupMessages(library(pbapply)) # v1.7
suppressPackageStartupMessages(library(lmerTest)) #v3.1
suppressPackageStartupMessages(library(WGCNA)) #v1.72-1
suppressPackageStartupMessages(library(anRichment)) #v1.22
suppressPackageStartupMessages(library(anRichmentMethods)) #v0.96-86
suppressPackageStartupMessages(library(biomaRt)) #v2.54.1
suppressPackageStartupMessages(library(Biobase)) #v2.58.0
suppressPackageStartupMessages(library(openxlsx)) #v4.2.8
suppressPackageStartupMessages(library(meta)) #v6.5
suppressPackageStartupMessages(library(rrvgo)) #v1.10.0
suppressPackageStartupMessages(library(matrixStats)) #v1.3.0

source("./Codes/fns.R") #source a list of frequently used functions

options(stringsAsFactors = FALSE)
enableWGCNAThreads(6)

indir <- ""
outdir <- ""

DATE = Sys.Date()

ResList <- readRDS(paste0(outdir, "/resids/ResList.rds")) # from 1_CreateResids.R
comRes <- readRDS(paste0(outdir, "/resids/comRes.rds")) # from 1_CreateResids.R
comTrait <- readRDS(paste0(indir, "/comTrait.rds")) # a file having cross-sectional endophenotypes for MCSA (from row 1-105) and ADNI (from row 106-196); row names are donorID
adni_mc_withAnno <- readRDS(paste0(indir, "/ADNI_s91_microarray_withAnno.rds")) #adni microarray file with gene annotation; sample ID as column names
bl_markers <- read.csv(paste0(indir,"/bl_markers.csv")) # Please refer to Supplementary Table 2 of the manuscript

# Calculate soft power ----
ind_sft <- ResList %>%
  lapply(function(set) WGCNA::pickSoftThreshold(
    set$data, 
    corFnc = "bicor",
    corOptions = list(use = 'p', maxPOutliers = 0.1),
    networkType = "signed", 
    RsquaredCut = 0.8,
    verbose = 3)
    )
save(ind_sft, file=paste0(outdir,"/bl_wgcna/ind_sft_",DATE,".RData"))
## organize the results into a dataframe
powers_df <- ind_sft %>% 
  names() %>%
  lapply(function(n) ind_sft[[n]]$fitIndices %>% 
           tibble::add_column(set = n, 
                              powerEstimate = ind_sft[[n]]$powerEstimate)) %>%
  dplyr::bind_rows()
write.csv(powers_df, paste0(outdir, "/bl_wgcna/powers_df_",DATE,".csv"), row.names = F) #soft power of 12 was chosen; the table is in supplementary table 19

# Check adjacency ----
adj_mcsa <- adjacency(comRes[1:105,], power = 12, type="signed", corFnc = "bicor")
adj_adni <- adjacency(comRes[106:196,], power = 12, type="signed", corFnc = "bicor")
## transform to long format
adj_mcsa_l <- adj_mcsa %>% as.data.frame() %>% tibble::rownames_to_column("Gene1") %>% pivot_longer(cols = 2:9744, names_to = "Gene2", values_to = "Corr") %>% mutate(Pair = paste0(Gene1, '_', Gene2)) %>% distinct(Pair, .keep_all = T)
adj_adni_l <- adj_adni %>% as.data.frame() %>% tibble::rownames_to_column("Gene1") %>% pivot_longer(cols = 2:9744, names_to = "Gene2", values_to = "Corr") %>% mutate(Pair = paste0(Gene1, '_', Gene2)) %>% distinct(Pair, .keep_all = T)
save(adj_mcsa_l, adj_adni_l, file=paste0(outdir, "/bl_wgcna/AdjacencyMx_MCSA_ADNI.RData"))

# Build WGCNA ----
## WGCNA parameters
params <- list(maxBlockSize = 30000, corType = "bicor", networkType = "signed", TOMType = "signed", minModuleSize = 30, 
               reassignThreshold = 0, mergeCutHeight = 0.25, maxPOutliers = 0.1, deepSplit=2,numericLabels = TRUE, 
               pamRespectsDendro = FALSE, verbose = 3, nThreads = 6, useMean = FALSE, power = 12,saveIndividualTOMs = FALSE, 
               cacheDir = paste0(outdir, "/bl_wgcna"),saveConsensusTOMs = TRUE, consensusTOMFilePattern = stringr::str_c(outdir, "/bl_wgcna/consensusTOM-block.%b_", DATE, ".RData"))

## Create individual topological overlap matrix (TOM)
sample_size = lapply(ResList, function(df) df$data %>% dim())
params_TOM = list(multiExpr = ResList, maxBlockSize = params$maxBlockSize, power= params$power, corType = params$corType, 
                  networkType = params$networkType,verbose = params$verbose, TOMType = params$TOMType, 
                  nThreads = params$nThreads, individualTOMFileNames = paste0(outdir, "/bl_wgcna/indvTOM_for_set%s(%N)_block%b.RData"))
## Build topological overlap matrix 
indvTOM = do.call(what = WGCNA::blockwiseIndividualTOMs,
                  args = params_TOM)
indvTOM_object = list(indvTOM = indvTOM, sessionInfo = sessionInfo(),
                      params = params_TOM, sample_size = sample_size)
saveRDS(indvTOM_object, paste0(outdir, "/bl_wgcna/indvTOM_object_", DATE, ".rds"))

## Build consensus networks
params_consNet = list(multiExpr = ResList, individualTOMInfo = indvTOM_object$indvTOM)

sample_size = lapply(ResList, function(df) df$data %>% dim()) %>% Reduce(f = "+", x = .) 

consNet = do.call(what = WGCNA::blockwiseConsensusModules,
                  args = c(params_consNet, params))
consNet_object = list(consNet = consNet, sessionInfo = sessionInfo(), 
                      params = c(params_consNet, params), sample_size = sample_size)
saveRDS(consNet_object, paste0(outdir, "/bl_wgcna/consNet_object_", DATE, ".rds"))

# Module preservation within blood consensus networks ----
## Module preservation tests whether consensus modules (derived from comRes) are preserved within each individual cohort.
## consNet uses the full combined expression matrix as the reference network; Discovery (MCSA) and Replication (ADNI) use cohort-specific matrices.
## All three share the same module color assignments from the consensus network.
mp_sp12 = WGCNA::modulePreservation(
  multiData = list(
    consNet = list(data = comRes), # combined MCSA + ADNI
    Discovery = ResList$MCSA, # MCSA only
    Replication = ResList$ADNI # ADNI only
  ),
  multiColor = list(
    consNet = consNet_object$consNet$colors,
    Discovery = consNet_object$consNet$colors,
    Replication = consNet_object$consNet$colors
  ),
  networkType = "signed", corFnc = "bicor", randomSeed = 1, nPermutations=100, verbose=3
)

saveRDS(mp_sp12, paste0(outdir, "/bl_wgcna/mp_SP12_nPerm100_2023-12-06.rds"))

# Module-trait association ----
## Calculate ME, separately for MCSA and ADNI
mcsa_me <- moduleEigengenes(comRes[1:105,], consNet_object$consNet$colors)$eigengenes
adni_me <- moduleEigengenes(comRes[106:196,], consNet_object$consNet$colors)$eigengenes
consMEs1_SP12 <- rbind(mcsa_me, adni_me)
write.csv(consMEs1_SP12 %>% rownames_to_column("SubjectID"), paste0(outdir, "/bl_wgcna/MEs1_SP12.csv"),row.names = F)

## Prepare files for association analysis
### MCSA
comTrait$Mag_HippICV <- factor(comTrait$Mag_HippICV)
consMEs1_mcsa_long_SP12_noHipp <- consMEs1_SP12[1:105,] %>%
  tibble::rownames_to_column("SubjectID") %>%
  pivot_longer(-SubjectID, names_to = "Module", values_to = "ME") %>%
  inner_join(comTrait %>% tibble::rownames_to_column("SubjectID") %>% dplyr::select(SubjectID, Sex, Educ, Age, Memory:AVDEL), by="SubjectID") %>%
  pivot_longer(-c("SubjectID", "Module","ME","Educ", "Age", "Sex"), names_to = "Pheno", values_to = "Measures") %>%
  named_group_split(Module, Pheno)
consMEs1_mcsa_long_SP12_Hipp <- consMEs1_SP12[1:105,] %>%
  tibble::rownames_to_column("SubjectID") %>%
  pivot_longer(-SubjectID, names_to = "Module", values_to = "ME") %>%
  inner_join(comTrait %>% tibble::rownames_to_column("SubjectID") %>% dplyr::select(SubjectID,Sex,Educ,Age,HippVol,Hipp_ICV, Mag_HippICV), by="SubjectID") %>%
  pivot_longer(-c("SubjectID", "Module","ME","Educ","Age", "Sex","Hipp_ICV", "Mag_HippICV"), names_to = "Pheno", values_to = "Measures") %>%
  named_group_split(Module)
### ADNI
consMEs1_adni_long_SP12_noHipp <- consMEs1_SP12[106:196,] %>%
  tibble::rownames_to_column("SubjectID") %>%
  pivot_longer(-SubjectID, names_to = "Module", values_to = "ME") %>%
  inner_join(comTrait %>% tibble::rownames_to_column("SubjectID") %>% dplyr::select(SubjectID, Sex, Educ, Age, Memory:AVDEL), by="SubjectID") %>%
  pivot_longer(-c("SubjectID", "Module","ME","Educ", "Age", "Sex"), names_to = "Pheno", values_to = "Measures") %>%
  named_group_split(Module, Pheno)
consMEs1_adni_long_SP12_Hipp <- consMEs1_SP12[106:196,] %>%
  tibble::rownames_to_column("SubjectID") %>%
  pivot_longer(-SubjectID, names_to = "Module", values_to = "ME") %>%
  inner_join(comTrait %>% tibble::rownames_to_column("SubjectID") %>% dplyr::select(SubjectID,Sex,Educ,Age,HippVol,Hipp_ICV, Mag_HippICV), by="SubjectID") %>%
  pivot_longer(-c("SubjectID", "Module","ME","Educ","Age", "Sex","Hipp_ICV", "Mag_HippICV"), names_to = "Pheno", values_to = "Measures") %>%
  named_group_split(Module)

### run ME-trait association analysis, separately for MCSA and ADNI, then bind the results
fit_mcsa_noHipp_models <- consMEs1_mcsa_long_SP12_noHipp %>% map(~ lm(ME ~ Measures + Educ + Age + Sex, data = .x))
fit_mcsa_noHipp <- fit_mcsa_noHipp_models %>% 
  map_dfr(broom::tidy, conf.int = TRUE, .id = "Module_Pheno") %>%
  filter(term == "Measures") %>%
  separate(Module_Pheno, into=c("Module", "Pheno"), sep = " / ") 
fit_mcsa_hipp_models <- consMEs1_mcsa_long_SP12_Hipp %>% map(~ lm(ME ~ Measures + Educ + Hipp_ICV + Age + Sex, data = .x)) 
fit_mcsa_hipp <- fit_mcsa_hipp_models %>% 
  map_dfr(broom::tidy, conf.int = TRUE, .id = "Module") %>%
  filter(term == "Measures") %>%
  mutate(Pheno="HippVol")
fit_mcsa <- bind_rows(fit_mcsa_noHipp, fit_mcsa_hipp) %>% 
  dplyr::select(-term) %>%
  group_by(Pheno) %>% 
  mutate(q.value = p.adjust(p.value, method = "fdr"))

fit_adni_noHipp_models <- consMEs1_adni_long_SP12_noHipp %>% map(~ lm(ME ~ Measures + Educ + Age + Sex, data = .x))
fit_adni_noHipp <- fit_adni_noHipp_models  %>% 
  map_dfr(broom::tidy, conf.int = TRUE, .id = "Module_Pheno") %>%
  filter(term == "Measures") %>%
  separate(Module_Pheno, into=c("Module", "Pheno"), sep = " / ") 
fit_adni_hipp_models <- consMEs1_adni_long_SP12_Hipp %>% map(~ lm(ME ~ Measures + Educ + Hipp_ICV + Mag_HippICV + Age + Sex, data = .x))
fit_adni_hipp <- fit_adni_hipp_models %>% 
  map_dfr(broom::tidy, conf.int = TRUE, .id = "Module") %>%
  filter(term == "Measures") %>%
  mutate(Pheno="HippVol")
fit_adni <- bind_rows(fit_adni_noHipp, fit_adni_hipp) %>% 
  dplyr::select(-term) %>%
  group_by(Pheno) %>% 
  mutate(q.value = p.adjust(p.value, method = "fdr"))
# save the models for later use
save(list= ls(pattern = "_models"), file = paste0(outdir, "/bl_wgcna/metr_models.RData"))

## Conduct meta-analysis of module eigengene associations
me_trait_joint <- fit_mcsa %>% inner_join(fit_adni, by = c("Module", "Pheno"), suffix=c("_MCSA", "_ADNI"))
### get column indices for estimates and standard errors
eff_idx1 = which(colnames(me_trait_joint) == "estimate_MCSA")
eff_idx2 = which(colnames(me_trait_joint) == "estimate_ADNI")
err_idx1 = which(colnames(me_trait_joint) == "std.error_MCSA")
err_idx2 = which(colnames(me_trait_joint) == "std.error_ADNI")
### run meta-analysis using the defined meta_gen_fn function (meta_gen_fn is sourced from above)
meta_raw = me_trait_joint %>% 
  pbapply::pbapply(1, meta_gen_fn, eff_idx1 = eff_idx1, eff_idx2 = eff_idx2, err_idx1 = err_idx1, err_idx2 = err_idx2)
names(meta_raw) <- paste0(me_trait_joint$Module, "-", me_trait_joint$Pheno)

### Get all fixed-effect results
fix_res = meta_raw %>%
  pbsapply(function(x) x$fixed) %>%
  t() %>%
  as.data.frame() %>%
  rename_with(~paste0(.x, "_fixed")) %>%
  rownames_to_column(var = "Module_Pheno") %>%
  mutate_all(.funs = ~unlist(.x)) %>%
  separate(col=Module_Pheno, into=c('Module', 'Pheno'), sep='-')

### Get all random-effect results
rand_res = meta_raw %>%
  pbsapply(function(x) x$random) %>%
  t() %>%
  data.frame() %>%
  rename_with(~paste0(.x, "_random")) %>%
  rownames_to_column(var = "Module_Pheno") %>%
  mutate_all(.funs = ~unlist(.x)) %>%
  separate(col=Module_Pheno, into=c('Module', 'Pheno'), sep='-')


### Get heterogeneity results
het_res = meta_raw %>%
  map(~.x[c("Q","df.Q","pval.Q","I2", "lower.I2", "upper.I2", "tau2", "se.tau2", "H", "lower.H", "upper.H")]) %>%
  bind_rows(.id = "Module_Pheno") %>%
  separate(col=Module_Pheno, into=c('Module', 'Pheno'), sep='-')

### Join the results to the original regression coeffcients
meta_df = me_trait_joint %>%
  inner_join(fix_res, by=c("Module", "Pheno")) %>%
  inner_join(rand_res, by=c("Module", "Pheno")) %>% 
  rename_with(~sub("^TE", "beta", .x), everything()) %>%
  rename_with(~sub("^seTE", "se", .x), everything()) %>%
  inner_join(het_res, by=c("Module", "Pheno"))

### final meta-analysis file
meta_final = meta_df %>%
  group_by(Pheno) %>%
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
write.csv(meta_final, paste0(outdir, "/bl_wgcna/module_trait.csv"), row.names=F)

# Module membership (MM) ----
## Module membership calculation, this creates columns called kME# (e.g. kME3 = MM of Bl_M3) for each module, while each row is the gene ID
MM_tbl <- consMEs1_SP12 %>% signedKME(datExpr = comRes, datME = ., corFnc = "bicor")
## Pvalues of MM for each calculated module, creates columns called p.MM# (e.g. p.MM3 = p value of MM of Bl_M3),  while each row is the gene ID
MMPvalue <-  as.data.frame(corPvalueStudent(as.matrix(MM_tbl), nrow(comRes)))
names(MMPvalue) <- gsub("^kME", "p.MM", names(MMPvalue))
## Organize the module membership table
MM_tbl <- MM_tbl %>% tibble::rownames_to_column("gene_id") %>% 
  ## join with adni_mc_withAnno, which has mapping between gene name and ensg ID
  inner_join(adni_mc_withAnno %>% dplyr::select(EnsgID_NCBI, GeneSymbol_NCBI), by=c("gene_id"="EnsgID_NCBI")) %>%
  dplyr::rename(gene_name = GeneSymbol_NCBI) %>%
  inner_join(MMPvalue %>% tibble::rownames_to_column("gene_id"), by="gene_id") %>%
  ## create a colume call Module that determine which module a gene belongs to
  mutate(Module=unname(consNet_object$consNet$colors)[match(names(consNet_object$consNet$colors), gene_id)]) %>%
  dplyr::select(gene_id, gene_name, Module, starts_with("kME"), starts_with("p.MM"))
## create a column called MM (i.e. module membership) and fill for each row the module membership number from the kME# columns
MM_tbl$MM = NA 
for(j in 1:nrow(MM_tbl)) {
  ## what this does is that for each row, fetch MM number from kME# columns based on the module number in the Module column
  ## so, e.g. if a gene is in Bl_M3, in the Module column, value = 3, and in the MM column, value = what's in kME3
  MM_tbl[j, "MM"] = MM_tbl[j, paste0("kME", MM_tbl[j,"Module"])]
}
MM_tbl <- MM_tbl %>% dplyr::select(gene_id, gene_name, Module, MM, any_of(paste0("kME", seq_along(0:36)-1)), any_of(paste0("p.MM", seq_along(0:36)-1)))
write.csv(MM_tbl, paste0(outdir, "/bl_wgcna/MM_tbl_", DATE, ".csv"), row.names = F)

# Gene ontology ----
go_collection = buildGOcollection(organism = "human")
mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl") 
bm = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"), bmHeader = TRUE, mart = mart) 
colnames(bm) = c("ensembl_gene_id", "hgnc_symbol", "entrez_gene_id")

wgcna_bm = bm[match(dimnames(comRes)[[2]], bm$ensembl_gene_id), ]
entrez_ids = wgcna_bm$entrez_gene_id

go_enrichment = enrichmentAnalysis(classLabels = consNet_object$consNet$colors,
                                   identifiers = entrez_ids,
                                   refCollection = go_collection,
                                   useBackground = "intersection",
                                   maxReportedOverlapGenes=10000,
                                   threshold = 0.05,
                                   thresholdType = "nominal",
                                   getOverlapEntrez = FALSE,
                                   getOverlapSymbols = TRUE)
collectGarbage()

GOEnrichmentTable_sig = go_enrichment$enrichmentTable %>% filter(FDR < 0.05)   
write.csv(GOEnrichmentTable_sig, paste0(outdir,"/bl_wgcna/GOEnrichmentTable_AllONT_", DATE, ".csv"), row.names = FALSE) # this is in supplementary table 4

## identify genes in Bl_M3 that are also mapped to at least one GO term
m3_gobp_genes <- GOEnrichmentTable_sig %>% 
  filter(class==3 & grepl("BP", inGroups)) %>% 
  dplyr::select(overlapGenes) %>% 
  splitstackshape::cSplit("overlapGenes", "|", direction = "long") %>% 
  distinct(overlapGenes) %>% 
  inner_join(MM_tbl %>% 
               filter(Module ==3) %>% 
               dplyr::select(gene_name, MM), by=c("overlapGenes" = "gene_name"))
nrow(m3_gobp_genes) #265 genes mapped to BP in Bl_M3
length(intersect(subset(MM_tbl, Module == 3 & MM>=0.7)$gene_name, m3_gobp_genes$overlapGenes)) #38 M3 hub genes mapped to a GO BP term
nrow(subset(MM_tbl, Module == 3 & grepl("NDUF", gene_name))) # 20 NDUFs in Bl_M3
nrow(subset(MM_tbl, Module == 3 & MM>=0.7 & grepl("NDUF", gene_name))) #8 NDUFs are hub genes in Bl_M3

# to test enrichment of NDUFs in Bl_M3
## we test: among the 38 hub genes that mapped to a GO BP in Bl_M3, whether it's enriched for NDUFs?
## there are 265 background genes mapped to a GO BP, 20 of them are NDUFs and 38 of them are hub genes; 8 of the 20 NDUFs are "hub NDUF genes mapped to a GO BP"
fisher.test(cbind(c(8, 12),  c(30, 215)), alternative = "greater") # p-value = 0.003008


# Blood cell type enrichment----
cell_types = c("NK", "Monocyte", "CD4T", "CD8T","Neutrophil", "Basophil", "Bcell","DC","Megakaryocytes","Erythroblasts")
ct_enrich_ls <- vector("list", 36)
names(ct_enrich_ls) <- paste0("M", seq_along(1:36))  

for(m in 1:36) { # iterate through each module
  module_genes = names(consNet_object$consNet$colors)[consNet_object$consNet$colors == m] # fetch the genes in the current module
  module_test_df = data.frame(Module = paste0("M",m), ModuleSize = length(module_genes)) # fetch the size of the current module
  
  for(n in 1: length(cell_types)) { # iterate through each cell type
    if (n==1){
      current_cell_type = cell_types[n] # fetch current cell type
      ct_genes_in_bg = intersect(MM_tbl$gene_id, subset(bl_markers, Celltype==current_cell_type)$EnsgID) # fetch the cell type genes that are in the background genes
      shared_genes = intersect(module_genes, ct_genes_in_bg) # fetch genes shared between module genes and cell type genes in the background
      # a contingency table for fisher's test
      contingency = cbind(c(length(shared_genes), (length(module_genes) - length(shared_genes))),
                          c((length(ct_genes_in_bg) - length(shared_genes)), (9743 - length(module_genes) - (length(ct_genes_in_bg) - length(shared_genes))))) #9743 is the total number of genes I've in my data and is set as background gene set
      # intermediate dataframe
      int_results = data.frame(CT = current_cell_type,
                               CTSizeInWGCNA = length(ct_genes_in_bg),
                               ModuleCTSharedSize = length(shared_genes),
                               SharedGenes = paste0(subset(MM_tbl, gene_id %in% shared_genes)$gene_name, collapse = ","),
                               Pvalue = fisher.test(contingency, alternative = "greater")$p.value) 
      module_test_df = cbind(module_test_df, int_results) # has to do cbind for first cell type
    } else{
      current_cell_type = cell_types[n]
      ct_genes_in_bg = intersect(MM_tbl$gene_id, subset(bl_markers, Celltype==current_cell_type)$EnsgID)
      shared_genes = intersect(module_genes, ct_genes_in_bg)
      contingency = cbind(c(length(shared_genes), (length(module_genes) - length(shared_genes))),
                          c((length(ct_genes_in_bg) - length(shared_genes)), (9743 - length(module_genes) - (length(ct_genes_in_bg) - length(shared_genes)))))
      int_results = data.frame(Module = paste0("M",m),
                               ModuleSize = length(module_genes),
                               CT = current_cell_type,
                               CTSizeInWGCNA = length(ct_genes_in_bg),
                               ModuleCTSharedSize = length(shared_genes),
                               SharedGenes = paste0(subset(MM_tbl, gene_id %in% shared_genes)$gene_name, collapse = ","),
                               Pvalue = fisher.test(contingency, alternative = "greater")$p.value) 
      module_test_df = rbind(module_test_df, int_results) # do rbind for the rest
    }
  }
  ct_enrich_ls[[paste0("M",m)]] = module_test_df
}

ct_enrich_df <- bind_rows(ct_enrich_ls) %>% group_by(CT) %>% mutate(FDR = p.adjust(Pvalue, "fdr"))
write.csv(ct_enrich_df, paste0(outdir,"/bl_wgcna/CT.csv"), row.names = FALSE) # this is in supplementary table 3