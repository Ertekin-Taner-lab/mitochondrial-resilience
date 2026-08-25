# Author: Wei (Adelyn) Tsai; tsai.wei@mayo.edu
# If you are using this script, please cite our study

suppressPackageStartupMessages(library(tidyverse)) #v2.0.0
suppressPackageStartupMessages(library(Biobase)) #v2.58.0
suppressPackageStartupMessages(library(pbapply)) # v1.7
suppressPackageStartupMessages(library(biomaRt)) #v2.54.1
suppressPackageStartupMessages(library(anRichment)) #v1.22
suppressPackageStartupMessages(library(anRichmentMethods)) #v0.96-86
suppressPackageStartupMessages(library(BRETIGEA))
suppressPackageStartupMessages(library(broom)) #v1.0.12
suppressPackageStartupMessages(library(lmerTest)) #v3.1
suppressPackageStartupMessages(library(WGCNA)) #v1.72
suppressPackageStartupMessages(library(matrixStats)) #4.2.2
suppressPackageStartupMessages(library(openxlsx)) #v4.2.8
suppressPackageStartupMessages(library(meta)) #v6.5
suppressPackageStartupMessages(library(rrvgo)) #v1.10

options(stringsAsFactors = FALSE)
allowWGCNAThreads(64)

indir <- ""
outdir <- ""

source("./Codes/fns.R") #source a list of frequently used functions

cohort= c("Mayo", "MSSM", "ROSMAP")
region= c("Mayo\nSTG", "CER", "FP", "MSSM\nSTG", "PHG", "IFG", "DLPFC")

comTrait_br <- readRDS(paste0(indir, "/comTrait_br.rds")) # a file containing organized phenotypes & covariates of brain datasets; rows are donors columns are phenotypes/covariates
cqn2use <- readRDS(paste0(indir, "/br_cqn2use.rds")) # a list containing cqn of brain datasets. In each list is a dataframe where rows are genes and columns are donors

# Generate residuals ----
cqn_ls_long <- cqn2use %>%
  map(., ~pivot_longer(.x, cols = 4:ncol(.), names_to = "SampleID", values_to = "cqn"))
br_covars_ls <- comTrait_br %>% named_group_split(Region)

for (i in 1:7){
  cqn_ls_long[[i]] <- cqn_ls_long[[i]] %>% 
    # join with covariates
    left_join(br_covars_ls[[i]] %>% dplyr::select(SampleID:Sex_M1F0), by="SampleID") %>%
    mutate(Batch = factor(Batch), Source = factor(Source))
}

## run residuals generation
resids_lmer = vector("list", 7 )
names(resids_lmer) = names(cqn_ls_long)
resids_lmer[c(1,2,7)] <- lapply(cqn_ls_long, function(df) {
  df %>% group_by(gene_id) %>% do(data.frame(., resid = residuals(lmer(cqn ~ AOD + Sex_M1F0 + RIN + Source + (1|Batch), data=.))))
}) # linear mixed model to get residuals
resids_lmer[c(3:6)] <- lapply(cqn_ls_long, function(df) {
  df %>% group_by(gene_id) %>% do(data.frame(., resid = residuals(lmer(cqn ~ AOD + Sex_M1F0 + RIN + (1|Batch), data=.))))
}) 

res2use_lmer <- resids_lmer %>% 
  map(., ~dplyr::select(.x, gene_id, SampleID, resid)) %>%
  map(., ~pivot_wider(.x, names_from = "SampleID", values_from = "resid")) %>%
  map(., ~tibble::column_to_rownames(.x, "gene_id")) %>%
  map(., ~t(.x))

for (i in 1:7){
  reg <- region[i]
  res2use_lmer[[i]] <- res2use_lmer[[i]][subset(comTrait_br, Region == reg)$SampleID,]
}

# create a list of expression residuals containing seven datasets
multiExpr_br_lmer = vector(mode = "list", length = 7)
names(multiExpr_br_lmer) <- region
multiExpr_br_lmer[[1]] <- list(data=(res2use_lmer[[1]]))
multiExpr_br_lmer[[2]] <- list(data=(res2use_lmer[[2]]))
multiExpr_br_lmer[[3]] <- list(data=(res2use_lmer[[3]])) 
multiExpr_br_lmer[[4]] <- list(data=(res2use_lmer[[4]]))
multiExpr_br_lmer[[5]] <- list(data=(res2use_lmer[[5]]))
multiExpr_br_lmer[[6]] <- list(data=(res2use_lmer[[6]]))
multiExpr_br_lmer[[7]] <- list(data=(res2use_lmer[[7]])) 

#check 
checkSets(multiExpr_br_lmer , checkStructure = T)
saveRDS(multiExpr_br_lmer, paste0(outdir,"/resids/multiExpr_br_lmer.rds"))

comRes_br_lmer <- lapply(multiExpr_br_lmer, function(set) set$data) %>% do.call(what = rbind, args = .)
saveRDS(comRes_br_lmer, paste0(outdir,"/resids/comRes_br_lmer.rds"))

# WGCNA ----
## Calculate soft power ----
ind_sft <- multiExpr_br_lmer %>%
  lapply(function(set) WGCNA::pickSoftThreshold(
    set$data, 
    corFnc = "bicor",
    corOptions = list(use = 'p', maxPOutliers = 0.1),
    networkType = "signed", 
    RsquaredCut = 0.8,
    verbose = 3))
## organize the results into a dataframe
powers_df <- ind_sft %>% 
  names() %>%
  lapply(function(n) ind_sft[[n]]$fitIndices %>% 
           tibble::add_column(set = n, 
                              powerEstimate = ind_sft[[n]]$powerEstimate)) %>%
  dplyr::bind_rows()
write.csv(powers_df, paste0(outdir, "/br_wgcna/powers_df.csv"), row.names = F) # this is in supplementary table 19
collectGarbage()

## Get consensus networks----
power = 12 #based on sft above
### parameters throughout
params <- list(maxBlockSize = 30000, corType = "bicor", networkType = "signed", TOMType = "signed",
               minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, maxPOutliers = 0.1,
               deepSplit=2,numericLabels = TRUE, pamRespectsDendro = FALSE, verbose = 3, nThreads = 6, useMean = FALSE, power = power,saveIndividualTOMs = FALSE,
               cacheDir = paste0(outdir,"/br_wgcna"),saveConsensusTOMs = TRUE, consensusTOMFilePattern = stringr::str_c(outdir, "/br_wgcna/consTOM_block.%b_SP", power, ".RData"))

### create individual TOMs
indvTOM_fname = stringr::str_c(outdir, "/br_wgcna/indvTOM.rds")

### Calculate and Save individual TOM
sample_size = lapply(multiExpr_br_lmer, function(df) df$data %>% dim())
TOM_fname_pattern = stringr::str_c(outdir, "/br_wgcna/indTOM_set%s(%N)_block%b_SP", power, ".RData")
params_TOM = list(multiExpr = multiExpr_br_lmer, maxBlockSize = params$maxBlockSize, power= params$power,
                  corType = params$corType, networkType = params$networkType,verbose = params$verbose, TOMType = params$TOMType, nThreads = params$nThreads,individualTOMFileNames = TOM_fname_pattern)
### Build topological overlap matrix 
indvTOM = do.call(what = WGCNA::blockwiseIndividualTOMs,
                  args = params_TOM)
indvTOM_object = list(indvTOM = indvTOM, sessionInfo = sessionInfo(),
                      params = params_TOM, sample_size = sample_size)
saveRDS(indvTOM_object, indvTOM_fname)

### consensus network
consNet_fname = stringr::str_c(outdir, "/br_wgcna/consNet_SP",power, ".rds")
params_consNet = list(multiExpr = multiExpr_br_lmer, individualTOMInfo = indvTOM_object$indvTOM)

sample_size = lapply(multiExpr_br_lmer, function(df) df$data %>% dim()) %>%
  Reduce(f = "+", x = .)

consNet = do.call(what = WGCNA::blockwiseConsensusModules,
                  args = c(params_consNet, params))
consNet_object = list(consNet = consNet, sessionInfo = sessionInfo(),
                      params = c(params_consNet, params), sample_size = sample_size)
saveRDS(consNet_object, consNet_fname)

## ME & MM----
ncbi_geneDB_full <- read_delim(paste0(indir, "/NCBI_Homo_sapiens_gene_info_20221020.txt"), delim = "\t") # a table of genes information downloaded from NCBI
ncbi_geneDB <- ncbi_geneDB_full %>% dplyr::mutate(ensgID = str_extract(dbXrefs, "ENSG\\d{11}"), .after = dbXrefs) %>% dplyr::mutate(GeneID = as.character(GeneID)) %>% dplyr::select(Symbol, ensgID, type_of_gene)  
br_genes_ncbi_protcoding <- ncbi_geneDB %>% filter(ensgID %in% colnames(comRes_br_lmer)) %>% group_by(ensgID) %>% slice_head()  # to get the Gene symbol based on ensgID


consNet_SP12 <- readRDS(paste0(outdir,"/br_wgcna/consNet_SP12.rds"))
### calculate module eigengenes for each dataset
consMEs_br=data.frame()
for (i in region){
  tmp=moduleEigengenes(comRes_br_lmer[rownames(comRes_br_lmer)%in%subset(comTrait_br, Region==i)$SampleID,], consNet_SP12$consNet$colors)$eigengenes
  
  consMEs_br=rbind(consMEs_br,tmp)
}
write.csv(consMEs_br %>% rownames_to_column("SubjectID"), paste0(outdir, "/br_wgcna/MEs1_SP12_Br.csv"), row.names = F)

### Module membership calculation, this creates columns called kME# (e.g. kME1 = MM of Br_M1) for each module, while each row is the gene ID
MM_tbl_br <- consMEs_br %>% signedKME(datExpr = comRes_br_lmer, datME = ., corFnc = "bicor")
### Pvalues of MM for each calculated module, creates columns called p.MM# (e.g. p.MM1 = p value of MM of Br_M1),  while each row is the gene ID
MMPvalue_br <-  as.data.frame(corPvalueStudent(as.matrix(MM_tbl_br), nrow(comRes_br_lmer)))
names(MMPvalue_br) <- gsub("^kME", "p.MM", names(MMPvalue_br))
### Organize the module membership table
MM_tbl_br <- MM_tbl_br %>% tibble::rownames_to_column("gene_id") %>% 
  ### join with gene annotation
  inner_join(br_genes_ncbi_protcoding %>% dplyr::select(ensgID, Symbol), by=c("gene_id"="ensgID")) %>%
  inner_join(MMPvalue_br %>% tibble::rownames_to_column("gene_id"), by="gene_id") %>%
  ### create a colume call Module that determine which module a gene belongs to
  mutate(Module=unname(consNet_SP12$consNet$colors)[match(names(consNet_SP12$consNet$colors), gene_id)]) %>%
  dplyr::select(gene_id, Symbol, Module, starts_with("kME"), starts_with("p.MM"))
### create a column called MM (i.e. module membership) and fill for each row the module membership number from the kME# columns
MM_tbl_br$MM = NA 
for(j in 1:nrow(MM_tbl_br)) {
  ## what this does is that for each row, fetch MM number from kME# columns based on the module number in the Module column
  ## so, e.g. if a gene is in Br_M1, in the Module column, value = 1, and in the MM column, value = what's in kME1
  MM_tbl_br[j, "MM"] = MM_tbl_br[j, paste0("kME", MM_tbl_br[j,"Module"])]
}
MM_tbl_br <- MM_tbl_br %>% dplyr::select(gene_id, Symbol, Module, MM, any_of(paste0("kME", seq_along(0:max(consNet_SP12$consNet$colors))-1)), any_of(paste0("p.MM", seq_along(0:max(consNet_SP12$consNet$colors))-1)))
write.csv(MM_tbl_br, paste0(outdir, "/br_wgcna/MM_tbl_br.csv"), row.names = F)

## Gene ontology ----
go_collection = buildGOcollection(organism = "human")
mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl") 
bm = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"), bmHeader = TRUE, mart = mart) 
colnames(bm) = c("ensembl_gene_id", "hgnc_symbol", "entrez_gene_id")

wgcna_bm = bm[match(dimnames(comRes_br_lmer)[[2]], bm$ensembl_gene_id), ]

go_enrichment = enrichmentAnalysis(classLabels = consNet_SP12$consNet$colors,
                                   identifiers = wgcna_bm$`NCBI gene (formerly Entrezgene) ID`,
                                   refCollection = go_collection,
                                   useBackground = "intersection",
                                   maxReportedOverlapGenes=10000,
                                   threshold = 0.05,
                                   thresholdType = "nominal",
                                   getOverlapEntrez = FALSE,
                                   getOverlapSymbols = TRUE)
collectGarbage()

GOEnrichmentTable = go_enrichment$enrichmentTable %>% filter(FDR < 0.05) 
write.csv(as_tibble(GOEnrichmentTable), paste0(outdir,"/br_wgcna/GO_br.csv"), row.names = FALSE) # this is in supplementary table 4

## Cell type enrichment ----
bl_markers <- read.csv(paste0(indir,"/bl_markers.csv")) # Please refer to Supplementary Table 2 of the manuscript
data("markers_df_brain") # calling from BRETIGEA

cell_types = c("NK", "Monocyte", "CD4T", "CD8T","Neutrophil", "Basophil", "Bcell","DC","Ast", "End", "Mic", "Neu", "Oli")
ct_enrich_ls <- vector("list", max(consNet_SP12$consNet$colors))
names(ct_enrich_ls) <- paste0("M", seq_along(1:max(consNet_SP12$consNet$colors)))  

for(m in 1:max(consNet_SP12$consNet$colors)) {
  module_genes = names(consNet_SP12$consNet$colors)[consNet_SP12$consNet$colors == m] # fetch the genes in the current module
  module_test_df = data.frame(Module = paste0("M",m), ModuleSize = length(module_genes)) # fetch the size of the current module
  
  for(n in 1: length(cell_types)) { # iterate through each cell type
    if (n==1){
      current_cell_type = cell_types[n] # fetch current cell type
      ct_genes_in_bg = intersect(MM_tbl_br$gene_id, get(current_cell_type)) # fetch the cell type genes that are in the background genes
      shared_genes = intersect(module_genes, ct_genes_in_bg) # fetch genes shared between module genes and cell type genes in the background
      # a contingency table for fisher's test
      contingency = cbind(c(length(shared_genes), (length(module_genes) - length(shared_genes))),
                          c((length(ct_genes_in_bg) - length(shared_genes)), (nrow(MM_tbl_br) - length(module_genes) - (length(ct_genes_in_bg) - length(shared_genes)))))
      # intermediate dataframe
      int_results = data.frame(CT = current_cell_type,
                               CTSizeInWGCNA = length(ct_genes_in_bg),
                               ModuleCTSharedSize = length(shared_genes),
                               SharedGenes = paste0(subset(MM_tbl_br, gene_id %in% shared_genes)$Symbol, collapse = ","),
                               Pvalue = fisher.test(contingency, alternative = "greater")$p.value) 
      module_test_df = cbind(module_test_df, int_results) # has to cbind for the first cell type
    } else{
      current_cell_type = cell_types[n]
      ct_genes_in_bg = intersect(MM_tbl_br$gene_id, get(current_cell_type))
      shared_genes = intersect(module_genes, ct_genes_in_bg)
      contingency = cbind(c(length(shared_genes), (length(module_genes) - length(shared_genes))),
                          c((length(ct_genes_in_bg) - length(shared_genes)), (nrow(MM_tbl_br) - length(module_genes) - (length(ct_genes_in_bg) - length(shared_genes)))))
      int_results = data.frame(Module = paste0("M",m),
                               ModuleSize = length(module_genes),
                               CT = current_cell_type,
                               CTSizeInWGCNA = length(ct_genes_in_bg),
                               ModuleCTSharedSize = length(shared_genes),
                               SharedGenes = paste0(subset(MM_tbl_br, gene_id %in% shared_genes)$Symbol, collapse = ","),
                               Pvalue = fisher.test(contingency, alternative = "greater")$p.value) 
      module_test_df = rbind(module_test_df, int_results) # do rbind for the rest
    }
  }
  ct_enrich_ls[[paste0("M",m)]] = module_test_df
}

ct_enrich_df <- bind_rows(ct_enrich_ls) %>% group_by(CT) %>% mutate(FDR = p.adjust(Pvalue, "fdr"))
write.csv(ct_enrich_df, paste0(outdir,"/br_wgcna/CT_42modules.csv"), row.names = FALSE) # this is in supplementary table 3

## ME-trait ----
consMEs_br_long_sub <- consMEs_br %>%
  dplyr::select(ME1, ME17, ME26) %>%
  tibble::rownames_to_column("SampleID") %>%
  pivot_longer(-SampleID, names_to = "Module", values_to = "ME") %>%
  inner_join(comTrait_br, by="SampleID") %>%
  named_group_split(Region) %>% 
  map(., ~mutate(.x, 
                 DX=if_else(DX=="Control", 1L, 0L), 
                 Batch = factor(Batch, levels=unique(Batch)),
                 Source = factor(Source, levels=unique(Source)))) %>%
  map(., ~named_group_split(.x, Module)) # this creates a list containing 7x3 dataframes

# function to run ME-trait association in brain data
me_tr_br <- function(region){
  if (region %in% c("Mayo\nSTG", "CER", "DLPFC")){
    fit_ls <- vector("list", 3)
    names(fit_ls) <- c("DX", "Braak", "Thal") 
    for (pheno in c("DX", "Braak", "Thal")){ 
      print(pheno)
      formula <- paste0("ME ~ ", pheno, "+ AOD + Sex_M1F0")
      fit <- consMEs_br_long_sub[[region]] %>%
        map(~ lm(formula, data = .x)) %>% 
        map_dfr(broom::tidy, conf.int = TRUE, .id = "Module") %>%
        filter(term == pheno) %>%
        mutate(q.value = p.adjust(p.value, method = "fdr"))
      fit_ls[[pheno]] <- fit
    }
  } else{
    fit_ls <- vector("list", 3)
    names(fit_ls) <- c("DX", "Braak", "CERAD")
    for (pheno in c("DX", "Braak", "CERAD")){
      print(pheno)
      formula <- paste0("ME ~ ", pheno, "+ AOD + Sex_M1F0")
      fit <- consMEs_br_long_sub[[region]] %>%
        map(~ lm(formula, data = .x)) %>% 
        map_dfr(broom::tidy, conf.int = TRUE, .id = "Module") %>%
        filter(term == pheno) %>%
        mutate(q.value = p.adjust(p.value, method = "fdr"))
      fit_ls[[pheno]] <- fit
    }
  } 
  return(fit_ls)
}  


METr_MayoSTG <- me_tr_br("Mayo\nSTG")
METr_CER <- me_tr_br("CER")
METr_FP <- me_tr_br("FP")
METr_MSSMSTG <- me_tr_br("MSSM\nSTG")
METr_PHG <- me_tr_br("PHG")
METr_IFG <- me_tr_br("IFG")
METr_DLPFC <- me_tr_br("DLPFC")


# Join the results together
## Create empty dfs so all METr lists have the same length for binding later
METr_MayoSTG$CERAD <- tibble(Module=c("ME26","ME17", "ME1"), term="CERAD", estimate=NA, std.error=NA, statistic=NA, p.value=NA, conf.low=NA, conf.high=NA, q.value=NA)
METr_CER$CERAD <- tibble(Module=c("ME26","ME17", "ME1"), term="CERAD", estimate=NA, std.error=NA, statistic=NA, p.value=NA, conf.low=NA, conf.high=NA, q.value=NA)
METr_FP$Thal <- tibble(Module=c("ME26","ME17", "ME1"), term="Thal", estimate=NA, std.error=NA, statistic=NA, p.value=NA, conf.low=NA, conf.high=NA, q.value=NA)
METr_MSSMSTG$Thal <- tibble(Module=c("ME26","ME17", "ME1"), term="Thal", estimate=NA, std.error=NA, statistic=NA, p.value=NA, conf.low=NA, conf.high=NA, q.value=NA)
METr_PHG$Thal <- tibble(Module=c("ME26","ME17", "ME1"), term="Thal", estimate=NA, std.error=NA, statistic=NA, p.value=NA, conf.low=NA, conf.high=NA, q.value=NA)
METr_IFG$Thal <- tibble(Module=c("ME26","ME17", "ME1"), term="Thal", estimate=NA, std.error=NA, statistic=NA, p.value=NA, conf.low=NA, conf.high=NA, q.value=NA)
METr_DLPFC$CERAD <- tibble(Module=c("ME26","ME17", "ME1"), term="CERAD", estimate=NA, std.error=NA, statistic=NA, p.value=NA, conf.low=NA, conf.high=NA, q.value=NA)

METr_all <- list(METr_MayoSTG, METr_CER, METr_FP, METr_MSSMSTG, METr_PHG, METr_IFG, METr_DLPFC)
names(METr_all) <- c("Mayo\nSTG", "CER", "FP", "MSSM\nSTG", "PHG", "IFG", "DLPFC")
METr_all <- map_dfr(METr_all, ~ bind_rows(.x),.id = 'Region')
METr_allNoNA <- drop_na(METr_all)

save(METr_all, METr_allNoNA, file=paste0(outdir,"/br_wgcna/METr.RData"))
rm(list=ls(pattern = "^METr_"))

## Preservation ----
### Preservation within brain consensus networks
### Module preservation tests whether consensus modules (derived from comRes_br_lmer) are preserved within each individual cohort.
### consNet uses the full combined expression matrix as the reference network; the rest use cohort-specific matrices.
### All share the same module color assignments from the consensus network.
mp = WGCNA::modulePreservation(
  multiData = list(
    consNet = list(data = comRes_br_lmer),
    MayoSTG = multiExpr_br_lmer[[1]],
    CER = multiExpr_br_lmer[[2]],
    FP = multiExpr_br_lmer[[3]],
    MSSMSTG = multiExpr_br_lmer[[4]],
    PHG = multiExpr_br_lmer[[5]],
    IFG = multiExpr_br_lmer[[6]],
    DLPFC = multiExpr_br_lmer[[7]]
  ),
  multiColor = list(
    consNet = consNet_SP12$consNet$colors
  ),
  networkType = "signed", corFnc = "bicor", randomSeed = 1, nPermutations=100
)

saveRDS(mp, paste0(outdir, "/br_wgcna/mp_br.rds"))

### Preservation of blood modules in the brain networks calculated by WGCNA::modulePreservation
comRes <- readRDS(paste0(outdir,"/resids/comRes.rds"))
consNet_bl <- readRDS(paste0(outdir,"/bl_wgcna/consNet_object_2023-08-15.rds"))

mp_bb = WGCNA::modulePreservation(
  multiData = list(
    Blood = list(data = comRes),
    Brain = list(data = comRes_br_lmer)
  ),
  multiColor = list(
    Blood = consNet_bl$consNet$colors
  ),
  networkType = "signed", corFnc = "bicor", randomSeed = 1, nPermutations=100, verbose=3
)
saveRDS(mp_bb, paste0(outdir, "/br_wgcna/mp_bb.rds"))

### Preservation of blood modules in brain networks calculated by fisher.test
MM_tbl <- read.csv(paste0(outdir, "/bl_wgcna/MM_tbl_2023-08-15.csv"))
brain_blood_shared <- intersect(MM_tbl$gene_id, MM_tbl_br$gene_id) #8957, background genes set

brain_blood_df = data.frame()
for(m in seq_along(1:max(MM_tbl$Module))) { # iterate through each blood module
  current_blood_gene = MM_tbl$gene_id[MM_tbl$Module == m] # fetch current blood module genes
  blood_in_shared = intersect(current_blood_gene, brain_blood_shared)# blood module genes in the background
  blood_in_shared_count = length(blood_in_shared)
  all_gene_count = length(brain_blood_shared) 
  
  module_test_df = data.frame(BloodModule = paste0("M",m),
                              BloodModuleSize = length(current_blood_gene),
                              BloodModuleSharedSize = blood_in_shared_count)
  
  # Iterate through each brain modules
  for(n in seq_along(1:max(MM_tbl_br$Module))) {
    current_brain_gene =  MM_tbl_br$gene_id[MM_tbl_br$Module == n] # fetch current brain module genes
    brain_in_shared = intersect(current_brain_gene, brain_blood_shared) # brain module genes in the background
    brain_in_shared_count = length(brain_in_shared)
    both_modules = length(intersect(blood_in_shared, brain_in_shared)) # number of genes in both modules
    
    # Hypergeometric contigency table
    Hypergeometric = cbind(c(both_modules, (blood_in_shared_count - both_modules)),
                           c((brain_in_shared_count - both_modules), (all_gene_count - blood_in_shared_count - (brain_in_shared_count - both_modules))))
    
    # intermediate dataframe
    int_results = data.frame(BrainModule=paste0("M",n),
                             BrainModuleSize = length(current_brain_gene),
                             BrainModuleSharedSize = brain_in_shared_count,
                             BloodBrainSharedSize = both_modules,
                             BloodBrainSharedGenes = paste0(subset(MM_tbl, gene_id %in% intersect(blood_in_shared, brain_in_shared))$gene_name, collapse = ","),
                             Pvalue = fisher.test(Hypergeometric, alternative = "greater")$p.value)
    if (n==1){
      int_results2 = cbind(module_test_df, int_results)
    }else{
      int_results3 = cbind(module_test_df, int_results)
      int_results2 = rbind(int_results2, int_results3)
      # int_results2 accumulates per-brain-module results for the current blood module
      # int_results3 is a single-row intermediate before binding
    }
    
  }
  
  brain_blood_df = rbind(brain_blood_df, int_results2)
}

brain_blood_df <- brain_blood_df %>%
  group_by(BloodModule) %>% #this is correcting for number of brain modules 
  mutate(FDR=p.adjust(Pvalue, method="fdr"))

write.csv(brain_blood_df, paste0(outdir,"/br_wgcna/ConsBBComp.csv"), row.names = F) # this is in supplementary table 9