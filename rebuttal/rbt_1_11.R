# Author: Wei (Adelyn) Tsai; tsai.wei@mayo.edu
# If you are using this script, please cite our study
# This script addresses comment 11 about cross-platform harmonization from reviewer 1

suppressPackageStartupMessages(library(tidyverse)) #v2.0.0
suppressPackageStartupMessages(library(sva)) #v3.46.0
suppressPackageStartupMessages(library(patchwork)) #v1.3.2
suppressPackageStartupMessages(library(WGCNA)) #v1.72
suppressPackageStartupMessages(library(ggVennDiagram)) #v1.5.2
suppressPackageStartupMessages(library(broom)) #v1.0.12
suppressPackageStartupMessages(library(meta)) #v6.5
suppressPackageStartupMessages(library(pbapply)) #v1.7

source("./Codes/fns.R") #source a list of frequently used functions

#######################################
# list of input/output directories
indir <- "./indir"
residir <- "./outdir/resids"
blwgcna_dir <- "./outdir/bl_wgcna"

outdir <- "./Rebuttal/outdir"

#######################################
# prepare files for MCSA and ADNI for Combat adjustment
cqn_raw <- read_delim(paste0(indir, "/PaX108_R01resilience_gene_CQN_neg3_postQC_s105.txt"), delim = "\t") # MCSA cqn file for 18046 genes. Each row is a gene and each column is a donor
adni_mc <- readRDS(paste0(indir,"/ADNI_s91_microarray.rds")) #microarray file for 91 donors from ADNI. This has 10116 protein-coding genes, gene id as column names, donorID as row names

# get the intersecting genes of two expression files
int_genes <- intersect(subset(cqn_raw, GeneBiotype == "protein_coding")$GeneId, colnames(adni_mc))

# prepare combat input expression files and batch variable
combat_expr <- cqn_raw %>%
  filter(GeneId %in% int_genes) %>% 
  select(-Chromosome, -Start, -End, -Length, -GeneName, -GeneBiotype) %>%
  left_join(
    adni_mc %>% 
      t() %>% 
      as.data.frame() %>% 
      rownames_to_column("GeneId")
  ) %>% 
  column_to_rownames("GeneId")
batch <- c(rep(1, 105), rep(2, 91))  

# combat adjustment
combat_adj <- ComBat(dat = combat_expr, batch = batch)
# save combat adjusted expression file
saveRDS(combat_adj, paste0(outdir,"/1_11/combat_adj_expr.rds"))
# plot to check the distribution of expression before and after adjustment 
p1 <- combat_expr %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(Cohort= c(rep("MCSA", 105), rep("ADNI", 91))) %>% 
  mutate(Cohort = factor(Cohort, levels = c("MCSA", "ADNI"))) %>% 
  pivot_longer(-Cohort, names_to = "gene", values_to = "Expression") %>% 
  ggplot(aes(x=Expression)) + 
  geom_histogram(aes(fill = Cohort), alpha = 0.2, bins = 50) +
  theme_classic() +
  theme(
    text = element_text(size=7),
    legend.key.height = unit(0.2, 'in'),
    legend.key.width = unit(0.1, 'in')
  ) +
  labs(title = "Before Combat adjustment")
p2 <- combat_adj %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(Cohort= c(rep("MCSA", 105), rep("ADNI", 91))) %>% 
  mutate(Cohort = factor(Cohort, levels = c("MCSA", "ADNI"))) %>% 
  pivot_longer(-Cohort, names_to = "gene", values_to = "Expression") %>% 
  ggplot(aes(x=Expression)) + 
  geom_histogram(aes(fill = Cohort), alpha = 0.2, bins = 50) +
  theme_classic() +
  theme(
    text = element_text(size=7),
    legend.key.height = unit(0.2, 'in'),
    legend.key.width = unit(0.1, 'in')
  ) +
  labs(title = "After Combat adjustment")
pdf(paste0(outdir,"/1_11/ExprDistBe4AftCombat.pdf"), width = 7, height = 3)
p1 + p2 + plot_layout(guides = "collect")
dev.off() # this is in the first rebuttal Rebuttal Fig.16
rm(p1, p2)
#######################################

#######################################
# generate new residuals adjusting for covariates in each dataset

# MCSA
mcsa_covars <- read_delim(paste0(indir,"/R01_Resilience_PAXgene_Covars.txt"), delim = "\t") %>% mutate(ptnum=as.character(ptnum)) #covariates file for MCSA, rows are donor columns are covariates
mcsa_covars <- mcsa_covars %>% 
  mutate(
    Sex = factor(Sex),
    PAXgene_flowcell = factor(PAXgene_flowcell), 
    PAXgene_Batch = factor(PAXgene_Batch) # factorize flowcell and batch
  )
mcsa_res <- combat_adj[,1:105] %>% 
  t() %>% 
  as.data.frame() %>%
  rownames_to_column("SubjectID") %>%
  pivot_longer(!SubjectID, names_to = "Gene", values_to = "GenExpr") %>% 
  inner_join(mcsa_covars, by=c("SubjectID" = "ptnum")) %>% 
  group_by(Gene) %>%
  do(data.frame(., resid = residuals(lm(GenExpr ~ PAXgene_Age + Sex + RIN + PAXgene_flowcell + PAXgene_Batch, data=.)))) %>% 
  dplyr::select(SubjectID, Gene, resid) %>% pivot_wider(names_from = Gene, values_from = resid)

adni_covars <- read_delim(paste0(residir, "/ADNI_s91_covars.txt"), delim = "\t") # a covariate file for ADNI, generated in 1_CreateResids.R
adni_covars <- adni_covars %>% mutate(PTGENDER = factor(PTGENDER))
adni_res <- combat_adj[,106:196] %>% 
  t() %>% 
  as.data.frame() %>%
  rownames_to_column("SubjectID") %>%
  pivot_longer(!SubjectID, names_to = "Gene", values_to = "GenExpr") %>% 
  inner_join(adni_covars) %>% 
  drop_na() %>% 
  group_by(Gene) %>%
  do(data.frame(., resid = residuals(lm(GenExpr ~ AGE_EXAM + PTGENDER + RIN + AffyPlate.meanREE + SITE.meanREE, data=.)))) %>% 
  dplyr::select(SubjectID, Gene, resid) %>% pivot_wider(names_from = Gene, values_from = resid)

# save the residuals for later use
mcsa_res <- mcsa_res %>% column_to_rownames("SubjectID") %>% data.matrix() 
adni_res <- adni_res %>% column_to_rownames("SubjectID") %>% data.matrix()
ResList_Combat <- list(MCSA = list(data = (mcsa_res)),
                ADNI = list(data = (adni_res)))
saveRDS(ResList_Combat, paste0(outdir, "/1_11/ResList_Combat.rds"))
checkSets(ResList_Combat)
comRes_Combat <- lapply(ResList_Combat, function(set) set$data) %>% do.call(what = rbind, args = .)
saveRDS(comRes_Combat, paste0(outdir, "/1_11/comRes_Combat.rds"))

# compare the residuals before and after combat
comRes <- readRDS(paste0(residir, "/comRes.rds")) # from 1_CreateResids.R
pdf(paste0(outdir,"/1_11/ResidsDistBe4AftCombat.pdf"), width = 5, height = 3)
comRes %>% 
  as.data.frame() %>% 
  rownames_to_column("SubjectID") %>% 
  mutate(Type = "Original", Cohort = c(rep("MCSA", 105), rep("ADNI", 91))) %>% 
  bind_rows(
    comRes_Combat %>% 
      as.data.frame() %>% 
      rownames_to_column("SubjectID") %>% 
      mutate(Type = "Combat", Cohort = c(rep("MCSA", 105), rep("ADNI", 91)))
  ) %>% 
  pivot_longer(
    cols = starts_with("ENSG"), names_to = "Gene", values_to = "resids"
  ) %>% 
  mutate(Cohort = factor(Cohort, levels = c("MCSA", "ADNI"))) %>% 
  ggplot(aes(x=resids)) + 
  geom_histogram(aes(fill = Type), alpha = 0.2, bins = 50) +
  scale_x_continuous(limits = c(-2, 2)) +
  facet_grid(~Cohort, scales = "free", space = "free") +
  theme_classic() +
  theme(
    text = element_text(size=7),
    legend.key.height = unit(0.2, 'in'),
    legend.key.width = unit(0.1, 'in')
  )
dev.off() # this is in the first rebuttal Rebuttal Fig.17
#######################################

#######################################
# Conduct new WGCNA using the combat adjust
options(stringsAsFactors = FALSE)
enableWGCNAThreads(6)

# Calculate soft power
ind_sft_Combat <- ResList_Combat %>%
  lapply(function(set) WGCNA::pickSoftThreshold(
    set$data, 
    corFnc = "bicor",
    corOptions = list(use = 'p', maxPOutliers = 0.1),
    networkType = "signed", 
    RsquaredCut = 0.8,
    verbose = 3)
  )
powers_df_Combat <- ind_sft_Combat %>% 
  names() %>%
  lapply(function(n) ind_sft_Combat[[n]]$fitIndices %>% 
           tibble::add_column(set = n, 
                              powerEstimate = ind_sft_Combat[[n]]$powerEstimate)) %>%
  dplyr::bind_rows()
write.csv(powers_df_Combat, paste0(outdir, "/1_11/powers_df_Combat.csv"), row.names = F) 
# plot the soft power for each cohort
p1 <- powers_df_Combat %>% 
  mutate(set = factor(set, levels=c("MCSA", "ADNI"))) %>% 
  ggplot(aes(x=Power, y= SFT.R.sq)) +
  geom_text(aes(label=Power, color=set), size=6/.pt) +
  scale_y_continuous(breaks = seq(0, 1, by=0.2), limits=c(0,1)) +
  labs(title="Scale Free Topology Model Fit", x="Scale Free Topology Model Fit", y="Soft Threshold Power") +
  theme_classic() +
  theme(text=element_text(size=7),
        legend.position = "none",
        plot.title = element_text(hjust=0.5, size=7),
        panel.grid.major.x = element_blank() ,
        panel.grid.major.y = element_line(linewidth=.1, color="black"), 
        panel.background = element_rect(color = "black", linewidth = 0.5))
p2 <- powers_df_Combat %>% 
  ggplot(aes(x=Power, y= mean.k.)) +
  geom_text(aes(label=Power, color=set), size=6/.pt) +
  labs(title="Mean Connectivity", x="Soft Threshold Power", y="Mean Connectivity") +
  theme_classic() +
  theme(text=element_text(size=7),
        plot.title = element_text(hjust=0.5, size=7),
        panel.grid.major.x = element_blank() , #set vertical line blank
        panel.grid.major.y = element_line(linewidth=.1, color="black"), #explicitly set horizontal line
        panel.background = element_rect(color = "black", linewidth = 0.5))
pdf(paste0(outdir,"/1_11/sft_combat.pdf"), width = 6, height = 2.5)
p1 + p2 # this is now in Extended Data Fig.9b
dev.off()
# will use soft power of 12

#######################################
# Build WGCNA
## WGCNA parameters
params <- list(maxBlockSize = 30000, corType = "bicor", networkType = "signed", TOMType = "signed", minModuleSize = 30, 
               reassignThreshold = 0, mergeCutHeight = 0.25, maxPOutliers = 0.1, deepSplit=2,numericLabels = TRUE, 
               pamRespectsDendro = FALSE, verbose = 3, nThreads = 6, useMean = FALSE, power = 12,saveIndividualTOMs = FALSE, 
               cacheDir = paste0(outdir, "/1_11"),saveConsensusTOMs = TRUE, consensusTOMFilePattern = stringr::str_c(outdir, "/1_11/consensusTOM-block.%b_20251018.RData"))

## Create individual topological overlap matrix (TOM)
sample_size = lapply(ResList_Combat, function(df) df$data %>% dim())
params_TOM = list(multiExpr = ResList_Combat, maxBlockSize = params$maxBlockSize, power= params$power, corType = params$corType, 
                  networkType = params$networkType,verbose = params$verbose, TOMType = params$TOMType, 
                  nThreads = params$nThreads, individualTOMFileNames = paste0(outdir, "/1_11/indvTOM_for_set%s(%N)_block%b_20251018.RData"))

indvTOM = do.call(what = WGCNA::blockwiseIndividualTOMs,
                  args = params_TOM)
indvTOM_object = list(indvTOM = indvTOM, sessionInfo = sessionInfo(),
                      params = params_TOM, sample_size = sample_size)
saveRDS(indvTOM_object, paste0(outdir, "/1_11/indvTOM_object_20251018.rds"))

## Build consensus networks
params_consNet = list(multiExpr = ResList_Combat, individualTOMInfo = indvTOM_object$indvTOM)

sample_size = lapply(ResList_Combat, function(df) df$data %>% dim()) %>% Reduce(f = "+", x = .) 

consNet = do.call(what = WGCNA::blockwiseConsensusModules,
                  args = c(params_consNet, params))
consNet_object = list(consNet = consNet, sessionInfo = sessionInfo(), 
                      params = c(params_consNet, params), sample_size = sample_size)
saveRDS(consNet_object, paste0(outdir, "/1_11/consNet_object_Combat_20251018.rds"))
#######################################

#######################################
# identify which module in the new network is a mitochondrial module and most similar to the original M3
# first, look at which new module has the greatest overlap with original M3
consNet_object_og <- readRDS(paste0(blwgcna_dir, "/consNet_object_2023-08-15.rds")) # this is created from 2_bl_wgcna.R
OG_M3_genes <- names(consNet_object_og$consNet$colors)[which(consNet_object_og$consNet$colors == 3)] # original M3 genes
new_mt_module_color <- consNet_object$consNet$colors[which(names(consNet_object$consNet$colors) %in% OG_M3_genes)]
table(new_mt_module_color) # all genes in the original M3 are in the new M3
length(which(consNet_object$consNet$colors == 3)) # there are 861 genes in the new M3, the additional gene is ENSG00000006451 (RALA)

# draw a Venn diagram to show the overlapping number
m3_genes <- list(Original = names(consNet_object_og$consNet$colors)[which(consNet_object_og$consNet$colors == 3)], 
                 Combat = names(consNet_object$consNet$colors)[which(consNet_object$consNet$colors == 3)])
pdf(paste0(outdir, "/1_11/overlap_genes_venn.pdf"), width = 3.5, height = 2)
ggVennDiagram::ggVennDiagram(m3_genes, color = "black", lwd = 0.8, lty = 1, label_size = 7/.pt) + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(text = element_text(size = 7),
        legend.key.height = unit(0.1,"in"),
        legend.key.width = unit(0.08, "in")) +
  coord_flip() 
dev.off() # this is Extended Data Fig.9c
#######################################
# We then calculate new M3 association with memory, LMDR, and hippocampal volume
## Calculate ME, separately for MCSA and ADNI
mcsa_me <- moduleEigengenes(comRes_Combat[1:105,], consNet_object$consNet$colors)$eigengenes
adni_me <- moduleEigengenes(comRes_Combat[106:196,], consNet_object$consNet$colors)$eigengenes
consMEs1_Combat <- rbind(mcsa_me, adni_me)
write.csv(consMEs1_Combat %>% rownames_to_column("SubjectID"), paste0(outdir, "/1_11/MEs1_Combat.csv"),row.names = F)

## conduct module-trait association, followed by plotting
comTrait <- readRDS(paste0(indir,"/comTrait.rds")) # a file having cross-sectional endophenotypes for MCSA (from row 1-105) and ADNI (from row 106-196); row names are donorID
me_trait_assoc_df <-  consMEs1_Combat %>%
  rownames_to_column("SubjectID") %>%
  left_join(comTrait %>% rownames_to_column("SubjectID"), by="SubjectID")
## this function can conduct association analysis and plot the correlation dot plot with association estimates and p-values underneath
me_trait_assoc_p_fn <- function(ME, Pheno, module, cov){
  formula <- paste0(ME, " ~ ", Pheno, cov)
  fit_mcsa <- lm(formula, data = me_trait_assoc_df[1:105,])
  mcsa_assoc_p <- ggplot(fit_mcsa$model, aes_string(x = names(fit_mcsa$model)[2], y = names(fit_mcsa$model)[1])) + 
    geom_point(shape=21, color="black", fill="brown1", size=2) +
    stat_smooth(method = "lm", col = "black") +
    labs(caption = paste0("Beta=",formatC(fit_mcsa$coef[[2]], format = "e", digits = 2),"\nP=",formatC(summary(fit_mcsa)$coef[2,4], format = "e", digits = 2)),subtitle = paste0("MCSA: ", str_remove(ME, "E"))) +
    theme_classic() +
    theme(plot.caption = element_text(hjust = 0, size=14),
          plot.subtitle = element_text(hjust=0.5, size=14),
          plot.title = element_blank(),
          legend.position = "none",
          plot.margin = margin(0,0,0,0,unit="mm"),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14)) #trbl
  
  fit_adni <- lm(formula, data = me_trait_assoc_df[106:196,])
  adni_assoc_p <- ggplot(fit_adni$model, aes_string(x = names(fit_adni$model)[2], y = names(fit_adni$model)[1])) + 
    geom_point(shape=21, color="black", fill="brown1", size=2) +
    stat_smooth(method = "lm", col = "black") +
    labs(caption = paste0("Beta=",formatC(fit_adni$coef[[2]], format = "e", digits = 2),"\nP=",formatC(summary(fit_adni)$coef[2,4], format = "e", digits = 2)),subtitle = paste0("ADNI: ", str_remove(ME, "E"))) +
    theme_classic() +
    theme(plot.caption = element_text(hjust = 0, size=14),
          plot.subtitle = element_text(hjust=0.5, size=14),
          plot.title = element_blank(),
          legend.position = "none",
          plot.margin = margin(0,0,3,0,unit="mm"),
          axis.text = element_text(size=14),
          axis.title = element_text(size=14))
  
  assoc_p <- mcsa_assoc_p + adni_assoc_p +  plot_layout(widths = c(6,6)) + plot_annotation(title = ME, theme = theme(plot.title = element_text(hjust=0.5, face="bold", size=14)))
  return(assoc_p)
}
assoc_p_ls <- list()
assoc_p_ls[["Memory"]] <- me_trait_assoc_p_fn(ME="ME3", Pheno="Memory", module = "M3", cov = "+ Sex + Educ + Age")
assoc_p_ls[["LMDR"]] <- me_trait_assoc_p_fn(ME="ME3", Pheno="LMDR", module = "M3", cov = "+ Sex + Educ + Age")
assoc_p_ls[["HippVol"]] <- me_trait_assoc_p_fn(ME="ME3", Pheno="HippVol", module = "M3", cov = "+ Sex + Educ + Age + Hipp_ICV + Mag_HippICV")
pdf(paste0(outdir, "/1_11/me_tr_combat.pdf"), width = 6, height = 10)
assoc_p_ls[["Memory"]]/assoc_p_ls[["LMDR"]]/assoc_p_ls[["HippVol"]]
dev.off() # this is Extended Data Fig.9d
# also save a table of association results
broom::tidy(lm(ME3 ~ Memory + Sex + Educ + Age, data = me_trait_assoc_df[1:105,]), conf.int =T) %>% 
  filter(term == "Memory") %>% 
  mutate(Cohort = "MCSA", .before = "term") %>% 
  bind_rows(
    broom::tidy(lm(ME3 ~ LMDR + Sex + Educ + Age, data = me_trait_assoc_df[1:105, ]), conf.int =T) %>% 
      filter(term == "LMDR") %>% 
      mutate(Cohort = "MCSA", .before = "term")
  ) %>% 
  bind_rows(
    broom::tidy(lm(ME3 ~ HippVol + Sex + Educ + Age + Hipp_ICV + Mag_HippICV, data = me_trait_assoc_df[1:105,]), conf.int =T) %>% 
      filter(term == "HippVol") %>% 
      mutate(Cohort = "MCSA", .before = "term")
  ) %>% 
  bind_rows(
    broom::tidy(lm(ME3 ~ Memory + Sex + Educ + Age, data = me_trait_assoc_df[106:196, ]), conf.int =T) %>% 
      filter(term == "Memory") %>% 
      mutate(Cohort = "ADNI", .before = "term")
  ) %>% 
  bind_rows(
    broom::tidy(lm(ME3 ~ LMDR + Sex + Educ + Age, data = me_trait_assoc_df[106:196, ]), conf.int =T) %>% 
      filter(term == "LMDR") %>% 
      mutate(Cohort = "ADNI", .before = "term")
  ) %>% 
  bind_rows(
    broom::tidy(lm(ME3 ~ HippVol + Sex + Educ + Age + Hipp_ICV + Mag_HippICV, data = me_trait_assoc_df[106:196,]), conf.int =T) %>% 
      filter(term == "HippVol") %>% 
      mutate(Cohort = "ADNI", .before = "term")
  ) %>% 
  write.csv(paste0(outdir, "/1_11/METr_Assoc_Combat.csv"), row.names = F)

# run meta-analysis
METr_Combat <- read.csv(paste0(outdir, "/1_11/METr_Assoc_Combat.csv"))
# prepare files for meta-analysis
meta_df <- inner_join(METr_Combat %>% filter(Cohort == "MCSA"), METr_Combat %>% filter(Cohort == "ADNI"), by = c("term"), suffix = c("_MCSA", "_ADNI")) %>% 
  column_to_rownames("term")
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
  dplyr::select(-c(ends_with("_random"), ends_with("_fixed")))
# save the results
write.csv(meta_final, paste0(outdir, "/1_11/METr_Assoc_meta.csv"), row.names = F) # this is in Supplementary Table 8b

rm(fix_res, rand_res, het_res, meta_df, meta_raw, eff_idx1, eff_idx2, err_idx1, err_idx2)