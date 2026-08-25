# Author: Wei (Adelyn) Tsai; tsai.wei@mayo.edu
# If you are using this script, please cite our study
# This script addresses comment 6 about network robustness from reviewer 1

suppressPackageStartupMessages(library(tidyverse)) #v2.0.0
suppressPackageStartupMessages(library(biomaRt)) #v2.54.1
suppressPackageStartupMessages(library(GenomeInfoDb)) #v1.34.9
suppressPackageStartupMessages(library(dbplyr)) #v2.3.4
suppressPackageStartupMessages(library(CEMiTool)) #v1.22.0
suppressPackageStartupMessages(library(anRichment)) #v1.22
suppressPackageStartupMessages(library(anRichmentMethods)) #v0.96-86
suppressPackageStartupMessages(library(openxlsx)) #v4.2.8
suppressPackageStartupMessages(library(patchwork)) #v1.3.2
suppressPackageStartupMessages(library(gground)) #v1.0.0
suppressPackageStartupMessages(library(ggprism)) #v1.0.7
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
# we first generate network separately for MCSA and ADNI
comRes <- readRDS(paste0(residir, "/comRes.rds")) # from 1_CreateResids.R

# generate mean-variance plot
# according to CEMiTool publication, VST should be done if there's strong linear dependence between the two
# Compute mean and variance per gene
df_mv_mcsa <- data.frame(
  Gene = colnames(comRes[1:105,]),
  Mean = colMeans(comRes[1:105,]),
  Variance = apply(comRes[1:105,], 2, var)
)
df_mv_adni <- data.frame(
  Gene = colnames(comRes[106:196,]),
  Mean = colMeans(comRes[106:196,]),
  Variance = apply(comRes[106:196,], 2, var)
)
# Plot mean-variance relationship
ggplot(df_mv_mcsa, aes(x = Mean, y = Variance)) +
  geom_point(alpha = 0.5, size = 1.2) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Mean–Variance Relationship of Gene Expression Residuals in MCSA",
    x = "Mean Residuals",
    y = "Variance"
  ) +
  geom_smooth(method = "loess", color = "red", se = FALSE, span = 0.5)
ggplot(df_mv_adni, aes(x = Mean, y = Variance)) +
  geom_point(alpha = 0.5, size = 1.2) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Mean–Variance Relationship of Gene Expression Residuals in ADNI",
    x = "Mean Residuals",
    y = "Variance"
  ) +
  geom_smooth(method = "loess", color = "red", se = FALSE, span = 0.5)
# there're none, as expected for residuals

#######################################
# generate MCSA network
mcsa_expr <- as.data.frame(t(comRes[1:105,])) #cemitool requires input data to be data.frame and rows as genes and columns as samples
mcsa_cem <- 
  cemitool(
    expr = mcsa_expr,
    filter = F,
    cor_function = "bicor",
    network_type = "signed"
  ) #about 13mins; 28 modules
generate_report(mcsa_cem, directory = paste0(outdir, "/1_6/Reports"))
write_files(mcsa_cem, directory = paste0(outdir, "/1_6/Tables/MCSA"))
mcsa_cem <- plot_beta_r2(mcsa_cem, plot_title = "Scale independence (beta selection): MCSA")
mcsa_cem <- plot_mean_k(mcsa_cem, title = "Mean connectivity: MCSA")
show_plot(mcsa_cem, "beta_r2")
show_plot(mcsa_cem, 'mean_k')
saveRDS(mcsa_cem, paste0(outdir, "/1_6/mcsa_cem.rds"))


mcsa_cem_spearman <- 
  cemitool(
    expr = mcsa_expr,
    filter = F,
    cor_method = "spearman",
    cor_function = "bicor",
    network_type = "signed"
  ) # this creates same set of MM as using pearson
saveRDS(mcsa_cem_spearman, paste0(outdir, "/1_6/mcsa_cem_spearman.rds"))

# generate ADNI network
adni_expr <- as.data.frame(t(comRes[106:196,])) #cemitool requires input data to be data.frame and rows as genes and columns as samples
adni_cem <- 
  cemitool(
    expr = adni_expr,
    filter = F,
    cor_function = "bicor",
    network_type = "signed"
  ) #about 13mins
generate_report(adni_cem, directory = paste0(outdir, "/1_6/Reports"), force = T)
write_files(adni_cem, directory = paste0(outdir, "/1_6/Tables/ADNI"))
adni_cem <- plot_beta_r2(adni_cem, plot_title = "Scale independence (beta selection): ADNI")
adni_cem <- plot_mean_k(adni_cem, title = "Mean connectivity: ADNI")
show_plot(adni_cem, "beta_r2")
show_plot(adni_cem, 'mean_k')
saveRDS(adni_cem, paste0(outdir, "/1_6/adni_cem.rds"))

#######################################
# save plots for scale-free topology and mean connectivity
pdf(paste0(outdir, "/1_6/beta_meanK_p.pdf"), width = 5, height = 4)
(mcsa_cem@beta_r2_plot$beta_r2_plot + geom_hline(yintercept = 0.8, color = "red", linetype = 2) + theme_classic(base_size = 6) + 
   mcsa_cem@mean_k_plot$mean_k_plot + theme_classic(base_size = 6))/(adni_cem@beta_r2_plot$beta_r2_plot + geom_hline(yintercept = 0.8, color = "red", linetype = 2) + theme_classic(base_size = 6) + 
                                                                       adni_cem@mean_k_plot$mean_k_plot + theme_classic(base_size = 6))
dev.off() # this is in first rebuttal Rebuttal Fig.8
#######################################
# create plots for number of modules and number of genes per module from MCSA and ADNI CEMitool network
mcsa_df4p <- as.data.frame(table(mcsa_cem@module$modules)) %>% # this creates a dataframe where Var1 is module name (M1, M2, etc.) and Freq is the number of genes in each module
  mutate(Var1 = factor(Var1, levels = paste0("M", seq_along(1:28)))) %>% 
  arrange(Var1) %>% 
  mutate(color = WGCNA::labels2colors(c(1:28)), 
         color = factor(color, levels = color)) 
p1 <- mcsa_df4p %>% 
  mutate(Var1 = str_c("CM_", Var1),
         Var1 = factor(Var1, levels = paste0("CM_M", seq_along(1:28)))
         ) %>% 
  ggplot(aes(x = Var1, y = Freq, fill = color)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  geom_text(aes(y = Freq, label = Freq), size = 5/.pt, color = "black", angle = 90, nudge_y = 3, hjust = -0.5, vjust = 0.5) +
  scale_y_continuous(expand = c(0,0), limits = c(-2, 1750)) +
  scale_fill_manual(values = levels(mcsa_df4p$color)) +
  labs(y = "Number of genes", x = "CEMiTool Module", title = "MCSA") +
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust =1))

adni_df4p <- as.data.frame(table(adni_cem@module$modules)) %>% # this creates a dataframe where Var1 is module name (M1, M2, etc.) and Freq is the number of genes in each module
  mutate(Var1 = factor(Var1, levels = paste0("M", seq_along(1:20)))) %>% 
  arrange(Var1) %>% 
  mutate(color = WGCNA::labels2colors(c(1:20)), 
         color = factor(color, levels = color)) 
p2 <- adni_df4p %>% 
  mutate(Var1 = str_c("CM_", Var1),
         Var1 = factor(Var1, levels = paste0("CM_M", seq_along(1:20)))
         ) %>% 
  ggplot(aes(x = Var1, y = Freq, fill = color)) +
  geom_bar(stat = "identity", position = position_dodge(), color = "black") +
  geom_text(aes(y = Freq, label = Freq), size = 5/.pt, color = "black", angle = 90, nudge_y = 3, hjust = -0.5, vjust = 0.5) +
  scale_y_continuous(expand = c(0,0), limits = c(-2, 2500)) +
  scale_fill_manual(values = levels(adni_df4p$color)) +
  labs(y = "Number of genes", x = "CEMiTool Module", title = "ADNI") +
  theme_classic(base_size = 6) +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.x = element_text(angle = 30, hjust =1))

pdf(paste0(outdir, "/1_6/module_p.pdf"), width = 5, height = 7)
p1/p2 + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size=7, face = "bold"))
dev.off() # this is in first rebuttal Rebuttal Fig.9

#######################################

#######################################
# to save module membership information 
mcsaMM <- get_hubs(mcsa_cem, n = "all", method = "kME") %>% 
  bind_rows(.id = "kME") %>% 
  column_to_rownames("kME") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("gene_id") %>% 
  right_join(mcsa_cem@module, by=c("gene_id"="genes")) %>% 
  right_join(MM_tbl %>% 
               dplyr::select(gene_id, gene_name, Module, MM) %>% 
               dplyr::rename(module_og = Module, MM_og = MM) %>% 
               mutate(module_og = paste0("M", module_og))) %>% 
  rowwise() %>%
  mutate(MM = get(modules)) %>%
  ungroup() %>% 
  dplyr::select(gene_id, gene_name, modules, MM, module_og, MM_og, paste0("M", seq_along(1:28)))
adniMM <- get_hubs(adni_cem, n = "all", method = "kME") %>% 
  bind_rows(.id = "kME") %>% 
  column_to_rownames("kME") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("gene_id") %>% 
  right_join(adni_cem@module, by=c("gene_id"="genes")) %>% 
  right_join(MM_tbl %>% 
               dplyr::select(gene_id, gene_name, Module, MM) %>% 
               dplyr::rename(module_og = Module, MM_og = MM) %>% 
               mutate(module_og = paste0("M", module_og))) %>% 
  rowwise() %>%
  mutate(MM = get(modules)) %>%
  ungroup() %>% 
  dplyr::select(gene_id, gene_name, modules, MM, module_og, MM_og, paste0("M", seq_along(1:20)))
write.xlsx(list("MCSA_MM" = mcsaMM, "ADNI_MM" = adniMM), paste0(outdir, "/1_6/MM.xlsx")) # this is in the first rebuttal Rebuttal Table 13

# calculate module eigengene and save results
mcsaME <- mod_summary(mcsa_cem, method = "eigengene")
adniME <- mod_summary(adni_cem, method = "eigengene")
write.xlsx(list("MCSA_ME" = mcsaME, "ADNI_ME" = adniME), paste0(outdir, "/1_6/ME.xlsx"))
#######################################
# next we looked at which modules from MCSA & ADNI Cemi-network are enriched with original M3 genes or mitochondria pathway
MM_tbl <-read.csv(paste0(blwgcna_dir, "/MM_tbl_2023-08-15.csv")) #generated in 2_bl_wgcna.R
m3_genes = subset(MM_tbl, Module == 3)$gene_id 

# enrichment test which MCSA cemi-modules are enriched with M3 genes
m3_mcsa_enrich <- vector('list', length = 28)
names(m3_mcsa_enrich) <- paste0("M", seq_along(1:28))  
for(m in names(m3_mcsa_enrich)) {
  module_genes = subset(mcsa_cem@module, modules == m)$genes 
  module_test_df = data.frame(Module = m, ModuleSize = length(module_genes))
  
  shared_genes = intersect(module_genes, m3_genes)
  contingency = cbind(c(length(shared_genes), (length(module_genes) - length(shared_genes))),
                      c((length(m3_genes) - length(shared_genes)), (9743 - length(module_genes) - (length(m3_genes) - length(shared_genes))))) #9743 is the total number of genes I've in my data
  int_results = data.frame(ModuleM3SharedSize = length(shared_genes),
                           SharedGenes = paste0(subset(MM_tbl, gene_id %in% shared_genes)$gene_name, collapse = ","),
                           Pvalue = fisher.test(contingency, alternative = "greater")$p.value) 
  module_test_df = cbind(module_test_df, int_results)
  
  m3_mcsa_enrich[[m]] = module_test_df
}
m3_mcsa_enrich_df <- bind_rows(m3_mcsa_enrich) %>% mutate(FDR = p.adjust(Pvalue, "fdr"))

# enrichment test which adni cemi-modules are enriched with M3 genes
m3_adni_enrich <- vector('list', length = 20)
names(m3_adni_enrich) <- paste0("M", seq_along(1:20))  
for(m in names(m3_adni_enrich)) {
  module_genes = subset(adni_cem@module, modules == m)$genes 
  module_test_df = data.frame(Module = m, ModuleSize = length(module_genes))
  
  shared_genes = intersect(module_genes, m3_genes)
  contingency = cbind(c(length(shared_genes), (length(module_genes) - length(shared_genes))),
                      c((length(m3_genes) - length(shared_genes)), (9743 - length(module_genes) - (length(m3_genes) - length(shared_genes))))) #9743 is the total number of genes I've in my data
  int_results = data.frame(ModuleM3SharedSize = length(shared_genes),
                           SharedGenes = paste0(subset(MM_tbl, gene_id %in% shared_genes)$gene_name, collapse = ","),
                           Pvalue = fisher.test(contingency, alternative = "greater")$p.value) 
  module_test_df = cbind(module_test_df, int_results)
  
  m3_adni_enrich[[m]] = module_test_df
}
m3_adni_enrich_df <- bind_rows(m3_adni_enrich) %>% mutate(FDR = p.adjust(Pvalue, "fdr"))
# save enrichment results for both cohorts
write.xlsx(list("MCSA" = m3_mcsa_enrich_df, "ADNI" = m3_adni_enrich_df), paste0(outdir, "/1_6/M3Enrich.xlsx")) # this is in the first rebuttal Rebuttal Table 11&12

# enrichment test which MCSA cemi-modules enriched with mitochondrial pathways
Sys.setenv(SQLITE_TMPDIR="/tmp/")
options(stringsAsFactors = FALSE)

go_collection = buildGOcollection(organism = "human")
mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl") 
bm = getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"), bmHeader = TRUE, mart = mart) 
colnames(bm) = c("ensembl_gene_id", "hgnc_symbol", "entrez_gene_id")

entrez_ids = bm[match(dimnames(comRes)[[2]], bm$ensembl_gene_id), ]$entrez_gene_id

mcsa_cem_genes <- mcsa_cem@module$modules
names(mcsa_cem_genes) <- mcsa_cem@module$genes
go_enrichment = enrichmentAnalysis(classLabels = mcsa_cem_genes,
                                   identifiers = entrez_ids,
                                   refCollection = go_collection,
                                   useBackground = "intersection",
                                   maxReportedOverlapGenes=10000,
                                   threshold = 0.05,
                                   thresholdType = "nominal",
                                   getOverlapEntrez = FALSE,
                                   getOverlapSymbols = TRUE)
collectGarbage()
go_mcsa = go_enrichment$enrichmentTable %>% filter(FDR < 0.05)   

# enrichment test which MCSA cemi-modules enriched with mitochondrial pathways
adni_cem_genes <- adni_cem@module$modules
names(adni_cem_genes) <- adni_cem@module$genes
go_enrichment = enrichmentAnalysis(classLabels = adni_cem_genes,
                                   identifiers = entrez_ids,
                                   refCollection = go_collection,
                                   useBackground = "intersection",
                                   maxReportedOverlapGenes=10000,
                                   threshold = 0.05,
                                   thresholdType = "nominal",
                                   getOverlapEntrez = FALSE,
                                   getOverlapSymbols = TRUE)
collectGarbage()
go_adni = go_enrichment$enrichmentTable %>% filter(FDR < 0.05) 
# save GO results for both cohorts
write.xlsx(list("MCSA" = go_mcsa, "ADNI" = go_adni), paste0(outdir, "/1_6/GO.xlsx")) 

# make plots for GO terms
# upload dataframes with GO terms that are cleaned up using REVIGO
# for MCSA, only M2, M3, M13, and M22 are organized; for ADNI, only M3 and M8 are organized
go_mcsa_organized <- read.xlsx(paste0(outdir, "/1_6/GO_organized.xlsx"), sheet = "MCSA_organized")
go_adni_organized <- read.xlsx(paste0(outdir, "/1_6/GO_organized.xlsx"), sheet = "ADNI_organized")

# add module color information 
go_mcsa_organized <- go_mcsa_organized %>% 
  mutate(
    color = case_when(
      class == "M2" ~ WGCNA::labels2colors(2),
      class == "M3" ~ WGCNA::labels2colors(3),
      class == "M13" ~ WGCNA::labels2colors(13),
      class == "M22" ~ WGCNA::labels2colors(22),
    )) %>% 
  mutate(color = factor(color, levels = unique(.$color))) %>% 
  mutate(class = factor(class, levels = c("M22", "M13", "M3", "M2")))
go_adni_organized <- go_adni_organized %>% 
  mutate(
    color = case_when(
      class == "M3" ~ WGCNA::labels2colors(3),
      class == "M8" ~ WGCNA::labels2colors(8),
    )) %>% 
  mutate(color = factor(color, levels = unique(.$color))) %>% 
  mutate(class = factor(class, levels = c("M8", "M3")))

# create a function for GO term plotting
go_p <- function(df, title){
  # left side width
  width <- 0.8
  # x-axis length
  xaxis_max <- max(-log10(df$FDR)) + 1
  # left side category 
  rect.data <- group_by(df, class) %>% #class is module
    top_n(10, wt = -FDR) %>% # got the top 10 terms ranked by FDR
    reframe(n = n()) %>% 
    ungroup() %>% 
    mutate(
      xmin = -3*width,
      xmax = -2*width,
      ymax = cumsum(n),
      ymin = lag(ymax, default = 0) + 0.6,
      ymax = ymax + 0.4
    )
  
  p <- df %>% 
    group_by(class) %>% 
    top_n(10, wt = -FDR) %>% 
    mutate(dataSetName = if_else(grepl("CC", inGroups), paste0(dataSetName, " (CC)"), dataSetName)) %>% # add cellular compartment label if the term belongs to that category
    mutate(module_num = as.numeric(gsub("M", "", class))) %>% 
    mutate(yterm = paste0(class, "_", dataSetName)) %>% # dataSetName is GO term description
    mutate(yterm = factor(yterm, levels = yterm[rev(order(module_num))])) %>% # reorganize GO term order
    ggplot(aes(-log10(FDR), y = yterm, fill = class)) +
    geom_round_col(aes(y = yterm), width = 0.6, alpha = 0.5) +
    geom_text(aes(x = 0.05, label = dataSetName), hjust = 0, size = 7/.pt) +
    # for gene count
    geom_point(aes(x = -width, size = nCommonGenes), 
               stroke = 0, 
               shape = 21, 
               alpha = 0.5) + #nCommonGenes is the number of genes in the dataset that are also related to GO terms
    geom_text(aes(x = -width, label = nCommonGenes), size = 7/.pt) +
    scale_size_continuous(name = "Count") +
    # for module label
    geom_round_rect(
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = class),
      data = rect.data,
      radius = unit(1, "mm"),
      alpha = 0.5,
      inherit.aes = F
    ) +
    geom_text(
      aes(x = (xmin + xmax)/2, y= (ymin + ymax)/2, label = class),
      size = 5/.pt,
      data = rect.data
    ) +
    geom_segment(
      aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
      linewidth = 1.5,
      inherit.aes = F
    ) +
    labs(y = NULL, title = title) +
    scale_fill_manual(name = "Module", values = rev(levels(df$color))) +
    scale_color_manual(values = rev(levels(df$color))) +
    scale_x_continuous(
      breaks = seq(0, xaxis_max, 2),
      expand = expansion(c(0, 0))
    ) +
    theme_prism(base_size = 9) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_text()
    )
  return(p)
}

# save the plot
pdf(paste0(outdir, "/1_6/GO_p.pdf"), width = 7, height = 10)
go_p(go_mcsa_organized, "MCSA")/go_p(go_adni_organized, "ADNI") # this is in the first rebuttal Rebuttal Fig.10
dev.off()

#######################################


#######################################
# conduct module-trait association in MCSA and ADNI using modules from CEMitool
comTrait <- readRDS(paste0(indir,"/comTrait.rds")) # a file having cross-sectional endophenotypes for MCSA (from row 1-105) and ADNI (from row 106-196); row names are donorID


# MCSA
# prepare files for association
mcsa_assoc_df <-  mcsaME %>%
  column_to_rownames("modules") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("SubjectID") %>% 
  left_join(comTrait %>% 
              rownames_to_column("SubjectID") %>% 
              dplyr::select(SubjectID, Sex, Educ, Age, Memory:Mag_HippICV), by="SubjectID") %>% 
  dplyr::select(SubjectID, Sex:Mag_HippICV, paste0("M", seq_along(1:28))) %>% 
  pivot_longer(cols = c(17:ncol(.)), names_to = "Module", values_to = "ME") %>% 
  pivot_longer(-c("SubjectID", "Module","ME","Educ","Age", "Sex","Hipp_ICV", "Mag_HippICV"), names_to = "Pheno", values_to = "Measures")
# run association and bind results
mcsa_assoc_res <- mcsa_assoc_df %>% 
  filter(Pheno != "HippVol") %>% 
  named_group_split(Module, Pheno) %>%
  map(~ lm(ME ~ Measures + Educ + Age + Sex, data = .x)) %>% 
  map_dfr(broom::tidy, conf.int = TRUE, .id = "Module_Pheno") %>%
  filter(term == "Measures") %>%
  separate(Module_Pheno, into=c("Module", "Pheno"), sep = " / ") %>% 
  bind_rows(mcsa_assoc_df %>% 
              filter(Pheno == "HippVol") %>% 
              named_group_split(Module) %>%
              map(~ lm(ME ~ Measures + Educ + Age + Sex + Hipp_ICV, data = .x)) %>% 
              map_dfr(broom::tidy, conf.int = TRUE, .id = "Module") %>%
              filter(term == "Measures")) %>% 
  mutate(Pheno = if_else(is.na(Pheno), "HippVol", Pheno)) %>% 
  group_by(Pheno) %>% 
  mutate(q.value = p.adjust(p.value, "fdr"))

# ADNI
# prepare files for association
adni_assoc_df <-  adniME %>%
  column_to_rownames("modules") %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column("SubjectID") %>% 
  left_join(comTrait %>% 
              rownames_to_column("SubjectID") %>% 
              dplyr::select(SubjectID, Sex, Educ, Age, Memory:Mag_HippICV), by="SubjectID") %>% 
  dplyr::select(SubjectID, Sex:Mag_HippICV, paste0("M", seq_along(1:20))) %>% 
  pivot_longer(cols = c(17:ncol(.)), names_to = "Module", values_to = "ME") %>% 
  pivot_longer(-c("SubjectID", "Module","ME","Educ","Age", "Sex","Hipp_ICV", "Mag_HippICV"), names_to = "Pheno", values_to = "Measures")
# run association and bind results
adni_assoc_res <- adni_assoc_df %>% 
  filter(Pheno != "HippVol") %>% 
  named_group_split(Module, Pheno) %>%
  map(~ lm(ME ~ Measures + Educ + Age + Sex, data = .x)) %>% 
  map_dfr(broom::tidy, conf.int = TRUE, .id = "Module_Pheno") %>%
  filter(term == "Measures") %>%
  separate(Module_Pheno, into=c("Module", "Pheno"), sep = " / ") %>% 
  bind_rows(adni_assoc_df %>% 
              filter(Pheno == "HippVol") %>% 
              named_group_split(Module) %>%
              map(~ lm(ME ~ Measures + Educ + Age + Sex + Hipp_ICV, data = .x)) %>% 
              map_dfr(broom::tidy, conf.int = TRUE, .id = "Module") %>%
              filter(term == "Measures")) %>% 
  mutate(Pheno = if_else(is.na(Pheno), "HippVol", Pheno)) %>% 
  group_by(Pheno) %>% 
  mutate(q.value = p.adjust(p.value, "fdr"))

# we also performed downsampled association
# we randomly selected half of the samples each cohort for 10 times and performed associations 
# we did this for modules from CEMiTool, as well as for WGCNA M3 with memory, LMDR, hippocampal volume

# prepare files for WGCNA M3 association
MEs <- read.csv(paste0(blwgcna_dir, "/MEs1_SP12.csv")) # from 2_bl_wgcna.R

# prepare a list for saving the results in the for loop
mcsa_downRes_ls <- list()
adni_downRes_ls <- list()
wgcnaM3_downRes_ls <- list()
for (i in 1:10){
  seed <- set.seed(i)
  print(paste0("seed: ", i))
  
  downsamples <- c(
  sample(rownames(comTrait)[1:105], size = 53),
  sample(rownames(comTrait)[106:196], size = 46)
  )
  # MCSA - prepare files for downsampled association
  mcsa_down_assoc_df <-  mcsaME %>%
    column_to_rownames("modules") %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("SubjectID") %>% 
    left_join(comTrait %>% 
                rownames_to_column("SubjectID") %>% 
                dplyr::select(SubjectID, Sex, Educ, Age, Memory:Mag_HippICV), by="SubjectID") %>% 
    dplyr::select(SubjectID, Sex:Mag_HippICV, paste0("M", seq_along(1:28))) %>% 
    filter(SubjectID %in% downsamples) %>% 
    pivot_longer(cols = c(17:ncol(.)), names_to = "Module", values_to = "ME") %>% 
    pivot_longer(-c("SubjectID", "Module","ME","Educ","Age", "Sex","Hipp_ICV", "Mag_HippICV"), names_to = "Pheno", values_to = "Measures")
  # run association and bind results for downsampled association
  mcsa_down_assoc_res <- mcsa_down_assoc_df %>% 
    filter(Pheno != "HippVol") %>% 
    named_group_split(Module, Pheno) %>%
    map(~ lm(ME ~ Measures + Educ + Age + Sex, data = .x)) %>% 
    map_dfr(broom::tidy, conf.int = TRUE, .id = "Module_Pheno") %>%
    filter(term == "Measures") %>%
    separate(Module_Pheno, into=c("Module", "Pheno"), sep = " / ") %>% 
    bind_rows(mcsa_down_assoc_df %>% 
                filter(Pheno == "HippVol") %>% 
                named_group_split(Module) %>%
                map(~ lm(ME ~ Measures + Educ + Age + Sex + Hipp_ICV, data = .x)) %>% 
                map_dfr(broom::tidy, conf.int = TRUE, .id = "Module") %>%
                filter(term == "Measures")) %>% 
    mutate(Pheno = if_else(is.na(Pheno), "HippVol", Pheno)) %>% 
    group_by(Pheno) %>% 
    mutate(q.value = p.adjust(p.value, "fdr"))
  # save the results in list
  mcsa_downRes_ls[[i]] <- mcsa_down_assoc_res
  
  
  # ADNI - prepare files for downsampled association
  adni_down_assoc_df <-  adniME %>%
    column_to_rownames("modules") %>% 
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("SubjectID") %>% 
    left_join(comTrait %>% 
                rownames_to_column("SubjectID") %>% 
                dplyr::select(SubjectID, Sex, Educ, Age, Memory:Mag_HippICV), by="SubjectID") %>% 
    dplyr::select(SubjectID, Sex:Mag_HippICV, paste0("M", seq_along(1:20))) %>% 
    filter(SubjectID %in% downsamples) %>% 
    pivot_longer(cols = c(17:ncol(.)), names_to = "Module", values_to = "ME") %>% 
    pivot_longer(-c("SubjectID", "Module","ME","Educ","Age", "Sex","Hipp_ICV", "Mag_HippICV"), names_to = "Pheno", values_to = "Measures")
  # run association and bind results for downsampled association
  adni_down_assoc_res <- adni_down_assoc_df %>% 
    filter(Pheno != "HippVol") %>% 
    named_group_split(Module, Pheno) %>%
    map(~ lm(ME ~ Measures + Educ + Age + Sex, data = .x)) %>% 
    map_dfr(broom::tidy, conf.int = TRUE, .id = "Module_Pheno") %>%
    filter(term == "Measures") %>%
    separate(Module_Pheno, into=c("Module", "Pheno"), sep = " / ") %>% 
    bind_rows(adni_down_assoc_df %>% 
                filter(Pheno == "HippVol") %>% 
                named_group_split(Module) %>%
                map(~ lm(ME ~ Measures + Educ + Age + Sex + Hipp_ICV, data = .x)) %>% 
                map_dfr(broom::tidy, conf.int = TRUE, .id = "Module") %>%
                filter(term == "Measures")) %>% 
    mutate(Pheno = if_else(is.na(Pheno), "HippVol", Pheno)) %>% 
    group_by(Pheno) %>% 
    mutate(q.value = p.adjust(p.value, "fdr"))
  # save the results in list
  adni_downRes_ls[[i]] <- adni_down_assoc_res
  
  # prepare file for wgcna M3 downsampled association
  wgcnaM3_trait_down_assoc_df <-  MEs %>%
    left_join(comTrait %>% rownames_to_column("SubjectID"), by="SubjectID") %>% 
    filter(SubjectID %in% downsamples)
  # conduct association
  wgcnaM3_trait_down_assoc_res <- broom::tidy(lm(ME3 ~ Memory + Sex + Educ + Age, data = wgcnaM3_trait_down_assoc_df[1:53,]), conf.int =T) %>% 
    filter(term == "Memory") %>% 
    mutate(Cohort = "MCSA", .before = "term") %>% 
    bind_rows(
      broom::tidy(lm(ME3 ~ LMDR + Sex + Educ + Age, data = wgcnaM3_trait_down_assoc_df[1:53, ]), conf.int =T) %>% 
        filter(term == "LMDR") %>% 
        mutate(Cohort = "MCSA", .before = "term")
    ) %>% 
    bind_rows(
      broom::tidy(lm(ME3 ~ HippVol + Sex + Educ + Age + Hipp_ICV + Mag_HippICV, data = wgcnaM3_trait_down_assoc_df[1:53,]), conf.int =T) %>% 
        filter(term == "HippVol") %>% 
        mutate(Cohort = "MCSA", .before = "term")
    ) %>% 
    bind_rows(
      broom::tidy(lm(ME3 ~ Memory + Sex + Educ + Age, data = wgcnaM3_trait_down_assoc_df[54:99, ]), conf.int =T) %>% 
        filter(term == "Memory") %>% 
        mutate(Cohort = "ADNI", .before = "term")
    ) %>% 
    bind_rows(
      broom::tidy(lm(ME3 ~ LMDR + Sex + Educ + Age, data = wgcnaM3_trait_down_assoc_df[54:99, ]), conf.int =T) %>% 
        filter(term == "LMDR") %>% 
        mutate(Cohort = "ADNI", .before = "term")
    ) %>% 
    bind_rows(
      broom::tidy(lm(ME3 ~ HippVol + Sex + Educ + Age + Hipp_ICV + Mag_HippICV, data = wgcnaM3_trait_down_assoc_df[54:99,]), conf.int =T) %>% 
        filter(term == "HippVol") %>% 
        mutate(Cohort = "ADNI", .before = "term")
    ) 
  # run meta-analysis
  # prepare files for meta-analysis
  meta_df <- inner_join(wgcnaM3_trait_down_assoc_res %>% filter(Cohort == "MCSA"), 
                        wgcnaM3_trait_down_assoc_res %>% filter(Cohort == "ADNI"), 
                        by = c("term"), suffix = c("_MCSA", "_ADNI")) %>% 
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
  # final meta-analyze dataframe, determining meta-analyzed beta, p-values based on heterogeneity measures I2
  meta_final = meta_df %>%
    mutate(model_meta = ifelse(I2 < 0.25, "Fixed", "Random"),
           beta_meta = ifelse(model_meta =="Fixed", beta_fixed, beta_random),
           se_meta = ifelse(model_meta =="Fixed", se_fixed, se_random),
           p_meta = ifelse(model_meta =="Fixed", p_fixed, p_random),
           CI_H_meta = ifelse(model_meta =="Fixed", upper_fixed, upper_random),
           CI_L_meta = ifelse(model_meta == "Fixed", lower_fixed, lower_random)) %>%
    dplyr::select(-c(ends_with("_random"), ends_with("_fixed")))
  # save the results in list 
  wgcnaM3_downRes_ls[[i]] <- meta_final
  
  rm(i, seed, downsamples, mcsa_down_assoc_df, mcsa_down_assoc_res, adni_down_assoc_df, adni_down_assoc_res, wgcnaM3_trait_down_assoc_df, wgcnaM3_trait_down_assoc_res, 
     fix_res, rand_res, het_res, meta_df, meta_raw, eff_idx1, eff_idx2, err_idx1, err_idx2)
}

mcsa_downRes <- bind_rows(mcsa_downRes_ls, .id = "Seed")
adni_downRes <- bind_rows(adni_downRes_ls, .id = "Seed")
wgcnaM3_downRes <- bind_rows(wgcnaM3_downRes_ls, .id = "Seed")
# save the association results
write.xlsx(list("MCSA" = mcsa_assoc_res, 
                "ADNI" = adni_assoc_res, 
                "MCSA_down" = mcsa_downRes, 
                "ADNI_down" = adni_downRes,
                "WGCNA_M3_down" = wgcnaM3_downRes), 
           paste0(outdir, "/1_6/CEM_Down_Assoc.xlsx")) # these are organized in the first rebuttal Rebuttal Table 13 & 14



# make plots of association results, focusing on mitochondria modules and memory, LMDR and hippocampal volume
p1 <- mcsa_assoc_res %>% 
  filter(Module %in% c("M2", "M3", "M13", "M22") & Pheno %in% c("Memory", "LMDR", "HippVol")) %>% 
  mutate(Phenotype=factor(Pheno, levels=c("HippVol", "LMDR", "Memory")),
         Module = factor(Module, levels = c("M2", "M3", "M13", "M22")),
         lbl = if_else(q.value<0.05, "✱", ""), 
         lbl = if_else(p.value<0.05 & q.value>0.05,"△", lbl)) %>% 
  ggplot(aes(x = Module, y = Phenotype)) +
  geom_point(shape=21, color="black", aes(fill=estimate, size = -log(p.value))) +
  geom_text(aes(label=lbl), fontface = "bold", size=5/.pt) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0)+
  guides(size = guide_legend(order = 1)) +
  labs(title = "MCSA", subtitle = "n = 105", x = "Mitochondrial modules", y = "Endophenotype") +
  theme_bw(base_size = 6) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "black"),
        legend.key.size = unit(0.1,"in"),
        text = element_text(family = "ArielMT"))
p2 <- adni_assoc_res %>% 
  filter(Module %in% c("M3", "M8") & Pheno %in% c("Memory", "LMDR", "HippVol")) %>% 
  mutate(Phenotype=factor(Pheno, levels=c("HippVol", "LMDR", "Memory")),
         Module = factor(Module, levels = c("M3", "M8")),
         lbl = if_else(q.value<0.05, "✱", ""), 
         lbl = if_else(p.value<0.05 & q.value>0.05,"△", lbl)) %>% 
  ggplot(aes(x = Module, y = Phenotype)) +
  geom_point(shape=21, color="black", aes(fill=estimate, size = -log(p.value))) +
  geom_text(aes(label=lbl), fontface = "bold", size=5/.pt) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0)+
  guides(size = guide_legend(order = 1)) +
  labs(title = "ADNI", subtitle = "n = 91", x = "Mitochondrial modules", y = "Endophenotype")+
  theme_bw(base_size = 6) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.key.size = unit(0.1,"in"),
        text = element_text(family = "ArielMT"))
p3 <- mcsa_downRes %>% 
  filter(Module %in% c("M2", "M3", "M13", "M22") & Pheno %in% c("Memory", "LMDR", "HippVol")) %>% 
  mutate(
    Module = factor(Module, levels = c("M2", "M3", "M13", "M22")),
    Pheno = factor(Pheno, levels = c("Memory", "LMDR", "HippVol"))
  ) %>%
  mutate(Seed = factor(Seed, levels = rev(c(1:10)))) %>% 
  mutate(lbl = if_else(q.value<0.05, "✱", ""), 
         lbl = if_else(p.value<0.05 & q.value>0.05,"△", lbl)) %>% 
  ggplot(aes(x = Pheno, y = Seed, fill = estimate)) +
  geom_tile(color = "black") +
  geom_text(aes(label = lbl), color = "black", size = 5/.pt) +
  scale_fill_gradient2(high = "firebrick", low = "dodgerblue", mid = "white") +
  facet_grid(~ Module) +
  labs(title = "MCSA", subtitle = "n = 53", x = "Endophenotype") +
  theme_classic(base_size = 6) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.key.size = unit(0.1,"in"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5),
        text = element_text(family = "ArielMT"))
p4 <- adni_downRes %>% 
  filter(Module %in% c("M3", "M8") & Pheno %in% c("Memory", "LMDR", "HippVol")) %>% 
  mutate(
    Module = factor(Module, levels = c("M3", "M8")),
    Pheno = factor(Pheno, levels = c("Memory", "LMDR", "HippVol"))
  ) %>%
  mutate(Seed = factor(Seed, levels = rev(c(1:10)))) %>% 
  mutate(lbl = if_else(q.value<0.05, "✱", ""), 
         lbl = if_else(p.value<0.05 & q.value>0.05,"△", lbl)) %>% 
  ggplot(aes(x = Pheno, y = Seed, fill = estimate)) +
  geom_tile(color = "black") +
  geom_text(aes(label = lbl), color = "black", size = 5/.pt) +
  scale_fill_gradient2(high = "firebrick", low = "dodgerblue", mid = "white") +
  facet_grid(~ Module) +
  labs(title = "ADNI", subtitle = "n = 46", x = "Endophenotype", caption = "✱: FDR-adjusted q<0.05\n△: p<0.05")+
  theme_classic(base_size = 6) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.key.size = unit(0.1,"in"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5),
        text = element_text(family = "ArielMT"))

design <- 'AAAABB
           CCCCDD'
cairo_pdf(paste0(outdir, "/1_6/assoc.pdf"), width = 6, height = 5)
p1 + p2+ p3 + p4 + plot_layout(design = design)
dev.off() # these are in the first rebuttal Rebuttal Fig.11

# make plot for split-half results for WGCNA M3
cairo_pdf(paste0(outdir, "/1_6/WGCNA_M3DownAssoc.pdf"), width = 4, height = 3)
wgcnaM3_downRes %>% 
  dplyr::select(Seed, Phenotype, estimate_MCSA, p.value_MCSA) %>% 
  dplyr::rename(estimate = estimate_MCSA, p.value = p.value_MCSA) %>% 
  mutate(Cohort = "MCSA") %>% 
  bind_rows(
    wgcnaM3_downRes %>% 
      dplyr::select(Seed, Phenotype, estimate_ADNI, p.value_ADNI) %>% 
      dplyr::rename(estimate = estimate_ADNI, p.value = p.value_ADNI) %>% 
      mutate(Cohort = "ADNI")
  ) %>% 
  bind_rows(
    wgcnaM3_downRes %>% 
      dplyr::select(Seed, Phenotype, beta_meta, p_meta) %>% 
      dplyr::rename(estimate = beta_meta, p.value = p_meta) %>% 
      mutate(Cohort = "Meta-analyzed")
  ) %>% 
  mutate(Seed = factor(Seed, levels = rev(c(1:10)))) %>% 
  mutate(lbl = if_else(p.value<0.05,"△", "")) %>% 
  mutate(Cohort = factor(Cohort, levels = c("MCSA", "ADNI", "Meta-analyzed")),
         Phenotype = factor(Phenotype, levels = c("Memory", "LMDR", "HippVol"))) %>% 
  ggplot(aes(x = Phenotype, y = Seed, fill = estimate)) +
  geom_tile(color = "black") +
  geom_text(aes(label = lbl), color = "black", size = 5/.pt) +
  scale_fill_gradient2(high = "firebrick", low = "dodgerblue", mid = "white") +
  facet_grid(~ Cohort) +
  labs(x = "Endophenotype", caption = "△: p<0.05")+
  theme_classic(base_size = 6.5) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.key.size = unit(0.1,"in"),
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5),
        text = element_text(family = "ArielMT"))
dev.off() # this is in Extended Data Fig.7a

