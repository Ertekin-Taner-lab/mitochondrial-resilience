# Author: Wei (Adelyn) Tsai; tsai.wei@mayo.edu
# If you are using this script, please cite our study
# If you are using this script to generate your own figures, they may not be perfect and may need some tweaking using Adobe Illustrator
# Most of the input files for figures are generated from R scripts 1-5; some are from rebuttal and will be referred.

library(tidyverse) #v2.0.0
library(openxlsx) #v4.2.8
library(rrvgo) #v1.10
library(ggh4x) #v0.3.1
library(ggtext) #v0.1.2
library(patchwork) #v1.3.2
library(cowplot) #v1.2
library(ggrepel) #v0.9.8
library(igraph) #v1.5
library(ggpubr) #v0.6.3
library(sjPlot) #v2.8 
library(scales) #v1.4
library(circlize) #v0.4.16
library(ComplexHeatmap) #v2.24.1
library(png) #v0.1

resids_dir <- "./outdir/resids"
blwgcna_dir <- "./outdir/bl_wgcna"
brwgcna_dir <-"./outdir/br_wgcna"
nduf_assoc_dir <- "./outdir/nduf_assoc_de"
long_dir <- "./outdir/longitudinal"
outdir <-"./outdir/Figs"
rbtdir <- "./Rebuttal/outdir"
rbtdir_v2 <- "./Rebuttal_2nd/outdir"

get_max_contrast_color = function(background_color_vec, cut_off = 0.5) {
  .rgb_to_bw = function(color_chr) {
    rgb_vec = col2rgb(color_chr)
    luminance = 0
    # Some magic values from stackoverflow
    luminance = luminance + rgb_vec["red", 1, drop = TRUE] * 0.299
    luminance = luminance + rgb_vec["green", 1, drop = TRUE] * 0.587
    luminance = luminance + rgb_vec["blue", 1, drop = TRUE] * 0.114
    luminance = luminance / 255
    luminance = as.double(luminance)
    ifelse(luminance > cut_off, "black", "white")
  }
  background_color_vec %>%
    purrr::map_chr(.rgb_to_bw)
}

# Fig.1----
## Fig.1a----
.mapping_module_process = function(module) {
  module %>%
    purrr::map_chr(stringr::str_extract, pattern = "[0-9]+") %>%
    rlang::set_names() %>%
    purrr::map_chr(~ switch(.x,
                            "0" = "Background Module",
                            "1" = "Cell Communication",
                            "2" = "Signaling/Immune",
                            "3" = "Cellular Respiration",
                            "4" = "ncRNA Processing",
                            "5" = "Microtubule Organization",
                            "6" = "Mitochondrion (CC)",
                            "7" = "Transcription",
                            "8" = "Immune",
                            "9" = "RNA Processing",
                            "10" = "H2O2 Catabolism",
                            "11" = "Unknown",
                            "12" = "Cytoplasmic Translation",
                            "13" = "Unknown",
                            "14" = "ncRNA processing",
                            "15" = "Nuclear Export",
                            "16" = "Gene Expression",
                            "17" = "Unknown",
                            "18" = "mRNA stabilization",
                            "19" = "Cellular Component Organization",
                            "20" = "Signal Transduction",
                            "21" = "Nephron Morphogenesis",
                            "22" = "Coagulation",
                            "23" = "Unknown",
                            "24" = "Negative Regulation of Catabolism",
                            "25" = "B Cell Activation",
                            "26" = "mRNA Metabolism",
                            "27" = "Transcription",
                            "28" = "RNA Splicing",
                            "29" = "Enzyme Binding (MF)",
                            "30" = "Mitochondrion (CC)",
                            "31" = "Unknown",
                            "32" = "Viral Defense",
                            "33" = "Unknown", 
                            "34" = "Defense Response",
                            "35" = "Long-chain Fatty Acid Metabolism",
                            "36" = "Unknown")) %>%
    purrr::imap_chr(~ stringr::str_c("M", .y, ": ", .x)) %>%
    rlang::set_names(nm = NULL)
}
consMEs1_SP12 <- read.csv(paste0(blwgcna_dir,"/MEs1_SP12.csv")) %>% column_to_rownames("SubjectID")
consNet_object_SP12 <- readRDS(paste0(blwgcna_dir,"/consNet_object_2023-08-15.rds"))
ct_enrich_df <- read.csv(paste0(blwgcna_dir,"/CT.csv"))

palette <- WGCNA::labels2colors(c(1:36)) 
names(palette) <- paste0("M", seq_along(1:36))

# Module info bar
module_hclust = consMEs1_SP12 %>%
  dplyr::select(-ME0) %>%
  t() %>%
  dist() %>%
  hclust(method="ave")
module_hclust$labels <- str_remove(module_hclust$labels,"E")
module_order = str_remove(module_hclust$labels[module_hclust$order], "E")

module_info <- tibble(module_name = paste0("M",seq_along(1:36))) %>%
  dplyr::mutate(module_name = factor(module_name , levels = module_order)) %>%
  dplyr::mutate(module_to_show = paste0(module_name, "\n(", table(consNet_object_SP12$consNet$colors)[2:37],")")) %>%
  dplyr::mutate(module_to_show = factor(module_to_show, module_to_show)) %>%
  dplyr::mutate(module_color = palette) %>%
  dplyr::mutate(module_color = factor(module_color, module_color)) %>%
  dplyr::mutate(text_color = get_max_contrast_color(module_color)) %>%
  dplyr::mutate(text_color = factor(text_color)) %>%
  dplyr::mutate(module_legend_text = .mapping_module_process(module_name))

module_info_bar <- module_info %>%
  ggplot(aes(x = module_name, y = 1, fill = module_color, label = module_to_show)) +
  geom_tile(color = "black")  +
  ggh4x::scale_x_dendrogram(hclust = module_hclust, position = "top", expand = expansion(mult = 0, add = 0)) +
  scale_fill_manual(values = levels(module_info$module_color)) +
  geom_text(aes(color = text_color), size = 14/.pt, family="ArielMT")  +
  scale_color_manual(values = levels(module_info$text_color)) + 
  labs(y="Module\n(Gene #)") +
  theme_classic() + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0.025,"in"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        axis.text = element_blank(),
        legend.position = "none",
        panel.border = element_blank(),
        plot.margin = margin(0,0,2,0,unit="mm"),
        text = element_text(family = "ArielMT")
  )

# Module-trait correlation
module_trait <- read.csv(paste0(blwgcna_dir,"/module_trait.csv"))

cor_tbl <- module_trait[,c("Module", "Pheno", "beta_meta", "p_meta", "q_meta")] %>% 
  filter(Module != "ME0"& Pheno %in% c("LMDR", "BN", "CF", "AVDEL", "TRAB", "Memory", "Language", "VSP", "Attention", "HippVol")) %>%
  dplyr::mutate(Module=gsub("^ME", "M", Module),
                Module = factor(Module, levels = module_order),
                lbl=if_else(p_meta<0.05, "△",""),
                lbl=if_else(q_meta<0.05, "✱", lbl),
                Pheno = factor(Pheno, levels = c("HippVol", "TRAB", "BN", "CF", "AVDEL", "LMDR", "VSP", "Attention", "Language", "Memory")),
                pheno_level = case_when(
                  Pheno == "HippVol" ~ "HV",
                  Pheno %in% c("TRAB", "BN", "CF", "AVDEL", "LMDR") ~ "Individual",
                  Pheno %in% c("VSP", "Attention", "Language", "Memory") ~ "Composite", 
                  TRUE ~ NA_character_),
                pheno_level = factor(pheno_level, levels = c("Composite", "Individual", "HV")))

mod_trait <- ggplot(cor_tbl, aes(x = Module, y = Pheno, fill = beta_meta)) + 
  geom_tile(color="black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, name=expression(beta))+
  geom_text(aes(label=lbl), size=14.5/.pt, fontface="bold", #position = position_nudge(y = -0.2), 
            family="ArielMT") +
  facet_nested(vars(pheno_level), scales = "free_y", switch = "y", space = "free_y") +
  labs(y="Phenotypes", caption = "△: p<0.05\n✱: FDR-adjusted q<0.05") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.y = element_text(size=16),
    axis.title.y = element_text(size=16),
    axis.text.x = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    strip.text = element_text(size=15),
    strip.background = element_rect(color = "black", linewidth = 0.5),
    plot.caption = element_text(size = 15.5, hjust = 0),
    legend.key.height= unit(0.25, 'in'),
    legend.key.width= unit(0.1, 'in'),
    legend.title = element_text(size=16),
    legend.text = element_text(size=15.5),
    plot.margin = margin(0,0,2,0,unit="mm"),
    text = element_text(family = "ArielMT")
  ) 

# Cell-type enrichment
ct <- ct_enrich_df %>% 
  mutate(lbl=if_else(FDR<0.05, "✱", ""),
         Module = factor(Module, levels = module_order)) %>%
  ggplot(aes(x = Module, y = CT, fill = -log(FDR))) +
  geom_tile(color="black") +
  scale_fill_gradient(low="white", high = "green") + 
  geom_text(aes(label=lbl,), size=15/.pt, fontface="bold", #position = position_nudge(y = -0.2), 
            family="ArielMT") +
  labs(y="Blood Cell Types", caption = "✱: FDR-adjusted q<0.05") +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.y = element_text(size=16),
    axis.title.y = element_text(size=16),
    axis.text.x = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.key.size = unit(0.15,"in"),
    legend.title = element_text(size=16),
    legend.text = element_text(size=15.5),
    plot.caption = element_text(size = 15.5, hjust = 0),
    plot.margin = margin(0,0,2,0,unit="mm"),
    text = element_text(family = "ArielMT")
  )

# GO legend
module_legend = module_info %>%
  mutate(module_name=factor(module_name, levels=paste0("M", 1:36)),
         module_legend_text = factor(module_legend_text, module_legend_text)) 

go_anno <- module_legend %>%
  ggplot(aes(x=module_name, y=1, fill=module_color, label=module_name)) + 
  geom_tile(height=0, width=0) +
  scale_fill_manual(values=levels(module_legend$module_color), labels=levels(module_legend$module_legend_text)) +
  guides(fill = guide_legend(nrow = 6, byrow = FALSE),
         color = guide_none()) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom", 
        line = element_blank(), 
        rect = element_blank(), 
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1.5,"mm"),
        legend.margin = margin(-1, 30, 0, 0),
        legend.spacing.x = unit(0.2,"cm"),
        legend.box.background = ggplot2::element_rect(color = "white"),
        plot.margin = margin(0,0,0,0,unit="mm"),
        text = element_text(family = "ArielMT")
  ) 

fig1a <- module_info_bar + mod_trait + ct + go_anno + plot_layout(heights = c(0.8, 6, 2.8,0), ncol = 1)

## Fig.1b----
comTrait <- readRDS("./indir/comTrait.rds") # a file having cross-sectional endophenotypes for MCSA (from row 1-105) and ADNI (from row 106-196); row names are donorID
me_trait_assoc_df <-  consMEs1_SP12 %>%
  rownames_to_column("SubjectID") %>%
  left_join(comTrait %>% rownames_to_column("SubjectID"), by="SubjectID")
# this function creates the scatter plot and also annotate the plot with beta and nominal p-value from regression
me_trait_assoc_p_fn <- function(ME, Pheno, module, cov){
  formula <- paste0(ME, " ~ ", Pheno, cov)
  fit_mcsa <- lm(formula, data = me_trait_assoc_df[1:105,])
  mcsa_assoc_p <- ggplot(fit_mcsa$model, aes_string(x = names(fit_mcsa$model)[2], y = names(fit_mcsa$model)[1])) + 
    geom_point(shape=21, color="black", fill=palette[module], size=2) +
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
    geom_point(shape=21, color="black", fill=palette[module], size=2) +
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

sig_lmdr <- c("ME2_LMDR","ME3_LMDR","ME30_LMDR") 
sig_mem <- c("ME2_Memory","ME3_Memory","ME30_Memory") 
sig_hv <- c("ME2_HippVol","ME3_HippVol","ME30_HippVol")

lmdr_assoc_p_ls <- list()
mem_assoc_p_ls <- list()
hv_assoc_p_ls <- list()

for (n in sig_lmdr){
  print(n)
  ME <- gsub("_.*","",n)
  Pheno <- gsub(".*_","",n)
  module <- gsub("E","",ME)
  
  cov <- "+ Sex + Educ + Age"
  final_p_ls <- me_trait_assoc_p_fn(ME=ME, Pheno=Pheno, module = module, cov = cov)
  
  lmdr_assoc_p_ls[[n]] <- final_p_ls
}

for (n in sig_mem){
  print(n)
  ME <- gsub("_.*","",n)
  Pheno <- gsub(".*_","",n)
  module <- gsub("E","",ME)
  
  cov <- "+ Sex + Educ + Age"
  final_p_ls <- me_trait_assoc_p_fn(ME=ME, Pheno=Pheno, module = module, cov = cov)
  
  mem_assoc_p_ls[[n]] <- final_p_ls
}

for (n in sig_hv){
  print(n)
  ME <- gsub("_.*","",n)
  Pheno <- gsub(".*_","",n)
  module <- gsub("E","",ME)
  
  cov <- "+ Sex + Educ + Age + Hipp_ICV + Mag_HippICV"
  final_p_ls <- me_trait_assoc_p_fn(ME=ME, Pheno=Pheno, module = module, cov = cov)
  
  hv_assoc_p_ls[[n]] <- final_p_ls
}

rm(sig_lmdr, sig_mem, sig_hv, ME, Pheno, module, cov, n, final_p_ls)

fig1b1 <- wrap_plots(mem_assoc_p_ls, nrow = 2) + plot_annotation(title = "Memory", theme=theme(plot.title = element_text(hjust=0.5, face="bold"),plot.margin = margin(0, 0.2, 0, 0.2, "pt")))
fig1b2 <- wrap_plots(lmdr_assoc_p_ls, nrow = 2)  + plot_annotation(title = "LMDR", theme=theme(plot.title = element_text(hjust=0.5, face="bold"), plot.margin = margin(0, 0.2, 0, 0.2, "pt"))) 
fig1b3 <- wrap_plots(hv_assoc_p_ls, nrow = 2) + plot_annotation(title = "Hippocampal Volume", theme=theme(plot.title = element_text(hjust=0.5, face="bold"),plot.margin = margin(0, 0.2, 0, 0.2, "pt")))

## Final Fig.1----
cairo_pdf(paste0(outdir,"/Fig1a.pdf"), width = 22, height = 11)
wrap_elements(fig1a) 
dev.off()


design = c(
  area(1,1,6,6), #t,l,b,r
  area(1,7,6,12),
  area(7,4,12,9)
)

pdf(paste0(outdir,"/Fig1b.pdf"), width = 15, height = 10)
wrap_elements(fig1b1) + wrap_elements(fig1b2) + wrap_elements(fig1b3) + plot_layout(design = design)
dev.off()

rm(cor_tbl, ct, ct_enrich_df, design, fig1a, fig1b1, fig1b2, fig1b3, go_anno, hv_assoc_p_ls, lmdr_assoc_p_ls, me_trait_assoc_df, mem_assoc_p_ls, mod_trait, module_legend)

# Fig.2 ----
cohort= c("Mayo", "MSSM", "ROSMAP")
region= c("Mayo\nSTG", "CER", "FP", "MSSM\nSTG", "PHG", "IFG", "DLPFC")

mpbb <- readRDS(paste0(brwgcna_dir,"/mp_bb.rds"))

mpbb_df <- mpbb$preservation$Z$ref.Blood$inColumnsAlsoPresentIn.Brain %>%
  rownames_to_column(var = "Module") %>%
  filter(!Module %in% c("0", "0.1")) %>%
  mutate(Module = as.numeric(Module)) %>%
  arrange(Module) %>%
  mutate(
    Module = str_c("M", as.character(Module)),
    Module = factor(Module, levels = paste0("M",seq_along(1:max(consNet_object_SP12$consNet$colors)))),
    moduleSize = as.numeric(table(consNet_object_SP12$consNet$colors)[2:37]),
    mod_fill_color = palette[Module],
    mod_txt_color = get_max_contrast_color(mod_fill_color))


set.seed(42)
fig2a <- mpbb_df %>%
  ggplot(aes(x = moduleSize, y = Zsummary.pres)) +
  geom_point(aes(fill = Module), shape = 21, color = "black", stroke = 1, size = 2) +
  geom_label_repel(aes(label = Module, fill = Module, color = mod_txt_color), size=6/.pt, max.overlaps = Inf, family="ArielMT") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 5, linetype = "dashed", color = "darkgreen") +
  geom_hline(yintercept = 2, linetype = "dashed", color = "blue") +
  scale_color_identity() +
  scale_fill_manual(values = palette) +
  scale_y_continuous(breaks=seq(0,10, by=2.5), limits =  c(-0.4,12)) +
  labs(x= "Module Size",y ="Preservation Z Summary") +
  theme_classic() +
  theme(
    legend.position = "none",
    panel.background =  element_rect(color = "black"),
    panel.border = element_blank(),
    axis.line =  element_blank(),
    axis.text = element_text(size=6),
    axis.title = element_text(size=7),
    text=element_text(family = "ArielMT")) 


load(paste0(brwgcna_dir, "/METr.RData"), verbose = T) 
fig2b <- METr_allNoNA %>%
  mutate(Cohort=case_when(Region %in% c("Mayo\nSTG", "CER")~"Mayo",
                          Region %in% c("FP","MSSM\nSTG","PHG","IFG")~"MSSM",
                          TRUE ~ "ROS\nMAP"),
         term=case_when(term=="DX"~"CTRL DX", 
                        term=="Thal" ~ "Amyloid\n(Thal in Mayo\n& ROSMAP)\n(CERAD in MSSM)",
                        term=="CERAD" ~ "Amyloid\n(Thal in Mayo\n& ROSMAP)\n(CERAD in MSSM)",
                        TRUE ~ term)) %>%
  mutate(Cohort=factor(Cohort, levels=c("Mayo", "MSSM", "ROS\nMAP")),
         Region=factor(Region, levels = c("Mayo\nSTG", "CER", "FP", "MSSM\nSTG", "PHG", "IFG", "DLPFC")),
         term=factor(term, levels = c("Braak", "Amyloid\n(Thal in Mayo\n& ROSMAP)\n(CERAD in MSSM)", "CTRL DX"))) %>%
  mutate(Module=gsub("^ME","Br_M", Module)) %>%
  mutate(Module=factor(Module, levels = c("Br_M1", "Br_M17", "Br_M26"))) %>%
  mutate(lbl=if_else(p.value<0.05, "△", "")) %>% #no one meets the critiera of p<0.05 & q>0.05
  mutate(lbl=if_else(q.value<0.05, "✱", lbl)) %>%
  ggplot(aes(y = term, x = Module, fill = estimate)) +
  geom_tile(color="black", linewidth=0.5) +
  scale_fill_gradient2(low = "blue",high = "red", mid = "white",midpoint = 0, space = "Lab",name="β") +
  geom_text(aes(label=lbl), size=6/.pt) +
  facet_nested(~ Cohort + Region, nest_line = element_line(linetype = 2)) +
  labs(caption = "✱: FDR-adjusted q<0.05", y="Phenotype", x="Brain Module")+
  theme_classic() +
  theme(
    axis.text = element_text(size=6),
    axis.title = element_text(size=7),
    axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    strip.text = element_text(face="bold", size=6),
    strip.background = element_blank(),
    plot.caption = element_text(size = 6, hjust = 0),
    legend.key.size = unit(3,"mm"),
    legend.title = element_text(size=7),
    legend.text = element_text(size=6),
    text=element_text(family = "ArielMT"))


design = c(
  area(1.5,1,4,4),
  area(2,5,3.5,10)
)

cairo_pdf(paste0(outdir,"/Fig2.pdf"), width = 7.08661, height = 4.5)
fig2a + fig2b + plot_layout(design = design)+  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size=7, face = "bold", family = "ArielMT"))
dev.off()

rm(mpbb, mpbb_df, METr_allNoNA, METr_all, fig2a, fig2b, design)

# Please refer to rbt_1_1_new for Fig2.c

# Fig.3 ----
comRes <- readRDS(paste0(resids_dir, "/comRes.rds")) # from 1_CreateResids.R
MM_tbl <- read.csv(paste0(blwgcna_dir,"/MM_tbl_2023-08-15.csv")) # from 2_bl_wgcna.R
nduf_bb <- read.csv(paste0(nduf_assoc_dir,"/nduf_bb.csv")) # from 4_nduf_assoc_de.R
go_df <- read.csv(paste0(blwgcna_dir,"/GOEnrichmentTable_AllONT_2023-08-16.csv")) %>% filter(inGroups == "GO|GO.BP|GO") # from 2_bl_wgcna.R
mapped_genes <- go_df %>% filter(class==3) %>%separate_rows(overlapGenes, sep = "\\|") # genes in Bl_M3

## Fig.3a ----
## we cleaned up the GO terms for visualization using rrvgo
rrvgo_ls <- list()
## clean up if a module has more than 1 term
tmp1 <- go_df %>% group_by(class) %>% mutate(duplicated = n() > 2) # class = module

for (i in unique(subset(tmp1, duplicated==TRUE)$class)){
  tmp2 <- filter(tmp1, class==i)
  simMX <- calculateSimMatrix(tmp2$dataSetID,
                              orgdb="org.Hs.eg.db",
                              ont="BP",
                              method="Rel")
  scores <- setNames(-log10(tmp2$FDR), tmp2$dataSetID) # score terms based on their FDR significance
  
  rt <- reduceSimMatrix(simMX,
                        scores,
                        threshold=0.5,
                        orgdb="org.Hs.eg.db") # rrvgo groups terms into cluster based on their semantic similarity and give score
  
  final <- rt %>% 
    arrange(cluster, desc(score)) %>% 
    group_by(cluster) %>% 
    filter(row_number()==1) # get the highest scored term per group (i.e. cluster)
  
  rrvgo_ls[[paste0("M", i)]] <- final
}
rm(tmp1, tmp2, simMX, scores, rt, final)

go_p_df <- rrvgo_ls$M3 %>% inner_join(go_df %>% filter(class==3), by=c("parent"="dataSetID")) # filter to keep cleaned terms of Bl_M3

lolliplot <- ggplot(go_p_df, aes(x = reorder(term, -log(FDR)), y = -log(FDR))) +
  geom_segment(aes(x=reorder(term, -log(FDR)), xend=reorder(term, -log(FDR)), y=0, yend=-log(FDR))) +
  geom_point(aes(fill=nCommonGenes), color="black",shape=21, size=6) + 
  scale_fill_gradient(low = "white", high="brown", name = "Number\nOverlap\nGenes", breaks=seq(min(go_p_df$nCommonGenes), max(go_p_df$nCommonGenes), by = 25), guide = guide_colorbar(label.hjust = 0.5, label.vjust = -0.02, label.theme = element_text(angle = 90, size=14))) +  
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  theme_classic() +
  coord_flip() +
  labs(x="Pathways", y="-log(FDR)")+ 
  geom_hline(yintercept = -log(0.05), linetype="dotted", color = "red", linewidth=1)+
  theme(panel.grid.major.y = element_blank(),
        text = element_text(family = "ArielMT"),
        axis.text = element_text(size=14),
        axis.title = element_text(size=15),
        legend.key.size = unit(0.25, 'in'),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 15),
        legend.position = "top") 

## fetch genes that are mapped to cleaned-up terms and many times they are mapped to those terms
lbl_genes <- go_p_df %>% dplyr::select(go, term, FDR, overlapGenes) %>% tidyr::separate_rows(overlapGenes, sep="\\|")
lbl_genes_count <- lbl_genes %>% group_by(overlapGenes) %>% dplyr::summarise(n=n())
lbl_genes <- lbl_genes %>% inner_join(lbl_genes_count) %>% inner_join(MM_tbl %>% dplyr::select(gene_name, MM), by=c("overlapGenes"="gene_name")) %>% filter(MM>=0.7)

go_genes <- ggplot(lbl_genes, aes(x=reorder(overlapGenes,-n),y=reorder(term, -log(FDR)))) + 
  geom_tile(aes(fill=n), color="white", linewidth=2) + 
  scale_fill_continuous(low="tan",high="brown", name="Occurrence #", 
                        breaks=seq(min(lbl_genes$n), max(lbl_genes$n)))+
  labs(x="Gene") +
  theme_classic() +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust=0.5, size=14, face = "italic"),
        text = element_text(family = "ArielMT"),
        legend.key.size = unit(0.2,"in"),
        legend.text = element_text(size=14),
        legend.title = element_text(size = 15),
        legend.position = "top")

fig3b <- lolliplot + plot_spacer() + go_genes + plot_layout(widths = c(0.6, 0.02, 1.5)) & theme(plot.margin = margin(0,0,0,0)) # originally Fig3b in initial submission; now is Fig3a

## Fig.3b----
de <- read.csv(paste0(nduf_assoc_dir,"/Bl_Assoc.csv")) # from 4_nduf_assoc_de.R
## organize the dataframe
blood_nduf_deg <- de %>% 
  dplyr::select(Phenotype, gene_name, estimate_MCSA, estimate_ADNI, p.value_MCSA, p.value_ADNI, q.value_MCSA, q.value_ADNI, beta_meta, p_meta, q_meta) %>%
  group_by(Phenotype, gene_name) %>%
  gather(key="Cohort", value="Beta", c("estimate_MCSA", "estimate_ADNI", "beta_meta")) %>%
  gather(key="Cohort2", value="P", c("p.value_MCSA", "p.value_ADNI", "p_meta")) %>%
  gather(key="Cohort3", value="Q", c("q.value_MCSA", "q.value_ADNI", "q_meta")) %>%
  dplyr::select(-Cohort2, -Cohort3) %>%
  dplyr::slice(c(1,14,27)) %>%
  ungroup() %>%
  mutate(Phenotype=factor(Phenotype, levels=c("Memory", "LMDR", "HippVol")),
         Cohort = gsub(".*_","", Cohort),
         Cohort = if_else(Cohort=="meta", "Meta-analyzed", Cohort),
         Cohort = factor(Cohort, levels=c("MCSA", "ADNI","Meta-analyzed")),
         gene_name = factor(gene_name, levels=unique(gene_name)),
         lbl = if_else(Q<0.05, "✱", ""), 
         lbl = if_else(P<0.05 & Q>0.05,"△", lbl),
         xpos = as.integer(gene_name),
         ypos = as.integer(Cohort))

fig3c <- blood_nduf_deg %>%  
  ggplot(aes(x=gene_name, y=Cohort)) +
  geom_point(shape=21, color="black", aes(fill=Beta, size = -log(P))) +
  geom_text(aes(label=lbl, x=xpos, y=ypos), fontface = "bold", size=14/.pt, position = position_nudge(y=0.1)
  ) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, name=expression(beta), 
                       guide = guide_colorbar(title.position = "top", label.position = "bottom", label.hjust = 0.5, label.vjust = -0.02, label.theme = element_text(angle = 90, size=14)))+
  scale_size_continuous(range = c(2,10), 
                        guide = guide_legend(title.position = "top")) + 
  labs(caption = "△:P<0.05\n✱:FDR-adjusted q<0.05") +
  facet_grid(Phenotype~.) +
  theme_classic() +
  theme(axis.text = element_text(size=15),
        axis.text.x = element_text(angle = 45, vjust=1, hjust=1, face = "italic"),
        axis.title= element_blank(),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14.5), 
        legend.key.size = unit(0.2,"in"),
        legend.position = "top",
        plot.caption = element_text(size = 14, hjust = 0),
        strip.background = element_rect(colour = "black", linewidth = 1),
        strip.text = element_text(face="bold", size=15),
        text = element_text(family = "ArielMT"))

design=c(
  area(2,1,7.5,4),
  area(1,4.5,8,12),
  area(9,1,12,12)
)

cairo_pdf(paste0(outdir,"/Fig3bc.pdf"), width = 18, height = 14)
plot_spacer() + wrap_elements(fig3c) + fig3b + plot_layout(design = design)
dev.off() # these are originall Fig3b & c in the initial submission; now is Fig3a&b

rm(rrvgo_ls, go_p_df, lolliplot,lbl_genes, lbl_genes_count, go_genes, de, fig3b, fig3c)

# Fig.3c ----
## Please refer rbt_2_1_3_BF_TM.Rmd for Fig.3c

# Fig.4----
## Fig.4a
br_de <- read.xlsx(paste0(nduf_assoc_dir,"/NDUF_Br_Assoc.xlsx"), sheet = "noNA") # from 4_nduf_assoc_de.R
brain_nduf_deg <- br_de %>% 
  left_join(nduf_bb %>% dplyr::select(gene_id, gene_name, Module_consbrain, MM_consbrain), by="gene_id")
## organize dataframe
fig4a_wgcna <- brain_nduf_deg %>%
  mutate(lbl=paste0("M",Module_consbrain,"\n(",round(MM_consbrain,2),")"),
         header="WGCNA") %>%
  ggplot(aes(x = 1, y = gene_name, fill = MM_consbrain)) +
  geom_tile(color="black", linewidth=0.5) +
  #geom_point(shape=21, color="black", size=8) +
  scale_fill_gradient(low = "white",high = "red",space = "Lab", na.value = "lightgrey") +
  geom_text(aes(label=lbl), size=14/.pt, family = "ArielMT") +
  facet_grid(.~header) +
  theme_classic() +
  theme(
    text = element_text(family = "ArielMT"),
    axis.text.y = element_text(size=14, face="italic"),
    axis.text.x = element_blank(),
    axis.title = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_line(linewidth = 0.5),
    legend.position = "none",
    strip.text = element_text(face="bold", size=14.5),
    strip.background = element_rect(colour = "black"))


fig4a_assoc <- brain_nduf_deg %>%
  mutate(lbl=if_else(q.value<0.05, "✱",""),
         lbl = if_else(p.value<0.05 & q.value>0.05,"△", lbl),
         term=case_when(term=="DX"~"CTRL DX", 
                        term=="Thal" ~ "Amyloid\n(Thal in Mayo & ROSMAP)\n(CERAD in MSSM)",
                        term=="CERAD" ~ "Amyloid\n(Thal in Mayo & ROSMAP)\n(CERAD in MSSM)",
                        TRUE ~ term),
         term=factor(term, levels=c("CTRL DX", "Amyloid\n(Thal in Mayo & ROSMAP)\n(CERAD in MSSM)", "Braak")),
         Region=factor(Region, levels=region),
         Cohort=case_when(Region %in% c("Mayo\nSTG", "CER")~"Mayo",
                          Region %in% c("FP","MSSM\nSTG","PHG","IFG")~"MSSM",
                          TRUE ~ "ROS\nMAP"),
         Cohort=factor(Cohort, levels=c("Mayo", "MSSM", "ROS\nMAP"))) %>%
  ggplot(aes(x = Region, y = gene_name, fill = estimate)) +
  geom_point(shape=21, color="black", aes(size = -log(q.value))) +
  scale_fill_gradient2(low = "blue",high = "red", mid = "white",midpoint = 0, space = "Lab",name="β",guide = guide_colorbar(title.position = "top", label.position = "bottom", label.hjust = 0.5, label.vjust = -0.02, label.theme = element_text(angle = 90, size=14))) +
  scale_size_continuous(range = c(5,20),guide = guide_legend(title.position = "top")) +
  guides(size = guide_legend(order = 1)) +
  geom_text(aes(label=lbl), size=14/.pt) +
  labs(caption = "△:P<0.05\n✱:FDR-adjusted q<0.05", x="Region") + 
  facet_nested(~ term + Cohort, nest_line = element_line(linetype = 2), scales = "free", space = "free") +
  theme_classic() +
  theme(
    text = element_text(family = "ArielMT"),
    axis.text.x = element_text(angle = 45, vjust=1, hjust=1, size=14),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(size=14.5),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_line(linewidth = 0.5),
    legend.key.size = unit(0.25,"in"),
    legend.position = "bottom",
    legend.text = element_text(size=14),
    legend.title = element_text(size=14.5),
    plot.caption = element_text(size = 14, hjust = 0),
    strip.text = element_text(face="bold", size=14.5),
    strip.background = element_blank()
  )

fig4a <- fig4a_wgcna+fig4a_assoc + plot_layout(widths = c(0.15,2))

# Please refer rbt_1_1_new.Rmd for Fig.4b

## Fig.4c
nduf_prot <- read.xlsx("./indir/NDUF_NeuroPro.xlsx") # file that has brain regions & study in which NDUFs are significantly upregulated in controls 
nduf_prot <- nduf_prot %>% 
  group_by(Gene, Locations) %>% 
  mutate(Count = n()) %>%
  ungroup() %>%
  mutate(Region =
           case_when(Locations=="Precuneus" ~ "PRE",
                     Locations=="Cingulate Gyrus" ~ "CG",
                     Locations=="Sensory Cortex" ~ "SC",
                     Locations=="Cerebellum" ~ "CER",
                     Locations=="Frontal Cortex"~"FC",
                     Locations=="Entorhinal Cortex"~"EC",
                     Locations=="Hippocampus"~"HP",
                     Locations=="Parahippocampal Cortex"~"PHC",
                     Locations=="Temporal Cortex"~"TC",
                     Locations=="Parietal Cortex"~"PC",
                     TRUE ~ Locations)) %>%
  mutate(Region = factor(Region, levels=c("FC", "SC", "PC", "PRE", "CG", "HP", "PHC", "TC", "EC", "CER")))

fig4b <- nduf_prot %>% 
  ggplot(aes(x=Region, y=Gene, fill=Count)) +
  geom_point(shape=21, color="black", size = 8) +
  scale_fill_continuous(low="#FFE3D9", high="#FF0E04", name="Number\nof Studies") + 
  geom_text(aes(label=Count), size=14/.pt, family="ArielMT") +
  geom_vline(xintercept = seq(1.5, 9.5, by=1), linetype="dotted") +
  theme_classic() +
  theme(axis.text.y = element_text(face="italic", size=15),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=15),
        axis.title.x= element_text(size=15.5),
        axis.line = element_line(linewidth = 0.5),
        text = element_text(family = "ArielMT")) # this is initially fig4b in the initial submission, now is fig.4c

# S15 ----
mice_deg <- read.csv(paste0(nduf_assoc_dir,"/mice_deg.csv")) # from 4_nduf_assoc_de.R

fig4c <- mice_deg %>% 
  mutate(Region=if_else(Region=="cortex","CTX","HP"),
         Model=factor(Model, levels=c("ADBXD", "5xFAD", "P301S", "APOE-TR")),
         Gene=factor(Gene, levels=unique(Gene)),
         lbl = if_else(q.value < 0.05, "✱",""),
         lbl = if_else(p.value<0.05 & q.value>0.05,"△", lbl),
         ypos = as.integer(Gene)) %>%
  ggplot(aes(x=1, y=Gene)) +
  geom_point(shape=21, color="black", aes(fill=-estimate, size = -log(q.value))) +
  geom_text(aes(label=lbl, x=1, y=ypos), size=14/.pt) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, name="log2FC", guide = guide_colorbar(title.position = "top", label.position = "bottom", label.hjust = 0.5, label.vjust = -0.02, label.theme = element_text(angle = 90, size=14)))+
  scale_size_continuous(
    range = c(3,15),
    guide = guide_legend(title.position = "top")) +
  guides(size = guide_legend(order = 1)) +
  labs(caption = "△:P<0.05\n✱:FDR-adjusted q<0.05") +
  facet_nested(~ Model + Region) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=14, face = "italic"),
        axis.title.y= element_blank(),
        axis.line = element_line(linewidth = 0.5),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14.5), 
        legend.key.size = unit(0.25,"in"),
        legend.position = "bottom",
        plot.caption = element_text(size = 14, hjust = 0),
        strip.text = element_text(face="bold", size=14.5),
        strip.background = element_rect(colour = "black"),
        text = element_text(family = "ArielMT")) # this is initially Fig.4c in the initial submission, now is Extended Data Fig.12

design = c(
  area(1,1,10,8),
  area(1,9,5,11),
  area(6,9,10,12)
)

cairo_pdf(paste0(outdir,"/Fig4.pdf"), width = 25, height = 16)
wrap_elements(fig4a) + fig4b + fig4c + plot_layout(design = design)
dev.off() # this is how Fig4 were generated originally; Illustrator was used to remove mouse DEG results and rearrange Fig4 to what it's now

rm(br_de, fig4a, fig4a_wgcna, fig4a_assoc, fig4b, fig4c)



# Fig.5 ----
longitudinal_res <- read.xlsx(paste0(long_dir, "/longitudinal_new.xlsx"), sheet = "Meta") # from 5_longitudinal.R
# first, create a dotted heatmap
# prepare files
df_p <- longitudinal_res %>% 
  dplyr::select(Phenotype, NDUF, Estimate_MCSA, Estimate_ADNI, `Pr(>|t|)_MCSA`, `Pr(>|t|)_ADNI`, qVal_MCSA, qVal_ADNI, beta_meta, p_meta, q_meta) %>%
  group_by(Phenotype, NDUF) %>%
  gather(key="Cohort", value="Beta", c("Estimate_MCSA", "Estimate_ADNI", "beta_meta")) %>%
  gather(key="Cohort2", value="P", c("Pr(>|t|)_MCSA", "Pr(>|t|)_ADNI", "p_meta")) %>%
  gather(key="Cohort3", value="Q", c("qVal_MCSA", "qVal_ADNI", "q_meta")) %>%
  dplyr::select(-Cohort2, -Cohort3) %>%
  dplyr::slice(c(1,14,27)) %>%
  ungroup() %>%
  mutate(Phenotype=factor(Phenotype, levels=c("Memory", "LMDR", "HippVol")),
         Cohort = gsub(".*_","", Cohort),
         Cohort = if_else(Cohort=="meta", "Meta-analyzed", Cohort),
         Cohort = factor(Cohort, levels=c("MCSA", "ADNI","Meta-analyzed")),
         NDUF = factor(NDUF, levels=unique(NDUF)),
         lbl = if_else(Q<0.05, "✱", ""), 
         lbl = if_else(P<0.05 & Q>0.05,"△", lbl),
         xpos = as.integer(NDUF),
         ypos = as.integer(Cohort))

# create plot
dot_hm <- df_p  %>%  
  ggplot(aes(x=NDUF, y=Cohort)) +
  geom_point(shape=21, color="black", aes(fill=Beta, size = -log(P))) +
  geom_text(aes(label=lbl, x=xpos, y=ypos), fontface = "bold", size=14/.pt, position = position_nudge(y=0.1)
  ) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, name=expression(beta))+
  scale_size_continuous(range = c(2,10)) + 
  labs(caption = "△:P<0.05\n✱:FDR-adjusted q<0.05") +
  facet_grid(Phenotype~.) +
  #geom_hline(yintercept=c(1.5, 2.5), linetype="dotted") +
  theme_classic() +
  theme(axis.text = element_text(size=15),
        axis.text.x = element_text(angle = 45, vjust=1, hjust=1, face = "italic"),
        axis.title= element_blank(),
        legend.text = element_text(size=14),
        legend.title = element_text(size=14.5), 
        legend.key.size = unit(0.2,"in"),
        plot.caption = element_text(size = 14, hjust = 0),
        strip.background = element_rect(colour = "black", linewidth = 1),
        strip.text = element_text(face="bold", size=15),
        text = element_text(family = "ArielMT"))

# second, plot interaction model plot
# function for plotting the model
plot_model_fn <- function(model_type, model_name, cohort_name, rows){
  model <- get(model_type)
  # plot quantiles of NDUFB9 values on the plot
  tertiles <- quantile(model[[model_name]]@frame[["NDUFB9"]], probs = c(0.25, 0.5, 0.75))
  
  plot_model(model[[model_name]],
               type = "pred",
               terms = c("YrDiff", paste0("NDUFB9 [", paste(round(tertiles, 2), collapse = ","), "]"))) +
      scale_colour_brewer(palette = "Set1", aesthetics = c('colour', 'fill'),
                          labels = c("Low (25%)", "Medium (50%)", "High (75%)"), 
                          # Assign red to positive values 
                          # and blue to negative ones
                          direction = -1)  +
      guides(color = guide_legend(reverse = TRUE)) +
      theme_sjplot() +
      ## also annotate the plot wit beta estimates and nominal p values from regression
      labs(caption = paste0("Beta=",formatC(lmerTest:::get_coefmat(model[[model_name]])[rows,1], format = "e", digits = 2),"\nP=",formatC(lmerTest:::get_coefmat(model[[model_name]])[rows,5], format = "e", digits = 2)),subtitle = cohort_name) +
      theme(plot.caption = element_text(hjust = 0, size=14),
            plot.subtitle = element_text(hjust=0.5, size=14),
            plot.title = element_blank(),
            plot.margin = margin(0,0,0,0,unit="mm"),
            axis.text.y = element_text(size=14),
            axis.text.x = element_text(size=14, vjust = 1, hjust = 1),
            legend.title = element_text(face = "italic", size=14),
            axis.title.y = element_text(size=14),
            text = element_text(family = "ArielMT"))
}

p1 <- plot_model_fn("mcsa_models", "Memory_NDUFB9", "MCSA", 7) + labs(x="Years from baseline") + plot_model_fn("adni_models", "Memory_NDUFB9", "ADNI", 7) + labs(x="Years from baseline") 
p2 <- plot_model_fn("mcsa_models", "LMDR_NDUFB9", "MCSA", 7) + labs(x="Years from baseline") + plot_model_fn("adni_models", "LMDR_NDUFB9", "ADNI", 7) + labs(x="Years from baseline") 
p3 <- plot_model_fn("mcsa_models", "HippVol_NDUFB9", "MCSA", 8) + labs(y="Hippocampal Volume", x="Years from baseline") + plot_model_fn("adni_models", "HippVol_NDUFB9", "ADNI", 9) + labs(y="Hippocampal Volume", x="Years from baseline")
# final plots of both MCSA and ADNI for memory and hippocampal volume
p4 <- wrap_plots(list(p1, p2, p3), ncol=3) + plot_layout(guides = "collect")

cairo_pdf(paste0(outdir,"/longitudinal_new_.pdf"), width = 16, height = 10)
dot_hm/p4 
dev.off()

# Fig.6----
## Data prep----
tick = as.raster(readPNG("./Codes/tick_mark2.png"))
## Please refer to supplementary table 18 for scoring
score_circ_df <- data.frame(
  Gene = sort(nduf_bb$gene_name),
  Score = c(50,	29,	37,	16,	8,	41,	34,	4,	12,	17,	38,	39,	31,	45,	33,	25,	25,	22,	25,	59),
  Rank = c(2,	11,	7,	17,	19,	4,	8,	20,	18,	16,	6,	5,	10,	3,	9,	12,	12,	15,	12,	1)
)
## blood circos heatmap dataframe
bl_circ_df <- blood_nduf_deg %>%
  filter(Cohort == "Meta-analyzed") %>%
  left_join(score_circ_df, by=c("gene_name"="Gene")) %>%
  mutate(gene_name= factor(gene_name, levels=unique(gene_name[order(Score)]))) %>%
  dplyr::select(gene_name, Phenotype, Beta) %>%
  pivot_wider(names_from = "Phenotype", values_from = "Beta") %>%
  dplyr::rename(`HV (Cross)` = HippVol, `LMDR (Cross)` = LMDR, `Memory (Cross)`=Memory)
## blood circos heatmap label dataframe
bl_circ_df_lbl <- nduf_bb %>%
  left_join(blood_nduf_deg %>% filter(Cohort == "Meta-analyzed")
  )  %>%
  left_join(score_circ_df, by=c("gene_name"="Gene")) %>%
  mutate(
    gene_name= factor(gene_name, levels=unique(gene_name[order(Score)])),
    Hub_genes_bl=if_else(MM_blood>=0.7, "V", "")) %>%
  dplyr::select(gene_name, Phenotype, Hub_genes_bl, lbl) %>%
  pivot_wider(names_from = "Phenotype", values_from = "lbl") %>%
  dplyr::rename(`HV (Cross)` = HippVol, `LMDR (Cross)` = LMDR, `Memory (Cross)`=Memory)
## blood circos AUC dataframe
bl_auc_circ_df <- data.frame(
  gene_name = sort(nduf_bb$gene_name),
  AUC_count = c(1, 0, 0, 0, 1, 2, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 2)
) %>%
  left_join(score_circ_df, by=c("gene_name"="Gene")) %>%
  mutate(
    gene_name= factor(gene_name, levels=unique(gene_name[order(Score)]))
  ) %>% 
  dplyr::select(-Score, -Rank)

cr_assoc_br_ndufs <- read.xlsx(paste0(rbtdir_v2, "/1_2/ROSMAP_resilience_Association.xlsx"), sheet = "NDUFs_needed") %>% # from /rebuttal_2nd/rbt_1_2.Rmd
  filter(X != "CERAD") %>%
  filter(X2 == term) %>% 
  rename(gene_name = X2) %>% 
  mutate(term = Y) %>% 
  dplyr::select(-X, - Y) %>% 
  left_join(nduf_bb %>% dplyr::select(gene_id, gene_name, MM_consbrain, Module_consbrain), by=c("gene_name")) %>% 
  mutate(Region = "DLPFC") 
## brain circos heatmap dataframe
br_circ_df <- brain_nduf_deg %>% 
  bind_rows(cr_assoc_br_ndufs) %>% 
  group_by(gene_name) %>% 
  summarise(score_dx = sum(term=="DX" & estimate >0 & q.value <0.05)*2 + sum(term=="DX" & estimate >0 & p.value <0.05 & q.value >0.05),
            score_braak = sum(term=="Braak" & estimate <0 & q.value <0.05)*2 + sum(term=="Braak" & estimate <0 & p.value <0.05 & q.value >0.05),
            score_cerad = sum(term=="CERAD" & estimate <0 & q.value <0.05)*2 + sum(term=="CERAD" & estimate <0 & p.value <0.05 & q.value >0.05),
            score_thal = sum(term=="Thal" & estimate <0 & q.value <0.05)*2 + sum(term=="Thal" & estimate <0 & p.value <0.05 & q.value >0.05),
            score_amyloid = sum(score_cerad, score_thal),
            score_ep = sum(grepl("cogn_ep", term) & estimate >0 & q.value <0.05)*2 + sum(grepl("cogn_ep", term) & estimate >0 & p.value <0.05 & q.value >0.05),
            score_ep_slope = sum(grepl("cognep_demog_slope", term) & estimate >0 & q.value <0.05)*2 + sum(grepl("cognep_demog_slope", term) & estimate >0 & p.value <0.05 & q.value >0.05),
            score_resilience = sum(score_ep, score_ep_slope)
  ) %>%
  mutate(across(score_braak:score_amyloid, ~.x *-1)) %>%
  full_join(
    nduf_prot %>% 
      group_by(Gene) %>% 
      distinct(Region, .keep_all = T) %>%
      mutate(Count_Sum_Prot = sum(Count)) %>%
      slice_head() %>%
      dplyr::select(Gene, Count_Sum_Prot), by=c("gene_name"="Gene")
  ) %>%
  left_join(score_circ_df, by=c("gene_name"="Gene")) %>%
  mutate(gene_name= factor(gene_name, levels=unique(gene_name[order(Score)]))) %>%
  left_join(nduf_bb %>% dplyr::select(gene_name, MM_consbrain, Module_consbrain), by=c("gene_name")) %>% 
  dplyr::select(-score_cerad, -score_thal, -score_ep, -score_ep_slope, -Score, -Rank) %>%
  dplyr::rename(
    `CTRL DX`=score_dx,
    Braak=score_braak,
    Amyloid=score_amyloid,
    Resilience = score_resilience,
    `CTRL DX (Prot)`= Count_Sum_Prot
  )
## brain circos heatmap label dataframe
br_circ_df_lbl <- br_circ_df %>%
  mutate_at(.vars=vars(Braak, Amyloid), ~.x*-1) %>%
  mutate_at(.vars=vars(`CTRL DX`:Amyloid, Resilience, `CTRL DX (Prot)`), ~as.character(.x)) %>%
  mutate(Hub_genes_br=if_else(MM_consbrain>=0.7 & Module_consbrain %in% c(1, 17, 26), "V", "")) %>%
  dplyr::select(gene_name, Hub_genes_br, `CTRL DX`:Resilience, `CTRL DX (Prot)`)


longitudinal <- read.xlsx(paste0(rbtdir,"/1_8/longitudinal_new.xlsx"), sheet = "Meta") #from rebuttal/rbt_1_8.R
## blood longitudinal circos heatmap dataframe
bl_circ_long_df <- longitudinal %>%
  dplyr::select(NDUF, Phenotype, beta_meta) %>% 
  rename(gene_name = NDUF, Beta = beta_meta) %>% 
  left_join(score_circ_df, c("gene_name" = "Gene")) %>%
  mutate(gene_name= factor(gene_name, levels=unique(gene_name[order(Score)]))) %>%
  dplyr::select(gene_name, Phenotype, Beta) %>%
  pivot_wider(names_from = "Phenotype", values_from = "Beta") %>%
  dplyr::rename(`HV (Long)` = HippVol, `LMDR (Long)` = LMDR, `Memory (Long)`=Memory)
## blood longitudinal circos labelling dataframe
bl_circ_long_df_lbl <- nduf_bb %>%
  left_join(longitudinal %>% 
              mutate(lbl = if_else(q_meta<0.05, "✱", ""), 
                     lbl = if_else(p_meta<0.05 & q_meta>0.05,"△", lbl)) %>% 
              dplyr::select(NDUF, Phenotype, lbl) %>% 
              rename(gene_name = NDUF), 
            by=c("gene_name")
  )  %>%
  left_join(score_circ_df, by=c("gene_name"="Gene")) %>%
  mutate(
    gene_name= factor(gene_name, levels=unique(gene_name[order(Score)]))) %>%
  dplyr::select(gene_name, Phenotype, lbl) %>%
  pivot_wider(names_from = "Phenotype", values_from = "lbl") %>%
  dplyr::rename(`HV (Long)` = HippVol, `LMDR (Long)` = LMDR, `Memory (Long)`=Memory)

deg_df <- bl_circ_df %>% 
  left_join(bl_circ_long_df, by="gene_name") %>% 
  left_join(br_circ_df) %>% 
  left_join(score_circ_df %>% dplyr::select(Gene, Score), by=c("gene_name"="Gene")) %>%
  mutate(gene_name= factor(gene_name, levels=unique(gene_name[order(Score)]))) 
lbl_df <- bl_circ_df_lbl %>% left_join(bl_circ_long_df_lbl, by="gene_name") %>% left_join(bl_auc_circ_df) %>% left_join(br_circ_df_lbl)
lbl_df$Hub_genes_br[which(lbl_df$gene_name=="NDUFA7")] <- "" # NDUFA7 is missing in AMPAD, so need to do this manually to add the gene to the dataframe


color <- colorRampPalette(c("#EDEEC6", "#F4985A"))(20)
color <- c(color[1:11], rep(color[12],3), color[15:20])
score_circ_df = score_circ_df %>% 
  left_join(data.frame(Color = color, Rank = rev(c(1:11, rep(12,3), 15, 16, 17, 18, 19, 20))), relationship = "many-to-many") %>% 
  distinct(Gene, .keep_all = T) %>%
  mutate(Gene = factor(Gene, levels = Gene[order(Score)]))

## Fig plot ----
cairo_pdf(paste0(outdir,"/Fig6_new_.pdf"), width = 24, height = 24)
#pushViewport(viewport(gp = gpar(fontfamily = "ArialMT")))
circos.clear()
# determines how 'tall' the tracks should be
circos.par(start.degree = 90, canvas.ylim = c(-1.1, 1.1), clock.wise=F, points.overflow.warning=F, gap.after=c(rep(0,19),45)) 

# initializes 
circos.heatmap.initialize(deg_df, split = deg_df$gene_name)

# 1.rank
circos.track(score_circ_df$Gene, ylim=c(0,1), track.height=0.1,
             panel.fun = function(x, y) { 
               circos.axis(labels = F, major.tick = F)
               circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2]+1, niceFacing = T, facing = "clockwise", 
                           labels=CELL_META$sector.index, font=3, adj = c(0.5, 0), cex=2, family="ArialMT")
             })

for(i in 1:nrow(score_circ_df)){
  circos.barplot(
    value = score_circ_df$Score[i]/max(score_circ_df$Score),
    pos = 0.5, 
    col = score_circ_df$Color[i], 
    border = NA, 
    bar_width = 1,
    sector.index = score_circ_df$Gene[i], 
    track.index = 1
  )}

for (i in 1:nrow(score_circ_df)){
  circos.text(0.5, 0.5, sector.index = score_circ_df$Gene[i], labels= as.character(score_circ_df$Score[i]), adj=c(0.5,0.5), facing="bending.inside", niceFacing = T, track.index=1, cex=2, family="ArialMT"
  )
}

# 2.blood_MM
circos.track(lbl_df$gene_name, ylim=c(0,1), track.height= uh(1, "cm"),
             panel.fun = function(x, y) { 
               #circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2]+1.5, niceFacing = T, facing = "clockwise", labels=CELL_META$sector.index, cex = 1, ps = 5, font=3, adj = c(0.5, 0))
               circos.axis(labels = F, major.tick = F)
             })

for (i in 1:nrow(lbl_df)){
  if (lbl_df$Hub_genes_bl[i]=="V"){
    circos.raster(tick, CELL_META$xcenter, CELL_META$ycenter, width = "0.7cm", height="0.7cm",
                  sector.index = lbl_df$gene_name[i], track.index = 2,
                  facing = "inside")
  }}

# 3.blood DEG
col_fun1 <- colorRamp2(c(-0.01, 0, 0.125), c("#7C50FF", "white", "firebrick2"))
deg_df %>%
  dplyr::select("Memory (Cross)", "LMDR (Cross)", "HV (Cross)", "Memory (Long)", "LMDR (Long)", "HV (Long)") %>%
  circlize::circos.heatmap(
    track.height = 0.2,
    bg.border = "black", cell.border = "black",
    col = col_fun1#, show.sector.labels = T
  )
for (i in 1:nrow(lbl_df)){
  circos.text(0.5, 1:6-0.5, sector.index = lbl_df$gene_name[i], 
              labels= c(lbl_df$`HV (Long)`[i], lbl_df$`LMDR (Long)`[i], lbl_df$`Memory (Long)`[i],lbl_df$`HV (Cross)`[i], lbl_df$`LMDR (Cross)`[i], lbl_df$`Memory (Cross)`[i]), 
              track.index=3, adj=c(0.5,0.5), facing="bending.inside", niceFacing = T, cex=2, family="ArialMT")
}

# 4.blood AUC
circos.track(lbl_df$gene_name, ylim=c(0,1), track.height= uh(1, "cm"),
             panel.fun = function(x, y) { 
               #circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2]+1.5, niceFacing = T, facing = "clockwise", labels=CELL_META$sector.index, cex = 1, ps = 5, font=3, adj = c(0.5, 0))
               circos.axis(labels = F, major.tick = F)
             })

for (i in 1:nrow(lbl_df)){
  circos.text(0.5, 0.5, sector.index = lbl_df$gene_name[i], 
              labels= as.character(lbl_df$AUC_count[i]), adj=c(0.5,0.5), 
              facing="bending.inside", niceFacing = T, track.index=4, cex=2, family="ArialMT"
  )}

# 5.brain MM
circos.track(lbl_df$gene_name, ylim=c(0,1), track.height=uh(1,"cm"),
             panel.fun = function(x, y) { 
               circos.axis(labels = F, major.tick = F)
             })

for (i in 1:nrow(lbl_df)){
  if (lbl_df$Hub_genes_br[i]=="V"){
    circos.raster(tick, CELL_META$xcenter, CELL_META$ycenter, width = "0.7cm", height="0.7cm",
                  sector.index = lbl_df$gene_name[i], track.index = 5,
                  facing = "inside")
  }}

# 6.brain DEG
col_fun2 <- colorRamp2(c(-12, 0, 12), c("#7C50FF", "white", "firebrick2"))
deg_df %>%
  dplyr::select(`CTRL DX`:`CTRL DX (Prot)`) %>%
  circlize::circos.heatmap(
    track.height = 0.15,
    bg.border = "black", cell.border = "black",
    col = col_fun2
  )
for (i in 1:nrow(br_circ_df_lbl)){
  circos.text(0.5, 1:5-0.5, sector.index = br_circ_df_lbl$gene_name[i], 
              labels= c(br_circ_df_lbl$`CTRL DX (Prot)`[i],
                        br_circ_df_lbl$Resilience[i],
                        br_circ_df_lbl$Amyloid[i],
                        br_circ_df_lbl$Braak[i],
                        br_circ_df_lbl$`CTRL DX`[i]), 
              track.index=6, adj=c(0.5,0.5), facing="bending.inside", niceFacing = T, cex=2, family="ArialMT")
}



# group labels
circos.track(track.index = 4, panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 1) { # the last sector
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), 0,
                CELL_META$cell.xlim[2] + convert_x(10, "mm"), 8.1,
                col = "darkred", border = "black")
    circos.text(CELL_META$cell.xlim[2] + convert_x(4.2, "mm"), 4,
                "Blood", facing = "clockwise", cex=2, family="ArialMT")
    
    circos.rect(CELL_META$cell.xlim[2] + convert_x(11, "mm"), 0,
                CELL_META$cell.xlim[2] + convert_x(19, "mm"), 6.5,
                col = "rosybrown1", border = "black")
    circos.text(CELL_META$cell.xlim[2] + convert_x(14.6, "mm"), 3,
                "Association", facing = "clockwise", cex=1.5, family="ArialMT")
  }
}, bg.border = NA)


circos.track(track.index = 6, panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 1) { # the last sector
    circos.rect(CELL_META$cell.xlim[2] + convert_x(1, "mm"), 0,
                CELL_META$cell.xlim[2] + convert_x(10, "mm"), 7,
                col = "plum1", border = "black")
    circos.text(CELL_META$cell.xlim[2] + convert_x(4.2, "mm"), 2.8,
                "Brain", facing = "clockwise", cex=2, family="ArialMT")
    
    circos.rect(CELL_META$cell.xlim[2] + convert_x(11, "mm"), 0,
                CELL_META$cell.xlim[2] + convert_x(19, "mm"), 5,
                col = "rosybrown1", border = "black")
    circos.text(CELL_META$cell.xlim[2] + convert_x(14.6, "mm"), 2.65,
                "Association", facing = "clockwise", cex=1.5, family="ArialMT")
  }
}, bg.border = NA)


# track y-labels
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 1) { # the last sector
    circos.text(CELL_META$cell.xlim[2] + convert_x(1, "mm"), 
                0.5, "Score", 
                adj = c(0, 0.5), facing = "inside", track.index=1, cex=1.5, family="ArialMT")
    circos.text(CELL_META$cell.xlim[2] + convert_x(8, "mm"), 
                0.5, "WGCNA M3 Hub?", 
                adj = c(0, 0.5), facing = "inside", track.index=2, cex=1.5)
    circos.text(rep(CELL_META$cell.xlim[2], 6) + convert_x(14, "mm"), 
                1:6 - 0.5, c("HV (Longitudinal)","LMDR (Longitudinal)","Memory (Longitudinal)", "HV (Cross-sectional)","LMDR (Cross-sectional)","Memory (Cross-sectional)"), 
                adj = c(0, 0.5), facing = "inside", track.index=3, cex=1.5, family="ArialMT")
    circos.text(CELL_META$cell.xlim[2] + convert_x(14, "mm"), 
                0.5, "Resilience status prediction", 
                adj = c(0, 0.5), facing = "inside", track.index=4, cex=1.5, family="ArialMT")
    circos.text(CELL_META$cell.xlim[2] + convert_x(10, "mm"), 
                0.5, "WGCNA M1/M17/M26 Hub?",
                adj = c(0, 0.5), facing = "inside", track.index=5, cex=1.5, family="ArialMT")
    circos.text(rep(CELL_META$cell.xlim[2], 5) + convert_x(19, "mm"), 
                1:5 - 0.5, c("CTRL DX (Protein)", "Resilience", "Amyloid","Braak","CTRL DX"),
                adj = c(0, 0.5), facing = "inside", track.index=6, cex=1.5, family="ArialMT")
      }
}, bg.border = NA)


# legend
lgd_deg <- Legend(at = c("Negative", "Positive", "No available data"),
                  grid_height = unit(12, "mm"),
                  grid_width = unit(8, "mm"),
                  legend_gp = gpar(fill=c("#7C50FF", "firebrick2", "gray")),
                  title_gp = gpar(fontsize = 20, fontface = "bold", fontfamily="ArialMT"),
                  labels_gp = gpar(fontsize = 20, fontfamily="ArialMT"),
                  title_position = "topleft", 
                  border="black",
                  title = "DEG", nrow = 3)

col_fun4 <- colorRamp2(sort(score_circ_df$Score), color)

lgd_score <- Legend(at = c(10, 20, 30, 40, 50,60), col_fun = col_fun4,
                    grid_height = unit(16.5, "mm"),
                    grid_width = unit(8, "mm"),
                    title_gp = gpar(fontsize = 20, fontface = "bold", fontfamily="ArialMT"),
                    labels_gp = gpar(fontsize = 20, fontfamily="ArialMT"),
                    title_position = "topleft", 
                    title = "Score", nrow = 1)
lgd_sig <- Legend(at = c("✱: FDR-adjusted q<0.05\n△: p<0.05"), grid_height = unit(12, "mm"),
                  grid_width = unit(0, "mm"),
                  labels_gp = gpar(fontsize = 20, fontfamily="ArialMT"),
                  title_position = "topleft", 
                  title = "", nrow = 1)

lgd_list <- packLegend(lgd_score, lgd_deg, lgd_sig)

draw(lgd_list, x = unit(25, "mm"), y = unit(25, "mm"), just = c("left", "bottom"))
#popViewport()
dev.off()

# S2----
ave1 <- comTrait[1:105,] %>%
  pivot_longer(-c(1:5,15:27,29:33), names_to = "Pheno", values_to = "Measures") %>% group_by(Pheno) %>% summarise(mean=mean(Measures, na.rm=T)) %>%
  mutate(Pheno=factor(Pheno, levels=c("Memory", "Language", "Attention", "VSP","LMDR", "AVDEL", "CF", "BN", "TRAB",  "HippVol"))) 

ave2 <- comTrait[106:196,] %>%
  mutate(TRAB=-TRAB) %>%
  pivot_longer(-c(1:5,15:27,29:33), names_to = "Pheno", values_to = "Measures") %>% group_by(Pheno) %>% summarise(mean=mean(Measures, na.rm=T)) %>%
  mutate(Pheno=factor(Pheno, levels=c("Memory", "Language", "Attention", "VSP","LMDR", "AVDEL", "CF", "BN", "TRAB",  "HippVol"))) 

s2a <- comTrait[1:105,] %>% 
  pivot_longer(-c(1:5,15:27,29:33), names_to = "Pheno", values_to = "Measures") %>%
  mutate(Pheno=factor(Pheno, levels=c("Memory", "Language", "Attention", "VSP","LMDR", "AVDEL", "CF", "BN", "TRAB",  "HippVol"))) %>%
  ggplot(aes(x=Measures)) +   
  geom_histogram(bins = 10, fill="lightblue", color="black") +
  facet_wrap(vars(Pheno), scales="free", nrow = 2) +
  geom_vline(data=ave1, aes(xintercept = mean),        # Add line for mean
             col = "red",
             lwd = 0.5) +
  labs(caption = "Red lines indicate mean", title="MCSA") +
  theme_classic() + 
  theme(text = element_text(size=7),
        strip.text = element_text(size=6.5),
        plot.caption = element_text(hjust=0, size=6.5),
        plot.title = element_text(hjust=0.5, size=7))

s2b <- comTrait[106:196,] %>% 
  mutate(TRAB=-TRAB) %>%
  pivot_longer(-c(1:5,15:27,29:33), names_to = "Pheno", values_to = "Measures") %>%
  mutate(Pheno=factor(Pheno, levels=c("Memory", "Language", "Attention", "VSP","LMDR", "AVDEL", "CF", "BN", "TRAB",  "HippVol"))) %>%
  ggplot(aes(x=Measures)) +   
  geom_histogram(bins = 10, fill="lightgreen", color="black") +
  facet_wrap(vars(Pheno), scales="free", nrow = 2) +
  geom_vline(data=ave2, aes(xintercept = mean),        # Add line for mean
             col = "red",
             lwd = 0.5) +
  labs(caption = "Red lines indicate mean", title="ADNI") +
  theme_classic() + 
  theme(text = element_text(size=7),
        strip.text = element_text(size=6.5),
        plot.caption = element_text(hjust=0, size=6.5),
        plot.title = element_text(hjust=0.5, size=7))

pdf(paste0(outdir,"/S2.pdf"), width = 7, height = 8.8)
s2a / s2b + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size=7, face="bold"))
dev.off()

rm(ave1, ave2, s2a, s2b)

# S3----
mp_blwgcna <- readRDS(paste0(blwgcna_dir, "/mp_SP12_nPerm100_2023-12-06.rds")) # from 2_bl_wgcna.R

mp_df_blwgcna <- list(
  `MCSA` = mp_blwgcna$preservation$Z$ref.consNet$inColumnsAlsoPresentIn.Discovery,
  `ADNI` = mp_blwgcna$preservation$Z$ref.consNet$inColumnsAlsoPresentIn.Replication
) %>%
  purrr::map(tibble::rownames_to_column, var = "Module") %>%
  bind_rows(.id = "Cohort") %>%
  filter(!Module %in% c("0", "0.1")) %>%
  mutate(Module = as.numeric(Module)) %>%
  arrange(desc(Cohort), Module) %>%
  mutate(
    Cohort = factor(Cohort, levels=c("MCSA", "ADNI")),
    Module = str_c("M", as.character(Module)),
    Module = factor(Module, levels = paste0("M",seq_along(1:36))),
    mod_fill_color = palette[Module],
    mod_txt_color = get_max_contrast_color(mod_fill_color))


s3 <- mp_df_blwgcna %>%
  mutate(Cohort=factor(Cohort, levels=c("MCSA", "ADNI")),
         moduleSize = rep(as.numeric(table(consNet_object_SP12$consNet$colors)[2:37]),2)) %>%
  ggplot(aes(x = moduleSize, y = Zsummary.pres)) +
  geom_point(aes(fill = Module), shape = 21, color = "black", stroke = 1, size = 4) +
  coord_cartesian(clip = "off") +
  geom_label_repel(aes(label = Module, fill = Module, color = mod_txt_color), max.overlaps = Inf, size=6/.pt, xlim = c(-Inf, Inf), ylim = c(-Inf, Inf)) +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 5, linetype = "dashed", color = "darkgreen") +
  geom_hline(yintercept = 2, linetype = "dashed", color = "blue") +
  scale_color_identity() +
  scale_fill_manual(values = palette) +
  facet_grid(cols =vars(Cohort)) +
  labs(
    x= "Module Size",
    y ="Preservation Z Summary"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 6.5),
    strip.background = element_rect(color = "black", fill = "white", linewidth = 1),
    panel.background = element_rect(color = "black"),
    axis.line = element_blank(),
    text = element_text(size=7)
  ) 

pdf(paste0(outdir,"/S3.pdf"), width = 5.6, height = 6)
s3
dev.off()

rm(mp_df, s3)

# S4 ----
# Please refer to rebuttal/rbt_1_12.R

# S5----
consNet_object_br <- readRDS(paste0(brwgcna_dir, "/consNet_SP12.rds"))
mpbr <- readRDS(paste0(brwgcna_dir,"/mp_br.rds"))
palette_br <- WGCNA::labels2colors(c(1:42))
names(palette_br) <- paste0("M", seq_along(1:42))

mpbr_df <- list(
  `Mayo STG` = mpbr$preservation$Z$ref.consNet$inColumnsAlsoPresentIn.MayoSTG,
  `CER` = mpbr$preservation$Z$ref.consNet$inColumnsAlsoPresentIn.CER,
  `FP` = mpbr$preservation$Z$ref.consNet$inColumnsAlsoPresentIn.FP,
  `IFG` = mpbr$preservation$Z$ref.consNet$inColumnsAlsoPresentIn.IFG,
  `PHG` = mpbr$preservation$Z$ref.consNet$inColumnsAlsoPresentIn.PHG,
  `MSSM STG` = mpbr$preservation$Z$ref.consNet$inColumnsAlsoPresentIn.MSSMSTG,
  `DLPFC` = mpbr$preservation$Z$ref.consNet$inColumnsAlsoPresentIn.DLPFC
) %>%
  purrr::map(tibble::rownames_to_column, var = "Module") %>%
  bind_rows(.id = "Region") %>%
  filter(!Module %in% c("0", "0.1")) %>%
  mutate(Module = as.numeric(Module),
         Cohort=case_when(Region %in% c("CER", "Mayo STG") ~ "Mayo",
                          Region %in% c("FP", "IFG", "PHG", "MSSM STG") ~ "MSSM",
                          TRUE ~ "ROSMAP")) %>%
  mutate(
    Cohort = factor(Cohort, levels=cohort),
    Region = factor(Region, levels=c("Mayo STG", "CER", "FP", "MSSM STG", "PHG", "IFG", "DLPFC")),
    Module = str_c("M", as.character(Module)),
    Module = factor(Module, levels = paste0("M",seq_along(1:42))),
    mod_fill_color = palette_br[Module],
    mod_txt_color = get_max_contrast_color(mod_fill_color),
    moduleSize = rep(as.numeric(table(consNet_object_br$consNet$colors)[2:43]),7)) %>%
  arrange(Module)

s4 <- ggplot(mpbr_df, aes(x = Zsummary.pres, y = Module, fill = Region)) +
  geom_bar(stat = "identity", position = position_dodge()) + 
  scale_y_discrete(limits=rev) +
  geom_vline(xintercept = 10, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 5, linetype = "dashed", color = "darkgreen") +
  geom_vline(xintercept = 2, linetype = "dashed", color = "blue") +
  labs(
    x ="Preservation Z Summary"
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 7),
    legend.key.size = unit(0.2, "in")
  ) 

pdf(paste0(outdir,"/S4.pdf"), width = 5, height = 7.08)
s4
dev.off() # originally extended data fig.4 in the initial submission now is extended data fig.5

rm(mpbr_df, s4)

# S6 ----
# Module GO terms
.mapping_module_process = function(module) {
  module %>%
    purrr::map_chr(stringr::str_extract, pattern = "[0-9]+") %>%
    rlang::set_names() %>%
    purrr::map_chr(~ switch(.x,
                            "0" = "Background Module",
                            "1" = "Cellular Respiration",
                            "2" = "Cell Migration",
                            "3" = "Telomerase RNA localization",
                            "4" = "Synaptic Signaling",
                            "5" = "ncRNA Processing",
                            "6" = "RNA Metabolism",
                            "7" = "Unknown",
                            "8" = "Immune",
                            "9" = "Protein Folding",
                            "10" = "Small Molecule Metabolism",
                            "11" = "Synaptic Signaling",
                            "12" = "Developmental Process",
                            "13" = "Ca2+ Channel Activity",
                            "14" = "Nucleic Acid Metabolism",
                            "15" = "Unknown",
                            "16" = "Immune",
                            "17" = "ncRNA Processing/Mitochondrial Translation",
                            "18" = "Endosome (CC)",
                            "19" = "RNA Splicing",
                            "20" = "MAPK Cascade",
                            "21" = "Synaptic Signaling",
                            "22" = "Cytoplasmic Translation",
                            "23" = "Myelination",
                            "24" = "Ras Signaling",
                            "25" = "mRNA Metabolism",
                            "26" = "Cellular Respiration",
                            "27" = "Unknown",
                            "28" = "Nucleic Acid Metabolism",
                            "29" = "K+ Transmembrane Transport",
                            "30" = "Macromolecule Metabolism",
                            "31" = "Biosynthesis Process",
                            "32" = "Synaptic Membrane Potential",
                            "33" = "Immune", 
                            "34" = "ABC-type Transporter Activity (MF)",
                            "35" = "Nucleic Acid Metabolism",
                            "36" = "Synaptic Signaling",
                            "37" = "Immune",
                            "38" = "Ca2+ Transport",
                            "39" = "GTPase Signaling",
                            "40" = "Protein Localization to Cell Surface",
                            "41" = "Hepatocyte Proliferation",
                            "42" = "Postsynapse (CC)")) %>%
    purrr::imap_chr(~ stringr::str_c("M", .y, ": ", .x)) %>%
    rlang::set_names(nm = NULL)
}


consMEs1_br <- read.csv(paste0(brwgcna_dir,"/MEs1_SP12_Br.csv")) %>% column_to_rownames("SubjectID") # from 3_br_wgcna.R
ct_enrich_df_br <- read.csv(paste0(brwgcna_dir,"/CT_42modules.csv")) # from 3_br_wgcna.R

# Module info bar
module_hclust_br = consMEs1_br %>%
  dplyr::select(-ME0) %>%
  t() %>%
  dist() %>%
  hclust(method="ave")
module_hclust_br$labels <- str_remove(module_hclust_br$labels,"E")
module_order_br = str_remove(module_hclust_br$labels[module_hclust_br$order], "E")

module_info_br <- tibble(module_name = paste0("M",seq_along(1:42))) %>%
  dplyr::mutate(module_name = factor(module_name , levels = module_order_br)) %>%
  dplyr::mutate(module_to_show = paste0(module_name, " (", table(consNet_object_br$consNet$colors)[2:43],")")) %>%
  dplyr::mutate(module_to_show = factor(module_to_show, module_to_show)) %>%
  dplyr::mutate(module_color = palette_br) %>%
  dplyr::mutate(module_color = factor(module_color, module_color)) %>%
  dplyr::mutate(text_color = get_max_contrast_color(module_color)) %>%
  dplyr::mutate(text_color = factor(text_color)) %>%
  dplyr::mutate(module_legend_text = .mapping_module_process(module_name))

module_info_bar_br <- module_info_br %>%
  ggplot(aes(x = 1, y = module_name, fill = module_color, label = module_to_show)) +
  geom_tile(color = "black")  +
  ggh4x::scale_y_dendrogram(hclust = module_hclust_br, position = "left", expand = expansion(mult = 0, add = 0)) +
  scale_fill_manual(values = levels(module_info_br$module_color)) +
  geom_text(aes(color = text_color), size = 5/.pt, family="ArielMT")  +
  scale_color_manual(values = levels(module_info_br$text_color)) + 
  labs(x="Module\n(Gene #)") +
  theme_classic() + 
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0.025,"in"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=5.5),
        axis.text = element_blank(),
        legend.position = "none",
        panel.border = element_blank(),
        plot.margin = margin(0,0,2,0,unit="mm"),
        text = element_text(family = "ArielMT"))

# Cell-type enrichment
ct_br <- ct_enrich_df_br %>% 
  mutate(lbl=if_else(FDR<0.05, "✱", ""),
         Module = factor(Module, levels = module_order_br),
         CT=factor(CT, levels=c("Basophil", "Bcell", "CD4T","CD8T", "DC", "Monocyte", "NK", "Neutrophil", "Ast", "End", "Mic", "Neu", "Oli"))) %>%
  ggplot(aes(x = CT, y = Module, fill = -log(FDR))) +
  geom_tile(color="black") +
  scale_fill_gradient(low="white", high = "green", guide = guide_colorbar(title.position = "top", position = "bottom", label.position = "bottom", label.hjust = 0.5, label.vjust = -0.02, label.theme = element_text(angle = 90, size=5))) + 
  geom_text(aes(label=lbl,), size=6/.pt, family="ArielMT") +
  labs(x="Cell Types", caption = "∗: FDR-adjusted q<0.05") +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.text.x = element_text(size=5, angle = 45, vjust = 1, hjust=1),
    axis.title.x = element_text(size=6),
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.key.size = unit(0.12,"in"),
    legend.title = element_text(size=6),
    legend.text = element_text(size=5.5),
    plot.caption = element_text(size = 5, hjust = 1),
    plot.margin = margin(0,0,2,0,unit="mm"),
    text = element_text(family = "ArielMT"))

# GO legend
module_legend_br = module_info_br %>%
  mutate(module_name=factor(module_name, levels=paste0("M", 1:42)),
         module_legend_text = factor(module_legend_text, module_legend_text)) 

go_anno_br <- module_legend_br %>%
  ggplot(aes(x=1, y=module_name, fill=module_color, label=module_name)) + 
  geom_tile(height=0, width=0) +
  scale_fill_manual(values=levels(module_legend_br$module_color), labels=levels(module_legend_br$module_legend_text)) +
  guides(fill = guide_legend(nrow = 21, byrow = FALSE),
         color = guide_none()) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        legend.position = "bottom", 
        line = element_blank(), 
        rect = element_blank(), 
        legend.title = element_blank(),
        legend.text = element_text(size = 5),
        legend.key.size = unit(1.5,"mm"),
        legend.margin = margin(-1, 30, 0, 0),
        legend.box.background = ggplot2::element_rect(color = "white"),
        legend.spacing.x = unit(0.5, 'cm'),
        plot.margin = margin(0,0,0,0,unit="mm"),
        text = element_text(family = "ArielMT")) 

s5 <- module_info_bar_br + ct_br + plot_layout(widths = c(0.8, 4), nrow = 1, ncol = 2)


cairo_pdf(paste0(outdir,"/S5_modules.pdf"), width = 3, height = 6)
wrap_elements(s5)
dev.off()

cairo_pdf(paste0(outdir,"/S5_GO.pdf"), width = 4, height = 8.26)
go_anno_br
dev.off() # originally extended data fig.4 in the initial submission now is extended data fig.6

# S7 ----
# Please refer rebuttal/rbt_1_6.R for Extended Data Fig.7a
# Please refer rebuttal/rbt_2_1.Rmd for Extended Data Fig.7b-d
# Please refer rebuttal/rbt_1_13.R for Extended Data Fig.7e

# S8 ----
# Please refer rebuttal/rbt_2_3.Rmd for Extended Data Fig.8

# S9 ----
# Please refer rebuttal_2nd/rbt_1_5.Rmd for Extended Data Fig.9a
# Please refer rebuttal/rbt_1_11.R for Extended Data Fig.9b-d

# S10 ----
# Please refer rebuttal_2nd/rbt_1_3.Rmd for Extended Data Fig.10

# S11 ----
# Please refer rebuttal_2nd/rbt_1_2.Rmd for Extended Data Fig.11

# S12----
## S12a----
# adapt from hdWGCNA (https://smorabit.github.io/hdWGCNA/index.html)
# blood M3
## first calculate the topology matrix from expression data; all parameters are set the same as WGCNA networks generation parameters
TOM_blood = WGCNA::TOMsimilarityFromExpr(comRes, power = 12, corType = "bicor", networkType = "signed", TOMType = "signed", maxPOutliers = 0.1, nThreads = 6, verbose = 3)
## focus on Bl_M3 genes that are hub genes and mapped to at least one GO term
M3hub2plot <- subset(MM_tbl, Module ==3 & MM >=0.7 & gene_name %in% mapped_genes$overlapGenes)[,c("gene_id", "gene_name", "MM")] %>% arrange(desc(MM))
inModule = is.finite(match(names(consNet_object_SP12$consNet$colors), M3hub2plot$gene_id)) ## need to do this so we can determine which genes among all WGCNA included genes should be plotted 
## while preserving the correct index of genes, so when indexing TOM_blood, it will be correct
M3_TOM_blood = TOM_blood[inModule, inModule]

gA <- graph_from_adjacency_matrix(M3_TOM_blood[1:15, 1:15], mode = "undirected", weighted = TRUE, diag = FALSE)
gB <- graph_from_adjacency_matrix(M3_TOM_blood[16:38, 16:38], mode = "undirected", weighted = TRUE, diag = FALSE)
layoutCircle1 <- rbind(layout.circle(gA)/2, layout.circle(gB))
g1 <- graph_from_adjacency_matrix(M3_TOM_blood, mode = "undirected", weighted = TRUE, diag = FALSE)
V(g1)$name <- M3hub2plot$gene_name[1:38]
V(g1)$color <- case_when(V(g1)$name %in% nduf_bb$gene_name ~ "red",
                         TRUE~ "black")

cairo_pdf(paste0(outdir,"/Fig3a.pdf"), width = 6, height =6, pointsize = 15) # originally Fig3a in the initial submission; now is Extended Fig.10a
par(mar=c(0,0,1,0))
plot(g1, edge.color = adjustcolor(palette[3], alpha.f = 0.25), edge.alpha = 0.25, edge.width=0.2, vertex.color = palette[3], vertex.label = V(g1)$name, vertex.label.dist = 1, vertex.label.degree = -pi/4, vertex.label.color = V(g1)$color, vertex.label.family="ArielMT", vertex.label.font = 3, vertex.frame.color = "black", layout = jitter(layoutCircle1), vertex.size = 10)
title("Blood M3: Top 38 hub genes", line=0.2, family="ArielMT")
dev.off()

rm(M3hub2plot, inModule, M3_TOM_blood, gA, gB)

## S12b ----
# Please refer rebuttal/rbt_2_1_3_BF_TM.Rmd for Extended Data Fig.12b

## S12d&e ----
# Please refer rebuttal_2nd/rbt_1_6.Rmd for Extended Data Fig.12d&e

# S13 & S14 & S16 ----
# Please refer rebuttal/rbt_1_4_cbx_LM10.Rmd for Extended Data Fig.13a, Extended Data Fig.14b&c and Extended Data Fig.16
# Please refer rebuttal_2nd/rbt_1_4.Rmd for Extended Data Fig.13b and Extended Data Fig.14a

# S17 ----
cqn_rawraw <-  read_delim("./indir/PaX108_R01resilience_gene_CQN.txt", delim = "\t") #CQN before filtering based on -3 threshold, 43704 genes
cqn_rawraw <- cqn_rawraw[,c("GeneId", "Chromosome", "Start", "End", "Length", "GeneName", "GeneBiotype", rownames(comRes)[1:105])]


s6a <- cqn_rawraw %>% 
  dplyr::select(GeneId, GeneName, 8:112) %>%
  pivot_longer(cols=3:ncol(.), names_to = "SubjectID", values_to = "Expression") %>%
  ggplot(aes(x=Expression, group=SubjectID, color=SubjectID)) + 
  geom_density() + 
  geom_vline(aes(xintercept = -3, color="red"), linetype="dashed")+
  theme_classic() + labs(title=paste0("MCSA CQN RNAseq Expression Distribution"), y='Density') +
  theme(plot.title=element_text(hjust=0, size=7),text = element_text(size=7), legend.position = "none")


## distribution after CQN normalization and distribution of microarray
cqn_raw <- read_delim("./indir/PaX108_R01resilience_gene_CQN_neg3_postQC_s105.txt", delim = "\t") # MCSA cqn file for 18046 genes. Each row is a gene and each column is a donor

adni_mc <- readRDS("./ADNI_s91_microarray_withAnno.rds") #microarray file for 91 donors from ADNI. This has 10116 protein-coding genes, gene id as column names, donorID as row names

mcsa_expr <- cqn_raw[cqn_raw$GeneId %in% colnames(comRes), rownames(comRes)[1:105]]
adni_expr <- adni_mc[adni_mc$EnsgID_NCBI %in% colnames(comRes), rownames(comRes)[106:196]]

s6b1 <- mcsa_expr %>%
  pivot_longer(cols=1:105, names_to = "SubjectID", values_to = "Expression") %>%
  ggplot(aes(x=Expression, group=SubjectID, color=SubjectID)) + 
  geom_density() + 
  theme_classic() + labs(title="MCSA", y="Density") +
  theme(plot.title=element_text(hjust=0.5, size=7), plot.subtitle = element_text(hjust=0.5, size=5.5), text = element_text(size=7), legend.position = "none")

s6b2 <- adni_expr %>% 
  pivot_longer(cols=1:91, names_to = "SubjectID", values_to = "Expression") %>%
  ggplot(aes(x=Expression, group=SubjectID, color=SubjectID)) + 
  geom_density() + 
  theme_classic() + labs(title="ADNI", y="Density") +
  theme(plot.title=element_text(hjust=0.5, size=7),text = element_text(size=7), legend.position = "none")

s6b <- s6b1+s6b2


pdf(paste0(outdir,"/S6.pdf"), width = 7, height = 4)  
s6a + wrap_elements(s6b) + plot_layout(widths = c(0.8,1.5)) + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size=7, face = "bold"))
dev.off() # this is originally extended data fig.6 in the initial submission, now is extended data fig.17

rm(s6a, s6b, s6b1, s6b2, cqn_raw, cqn_rawraw, adni_mc, mcsa_expr, adni_expr)

# S18 ----
setLabels = c("Mayo STG", "Mayo CER", "MSSM FP","MSSM STG","MSSM PHG","MSSM IFG","ROSMAP DLPFC")
br_cqn_thresholds=c(0.5, -2, 0.5, 0.5, 0.5, 2.5, 1.5)
fileLabels = c("Mayo_STG", "Mayo_CER", "MSSM_FP","MSSM_STG","MSSM_PHG","MSSM_IFG","ROSMAP_DLPFC")

br_cqn_dfs <- readRDS("./indir/ampad_cqn_dfs_57602genes_postQC-Subj.rds") # a list of AMP-AD cqn before filtering based on threshold; rows are genes columns are donors
br_cqn_dfs_long <- br_cqn_dfs %>% purrr::map(., ~ pivot_longer(.x, -gene_id, names_to = "SampleID", values_to = "Expression"))

for (i in 1:7) {
  # Plot expression densities per sample
  p = br_cqn_dfs_long[[i]] %>% 
    ggplot(aes(x = Expression)) +
    geom_density(aes(group = SampleID , color = SampleID)) + 
    geom_vline(xintercept = br_cqn_thresholds[i] , color="red", linetype="dashed", linewidth=0.5) +
    labs(title = paste0("Dataset: ",setLabels[i]) , x = "Expression", y = "Density") +
    scale_x_continuous(breaks = seq(-12, 12, 1)) +
    theme_classic() +
    theme(legend.position = "none", 
          text = element_text(size=24),
          axis.text.x = element_text(angle=90, hjust = 1, vjust = 1))
  
  
  pdf(paste0(outdir,"/S7a_",fileLabels[i],".pdf"), width = 6.5, height = 6)  
  plot(p)
  dev.off()
}
rm(i, p)


br_cqn_filtered <- readRDS("./indir/br_cqn2use.rds") # a list containing cqn of brain datasets after filtering based on thresholds. 
# In each list is a dataframe where rows are genes and columns are donors
br_cqn_long_filtered <- br_cqn_filtered %>% purrr::map(., ~pivot_longer(.x, cols=4:ncol(.x), names_to = "SampleID", values_to = "CQN")) %>% bind_rows()

s7b <- br_cqn_long_filtered %>%
  mutate(Cohort=factor(Cohort, levels=c("Mayo", "MSSM", "ROSMAP")),
         Region=factor(Region, levels=c("Mayo\nSTG", "CER", "FP","MSSM\nSTG","PHG","IFG","DLPFC"))) %>%
  ggplot(., aes(x=CQN, group=SampleID, color=SampleID)) + 
  geom_density() + ylim(0, 0.35) +
  facet_nested(~ Cohort + Region, nest_line = element_line(linetype = 2)) +
  labs(y="Density") +
  theme_classic() + 
  theme(legend.position = "none",
        strip.background = element_rect(color = "black", fill = "white", linetype = "solid", linewidth = 1))

pdf(paste0(outdir,"/S7b.pdf"), width = 7, height = 4)  
s7b
dev.off() # originally extended data fig.7 in the initial submission, now is extended data fig.18

rm(setLabels, br_cqn_thresholds, fileLabels, br_cqn_dfs, br_cqn_dfs_long, br_cqn_filtered, br_cqn_long_filtered, s7b)

# S19 ----
powers_df_bl <- read.csv(paste0(blwgcna_dir,"/powers_df_2023-08-15.csv")) %>% mutate(set = factor(set, levels=c("MCSA", "ADNI"))) # from 2_bl_wgcna.R
s8a1 <- powers_df_bl %>% 
  ggplot(aes(x=Power, y= SFT.R.sq)) +
  geom_text(aes(label=Power, color=set), size=6/.pt) +
  scale_y_continuous(breaks = seq(0, 1, by=0.2), limits=c(0,1)) +
  labs(title="Scale Free Topology Model Fit", x="Scale Free Topology Model Fit", y="Soft Threshold Power") +
  theme_classic() +
  theme(text=element_text(size=7),
        plot.title = element_text(hjust=0.5, size=7),
        panel.grid.major.x = element_blank() , #set vertical line blank
        panel.grid.major.y = element_line(linewidth=.1, color="black"), #explicitly set horizontal line
        panel.background = element_rect(color = "black", linewidth = 0.5))

s8a2 <- powers_df_bl %>% 
  ggplot(aes(x=Power, y= mean.k.)) +
  geom_text(aes(label=Power, color=set), size=6/.pt) +
  labs(title="Mean Connectivity", x="Soft Threshold Power", y="Mean Connectivity") +
  theme_classic() +
  theme(text=element_text(size=7),
        plot.title = element_text(hjust=0.5, size=7),
        panel.grid.major.x = element_blank() , #set vertical line blank
        panel.grid.major.y = element_line(linewidth=.1, color="black"), #explicitly set horizontal line
        panel.background = element_rect(color = "black", linewidth = 0.5))

s8a <- s8a1 + s8a2 + plot_layout(guides = "collect") & theme(legend.key.size = unit(0.1, 'in'))

powers_df_br <- read.csv(paste0(brwgcna_dir,"/powers_df.csv")) %>% mutate(set = factor(set, levels=region)) # from 3_br_wgcna.R


s8b1 <- powers_df_br %>% 
  ggplot(aes(x=Power, y= SFT.R.sq)) +
  geom_text(aes(label=Power, color=set), size=6/.pt) +
  scale_y_continuous(breaks = seq(0, 1, by=0.2), limits=c(0,1)) +
  labs(title="Scale Free Topology Model Fit", x="Scale Free Topology Model Fit", y="Soft Threshold Power") +
  theme_classic() +
  theme(text=element_text(size=7),
        plot.title = element_text(hjust=0.5,size=7),
        panel.grid.major.x = element_blank() , #set vertical line blank
        panel.grid.major.y = element_line(linewidth=.1, color="black"), #explicitly set horizontal line
        panel.background = element_rect(color = "black", linewidth = 0.5))

s8b2 <- powers_df_br %>% 
  ggplot(aes(x=Power, y= mean.k.)) +
  geom_text(aes(label=Power, color=set), size=6/.pt) +
  scale_y_continuous(breaks = seq(0, 7000, by=1000), limits=c(0,7000)) +
  labs(title="Mean Connectivity", x="Mean Connectivity", y="Soft Threshold Power") +
  theme_classic() +
  theme(text=element_text(size=7),
        plot.title = element_text(hjust=0.5,size=7),
        panel.grid.major.x = element_blank() , #set vertical line blank
        panel.grid.major.y = element_line(linewidth=.1, color="black"), #explicitly set horizontal line
        panel.background = element_rect(color = "black", linewidth = 0.5),
        legend.key.size = unit(0.1, 'in'))  

s8b <- s8b1 + s8b2 + plot_layout(guides = "collect")  & theme(legend.key.size = unit(0.1, 'in'))

pdf(paste0(outdir,"/S8.pdf"), width = 7, height = 7)
wrap_elements(s8a)/wrap_elements(s8b) + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size=7, face = "bold"))
dev.off() # originally extended data fig.8 in the initial submission, now is extended data fig.19

# S20 ----
load(paste0(blwgcna_dir,"/AdjacencyMx_MCSA_ADNI.RData"), verbose=T) # from 2_bl_wgcna.R

## randomly select 100 genes and check their adjacency. Do this 10 times.
set.seed(123)
mcsa_ls <- lapply(1:10, function(x) adj_mcsa_l[sample(nrow(adj_mcsa_l), size = 100),])
mcsa_ls <- mcsa_ls %>% map(~ .x %>% mutate(Cohort = "MCSA"))

adni_ls <- list()
for (i in 1:10){
  adni_ls[[i]] <- adj_adni_l %>% filter(Pair %in% mcsa_ls[[i]]$Pair)
}
adni_ls <- adni_ls %>% map(~ .x %>% mutate(Cohort = "ADNI"))

adj_ls <- map2(mcsa_ls, adni_ls, ~ rbind(.x, .y)) %>% map(~.x %>% mutate(Cohort=factor(Cohort, levels=c("MCSA", "ADNI"))))


adj_hist <- list()
for (i in 1:10){
  tmp <- ggplot(adj_ls[[i]], aes(x=log(Corr), fill=Cohort)) + 
    geom_histogram(alpha=0.2, position="identity") +
    labs(y="Count", x="log(Adjacency)") + 
    theme_classic() +
    theme(text = element_text(size=7),
          legend.key.height = unit(0.2, 'in'),
          legend.key.width = unit(0.1, 'in'),
          axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1))
  
  adj_hist[[i]] <- tmp
  rm(tmp)
}
# adjacency for only Bl_M3 genes
m3_pairs <- expand.grid(MM_tbl$gene_id[which(MM_tbl$Module==3)], MM_tbl$gene_id[which(MM_tbl$Module==3)])
m3_pairs$pair = paste0(m3_pairs$Var1, "_", m3_pairs$Var2)
m3_adj <- adj_mcsa_l %>% 
  filter(Pair %in% m3_pairs$pair) %>% 
  mutate(Cohort = "MCSA") %>% 
  bind_rows(
    adj_adni_l %>%
      filter(Pair %in% m3_pairs$pair) %>%
      mutate(Cohort="ADNI")
  ) %>%
  mutate(Cohort=factor(Cohort, levels=c("MCSA", "ADNI")))

s9a <- wrap_plots(adj_hist, nrow = 2, ncol = 5) + plot_layout(guides = "collect") + plot_annotation(title = "Random 100 Genes Adjacency", theme = theme(plot.title = element_text(size=7, hjust=0.5)))
s9b <- ggplot(m3_adj, aes(x=log(Corr), fill=Cohort)) + 
  geom_histogram(alpha=0.2, position="identity", bins = 50) +
  labs(y="Count", x="log(Adjacency)", title = "Blood M3 Genes Adjacency") + 
  theme_classic() +
  theme(text = element_text(size=7),
        legend.key.height = unit(0.2, 'in'),
        legend.key.width = unit(0.1, 'in'),
        plot.title = element_text(hjust=0.5))

design = c(
  area(1,1,4,6),
  area(5,3,7,5)
)

pdf(paste0(outdir,"/S9.pdf"), width = 6.5, height = 6)
wrap_elements(s9a) + s9b + plot_layout(design = design) + plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size=7, face = "bold"))
dev.off() # originally extended data fig.9 in the initial submission, now is extended data fig.20.