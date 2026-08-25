# Author: Wei (Adelyn) Tsai; tsai.wei@mayo.edu
# If you are using this script, please cite our study
# This script addresses comment 7 about competitive gene testing from reviewer 1

suppressPackageStartupMessages(library(tidyverse)) #v2.0.0
suppressPackageStartupMessages(library(fgsea)) #v1.34
suppressPackageStartupMessages(library(openxlsx)) #v4.2.8

# list of input/output directories
indir <- "./indir"
blwgcna_dir <- "./outdir/bl_wgcna"

outdir <- "./Rebuttal/outdir"

# load in mitochondrial pathway associated genes 
mitopaths <- read_delim(paste0(indir, "/MitoPathways3.0.gmx")) #downloaded from MitoCarta 3.0 (https://www.broadinstitute.org/mitocarta/mitocarta30-inventory-mammalian-mitochondrial-proteins-and-pathways)
mitopaths <- mitopaths[-1,]

# load in module membership information 
MM_tbl <- read.csv(paste0(blwgcna_dir, "/MM_tbl_2023-08-15.csv")) #from 2_bl_wgcna.R
# sor the M3 genes based on their module membership
M3_sort <- MM_tbl %>% filter(Module == 3) %>% arrange(desc(MM))
ranks <- M3_sort$MM
names(ranks) <- M3_sort$gene_name

# Focus on all broad pathways in mitochondria ----
# Broad institute organize MitoPathways into following bigger categories: (see C MitoPathways in https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Human.MitoCarta3.0.xls)
# mtDNA maintenance
# mtRNA metabolism
# Translation
# Protein import, sorting and homeostasis
# OXPHOS
# Metabolism
# Small molecule transport
# Signaling
# Mitochondrial dynamics and surveillance

# turn the dataframe into a list of mitochondrial genes, each list correspond to a pathway
mt_ls_allpaths <- lapply(mitopaths, function(x) {
  x[!is.na(x)]
})
mt_ls_testedpaths <- mt_ls_allpaths[c("mtDNA_maintenance", "mtRNA_metabolism", "Translation", 
                                      "Protein_import_sorting_and_homeostasis", "OXPHOS", "Metabolism", "Small_molecule_transport", 
                                      "Signaling", "Mitochondrial_dynamics_and_surveillance")]
# look at number of overlapping between mt_ls_testedpaths and genes in M3 genes
for (n in names(mt_ls_testedpaths)){
  print(paste0("Number of overlap of ", n, ": ", length(intersect(subset(MM_tbl, Module == 3)$gene_name, mt_ls_testedpaths[[n]]))))
}
# "Number of overlap of mtDNA_maintenance: 5"
# "Number of overlap of mtRNA_metabolism: 13"
# "Number of overlap of Translation: 40"
# "Number of overlap of Protein_import_sorting_and_homeostasis: 16"
# "Number of overlap of OXPHOS: 45"
# "Number of overlap of Metabolism: 65"
# "Number of overlap of Small_molecule_transport: 9"
# "Number of overlap of Signaling: 5"
# "Number of overlap of Mitochondrial_dynamics_and_surveillance: 11"

set.seed(123) # need to set.seed right before running to get the same results
fgseaRes_broadpaths <- fgsea(pathways = mt_ls_testedpaths, 
                             stats    = ranks,
                             scoreType = "pos",
                             minSize  = 5,
                             maxSize  = 500)
fgseaRes_broadpaths <- as.data.frame(fgseaRes_broadpaths) # this shows Translation is also enriched, but not as much as OXPHOS


# Focus on OXPHOS complexes & translation subcategories ----
# OXPHOS is subcategorized into five different complexes. 
# Translation is subcategorized into: Mitochondrial ribosome, Mitochondrial ribosome assembly, Translation factors, mt-tRNA synthetases, fMet processing

mt_ls_oxphos_txl <- mt_ls_allpaths[c("Complex_I", "Complex_II", "Complex_III", "Complex_IV", "Complex_V",
                                     "Mitochondrial_ribosome", "Mitochondrial_ribosome_assembly", "Translation_factors", "mt-tRNA_synthetases", "fMet_processing")]

# look at number of overlapping between mt_ls_oxphos_txl and genes in M3 genes
for (n in names(mt_ls_oxphos_txl)){
  print(paste0("Number of overlap of ", n, ": ", length(intersect(subset(MM_tbl, Module == 3)$gene_name, mt_ls_oxphos_txl[[n]]))))
}
# "Number of overlap of Complex_I: 21"
# "Number of overlap of Complex_II: 0"
# "Number of overlap of Complex_III: 5"
# "Number of overlap of Complex_IV: 13"
# "Number of overlap of Complex_V: 4"
# "Number of overlap of Mitochondrial_ribosome: 24"
# "Number of overlap of Mitochondrial_ribosome_assembly: 8"
# "Number of overlap of Translation_factors: 4"
# "Number of overlap of mt-tRNA_synthetases: 1"
# "Number of overlap of fMet_processing: 0"

set.seed(123) # need to set.seed right before running to get the same results
fgseaRes_oxphos_txl <- fgsea(pathways = mt_ls_oxphos_txl, 
                  stats    = ranks,
                  scoreType = "pos",
                  minSize  = 5,
                  maxSize  = 500)
fgseaRes_oxphos_txl <-as.data.frame(fgseaRes_oxphos_txl) # this shows ComplexI genes (NDUFs) and mitochondrial ribosomal genes are significantly enriched

# our initial reason of choosing NDUFs was that 
# there's significant enrichment of "NDUF hub genes" among all the hub genes mapped to a GO biologcal process pathways (p=0.003008)
# this was done in 2_bl_wgcna.R
# now we also test if the is true for the mitochondrial ribosomal genes

ribosomal_genes <- intersect(M3_sort$gene_name, mt_ls_oxphos_txl$Mitochondrial_ribosome) #24 genes
length(intersect(subset(M3_sort, MM>=0.7)$gene_name, ribosomal_genes)) #4 hub ribosomal genes

## there are 265 background genes mapped to a GO BP, 22 of them are ribosomal genes and 38 of them are hub genes; 4 of the 22 ribosomal genes are "hub ribosomal genes mapped to a GO BP"
fisher.test(cbind(c(4, 18),  c(34, 211)), alternative = "greater") # p = 0.3838


write.xlsx(list(
  "BroadPaths" = fgseaRes_broadpaths,
  "OXPHOS_TXL" = fgseaRes_oxphos_txl),
  paste0(outdir, "/1_7/fgseaRes.xlsx")) # this is in Supplementary Table 11



