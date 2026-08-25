#!/usr/bin/env Rscript

source("~/bin/myRpackages.R")
projDir <- "./R01_Resilience_WGS_022321/"
inDir <- paste0(projDir, "input/from_Adelyn/100125/")
outDir <- paste0(projDir, "output/")
analDir <- paste0(outDir, "analyses/")
annoDir <- paste0(outDir, "annotation/")
qcDir <- paste0(outDir, "QC/Blood/")
qtlDir <- paste0(analDir, "NDUFs_eQTL/")

# read in a file containing NDUF gene names and IDs
nduf.bb <- read_csv(paste0(inDir, "nduf_bb.csv")) %>% rename(GeneID=gene_id, GeneSymbol=gene_name)
# read in a file contain gene loci information
geneLoci <- read.delim("./Ensembl/hg38.Gene_features.GC_Content.txt", header = T, stringsAsFactors = F, quote = "") %>% right_join(nduf.bb)
geneLoci$Chr <- as.numeric(gsub("chr","",geneLoci$Chr))
# load dataframe containing NDUFs residuals correcting for only technical covariates; these are prepped using rebuttal/rbt_1_5.Rmd
load(paste0(inDir, "Resilience_Residuals_eQTL.RData"), verbose = T)
# load dataframe for covariates including sex and age
covar.rnaseq.R01 <- read.delim(paste0(inDir, "mcsa_covars_eqtl.txt"), header = T, stringsAsFactors = F, quote = "")
covar.array.adni <- read.delim(paste0(inDir, "adni_covars_eqtl.txt"), header = T, stringsAsFactors = F, quote = "")
# read in files containing genetic PCs
pcs.wgs.adni <- read.table("./ADNI_082619/output/QC/ADNI_WGS/ADNI_WGS_CSHR_GBHNFL.pca.evec", skip = 1, header = F, stringsAsFactors = F)
colnames(pcs.wgs.adni)<- c("IID","EV1", "EV2", "EV3", "EV4", "EV5", "EV6", "EV7", "EV8", "EV9", "EV10", "Dx")

pcs.wgs.R01 <- read.delim(paste0(qcDir, "R01_Resilience_WGS_CSHRP_VNGBHF_Trimmed_Covariates_100721.tsv"), header = T, stringsAsFactors = F, quote = "")%>% select(FID, IID, EV1:EV3)

#############################################################
#prepare Pheno files for cisQTL

resid.ndufs.R01 <- mcsa_res_eqtl %>% 
  select(SubjectID, all_of(geneLoci$GeneID)) %>% 
  mutate(ptnum=as.numeric(SubjectID)) %>% 
  left_join(covar.rnaseq.R01 %>% select(FID, IID, ptnum)) %>% 
  select(FID, IID, matches("ENSG"))
resid.ndufs.adni <- adni_res_eqtl %>% 
  select(SubjectID, all_of(geneLoci$GeneID)) %>% 
  mutate(FID=0, IID=toupper(SubjectID)) %>% 
  select(FID, IID, matches("ENSG")) 

path.resid.adni <- paste0(qtlDir, "NDUFs_ADNI.pheno")
path.resid.R01 <- paste0(qtlDir, "NDUFs_R01_Blood.pheno")
write.table(resid.ndufs.adni, path.resid.adni,row.names = F, col.names = T, sep = "\t", quote = F)
write.table(resid.ndufs.R01, path.resid.R01,row.names = F, col.names = T, sep = "\t", quote = F)

#############################################################
#prepare Covar files for cisQTL


covar.wgs.R01 <- covar.rnaseq.R01 %>% mutate(Sex=ifelse(Sex=="Male", 0, 1)) %>% left_join(pcs.wgs.R01) %>% select(FID, IID, Age, Sex, matches("EV"))
covar.wgs.adni <- covar.array.adni %>% mutate(FID=0, IID=toupper(SubjectID), Sex=ifelse(Sex=="Male", 0, 1)) %>% select(FID, IID, Age, Sex) %>% left_join(pcs.wgs.adni %>% select(IID, EV1, EV2, EV3) %>% mutate(IID=toupper(IID))) 

path.covar.adni <- paste0(qtlDir, "NDUFs_ADNI.covar")
path.covar.R01 <- paste0(qtlDir, "NDUFs_R01_Blood.covar")
write.table(covar.wgs.adni, path.covar.adni, row.names = F, col.names = T, sep = "\t", quote = F)
write.table(covar.wgs.R01, path.covar.R01, row.names = F, col.names = T, sep = "\t", quote = F)

#############################################################
#Plink genetic files

inPlink.wgs.R01 <- "./R01_Resilience_WGS_022321/output/QC/Blood/R01_Resilience_WGS_CSHRP_VNGBHF"
inPlink.wgs.adni <- "./ADNI_082619/output/liftOver/ADNI_WGS_CSHRP_GBHNF_UpdatedIIDs_HG38"

#############################################################

rm(adni_res_eqtl, covar.array.adni, covar.rnaseq.R01, mcsa_res_eqtl, nduf.bb, pcs.wgs.R01, pcs.wgs.adni)


#############################################################

# run eQTL
for(i in 1:nrow(geneLoci)){
  #i <- 1
  chrom <- geneLoci[i, "Chr"]
  start <- geneLoci[i, "Start"]-1000000
  stop <- geneLoci[i, "Stop"]+1000000
  pheno <- geneLoci[i, "GeneID"]
  geneName <- geneLoci[i, "GeneSymbol"]
  
  system(paste0("plink2 --bfile ", inPlink.wgs.adni, " --chr ", chrom, " --from-bp ", start, " --to-bp ", stop, " --linear hide-covar --ci 0.95 --keep ", path.resid.adni, " --covar-variance-standardize --set-all-var-ids @:# --pheno ", path.resid.adni, " --pheno-name ", pheno, " --covar ", path.covar.adni, " --out ", qtlDir, "ADNI"))
  
  system(paste0("plink2 --bfile ", inPlink.wgs.R01, " --chr ", chrom, " --from-bp ", start, " --to-bp ", stop, " --linear hide-covar --ci 0.95 --keep ", path.resid.R01, " --covar-variance-standardize --set-all-var-ids @:# --pheno ", path.resid.R01, " --pheno-name ", pheno, " --covar ", path.covar.R01, " --out ", qtlDir, "R01_Blood"))
  
  res.adni <- read.delim(paste0(qtlDir, "ADNI.", pheno,".glm.linear"), header = T, stringsAsFactors = F, quote = "", check.names = F) %>% filter(!is.na(P)) %>%  rename(CHR=`#CHROM`, BP=POS, SNP=ID, A2=OMITTED, MAF=A1_FREQ, NMISS=OBS_CT) %>% select(CHR, BP, SNP, A1, A2, MAF, NMISS, BETA, SE, P)
  
  res.R01 <- read.delim(paste0(qtlDir, "R01_Blood.", pheno,".glm.linear"), header = T, stringsAsFactors = F, quote = "", check.names = F) %>% filter(!is.na(P)) %>%  rename(CHR=`#CHROM`, BP=POS, SNP=ID, A2=OMITTED, MAF=A1_FREQ, NMISS=OBS_CT) %>% select(CHR, BP, SNP, A1, A2, MAF, NMISS, BETA, SE, P)
  
  path.res.adni <- paste0(qtlDir, "ADNI_", pheno, "_metaInput.tsv")
  path.res.R01 <- paste0(qtlDir, "R01_Blood_", pheno, "_metaInput.tsv")
  
  write.table(res.adni, path.res.adni, row.names = F, col.names = T, sep = "\t", quote = F)
  write.table(res.R01, path.res.R01, row.names = F, col.names = T, sep = "\t", quote = F)
  
  system(paste0("plink --meta-analysis ", path.res.adni, " ", path.res.R01, " + qt weighted-z --out ",qtlDir, pheno))
  
  res.meta <-  read.table(paste0(qtlDir, pheno, ".meta"), header = T, stringsAsFactors = F, quote = "", check.names = F) %>% mutate(gene_id=pheno, gene_name=geneName) %>% select(gene_id, gene_name, everything())
  
  if(i==1){
    comb.res <- res.meta
  }else{
    comb.res <- bind_rows(comb.res, res.meta)
  }
}
system(paste0("rm ", qtlDir,"*.{meta,prob,log,linear,tsv}"))
# annotated the results 
anno <- readRDS(paste0(annoDir, "R01_Resilience_Blood_WGS_CSHRP_VNGBHF.rds")) %>% select(CHR, POS, REF, ALT, rsid, matches("refGene"), SIFT_pred:CADD_phred)

anno.comb.res <- comb.res %>% left_join(anno %>% rename(BP=POS))
saveRDS(anno.comb.res, file=paste0(qtlDir, "ADNI_R01_Resilience_Blood_cisQTL_metaAnalysis_021726.rds"))
