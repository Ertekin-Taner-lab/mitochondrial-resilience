# Mitochondrial-Resilience
This repository cotains analysis and visualization scripts used in the manuscript **Discovery of Genes Underlying Cognitive Resilience in Individuals Predisposed to Alzheimer’s Risk (Tsai et al., 2026)**.  
Dataset page for the manuscript is located at: https://www.synapse.org/Synapse:syn68156041/datasets/.  

_If you use these codes in this Github, the following manuscript must be cited and the grants need to be included in the Acknowledgement:_ \
**Paper citation**: Discovery of Genes Underlying Cognitive Resilience in Individuals Predisposed to Alzheimer’s Risk (Tsai et al., 2026)\
**Grant Acknowledgement**: This work was supported by the National Institutes of Health and National Institute on Aging [RF1 AG051504, U01 AG046139, R01 AG061796, U19 AG074879] and Alzheimer’s Association Zenith Fellows Award [ZEN-22-969810] to Dr. Nilufer Ertekin-Taner.

# Repository structure
**Main**\
`1_CreatResids.R` : create residuals of gene expression data from MCSA and ADNI.\
`2_bl_wgcna.R` : blood (MCSA and ADNI) consensus WGCNA analysis.\
`3_br_wgcna.R` : brain (AMP-AD) consensus WGCNA analysis.\
`4_nduf_assoc_de.R` : DEG analysis of _NDUFs_ in blood and brain cohorts.\
`5_longitudinal.R` : longitudinal association of _NDUFs_ in MCSA and ADNI.\
`visualization.R` : create all main and supplementary figures to visualize results.\
`fns.R` : a list of frequently used functions

**rebuttal**: Scripts generated in response to reviewers' comments.\
`rbt_1_1_new.Rmd` : This script addresses comment 1 about cognitive resilience in the AMP-AD datasets from reviewer 1.\
`rbt_1_2.R` : This script addresses comment 2 about hippocampus and entorhinal cortex specificity from reviewer 1.\
`rbt_1_3.R` : This script addresses comment 3 about mediation from reviewer 1.\
`rbt_1_4_cbx_LM10.Rmd` : This script address comment 4 about cell type proportion from reviewer 1.\
`rbt_1_5.Rmd` :  This script addresses comment 5 about genetic support from reviewer 1.\
`rbt_1_5_run_cisQTL_NDUFs_021726.R` : This script is used for running cis-eQTL and addresses comment 5 about genetic support from reviewer 1. 
`rbt_1_6.R` : This script addresses comment 6 about network robustness from reviewer 1.\
`rbt_1_7.R` : This script addresses comment 7 about competitive gene testing from reviewer 1.\
`rbt_1_10_PRSv3.Rmd` : This script addresses comment 10 about genetic risk interaction from reviewer 1.\
`rbt_1_10.R` : This script addresses comment 10 about genetic risk interaction from reviewer 1. This script in particular runs interaction with APOE4 dose.\
`rbt_1_11.R` : This script addresses comment 11 about cross-platform harmonization from reviewer 1.\
`rbt_1_12.R` : This script addresses comment 12 about effect size reporting and forest plots from reviewer 1.\
`rbt_1_13.R` : This script addresses comment 13 about outlier diagnostics from reviewer 1.\
`rbt_2_1.Rmd` : This script addresses comment 1 about analysis correcting for diagnosis from Reviewer 2.\
`rbt_2_1_3_BF_TM.Rmd` : This script addresses comment 1 about number of resilience per diagnosis from Reviewer 2. This script also addresses comment B and F from Reviewer 3 regarding prediction on resilience status.\
`rbt_2_3.Rmd` : This script addresses comment 3 about APOE4 sensitivity analysis Reviewer 2.\
`rbt_2_9.R` : This script addresses comment 9 about hippocampal volume harmonization from reviewer 2.

**rebuttal_2nd**: Scripts generated in response to the second round of reviewers' comments.\
`rbt_1_2.Rmd` : This script addresses comment 2 about resilience in ROSMAP from reviewer 1.\
`rbt_1_3.Rmd` : This script addresses comment 3 about regional specificity from reviewer 1.\
`rbt_1_4.Rmd` : This script address comment 4 about cell type proportion from reviewer 1.\
`rbt_1_5.Rmd` :  This script addresses comment 5 about Combat harmonization from reviewer 1.\
`rbt_1_6.Rmd` :  This script addresses comment 6 about prediction analysis from reviewer 1.

Below summarizes the corresponding codes related to supplementary tables that contain analytical results:
|        **Table**       |                                     **Codes**                                    |
|:----------------------:|:--------------------------------------------------------------------------------:|
|  Supplementary Table 3 |                            `2_bl_wgcna.R` & `3_br_wgcna.R`                       |
|  Supplementary Table 4 |                            `2_bl_wgcna.R` & `3_br_wgcna.R`                       |
|  Supplementary Table 5 |                                `rebuttal/rbt_1_12.R`                              |
|  Supplementary Table 7 |                                `rebuttal/rbt_1_13.R`                              |
|  Supplementary Table 8 |                                `rebuttal/rbt_1_11.R`                              |
|  Supplementary Table 9 |                                   `3_br_wgcna.R`                                   |
| Supplementary Table 10 |                             `rebuttal_2nd/rbt_1_3.Rmd`                           |
| Supplementary Table 11 |                                `rebuttal/rbt_1_7.R`                               |
| Supplementary Table 12 |                                 `4_nduf_assoc_de.R`                               |
| Supplementary Table 13 |                             `rebuttal_2nd/rbt_1_4.Rmd`                           |
| Supplementary Table 14 |                             `rebuttal_2nd/rbt_1_6.Rmd`                           |
| Supplementary Table 15 |                                `rebuttal/rbt_1_3.R`                               |
| Supplementary Table 16 | `rebuttal/rbt_1_5_run_cisQTL_NDUFs_021726.R`;<br>`rebuttal/rbt_1_5.Rmd` (MR section) |
| Supplementary Table 19 |                            `2_bl_wgcna.R` & `3_br_wgcna.R`                       |


# Data availability
All human de-identified data in this manuscript is available via the [AD Knowledge Portal](https://adknowledgeportal.synapse.org). The AD Knowledge Portal is a platform for accessing data, analyses and tools generated by the Accelerating Medicines Partnership (AMP-AD) Target Discovery Program and other National Institute on Aging (NIA)-supported programs to enable open-science practices and accelerate translational learning. Data is available for general research use according to the following requirements for data access and data attribution (https://adknowledgeportal.synapse.org/DataAccess/Instructions). The bulk blood RNAseq data generated in this study will be deposited and made available in the AD Knowledge Portal upon publication. Other datasets used in this study were downloaded as described in the Methods or provided by collaborators. 

# Acknowledgement
We would like to thank the patients and their families for their participation, without whom these studies would not have been possible. We thank Thuy Nguyen and Kimberly Malphrus for extracting RNA from MCSA blood samples. We thank Jeremy Syrjanen and Steve Smith for their help with providing MCSA data and answering related questions.  We thank the Mayo Clinic Genome Analysis Core (GAC) director, Stephen Murphy, PhD, and supervisor Julie Lau, for their collaboration in collection of omics data.

## Funding 
This work was supported by the National Institutes of Health and National Institute on Aging [RF1 AG051504, U01 AG046139, R01 AG061796, U19 AG074879 to NET; P30 AG062677, U01 AG006786 to RCP; R01 AG075959, R01 AG082348, R01 AG081951, R01 AG057739, R01 AG070883, U01 AG024904, R01 LM013463, T32 AG071444, U24 AG074855, U01 AG068057, U01 AG072177, U19 AG074879, P30 AG010133, P30 AG072976, R01 AG019771, R01 AG057739, U19 AG024904, R01 LM013463, R01 AG068193 to AJS; R01 AG081951, R01 LM012535, U01 AG072177, and U19 AG074879 to KN] and Alzheimer’s Association Zenith Fellows Award [ZEN-22-969810 to NET].

## Datasets
**MayoRNAseq** :  The Mayo RNAseq study data was led by Dr. Nilüfer Ertekin-Taner, Mayo Clinic, Jacksonville, FL as part of the multi-PI U01 AG046139 (MPIs Golde, Ertekin-Taner, Younkin, Price). Samples were provided from the following sources: The Mayo Clinic Brain Bank. Data collection was supported through funding by NIA grants P50 AG016574, R01 AG032990, U01 AG046139, R01 AG018023, U01 AG006576, U01 AG006786, R01 AG025711, R01 AG017216, R01 AG003949, NINDS grant R01 NS080820, CurePSP Foundation, and support from Mayo Foundation. Study data includes samples collected through the Sun Health Research Institute Brain and Body Donation Program of Sun City, Arizona. The Brain and Body Donation Program is supported by the National Institute of Neurological Disorders and Stroke (U24 NS072026 National Brain and Tissue Resource for Parkinsons Disease and Related Disorders), the National Institute on Aging (P30 AG19610 Arizona Alzheimers Disease Core Center), the Arizona Department of Health Services (contract 211002, Arizona Alzheimers Research Center), the Arizona Biomedical Research Commission (contracts 4001, 0011, 05-901 and 1001 to the Arizona Parkinson's Disease Consortium) and the Michael J. Fox Foundation for Parkinsons Research. 

**Mount Sinai Brain Bank** : These data were generated from postmortem brain tissue collected through the Mount Sinai VA Medical Center Brain Bank and were provided by Dr. Eric Schadt from Mount Sinai School of Medicine.

**ROSMAP** : Study data were provided by the Rush Alzheimer’s Disease Center, Rush University Medical Center, Chicago. Data collection was supported through funding by NIA grants P30AG10161 (ROS), R01AG15819 (ROSMAP; genomics and RNAseq), R01AG17917 (MAP), R01AG30146, R01AG36042 (5hC methylation, ATACseq), RC2AG036547 (H3K9Ac), R01AG36836 (RNAseq), R01AG48015 (monocyte RNAseq) RF1AG57473 (single nucleus RNAseq), U01AG32984 (genomic and whole exome sequencing), U01AG46152 (ROSMAP AMP-AD, targeted proteomics), U01AG46161(TMT proteomics), U01AG61356 (whole genome sequencing, targeted proteomics, ROSMAP AMP-AD), the Illinois Department of Public Health (ROSMAP), and the Translational Genomics Research Institute (genomic). Additional phenotypic data can be requested at www.radc.rush.edu.

**RNAseq Harmonization Study (rnaSeqReprocessing)** : Data generation was supported by the following NIH grants: P30AG10161, P30AG72975, R01AG15819, R01AG17917, R01AG036836, U01AG46152, U01AG61356, U01AG046139, P50 AG016574, R01 AG032990, U01AG046139, R01AG018023, U01AG006576, U01AG006786, R01AG025711, R01AG017216, R01AG003949, R01NS080820, U24NS072026, P30AG19610, U01AG046170, RF1AG057440, and U24AG061340, and the Cure PSP, Mayo and Michael J Fox foundations, Arizona Department of Health Services and the Arizona Biomedical Research Commission. We thank the participants of the Religious Order Study and Memory and Aging projects for the generous donation, the Sun Health Research Institute Brain and Body Donation Program, the Mayo Clinic Brain Bank, and the Mount Sinai/JJ Peters VA Medical Center NIH Brain and Tissue Repository. Data and analysis contributing investigators include Nilüfer Ertekin-Taner, Steven Younkin (Mayo Clinic, Jacksonville, FL), Todd Golde (University of Florida), Nathan Price (Institute for Systems Biology), David Bennett, Christopher Gaiteri (Rush University), Philip De Jager (Columbia University), Bin Zhang, Eric Schadt, Michelle Ehrlich, Vahram Haroutunian, Sam Gandy (Icahn School of Medicine at Mount Sinai), Koichi Iijima (National Center for Geriatrics and Gerontology, Japan), Scott Noggle (New York Stem Cell Foundation), Lara Mangravite (Sage Bionetworks).

**UCI_5XFAD Study** : The IU/JAX/UCI MODEL-AD Center was established with funding from The National Institute on Aging (U54 AG054345-01 and AG054349). Aging studies are also supported by the Nathan Shock Center of Excellence in the Basic Biology of Aging (NIH P30 AG0380770).

**RADC** : We thank the study participants and staff of the Rush Alzheimer’s Disease Center. ROSMAP is supported by P30AG10161, P30AG72975, R01AG15819, R01AG17917, U01AG46152, and U01AG61356. ROSMAP resources can be requested at https://www.radc.rush.edu and www.synpase.org. All ROSMAP participants are enrolled without known dementia and agreed to detailed clinical evaluation and brain donation at death. Both ROS and MAP studies were approved by an Institutional Review Board of Rush University Medical Center. Each participant signed an informed consent, Anatomic Gift Act, and an RADC Repository consent allowing their data and biospecimens to be repurposed.

**ADNI dataset**: Data collection and sharing for this project was funded by the Alzheimer’s Disease Neuroimaging Initiative (ADNI) (National Institutes of Health Grant U01 AG024904) and DOD ADNI (Department of Defense award number W81XWH-12-2-0012). ADNI is funded by the National Institute on Aging, the National Institute of Biomedical Imaging and Bioengineering, and through generous contributions from the following: AbbVie, Alzheimer’s Association; Alzheimer’s Drug Discovery Foundation; Araclon Biotech; BioClinica, Inc.; Biogen; Bristol-Myers Squibb Company; CereSpir, Inc.; Cogstate; Eisai Inc.; ElanPharmaceuticals, Inc.; Eli Lilly and Company; EuroImmun; F. Hoffmann-La Roche Ltd. and its affiliated company Genentech, Inc.; Fujirebio; GE Healthcare; IXICO Ltd.; Janssen Alzheimer Immunotherapy Research & Development, LLC.; Johnson & Johnson Pharmaceutical Research & Development LLC.; Lumosity; Lundbeck; Merck & Co., Inc.; Meso Scale Diagnostics, LLC.; NeuroRx Research; Neurotrack Technologies; Novartis Pharmaceuticals Corporation; Pfizer Inc.; Piramal Imaging; Servier; Takeda Pharmaceutical Company; and Transition Therapeutics. The Canadian Institutes of Health Research is providing funds to support ADNI clinical sites in Canada. Private sector contributions are facilitated by the Foundation for the National Institutes of Health (www.fnih.org). The grantee organization is the Northern California Institute for Research and Education, and the study is coordinated by the Alzheimer’s Therapeutic Research Institute at the University of Southern California. ADNI data are disseminated by the Laboratory for Neuro Imaging at the University of Southern California.

# Contact
For further information about the codes or the project, please contact tsai.wei@mayo.edu or Taner.Nilufer@mayo.edu
