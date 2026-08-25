# mitochondrial-resilience
This repository cotains analysis and visualization scripts used in the manuscript **Discovery of Genes Underlying Cognitive Resilience in Individuals Predisposed to Alzheimer’s Risk (Tsai et al., 2026)**.  
Dataset page for the manuscript is located at: https://www.synapse.org/Synapse:syn68156041/datasets/.  
Availability of input data used in this manuscript is outlined in the methods and data availability sections of the manuscript.\
_If you use these codes in your research, please cite our manuscript._ 

# Repository structure
**Main**\
`1_CreatResids.R` : create residuals of gene expression data from MCSA and ADNI.\
`2_bl_wgcna.R` : blood (MCSA and ADNI) consensus WGCNA analysis.\
`3_br_wgcna.R` : brain (AMP-AD) consensus WGCNA analysis.\
`4_nduf_assoc_de.R` : DEG analysis of _NDUFs _in blood and brain cohorts.\
`5_longitudinal.R` : longitudinal association of _NDUFs_ in MCSA and ADNI.\
`visualization.R` : create all main and supplementary figures to visualize results.\
`fns.R` : a list of frequently used functions

**rebuttal**: Scripts generated in response to reviewers' comments.\
`rbt_1_1_new.Rmd` : This script addresses comment 1 about cognitive resilience in the AMP-AD datasets from reviewer 1.\
`rbt_1_2.R` : This script addresses comment 2 about hippocampus and entorhinal cortex specificity from reviewer 1.\
`rbt_1_3.R` : This script addresses comment 3 about mediation from reviewer 1.\
`rbt_1_4_cbx_LM10.Rmd` : This script address comment 4 about cell type proportion from reviewer 1.\
`rbt_1_5.Rmd` :  This script addresses comment 5 about genetic support from reviewer 1.\
`rbt_1_5_run_cisQTL_NDUFs_021726.R` : This script is used for running cis-eQTL and addresses comment 5 about genetic support from reviewer 1.\ 
`rbt_1_6.R` : This script addresses comment 6 about network robustness from reviewer 1.\
`rbt_1_7.R` : This script addresses comment 7 about competitive gene testing from reviewer 1.\
`rbt_1_10_PRSv3.Rmd` : This script addresses comment 10 about genetic risk interaction from reviewer 1.\
`rbt_1_10.R` : This script addresses comment 10 about genetic risk interaction from reviewer 1. This script in particular runs interaction with _APOE_4 dose.\
`rbt_1_11.R` : This script addresses comment 11 about cross-platform harmonization from reviewer 1.\
`rbt_1_12.R` : This script addresses comment 12 about effect size reporting and forest plots from reviewer 1.\
`rbt_1_13.R` : This script addresses comment 13 about outlier diagnostics from reviewer 1.\
`rbt_2_1.Rmd` : This script addresses comment 1 about analysis correcting for diagnosis from Reviewer 2.\
`rbt_2_1_3_BF_TM.Rmd` : This script addresses comment 1 about number of resilience per diagnosis from Reviewer 2. This script also addresses comment B and F from Reviewer 3 regarding prediction on resilience status.\
`rbt_2_3.Rmd` : This script addresses comment 3 about _APOE_4 sensitivity analysis Reviewer 2.\
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


# Contact
For further information about the codes or the project, please contact tsai.wei@mayo.edu or Taner.Nilufer@mayo.edu
