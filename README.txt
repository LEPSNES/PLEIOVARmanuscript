#----------------------------------------------------------------------------------------------------#
#--------------- This document details the main steps for running PLEIOVAR  -------------------------# 
#----------------------------------------------------------------------------------------------------#

#---- MAIN FOLDERS:
#     PLEIOVAR COMPLETE - Code and input files including all pre-processing steps.  
#          This folder also has a sub-folder \Extra with the complete gene set scoring script and a script (still in development) to identify an original trait and original SNP given a PCT and a PCS 
#     PLEIOVAR COMPLETE DOCUMENTATION - Main documentation file "PLEIOVAR DOCUMENTATION  2023_05_07.pdf (update)"
#     PLEIOVAR TEST     - code and input files for running a PLEIOVAR test dataset 
#     PLEIOVAR test Gene Network - code and input files for scoring a pre-specified set of genes (gene network) 
#     Supplemental Data - Supplemental Tables that were too large to fit on a single page
#     Supplemental Materials - Supplemental Tables and Figures


#     Main code, input files and output files -----#

#---- PLEIOVAR COMPLETE - At least 80% of all the code, includes multiple input files, R scripts and Perl scripts
#     designed for using multiple processors.
#     Also, includes all pre-processing for (e.g. 11.8 million SNP's in SardiNIA) generating: 
#     gene region files (e.g. 22,397 genes in SardiNIA),
#     generation of PC-SNPs with dimension reduction for each gene (22,397 PC-SNP files)
#     generation of PC-Traits with dimension reducion
#     association (p-values) between each PC-Trait file and each of the 22,397 PC-SNP files  
#     Z2 Tables for each gene (22,397 files) adjusted or not for variance inflation
#     Geeration of joint score for a set of genes

#---- PLEIOVAR TEST -----#
#-----Less than 20% of all code, includes three input files and an R script that generates
#     PC-Traits with dimension reduction (6 PC-Traits out of the initial 8 phenotypes based on Glycemia,
#     Insulinemia",HbA1c","HDL","LDL","Triglycerides","SBP","DBP"), PC-SNPs with dimension reduction for CETP gene (chosen as an example) and 
#     the association between PC-Traits and the set of 10 PC-SNPs generated from CETP SNP data 


# ---- PLEIOVAR test Gene Network - code and input files ------#
# ---- Used for scoring a pre-specified set of genes (gene network) 
# ---- Uses a predefined set of 14 genes (from the Gene Ontology biological process "oxygen transport" 
#      based on PLEIOVAR results of a set of 11 inflamatory biomarkers)
#      Outputs p-value for the gene network before and after adjusting for variance inflation (VIF)  









