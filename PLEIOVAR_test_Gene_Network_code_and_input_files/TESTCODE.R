
  #setwd("T:\\LEPS\\NES\\newStaff\\OSORIO\\PLEIOVAR\\ASHG Submission 20230711\\PLEIOVAR test Gene Network - code and input files")

  #***********************************************************************************************************************************************#
  #***********************************************************************************************************************************************#
  #--------------- THIS SCRIPT AIMS AT RUNNING ALL MAIN STEPS FOR THE SCORING OF A GIVEN GENE SET USING A GIVEN SET OF GENES ---------------------#
  #***********************************************************************************************************************************************#
  #***********************************************************************************************************************************************#

  #------ INPUT (1): results_inflamation_biomarkers
  #------ this file has one line per gene for all 22397 genes and it has the PLEIOVAR results for 
  #------ inflamation biomarkers as the phenotype dataset (please see Excel Workbook in same folder)
  #------ INPUT (2): oxygen_transport (this file has a list of 14 genes from the Gene Ontology 
  #------ biological process and represents the most significant output from submitting 
  #------ all 92 genes that had PLEIOVAR false discovery rate <= 0.01  
 
  #------ use the program TESTCODE to run it by pasting it into the R console -----#

  #------ the output will be the p-value after adjusting for variance inflation factor (VIF)---#
  #------ Due to the very large files using all SNPS in each gene for the estimation of VIF
  #------ we estimated VIF in our server and then adjusted for it to get the final p-value
  #------ we also output the p-value without VIF adjustment which is expected to be 
  #------ much smaller, specially if VIF is much larger than one (in our case it is 2.643) ----#

  #------ at the end of this script we present a way to estimate the p-value for the gene set 
  #------ for situations when the R package is not capable of calculating it since it 
  #------ is too small for using the R pnorm function

 #***********************************************************#
 #---------------------- TEST CODE --------------------------#
 #***********************************************************#

 #---- results dataset for the 11 inflamatory traits for SardiNIA: 
 # CRP
 # IL6
 # erythrocyte sedimentation rate 
 # White blood cells
 # Platelets
 # Red blood cells
 # eosinophil count
 # basophil count
 # lymphocytes count
 # monocyte count
 # neutrophil level

 B=read.table("results_inflamation_biomarkers",header=FALSE)
 colnames(B)=c("chr","start","end","gene","start_reg","end_reg","n_PCT","n_PCS","SSQ","DF","p-value")
 dim(B)

 #-----------------------------------------------------------------------------------------------#
 #---- need to define the set of genes in top GO network ----------------------------------------#
 #---- use Excel sheet Network Scoring Example to get GO sets of genes (gene networks)-----------#
 #---- get output from GO -----------------------------------------------------------------------#
 #---- top biological process was "oxygen transport" with 14 genes (fdr = 2.2E-14) --------------#
 #-----------------------------------------------------------------------------------------------#
 
 C = read.table("oxygen_transport.txt",header=FALSE); colnames(C)[1] = "genes"
 C = as.matrix(C) ; C= as.vector(C) ;  NROW = length(C)
 
 #---- need to estimate the variance inflation factor (VIF)
 #---- for top GO network set of genes ("oxygen transport") 
 #---- For this, we used SNPs within each gene in the set (large files which are omitted in the  
 #---- supplemental materials) but were run on our server to estimate the VIF (2.643) and 
 #---- then the adjustment factor which is 1.626, where the adjustment factor is the sqrt of VIF 
 
 adjustment_factor = 1.626

 #-- use C genes (gene set) and locate in B (PLEIOVAR results) to estimate TSSQ and TDF 
 #-- where TSSQ is the sum of SSQ for each gene and TDF is the sum of DF for each gene 

 TSSQ = 0; TDF = 0
 for(i in 1:NROW)
 {
     gene <- as.character(C[i]) 
     POS = which( B$gene == gene) 
     if( length(POS) > 0)
     {
        SSQ  = B$SSQ[POS] ; DF  = B$DF[POS]
        TSSQ = TSSQ + SSQ ; TDF = TDF + DF
     }
 }
 TSSQ
 TDF

 #----------------------------- TSSQ is the total sum of squares ----------------------------------------------------------#
 #----  How do we use the adjusmtnt factor?  We first need to estimate the z-score corresponding to TSSQ 
 #----  To do this we would calculate the p-value for TSSQ and then get the z-score and the divide the zscore  
 #----  by the adjustment factor.  However, the p-value for TSSQ under the Chi-square might overflow if TSSQ is large    
 #----  Canal approximation of the Chi-sqaure if a way to obtain the z-score with having to first obtain the p-value  

 Q = TSSQ/TDF
 L = Q^(1/6) - 0.5*Q^(1/3) + (1/3)*Q^(1/2)
 EL = (5/6) - (1/9)*(1/TDF) - (7/648)*(1/TDF)^2 
 VL = (1/18)*(1/TDF) + (1/162)*(1/TDF)^2 - (37/11664)*(1/TDF)^3  

 Z = -(L - EL)/sqrt(VL) ;  pnorm(Z,0,1)
 Z_adj = Z/adjustment_factor ; pnorm(Z_adj,0,1) 

 #--------------------------- Here we have the adjusted z-score -----------------------------------------------------------#
 #--- however, we can still run into the problem of overflow when Z_adj is too large in absolute value 
 #--- To solve this we use a Q-function to estimate the log(Q-function(Z_adj))
 #--- as the Q-function gives a very good approximation for the p-value of a Normal c.d.f for large |Z| scores
 #--- when Z_adj is large in magnitude. Thus, PVALUE = (1/sqrt(2*3.14159))*(abs(Z_adj)/(1+Z_adj^2))*exp( -(Z_adj^2)/2 ) 
 
 minus_Log_PVALUE = -log((1/sqrt(2*3.14159))) -log(abs(Z_adj)/(1+Z_adj^2)) + (Z_adj^2)/2
 minus_Log10_PVALUE = log10(2.718281828)*minus_Log_PVALUE
 power = floor(minus_Log10_PVALUE)
 frac = minus_Log10_PVALUE-power
 power
 base = paste("E-",power+1,sep="")
 mantissa = 10^(1-frac)
 mantissa = floor(100*mantissa)/100
 final = paste("Extreme p-value = ",mantissa,base,sep="")
 final
 pnorm(Z_adj,0,1)

 #------ as a final summary, our final p-value for the gene set "oxygen transport" is 8.26E-32
 #------ without VIF correction our p-value would have been incorrectly 
 #------ overly optimistic (1.066E-80). This example illustrates that importance for correction 
 #------ for VIF, primarily caused by neighboring genes in the same chromossome 
 #------ FOr example, of the set of 14 genes in oxygen transport, there were 5 neighboring genes in chr 11
 #------ and 4 neighboring genes in chr 16 which explains why VIF was equal to 2.643 when 
 #------ it should be close to one if genes where far away from one another (if in the same chromossome)
 #------ or in different chromossomes (and therefore their corresponding PC-SNPs 
 #------ would be independent from one another).  

