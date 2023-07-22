
 #-----------------------------------------------------------------------------------#
 #--------------- This program scores a set of gene speficied by the user -----------#
 #-- it will eventully become a function which can be called tousands of times ------# 
 #-- where each loop is a set of genes used to call the function --------------------#  
 


 #---- folder where the gene dosage files are located (one file per each gene) -----#
 infolder = "/data/Osorio/VALIDATION_PLEIOVAR_SARDINIA_NEW"
 #----- simplified hg file (containing CETP) and some of its neighbors - total 10 genes -----#
 #---- for complete file please see All_genes_hg19 in this same outfolder -------------------#
 infile <- "Sample_hg19"; AA <- read.table(infile,header=FALSE); colnames(AA) = c("chr","start","end","gene"); nrowA <- dim(AA)[1]; ncolA <- dim(AA)[2]; extension = 50000
 #---- read PC-TRAITS file -------#
 PCT <- read.table("PC-Traits.txt",header=FALSE); colnames(PCT)[1] <- "ID"; nPCT <- dim(PCT)[2]-1 ; #nrow = dim(PCT)[1]
 for(m in 1:nPCT) { colnames(PCT)[m+1] <- paste("PC_TRAIT-",m,sep="")} ; count = 0


 #----- The loop below merges all PC-SNPs from the specified set of genes in Sample_hg19 -----#
 #----- the set of genes will be used to obtain a correlation matrix -------------------------#
 #----- which will then be used to adjust the inflation of the sum of SSS --------------------#
 #----- from all genes in the network due to gene proximity of some of the -------------------#
 #----- genes in the network that are in the same chromossome or genes -----------------------#
 #----- in the same chromossome that are correlated with one another -------------------------# 

 for( i in 1:nrowA)
 {
     chr <- AA[i,1]; gene <- AA[i,4]; SNPfile <- paste(infolder,"/PC-SNPS_B/PC-SNP_",gene,sep="")
     if(file.exists(SNPfile))
     {
         count = count + 1
         if(count == 1)
         { 
            OSNP = read.table(SNPfile,header=TRUE);  colnames(OSNP)[1] <- "ID"
         } else {
            TEMP = read.table(SNPfile,header=TRUE); colnames(OSNP)[1] <- "ID" ; nc1 = ncol(TEMP) -1 ; nc0 = ncol(OSNP)-1
            temp = NULL
            for(j in 1: nc1)
            {
              temp[j] = paste("PC_SNP.",j+nc0,sep="")
            } 
            colnames(TEMP)[ c(2:(nc1+1))]  = temp  ;   OSNP = merge(OSNP,TEMP)                 
        } 
     }
 }

 #-------------------------------------------------------------------------------------------------#
 #---- adjustment factor is obtained by obtaining the ratio between -------------------------------#
 #----- NUMERATOR  :  Twice the Sum of the squares of the correlations ----------------------------# 
 #----- DENOMINATOR:  Twice the snumber of degrees of freedom -------------------------------------#
 #---- under NULL, the expected value of the NUMERATOR is 2*DF (the DENOMINATOR) ------------------# 
 #---- adjustment factor is the square root of (NUMERATOR/DENOMINATOR) ----------------------------#

 DF = ncol(OSNP)-1; TEMP = OSNP[,-1];CORREL_MATRIX = cor(TEMP)^2;n=nrow(OSNP);NUM=sum(2*CORREL_MATRIX);DEN=2*DF+2*DF*(DF-1)/(n-1)
 var_sum_inflation= NUM/DEN; adjustment_factor = sqrt( var_sum_inflation )

 #--- getting total SSQ (sum across all SSQ across genes) with real data -----------------------------------------#
 #--- after looping each gene in the set, finding it on results-file and getting the corresponding SSQ and DF ----# 
 MAINSSQ = read.table("results_LIPIDS",header=FALSE); MAINSSQ_gene = MAINSSQ$V4 ;  TSSQ = 0; TDF = 0
 for(i in 1:nrowA)
 {
     gene <- as.character(AA[i,4]) ; print(" i = ");print(i)
     if( nchar(gene) > 0)
     {
        POS = which( gene == MAINSSQ_gene ) 
        SSQ = MAINSSQ$V9[POS]; DF  = MAINSSQ$V10[POS]
        TSSQ = TSSQ + SSQ ; TDF = TDF + DF
     }
 }

 #----------------------------------- TSSQ is the total sum of squares ------------------------------------------------------------------#
 #-----  How do we use the adjustment factor?  We first need to estimate the z-score corresponding to TSSQ -------------------------------#
 #-----  To do this we would calculate the p-value for TSSQ and then get the z-score and the divide the zscore --------------------------# 
 #-----  by the adjsutment factor.  However, the p-value for TSSQ under the Chi-square might overflow if TSSQ is large ------------------#   
 #------ Canal approximation of the Chi-sqaure if a way to obtain the z-score with having to first obtain the p-value  ------------------#

 Q = TSSQ/TDF
 L = Q^(1/6) - 0.5*Q^(1/3) + (1/3)*Q^(1/2)
 EL = (5/6) - (1/9)*(1/TDF) - (7/648)*(1/TDF)^2 
 VL = (1/18)*(1/TDF) + (1/162)*(1/TDF)^2 - (37/11664)*(1/TDF)^3  

 Z = -(L - EL)/sqrt(VL) ;  pnorm(Z,0,1)
 Z_adj = Z/adjustment_factor ; pnorm(Z_adj,0,1) 

 #--------------------------- Here we have the adjusted z-score -------------------------------------------------------------------------#
 #---- however, we can still run into overflow of Z_adj is too large in absolute value --------------------------------------------------#
 #--- To solve this we use a Q-function to estimate the log(Q-function(Z_adj)) ----------------------------------------------------------# 
 #--- as the Q-function give a very good approximation for the p-value of a Normal c.d.f ------------------------------------------------#
 #--- when Z_adj is large in magnitude. Thus, PVALUE = (1/sqrt(2*3.14159))*(abs(Z_adj)/(1+Z_adj^2))*exp( -(Z_adj^2)/2 ) -----------------#  
 #--- code below should also be a fucntion ----------------------------------------------------------------------------------------------# 

  #---- double check from manuscript -----#

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


 

