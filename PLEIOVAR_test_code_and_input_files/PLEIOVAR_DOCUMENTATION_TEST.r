
  #***********************************************************************************************************************************************#
  #***********************************************************************************************************************************************#
  #--------------- THIS SCRIPT AIMS AT RUNNING ALL MAIN STEPS FOR THE ASSOCIATION INVOLVING AN INITIAL TRAIT FILE AND SINGLE GENE REGION ---------#
  #***********************************************************************************************************************************************#
  #***********************************************************************************************************************************************#

  #------ INPUT (1): mainpheno.dat  (this file has one trait per ID) it can be obtained by either using the residual after controlling for covariates of a cross-sectional dataset 
  #------ or by estimating a random intercept estimate if using a longitudinal dataset also after controlling for covariates 
  #------ INPUT (2): CETP_assembled (this is the gene region file for CETP containing allele dosages for each SNP 
  #------ INPUT (3): small file with the variance inflation for each PC-trait (assumes previously that all genes were   
  #------ Note that to run this script, all inuput files should be located at the same folder of the R console
  #------ where this script is located.  This, all we need to do is paste this script to get the results  

  #------------ KEY PARAMTERS ---------#
  #---- par1: proportion of variance explained for using dimension reduction when generating the PC-traits file
  #---- par2: proportion of variance explained for using dimension redcution when generating the PC-SNPs file 
  #---- par3: minimum MAF for a SNP to be considered for running the PC-SNP step
  #---- par4: flag for adjusting for variance inflation 


  #------------ MAIN OUTPUT -----------#
  #---- 1-PC-Traits file 
  #---- 2-Loadings for PC-Traits 
  #---- 3-variance explained for each PC-Trait 
  #---- 4-PC-SNPs file
  #---- 5-loadgings for PC-SNPs 
  #---- 6-variance explained for each PC-SNP 
  #---- 7-Association results (SSQ, DF, p-value) between the PC-Traits file and the PC-SNP file for CETP gene 
  #---- 8-Z2 Table (adjusted or unadjusted for genomic control type of variance inflation)

  
  #----- reading main input files -----------#
  phenofile =  read.table("TRAIT_TEST.dat",header=TRUE)
  SNPfile   =  read.table("CETP_TEST.dat",header=TRUE) 
  #--- varaince inflation file assumes that PLEIVOAR association has been run for all genes (vs. PC-Traits file) and that the VIF was estimated for each PC-Trait -------#
  VIF       =  read.table("PC-Traits_TEST_VIF.dat",header=FALSE)

  #----- setting main paramters ---------------------------------------------------------------#
  #----- variance explained for PC-Traits dimension reduction ---------------------------------#  
  par1 = 0.75
  #----- variance explained for PC-SNPs dimension reduction -----------------------------------#
  par2 = 0.75
  #----- Minimum Minor Allele Frequency to be used as a SNP in the generation of PC-SNPs ------#
  par3 = 0.005
  #----- Flag to perform the variance inflation adjustment ------------------------------------#
  par4 = 1


  #*******************************************#
  #----- generation of PC-Traits file --------#
  #*******************************************#
 
  OPF_trait = par1 ;  ID <- phenofile[,1]
  #--- take out the ID's from traits file ----#
  TRAITS <-  phenofile[ , c(-1)]
  #--- get number of columns and rows from trait file without the ID ----#
  ncolT <- ncol(phenofile)-1 ; nrowT <- nrow(phenofile)
  #----- if only one trait, then standardize and assign 1 to its loadings -----#
  if (ncolT == 1)
  {
     M <- mean(phenofile[,2]);  S <- sd(phenofile[,2]);  OTrait <- (phenofile[,2]-M)/S ;   loadings = "PC-Trait-1 1.00"
  } else {
    #--- Generating PC-Traits ----------------####
    prW <- prcomp(TRAITS);  OTrait <- prW$x; loadings <- prW$rotation ; loadings = floor(loadings*1000+0.5)/1000
    #--- Get % variance from each PC-TRAIT  ---#
    OPF_s <- apply(OTrait,2,var);  OPF_s <- OPF_s/sum(OPF_s); lambda = floor(OPF_s*1000+0.5)/1000 
    #----- write the lambdas to a file -------#
    write.table(lambda,"PC-trait_lambdas",col.names=FALSE)
    #--- defining the vector of cummulative % variance for PC-TRAITS ---#
    OPFcum_s <- matrix(0,1,ncolT); CUM <- 0;
    for(i in 1:ncolT) {  CUM <- CUM + OPF_s[i]; OPFcum_s[i] <- CUM }
    #----- select PC-traits that meet variance explained cutoff -----#
    tflag <- 1*(OPFcum_s >= OPF_trait); select_trait <- min((1:ncolT)[tflag>0]); OTrait <- OTrait[,1:select_trait]; FLAG = 2
  }
  
  #---- roundoff PC-Traits to neighrest hundreds -------------------------------------------------#
  OTrait = sign(OTrait)*floor(100*abs(OTrait) + 0.5)/100 
  #---- join the ID's with the PC-Traits ---------------------------------------------------------#
  ID <- as.matrix(ID);  OTraitB <- cbind(ID,OTrait); colnames(OTraitB)[1] = "ID" 
  write.table(OTraitB,"PC-Traits.txt",col.names = FALSE , row.names = FALSE, quote = FALSE)
  #----- write PC-Trait loadings to file. Here we differentiated between having 1 or more PC-traits saving rownames for the single PC-trait ------#
  #----- although I am not sure this differentiation will be needed when using more powerfull software tools -------------------------------------#
  write.table(loadings,"PC-Trait_loadings",col.names=FALSE,row.names=TRUE,quote=FALSE)

  #*******************************************#
  #----- generation of PC-SNPs file ----------#
  #*******************************************#

   OPF_snp = par2 ; freq_cut = par3 ;  X  <- SNPfile 
   #---- condition if there is at least one SNP (the first column is the ID, so this is why we use > 1) ----------#
   if (dim(X)[2] > 1)
   {
       #--- if more than one SNP, eliminate the first columnm otherwise just calculate the MAF of the single SNP ----#
       if (dim(X)[2] > 2)
       { 
           freq <- apply(X[,-1],2,mean)/2; MAF = freq*(freq <= 0.5) + (freq-0.5)*(freq > 0.5) 
       } else { 
           freq = mean(X[,2])/2 ;  MAF = freq*(freq <= 0.5) + (freq-0.5)*(freq > 0.5)
       }
       #---- elim is the vector that that has the column number of SNPS that did not meet the minimum varance Threshold
       elim <- which(MAF < freq_cut)
       #----- eliminate the SNPs that did not meet the threshold ----#
       if(length(elim) > 0 )
       { 
           cutpos <- elim + 1;  X <- X[ ,-cutpos] 
       }
       X = as.matrix(X)
   }

   if (dim(X)[2] > 1)
   {
          #----- only do if at least two SNPs at this stage -----#
          if (dim(X)[2] > 2)
          {
             #---- assinging X to S and eliminate the ID column ---------#
             S <- X ; ID <- S[,1] ; S <-  S[ , c(-1)]; S<- as.matrix(S); rownames(S) <- NULL;  colnames(S) <- NULL
             #-- getting the dimensions of SNP dataset --#
             d <- dim(S);  nrows <- d[1] ;  ncols <- d[2]
             #-- insert small noise to overcome LD as SNPs with 100% correlation will crash when running PCA ---#
             set.seed(1); noise = runif(nrows*ncols,-0.0001,0.0001); dim(noise) = c(nrows,ncols); S = S + noise
             #--- principal components of SNPs ---#
             pr <- prcomp(S); SPC <- pr$x;  colSPC <- dim(SPC)[2]; loadings <- pr$rotation
             #--- roundoff loadings to the fourth decimal and get lambdas (the varaince explained for each PC-SNP)-----#
             loadings <- floor(10000*loadings +0.5)/10000; lambda <- (pr$sdev)^2
             #--- getting the cummulative variance of PC-SNPs -------------------------------#
             OPF_s <- lambda/(sum(lambda)); OPFcum_s <- matrix(0,1,colSPC); CUM <- 0
             for( j in 1:colSPC)
             {
                   CUM <- CUM + OPF_s[j];  OPFcum_s[j] <- CUM
             }
             #--- number of PC-SNPs needed to meet the varaince explained cutoff (#PC-SNP's with cummulative var explained > cutoff (0.75 in our case) ----#
             select_snp <- min(which(OPFcum_s >= OPF_snp))
             #----- rounding off varaince explained  ----#
             OPFcum_s <- floor(1000*OPFcum_s +0.5)/1000
             #----- OSNP is the set of PC-SNPS after dimension reduction
             OSNP   <- SPC[,1:select_snp]; OSNP <- as.matrix(OSNP)
             #---- rounding to the neighrest hundreth ------
             OSNP <- floor(100*OSNP+0.5)/100;  colnames(OSNP)=rep("",select_snp);
             #--- lambda now is the variance explained  ----#
             lambda = OPF_s ; lambda = floor(1000*lambda+0.5)/1000
              #---- generatnig the names for the PC-SNPS (PC-SNP1, PC-SNP2, ... )
             for(j in 1:select_snp)
             {
               colnames(OSNP)[j] <- paste("PC_SNP-",j,sep="")
             }
             OSNP <- cbind(ID,OSNP)
             loadings <- loadings[1:select_snp,]
          } else {
             #--- if there was only a single SNP, then the PC-SNP is the SNP itself and lambda and the loadings is equal to 1 -----#
             OSNP <- X; OPFcum_s <- 1; loadings <- 1 ; lambda <- 1
          }
          #---- write OSNP (the PC-SNPs) to the PC-SNP file -----#
          write.table(OSNP,"PC-SNP_CETP",col.names=TRUE,row.names=FALSE,quote=FALSE)
          #------  write the loadings to the loadings file -----------------------#
          write.table(loadings,"PC-SNP_CETP_loadings",row.names=TRUE,col.names=FALSE,quote=FALSE)
          #---- write the variance explained (lambdas) ---------------------------#
          write.table(OPFcum_s,"PC-SNP_CETP_varexp",col.names=FALSE,row.names=FALSE,quote=FALSE)
  }  



  #*******************************************#
  #--------- Association Step ----------------#
  #*******************************************#

  #--- merge PC-Traits (OtraitB) with PC-SNPs (OSNP) , calculate the number of degrees of freedom usind in correlation --#   
    
  n_pctraits = dim(OTraitB)[2]-1 ; n_pcsnps = dim(OSNP)[2] - 1 ;   JOIN  = merge(OTraitB,OSNP) ; JOIN = JOIN[,-1] ;  NDF = dim(JOIN)[1] - 2;

  PCT = JOIN[ , c(1:n_pctraits)] ; JOIN = JOIN[ , -c(1:n_pctraits)  ] ; PCS = JOIN ; PCT = as.matrix(PCT) ; PCS = as.matrix(PCS)  

  C <- cor(PCT,PCS); Z2 <- (C^2)/(1-C^2+0.0001)*NDF  

  #**************************************************
  #------- variance inflation adjustment of Z2 -----#
  #**************************************************   
  
  flag_VIF = par4  ;  lambda = VIF[,2] ; Z2B = matrix(0,n_pctraits,n_pcsnps)   

  if(flag_VIF)
  {
     #----- correcting Z2 associations using stored variance inflation factors and generating the Z2PRE_adj which will be used next to calculate SSQ and its p-value -----#
     #----- here the Z2 for each PC-Trait (row) is divided by its corresponding VIF 
     for(i in 1:n_pctraits) 
     { 
         Z2B[i,] = Z2[i,]/lambda[i]  
     }      
  } else {
     Z2B = Z2
  }

  Z2B = floor(100*Z2B+0.5)/100    
  DF = n_pctraits*n_pcsnps ; SSQ = sum(Z2B) ;  pvalue <- pchisq(SSQ,DF,lower.tail=FALSE) ; pvalue <- formatC(pvalue,digits = 2,format="E")
  
  #*******************************
  #----- saving main results ----#
  #*******************************

  OUTPUT = c("CETP",n_pctraits,n_pcsnps,SSQ,DF,pvalue)  
  OUTPUT = rbind( c("gene","n_PCT","n_PCS","SSQ","DF","pvalue"),OUTPUT) 
  OUTPUT

  write.table(OUTPUT,"ASSOCIATION_RESULTS",col.names=FALSE,row.names=FALSE,quote=FALSE) 
 
  #***********************************#
  #---- saving final Z2 Table --------#
  #***********************************#
  
  name_col = rep("",n_pcsnps) ; name_row = rep("",n_pctraits)
  for(i in 1:n_pcsnps)
  {
    name_col[i] = paste("PCS-",i,sep="")
  }
  for(i in 1:n_pctraits)
  {
    name_row[i] = paste("PCT-",i,sep="")
  }

  colnames(Z2B) = name_col ; rownames(Z2B) = name_row  
  write.table(Z2B,"CETP_Z2_Table",col.names=TRUE,row.names=TRUE) 


 
  








 
