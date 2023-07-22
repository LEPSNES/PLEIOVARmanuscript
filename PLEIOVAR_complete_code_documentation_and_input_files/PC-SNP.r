
#---- This programs aims at selecting the SNPS from the SardiNIA vcf file which correspond to a specific gene region 
#---- After preprocessing was performed fom previous steps, we have a gene region dosage file for each gene 
#---- block files (one for each CPU when using multiple jobs in parallel) contain a list of all the genes in the block in UCSC format (chr,start,end,gene)  
#---- So, in summary, these block files are really subsets of the file All_Genes_hg19 which contains all the genes 
#---- Here is an example for the block file number 27 (out of a total of 64 blocks) 
#---- 7 137074384 137531609 DGKI
#-----7 137597556 137686847 CREB3L2
#-----7 137638093 137642712 LOC100130880
#-----7 137761177 137803050 AKR1D1

#---- Of the block file described above, we use only the gene name so we can locate the corresponding gene region file {gene name}_assembled 
#---- Next, for each gene region (in the block) , read the "assembled" file for that gene and then we generate the PC-SNPs for each gene.  
#-----Next, we perform a dimension reduction which results from a much smaller set of PC-SNPs than the number of SNPs in that gene region 

#---- Now we will show the main values from the CONFIG file and explain its current settings --------#

#---- value 1 (not used in this program) folder address where vcf files are located (/data/Osorio/VALIDATION_PLEIOVAR_EPACTS/VCF_files)
#---- value 2 This is the folder where the output results will be written to (/data/Osorio/VALIDATION_PLEIOVAR_SARDINIA_NEW) 
#---- value 3 (not used in this program) Number of digits (2 in this case) after rounding off the dosages (this is written to the SAMPLE files
#---- value 4 (not used in this program) Number of SNPS (500 in this case) in the MICROFILES (subsets of the much larger SAMPLE_{chr} files)
#---- value 5 (not used in this program) Number of CPU's (64 in this case) used in the parallel processing (this is also the same as the number of block files) 
#---- value 6 (not used in this program) Number of Kb (50 in our case) to be used to define a gene regions (start - 50kb, end + 50kb in our case)   
#---- value 7 variance explained parameter (0.75 in our case) for PC-SNP dimension reduction 
#---- value 8 minor allele frequency cutoff for a SNP to be used for PC-SNP (in our case this is 0.005)
#---- value 9 variance explained parameter (0.99 in our case) for PC-Trait dimension reduction 
#---- value 10 variance inflation corretion flag (to correct or not for variance inflation) In our case it is a 1 (since SardiNIA uses it as some members are from the same family). 
#---- value 11 is the sub-sample size for running PC-SNP-FAST   

 #------- open the CONFIG file -----#
 A <- read.table("CONFIG",header=FALSE); A <- as.list(A)
 outfolder <- A$V1[[2]]; OPF_snp <- A$V1[[7]] ; OPF_snp <- as.character(OPF_snp)   ; OPF_snp <- as.numeric(OPF_snp)
 freq_cut <- A$V1[[8]]; freq_cut <- as.character(freq_cut) ; freq_cut <- as.numeric(freq_cut);
 #---- VARCUT is a cutoff that would exclude SNPs that might have a variance lower than  a threhold.  
 #---- The advanteage of using this rather than minor allele frequency is that we could have an allele frequency that is greater than the threshold 
 #---- but has very low variance (and consequently, low information).  Not sure though if this makes much difference compared to using a MAF cutoff instead.   
 VARCUT=2*freq_cut*(1-freq_cut)
 #----- main arguments ----#
 args   <- commandArgs(); block_name <- args[8]; 
 #block_name <- "block64" 
 #--- infile is the folder where the block files are located 
 infile <- paste(outfolder,"/Regions/",block_name,sep="")
 #--- AA is the block file 
 AA <- read.table(infile,header=FALSE); 
 #--- nrowA is the number if rows if AA and ncolA is the number of columns of AA ----# 
 nrowA <- dim(AA)[1]; ncolA <- dim(AA)[2]
 #--- here we loop for every gene in the block file ----#
 for (i in 1:nrowA)
 {
    #----- get chr and gene region file called {gene name}_assembled ---
    chr <- AA[i,1]; gene <- AA[i,4]; SNPfile <- paste(outfolder,"/Assemble/",gene,"_assembled",sep="") 
    #------ output folder for PC-SNP_{gene name} files and saved on the PC-SNP folder  
    outfile1 <- paste(outfolder,"/PC-SNPS/PC-SNP_",gene,sep="")
    #------ output folder for PC-SNP_load_{gene name}.  Loadings files (paramters that generate the PC-SNPs when combined with the original SNPs (dosages) ----#
    outfile2 <- paste(outfolder,"/PC-SNPS_loadings/PC-SNP_load_",gene,sep="")
    #------ output folder for PC-SNP_var_{gene name} files.  These contain the cummulative variance explained (var explained PC-1, var explained PC-1 and PC-2, ... ) 
    outfile3 <- paste(outfolder,"/PC-SNPS_var_exp/PC-SNP_var_",gene,sep=""); 
    #---- lambdas are just the variance explained (non-cummulative) such as var-explained PC-1, var explained PC-2, .... 
    outfile4 <- paste(outfolder,"/PC-SNPS_lambdas/PC-SNP_lambdas_",gene,sep="")
    # -----------  READ X MATRIX on FILE -------------------* 
    
    #---- the SNPfile (gene dosages for the gene) will only exist if there was at least one SNP in the gene from the SardiNIA dosage dataset  
    if (file.exists(SNPfile)) 
    {
       X  <- read.table(SNPfile,header=TRUE) 
       #---- condition if there is at least one SNP (the first column is the ID, so this is why we use > 1) ----------#
       if (dim(X)[2] > 1)
       {
          #--- if more than one SNP eliminate the first columnm otherwise 
          #---- just calculate the variance of the single SNP ----#   
          if (dim(X)[2] > 2){ VarSnps <- apply(X[,-1],2,var) } else { VarSnps <- var(X[,2]) }
          #---- elim is the vector that that has the column number of SNPS that did not meet the minimum varance Threshold 
          elim <- which(VarSnps < VARCUT)
          #----- eliminate the SNPs that did not meet the threshold ----#
          if(length(elim) > 0){ cutpos <- elim + 1;  X <- X[ ,-cutpos] }; X <-as.matrix(X)
       }
       #----- This next step looks like a repetion of the previous but it needed to be done
       #----- because at this point, it might not have 2 or more SNPs any longer
       #----- due to the exclusion by variance cutoof above.  I might not have even a single 
       #----- SNP at this point, so we need to ask again about these conditions. 
       
       #---- only do if at least one SNP at this stage ----------#
       if (dim(X)[2] > 1)
       {       
          #----- only do if at least two SNPs at this stage -----#
          if (dim(X)[2] > 2) 
          {  
             #---- assinging X to S ---------#
             S <- X ; ID <- S[,1] ; 
             #---- eliminate the ID column ----# 
             S <-  S[ , c(-1)]; S<- as.matrix(S); rownames(S) <- NULL;  colnames(S) <- NULL
             #-- getting the dimensions of SNP dataset --# 
             d <- dim(S);  nrowS <- d[1] ;  ncolS <- d[2]
           
             #-- insert small noise to overcome LD as SNPs with 100% correlation will crash when running PCA ---#
             #tot <- dim(S)[2]; a <- runif(tot,-0.000000001,0.000000001);  b <- seq(1:dim(S)[1]); c <- runif(dim(S)[1],0,1); b <- b[order(c)]
             #for(j in 1:tot) { S[b[j],j] <- S[b[j],j] + a[j]  }      
              
             #-- insert small noise to overcome LD as SNPs with 100% correlation will crash when running PCA ---#            
             noise = runif(nrowS*ncolS,-0.0001,0.0001); dim(noise) = c(nrowS,ncolS); S = S + noise
 
             #--- principal components of SNPs ---#
             pr <- prcomp(S); SPC <- pr$x;  colSPC <- dim(SPC)[2]; loadings <- pr$rotation 
             #--- roundoff loadings to the fourth decimal and get lambdas (the varaince explained for each PC-SNP)-----#
             loadings <- floor(10000*loadings +0.5)/10000; lambda <- (pr$sdev)^2
             #--- getting the cummulative variance of PC-SNPs -------------------------------#
             OPF_s <- lambda/(sum(lambda)); OPFcum_s <- matrix(0,1,colSPC); CUM <- 0      
             for( j in 1:colSPC) 
             {
                   CUM <- CUM + OPF_s[j]
                   OPFcum_s[j] <- CUM 
             }
             #---- GAOS's estimation of # of independent SNPs and # of PC-SNPS meeting cutoff ----#
             indep_SNPs <- min(which(OPFcum_s >= 0.995)) 
             #--- number of PC-SNPs needed to meet the varaince explained cutoff (#PC-SNP's with cummulative var explained > cutoff (0.75 in our case) ----#
             select_snp <- min(which(OPFcum_s >= OPF_snp)) 
             #----- calculating cummulative variance explained to save into file later ----# 
             OPFcum_s <- OPFcum_s[1:indep_SNPs]; OPFcum_s <- floor(1000*OPFcum_s +0.5)/1000
             #----- OSNP is the set of PC-SNPS after dimension reduction 
             OSNP   <- SPC[,1:select_snp]; OSNP <- as.matrix(OSNP)
             #---- rounding to the neighrest hundreth ------
             OSNP <- floor(100*OSNP+0.5)/100;  colnames(OSNP)=rep("",select_snp); 
             #--- lambda now is the variance explained only for the PC-SNPS that met the threshold ----#
             lambda = OPF_s[1:select_snp]
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
          #----- roundoff lambda ------#
          lambda = floor(1000*lambda+0.5)/1000
          #---- write OSNP (the PC-SNPs) to the PC-SNP file -----#
          #write.table(OSNP,outfile1,col.names=TRUE,row.names=FALSE,quote=FALSE)
          saveRDS(OSNP,file=outfile1) 
          #------  write the loadings to the loadings file -----------------------# 
          #write.table(loadings,outfile2,row.names=TRUE,col.names=FALSE,quote=FALSE)
          saveRDS(loadings,file=outfile2)
          #---- write the cummulative lambdas into the PC-SNP_var files  ---------#
          #write.table(OPFcum_s,outfile3,col.names=FALSE,row.names=FALSE,quote=FALSE)
          saveRDS(OPFcum_s,file=outfile3)
          #------  write the lambdas to the lambda file PC-SNP_lambda files-------# 
          #write.table(lambda,outfile4,col.names=FALSE,row.names=FALSE,quote=FALSE)
          saveRDS(lambda,file=outfile4)
       }
    }             
}





            
