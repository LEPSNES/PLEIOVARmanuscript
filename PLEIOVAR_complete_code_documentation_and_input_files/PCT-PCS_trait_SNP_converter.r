

#-------- This program takes a specific PC-Trait and a specific PC-SNP and ---------#
#-------- find the most correlated Trait and most correlated SNP -------------------#
#-------- to get a bette interpretation --------------------------------------------#

 A <- read.table("CONFIG",header=FALSE); A <- as.list(A)
 #---- output folder ; gene region extension , flag for related related individuals inflation correcttion ----#
 outfolder <- A$V1[[2]]; results_folder = A$V1[[11]]
 
 #GENE = "CLIC1" 
 #GENE = "CPD"
 GENE = "COL4A1"

 #----- Z2 Table ------#
 Z2_filename = paste("/data/Osorio/PLEIOVAR_TEST/OUTPUT/Z2_set0/Z2-file_",GENE,sep="")
 Z2_TABLE = read.table(Z2_filename,header=TRUE); Z2_TABLE
 if( ncol(Z2_TABLE) == 1) { Z2_TABLE = t(Z2_TABLE) }
 
 #--------- getting trait file (random intercept estimates of mixed-model------------#
 TRAIT_filename <- "/data/Osorio/PLEIOVAR_TEST/TRAIT_FILES/mainpheno_0.dat"
 TRAITS <- read.table(TRAIT_filename,header=TRUE)
 #--------- getting PC-TRAITS --------#
 PCT_filename <- "/data/Osorio/PLEIOVAR_TEST/TRAIT_FILES/PC-Traits_0.txt"
 PCT <- read.table(PCT_filename,header=FALSE); colnames(PCT)[1] = "ID"
 #------- getting loadings for PCT-Traits ---------------#
 PCT_loadings <- read.table(paste("/data/Osorio/PLEIOVAR_TEST/TRAIT_FILES/PC-Traits_0_loadings.txt",sep=""),header=TRUE)
 #------ getting loadings for PC-SNPs in the gene -------#
 PC_SNP_loadings_filename <- paste(outfolder,"/PC-SNPS_B_loadings/PC-SNP_load_",GENE,sep="")
 PC_SNP_loadings <- read.table(PC_SNP_loadings_filename,header=FALSE)
 #------ getting PC-SNPs --------------------------------# 
 PCSNP_data_filename     <- paste(outfolder,"/PC-SNPS_B/PC-SNP_",GENE,sep="")
 PCS <- read.table(PCSNP_data_filename,header = TRUE)
 #----- getting the gene region file -------------------#
 Assemble_filename  <- paste(outfolder,"/Assemble/",GENE,"_assembled",sep="")
 SNPFILE <- read.table(Assemble_filename,header=TRUE)
 #----- preparing to merge PC-Traits, PC-SNPs and SNPFILE to get common ID's ------#
 TEMP1 = cbind(SNPFILE$ID,c(1:nrow(SNPFILE))) ; colnames(TEMP1) = c("ID","SEQ1")
 TEMP2 = cbind(PCS$ID,c(1:nrow(PCS))) ; colnames(TEMP2) = c("ID","SEQ2")
 TEMP3 = cbind(PCT$ID,c(1:nrow(PCT))) ; colnames(TEMP3) = c("ID","SEQ3")
 TEMP4 = cbind(TRAITS$ID,c(1:nrow(TRAITS))) ; colnames(TEMP4) = c("ID","SEQ4")

 JOINT = merge(TEMP1,TEMP2) ; JOINT = merge(JOINT,TEMP3) ; JOINT = merge(JOINT,TEMP4)
 SEL1 = JOINT[,2]; SEL2 = JOINT[,3] ; SEL3 = JOINT[,4] ; SEL4 = JOINT[,5]  
 SNPFILE = SNPFILE[ SEL1,]; PCS = PCS[ SEL2,]; PCT = PCT[ SEL3, ];  TRAITS = TRAITS[ SEL4, ]

 #--- parameters for Trait and SNP search ---------#
 PCT_num = 1  ; PCS_num = 24

 #----- most correlated SNP ------#
 P1 = PCS[,PCS_num + 1] ; P0 = SNPFILE[,-1]; C = cor(P1,P0) 
 MAX = max(abs(C),na.rm=TRUE) ; MAX 
 POS = which(abs(C)==MAX) ; POS = min(POS) ; POS 
 C[POS]
 colnames(SNPFILE)[POS]  
 P2 = P0[,POS]  
 
 C = cor(TRAITS[,-1],P2); C ; nr = nrow(TRAITS) ; ZT = C*sqrt((nr-2)/(1-C^2)) ;  ZT 
 C = cor(PCT[,-1],P2) ; C ; nr = nrow(PCT) ; ZPCT = C*sqrt((nr-2)/(1-C^2)) ;  ZPCT
 colnames(SNPFILE)[POS]

 #----- correlated TRAITS with PCT -------------#
 
 C = cor(TRAITS[,-1],PCT[,-1]) ; C

 #----- best SNPs as predictors of PC_Trait ----#
 MODEL <- lm(PCT[,2]  ~ P0[,6] + P0[,84] + P0[,127]  ) 
 summary(MODEL) 







 










 




