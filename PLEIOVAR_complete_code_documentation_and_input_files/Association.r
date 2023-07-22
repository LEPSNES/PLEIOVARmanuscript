
 #------ NEWEST VERSION for ASHG journal --------------#


 #-----------------------------------------------------#
 #------- THis program is the association step --------#
 #--- whole association takes aprox. 3.5 minutes ------#
 #-----------------------------------------------------#

 #------- CONFIG file is open with key parameters for the association ----#
 A <- read.table("CONFIG",header=FALSE); A <- as.list(A)
 #---- outputfolder is the top folder where PC-SNP fiels are located -----#
 outfolder <- A$V1[[2]]
 #----- extension is the amount of Kb to subtract from the gene start position and to add to the gene end position ----------------------------#
 extension <- A$V1[[6]]; extension <- as.character(extension) ; extension <- 1000*as.numeric(extension)
 #----- flag_inflattion_correction is the flag to perform the variance inflation correction (generally should be set to one as the default ) ---#
 flag_inflation_correction <- A$V1[[10]]; flag_inflation_correction <- as.character(flag_inflation_correction); flag_inflation_correction <- as.numeric(flag_inflation_correction)
 #---- open UCSC gene file where chr, gene, start, end fields are captured -------#
 infile <- "All_Genes_hg19";  AA <- read.table(paste("./",infile,sep=""),header=FALSE); nrowA <- dim(AA)[1]; ncolA <- dim(AA)[2]
 #---- read PC-TRAITS file -------#
 PCT <- read.table(paste(outfolder,"/PC-Traits.txt",sep=""),header=FALSE); colnames(PCT)[1] <- "ID"; ncolT <- dim(PCT)[2]-1
 #---- assign columns names for PC-traits ------#
 for(j in 1:ncolT) { colnames(PCT)[j+1] <- paste("PC_TRAIT-",j,sep="")}
 #---- This version does not merge OSNP with PCT, but instead, uses the index files (two), containing the order of the OSNP 
 #---- and the order of PCT which were in the intersection between OSNP and PCT.  This step will be done once (on the first non-null PC-SNP file) since the 
 #---- overlap of ID's between each PC-SNP file and PCT file is the same      
 #---- initilizing output variables ----#
 OUTPUT <- NULL ; 
 #--- ZPRE is the matrix of Z2 before adjustment for variance inflation.  This ZPRE has a constant number of rows (number of PC-Traits) but 
 #--- will grow as we will append to it, the Z2 matrix for the latest gene, so in the end it will have a few PC-Traits but thousands of PC-SNPs  
 #--- POS1 is the beggining position in ZPRE of the Z2 that is appended
 #--- POS2 is the end position in ZPRE of the Z2 that is appended ( then next POS1 will be the previous POS2 + 1 )
 Z2PRE  <- NULL ; Z2PRE_adj <- NULL; count <- 0; POS1 = rep(0,50000); POS2 = rep(0,50000); GENE <- NULL
 for (i in 1:nrowA)
 {
   #---- get gene name and open PC-SNP file ----#
   chr <- AA[i,1]; gene <- AA[i,4]; ini <- AA[i,2]; fin <- AA[i,3] ;begpos <- max(1,ini-extension); endpos <- fin+extension
   SNPfile <- paste(outfolder,"/PC-SNPS/PC-SNP_",gene,sep="")
   # -----------  READ X MATRIX on FILE -------------------*
    if(file.exists(SNPfile))
    {
       if(count==0)
       { 
           #--- if first PC-SNP file read, merge with PC-trait and from the merged file, get the PC-SNP sequence number and PC-Traits sequence number -------------# 
           TEMP1 = readRDS(file=SNPfile); TEMP1 = as.data.frame(TEMP1); TEMP2 = as.data.frame(PCT); c1 = dim(TEMP1)[1] ; c2 = dim(TEMP2)[1] ; TEMP1$seqS = c(1:c1) ; TEMP2$seqT = c(1:c2)
           IDX=merge(TEMP1,TEMP2);P1 = which(colnames(IDX)== "seqS") ; P2=which(colnames(IDX)== "seqT");IDX = IDX[,c(P1,P2)]; 
           #---- PC-traits file (PCT) is now reduced for the rest of the program by using only the cases that merged with PC-SNP from PC-SNP sequecne number ------#
           PCT= PCT[IDX[,2],]; PCT= PCT[,-1]; PCT=as.data.frame(PCT) ;  select_trait <- dim(PCT)[2]; NDF <- nrow(IDX)-2
           #--- at this point we have the index file IDX which consists of "seqS" the row number of the original PC-SNP (before merge) and "seqT" the
           #--- original row number of PCT (before merge.  This step is only performed once (when count = 0)
           #--- PCT at this point has already selected only the cases numbers in seqT (in which ID matched the ID's in the first PC-SNP file), and excluded the ID 
       }       
       #--- here we have the sequence number of PC-SNP and then PC-Traits so we can just use reduced files (wihout merge) to get the correlation ----------------# 
       #--- since the SNPfile exists, increase the count by one 
       count = count + 1       
       #--- read the PC-SNP file, take out the ID, select only the cases with IDX[,1] row numbers 
       OSNP <- readRDS(file=SNPfile); OSNP = OSNP[,-1]; OSNP = as.matrix(OSNP); OSNP = OSNP[IDX[,1],];OSNP = as.matrix(OSNP)
       #--- generated the T2 matrix between PCT and OSNP which is approximately a Z2 matrix
       #--- select_snp is the number of PC-SNPs
       select_snp <- dim(OSNP)[2]; C <- cor(PCT,OSNP); Z2 <- (C^2)/(1-C^2+0.0001)*NDF
       #--- Z2B is just the Z2 matrix rounded to the neighrest hundreth ----#
       Z2B <- floor(1000*Z2+0.5)/1000; #if(ncolT==1) { Z2B <- t(Z2B) } 
       GENE[count] = gene ; OUTPUT <- rbind(OUTPUT,paste(chr,ini,fin,gene,begpos,endpos,select_trait,select_snp,sep=" "))
       #--- saving the position of the Z2B gene grid for the large Z2 matrix  which has k PC-traits (3 in this study) and their association vs. all PC-SNPs (thousands) ------
       if(count==1)  { POS1[count] = 1  } else { POS1[count] = POS2[count-1] + 1 }
       POS2[count] = POS1[count] + select_snp - 1
       #---- append Z2B to the long horizontal ZPRE matrix ------#
       Z2PRE = cbind(Z2PRE,Z2B) 
    }    
}


#------ getting the median of the Z2PRE association for each PC-trait across thousands of PC-SNPs  -----------#   
POS1 = POS1[1:count]; POS2 = POS2[1:count] ; lambda = apply(Z2PRE,1,median)/0.455; lambda <- floor(1000*lambda+0.5)/1000; names <- NULL
#----- writing the variance inflation factors (lambda) to adjust all Z2 associations for each PC-traits ---------------#
#outfile <- paste(results_folder,"/PLEIOVAR_variance_inflation_factors_nostd",sep="")
outfile <- paste(outfolder,"/PLEIOVAR_variance_inflation_factors_nostd",sep="")
write.table(lambda,outfile,col.names=FALSE,quote=FALSE)

#----- Now we will take each row of Z2PRE and divide it the the corresponding VIF (lambda) which in essence is correcting  each Z2 score of Z2PRE for variance inflation   
Z2PRE_adj = Z2PRE
for(j in 1:ncolT) 
{ 
    Z2PRE_adj[j,] = Z2PRE[j,]/lambda[j]
} 
Z2PRE_adj = floor(100*Z2PRE_adj+0.5)/100

#----- Now we will recover the Z2B for each gene and corrected for variance inflation by parsing Z2PRE_adj using the stored position POS1 and POS2 for the gene ------------------#
#----- Then we will use all the corrected Z2B scores to calculate SSQ, get the DF and generate its corresponding p-value ---------------------------------------------------------#
Z2_POS_MAT  <- rep(0,4*count) ; dim(Z2_POS_MAT) <- c(count,4) 
for(i in 1:count)
{
   TEMP = Z2PRE_adj[ , c(POS1[i]:POS2[i]) ]; TEMP = as.matrix(TEMP) 
   SSQ <- sum(TEMP) ; SSQ <- floor(100*SSQ+0.5)/100 ; npct <- dim(TEMP)[1]; npcs <- dim(TEMP)[2]
   DF <- npct*npcs  ;   pvalue <- pchisq(SSQ,DF,lower.tail=FALSE) ; pvalue <- formatC(pvalue,digits = 2,format="E") ; OUTPUT[i] <- paste(OUTPUT[i],SSQ,DF,pvalue,sep=" ")  
   #----- saving into Z2_POS_MAT -----#
   Z2_POS_MAT[i,1] = POS1[i] ; Z2_POS_MAT[i,2] = POS2[i] ; Z2_POS_MAT[i,3] = SSQ ; Z2_POS_MAT[i,4] = DF 
}    


TEMP = cbind(GENE,Z2_POS_MAT) ; TEMP = as.data.frame(TEMP); colnames(TEMP)[2:5] = c("POS1","POS2","SSQ","DF") 
#---- write gene, pos1, pos2, SSQ and DF to a file ----------#
#write.table(TEMP,"Z2_POSITIONS_BT_GENE.dat",col.names=TRUE,row.names=FALSE,quote=FALSE)
write.table(TEMP,paste0(outfolder, "/Z2_POSITIONS_BY_GENE.dat"),col.names=TRUE,row.names=FALSE,quote=FALSE)
#---- write Z2PRE_adj into a file ---------------------------#
TEMP = as.data.frame(t(Z2PRE_adj)); rownames(TEMP) = NULL ; write.table(TEMP,paste0(outfolder, "/Z2_adjusted.dat"),col.names=TRUE,row.names=FALSE)


#---- nsame of the OUTPUT variable (thousands of rows with resutls for each gene --------#
colnames(OUTPUT) <- paste("chr","ini","fin","gene","begpos","endpos","select_trait","select_snp","SSQ","DF","p-value",sep=" ")
write.table(OUTPUT,paste(outfolder,"/results-file",sep=""),quote=FALSE,col.names=FALSE,row.names=FALSE)


#**********************************************************************#
#**********************************************************************#
#------------------------ PROGRAM ENDS HERE ---------------------------# 
#**********************************************************************#
#**********************************************************************#


