
#---- This programs aims at selecting the SNPS from the SardiNIA vcf file ----#
#---- which correspond to a specific gene region -------#

#---- INPUT: chr begpos endpos gene 
#---- step 1: identify the chunk
#---- step 2: rowbind each chunk  file
#---- step 3: select rows that    

#     Include structure below when calling this program #####

 #---- output file version ( 1 for each block) ------#

 A <- read.table("CONFIG",header=FALSE)
 A <- as.list(A)
 outfolder <- A$V1[[2]]
 extension <- A$V1[[6]]; extension <- as.character(extension) ; extension <- 1000*as.numeric(extension)

 args <- commandArgs()
 block_name <- args[8]
 infile <- paste(outfolder,"/Regions/",block_name,sep="")

 AA <- read.table(infile,header=F)
 nrowA <- dim(AA)[1]; ncolA <- dim(AA)[2];


 MYOUT <- c(outfolder,extension,block_name,infile)
 #write.table(MYOUT,"temp",row.name=FALSE,quote=FALSE)
 
 #------ chromossome number of previous gene region ------#
 oldchr <- 0

 for (i in 1:nrowA)
 {
   chr <- AA[i,1]; ini <- AA[i,2]; fin <- AA[i,3];gene <- AA[i,4]
   begpos <- max(1,ini-extension); endpos <- fin+extension
   Table_file <- paste(outfolder,"/TABLE_CUT_",chr,sep="")
   A  <- read.table(Table_file,header=F)
   D  <- dim(A); nrows <- D[1]; ncols <- D[2]; C <- rep(0,2*nrows)
   dim(C) <- c(nrows,2)

   for(i0 in 1:nrows)
   {
     C[i0,1] <- 1*(A[i0,3] <= begpos)*(A[i0,4] >= begpos)
     C[i0,2] <- 1*(A[i0,3] <= endpos)*(A[i0,4] >= endpos)          
   }

   chunk0 <-  which(C[,1]==1) ; chunk1 <-  which(C[,2]==1)
    
   if (length(chunk0)==0)
   {
     chunk0 <- 0
   }
   if (length(chunk1)==0)
   {
     chunk1 <- 0
   }

   if (chunk0 > 0 && chunk1 > 0)
   {
     FILDATA <- paste(chr,begpos,endpos,gene,"chunk0 -> ",chunk0,"chunk1 -> ",chunk1,sep=" ")
   } else {
     next
   }
  
   myarray <- NULL
   for (ii in chunk0:chunk1)
   {
     filename <- paste(outfolder,"/MICROFILES/SAMPLE_",chr,"_",ii,sep="")
     B <- read.table(filename,header = TRUE)
     flag <- (B[,2] >= begpos)*(B[,2] <= endpos)
     loc <- which(flag==1)
     myarray <- rbind(myarray,B[loc,])
   }
   myarray <- as.matrix(myarray)
   if( dim(myarray)[1] > 0 )
   {
      C <- colnames(myarray); SNPPOS <- myarray[,2]; C <- C[-c(1,2)]; L <- nchar(C)
      C0 <- substr(C,2,L); C0 <- as.numeric(C0)
      myarray  <- myarray[,-c(1,2)]; colnames(myarray) <- NULL; rownames(myarray) <- NULL
      C0 <- as.matrix(C0); X <- cbind(C0,t(myarray)); colnames(X) <- c("ID",SNPPOS)
      FILENAME <- paste(outfolder,"/Assemble/",gene,"_assembled",sep="")
      write.table(X,FILENAME,col.names=TRUE,row.names=FALSE,quote=FALSE)    
   }   
}

