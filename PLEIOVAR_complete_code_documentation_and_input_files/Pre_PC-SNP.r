

 #---- Now that we got the files, our goal is to write L commands to the file to be run by runjobs.pl ----#

 A <- read.table("CONFIG",header=FALSE)
 A <- as.list(A)

 N <- A$V1[[5]]
 N <- as.character(N); N <- as.numeric(N)


 C <- rep(0,N)
 for (i in 1:N)
 {
   C[i] <- paste("R --slave --no-save --no-restore --no-environ --silent --args ","block",i," < PC-SNP.r",sep="")
 }

 write.table(C,"jobs_PC-SNPs.lst",col.names=F,row.names=F,quote=F)

