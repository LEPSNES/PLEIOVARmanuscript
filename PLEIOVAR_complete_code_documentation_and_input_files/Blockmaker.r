
 #------ break gene file into blocks -----#

A <- read.table("CONFIG",header=FALSE) 
A <- as.list(A)

outfolder <- A$V1[[2]]
N <- A$V1[[5]]
N <- as.character(N); N <- as.numeric(N)

 B <- read.table("All_Genes_hg19",header=F)
 cmd <- paste("mkdir -p ",outfolder,"/Regions",sep="")
 system(cmd)

 nrow <- dim(B)[1]

 max_blk <- floor(nrow/(N+0.000001))
 left <- nrow %% max_blk

 pos1 <- 0
 for(i in 1:N)
 {
   pos0 <- pos1 + 1
   temp <- 1*(left > 0)*(i <= left)
   pos1 <- pos0 + max_blk -1 + temp
   block <- B[pos0:pos1,]
   outfile <- paste(outfolder,"/Regions/block",i,sep="")
   write.table(block,outfile,row.names=F,col.names=F,quote=F)
 }

 #---- Now that we got the files, our goal is to write L commands to the file to be run by runjobs.pl ----#

 C <- rep(0,N)
 for (i in 1:N)
 {
   C[i] <- paste("R --slave --no-save --no-restore --no-environ --silent --args ","block",i," < Assemble.r",sep="")
 }

write.table(C,"jobs_assemble.lst",col.names=F,row.names=F,quote=F)

