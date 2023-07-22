
#---------------- goal is to save all PLEIOVAR results into Z2 Tables, one for each gene ----------------#
#------- we start by opening the config file and capturing the output folder for the association --------#
#--------------------------------------------------------------------------------------------------------#  


 A <- read.table("CONFIG",header=FALSE); A <- as.list(A)
 #---- outputfolder is the top folder where PC-SNP fiels are located -----#
 outfolder <- A$V1[[2]]
 #---- results_folder is where PC-Traits.txt is located and it is where most of the output for the association will be saved ------------------#
 #results_folder = A$V1[[11]]

 cmd = paste("mkdir -p ",outfolder,"/Z2_RDS",sep=""); system(cmd)
 cmd = paste("mkdir -p ",outfolder,"/Z2",sep=""); system(cmd)

 #---- next, we start by defining the files  
 #infile_pos =  paste(results_folder,"/Z2_POSITIONS_BY_GENE.dat",sep=""); infile_Z2  =  paste(results_folder,"/Z2_adjusted.dat",sep="")
 infile_pos =  paste(outfolder,"/Z2_POSITIONS_BY_GENE.dat",sep=""); infile_Z2  =  paste(outfolder,"/Z2_adjusted.dat",sep="")

 infile_pos
 infile_Z2 

 #---- next, we start by opening the positions file and the Z2_adj file --------#
 temp1 = read.table(infile_pos,header=TRUE); temp2 = read.table(infile_Z2,header=TRUE); nr = nrow(temp1) 
 
 for(i in 1:nr)
 {
     GENE = temp1[i, 1]
     POS1 = temp1[i,2]; POS2 = temp1[i,3]
     Z2   = t(temp2[POS1:POS2,]); nc = ncol(Z2)
     cnames = rep("A",nc)
     for(j in 1:nc)
     { 
        cnames[j]=  paste("PC-SNP.",j,sep="") 
     } 
     colnames(Z2) = cnames
     #fil_A <- paste(outfolder,"/Z2_RDS/Z2-file_",GENE[i],sep="");  Z2 ;   saveRDS(Z2B,file=fil_A)
     #fil_B <- paste(outfolder,"/Z2/Z2-file_",GENE[i],sep="") ;  write.table(Z2,fil_B,row.names=TRUE,col.names=TRUE,quote=FALSE)
     fil_A <- paste(outfolder,"/Z2_RDS/Z2-file_",GENE,sep="");  Z2 ;   saveRDS(Z2,file=fil_A)
     fil_B <- paste(outfolder,"/Z2/Z2-file_",GENE,sep="") ;  write.table(Z2,fil_B,row.names=TRUE,col.names=TRUE,quote=FALSE)
 } 
 

#------------------ testing for a given gene   -------------------#

#temp1 = read.table("Z2_POSITIONS_BT_GENE.dat",header=TRUE)
#temp2 = read.table("Z2_adjusted.dat",header=TRUE)
#mygene = "CETP" ; POS = which(temp1[,1] == mygene)
#tp1 = temp1[POS,2] ; tp2 = temp1[POS,3];Z2TEST = t(temp2[tp1:tp2,]); nc = ncol(Z2TEST)
#TEMP = NULL; for(i in 1:nc){ TEMP[i] = paste("PC-SNP_",i,sep="") } ; colnames(Z2TEST) = TEMP
#Z2TEST
 




