

  A <- read.table("CONFIG",header=FALSE); A <- as.list(A)

  outfolder <- A$V1[[2]]

  # variance based cutoff for PC-Traits
  OPF_trait <- A$V1[[9]]; OPF_trait <- as.character(OPF_trait)   ; OPF_trait <- as.numeric(OPF_trait)

  phenofile <- read.table("mainpheno.dat",header=T)
  colnames(phenofile)[1] <- "ID"; ID <- phenofile[,1]
  #--- take out the ID's from traits file ----#
  TRAITS <-  phenofile[ , c(-1)]

  #temp = TRAITS[,1] ; TRAITS[,1] = TRAITS[,2] ; TRAITS[,2] = temp
  #colnames(TRAITS)[1] = "LDL" ;  #colnames(TRAITS)[2] = "HDL" ;  #colnames(TRAITS)[3] = "Tri"


  #--- get number of columns and rows from trait file without the ID ----#
  ncolT <- ncol(phenofile)-1 ; nrowT <- nrow(phenofile)  

  #----- if only one trait, then standardize and assign 1 to its loadings -----#
  if (ncolT == 1)
  {
     M <- mean(phenofile[,2])
     S <- sd(phenofile[,2])
     OTrait <- (phenofile[,2]-M)/S 
     loadings = "PC-Trait-1 1.00"
     #write.table(loadings,"PC-Trait_loadings",row.names=FALSE,col.names=FALSE)
     FLAG = 1
  } else {

    #--- Generating PC-Traits ----------------####
    prW <- prcomp(TRAITS);  OTrait <- prW$x; loadings <- prW$rotation
    #----- write PC-Trait loadings to file ------#
    write.table(loadings,"PC-Trait_loadings")
    #--- Get % variance from each PC-TRAIT  ---#
    OPF_s <- apply(OTrait,2,var);  OPF_s <- OPF_s/sum(OPF_s);
    #--- defining the vector of cummulative % variance for PC-TRAITS ---#
    OPFcum_s <- matrix(0,1,ncolT); CUM <- 0;
    for(i in 1:ncolT) {  CUM <- CUM + OPF_s[i]; OPFcum_s[i] <- CUM }
    #----- select PC-traits that meet variance explained cutoff -----#
    tflag <- 1*(OPFcum_s >= OPF_trait); select_trait <- min((1:ncolT)[tflag>0]); OTrait <- OTrait[,1:select_trait];
    FLAG = 2
  }

  #---- roundoff PC-Traits to neighrest hundredths -----#
  OTrait <- (OTrait > 0)*floor(100*OTrait+0.5)/100 +(OTrait < 0)*floor(100*OTrait-0.5)/100
  #---- join the ID's with the PC-Traits -----#
  ID <- as.matrix(ID);  OTraitB <- cbind(ID,OTrait)
  path <- paste(outfolder,"/PC-Traits.txt",sep="")
  write.table(OTraitB,path,col.names = FALSE , row.names = FALSE, quote = FALSE)
  path <- paste(outfolder,"/PC-Traits_loadings.txt",sep="")
  #----- write PC-Trait loadings to file. Here we differentiated between having 1 or more PC-traits saving rownames for the single PC-trait ------#
  #---- although I am not sure this differentiation will be needed when using more powerfull software tools -----#
  if (FLAG == 2)
  {
    write.table(loadings,path,col.names=FALSE,row.names=TRUE,quote=FALSE)
  }
  if( FLAG == 1)
  {
    write.table(loadings,path,col.names=FALSE,row.names=FALSE,quote=FALSE)
  }
 

