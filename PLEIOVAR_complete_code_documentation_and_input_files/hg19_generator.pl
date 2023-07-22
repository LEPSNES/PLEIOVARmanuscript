
#-----------------------------------------------------------------------------------------------------
#---- goal is to get the information from UCSC website and structure it into chr GENE start end -----#
#-----------------------------------------------------------------------------------------------------

  my $command0 = "wget -nc http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz";
  system($command0);
  #my $command0 = "mv refGene.txt.gz $folder/refGene.txt.gz";
  #system($command0);
  $command0 = "gunzip ./refGene.txt.gz.txt";
  system($command0);
  open(IN,"./refGene.txt") ||  die "Unable to read ./refGene.txt\n";
  open(OUT,">./UCSC_gene_info") ||  die "Unable to read ./UCSC_gene_info\n";
  while(my $line = <IN>)
  {
     chomp($line);
     my @fox = split("\t",$line);
     my @fout = ();
     $fout[0] = $fox[2];
     $fout[1] = $fox[3];
     $fout[2] = $fox[4];
     $fout[3] = $fox[5];
     $fout[4] = $fox[12];
     $fout[5] = "\n";
     $, = "\t";
     print OUT @fout;
   }
  close(OUT);
  close(IN);

#---- processing the gene file ------------#
  my $outfolder = "./genefiles";
  my $infile = "UCSC_gene_info";
  my $infolder =  "./";
  $command0 = "mkdir -p $outfolder";
  system($command0);
  my $regfile = "$infolder/$infile";
  open(IN,$regfile) ||  die "Unable to read $regfile\n";
  while(my $line = <IN>)
  {
     chomp($line);
     #---- build the structure chr strnad beg end genename -----
     my ($chr,$strand,$beg,$end,$genename) = split(/\t/, $line);
     my $chrnum = substr($chr,3,2);
     next if(!looks_like_number($chrnum));
     if ( $chrnum >= 1 && $chrnum <= 22)
     {
        open(OUT,">>$outfolder/Rxxx$chrnum");
        print OUT "$chrnum $beg $end $genename\n";
     }
     close(IN);
     close(OUT);
     for (my $i = 1; $i <= 22; $i++)
     {
        my $nfile = "$outfolder/Rxxx$i";
        my $commandA = "sort -k1,1n -k2,2n $nfile | uniq > tmp";
        system($commandA);
        $commandA = "mv tmp $nfile";
        system($commandA);
     }
     my $commandA = "cat $outfolder/Rxxx* > All_Genes";
     system($commandA);
     $commandA = "rm -rf $outfolder";
     system($commandA);
     $commandA = "rm refGene.txt UCSC_gene_info";
     system($commandA);

  }



