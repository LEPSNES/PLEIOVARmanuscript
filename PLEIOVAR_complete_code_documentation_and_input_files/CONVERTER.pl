#!/usr/bin/perl -w
use warnings;
use strict;

#  ------ PROGRAM TO CONVERT vmf files to smaller files ------ # 

die "Usage: perl CONVERTER.pl [Chr]\n" unless ($#ARGV == 0);
my ($chromossome) = @ARGV;

#my $chromossome = 22; 

#--------- infile needs to read-pipeline the input directory -----------------#
#--------- number of rounding of dosages digits peipeline specified ----------#
#--------- outfile needs to resad-pipeline the output directory --------------#

open(MAIN_INPUT,"CONFIG") || die "Unable to create CONFIG file\n";

my $infolder  = <MAIN_INPUT>; chomp($infolder);
my $outfolder = <MAIN_INPUT>; chomp($outfolder);
my $rnd = <MAIN_INPUT>; chomp($rnd);
my $power = 10**$rnd;
    
print"\n $infolder \n";
print"\n $outfolder \n";
print"\n rounding to $rnd decimals\n";

my $infile = "$infolder/Sardinia.b37.ss2120.FAref.impv4.chr".$chromossome.".vcf";
my $outfile = "$outfolder/SAMPLE_".$chromossome; 

print"\n $infile \n";
print"\n $outfile \n\n\n";

open(IN, $infile) || die "Unable to create $infile\n";
open(OUT,">$outfile") || die "Unable to create $outfile\n";

my $count = 0;

while(my $line = <IN>) 
{
    $count++;
    if ($count > 6)
    {
       chomp($line);
       my(@MYLINE) = split(/\t/,$line);
       my $numfields = @MYLINE;
       if ($count == 7)
       { 
           my @HEADER = '';
           $HEADER[0]="CHROM";
           $HEADER[1] = "POS";
           for (my $i = 9; $i <= $numfields-1; $i++)
           {
              $HEADER[$i-7] = $MYLINE[$i];
           }  
           print OUT ("@HEADER\n"); 
       }
       if ($count > 7)
       {
           my @LINE = '';
           $LINE[0]=$MYLINE[0];
           $LINE[1]=$MYLINE[1];
           if ($MYLINE[6] eq "PASS")
           {
              my $flagprint = 1;
              for (my $i = 9; $i <= $numfields-1; $i++)
              {
                 my $temp = $MYLINE[$i];
                 my $pos = index($temp,":");
                 $temp = substr($temp,$pos+1,4); 
                 if ($temp eq ".")
                 { 
                   $flagprint = 0;
                 }
                 else 
                 { 
                   #$LINE[$i-7] = int($temp*(10**$rnd)+0.5)/(10**$rnd);
                   $LINE[$i-7] = int($temp*$power+0.5)/$power;
                 }
              }
              if ($flagprint == 1)
              {
                print OUT ("@LINE\n");
              }
           }
       } 
    }      
}

close MAIN_INPUT;
close IN;
close OUT;

