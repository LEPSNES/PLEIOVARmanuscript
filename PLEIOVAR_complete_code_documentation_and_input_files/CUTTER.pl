#!/usr/bin/perl -w
use warnings;
use strict;

#  ------ PROGRAM TO CONVERT vmf files to smaller files ------ # 

die "Usage: perl CUTTER.pl [Chr]\n" unless ($#ARGV == 0);

open(MAIN_INPUT,"CONFIG") || die "Unable to create CONFIG file\n";
my $infolder  = <MAIN_INPUT>;
my $outfolder = <MAIN_INPUT>;chomp($outfolder);
my $precision = <MAIN_INPUT>;chomp($outfolder);
my $microsize =  <MAIN_INPUT>;chomp($outfolder);
close MAIN_INPUT;



my $command = "mkdir -p ./MICROFILES";
system($command);
my ($chr) = @ARGV;
my $infile = "$outfolder/SAMPLE_"."$chr";

print"\n /$infile \n";

open(IN, "$infile") || die "Unable to create $infile\n";

my $header = <IN>;
my $count = 0;
my $oldchunk = 0;
my $chunk = 0;

while(my $line = <IN>) 
{
    $count++;
    $chunk = int($count/($microsize+0.00000001))+1;
    if ($chunk > $oldchunk) 
    {
       my $outfile = "$outfolder/MICROFILES/SAMPLE_"."$chr"."_$chunk";
       if ($count > 1) 
       {
         close(OUT);
       } 
       open(OUT,">$outfile") || die "Unable to create $outfile\n";
       print OUT "$header";
    } 
    print OUT "$line"; 
    $oldchunk = $chunk;
}
close(OUT);

