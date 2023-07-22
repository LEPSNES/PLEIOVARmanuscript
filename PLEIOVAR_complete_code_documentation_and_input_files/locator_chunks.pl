#!/usr/bin/perl -w
use warnings;
use strict;

#  ------ PROGRAM TO CONVERT vmf files to smaller files ------ # 

die "Usage: perl locator_chunks.pl [Chr]\n" unless ($#ARGV == 0);
my ($chr) = @ARGV;

open(MAIN_INPUT,"CONFIG") || die "Unable to create CONFIG file\n";
my $infolder  = <MAIN_INPUT>;
my $outfolder = <MAIN_INPUT>;chomp($outfolder);
my $precision = <MAIN_INPUT>;chomp($outfolder);
my $micro_size = <MAIN_INPUT>;chomp($micro_size);
close MAIN_INPUT;

my $infile = "$outfolder/CUT_"."$chr";
my $outfile = "$outfolder/TABLE_CUT_"."$chr";

print"\n $infile \n";

#--- determining the zise of the file in total lines - header ----#

my $command = "wc -l $infile > a.$chr.out";
system($command);
open(IN0, "a.$chr.out") || die "Unable to create a.$chr.out\n";
my $line = <IN0>;
chomp($line);
my ($numlines,@tmp) = split(/ /,$line);

#---- exclude header from count -----#

$numlines = $numlines - 1;
#my $maxcount = int($numlines/500.0000000001)*500;
my $maxcount = int($numlines/($micro_size+0.0000000001))*$micro_size;
close(IN0);

open(IN, $infile) || die "Unable to create $infile\n";
$line = <IN>;
my $header = $line;

my $chunk = 0;
my $remainder = 0;

my $begpos = 0;
my $endpos = 0;

open(OUT,">$outfile") || die "Unable to create $outfile\n";

for (my $i = 1; $i <= $numlines; $i++) 
{
  $line = <IN>;
  $chunk = int($i/($micro_size+0.0000000001))+1;
  my $maxcount = int($numlines/($micro_size+0.0000000001))*$micro_size;
  $remainder = $i - $micro_size*($chunk-1);
  if ($remainder == 1)
  {
    chomp($line);
    my ($tmp1,$tmp2) = split(/ /,$line);       
    $begpos = $tmp2;
  }
  if ($remainder == $micro_size)
  { 
    chomp($line);
    my ($tmp1,$tmp2) = split(/ /,$line);
    $endpos = $tmp2;
    print OUT "$chr $chunk $begpos $endpos\n";
  }
}
chomp($line);
my ($tmp1,$tmp2) = split(/ /,$line);
$endpos = $tmp2;
print OUT "$chr $chunk $begpos $endpos\n";
close(OUT);
