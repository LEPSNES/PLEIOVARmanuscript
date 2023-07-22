#!/usr/bin/perl -w
#use warnings;
#use strict;

#------ generating the cut file by cutting the first two columns ofthe vcf file -------#

die "Usage: perl cut_generator.pl [Chr]\n" unless ($#ARGV == 0);
my ($chr) = @ARGV;

open(MAIN_INPUT,"CONFIG") || die "Unable to create CONFIG file\n";
my $infolder  = <MAIN_INPUT>;
my $outfolder = <MAIN_INPUT>;chomp($outfolder);
close MAIN_INPUT;

my $infile = "CUT_$chr";
my $command = "cut -f1-2 -d' ' $outfolder/SAMPLE_$chr  > $outfolder/$infile";
system($command);

