#!/usr/bin/perl -w
#use warnings;
#use strict;

#************************************************************
#*********** average duration of process per gene region ****
#------- assemble takes 10sec per region -------------------#
#------- PC-SNP takes 16 sec per region --------------------#
#------- fast PC-SNP takes 4 sec per region ----------------#
#------- statistical analysis takes 0.25 sec per region ----#
#************************************************************
#************************************************************


open(MAIN_INPUT,"CONFIG") || die "Unable to create CONFIG file\n";
my $infolder  = <MAIN_INPUT>;
my $outfolder = <MAIN_INPUT>;chomp($outfolder);
close MAIN_INPUT;

#---------- creation of main folders and generation of sorted (by p-value) results file -----------#

# $command = "sort -k13 -g $outfolder/results-file > sorted_results-file";
# print "$command\n";
# system($command);

$command = "sort -k11 -g $outfolder/results-file > $outfolder/sorted_results-file";
print "$command\n";
system($command);


