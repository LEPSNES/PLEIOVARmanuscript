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


#open(MAIN_INPUT,"CONFIG") || die "Unable to create CONFIG file\n";
#my $infolder  = <MAIN_INPUT>;
#my $outfolder = <MAIN_INPUT>;chomp($outfolder);
#close MAIN_INPUT;


my $command = "R CMD BATCH --slave --no-restore ./Blockmaker.r";
system($command);

my $command = "runjobs.pl jobs_assemble.lst";
system($command);

