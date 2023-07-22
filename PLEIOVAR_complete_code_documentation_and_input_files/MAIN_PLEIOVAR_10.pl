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


#--------------------- running the association ----------------#
my $command = "R CMD BATCH --slave --no-restore ./Z2_scoring.r";
system($command);
#
