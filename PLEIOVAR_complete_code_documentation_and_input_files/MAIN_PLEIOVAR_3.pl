#!/usr/bin/perl -w
#use warnings;
#use strict;

#------ generating the cut file by cutting the first two columns ofthe vcf file -------#

my $command = "runjobs.pl jobs_cut_generator.lst"; 
system($command);

