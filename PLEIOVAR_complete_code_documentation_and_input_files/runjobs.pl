#!/usr/bin/perl

use strict;
my $usage = "usage: $0 joblist\n";
my $arg=0;
my $input_file = $ARGV[$arg++] or die $usage;

open(CMDS, $input_file) or die $!;

while(my $jobcommand=<CMDS>){
   chop($jobcommand);
   system("sbatch --partition=TMP --wrap=\"$jobcommand\"");
   #sleep 3600;
}
close CMDS;
exit(0);
