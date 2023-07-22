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

#---------- creation of main folders -----------#

my $command = "mkdir -p $outfolder";
print "$command\n";
system($command);

$command = "mkdir -p $outfolder/MICROFILES";
print "$command\n";
system($command);

$command = "mkdir -p $outfolder/Regions";
print "$command\n";
system($command);

$command = "mkdir -p $outfolder/Assemble";
print "$command\n";
system($command);

$command = "mkdir -p $outfolder/PC-SNPS";
print "$command\n";
system($command);


$command = "mkdir -p $outfolder/PC-SNPS_lambdas";
print "$command\n";
system($command);

$command = "mkdir -p $outfolder/PC-SNPS_loadings";
print "$command\n";
system($command);

$command = "mkdir -p $outfolder/PC-SNPS_var_exp";
print "$command\n";
system($command);

$command = "mkdir -p $outfolder/Z2";
print "$command\n";
system($command);

$command = "mkdir -p $outfolder/GWAS-Z";
print "$command\n";
system($command);

$command = "mkdir -p $outfolder/GWAS-Z-adj";
print "$command\n";
system($command);




#---------- extracting vcf files to the specified form (SAMPLE files) ------#


$command = "mkdir -p $outfolder/PC-SNPS_loadings";
print "$command\n";
system($command);

$command = "mkdir -p $outfolder/PC-SNPS_var_exp";
print "$command\n";
system($command);

$command = "mkdir -p $outfolder/Z2";
print "$command\n";
system($command);

$command = "mkdir -p $outfolder/GWAS-Z";
print "$command\n";
system($command);

$command = "mkdir -p $outfolder/GWAS-Z-adj";
print "$command\n";
system($command);


#---------- extracting vcf files to the specified form (SAMPLE files) ------#

$command = "runjobs.pl converter_jobs.lst";
system($command);





