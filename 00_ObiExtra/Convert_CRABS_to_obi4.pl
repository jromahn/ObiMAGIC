#!/usr/bin/env perl

use strict;
use warnings;
use Encode qw(encode);


###################################################
# Pipeline written by Juliane Romahn ( email: romahnjuliane@gmail.com)
# version: 25.01.2024
# aim: change crabs databases to be compatible with ObiTool4 for obitag
####### input file: is a tsv file from CRABS after the  6. dereplicate or 7, clean step (because of memory issues, db needs a insilico_Pcr step)
## important: written for the 10 columns tsv file and the 7 rank system (default): superkingdom, phylum, class, order, family, genus, and species
## input style as tab separated: accession taxID rank_1 rank_2 rank_3 rank_4 rank_5 rank_6 rank_7 sequence
######## output file: fasta file in obitools format
##################
# dependencies: perl
###################################################


###################################################
#### manage ARGV

# Help message
my $usage = "
Usage: Convert_CRABS_to_obi4.pl <input_file> 
Purpose: Converting CRABS tsv table including taxonomy into DB for ObiTools4 in the same path

Arguments:
  input_file\t\tCRABS table with NCBI ID, taxonomy and sequences in tsv format (full path)

Options:
  -h, --help\t\tShow this help message and exit

Example:
  Convert_CRABS_to_obi4.pl CRABS_taxonomy.tsv
";

# Check for help flag
if (@ARGV == 1 && ($ARGV[0] eq "-h" || $ARGV[0] eq "--help")) {
  print $usage;
  exit 0;
}

# Check number of arguments
if (@ARGV != 1) {
  print "ERROR: Wrong number of arguments. Expected 1, got " . scalar(@ARGV) . "\n";
  print $usage;
  exit;
}

# Assign arguments to named variables
my $input_file  = $ARGV[0];
(my $output_file = $input_file) =~ s/\.\w+$/_obi4transformed.fasta/;

# check if inputfile exists 
unless (-e $input_file) {
  print "ERROR: \t Input file not found: $input_file\n";
  exit;
}

############## old
#my $input_file="Euka02_CRABS_EMBL_new//05_cleaned_EMBL_Euka02_L500_tax.tsv"; #input filename and if necceassary path to file
#my $output_file="Euka02_CRABS_EMBL_new/EMBL_CRABS_database_Aug24_Euka02_L500_new.fasta"; #output filename and if necceassary path to file

####################################################################################################################################################
##################################################  PLEASE DO NOT CHANGE ANYTHING FROM HERE ON #####################################################
####################################################################################################################################################


open(my $FH, "<", $input_file) or die "Can not open $input_file: $! \n";
open(my $FH_new, ">", $output_file) or die "Can not create $output_file: $! \n";

while( defined (my $line = <$FH> )){
    chomp $line;
    next if $line =~ /^seqID/;
    #remove any possible problematic special characters
    $line = encode("UTF-8", $line);
    $line =~ s/'//g;
    $line =~ s/&//g;
    $line =~ s/\//_/g;
    $line =~ s/\+/_/g;
    $line =~ s/\-/_/g;
    $line =~ s/\|/_/g;
    $line =~ s/\[//g;
    $line =~ s/\#//g;
    $line =~ s/\]//g;
    
    
    ## extract important information from file
    my @a_line = split("\t", $line);
    my $tax_ID =$a_line[1] ;
    my $asseccion = $a_line[0];
    my $species_name = $a_line[8];
    $species_name =~ s/"//g;
    $species_name=~ s/^cf\._//;
    (my $scientific_name= $species_name) =~ s/([a-zA-Z]+)_([a-zA-Z\.]+).*/$1 $2/;
    my $seq = $a_line[9];

    ## create new output
    #print {$FH_new} ">$asseccion {\"direction\":\"forward\",\"forward_error\":0,",
    #                "\"forward_match\":\"$forward\",\"forward_primer\":\"$forward\",\"reverse_error\":0,",
    #                "\"reverse_match\":\"$reverse\",\"reverse_primer\":\"$reverse\",",
    #               "\"scientific_name\":\"$scientific_name\",\"taxid\":$tax_ID} ",
    #                "$species_name\n";

    print {$FH_new} ">$asseccion {\"species_name\":\"$scientific_name\",\"taxid\":$tax_ID} ",
                    "$species_name\n";
    print {$FH_new} $seq, "\n";
    #exit;
}