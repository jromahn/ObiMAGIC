#1/usr/bin/nev perl

use strict;
use warnings;


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

my $input_file="05_cleaned_EMBL_Euka02_L500_tax.tsv"; #input filename and if necceassary path to file
my $output_file="EMBL_CRABS_database_Aug24.fasta"; #output filename and if necceassary path to file

####################################################################################################################################################
##################################################  PLEASE DO NOT CHANGE ANYTHING FROM HERE ON #####################################################
####################################################################################################################################################


open(my $FH, "<", $input_file) or die "Can not open $input_file: $! \n";
open(my $FH_new, ">", $output_file) or die "Can not create $output_file: $! \n";

while( defined (my $line = <$FH> )){
    chomp $line;
    next if $line =~ /^seqID/;
    
    ## extract important information from file
    my @a_line = split("\t", $line);
    my $tax_ID =$a_line[1] ;
    my $asseccion = $a_line[0];
    my $species_name = $a_line[8];
    $species_name =~ s/"//g;
    (my $scientific_name= $species_name) =~ s/(\w+)_([a-zA-Z.]+).*/$1 $2/;
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