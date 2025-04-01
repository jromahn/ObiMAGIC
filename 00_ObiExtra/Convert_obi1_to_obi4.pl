#!/usr/bin/nev perl

use strict;
use warnings;
use utf8;  # Ensure the script itself is treated as UTF-8
use open ':std', ':encoding(UTF-8)';  # Ensure input and output are treated as UTF-8
use Text::Unidecode; # transform special letters into ascii friendly ncbi readable output


###################################################
# Pipeline written by Juliane Romahn ( email: romahnjuliane@gmail.com)
# version: 25.01.2024
# aim: change PhyloAlps/ObiTools1 fasta databases to be compatible with ObiTool4 for obitag
##################
# dependencies: perl
###################################################

my $input_file="db_v05_r117.fasta"; #input filename and if necceassary path to file
my $output_file="db_v05_r117__obi4.fasta"; #output filename and if necceassary path to file

####################################################################################################################################################
##################################################  PLEASE DO NOT CHANGE ANYTHING FROM HERE ON #####################################################
####################################################################################################################################################
## transform >RSZ-RSZAXPI000687-79 family_name=Pinaceae; species_name=Abies alba; family=3318; reverse_match=CCATTGAGTCTCTGCACCTATC; taxid=45372; rank=species; forward_error=1; forward_tm=48.78; genus_name=Abies; seq_length_ori=17228; forward_match=GGGTAATCCTGAGCCAA; reverse_tm=56.99; genus=3319; reverse_error=0; species=45372; strand=D;  PHA000002 assembled
### need to extract: accesion species_name & species and the last part

open(my $FH, "<", $input_file) or die "Can not open $input_file: $! \n";
open(my $FH_new, ">", $output_file) or die "Can not create $output_file: $! \n";

while( defined (my $line = <$FH> )){
    chomp $line;

    
    ## extract important information from file
    if( $line =~ /^>/){
        (my $accession= $line) =~s/ .*//;
        $line=~ /(species_name|species_sn)=([\[\]\p{L}\-\.0-9_ #\(\)',\/&\:]+);/; ### includes also weird letters 
        my $speciesname=$2;
        unless($speciesname){ 
            print $line, "\n"; 
            exit
        }
        $line=~ /species=(\d*);/;
        my $taxid=$1;
        $line=~ /; ([\[\]\p{L}\-\.0-9_ #\(\)',\/&\:]*)$/;
        my $last =$1;
        $line ="$accession {\"species_name\":\"$speciesname\",\"taxid\":$taxid} $last";
        #print $line, "\n";
        #exit;
    }

    print $FH_new $line, "\n";
    #exit;
}