#!/usr/bin/env perl


#################### Libraries ################
use strict;
use warnings;
use Text::Iconv;
my $converter = Text::Iconv->new("utf-8", "windows-1251");
use Spreadsheet::XLSX;


###################################################
# Pipeline written by Juliane Romahn 
# version: 05.09.2025
# aim: converting apscale results to assign data with ObiTools4 via ObiWizard
##################
###################################################

###################################################
#### manage ARGV

# Help message
my $usage = "
Usage: Convert_obi1_to_obi4.pl <input_file> 
Purpose: Converting APSCALE output into ObiTools VErsion 4 format, in the same path


Arguments:
  input_file\t\tExcel table of APSCALE (full path)

Options:
  -h, --help\t\tShow this help message and exit

Example:
  Convert_obi1_to_obi4.pl PROJECT_apscale/11_read_table/data/0_PROJECT_sequence_read_table_part_0.xlsx
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

#################### INPUT #########################
#### PATH to apscale results, shoul fit the following pattern: {PROJECT}_apscale/11_read_table/data/0_{PROJECT}_sequence_read_table_part_0.xlsx
my $file_excel="02_sedaDNA_apscale/11_read_table/data/0_02_sedaDNA_sequence_read_table_part_0.xlsx";
###################################################


####################################################################################################################################################
##################################################  PLEASE DO NOT CHANGE ANYTHING FROM HERE ON #####################################################
####################################################################################################################################################

(my $file_fasta=$file_excel)=~s/\.(\w+)$/_obitools.fasta/;

## open files
my $excel = Spreadsheet::XLSX->new($file_excel, $converter);
open(my $FH_new, ">", $file_fasta) or die "Can not create $file_fasta: $!\n";



### zero - asv name, one - sequence, rest is samples
my @a_samples; my $flag_first=1;


my $sheet = ${$excel->{Worksheet}}[0];
$sheet->{MaxRow} ||= $sheet->{MinRow};

foreach my $row ($sheet->{MinRow} .. $sheet->{MaxRow}) {    
    $sheet->{MaxCol} ||= $sheet->{MinCol};
    my $counter=0;
    my $samples="";


    #write into fasta
    my $cell_asv = $sheet->{Cells}[$row][0];
    my $cell_Seq = $sheet->{Cells}[$row][1];
    

    #### loop through counts
    foreach my $col ($sheet->{MinCol} ..  $sheet->{MaxCol}) {
        my $cell = $sheet->{Cells}[$row][$col];

        #save replicates if first row identified by a defined flag
        if ( $flag_first && $col > 1 ) {
            push @a_samples, $cell->{Val};
            $flag_first = 0 if $row > 0;
            next;
        }
        #else count reads per asv and combine the samples
        my $i= $col-2;
        if($col>1){ # skip asv name and sequence
            if($cell->{Val}>0){ # check if read number is bigger than 0
                if($counter>0){$samples.=","} # add ; to sample if something was stored already to it
                $samples.="\"".$a_samples[$i]."\"".":". $cell->{Val};
                $counter+=$cell->{Val};
            }
        }
        
    }
    my $extension = " {\"count\":$counter,\"merged_sample\":{$samples}}";

    #print "$extension \n";
    
    
    #only print into fasta if not first row
    unless($flag_first){
        print {$FH_new} ">",$cell_asv->{Val},$extension, "\n", $cell_Seq->{Val}, "\n"; 
    }
    
     
}

#print "@a_samples";

exit 0;
