#!/usr/bin/env perl

#### Script help analysing via sequences couldn't be demultiplexed via ObiTools4
#### execute via: perl Perl__reasons_unidentified.pl /path/to/file/filname.fastq.gz

use strict;
use warnings;
use IO::Zlib;  #requires Compress::Zlib, of course


## pattern can change sometimes
my $error_message_pattern="obimultiplex_error";

################################## DO NOT CHANGE FROM HERE ON ##################################
my $file   = $ARGV[0];
$file=~ /([\w\-\.]+)\/([\w\-]+)[a-zA-Z\.]+$/;
my $folder=$1;
my $output_file="$folder/$2__reasons.tsv";


open(my $FH_new, ">", $output_file) or die "Can not create $output_file: $! \n";
print {$FH_new} "#Error_message\tForward_Tag\tForward_Tag_Length\tForward_Primer\tForward_Primer_Mimatches\tReverse_Tag\tReverse_Tag_Length\tReverse_Primer\tReverse_Primer_Mimatches\tProportion_Ns\tSeq_length\n";

open(my $FH, "gunzip -c $file |") or die "gunzip $file: $!";
my $counter=0;
while (defined( my $line= <$FH>)) {
    chomp $line;
    
    ## skip non header lines by count
    my $modulo = $counter % 4;
    ++$counter;
    if ($modulo != 0){
        #print $line, "\n";
        next;
    }

    #get sequence name
    $line=~ /^@(.*)\{"/;
    my $name= $1;


    ## extract general error message
    $line=~ /"$error_message_pattern":"([a-zA-Z_ \(\):]+)",/;
    unless( $line =~ /"$error_message_pattern":"([a-zA-Z_ \(\):]+)",/){
        print "No error message \t $error_message_pattern \n", 
              "Did the pattern changed maybe? please check: \n",
              $counter, "\t", $line, "\n";
        exit;
    }
    #next;
    my $error_message = $1;
    if ($error_message eq $folder){
        print "Equal folder \n";
        print $line, "\n"; exit;
    }
    if ($error_message eq $name){
        print "Equal name \n";
        print $line, "\n"; exit;
    }
    unless($error_message){
        print $line, "\n"; exit;
    }
    $error_message=~ s/\(.*\)//;


    my ($demulti_error, $forward_primer, $forward_tag,$forward_error, $forward_tag_length,$reverse_primer, $reverse_tag, $reverse_error,$reverse_tag_length,$ali_length);
    if($line=~ /"\w+forward_error":(\d+),"/){   $forward_error = $1   } else{$forward_error = "NA"  }
    if($line=~ /"\w+forward_primer":"(\w+)",/){  $forward_primer = $1      } else{$forward_primer = "NA"  }
    if($line=~ /"\w+forward_tag":"(\w+)"/){  $forward_tag = $1;  $forward_tag_length = length($forward_tag) } else{$forward_tag = 0 ;  $forward_tag_length = 0   }
    #
    if($line=~ /"\w+reverse_error":(\d+),"/){ $reverse_error = $1       } else{$reverse_error = "NA"  }
    if($line=~ /"\w+reverse_primer":"(\w+)","/){  $reverse_primer = $1      } else{$reverse_primer = "NA"  }
    if($line=~ /"\w+reverse_tag":"(\w+)","/){ $reverse_tag = $1 ;  $reverse_tag_length = length($reverse_tag) } else{$reverse_tag = 0 ;  $reverse_tag_length = 0   }
    #
    # calculate sequence length
    ++$counter;
    $line= <$FH>; chomp $line;$ali_length = length($line) ;  

    # calculate proportion of n in primer pair
    my $primer_pair_seq="$forward_primer$forward_tag$reverse_primer$reverse_tag";
    my $number= $primer_pair_seq =~s/n/n/g;
    my$n_prop=  sprintf("%.2f", $number/$ali_length);
    print {$FH_new} "$error_message\t$forward_tag\t$forward_tag_length\t$forward_primer\t$forward_error\t$reverse_tag\t$reverse_tag_length\t$reverse_primer\t$reverse_error\t$n_prop\t$ali_length \n";
}

# execute R script to visualize the results
system("Rscript 02_diagnosis_unidentified.R $output_file")