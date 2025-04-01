#!/usr/bin/env perl

#### Script help analysing via sequences couldn't be demultiplexed via ObiTools4
#### execute via: perl 03_Perl_split_demultiplexed_files.pl /path/to/file/filname.fastq.gz /path/to/file/ngsfilname.txt

use strict;
use warnings;
use IO::Zlib;  #requires Compress::Zlib, of course
use List::MoreUtils qw(uniq);

my $file   =  $ARGV[0]; # "eDNA_traps_MP1_08_results/7_eDNA_traps_MP1_08_final.fasta"; 
my $ngs_file=  $ARGV[1]; # "eDNA_traps_MP1_08_results/eDNA_traps_MP1_08__ngsfile.tsv"; # 

### prepare output
my @file= split("/",$file );
my $filename= $file[-1];
(my $prefix=$filename) =~ s/\.\w+$//;
$prefix=~ s/^(\d+)_//;
my $prefix_no = $1 +1;
my $formatted = sprintf("%02d", $prefix_no);
$prefix= $formatted . "_". $prefix;
(my $path= $file) =~s/$filename//;
print $path, "\n";

### loop through ngs and get primer names and save them into list
my @a_primer;
open(my $FH_ngs, "<", $ngs_file) or die "Can not open $ngs_file: $! \n";
while (defined( my $line= <$FH_ngs>)) {
    chomp $line;
    next if $line =~ /^#/;
    my @a_line= split("\t", $line);
    $a_line[1]=~ /^(.*)__(.*)$/;
    push(@a_primer, $1);
}
close $FH_ngs;
@a_primer= uniq @a_primer;

# loop through fasta file for each primer identifier in fasta file
for( my $i=0; $i < scalar(@a_primer); ++$i){
    my $primer= $a_primer[$i];
    my $primer_filename = $path. $prefix. "__". $primer.  ".fasta";

    #open fasta
    open(my $FH, "<", $file) or die "Can not open $filename: $! \n";

    #create new primer file 
    open(my $FH_new, ">", $primer_filename) or die "Can not create $primer_filename: $! \n";

    my $flag_p=0;
    while (defined( my $line= <$FH>)) {
        chomp $line;
        #check and store in flag if fasta header contains primer idenifier
        if( $line =~ />/){
            if ($line =~ /$primer/){
                $flag_p=1;
            }else{
                $flag_p=0;
            }
        }
        #print header and seq as long as flag is TRUE(1)
        if($flag_p){
            print{ $FH_new} $line , "\n";
        }
    }
    close $FH_new;
    close $FH;
}
