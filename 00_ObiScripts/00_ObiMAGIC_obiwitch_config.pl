#!/usr/bin/env perl

use strict;
use warnings;

if(scalar(@ARGV) != 7){
	print "ERROR\tWrong number of arguments for create_obiwitch_config.pl\n";
	exit 1;
}
my ($config_file,$project,$read1,$read2,$ngs_file,$threads,$pipeline_path) = @ARGV;

open (IN,$config_file) or die "Could not open " . $config_file . " for reading!\n";

while (my $line = <IN>){
	chomp $line;
	my @parameter = split(/=/,$line);
	if (scalar(@parameter > 1)){
		if($parameter[0] eq "project"){
			$parameter[1] = "\"$project\"";
		}
		if($parameter[0] eq "read1"){
			$parameter[1] = "\"$read1\"";
		}
		if($parameter[0] eq "read2"){
			$parameter[1] = "\"$read2\"";
		}
		if($parameter[0] eq "ngs_file"){
			$parameter[1] = "\"$ngs_file\"";
		}
		if($parameter[0] eq "threads"){
			$parameter[1] = "$threads";
		}
		if($parameter[0] eq "FLAG_DEMULTI_FIRST"){
			if($pipeline_path eq "1"){
				$parameter[1] = "FALSE";
			}
			else{
				$parameter[1] = "TRUE";
			}
		}
		$line = join("=",@parameter);
	}
	print $line . "\n";
}

close IN;

exit 0;
