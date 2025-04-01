#!/usr/bin/env perl

use strict;
use warnings;

if(scalar(@ARGV) != 7){
	print "ERROR\tWrong number of arguments for 04_create_obiwizard_config.pl\n";
	exit 1;
}
my ($config_file,$input_assign,$identifier,$threads,$db_path_to_file,$db_path_new,$new_taxdump) = @ARGV;

open (IN,$config_file) or die "Could not open " . $config_file . " for reading!\n";

while (my $line = <IN>){
	chomp $line;
	my @parameter = split(/=/,$line);
	if (scalar(@parameter > 1)){
		if($parameter[0] eq "input_assign"){
			$parameter[1] = "\"$input_assign\"";
		}
		if($parameter[0] eq "identifier"){
			$parameter[1] = "\"$identifier\"";
		}
		if($parameter[0] eq "threads"){
			$parameter[1] = "$threads";
		}
		if($parameter[0] eq "db_path_to_file"){
			$parameter[1] = "\"$db_path_to_file\"";
		}
		if($parameter[0] eq "db_path_new"){
			$parameter[1] = "\"$db_path_new\"";
		}
		if($parameter[0] eq "tax_path_to_file"){
			$parameter[1] = "\"$new_taxdump\"";
		}
		$line = join("=",@parameter);
	}
	print $line . "\n";
}

close IN;

exit 0;
