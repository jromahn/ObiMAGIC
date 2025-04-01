#!/usr/bin/env perl

use strict;
use warnings;
use IPC::Cmd qw[can_run run];

#usage: ngs_rm-cols_placeholder.pl [-primer <primer_name>] <ngs-file1> <ngs-file2> ...
#this script will create an adjusted ngs file for obimagic and print to STDOUT.
#Adjustments include:
#- allow multiple input colum delimiters [,;\t]
#- removing windows line breaks
#- substitution of "I" nucleotides in primer sequences to "N"
#- show (but not exclude) if a tag pair is duplicated in the file and exit 1
#- show (but not exclude) if a Primer_Sample_ID is duplicated in the file and exit 1
#- adding a placeholder for empty samples
#- removing columns 2 and 3 from the original file
#- concatenating multiple ngs files if primers were multiplexed
#TODO: double underscore in Primer_Sample_ID: primer name exists; if not add primer name later (obimagic global config)

my $primer_id = "";

if($ARGV[0] eq "-primer"){
	$primer_id = $ARGV[1];
	splice(@ARGV,0,2);
}

my @header = ("#experiment","Primer_Sample_ID","tags","forward_primer","reverse_primer","extra_information");
print STDOUT join("\t",@header) . "\n";

my $tag_pair_error = 0;
my $primer_sample_id_error = 0;
my $primer_error = 0;
my @delimiter = (
	",",
	";",
	"  ",
	"\\t",
);

for (my $i=0; $i < scalar(@ARGV); $i++){ #process all ngs files specified in @ARGV
	my $input_file = $ARGV[$i];
	my %tag_pair;
	my %primer_sample_id;
	my $placeholder = "8001";
	my $current_delimiter = "";

	open (IN,"$input_file") or die "Could not open " . $input_file . " for reading!\n";

	#check format
	#which delimiter
	my $test_line = `head -n 2 $input_file | tail -n 1`;
	chomp $test_line;
	my $test_error = "";
	my @test;
	foreach(@delimiter){
		if($_ eq "  "){
			@test = split(/ +/,$test_line); #testing for >1 space but splitting with one space if odd number of spaces as delimiter
		}
		else{
			@test = split(/$_/,$test_line);
		}
		my $column_number = scalar(@test);
		$test_error = $test_error . "\"$_\"\t\t$column_number\n";
		if(scalar(@test) >= 6){
			if($_ eq "  "){
				$current_delimiter = " "; #passing one space as delimiter if odd number of spaces as delimiter
			}
			else{
				$current_delimiter = $_;
			}
		}
	}
	if($current_delimiter eq ""){
		print STDERR "ERROR\tDelimiter could not be guessed.\n";
		print STDERR "Delimiter\tTotal number of columns\n";
		print STDERR $test_error;
		exit 1;
	}

	#check format, e.g. which column contains the tag
	my $tag_column_index = "";
	if($current_delimiter eq " "){
		@test = split(/ +/,$test_line);
	}
	else{
		@test = split(/$current_delimiter/,$test_line);
	}
	for(my $i = 0; $i < scalar(@test); $i++){
		if($test[$i] =~ m/:/){
			$tag_column_index = $i;
		}
	}
	if($tag_column_index eq ""){
		print STDERR "ERROR\tTag column could not be guessed.\n";
		exit 1;
	}

	while (my $line = <IN>){
		my $primer_forward_column_index = $tag_column_index + 1;
		my $primer_reverse_column_index = $tag_column_index + 2;
		my $primer_sampleid_column_index = $tag_column_index - 1;
		
		if ($line =~ m/^#/){ #skip header
			next;
		}

		#remove win linebreaks
		$line=~s/\n//g;
		$line=~s/\r//g;
		
		#process columns
		my @cols = ();
		if($current_delimiter eq " "){
			@cols = split(/ +/,$line);
		}
		else{
			@cols = split(/$current_delimiter/,$line);
		}
	
		#substitute I to N in primer sequences
		if($cols[$primer_forward_column_index] =~ m/I/s){
			$cols[$primer_forward_column_index] =~ s/[Ii]/N/g;
		}
		if($cols[$primer_reverse_column_index] =~ m/I/s){
			$cols[$primer_reverse_column_index] =~ s/[Ii]/N/g;
		}
	
		#show warning when tag pair is used twice
		if (exists($tag_pair{$cols[$tag_column_index]})){
			print STDERR "ERROR\tTag pair $cols[$tag_column_index] not unique.\n";
			$tag_pair_error = 1;
		}
		else{
			$tag_pair{$cols[$tag_column_index]} = 1;
		}
		
		#check if all primer sample IDs are unique
		if($tag_column_index >= 2){ #if tags are in column 3 or higher, a Primer_Sample_ID column exists before.
			if($cols[$primer_sampleid_column_index] =~ m/__/){ #check if primer name is set in ngs file
				my @current_primer_sampleid_array = split(/__/,$cols[$primer_sampleid_column_index]);
				if(scalar(@current_primer_sampleid_array) > 2){
					print STDERR "ERROR\tMore than one double underscore in Primer_Sample_ID $cols[$primer_sampleid_column_index].\n";
					$primer_error = 1;
				}
				else{
					my $current_primer_id = $current_primer_sampleid_array[0];
					if($current_primer_id eq ""){
						if($primer_id eq ""){
							print STDERR "ERROR\tNo primer name set and none in Primer_Sample_ID $cols[$primer_sampleid_column_index].\n";
							$primer_error = 1;
						}
						else{
							$cols[$primer_sampleid_column_index] = $primer_id . $cols[$primer_sampleid_column_index];
						}
					}
				}
			}
			else{
				if($primer_id eq ""){
					print STDERR "ERROR\tNo primer name set and none in Primer_Sample_ID $cols[$primer_sampleid_column_index].\n";
					$primer_error = 1;
				}
				else{
					$cols[$primer_sampleid_column_index] = $primer_id . "__" . $cols[$primer_sampleid_column_index];
				}
			}

			if (exists($primer_sample_id{$cols[$primer_sampleid_column_index]})){ #check if all primer sample IDs are unique
				print STDERR "ERROR\tPrimer_Sample_ID $cols[$primer_sampleid_column_index] not unique.\n";
				$primer_sample_id_error = 1;
			}
			else{
				$primer_sample_id{$cols[$primer_sampleid_column_index]} = 1;
			}
		}

		#allow no plate and position info in the last column and add placeholder
		if($cols[-1] eq "F"){
			$cols[-1] = "F @ position=00_00Z;";
		}

		#add placeholder to empty samples
		if(scalar(@cols) > 6){
			my $sampleid_column_index = $tag_column_index - 2;
			if($cols[$sampleid_column_index] eq "" or $cols[$sampleid_column_index] =~ m/^__R/){
				my @extra = split(/=/,$cols[-1]); #extract plate info from extra_information column
				my $plate = $extra[1];
				$plate =~ s/_.*//;
		
				my @psid = split(/__/,$cols[$sampleid_column_index + 1]); #$sampleid_column_index + 1 = Primer_Sample_ID
				$cols[$primer_sampleid_column_index] = $psid[0] . "__TS24MN_" . $placeholder . "_P" . $plate;
				$placeholder++;
			}
		}
	
		#flexible remove columns like e.g. Primer, Sample, which are not allowed in the format processed by obitools
		if(scalar(@cols) > 6){
			my $splice_length = scalar(@cols) - 6;
			splice(@cols,1,$splice_length);
		}
		print STDOUT join("\t",@cols) . "\n";
	}
}

close IN;

if($tag_pair_error == 1){
	print STDERR "ERROR\tAt least one tag pair is not unique.\n";
}
if($primer_sample_id_error == 1){
	print STDERR "ERROR\tAt least one Primer_Sample_ID is not unique.\n";
}
if($primer_error == 1){
	print STDERR "ERROR\tPrimer name could not be set correctly.\n";
}
if($tag_pair_error == 1 or $primer_sample_id_error == 1 or $primer_error == 1){
	exit 1;
}

exit 0;
