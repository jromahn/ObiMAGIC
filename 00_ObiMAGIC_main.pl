#!/usr/bin/env perl

use strict;
use warnings;
use Cwd 'abs_path';
use IPC::Cmd qw[can_run run];
use File::Basename;
use Scalar::Util qw(looks_like_number);

#define variables
my $version = "0.1";
my @orig_ngs = ();
my $out_dir = abs_path("./");
my @orig_fq = ();
my $project = "";
my $threads = 10;
my $slurm = 0;
my $no_slurm = 0;
my $slurm_mem = "10G";
my $slurm_opts = "--cpus-per-task=$threads --mem=$slurm_mem --output %x-%N-%j.log --error %x-%N-%j.err --export=PATH=\$PATH";
my $verbose = 0;
my $keep_tmp = 0;
my $pipeline_path = "1";
my $dry = 0;
my $cmd;
my $home = `echo \$HOME`;
chomp $home;
my $db_config_default = abs_path(dirname($0) . "/config_db.tsv");
my $obiwitch_config_default = abs_path(dirname($0) . "/00_ObiScripts/config_ObiWitch.ini");
my $obiwizard_config_default = abs_path(dirname($0) . "/00_ObiScripts/config_ObiWizard.ini");
my $db_config = "";
my $obiwitch_config = "";
my $obiwizard_config = "";
my %db;
my $primer = "";
my $new_taxdump_default = `grep "tax_path_to_file" $obiwizard_config_default`;
chomp $new_taxdump_default;
$new_taxdump_default =~ s/^tax_path_to_file="//;
$new_taxdump_default =~ s/".*$//;
my $new_taxdump = "";
my @dependencies = (
	#external dependencies
	"Rscript",
	"pigz",
	"obi4_obiannotate",
	"obi4_obipairing",
	"obi4_obigrep",
	"obi4_obimultiplex",
	"obi4_obiuniq",
	"obi4_obicount",
	"obi4_obiannotate",
	"obi4_obiclean",
	"obi4_obipcr",
	"obi4_obirefidx",
	"obi4_obitag",
	"obi4_obidistribute",
	#ObiMAGIC main scripts
	"00_ObiHelp_NGS_control.pl",
	"00_ObiMAGIC_main.pl",
	"01_ObiWitch_diagnostic.R",
	"01_ObiWitch_main.sh",
	"02_ObiHelp_split_samples.pl",
	"03_ObiWizard_diagnostic.R",
	"03_ObiWizard_main.sh",
	#ObiScripts
	"00_OBIMAGIC_functions.R",
	"00_ObiMAGIC_obiwitch_config.pl",
	"00_ObiMAGIC_obiwizard_config.pl",
	"00_ObiMAGIC_obiwizard_jobarray.sh",
	"01_ObiWitch_unidentified_diagnostic.R",
	"01_ObiWitch_unidentified_main.pl",
	#ObiExtra
	"00_ObiHelp_perl_control.pl",
	"00_ObiHelp_R_control.R"
);

my $input_error = 0;

#subroutines
sub print_help{
	print STDOUT "\n";
	print STDOUT "00_ObiMAGIC_main.pl v$version\n";
	print STDOUT "\n";
	print STDOUT "Description:\n";
	print STDOUT "\tThis script will run a pipeline to filter and assign metabarcoding Illumina data with obitools.\n";
	print STDOUT "\tThe pipeline will run these steps:\n";
	print STDOUT "\t\t1) pre-process the specified NGS file for usage in obitools\n";
	print STDOUT "\t\t2) obiwitch filtering\n";
	print STDOUT "\t\t3) demultiplexing\n";
	print STDOUT "\t\t4) obiwizard assignment\n";
	print STDOUT "\t\t5) create tar archive containing relevant results\n";
	print STDOUT "\n";
	print STDOUT "Usage:\n";
	print STDOUT "\t00_ObiMAGIC_main.pl -ngs <NGS file> -fq <paired_1.fq.gz,paired_2.fq.gz> -pipeline-path {1,2}\n";
	print STDOUT "\t                   [-project <custom project name> -db-config <custom database conig>\n";
	print STDOUT "\t                    -o <output directory> -t <number of threads> -slurm-mem <slurm memory>\n";
	print STDOUT "\t                    -slurm-opts \"<option1> [<option2> ...]\"]\n";
	print STDOUT "\n";
	print STDOUT "Mandatory:\n";
	print STDOUT "\tMandatory options -ngs and -fq can be specified multiple times\n";
	print STDOUT "\t-ngs STR\t\tNGS file containing primer, sample ID and tag information\n";
	print STDOUT "\t-fq STR\t\t\tPaths to two fastq files containing paired Illumina reads comma sperated\n";
	print STDOUT "\n";
	print STDOUT "Options: [default]\n";
	print STDOUT "\t-pipeline-path {1,2}\tSpecify the order of merging and sample demultiplexing [$pipeline_path]\n";
	print STDOUT "\t\t\t\t1: Run first merging and second sample demultiplexing\n";
	print STDOUT "\t\t\t\t2: Run first sample demultiplexing and second merging\n";
	print STDOUT "\t-project STR\t\tSpecify a custom project name\n";
	print STDOUT "\t\t\t\t[string before first underscore in first specified NGS file]\n";
	print STDOUT "\t-obiwitch-config STR\tCustom config file for obiwitch\n";
	print STDOUT "\t\t\t\t[$obiwitch_config_default]\n";
	print STDOUT "\t-obiwizard-config STR\tCustom config file for obiwizard\n";
	print STDOUT "\t\t\t\t[$obiwizard_config_default]\n";
	print STDOUT "\t-db-config STR\t\tTSV file specifying primer database combination\n";
	print STDOUT "\t\t\t\t[$db_config_default]\n";
	print STDOUT "\t-new-taxdump STR\tCustom path to NCBI's new_taxdump directory. Overrides default from\n";
	print STDOUT "\t\t\t\t-obiwizard-config [$new_taxdump_default]\n";
	print STDOUT "\t-primer STR\t\tSet a primer name when not present in the NGS file. Does not override\n";
	print STDOUT "\t\t\t\tprimer name in NGS file. Only allowed when one NGS file is specified.\n";
	print STDOUT "\t\t\t\tPossible values have to be defined in config_db.tsv\n";
	print STDOUT "\t-o STR\t\t\tOutput directory [.]\n";
	print STDOUT "\t\t\t\tWill be created if not existing\n";
	print STDOUT "\t-keep-tmp\t\tDo not delete concatenated fastq files [off]\n";
	print STDOUT "\t-t INT\t\t\tNumber of threads [$threads]\n";
	print STDOUT "\t-slurm-mem STR\t\tMemory that will be allocated by sbatch (--mem) [$slurm_mem]\n";
	print STDOUT "\t-slurm-opts STR\t\tOptions that will be passed to sbatch\n";
	print STDOUT "\t\t\t\t--cpus-per-task is set according to -t\n";
	print STDOUT "\t\t\t\t--mem is set according to -slurm-mem\n";
	print STDOUT "\t\t\t\t[$slurm_opts]\n";
	print STDOUT "\t\t\t\tPass options with quotes: -slurm-opts \"<option1> <option2>\"\n";
	print STDOUT "\t-no-slurm\t\tDo not submit internal processes via slurm if sbatch is detected [off]\n";
	print STDOUT "\t-v\t\t\tPrint executed commands to STDERR [off]\n";
	print STDOUT "\t-dry-run\t\tOnly print commands to STDERR instead of executing [off]\n";
	print STDOUT "\n";
	print STDOUT "\t-h or --help\t\tPrint this help and exit\n";
	print STDOUT "\t-version\t\tPrint version number and exit\n";
	exit;
}

sub exe_cmd{
	my ($cmd,$verbose,$dry) = @_;
	if($verbose == 1){
		print STDERR "CMD\t$cmd\n";
	}
	if($dry == 0){
		system("$cmd") == 0 or die "ERROR\tsystem $cmd failed: $?";
	}
}

#print help if no argument is passed
if(scalar(@ARGV==0)){
	print_help;
}

#process @ARGV
for (my $i = 0; $i < scalar(@ARGV);$i++){
	if ($ARGV[$i] eq "-ngs"){
		push(@orig_ngs, abs_path($ARGV[$i+1]));
	}
	if ($ARGV[$i] eq "-fq"){
		push(@orig_fq,$ARGV[$i+1]);
	}
	if ($ARGV[$i] eq "-pipeline-path"){
		$pipeline_path = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "-new-taxdump"){
		$new_taxdump = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "-project"){
		$project = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "-obiwitch-config"){
		$obiwitch_config = abs_path($ARGV[$i+1]);
	}
	if ($ARGV[$i] eq "-obiwizard-config"){
		$obiwizard_config = abs_path($ARGV[$i+1]);
	}
	if ($ARGV[$i] eq "-db-config"){
		$db_config = abs_path($ARGV[$i+1]);
	}
	if ($ARGV[$i] eq "-primer"){
		$primer = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "-t"){
		$threads = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "-o"){
		$out_dir = abs_path($ARGV[$i+1]);
	}
	if ($ARGV[$i] eq "-keep-tmp"){
		$keep_tmp = 1;
	}
	if ($ARGV[$i] eq "-slurm-mem"){
		$slurm_mem = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "-slurm-opts"){
		$slurm_opts = $ARGV[$i+1];
	}
	if ($ARGV[$i] eq "-no-slurm"){
		$no_slurm = 1;
	}
	if ($ARGV[$i] eq "-v"){
		$verbose = 1;
	}
	if ($ARGV[$i] eq "-dry-run"){
		$dry = 1;
		$verbose = 1;
	}
	if ($ARGV[$i] eq "-h" or $ARGV[$i] eq "-help" or $ARGV[$i] eq "--help"){
		print_help;
	}
	if ($ARGV[$i] eq "-version"){
		print STDERR $version . "\n";
		exit;
	}
}

#print the original command
print STDERR "CMD\t" . $0 . " " . join(" ",@ARGV) . "\n";

#catch multiple input errors
if(scalar(@orig_ngs) == 0){
	print STDERR "ERROR\tSpecify -ngs at least once\n";
	$input_error = 1;
}
if(scalar(@orig_ngs) > 1){
	if($primer ne ""){
		print STDERR "ERROR\tSetting a primer name is only possible when a single NGS file is specified.\n";
		$input_error = 1;
	}
}

if(scalar(@orig_fq) == 0){
	print STDERR "ERROR\tSpecify -fq at least once\n";
	$input_error = 1;
}

if($pipeline_path eq ""){
	print STDERR "ERROR\tSpecify either 1 or 2 for -pipeline-path\n";
	$input_error = 1;
}
else{
	if(looks_like_number($pipeline_path) == 1){
		if($pipeline_path != 1 and $pipeline_path != 2){
			print STDERR "ERROR\tSpecify either 1 or 2 for -pipeline-path\n";
			$input_error = 1;
		}
		if($pipeline_path == 2){
			if(not defined(can_run("cutadapt"))){
				print STDERR "ERROR\tcutadapt is not in your \$PATH\n";
				$input_error = 1;
			}
		}
	}
	else{
		print STDERR "ERROR\tSpecify either 1 or 2 for -pipeline-path\n";
		$input_error = 1;
	}
}

foreach(@orig_ngs){
	if(not -f $_){
		print STDERR "ERROR\tNGS file $_ is not a file\n";
		$input_error = 1;
	}
}

#check if custom or default new taxdump exists
if($new_taxdump eq ""){
	if($new_taxdump_default ne ""){
		if(not -d $new_taxdump_default){
			print STDERR "ERROR\tDefault new taxdump directory $new_taxdump_default does not exist!\n";
			$input_error = 1;
		}
		else{
			$new_taxdump = $new_taxdump_default;
		}
	}
	else{
		print STDERR "ERROR\tDefault new taxdump directory not specified or empty!\n";
		$input_error = 1;
	}
}
else{
	if(not -d $new_taxdump){
		print STDERR "ERROR\tNew taxdump directory $new_taxdump does not exist!\n";
		$input_error = 1;
	}
}

if(-f $out_dir){
	print STDERR "ERROR\tOutput directory $out_dir is already a file!\n";
	$input_error = 1;
}

if($threads !~ m/^\d+$/ or $threads < 1){
	print STDERR "ERROR\tThreads is no integer >= 1!\n";
	$input_error = 1;
}

#update slrum options
$slurm_opts = "--cpus-per-task=$threads --mem=$slurm_mem --output %x-%N-%j.log --error %x-%N-%j.err --export=PATH=\$PATH";

#check if dependencies are in $PATH
foreach(@dependencies){
	if(not defined(can_run("$_"))){
		print STDERR "ERROR\t$_ is not in your \$PATH\n";
		$input_error = 1;
	}
}

#test if all needed perl modules of depending perl scripts are available
$cmd = "00_ObiHelp_perl_control.pl > /dev/null 2>&1";
if($verbose == 1){
	print STDERR "CMD\t$cmd\n";
}
if($dry == 0){
	system("$cmd");
	if($? != 0){
		print "ERROR\tsystem $cmd failed: $?\n";
		print "ERROR\tAt least one required perl module is not available!\n";
		$input_error = 1;
	}
}

#test if all needed R libraries are available
$cmd = "00_ObiHelp_R_control.R > /dev/null 2>&1";
if($verbose == 1){
	print STDERR "CMD\t$cmd\n";
}
if($dry == 0 and defined(can_run("Rscript"))){
	system("$cmd");
	if($? != 0){
		print "ERROR\tsystem $cmd failed: $?\n";
		print "ERROR\tAt least one required R library is not available!\n";
		$input_error = 1;
	}
}

if(defined(can_run("sbatch"))){
	if($no_slurm == 0){
		print STDERR "INFO\tSlurm detected. Will submit jobs via sbatch.\n";
		$slurm = 1;
	}
	else{
		print STDERR "INFO\tSlurm detected but -no-slurm set. Will not submit jobs via sbatch.\n";
	}
}

#define obiwitch config
if ($obiwitch_config eq ""){
	if(-f $obiwitch_config_default){
		print STDERR "INFO\tUsing default obiwitch config $obiwitch_config_default\n";
		$obiwitch_config = $obiwitch_config_default;
	}
	else{
		print STDERR "ERROR\tNo obiwitch config file specified and no default obiwitch config file found.\n";
		$input_error = 1;
	}
}

#define obiwizard config
if ($obiwizard_config eq ""){
	if(-f $obiwizard_config_default){
		print STDERR "INFO\tUsing default obiwizard config $obiwizard_config_default\n";
		$obiwizard_config = $obiwizard_config_default;
	}
	else{
		print STDERR "ERROR\tNo obiwizard config file specified and no default obiwitch config file found.\n";
		$input_error = 1;
	}
}

#define DB config
if($db_config eq ""){
	if(-f $db_config_default){
		print STDERR "INFO\tUsing default DB config $db_config_default\n";
		$db_config = $db_config_default;
	}
	else{
		print STDERR "ERROR\tNo default DB config found and no custom DB config specified\n";
		$input_error = 1;
	}
}

#read DB config
if($verbose == 1){
        print STDERR "INFO\tReading DB config $db_config\n";
}

open (DBCONF,'<',$db_config) or die "Could not open $db_config for reading!\n";

while (my $line = <DBCONF>){
	chomp $line;
	my ($primer,$db_path) = split(/\t/,$line);
	$db_path = abs_path(dirname($0)) . "/$db_path";
	if(exists($db{$primer})){
		print STDERR "WARNING\t$primer specified multiple times in $db_config. Using first mention only:\n";
		print STDERR "WARNING\t$primer\t$db{$primer}\n";
	}
	else{
		if(-f $db_path){
			$db{$primer} = $db_path;
		}
	}
}

close DBCONF;

if($primer ne ""){
 	my $primer_match = 0;
	foreach(keys(%db)){
		if($primer eq $_){
 			$primer_match = 1;
		}
	}
	if($primer_match == 0){
		print STDERR "ERROR\tPrimer name $primer not in the list of possible primers: {" . join(",",keys(%db)) . "}\n";
		$input_error = 1;
	}
}

if ($input_error == 1){
	print STDERR "ERROR\tInput error detected!\n";
	exit 1;
}

if(not -d "$out_dir"){
	print STDERR "INFO\tCreating output directory $out_dir\n";
	$cmd = "mkdir -p $out_dir";
	exe_cmd($cmd,$verbose,$dry);
}

#check if fq file pairs exist
my %fq_filter;

foreach(@orig_fq){
	my @pair = split(/,/,$_);
	foreach(@pair){
		if($_ =~ m/^~/){
			$_ =~ s/^~/$home/; #~ is translated by bash into $HOME. This does not work if there is no space infront. That means if the second file starts with "~" it will not be recognized even though it exists
		}
	}
	if(scalar(@pair) != 2){
		print STDERR "INFO\tNot a pair: $_ - skipping file(s)\n";
	}
	else{
		my $file_error = 0;
		if(not -f "$pair[0]"){
			print STDERR "INFO\tNo file $pair[0] - skipping pair $_\n";
			$file_error = 1;
		}
		if(not -f "$pair[1]"){
			print STDERR "INFO\tNo file $pair[1] - skipping pair $_\n";
			$file_error = 1;
		}
		if($file_error == 0){
			if(exists($fq_filter{abs_path($pair[0]) . "," . abs_path($pair[1])})){
				print STDERR "INFO\tPair " . abs_path($pair[0]) . "," . abs_path($pair[1]) . " already specified\n";
			}
			else{
				$fq_filter{abs_path($pair[0]) . "," . abs_path($pair[1])} = 1;
			}
		}
	}
}

#pre-process ngs file
#if not specified take project name from ngs file
if($project eq ""){
	my $ngs_line = `head -n 2 $orig_ngs[0] | tail -n 1`;
	chomp $ngs_line;
	$project = (split(/_/,$ngs_line))[0];
}

if($primer ne ""){
	$cmd = "00_ObiHelp_NGS_control.pl -primer $primer " . join(" ",@orig_ngs) . " > $out_dir/$project.ngs";
}
else{
	$cmd = "00_ObiHelp_NGS_control.pl " . join(" ",@orig_ngs) . " > $out_dir/$project.ngs";
}
exe_cmd($cmd,$verbose,$dry);

#concatenate or link fastq files
my @forward;
my @reverse;
foreach(keys(%fq_filter)){
	my ($for,$rev) = split(/,/,$_);
	push(@forward,$for);
	push(@reverse,$rev);
	if(scalar(keys(%fq_filter)) > 1){ #more than one fastq file pair will be concatenated into one fastq file pair
		$cmd = "cat " . join(" ",@forward) . " > $out_dir/$project\_1.fq.gz & cat " . join(" ",@reverse) . " > $out_dir/$project\_2.fq.gz";
		exe_cmd($cmd,$verbose,$dry);
	}
	else{ #if only one fastq file pair is specified create soft links
		if(-l "$out_dir/$project\_1.fq.gz"){
			$cmd = "rm $out_dir/$project\_1.fq.gz";
			exe_cmd($cmd,$verbose,$dry);
		}
		if(-l "$out_dir/$project\_2.fq.gz"){
			$cmd = "rm $out_dir/$project\_2.fq.gz";
			exe_cmd($cmd,$verbose,$dry);
		}
		$cmd = "ln -s $for $out_dir/$project\_1.fq.gz && ln -s $rev $out_dir/$project\_2.fq.gz";
		exe_cmd($cmd,$verbose,$dry);
	}
}

#create obiwitch config file according to pipeline path
$cmd = "00_ObiMAGIC_obiwitch_config.pl $obiwitch_config $project $out_dir/$project\_1.fq.gz $out_dir/$project\_2.fq.gz $out_dir/$project.ngs $threads $pipeline_path > $out_dir/$project\_obiwitch.ini";
exe_cmd($cmd,$verbose,$dry);

#run obiwitch
$cmd = "01_ObiWitch_main.sh -obiwitch-config $out_dir/$project\_obiwitch.ini";
if($slurm == 1){
	$cmd = "sbatch --job-name=$project\_obiwitch $slurm_opts --wait --wrap=\"$cmd\"";
}
exe_cmd($cmd,$verbose,$dry);

#remove concatenated fastq files
if($keep_tmp == 0){
	if(scalar(keys(%fq_filter)) > 1){
		$cmd = "rm $out_dir/$project\_1.fq.gz $out_dir/$project\_2.fq.gz";
		exe_cmd($cmd,$verbose,$dry);
	}
}

#split demultiplexed files
#not submitted via slurm
$cmd = "02_ObiHelp_split_samples.pl $project\_results/07_$project\_final.fasta $project\_results/00_$project\__ngsfile.tsv";
#if($slurm == 1){
#	$cmd = "sbatch --job-name=$project\_demux $slurm_opts --wrap=\"$cmd\"";
#}
exe_cmd($cmd,$verbose,$dry);

#find demultiplexed final fastas
my @files = ();
if($dry == 0){
	opendir (RESULTS, "$out_dir/$project\_results") or die "ERROR\tCould not open directory $out_dir/$project\_results\n";	#opens the results directory
	while (my $file = readdir(RESULTS)) {					#reads the results directory
		next unless (-f "$out_dir/$project\_results/$file");		#returns only files from the results directory
		if ($file =~ m/08_$project\_final__.*\.fasta/){			#ignores files that don't match the pattern of demultiplexed final fastas
			push (@files, "$out_dir/$project\_results/$file");
		}
	}
	closedir RESULTS;
}
else{
	print STDERR "INFO\tDry run: Would search for files matching 9_$project\_final__.*\\.fasta in $out_dir/$project\_results\n";
}

#create a obiwizard config file for each file resulting from demultiplexing e.g. for each primer
my %configs;
if($dry == 0){
	foreach(@files){
		my $input_assign = $_;
		my $primer = (split(/__/,basename($input_assign)))[-1];
		$primer =~ s/\.fasta$//;
		my $db_path_new = "$out_dir/00_EmysR_ReferenceDB";
		$cmd = "00_ObiMAGIC_obiwizard_config.pl $obiwizard_config $input_assign $primer $threads $db{$primer} $db_path_new $new_taxdump > $out_dir/$project\_obiwizard_$primer.ini";
		exe_cmd($cmd,$verbose,$dry);
		$configs{$primer} = "$out_dir/$project\_obiwizard_$primer.ini";
	}
}
else{
	print STDERR "INFO\tDry run: Would create a obiwizard config file for each 9_$project\_final__.*\\.fasta file\n";
}

#run obiwizard
if($slurm == 1){
	my $array_end = scalar(keys(%configs));
	$cmd = "sbatch --job-name=$project\_obiwizard $slurm_opts --array=1-$array_end --wait --wrap=\"00_ObiMAGIC_obiwizard_jobarray.sh $out_dir\"";
	exe_cmd($cmd,$verbose,$dry);
}
else{
	for(my $i=0; $i < scalar(keys(%configs)); $i++){
		my $primer = (keys(%configs))[$i];
		$cmd = "03_ObiWizard_main.sh -obiwizard-config $configs{$primer}";
		exe_cmd($cmd,$verbose,$dry);
	}
}

#create archives with relevant results

$cmd = "cd $out_dir && tar -czvf $project.raw-results.tar.gz \$(find -name \"*$project*fasta\" -or -name \"*$project*csv\" | sort)";
exe_cmd($cmd,$verbose,$dry);

$cmd = "cd $out_dir && tar -czvf $project.final-results.tar.gz \$(find $project\_results -type f -name \"*RData\" -or -name \"*pdf\" -or -name \"*tsv\" -or -name \"*csv\" -or -name \"*txt\" | sort)";
exe_cmd($cmd,$verbose,$dry);

exit 0;
