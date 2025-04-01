#!/usr/bin/env bash
work_dir=$(pwd)
echo "$(realpath $0) $*"
########################################################################################################
# ObiTools4 Script for ASSIGNMEN    
# 
###################################################
# Pipeline written by Juliane Romahn ( email: romahnjuliane@gmail.com)
# version: "0.1" - 05.03.2025
version="0.1"
###################################################


######################################################################################################################################
###### Prepare handling ARGV input
######################################################################################################################################

#declare ARGV variables saved as default to overwrite config
input_DF=""    
identifiers_DF=""
ref_DF="" 
taxdump_DF="" 

#declare pipeline used variables
input_assign=""    
identifier="" 
db_path_to_file="" 
db_path_new=""
tax_path_to_file="" 

#set default config path
my_path=$(dirname $0)
config_file="$my_path/00_ObiScripts/config_ObiWizard.ini"

# Function to display help message
usage() {
    execution=$(basename $0)
    echo "Usage: $execution -input <obitools_final.fasta> -identifier <primer> -obiwizard-config <config_file> -ref <reference_database> -new-taxdump <NCBI_taxdump>"
    echo ""
    echo "Options:"
    echo "  -input <obitools_final.fasta>   Paths to obitools formated fasta file"
    echo "  -identifier <primer>   Specifiy the primer name"
    echo "  -obiwizard-config <config_file>	Custom config file for obiwitch if you want to modify the settings"
    echo "  -ref <reference_database>   Path to reference database for assigning"
    echo "  -new-taxdump <NCBI_taxdump>   Custom path to NCBI's new_taxdump directory"
    echo "  -version    Show ObiWitch version and exit"
    echo "  -h, --help    Show this help message and exit"
    exit 0
}

if [[ $# -eq 0 ]];then
    echo "Unknown option: $1"
    usage
fi

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help)
            usage
            ;;
        -version)
            echo $version
            exit 0
            ;;
        -input)
            input_DF="$2"
            shift 2
            ;;
        -identifier)
            identifiers_DF="$2"
            shift 2
            ;;
        -obiwizard-config)
            config_file="$2"
            shift 2
            ;;
        -ref)
            ref_DF="$2"
            shift 2
            ;;
        -new-taxdump)
            taxdump_DF="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# check if config starts with './' or is a relative path (not starting with '/') and add full path
if [[ "$config_file" == ./* ]]; then
    # Remove the leading './' from the file path
    config_file="${config_file#./}"
fi
if [[ "$config_file" != /* ]]; then
    # If it's a relative path (doesn't start with '/'), prepend the working directory
    config_file="$work_dir/$config_file"
fi

#check if config file exists
if [ ! -f $config_file ]; then
    echo "Sorry but your config file doesn't exist ......."
    echo $config_file
    echo "Check again"
    exit 42
fi

#execute config
source $config_file

# overwrite config variables if they are set in ARGV
[[ -n "$input_DF" ]] && input_assign="$input_DF"
[[ -n "$identifiers_DF" ]] && identifier="$identifiers_DF"
[[ -n "$ref_DF" ]] && db_path_to_file="$ref_DF"&& db_path_new="$work_dir/00_indexed_ReferenceDB"
[[ -n "$taxdump_DF" ]] && tax_path_to_file="$taxdump_DF"

####################################################################################################################################################
##################################################  PLEASE DO NOT CHANGE ANYTHING FROM HERE ON #####################################################
####################################################################################################################################################

######################################################################################################################################
###### Make files to absolute paths
######################################################################################################################################


######### Check if the file path starts with './' or is a relative path (not starting with '/')
if [[ "$input_assign" == ./* ]]; then
    # Remove the leading './' from the file path
    input_assign="${input_assign#./}"
fi
if [[ "$db_path_to_file" == ./* ]]; then
    # Remove the leading './' from the file path
    db_path_to_file="${db_path_to_file#./}"
fi
if [[ "$db_path_new" == ./* ]]; then
    # Remove the leading './' from the file path
    db_path_new="${db_path_new#./}"
fi
if [[ "$tax_path_to_file" == ./* ]]; then
    # Remove the leading './' from the file path
    tax_path_to_file="${tax_path_to_file#./}"
fi

######### add full path to file
# input file
if [[ "$input_assign" != /* ]]; then
    # If it's a relative path (doesn't start with '/'), prepend the working directory
    input_assign="$work_dir/$input_assign"
fi
# reference db file
if [[ "$db_path_to_file" != /* ]]; then
    # If it's a relative path (doesn't start with '/'), prepend the working directory
    db_path_to_file="$work_dir/$db_path_to_file"
fi
# path of indexed reference db folder
if [[ "$db_path_new" != /* ]]; then
    # If it's a relative path (doesn't start with '/'), prepend the working directory
    db_path_new="$work_dir/$db_path_new"
fi
# taxonomy file
if [[ "$tax_path_to_file" != /* ]]; then
    # If it's a relative path (doesn't start with '/'), prepend the working directory
    tax_path_to_file="$work_dir/$tax_path_to_file"
fi

######################################################################################################################################
###### Set output variables
######################################################################################################################################


## from input file
output=$( dirname $input_assign )
project="${output##*/}"        # Get the last folder name
project=${project/%_results} # remove result to get project name
echo $project

now=$(date)
#output=$output/${project}_results


now=$(date)
date_log=$(date '+%Y-%m-%d')


file3=""
readme_file="$output/read_me__${project}_${identifier}_${date_log}_assignment.txt" 
stats_file="$output/stats_files__${project}_${identifier}_${date_log}_assignment.txt"
input_db=$db_path_to_file


######################################################################################################################################
################################################ Preprare Reference Database #########################################################
######################################################################################################################################
mkdir -p $work_dir
mkdir -p $db_path_new

# check if database was already indexed otherwise index it
if [ $( grep "obitag_ref_index" $db_path_to_file | wc -l ) -gt 0 ] ; then 
            input_db=$db_path_to_file
else
            #ensure that sequences each have a unique identification & define new db
            input_db=$( basename $db_path_to_file .fasta )
            input_db="$db_path_new/$input_db.indexed.fasta"
            obi4_obirefidx --max-cpu $threads -t $tax_path_to_file $db_path_to_file  > $input_db
fi

#exit
######################################################################################################################################
############################################################ Assignment ##############################################################
######################################################################################################################################

printf "Project name: \t $project \n" >> $readme_file
printf "Date: \t $now \n" >> $readme_file
printf "Working directory: \t $work_dir \n" >> $readme_file
printf "Thread number : \t $threads\n" >> $readme_file
printf "OBiWizard Version: \t $version\n" >> $readme_file
printf "OBiTools4 Version: \t" >> $readme_file
obi4_obiannotate --version &>> $readme_file
echo ""
echo "Start assignment with reference DB: $db_path_to_file "

# created header cleaned files ( 1.step)
input_filename=$( basename $input_assign )
input_filename2=$( basename $input_assign .fasta)
number=$(echo "$input_filename" | grep -o '^[0-9]*')
number=$((10#$number + 1)) # Remove leading zeros to ensure the number is treated as decimal
echo $input_filename
echo $number

output_first=$output/${input_filename2}_clean.fasta

#create length filtered files (2.step)
formatted_number=$(printf "%02d" "$number")
input_filename3="${formatted_number}_${project}_final_${identifier}"
output_second=$output/${input_filename3}_clean_5bp.fasta
short_amplicons=$output/${input_filename3}_clean_tooshort.fasta

# create files for final output ( 3. step)
number=$(( $number + 1 ))
formatted_number=$(printf "%02d" "$number")
file3="${formatted_number}_${project}_final_${identifier}__assigned.fasta"

#################################### 1. Step - remove unnessesary information ####################################
# stats about tax dump file to know when downloaded (which version)
printf  "ObiTools4 \t $obitools_version \n\n" >> $readme_file
printf  "$(basename $input_assign) \t Inputfile for assignment \n" >> $readme_file

## only keep count and sample information
obi4_obiannotate -k count -k merged_sample  $input_assign > $output_first


# stats about tax dump file to know when downloaded (which version)
printf  "$(basename $output_first) \t Cleaned Inputfile for assignment containing only count & sample information \n" >> $readme_file
printf  "$(basename $output_second) \t Cleaned Inputfile for assignment containing sequences longer than 5bp \n" >> $readme_file
printf  "$(basename $short_amplicons) \t Files containing removed sequences  which are shorter than 5 bp \n" >> $readme_file

#################################### 2. Step - remove too short sequences  ####################################
#### remove sequences which are shorter than 5 bp
obi4_obigrep -l 5 $output_first > $output_second
obi4_obigrep -L 4 $output_first > $short_amplicons


#################################### 3. Step - start assignment  ####################################
# start assignment
input_assign=$output_second
LOG_assignment=LOG_assignment_errors_${identifier}_$now.txt
#if small no problems occur
if [ $( grep ">" $input_assign | wc -l ) -lt 1000 ] ; then
        echo "obi4_obitag --max-cpu $threads --no-order -t $tax_path_to_file -R $input_db $input_assign > $output/08_${project}_final_${identifier}__assigned.fasta"
        obi4_obitag --max-cpu 50 --no-order -t $tax_path_to_file -R $input_db $input_assign > $output/08_${project}_final_${identifier}__assigned.fasta 2>>"$LOG_assignment"

else #split into several fasta files due to software struggles
        mkdir -p $work_dir/$identifier
        cd $work_dir/$identifier
        touch LOG_file_obi4_assignment.txt
        headers=$( grep '>' $input_assign | wc -l  )
        line_no=$((( $headers + 999) / 1000 ))
        echo "obi4_obidistribute --pattern toto__%s.fasta --batches $line_no $input_assign" 
        obi4_obidistribute --pattern toto__%s.fasta --batches $line_no $input_assign

        counter=1
        for file in $( ls toto__*.fasta ); do
            echo "obi4_obitag --max-cpu 50  --no-order -t $tax_path_to_file -R $input_db $file >> $output/08_${project}_final_${identifier}__assigned_$counter.fasta">> LOG_file_obi4_assignment.txt
            obi4_obitag --max-cpu 50  --no-order -t $tax_path_to_file -R $input_db $file >> $output/08_${project}_final_${identifier}__assigned_$counter.fasta 2>>"$LOG_assignment"
            counter=$( expr $counter + 1 )
        done
        cd $work_dir
        rm -r $work_dir/$identifier

        
        cat $output/08_${project}_final_${identifier}__assigned_*.fasta > $output/$file3
        rm $output/08_${project}_final_${identifier}__assigned_*.fasta
fi
#exit

###### Check if there was sequence loss due to ObitTools4
if [ ! $( grep ">" $output/$file3 | wc -l ) -eq $( grep ">" $input_assign | wc -l ) ] ; then
    echo  "Lovell, we have a problem!!! - there is a problem : we losing sequences"
    filename_first=$(basename $input_assign )
    seqno_first=$( grep ">" $input_assign | wc -l )
    seqno_second=$( grep ">" $output/$file3 | wc -l )
    printf "$filename_first \t $seqno_first sequences \n $file3 \t $seqno_second sequences \n"
    exit 666
fi


############ Create Read me input
echo >> $readme_file
printf  "$tax_path_to_file \t ">> $readme_file
stat -c '%w' $tax_path_to_file/names.dmp  >> $readme_file

#
printf  "${input_db} \t Reference db for ${file3}  \n" >> $readme_file
printf  "${file3} \t Assigned sequences with ObiTools ASV names \n" >> $readme_file

######################################################################################################################################
############################################################ Diagnosis plot ##########################################################
######################################################################################################################################

if [ -f $output/$file3 ] && [ $( stat -c%s $output/$file3 ) -gt 0 ]; then
    echo "03_ObiWizard_diagnostic.R $output $file3 $tax_path_to_file $threads $(which 00_OBIMAGIC_functions.R)"
    03_ObiWizard_diagnostic.R $output $file3 $tax_path_to_file $threads $(which 00_OBIMAGIC_functions.R)

    #exit
    file3_name=$( basename $file3 .fasta )

    # remove created file after it was needed
    rm $output/${file3_name}_accessionTaxa.sql

    printf  "\n" >> $readme_file
    printf  "${file3_name} \t ASVs are renamed in this step and shorten   \n" >> $readme_file
    printf  "${file3_name}__renamed.fasta \t Fasta table with renamed ASV names    \n" >> $readme_file
    printf  "${file3_name}__community.RData \t RData community table, rows - sample, columns - ASVs  \n" >> $readme_file
    printf  "${file3_name}__diagnosis_plots_assignment.pdf \t Diagnositic plots to check controls and assignments \n" >> $readme_file
    printf  "${file3_name}__taxonomy_info.tsv \t Taxonomic feedback from ObiTools4 for each ASV \n" >> $readme_file
    printf  "${file3_name}__taxonomy_info_COMPLETE.tsv \t Taxonomic feedback from ObiTools4 including NCBI taxdump information with new & old ID names\n" >> $readme_file
    printf  "${file3_name}__taxonomy_SPECIES_summary.tsv \t Taxonomic feedback to each NCBI taxID containing NCBI taxdump information\n" >> $readme_file   
fi
now=$(date)
printf "\nDate stoped: \t $now \n" >> $readme_file


######################################################################################################################################
######################################################### Produce File statistics ####################################################
######################################################################################################################################


echo "Calulcate File statistics"
printf "#ObiTools4 Stats \t \n" > $stats_file
printf "#Filename\tAverage.Amplicon.Length \tTOTAL.ASVs TOTAL.Reads\n" >> $stats_file
outputfiles=$(find "$output" -type f -maxdepth 1 -name "*.fasta" ! -name "*renamed*" | sort)
for file in $( ls $output/*.fasta ); do 
    filename=$(basename $file)

    echo "$filename"
    # Extract the file extension
    extension="${filename##*.}"
    
    #Calulcate seq length
    seq_length=$( awk '/^>/ {if (seqlen){print seqlen; seqlen=0}} !/^>/ {seqlen += length($0)} END {if (seqlen) print seqlen}'\
                        $file | awk '{total += $1; count++} END {print total/count}' )
    #Create Output
    printf "$filename \t $seq_length \t" >> $stats_file
    obi4_obicount -v -r $file >> $stats_file
done
printf  "$stats_file \t Stats for every ObiTools Step to ASV & Reads   \n" >> $readme_file

