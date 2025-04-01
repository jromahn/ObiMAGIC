#!/usr/bin/env bash
work_dir=$(pwd)
echo "$(realpath $0) $*"

########################################################################################################
# ObiTools4 Script for Metabarcoding --> ObiWitch
# Pipeline written by Juliane Romahn ( email: romahnjuliane@gmail.com)
# version: "0.1" - 05.03.2025
version="0.1"
###################################################

#declare ARGV variables saved as default to overwrite config
project_DF=""    
fq_files_DF=""
read1_DF="" 
read2_DF="" 
ngs_file_DF="" 

#declare pipeline used variables
project=""    
read1="" 
read2="" 
ngs_file="" 

#set default config path
my_path=$(dirname $0)
config_file="$my_path/00_ObiScripts/config_ObiWitch.ini"

######################################################################################################################################
###### Prepare handling ARGV input
######################################################################################################################################


# Function to display help message
usage() {
    execution=$(basename $0)
    echo "Usage: $execution -ngs <demultiplexing_file> -fq <paired_1.fq.gz,paired_2.fq.gz> -project <custom project name> -obiwitch-config <config_file> -o <path_to_result>"
    echo ""
    echo "Options:"
    echo "  -fq <file1>   			Paths to two fastq files containing paired Illumina reads comma sperated"
    echo "  -ngs <demultiplexing_file>   	Specify the demultiplexing file containing primer, sample ID and tag information"
    echo "  -obiwitch-config <config_file>	Custom config file for obiwitch if you want to modify the settings"
    echo "  -project <project_name>		Specify projects name"
    echo "  -o <path_to_result>   		Path where all results should be stored, Default: working directory"
    echo "  -version    			Show ObiWitch version and exit"
    echo "  -h, --help    			Show this help message and exit"
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
        -fq)
            fq_files="$2"
            shift 2
            # Set IFS to comma and read into two variables
            IFS=',' read -r read1_DF read2_DF <<< "$fq_files"
            ;;
        -project)
            project_DF="$2"
            shift 2
            ;;
        -ngs)
            ngs_file_DF="$2"
            echo $ngs_file_DF
            shift 2
            ;;
        -obiwitch-config)
            config_file="$2"
            shift 2
            ;;
        -o)
            work_dir="$2"
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
[[ -n "$read1_DF" ]] && read1="$read1_DF"
[[ -n "$read2_DF" ]] && read2="$read2_DF"
[[ -n "$ngs_file_DF" ]] && ngs_file="$ngs_file_DF"
[[ -n "$project_DF" ]] && project="$project_DF"


######################################################################################################################################
###### Set output variables
######################################################################################################################################

######### Input and results folder
start=$(date)
now=$(date)

#work_dir=$(pwd)
echo $work_dir
output=$work_dir/${project}_results
mkdir -p  $work_dir
mkdir -p  $output

##############################
##change tmp  because of storage issues otherwise
new_temp=$( echo "$work_dir/${project}_tmp/")
mkdir -p $new_temp
# Export TMPDIR to point to the new directory
export TMPDIR=$new_temp
# Your script logic goes here
echo "Using new temporary directory: $TMPDIR" 
##############################

#### Output variables
date_log=$(date '+%Y-%m-%d')
readme_file_stats="$output/stats_ASV_counts__${project}_${date_log}_bioinformatics.txt" 
readme_file="$output/read_me__${project}_${date_log}_bioinformatics.txt" 
stats_file="$output/stats_files__${project}_${date_log}_bioinformatics.txt"
new_ngsfile="$output/00_${project}__ngsfile.tsv"

##### Create README
printf "Project name: \t $project \n" >> $readme_file
printf "Date: \t $now \n" >> $readme_file
printf "Working directory: \t $work_dir \n" >> $readme_file
printf "Forward reads : \t $read1 \n" >> $readme_file
printf "Reverse reads : \t $read2 \n" >> $readme_file
printf "Demultiplex file: \t $ngs_file \n" >> $readme_file
printf "Thread number : \t $threads\n" >> $readme_file
printf "OBiWitch Version: \t $version\n" >> $readme_file
printf "OBiTools4 Version: \t" >> $readme_file
obi4_obiannotate --version &>> $readme_file
printf "Using new temporary directory: \t $TMPDIR \n" >> $readme_file
echo "" >> $readme_file

echo "cp $ngs_file $new_ngsfile"
cp $ngs_file $new_ngsfile
echo ""

## default files created in intermediate steps 
### has to be changed depending on the acutal biinformatic steps
merged_file="$output/01_merged.fastq.gz"
step3_file="" #spaceholder
step4_file="" #spaceholder


###########################################
## test if input files exist
error=0
if  [ ! -f $read1 ] || [ ! -f $read2 ] || [ ! -f $ngs_file ];then
        if  [ ! -f $ngs_file ];then
            error=$(($error+1))
        fi
        if  [ ! -f $read1 ]|| [ ! -f $read2 ];then
            error=$(($error+2))
        fi 
fi
if  [ $error -gt 0 ];then  
        if  [ $error -eq 1 ];then    printf "CAUTION! \t $ngs_file  not exist \n"; exit 1
        elif  [ $error -eq 2 ];then  printf "CAUTION! \t $read1 &/or $read2  not exist \n"; exit 2
        elif  [ $error -eq 3 ];then  printf "CAUTION! \t $ngs_file & ($read1 &/or $read2)  not exist \n"; exit 3 
        fi
fi
###########################################


######################################################################################################################################
###### Check  FLAG_DEMULTI_FIRST -- demultiplex if true
######################################################################################################################################
if  [ "$FLAG_DEMULTI_FIRST" == "TRUE" ]; then

    echo "Decision : demultiplex before merging" >> $readme_file

    #declare variables
    output_path="$output/01_demultiplexed_all"
    output_unknown="$output/01_unidentified"
    base_name=$(basename "$new_ngsfile" .tsv)
    output_file1="$output/${base_name}_forward.fasta"
    output_file2="$output/${base_name}_reverse.fasta"

    # Ensure output files are empty
    > "$output_file1"
    > "$output_file2"

    ## create new folders
    mkdir -p $output_path
    mkdir -p $output_unknown

    #declare boolean "hash"
    declare -A seen_forward
    declare -A seen_reverse

    # Read the ngs file line by line, save unique tags and primer combination for cutadapt
    while IFS= read -r line; do
        # Skip lines starting with #
        if  [[ "$line" =~ ^# ]]; then
            continue
        fi

        # Extract the third column (before and after the colon) -- tags , the fourth -- forward primer, and the fifth columns -- reverse primer
        sample=$(echo "$line" | awk '{print $2}')
        third_col_before=$(echo "$line" | awk '{split($3, a, ":"); print a[1]}')
        third_col_after=$(echo "$line" | awk '{split($3, a, ":"); print a[2]}')
        fourth_col=$(echo "$line" | awk '{print $4}')
        fifth_col=$(echo "$line" | awk '{print $5}')
        
        #combine tag and primer
        forward="${third_col_before}${fourth_col}"
        reverse="${third_col_after}${fifth_col}"
        
        # append to respective output files if not already known
        if  [[ -z ${seen_forward[$forward]} ]];then
            #echo $forward
            printf ">$forward\n$forward\n" >> "$output_file1"  
            seen_forward[$forward]=1
        fi
        if  [[ -z ${seen_reverse[$reverse]} ]];then
            printf ">$reverse\n$reverse\n" >> "$output_file2"    
            seen_reverse[$reverse]=1
        fi
    done < "$new_ngsfile"

    # extract cutadapt version
    printf "cutadapt Version: \t" >> $readme_file
    cutadapt --version &>> $readme_file

    #demultiplex with cutadapt ( minimum length is 10!!!)
    command="cutadapt -e $mismatches --no-indels --minimum-length 10 --cores $threads -g file:${output_file1} -G file:${output_file2}  \
        -o $output_path/{name1}-{name2}.1.fastq.gz -p $output_path/{name1}-{name2}.2.fastq.gz \
        $read1 $read2"

    echo $command
    $command

    # remove empty files
    for file in $( ls $output_path/*fastq.gz ); do if [[ $(zcat $file | head | wc -l) -eq 0  ]]; then  rm $file ; fi ; done
    sleep 1

    if [[ $(ls $output_path/ | wc -l) -eq 0  ]]; then
	echo "No files! Sample demultiplexing was not successful. Check the NGS demultiplexing file."
	exit 4
    fi

    #mv unidentified tags
    mv $file $output_path/*unknown*.fastq.gz $output_unknown
fi 

######################################################################################################################################
###### Merge/ Recover full sequence reads from forward and reverse partial reads
######################################################################################################################################
echo "### Recover full sequence with obipairing via merging R1 & R2 files"
if  [ "$FLAG_DEMULTI_FIRST" == "TRUE" ]; then

    merged_file="$output/02_merged.fastq.gz"
    output_known="$output/01_demultiplexed"
    mkdir -p $output_known

    LOG_mergin="$output/01_demultiplexed_LOG.txt"

    #### use ngs file again to merge each  sample independently (easier to reformat the files)
    while IFS= read -r line  ; do

        #extract important information
        sample=$(echo "$line" | awk '{print $2}')
        forward_tag=$(echo "$line" | awk '{split($3, a, ":"); print a[1]}')
        reverse_tag=$(echo "$line" | awk '{split($3, a, ":"); print a[2]}')
        forward_primer=$(echo "$line" | awk '{print $4}')
        reverse_primer=$(echo "$line" | awk '{print $5}')
    
        #test if both files exists
        if [[ -f $output_path/$forward_tag$forward_primer-$reverse_tag$reverse_primer.1.fastq.gz && -f $output_path/$forward_tag$forward_primer-$reverse_tag$reverse_primer.2.fastq.gz ]]; then 
        echo $sample 
        #execute the merging
        echo "$output_path - $sample" 
        printf "$forward_tag$forward_primer\t$reverse_tag$reverse_primer\t$sample\n" >> $LOG_mergin

        obi4_obipairing --max-cpu $threads --compress --min-overlap=10 -F "$output_path/$forward_tag$forward_primer-$reverse_tag$reverse_primer.1.fastq.gz" \
            -R "$output_path/$forward_tag$forward_primer-$reverse_tag$reverse_primer.2.fastq.gz" >  "$output_known/$sample.fastq.gz" 
        
        sleep 1

        #remove existing files
        rm  "$output_path/$forward_tag$forward_primer-$reverse_tag$reverse_primer.1.fastq.gz"
        rm  "$output_path/$forward_tag$forward_primer-$reverse_tag$reverse_primer.2.fastq.gz"

        #add new line character
        echo "" | gzip - | cat - >> $output_known/$sample.fastq.gz
          
        # unzipped change header and add sample name
        zcat $output_known/$sample.fastq.gz  | awk -v name="$sample" 'NR%4==1 {sub(/\}$/, ",\"sample\":\""name"\"}");} {print}'| gzip  > $output_known/$sample.merged.obi4.fastq.gz

        fi 
    done < $new_ngsfile
    #exit

    #### combine the merged and demultiplexed files
    touch $merged_file
    for file in $( ls $output_known/*.merged.obi4.fastq.gz ); do

        #add new line character
        echo "" | gzip - | cat - >> $$file
        #echo $file
        cat $file >> $merged_file
        sleep 1
        rm $file
    done
    #exit


    #### move all files left because they are also unknown samples ( leftovers can exist if primers were multiplexed)
    echo "move all unknown tag combination"

    #forward read
    touch $output/01_unidentified.R1.fastq.gz
    #for file in $( ls $output_unknown/*.fastq.gz & ls $output_path/*.fastq.gz  ); do
    for file in "$output_unknown"/*.1.fastq.gz "$output_path"/*.1.fastq.gz; do
        tag=$( basename $file .1.fastq.gz)
        zcat $file  | awk -v name="$tag" 'NR%4==1 {sub(/\}$/, ",\"sample\":\""name"\"}");} {print}'| gzip  >> $output/01_unidentified.R1.fastq.gz
    done 
    #reverse read
    touch $output/01_unidentified.R2.fastq.gz
    #for file in $( ls $output_unknown/*.fastq.gz & ls $output_path/*.fastq.gz  ); do
    for file in "$output_unknown"/*.2.fastq.gz "$output_path"/*.2.fastq.gz; do
        tag=$( basename $file .2.fastq.gz)
        zcat $file  | awk -v name="$tag" 'NR%4==1 {sub(/\}$/, ",\"sample\":\""name"\"}");} {print}'| gzip  >> $output/01_unidentified.R2.fastq.gz
    done 
    rm -r $output_unknown

    #mv $output_path/*.fastq.gz $output_unknown

    echo "Remove  $output_path"
    #ls $output_path/*.fastq.gz | wc -l
    rm -r $output_path

    ## create readme output of cutadapt
    printf  "#1.Step \t demultiplexing of the samples with cutadapt(!) with length filter of minimum 10bp \n " >> $readme_file
    printf  "01_demultiplexed \t folder contains demultiplexed files of the samples \n" >> $readme_file
    printf  "01_unidentified \t folder contains files with unknown tag combination or missing tag or primers sequences \n\n" >> $readme_file

    ## create readme output of merged samples
    filename=$(basename $merged_file)
    printf  "#2.Step \t merging of the samples again with ObiTools4 for each sample file\n" >> $readme_file
    printf  "$filename \t file of demultiplexed and concatinated/combined sample files \n\n" >> $readme_file
elif  [[ "$FLAG_DEMULTI_FIRST" == "FALSE" && ! -f "$merged_file" ]]; then

    echo "Decision : merging before demultiplexing" >> $readme_file
    #define options
    options="--max-cpu $threads --compress"

    ########################
    if [ "$FLAG_remove" == "Y" ]; then
            echo "obi4_obipairing $options --min-identity=$minimum_alignment_score --min-overlap=10 -F $read1 -R $read2 > $merged_file"
            obi4_obipairing $options --min-identity=$minimum_alignment_score --min-overlap=10 -F $read1 -R $read2 > $merged_file
            ## 
            printf  "#1.Step \t recover full sequence read \n" >> $readme_file
            printf  "1_merged.fastq.gz \t Recovered sequence of library 1 with alignmentscore $minimum_alignment_score  and overlap of 10 \n" >> $readme_file
    else
            echo "obi4_obipairing $options --min-overlap=10 -F $read1 -R $read2 > $merged_file"
            obi4_obipairing $options --min-overlap=10 -F $read1 -R $read2 > $merged_file

            ## readme
            printf  "#1.Step \t recover full sequence read \n" >> $readme_file
            printf  "1_merged.fastq.gz \t Recovered sequence of library 1 without alignmentscore and overlap of 10  \n" >> $readme_file
    fi
fi
echo "" >> $readme_file

######################################################################################################################################
###### remove unaligned sequence records
######################################################################################################################################
if  [ "$FLAG_concat" == "Y"  ]; then

    #change outputfiles depending on bioinformatic steps
    if  [ "$FLAG_DEMULTI_FIRST" == "TRUE" ]; then
        step3_file="$output/03_aligned.fastq.gz"
        step4_file=$step3_file # define input for dereplication
        step="3"
    else
        step3_file="$output/02_aligned.fastq.gz"
        step="2"
    fi
    #execute script if file is not existing
    if  [ ! -f $step3_file ];then
        if  [ -f $merged_file ] && [ $( stat -c%s $merged_file ) -gt 0 ]; then
            echo "Remove unaligned sequences" 
            obi4_obigrep --max-cpu $threads  --compress -p 'annotations.mode != "join"' $merged_file > $step3_file

            printf  "#$step.Step \t remove unaligned/just concatenated sequence \n" >> $readme_file
            printf  "02_aligned.fastq.gz \t Removed unaligned sequences from dataset \n\n" >> $readme_file
        else
            printf "CAUTION! \t $merged_file \t doesn't exist.\n"
            exit 5
        fi
    fi 
else 
    step3_file=$merged_fil
    step4_file=$merged_file
fi
######################################################################################################################################
###### Assign each sequence record to the corresponding sample/marker combination if not done before
######################################################################################################################################
if  [ "$FLAG_DEMULTI_FIRST" == "FALSE" ]; then
    if  [ ! -f  $output/03_demultiplexed.fastq.gz ] || [ $( stat -c%s "$output/03_demultiplexed.fastq.gz" ) -eq 0 ]; then
        # define output file
        step4_file=$output/03_demultiplexed.fastq.gz
        #sleep 5
        if  [ -f "$step3_file" ] && [ $( stat -c%s "$step3_file" ) -gt 0 ]; then
            #obi4_obimultiplex --max-cpu $threads --compress -t $new_ngsfile -e $mismatches -u $output/03_unidentified.fastq.gz $step3_file > $step4_file
	    obi4_obimultiplex --max-cpu $threads --compress --tag-list $new_ngsfile -e $mismatches -u $output/03_unidentified.fastq.gz $step3_file > $step4_file
            
	    
            printf  "#3.Step \t demultiplex sequences with $$mismatches mismatches allowed \n" >> $readme_file
            printf  "03_unidentified.fastq.gz \t Sequences with unidentified tag combination  \n" >> $readme_file
            printf  "03_demultiplexed.fastq.gz \t Demultiplexed Sequences with removed unidentified tag combination  \n" >> $readme_file
            printf  "${project}__ngsfile.tsv \t Cleaned and reformated demultiplex file to change most special characters \(just in case\) used for demultiplexing \n" >> $readme_file      
            
        else
            printf "CAUTION! \t $step3_file  \t doesn't exist.\n"
            exit 6
        fi
    fi
fi

### stop analysis if demultiplexed file is empty
if  [ "$FLAG_DEMULTI_FIRST" == "FALSE"  ]; then 

    if  [  $( stat -c%s "$output/03_demultiplexed.fastq.gz" ) -eq 0  ];then 
        printf "CAUTION! \t $output/03_demultiplexed.fastq.gz  \t is empty Something went wrong while demultiplexing.\n"
        printf "CAUTION! \t Do you have the correct primer sequences in the demultiplex file?!?! \n"

        # create file for error searching
        echo "01_ObiWitch_unidentified_main.pl $output/03_unidentified.fastq.gz"
        01_ObiWitch_unidentified_main.pl $output/03_unidentified.fastq.gz
        printf  "03_unidentified__reasons.tsv \t Listing the reasons why sequences could not be demultiplexed  \n\n" >> $readme_file
        printf  "03_unidentified__diagnosis_plot.pdf \t Diagnostic plot of common errors, tags and other issues in unidentified \n\n" >> $readme_file

        exit 7
    fi 
    ## if unidentified is larger than 15% of demultiplexed file
    # Get the sizes of the files in bytes
    size1=$(stat -c%s "$output/03_demultiplexed.fastq.gz")
    size2=$(stat -c%s "$output/03_unidentified.fastq.gz")

    # Calculate 15% of  multiplexed file size
    fifteen_percent_of_size1=$((size1 * 15 / 100))

    if  [ $size2 -ge $fifteen_percent_of_size1 ]; then 
        printf "$output/03_demultiplexed.fastq.gz  \t is smaller than 03_unidentified.fastq.gz -> Something might went wrong.\n" >> $readme_file

        # create file for error searching
        echo "01_ObiWitch_unidentified_main.pl $output/03_unidentified.fastq.gz"
        01_ObiWitch_unidentified_main.pl $output/03_unidentified.fastq.gz
        printf  "03_unidentified__reasons.tsv \t Listing the reasons why sequences could not be demultiplexed  \n\n" >> $readme_file
        printf  "03_unidentified__diagnosis_plot.pdf \t Diagnostic plot of common errors, tags and other issues in unidentified \n\n" >> $readme_file
        #exit
    fi
    #exit
fi


######################################################################################################################################
###### dereplicate reads into uniq sequences --> ASVs
######################################################################################################################################
if  [ ! -f $output/04_dereplicated.fasta.gz ]; then 
    if  [ -f $step4_file ] && [ $( stat -c%s $step4_file  ) -gt 0 ]; then

        #execute command
        echo "obi4_obiuniq --max-cpu $threads  -O --compress -m sample $step4_file  > $output/04_dereplicated.fasta.gz"
        obi4_obiuniq --max-cpu $threads  -O --compress -m sample $step4_file  > $output/04_dereplicated.fasta.gz


        #### test for possible issues during dereplciation ( known to be problematic)
        reads_demultiplexed=$(obi4_obicount -r $step4_file )
        reads_dereplicated=$(obi4_obicount -r $output/04_dereplicated.fasta.gz)

        # Check if the two strings are identical to check if something went wrong
        if [ "$reads_demultiplexed" == "$reads_dereplicated" ]; then
            echo "---> Dereplication was smooth " # :)"
        else
            echo "Huston, we have a problem!!!"
            echo "Something went wrong during the dereplication"
            echo "Two possible problems: RAM was too small or storage was to small"
            exit 66
        fi

        # create output for the readme file
        printf  "#4.Step \t dereplicate sequences \n" >> $readme_file
        printf  "04_dereplicated.fasta.gz \t Dereplicate reads into uniq sequences  \n" >> $readme_file
    else
        printf "CAUTION! \t $step4_file  \t doesn't exist.\n"
        exit 8
    fi
fi
######################################################################################################################################
###### keep only count and merged_sample information
######################################################################################################################################

if  [ ! -f $output/05_dereplicated.fasta.gz ]; then 
    if  [ -f $output/04_dereplicated.fasta.gz ] && [ $( stat -c%s $output/04_dereplicated.fasta.gz ) -gt 0 ]; then
        echo "obi4_obiannotate -k count -k merged_sample $output/04_dereplicated.fasta.gz > $output/05_dereplicated.fasta.gz"
        obi4_obiannotate -k count -k merged_sample  --compress $output/04_dereplicated.fasta.gz > $output/05_dereplicated.fasta.gz

        printf  "#5.Step \t Clean up information from dereplicated sequences \n" >> $readme_file
        printf  "05_dereplicated.fasta.gz \t Dereplicate reads into uniq sequences and reduced information just to count and merged_sample \n" >> $readme_file
    else
        printf "CAUTION! \t $output/04_dereplicated.fasta.gz \t doesn't exist.\n"
        exit 9
    fi
fi
#exit


######################################################################################################################################
################################################################ Final cleaning ######################################################
######################################################################################################################################

####### Create stats about toal ASV and read number
echo "Create stats about toal ASV and read number - after dereplication"
printf "#Stats of singletons, low and high abundant ASVs\t\n" > $readme_file_stats
printf  "05_dereplicated.fasta.gz \t Total reads: \t">> $readme_file_stats
obi4_obicount -r $output/05_dereplicated.fasta.gz >> $readme_file_stats
printf  "05_dereplicated.fasta.gz \t Total ASVs: \t">> $readme_file_stats
obi4_obicount -v $output/05_dereplicated.fasta.gz >> $readme_file_stats
echo "" >> $readme_file_stats
echo ""

file1="" # denoised
file2=6  # clean
file3=6  # final

###### Denoise the sequence dataset
##Tag the sequences for PCR errors (sequence variants)
if  [ "$FLAG_denoise" == "Y" ]; then
    if  [ ! -f $output/06_denoised.fasta.gz ]; then 
        if  [ -f $output/05_dereplicated.fasta.gz ] && [ $( stat -c%s $output/05_dereplicated.fasta.gz) -gt 0 ]; then
            echo "obi4_obiclean -s sample -r $pcr_threshold -H --compress $output/05_dereplicated.fasta.gz > $output/06_denoised.fasta.gz"
            obi4_obiclean -s sample -r $pcr_threshold -H --compress $output/05_dereplicated.fasta.gz > $output/06_denoised.fasta.gz

            percent=$(echo $pcr_threshold + 100 | bc -l)
            
            printf  "#6.Step \t Denoise the sequence dataset - excluding head seq sequences which are not variants of another sequence with a count greater than 5% of their own count \n" >> $readme_file
            printf  "06_denoised.fasta \t Denoise the sequence dataset from PCR errors (sequence variants) $pcr_threshold \n" >> $readme_file
            printf  " \t Only keep sequences which are not variants of another sequence with a count greater than $percent% of their own count  \n\n" >> $readme_file

        else
            printf "CAUTION! \t $output/05_dereplicated.fasta.gz \t doesn't exist.\n"
            exit 10
        fi
    fi
    file1="06_denoised.fasta.gz"
    file2=$(($file2+1))
    file3=$(($file3+1))

else
    file1="05_dereplicated.fasta.gz"
fi

################################################# Produce statistics after denoicing #################################################
############################################

####### Create stats about toal ASV and read number
echo "Create stats about toal ASV and read number - after denoising"

#general
#printf "#Stats of singletons, low and high abundant ASVs\t$file1\n" >> $readme_file_stats
printf  "$file1 \t Total reads: \t">> $readme_file_stats
obi4_obicount -r $output/$file1 >> $readme_file_stats
printf  "$file1 \t Total ASVs: \t">> $readme_file_stats
obi4_obicount -v $output/$file1 >> $readme_file_stats

## more than 10
echo "" >> $readme_file_stats
printf  "$file1 \t Count >=10 Reads: \t">> $readme_file_stats
obi4_obigrep -p 'sequence.Count() >= 10' $output/$file1 | obi4_obicount -r >> $readme_file_stats
printf  "$file1 \t Count >=10 ASVs: \t">> $readme_file_stats
obi4_obigrep -p 'sequence.Count() >= 10' $output/$file1 | obi4_obicount -v >> $readme_file_stats

## less than 10
printf  "$file1 \t Count <10 Reads: \t">> $readme_file_stats
obi4_obigrep -p 'sequence.Count() < 10' $output/$file1 | obi4_obicount -r >> $readme_file_stats
printf  "$file1 \t Count <10 ASVs: \t">> $readme_file_stats
obi4_obigrep -p 'sequence.Count() < 10' $output/$file1 | obi4_obicount -v >> $readme_file_stats

## Singletons
echo "" >> $readme_file_stats
printf  "$file1 \t Singleton Reads: \t">> $readme_file_stats
obi4_obigrep -p 'sequence.Count() == 1' $output/$file1 | obi4_obicount -r >> $readme_file_stats
printf  "$file1 \t Singleton ASVs: \t">> $readme_file_stats
obi4_obigrep -p 'sequence.Count() == 1' $output/$file1 | obi4_obicount -v >> $readme_file_stats


spaceholder=$( basename $readme_file_stats )
echo "" >> $readme_file
printf  "$spaceholder \t Stats for $file1 about amount of Singletons, ASV count higher and lower than 10 reads   \n" >> $readme_file


##################################################### Count and length filtering #####################################################
####### Keep only the sequences having a count greater or equal to X and a length shorter than Y bp
echo "-------> Prefilter sequences if asked for"

if [ "$FLAG_length" == "Y" ] || [ "$FLAG_count" == "Y"  ]; then
    file3=$(($file3+1))
    if [ -f $output/$file1 ] && [ $( stat -c%s  $output/$file1 ) -gt 0 ]; then
        if [ "$FLAG_length" == "Y" ] && [ "$FLAG_count" == "Y"  ] ; then

            formatted_number=$(printf "%02d" "$file2")
            file2="${formatted_number}_sequences_${min_length}.${max_length}_c${count}_cleaned.fasta.gz"
            
            if [  ! -f $output/$file2 ]; then
                obi4_obigrep --compress -l $min_length -L $max_length -p  "sequence.Count() >= $count" $output/$file1 > $output/$file2
                printf  "#7.Step \t Clean up sequences \n" >> $readme_file
                printf  "$file2 \t Cleaned sequences with length beween $min_length & $max_length  and count of $count  \n\n" >> $readme_file
            fi

        elif [ "$FLAG_length" == "Y" ] && [ "$FLAG_count" == "N" ]; then

            formatted_number=$(printf "%02d" "$file2")
            file2="${formatted_number}_sequences_${min_length}.${max_length}_cleaned.fasta.gz"

            if [  ! -f $output/$file2 ]; then
                obi4_obigrep --compress -l $min_length -L $max_length  $output/$file1 > $output/$file2
                printf  "#7.Step \t Clean up sequences \n" >> $readme_file
                printf  "$file2 \t Cleaned sequences with length beween $min_length & $max_length  and not after count \n\n" >> $readme_file
            fi

        elif [ "$FLAG_length" == "N" ] && [ "$FLAG_count" == "Y" ] ; then

            formatted_number=$(printf "%02d" "$file2")
            file2="${formatted_number}_sequences_c${count}_cleaned.fasta.gz"

            if [  ! -f $output/$file2 ]; then
                obi4_obigrep --compress -p  "sequence.Count() >= $count" $output/$file1 > $output/$file2
                printf  "#7.Step \t Clean up sequences \n" >> $readme_file
                printf  "$file2 \t Cleaned sequences count of $count and not after length  \n\n" >> $readme_file
            fi
        else
            printf "CAUTION! \t Something is really weird :D .\n"
            exit 11
        fi
    else
        printf "CAUTION! \t $output/$file1 \t doesn't exist.\n"
        exit 12
    fi
else
    file2=$file1
fi



echo "Generate the final fasta file"
############################################### Generate the final fasta file #########################################################
#if [ ! -f $output/$file3 ]; then
    if [ -f $output/$file2 ] && [ $( stat -c%s $output/$file2 ) -gt 0 ]; then

        formatted_number=$(printf "%02d" "$file3")
        file3="${formatted_number}_${project}_final.fasta"

        #execute command
        echo "obi4_obiannotate -k count -k merged_sample $output/$file2 > $output/$file3"
        obi4_obiannotate -k count -k merged_sample $output/$file2 > $output/$file3

        printf  "#Last Step \t  Generate the final fasta file \n" >> $readme_file
        printf  "$file3 \t Final fasta file demultiplexed, deprelicated and cleaned if wanted from ObiTols4 with ObiTools ASV names \n" >> $readme_file

    else
        printf "CAUTION! \t $output/$file2 \t doesn't exist. \n"
        exit 13
    fi
#fi

printf  "\n#Extra Step \t  Additional Files & statistics \n" >> $readme_file


######################################################################################################################################
############################################################ Diagnosis plot ##########################################################
######################################################################################################################################

if [ -f $output/$file3 ] && [ $( stat -c%s $output/$file3 ) -gt 0 ]; then
    file3=$( basename $file3 )
    new_ngsfile=$( basename $new_ngsfile )

    echo "TMPDIR=\"$new_temp\" 01_ObiWitch_diagnostic.R $output $file3 $new_ngsfile $threads $(which 00_OBIMAGIC_functions.R)"
    TMPDIR="$new_temp" 01_ObiWitch_diagnostic.R $output $file3 $new_ngsfile $threads $(which 00_OBIMAGIC_functions.R)

    file3_name=$( basename $file3 .fasta )


    echo  "" >> $readme_file
    printf  "${file3_name} \t ASVs are renamed in this step and shorten   \n" >> $readme_file
    printf  "${file3_name}__renamed_sequences.fasta \t Fasta file with renamed ASVs for further assignment in mothur & BOLDdigger \n" >> $readme_file
    printf  "${file3_name}__mothur_counttable.csv \t Count table with renamed ASVs for further assignment in mothur  \n" >> $readme_file
    printf  "${file3_name}__general_infos.csv, \t General information of sequences/ASVs provided by ObiTools4 including original and new ASV names\n" >> $readme_file
    printf  "${file3_name}__community.RData \t RData community table, rows - sample, columns - ASVs (new ASV names!)  \n" >> $readme_file
    printf  "${file3_name}__diagnosis_plots_bioinformatics.pdf \t First diasgnosis plots to check samples and controls \n" >> $readme_file

    if [ ! -f "$output/${file3_name}__diagnosis_plots_bioinformatics.pdf" ] || [ $( stat -c%s "$output/${file3_name}__diagnosis_plots_bioinformatics.pdf" ) -eq 0 ]; then
        echo "$output/${file3_name}__diagnosis_plots_bioinformatics.pdf can not be found"
        exit 14
    fi

fi
now=$(date)
printf "\nDate stoped: \t $now \n" >> $readme_file

#exit

######################################################################################################################################
############################################################ Produce FILE statistics #################################################
######################################################################################################################################


echo "Calulcate File statistics"
printf "#ObiTools4 Stats \t \n" > $stats_file
printf "#Filename\tAverage.Amplicon.Length \tTOTAL.ASVs TOTAL.Reads\n" >> $stats_file

#create list of all files of interest
outputfiles=$(find "$output" -type f -maxdepth 1 \( -name "*.fastq.gz" -o -name "*.fasta.gz" -o -name "*.fasta" \) ! -name "*ngsfile*" ! -name "*renamed*" | sort)
# Add the two additional files to the list
outputfiles="$read1"$'\n'"$read2"$'\n'"$outputfiles"

#for file in $( ls $output/*.fastq.gz && ls $output/*.fasta.gz && ls $output/*.fasta ); do 
for file in $outputfiles; do
    filename=$(basename $file)

    echo "$filename"
    # Extract the file extension
    extension="${filename##*.}"
    
    #Calulcate sequence length
    seq_length=0
    if [[ "$extension" == "fasta" || "$extension" == "fa" ]]; then
        seq_length=$( awk '/^>/ {if (seqlen){print seqlen; seqlen=0}} !/^>/ {seqlen += length($0)} END {if (seqlen) print seqlen}'\
                        $file | awk '{total += $1; count++} END {print total/count}' )
    else
        # get extension without gz
        filename2=$(basename $file .gz)
        extension="${filename2##*.}"

        if [[ "$extension" == "fasta" || "$extension" == "fa" ]]; then
            seq_length=$( zcat $file | awk '/^>/ {if (seqlen){print seqlen; seqlen=0}} !/^>/ {seqlen += length($0)} END {if (seqlen) print seqlen}'\
                         | awk '{total += $1; count++} END {print total/count}' )
        else
            seq_length=$( zcat $file | awk 'NR%4==2 { total += length($0); count++ } END { print total/count }')
        fi
        
    fi
    #Create Output
    printf "$filename \t $seq_length \t" >> $stats_file
    obi4_obicount -v -r $file >> $stats_file
done
printf  "$stats_file \t Stats for every ObiTools Step to ASV & Reads   \n" >> $readme_file


######################################################################################################################################
## remove tempory directory
rm -r $new_temp
######################################################################################################################################

end=$(date)

echo "Finished: the magic is done"
echo "The analysis started : $start"
echo "The analysis started : $start">> $readme_file
echo "The analysis finished: $end"
echo "The analysis finished: $end" >> $readme_file
