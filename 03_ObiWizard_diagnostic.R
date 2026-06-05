#!/usr/bin/env Rscript

########################################
# Diagnosis script for ObiTools4
# Juliane Romahn
# January 14, 2025
# html_document
# execute via: 03_ObiWizard_diagnostic.R Resuults_folder FINAL_fasta_file PATH/to/taxdump threads full/path/to/00_ObiScripts/00_OBIMAGIC_functions.R
########################################

#rm(list=ls())  

############################
## Handle Input
args = commandArgs(trailingOnly=TRUE)
# Help message
usage <- "
Usage: 03_ObiWizard_diagnostic.R <path> <input_file> <taxdump_path> <threads> <functions_path>

Arguments:
  path \t Path to the ObiTools matrix file
  input_file \t Filename of the ObiTools matrix file
  taxdump_path \t Path to Taxdump File (used for ObiMAGIC)
  threads \t Number of threads, necessary for ordinations
  functions_path  \t Path to 00_OBIMAGIC_functions.R file

Options:
  -h, --help \t Show this help message and exit

Example:
  03_ObiWizard_diagnostic.R ObiMagic_results communtiy_matrix.csv /PATH/TO/TAXDUMPFILE 10 00_OBIMAGIC_functions.R
"

# Check for help flag
if (any(args %in% c("-h", "--help"))) {
  cat(usage)
  quit(status = 0)
}

# Check number of arguments
if (length(args) != 5) {
  cat("ERROR: Wrong number of arguments. Expected 5, got", length(args), "\n")
  cat(usage)
  quit(status = 1)
}

## read input if check was successfull
path=args[1] #"Malawi_2023_eDNA_results/"
input_file =args[2] #"7_Malawi_2023_eDNA_final.fasta" #
taxdump_path=args[3] 
threads= args[4]
functions_path=args[5]

## if you execute it in RStudio or similar
#path="BiodivSoup_trail_results"
#input_file ="10_BiodivSoup_final_FolDegenRev__assigned__matrix.csv"
#setwd("/Users/juliane/Documents/00_Work_SGN/00_PhytoArk/XX_PAPERS/2025sub_ObiMAGIC/01_Benchmarking_resubmission/testing")
#path="01_BiodivSoup_apscale/11_read_table/data/"
#input_file ="02_data_final_FolDegenRev__assigned__matrix.csv"
#taxdump_path="."
#threads=1
#functions_path="00_OBIMAGIC_functions.R"
############################

###### Library
#install.packages("jsonify")
#install.packages("devtools")
 # install.packages("ghibli")
#devtools::install_git("https://git.metabarcoding.org/obitools/obitools4/robireadfasta.git")
#install.packages('taxize', repos = c('https://ropensci.r-universe.dev', 'https://cloud.r-project.org'))
#library(ROBIFastread)
library(tidyverse)
library(vegan)
library(magrittr)
library(seqinr)
library(ggpmisc)
library(gapminder) # allow commas in numbers --> +scale_y_continuous(labels = scales::comma)
library(taxonomizr) ## NCBI taxonomy
library(taxize) ## NCBI taxonomy
library(ggpubr) # ggplot stats
library(ghibli) ## ggplot colors
library(cowplot) ## add text under plot
library(plotly) ##ggtreemap
library(treemapify) ##ggtreemap

options(scipen = 999) # stop scientific notation


#### read in all functions
source(functions_path)


############################ prepare output
input_file_pattern = gsub( "assigned__matrix.csv","assigned",input_file)
input_file2= paste(input_file_pattern,"__taxonomy_info.csv", sep="")

output=input_file_pattern
#output=gsub(".fasta", "", input_file)


### create folder for the figures
figure_folder=paste(input_file_pattern, "diagnostic_figures", sep="_")
output_figures=file.path(path,figure_folder)
print("Create output folder")
if (!dir.create(output_figures)){
  dir.create(output_figures)
}

############################ read in input

#read it file with NCBI IDs l
print("Load NCBI IDs")
print(file.path(path,input_file2))
data_ids <- read.table(file.path(path,input_file2), header=T, sep=",", fill = TRUE, quote = "\"")

# read in community matrix
print("Load community matrix")
community <- read.table(file.path(path,input_file), header=T, sep=",", check.names = FALSE) # check.names = FALSE otherwise : will be replaced, okay since it is changes soon
community <- unique(community)
row.names(community) <- community$id
community$id <- NULL




###### check if there is a renaming file, otherwise create new
print("Change ASV names")
### search for renaming file
files <- list.files(path = path, 
                    pattern = "^[0-9]+_.*__renamePattern\\.tsv$", 
                    full.names = TRUE)

# rename with existing pattern or create new
if (length(files) > 0) {
  renaming <- read.table(file=files[1], header = T, sep ="\t")
  for (i in 1:nrow(renaming)) {
    names(community)[names(community) == renaming$id[i]] <- renaming$id_new[i]
  }
}else{
  old_asv_names <- colnames(community)
  new_asv_names <- paste(asv_start,sprintf("%06d", c(1:length(colnames(community)))), sep="_")
  colnames(community) <- new_asv_names
  
  #overview how things are renames
  renaming <- data.frame(id= old_asv_names, id_new=new_asv_names )
}
data_ids <- data_ids %>% left_join(renaming, by="id")

########### prepare taxonomy

### get taxonomy & save it in same folder then the taxdump
print("Load NCBI")
loc_ncbi <- c(paste(taxdump_path,  "names.dmp", sep="/"),
              paste(taxdump_path,  "nodes.dmp", sep="/")) 
sqlFile <- file.path(path, paste(output,'accessionTaxa.sql', sep="_"))
read.accession2taxid(list.files(loc_ncbi),sqlFile,indexTaxa=TRUE,overwrite=TRUE)
taxaNames<-read.names.sql(loc_ncbi[1],sqlFile)
taxaNodes<-read.nodes.sql(loc_ncbi[2],sqlFile)


########  get taxonomy
taxonomy <- data_ids %>%
      mutate(obitag_bestid100 = obitag_bestid *100)%>%
      mutate(SeqLength=nchar(sequence))%>%
      relocate(SeqLength, .after=sequence)

#depending on the obitools version extract scientific name and taxid
if (grepl("^taxon", taxonomy$taxid[1])) {
  #scientific_name looks like this: taxon:2759 [Eukaryota]@superkingdom
  taxonomy <- taxonomy %>% 
    mutate(save=taxonomy$taxid ) %>%
    mutate(taxid= gsub("taxon:(\\d+).*", "\\1", save, perl=T))%>% 
    mutate(scientific_name= gsub("taxon:\\d+.*\\[(.*)\\]?.*$", "\\1", save, perl=T))%>%
    mutate(save=NULL)
}

# feedback by obitools 4  
write.table(taxonomy, file=file.path(path, paste(output, "taxonomy_info.tsv", sep="__")), quote=FALSE, sep = "\t", row.names = F)
 
# taxonomy of ncbi taxon ids  
detailed_taxonomy <- ncbi_taxonomy(taxonomy$taxid)
write.table(detailed_taxonomy, file=file.path(path, paste(output, "taxonomy_SPECIES_summary.tsv", sep="__")), quote=FALSE, sep = "\t", row.names = F)
  
# combined overview of sequences and their taxonomy
detailed_taxonomy <- merge( taxonomy, detailed_taxonomy,by.x = "taxid", by.y="taxID")
write.table(detailed_taxonomy, file=file.path(path, paste(output, "taxonomy_info_COMPLETE.tsv", sep="__")), quote=FALSE, sep = "\t", row.names = F)


##################

##############
save(community, file=file.path(path, paste(output, "community.RData", sep="__")))



################################## PLOTTING ################################## 
## prepare data and sample identification

#extract replicate and sample information from communtiy rownames
community <- extra_sample_names(community, type_identifier)

# get index from asv and samples
index <- grep(asv_start, colnames(community))
index_samples <-which(!grepl(grep_controls,community$type))
pr <- unique(community$primer)

community$total.reads <- rowSums(community[,index])
community$total.asvs <- rowSums(1*(community[,index] >0))


##############

pdf(file = file.path(path, paste(output, "diagnosis_plots_assignment.pdf", sep="__")), paper = "a4r")

# secure that sampel types are always named the same by sorting them before
type_List <-community %>% arrange(type)%>% pull(type) %>% unique()
## define colors
selected <- colorRampPalette(colors)(length(unique(community$type)))
color_df<- data.frame(type= type_List, color=selected)

###### project overview
community_sum <- categorize_community(community, asv_start)



## ggtree map: sample type overview
plot <- ggplot(community_sum, aes(area = replicates, fill = type, label= type, subgroup=type_long)) +
  geom_treemap()+
  geom_treemap_subgroup_border(color="white", size=4)+
  geom_treemap_subgroup_text(place ="center",color="white",reflow = T)+
  geom_treemap_text(color="black", alpha = 0.7)+
  scale_fill_manual(values=color_df$color, breaks = color_df$type)+
  labs(title= paste("Project overview:",output , "\nPrimer:", unique(community$primer)),
       subtitle = "Size represents replicate number\nColor coding is consistent within the document")+
  plot_theme  
print(plot)
plot_function(plot, file.path(output_figures, "01_project_overview__treemap"))


#print it as table
community_sum_reads <- community_sum %>% select(-total_asv, -average_asv, -median_asv, -replicates, -samples )
community_sum_asvs <- community_sum %>% select(-total_reads, -average_reads, -median_reads, -replicates, -samples )
community_sum <- community_sum %>% select(-average_reads, -median_reads,-average_asv, -median_asv )

#print it as table
plot<- ggplot() + theme_void() + annotate(geom="table",x=1,y=1,label=list(community_sum), size=3)+
  labs(title= paste("Project overview:",output ),
       subtitle = "Total includes multiplexed primers (if existing)")
print(plot)
plot_function(plot, file.path(output_figures, "02_project_overview_table"))


plot<- ggplot() + theme_void() + annotate(geom="table",x=1,y=1,label=list(community_sum_reads), size=3)+
  labs(title= paste("Project overview:",output ),
       subtitle = "Total includes multiplexed primers (if existing)")
print(plot)
plot_function(plot, file.path(output_figures, "02_project_overview_table__reads"))

plot<- ggplot() + theme_void() + annotate(geom="table",x=1,y=1,label=list(community_sum_asvs), size=3)+
  labs(title= paste("Project overview:",output ),
       subtitle = "Total includes multiplexed primers (if existing)")
print(plot)
plot_function(plot, file.path(output_figures, "02_project_overview_table__Asv"))


# colors assignment
pal <- as.vector(ghibli_palette(name = "MarnieMedium2", n = 10, type = "continuous"))
color_df_assign <- data.frame(type=seq(10, 100, by = 10), color=pal)


# Distribution of Similarity to closest Database Hit
plot <- ggplot(detailed_taxonomy, aes(x=obitag_bestid100)) + 
  geom_histogram(fill="deeppink4", alpha=0.8)+
  geom_vline(xintercept = mean(detailed_taxonomy$obitag_bestid100), linetype="dotted", color = "black", linewidth=1.5)+
  plot_theme+
  labs(title= "Distribution of Similarity to closest Database Hit", 
                     subtitle= paste("Dotted line represents mean -", round(mean(detailed_taxonomy$obitag_bestid100),2)),
                     x= "Similarity to Closest Database sequence (%)")+ 
  scale_y_continuous(labels = scales::comma)
print(plot)
plot_function(plot, file.path(output_figures, "03_Sequencing_Length_Distribution"))

#Correlation between Amplicon length & Similarity to closest Database Hit
plot <- ggplot(detailed_taxonomy,aes(x= SeqLength, y= obitag_bestid100)) + 
  geom_point( color= "deeppink4", alpha=0.5)+
  geom_smooth(method=lm , color="grey25", fill="#8B8989", se=TRUE, alpha=0.9) +
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01, geom = "label",  
           label.x.npc = "middle",label.y.npc = "bottom")+
  ylim(min(detailed_taxonomy$obitag_bestid100), max(detailed_taxonomy$obitag_bestid100))+
  plot_theme+
  theme(plot.title = element_text(face = "bold"))+
  labs(title= "Correlation between Amplicon length & Similarity to closest Database Hit",
        x= "Sequence length (bp)" ,
        y= "Similarity to Closest Database sequence (%)")
print(plot)
plot_function(plot, file.path(output_figures, "04_Correlation_SeqLength_AssignSuccess"))

#Disturibution of Assignment Success stacked
plot <- detailed_taxonomy %>% mutate(identity= round(obitag_bestid100, -1))%>%
  ggplot( aes(x=factor(obitag_rank,levels=obirank_level),fill=as.factor(identity),group=identity)) + 
  geom_bar(color = "black", alpha = 0.7) +
  plot_theme+ 
  labs(title= "Distribution of Assignment Success", x= "Taxonomic rank", fill="~% Similarity to RefDB")+
  theme(axis.text.x = element_text(angle = 45,hjust=1))+
  scale_y_continuous(labels = scales::comma)+
  scale_fill_manual(values=color_df_assign$color, breaks = color_df_assign$type)
print(plot)
plot_function(plot, file.path(output_figures, "05_Distribution_AssignmentSuccess_detailed"))


#"Disturibution of Assignment Success
plot <- ggplot(detailed_taxonomy, aes(x=factor(obitag_rank,levels=obirank_level))) + 
  geom_bar(color = "black", alpha = 0.7) +
  theme_light()+ 
  labs(title= "Distribution of Assignment Success", x= "Taxonomic rank")+
  theme(axis.text.x = element_text(angle = 45,hjust=1))+
  scale_y_continuous(labels = scales::comma)
print(plot)
plot_function(plot, file.path(output_figures, "05_Distribution_AssignmentSuccess"))


#Disturibution of Amplicon Length without primer pair sequence
plot <- ggplot(detailed_taxonomy, aes(x=SeqLength)) + 
  geom_histogram(fill="deeppink4", alpha=0.6) +
  plot_theme+
  labs(title= "Distribution of Amplicon Length without primer pair sequence", x= "Amplicon length (bp)")+
  theme(axis.text.x = element_text(angle = 45,hjust=1))+
  scale_y_continuous(labels = scales::comma)
print(plot)
plot_function(plot, file.path(output_figures, "06_Distribution_SeqLength"))


##############
## First plots about Metabarcoding success and quality
## Reads
plot <- ggplot(community, aes(x=as.factor(type), y=total.reads, fill=type)) + 
  geom_boxplot( alpha=0.5)+
  #stats
  stat_compare_means(label = "p.signif",  method="t.test",ref.group = ".all.", geom="label")  +
  #theme
  labs(x = "Sample Type", y="Total Read no.", 
       title = "Control Overview - Reads per replicate") +
  plot_theme+
  theme(axis.text.x = element_text(angle = 45,hjust=1))+ 
  scale_y_continuous(labels = scales::comma)+
  scale_fill_manual(values=color_df$color, breaks = color_df$type)
plot <- ggdraw(add_sub(plot, "Asteriks represent mean multiple pairwise tests against all", size=10))
print(plot)
plot_function(plot, file.path(output_figures, "07_ReadNo_control_overview_sampleType"))


plot <- community%>%filter( mapply(grepl, grep_negative_controls, type, perl=T))%>%
  ggplot(aes(x=as.factor(type), y=total.reads, fill=type)) + 
  geom_boxplot( alpha=0.5)+
  #stats
  stat_compare_means(label = "p.signif",  method="t.test",ref.group = ".all.", geom="label")+
  labs(x = "Type", y="Total Read no.", 
       title = "Negative Control Overview - Reads",
       subtitle="Reads per relicate splitted after Sample Type")+
  plot_theme+ 
  theme(axis.text.x = element_text(angle = 45,hjust=1))+ 
  scale_y_continuous(labels = scales::comma)+
  scale_fill_manual(values=color_df$color, breaks = color_df$type)
print(plot)
plot_function(plot, file.path(output_figures, "07_ReadNo_negative_control_overview_sampleType"))


### ASVS
plot <- ggplot(community, aes(x=as.factor(type), y=total.asvs, fill=type)) + 
  geom_boxplot( alpha=0.5) +
  stat_compare_means(label = "p.signif",  method="t.test",ref.group = ".all.", geom="label")  +
  labs(x = "Sample Type", y="Total ASV no.", 
       title = "Control Overview - ASVs ",
       subtitle = "Asteriks represent T-Test against all others") +
  plot_theme+ 
  theme(axis.text.x = element_text(angle = 45,hjust=1))+
  scale_y_continuous(labels = scales::comma)+
  scale_fill_manual(values=color_df$color, breaks = color_df$type)
print(plot)
plot_function(plot, file.path(output_figures, "07_ASVNo_control_overview_sampleType"))

plot <- community%>%filter( mapply(grepl, grep_negative_controls, type, perl=T))%>%
  ggplot(aes(x=as.factor(type), y=total.asvs, fill=type)) + 
  geom_boxplot( alpha=0.5)+
  #stats
  stat_compare_means(label = "p.signif",  method="t.test",ref.group = ".all.", geom="label")+
  labs(x = "Type", y="Total Read no.", 
       title = "Negative Control Overview - Reads",
       subtitle="Reads per relicate splitted after Sample Type")+
  plot_theme+ 
  theme(axis.text.x = element_text(angle = 45,hjust=1))+ 
  scale_y_continuous(labels = scales::comma)+
  scale_fill_manual(values=color_df$color, breaks = color_df$type)
print(plot)
plot_function(plot, file.path(output_figures, "07_ASVNo_negative_control_overview_sampleType"))



## Checking replication control
plot <- ggplot(community, aes(x=as.factor(repl_type), y=total.reads)) + 
  geom_boxplot(fill="deeppink4", alpha=0.5)+
  labs(x = "Type", y="Total Read no.", 
       title = "Plate Control Overview - Reads",
       subtitle = "Total Read number distribution of all samples for every replicate") +
  plot_theme+
  theme(axis.text.x = element_text(angle = 45,hjust=1))+ 
  scale_y_continuous(labels = scales::comma)
print(plot)
plot_function(plot, file.path(output_figures, "07_ReadNo_control_overview_plate"))


plot <- community%>%filter( mapply(grepl, grep_negative_controls, type, perl=T))%>%
  ggplot(aes(x=as.factor(repl_type), y=total.reads, fill= type)) + 
  geom_boxplot( alpha=0.5)+
  labs(x = "Type", y="Total Read no.", 
       title = "Negative Control Overview - Total Read number distribution for every replicate") +
  plot_theme+
  theme(axis.text.x = element_text(angle = 45,hjust=1))+ 
  scale_y_continuous(labels = scales::comma)+
  scale_fill_manual(values=color_df$color, breaks = color_df$type)
plot <- ggdraw(add_sub(plot, "P represents Plates of PCR & Multiplex controls, R replciates of all other sample types", size=10))
print(plot)
plot_function(plot, file.path(output_figures, "07_ReadNo_negative_control_overview_plate"))


plot <- ggplot(community, aes(x=as.factor(repl_type), y=total.asvs)) + 
  geom_boxplot(fill="deeppink4", alpha=0.5) +
  labs(x = "Type", y="Total ASVs. no.", 
       title = "Plate Control Overview - Total ASV number distribution for every replicate",
      subtitle = "Total ASV number distribution of all samples for every replicate") +
  plot_theme+
  scale_y_continuous(labels = scales::comma)
plot <-ggdraw(add_sub(plot, "P represents Plates of PCR & Multiplex controls, R replciates of all other sample types", size=10))
print(plot)
plot_function(plot, file.path(output_figures, "07_ASVNo_control_overview_plate"))

plot <- community%>%filter( mapply(grepl, grep_negative_controls, type, perl=T))%>%
  ggplot( aes(x=as.factor(repl_type), y=total.asvs)) + 
  geom_boxplot(fill="deeppink4", alpha=0.5) +
  labs(x = "Type", y="Total ASVs. no.", 
       title = "Plate Control Overview - Total ASV number distribution for every replicate",
       subtitle = "Total ASV number distribution of all samples for every replicate") +
  plot_theme+
  scale_y_continuous(labels = scales::comma)
plot <-ggdraw(add_sub(plot, "P represents Plates of PCR & Multiplex controls, R replciates of all other sample types", size=10))
print(plot)
plot_function(plot, file.path(output_figures, "07_ASVNo_negative_control_overview_plate"))


##############
#Frequency distribution of read numbers in **ASVs**

# 1) Summed reads per ASV (across selected samples/columns)
reads_per_seq_variant <- apply(community[, index], 2, sum)

# 2) Thresholds (as proportions of the total reads)
summ_of_reads <- sum(reads_per_seq_variant)
thr_props  <- c(1e-5, 5e-6, 2e-6, 1e-6, 5e-7)
thr_vals   <- summ_of_reads * thr_props
thr_colors <- c("red", "darkred", "blue", "darkblue", "darkgreen")
thr_labs   <- format(thr_props, scientific = FALSE)

thr_df <- data.frame(x = thr_vals, label = thr_labs, thr_colors = thr_colors, props=thr_props)

# 3) Frequency table (how many ASVs have N reads)
freq_table <- as.data.frame(table(reads_per_seq_variant))
names(freq_table) <- c("reads", "freq")
freq_table$reads <- as.numeric(as.character(freq_table$reads))

# Remove zeros for log scales
freq_plot <- subset(freq_table, reads > 0 & freq > 0)

# 4) Plot

plot <- ggplot(freq_plot, aes(x = reads, y = freq)) +
  geom_point(alpha = 0.6, size = 2) +
  # vertical threshold lines
  geom_vline(data = thr_df, aes(xintercept = x, colour = thr_labs), linewidth = 0.7, show.legend = FALSE) +
  # text labels for thresholds
  annotate("text",
           x = mean(freq_plot$reads) * 3, y =  max(freq_plot$freq) * c(0.90, 0.70, 0.55, 0.43, 0.30),
           label = paste(thr_df$props,round(thr_df$x, 1), sep="/~") , hjust = 1, vjust = 1, size = 3.2,
           color = rev(thr_df$thr_colors)) +
  scale_color_manual(values = setNames(thr_df$thr_colors, thr_df$thr_labs)) +
  scale_x_log10() +  scale_y_log10() +plot_theme+
  labs(
    title = "Frequency distribution of read numbers in ASVs",
    subtitle = "Lines & colors represent different proportions of total read number/actual read number",
    x = "Read count of an ASV",
    y = "Frequency"
  )+theme(legend.position = "bottom")

print(plot)
plot_function(plot, file.path(output_figures, "07_Low_frequency_ASV_detection"))


##############
#Occupancy of ASVs in **replicates**
locations_per_seq_variant <- apply(community[,index]!= 0 , 2, sum)  #remove rows with any zero
freq_loc_table <- data.frame(table(locations_per_seq_variant))
names(freq_loc_table) <- c("occ", "freq")
freq_loc_table$occ  <- as.numeric(as.character(freq_loc_table$occ))
freq_loc_table$freq <- as.numeric(freq_loc_table$freq)

# Remove zeros for log scales
freq_plot2 <- subset(freq_loc_table, occ > 0 & freq > 0)

# 3) Plot
plot <- ggplot(freq_plot2, aes(x = occ, y = freq)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_x_log10() +  scale_y_log10() + plot_theme+
  labs(
    title = "Occupancy distribution of ASVs across replicates",
    subtitle = "How often each ASV is represented",
    x = "Number of replicates in which an ASV occurs",
    y = "Frequency of ASVs"
  )
print(plot)
plot_function(plot, file.path(output_figures, "07_Replicate_number_ASV"))

rm(locations_per_seq_variant,freq_loc_table)

##############
#Frequency distribution 
plot <-ggplot(community, aes(x=as.numeric(total.reads), fill=type)) + #
  geom_histogram(  color="#e9ecef", alpha=0.8, position = 'identity') + # binwidth=0.1,
  labs(fill="", title= "Frequency distribution of Total Read no. across replicates", 
       x= "Total read no.", y="Freq")+
  plot_theme+
  scale_x_log10(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma)+
  scale_fill_manual(values=color_df$color, breaks = color_df$type)
print(plot)
plot_function(plot, file.path(output_figures, "08_ReadNo_frequency_distribution"))

medians <- community %>%
  group_by(type) %>%
  summarise(median_reads = median(as.numeric(total.reads), na.rm = TRUE))

plot <- plot +
  geom_vline(data = medians, 
             aes(xintercept = median_reads, color = type),
             linetype = "dashed", linewidth = 1) +
  scale_color_manual(values = color_df$color, breaks = color_df$type)+
  labs(subtitle ="Line represents median of the reads groupd after reads", color="")

print(plot)
plot_function(plot, file.path(output_figures, "08_ReadNo_frequency_replicates_median"))

plot <- ggplot(community, aes(x=total.asvs, fill=type)) + #
  geom_histogram(  color="#e9ecef", alpha=0.8, position = 'identity') + # binwidth=30,
  labs(fill="", title= "Frequency distribution of ASV no. across replicates", 
       x= "Total ASV no.", y="Freq")+
  plot_theme+
  scale_x_log10(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma)+
  scale_fill_manual(values=color_df$color, breaks = color_df$type)
print(plot)
plot_function(plot, file.path(output_figures, "08_ASVNo_frequency_distribution"))

#Frequency distribution  splitted after Sample Type
plot <- ggplot(community, aes(x=as.numeric(total.reads), fill=type)) + #
  geom_histogram(  color="#e9ecef", alpha=0.8, position = 'identity') + # binwidth=0.1,
  labs(fill="", title= "Frequency distribution of Read no. across replicates after Sample Type", 
       x= "Total read no.", y="Freq")+
  facet_wrap(~type, scales = "free_y")+theme_light()+ 
  theme(legend.position="none")+ 
  scale_x_log10(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma)+
  scale_fill_manual(values=color_df$color, breaks = color_df$type)
print(plot)
plot_function(plot, file.path(output_figures, "08_ReadNo_frequency_distribution_splitted"))

plot <- ggplot(community, aes(x=total.asvs, fill=type)) + #
  geom_histogram(  color="#e9ecef", alpha=0.8, position = 'identity') + # binwidth=0.1,
  labs(fill="", title= "Frequency distribution of ASV no. across replicates after Sample Type", 
       x= "Total read no.", y="Freq")+
  facet_wrap(~type, scales = "free_y")+
  plot_theme+
  theme(legend.position="none")+ 
  scale_y_continuous(labels = scales::comma)+
  scale_x_log10(labels = scales::comma) +
  scale_fill_manual(values=color_df$color, breaks = color_df$type)
print(plot)
plot_function(plot, file.path(output_figures, "08_ASVNo_frequency_distribution_splitted"))

# Richness ~ reads
plot <- ggplot(community, aes(x= total.reads, y= total.asvs, color= type)) + 
  geom_point() + 
  coord_trans(x="log2", y="log2") + theme_classic()+ 
  labs(title = "Total ASVs vs. Total Read no. per replicates grouped by Sample Type",
       x="Total Read Number", y= "Total ASV number",color="Sample Type")+ 
  plot_theme+
  theme(legend.position="bottom")+ 
  scale_y_continuous(labels = scales::comma)+ 
  scale_x_continuous(labels = scales::comma,
                     guide = guide_axis(angle = 45))+
  scale_color_manual(values=color_df$color, breaks = color_df$type)
print(plot)
plot_function(plot, file.path(output_figures, "09_ASV_vs_readno_per_replicate"))

##############
# Ordination - 
set.seed(25)

### caluclate distace
rownames(community) <- community$full_replicate

bray_data <- metaMDS(community[,index], distance ="bray", parallel =threads)
nmds <- plot_simpleNDMS_symbol(bray_data, community, "type")

nmds[[2]] <- nmds[[2]] + 
  labs(title = paste(pr, "Similarity between all replicates ", sep=" - "),
       subtitle=" Bray Curtis Distance grouped after Sample Type (seed no. 25)",
       color="Sample Type")+ 
  plot_theme+
  theme(legend.position="bottom")+
  scale_color_manual(values=color_df$color, breaks = color_df$type)
print(nmds[[2]])
plot_function(nmds[[2]], file.path(output_figures, "10_NMDS_ordination_after_SampleType"))


##############
## remove low read samples
sample_reads <- rowSums(community[,index])
sample_reads <- sample_reads > max(sample_reads)*0.01

bray_data <- metaMDS(community[sample_reads,index], distance ="bray", parallel =threads)

nmds <- plot_simpleNDMS_symbol(bray_data, community[sample_reads,], "type")

##normal ordination plot
nmds[[2]] <- nmds[[2]] + 
  labs(title = "Community Similarity between all replicates (seed no. 25 ) ",
       subtitle= "Low abundant replicates excluded, Bray Curtis Distance")+ 
  plot_theme+
  scale_color_manual(values=color_df$color, breaks = color_df$type)

nmds[[2]] <- ggdraw(add_sub(nmds[[2]], "Low abundant replicates: replicates with less than 1% reads of max sample read number\n Colors represent samples", size=10))
plot_function(nmds[[2]], file.path(output_figures, "10_NMDS_ordination_after_SampleType_filtered"))

##### ordination with spiders for samples not controls
#prepare color for ghibli
sample_no <- length(unique(community[index_samples,"sample"]))
color_samples <-as.vector(ghibli_palette(name = "PonyoMedium", n = sample_no, type = "continuous"))

bray_data <- metaMDS(community[index_samples,index], distance ="bray")

nmds <- plot_simpleNDMS_symbol(bray_data, community[index_samples,], "sample")
##normal ordination plot
nmds[[2]] <- nmds[[2]] + 
  labs(title = "Community Similarity between the Sample replicates (seed no. 25)",
       subtitle= "Controls were excluded, Bray Curtis Distance")+ 
  plot_theme+
  theme(legend.position="none")+
  scale_color_manual(values=color_samples)
print(nmds[[2]])
plot_function(nmds[[2]], file.path(output_figures, "10_NMDS_ordination_after_SampleType_noControls"))

##spider plot
nmds <- plot_simpleNDMS_spider(bray_data, community[index_samples,], "sample")

nmds[[2]] <- nmds[[2]] + 
  labs(title = "Community Similarity between the Sample replicates (seed no. 25) ",
       subtitle= "Controls were excluded, Bray Curtis Distance")+ 
  plot_theme+
  theme(legend.position="none")+
  scale_color_manual(values=color_samples)
print(nmds[[2]])
plot_function(nmds[[2]], file.path(output_figures, "10_NMDS_spiderplot_after_SampleType_noControls"))

##############
## remove low read samples
sample_reads <- rowSums(community[index_samples,index])
sample_reads <- sample_reads > max(sample_reads)*0.01

bray_data <- metaMDS(community[index_samples[sample_reads],index], distance ="bray", parallel =threads)

nmds <- plot_simpleNDMS_symbol(bray_data, community[index_samples,], "sample")

##normal ordination plot
nmds[[2]] <- nmds[[2]] + 
  labs(title = "Community Similarity between the Sample replicates (seed no. 25 ) ",
       subtitle= "Controls were excluded & low abundant replicates, Bray Curtis Distance")+ 
  plot_theme+
  theme(legend.position="none")+
  scale_color_manual(values=color_samples)

nmds[[2]] <- ggdraw(add_sub(nmds[[2]], "Low abundant replicates: replicates with less than 1% reads of max sample read number\n Colors represent samples", size=10))
plot_function(nmds[[2]], file.path(output_figures, "10_NMDS_ordination_after_SampleType_noControls_filtered"))

##spider plot
nmds <- plot_simpleNDMS_spider(bray_data, community[index_samples,], "sample")

nmds[[2]] <- nmds[[2]] + 
  labs(title = "Community Similarity between the Sample replicates (seed no. 25) ",
       subtitle= "Controls were excluded & low abundant replicates, Bray Curtis Distance")+ 
  plot_theme+
  theme(legend.position="none")+
  scale_color_manual(values=color_samples)
nmds[[2]]<- ggdraw(add_sub(nmds[[2]], "Low abundant replicates: replicates with less than 1% reads of max sample read number\n Colors represent samples", size=10))
plot_function(nmds[[2]], file.path(output_figures, "10_NMDS_spiderplot_after_SampleType_noControls_noControls_filtered"))

rm(bray_data, nmds)

dev.off()
