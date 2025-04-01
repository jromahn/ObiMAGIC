#!/usr/bin/env Rscript

########################################
# Diagnosis script fro ObiTools4
# Juliane Romahn
# January 14, 2025
# html_document
# execute via: 03_ObiWizard_diagnostic.R Resuults_folder FINAL_fasta_file PATH/to/taxdump threads full/path/to/00_ObiScripts/00_OBIMAGIC_functions.R
########################################

#rm(list=ls())  
#install.packages("jsonify")
#install.packages("devtools")
 # install.packages("ghibli")
#devtools::install_git("https://git.metabarcoding.org/obitools/obitools4/robireadfasta.git")
#install.packages('taxize', repos = c('https://ropensci.r-universe.dev', 'https://cloud.r-project.org'))
library(ROBIFastread)
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


##############
## Input
args = commandArgs(trailingOnly=TRUE)
path=args[1] #"Malawi_2023_eDNA_results/"
input_file =args[2] #"7_Malawi_2023_eDNA_final.fasta" #
taxdump_path=args[3] 
threads= args[4]
functions_path=args[5]

#### read in all functions
source(functions_path)


output=gsub(".fasta", "", input_file)

##############
# convert ObiTools Output for Mothur
data <- read_obifasta(file.path(path,input_file))
data <- unique(data)
Large <- extract_features(data,key="count")
tab <- extract_readcount(data,key="merged_sample")
community <- as.data.frame(as.matrix(tab))

old_asv_names <- colnames(community)
new_asv_names <- paste(asv_start,sprintf("%06d", c(1:length(colnames(community)))), sep="_")
colnames(community) <- new_asv_names


########### prepare taxonomy

### get taxonomy & save it in same folder then the taxdump
loc_ncbi <- c(paste(taxdump_path,  "names.dmp", sep="/"),
              paste(taxdump_path,  "nodes.dmp", sep="/")) 
sqlFile <- file.path(path, paste(output,'accessionTaxa.sql', sep="_"))
read.accession2taxid(list.files(loc_ncbi),sqlFile,indexTaxa=TRUE,overwrite=TRUE)
taxaNames<-read.names.sql(loc_ncbi[1],sqlFile)
taxaNodes<-read.nodes.sql(loc_ncbi[2],sqlFile)


########  get taxonomy
taxonomy <- extract_features(data,key=c("obitag_bestmatch","obitag_bestid", "obitag_rank", "obitag_similarity_method","scientific_name", "taxid"))
taxonomy <- taxonomy %>%
      select(-features, - definition)%>%
      mutate(obitag_bestid100 = obitag_bestid *100)%>%
      mutate(SeqLength=nchar(sequence))%>%
      relocate(SeqLength, .after=sequence)

#depending on the obitools version extract scientific name and taxid
if (!"taxid" %in% colnames(taxonomy)) {
  #scientific_name looks like this: taxon:2759 [Eukaryota]@superkingdom
  taxonomy <- taxonomy %>% 
    mutate(taxid= gsub("taxon:(\\d+).*", "\\1", taxonomy$scientific_name, perl=T))%>% 
    mutate(scientific_name= gsub("taxon:\\d+.*\\[(.*)\\].*", "\\1", taxonomy$scientific_name, perl=T))
} 


# feedback by obitools 4  
taxonomy <- merge(taxonomy, data.frame(id=old_asv_names, id_new = new_asv_names), by = "id")
write.table(taxonomy, file=file.path(path, paste(output, "taxonomy_info.tsv", sep="__")), quote=FALSE, sep = "\t", row.names = F)
  
# taxonomy of ncbi taxon ids  
detailed_taxonomy <- ncbi_taxonomy(taxonomy$taxid)
write.table(detailed_taxonomy, file=file.path(path, paste(output, "taxonomy_SPECIES_summary.tsv", sep="__")), quote=FALSE, sep = "\t", row.names = F)
  
# combined overview of sequences and their taxonomy
detailed_taxonomy <- merge( taxonomy, detailed_taxonomy,by.x = "taxid", by.y="taxID")
write.table(detailed_taxonomy, file=file.path(path, paste(output, "taxonomy_info_COMPLETE.tsv", sep="__")), quote=FALSE, sep = "\t", row.names = F)


##################

##############
## reformat data and produce files for follow up steps
Large <- merge(Large, data.frame(id=old_asv_names, id_new = new_asv_names), by = "id")
Large$definition <-NULL
counttable <- as.data.frame(t(community))

total <- rowSums(counttable)
counttable <- cbind(total, counttable)
Representative_Sequence <- rownames(counttable)
counttable <- cbind(Representative_Sequence, counttable)

write.fasta(sequences = as.list(Large$sequence), names= Large$id_new, file.path(path, paste(output, "renamed.fasta", sep="__")), 
            open = "w", nbchar = nchar(max(Large$sequence)), as.string = FALSE)
#write.table(Large, file=file.path(path, paste(output, "general_infos.tsv", sep="__")), quote=FALSE, sep = "\t", row.names = F)

save(community, file=file.path(path, paste(output, "community.RData", sep="__")))

rm(Large,data,counttable)

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
ggplot(community_sum, aes(area = replicates, fill = type, label= type, subgroup=type_long)) +
  geom_treemap()+
  geom_treemap_subgroup_border(color="white", size=4)+
  geom_treemap_subgroup_text(place ="center",color="white",reflow = T)+
  geom_treemap_text(color="black", alpha = 0.7)+
  scale_fill_manual(values=color_df$color, breaks = color_df$type)+
  labs(title= paste("Project overview:",output  , "\nPrimer:", pr),
       subtitle = "Size represents replicate number")+
  plot_theme  


#print it as table
ggplot() + theme_void() + annotate(geom="table",x=1,y=1,label=list(community_sum),size=3)+
  labs(title= paste("Project overview:",output ),
       subtitle = "Total includes multiplexed primers (if existing)")

# colors assignment
pal <- as.vector(ghibli_palette(name = "MarnieMedium2", n = 10, type = "continuous"))
color_df_assign <- data.frame(type=seq(10, 100, by = 10), color=pal)


# Distribution of Similarity to closest Database Hit
ggplot(detailed_taxonomy, aes(x=obitag_bestid100)) + 
  geom_histogram(fill="deeppink4", alpha=0.8)+
  geom_vline(xintercept = mean(detailed_taxonomy$obitag_bestid100), linetype="dotted", color = "black", linewidth=1.5)+
  plot_theme+
  labs(title= "Distribution of Similarity to closest Database Hit", 
                     subtitle= paste("Dotted line represents mean -", round(mean(detailed_taxonomy$obitag_bestid100),2)),
                     x= "Similarity to Closest Database sequence (%)")+ 
  scale_y_continuous(labels = scales::comma)

#Correlation between Amplicon length & Similarity to closest Database Hit
ggplot(detailed_taxonomy,aes(x= SeqLength, y= obitag_bestid100)) + 
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

#Disturibution of Assignment Success stacked
detailed_taxonomy %>% mutate(identity= round(obitag_bestid100, -1))%>%
  ggplot( aes(x=factor(obitag_rank,levels=obirank_level),fill=as.factor(identity),group=identity)) + 
  geom_histogram(color="black", alpha=0.7,stat="count") +
  plot_theme+ 
  labs(title= "Distribution of Assignment Success", x= "Taxonomic rank", fill="~% Similarity to RefDB")+
  theme(axis.text.x = element_text(angle = 45,hjust=1))+
  scale_y_continuous(labels = scales::comma)+
  scale_fill_manual(values=color_df_assign$color, breaks = color_df_assign$type)

#"Disturibution of Assignment Success
ggplot(detailed_taxonomy, aes(x=factor(obitag_rank,levels=obirank_level))) + 
  geom_histogram(fill="deeppink4", alpha=0.6,stat="count") +
  theme_light()+ 
  labs(title= "Distribution of Assignment Success", x= "Taxonomic rank")+
  theme(axis.text.x = element_text(angle = 45,hjust=1))+
  scale_y_continuous(labels = scales::comma)


#Disturibution of Amplicon Length without primer pair sequence
ggplot(detailed_taxonomy, aes(x=SeqLength)) + 
  geom_histogram(fill="deeppink4", alpha=0.6) +
  plot_theme+
  labs(title= "Distribution of Amplicon Length without primer pair sequence", x= "Amplicon length (bp)")+
  theme(axis.text.x = element_text(angle = 45,hjust=1))+
  scale_y_continuous(labels = scales::comma)


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
ggdraw(add_sub(plot, "Asteriks represent mean multiple pairwise tests against all", size=10))


community%>%filter( mapply(grepl, grep_negative_controls, type, perl=T))%>%
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

### ASVS
ggplot(community, aes(x=as.factor(type), y=total.asvs, fill=type)) + 
  geom_boxplot( alpha=0.5) +
  stat_compare_means(label = "p.signif",  method="t.test",ref.group = ".all.", geom="label")  +
  labs(x = "Sample Type", y="Total ASV no.", 
       title = "Control Overview - ASVs ",
       subtitle = "Asteriks represent T-Test against all others") +
  plot_theme+ 
  theme(axis.text.x = element_text(angle = 45,hjust=1))+
  scale_y_continuous(labels = scales::comma)+
  scale_fill_manual(values=color_df$color, breaks = color_df$type)


## Checking replication control
ggplot(community, aes(x=as.factor(repl_type), y=total.reads)) + 
  geom_boxplot(fill="deeppink4", alpha=0.5)+
  labs(x = "Type", y="Total Read no.", 
       title = "Plate Control Overview - Reads",
       subtitle = "Total Read number distribution of all samples for every replicate") +
  plot_theme+
  theme(axis.text.x = element_text(angle = 45,hjust=1))+ 
  scale_y_continuous(labels = scales::comma)

plot <- community%>%filter( mapply(grepl, grep_negative_controls, type, perl=T))%>%
  ggplot(aes(x=as.factor(repl_type), y=total.reads, fill= type)) + 
  geom_boxplot( alpha=0.5)+
  labs(x = "Type", y="Total Read no.", 
       title = "Negative Control Overview - Total Read number distribution for every replicate") +
  plot_theme+
  theme(axis.text.x = element_text(angle = 45,hjust=1))+ 
  scale_y_continuous(labels = scales::comma)+
  scale_fill_manual(values=color_df$color, breaks = color_df$type)
ggdraw(add_sub(plot, "P represents Plates of PCR & Multiplex controls, R replciates of all other sample types", size=10))


plot <- ggplot(community, aes(x=as.factor(repl_type), y=total.asvs)) + 
  geom_boxplot(fill="deeppink4", alpha=0.5) +
  labs(x = "Type", y="Total ASVs. no.", 
       title = "Plate Control Overview - Total ASV number distribution for every replicate",
      subtitle = "Total ASV number distribution of all samples for every replicate") +
  plot_theme+
  scale_y_continuous(labels = scales::comma)
ggdraw(add_sub(plot, "P represents Plates of PCR & Multiplex controls, R replciates of all other sample types", size=10))


##############
#Frequency distribution of read numbers in **ASVs**

reads_per_seq_variant <- apply(community[,index],2,sum)
# caluclate threshold:
summ_of_reads <- sum(reads_per_seq_variant)
rarity_threshold1 <- summ_of_reads*0.00001
rarity_threshold2 <- summ_of_reads*0.000005
rarity_threshold3 <- summ_of_reads*0.000002
rarity_threshold4 <- summ_of_reads*0.000001
rarity_threshold5 <- summ_of_reads*0.0000005

#plot frequen
freq_table = data.frame(table(reads_per_seq_variant))
x_max<-mean(as.numeric(freq_table$reads_per_seq_variant))*3
y_max<-max(as.numeric(freq_table$Freq))


freq_table$reads_per_seq_variant <- 
  as.vector(freq_table$reads_per_seq_variant)
plot(freq_table, log = c("xy"),
     main="Frequency distribution of read numbers in ASVs\n How abundanted is each ASV?",
     sub= "Note: Lines & colors represent different proportion of total read number",
     xlab = c("Read count in ASV"),
     ylab= "Freq/Abu") +
  abline(v =rarity_threshold1, col="red")+
  text(x_max,y_max*0.9, "0.00001", col="red")+
  abline(v =rarity_threshold2, col="darkred")+
  text(x_max,y_max*0.70, "0.000005", col="darkred")+
  abline(v =rarity_threshold3, col="blue")+
  text(x_max,y_max*0.55, "0.000002", col="blue")+
  abline(v =rarity_threshold4, col="darkblue")+
  text(x_max,y_max*0.43, "0.000001", col="darkblue")+
  abline(v =rarity_threshold5, col="darkgreen")+
  text(x_max,y_max*0.30, "0.0000005", col="darkgreen")

##############
#Occupancy of ASVs in **replicates**

locations_per_seq_variant <- apply(community[,index]!= 0 , 2, sum)  #remove rows with any zero
freq_loc_table <- data.frame(table(locations_per_seq_variant))
freq_loc_table$locations_per_seq_variant <- as.vector(freq_loc_table$locations_per_seq_variant)
plot(freq_loc_table,log = c("xy"), 
     main="Frequency distribution of read numbers in ASVs\n How abundanted is each ASV?",
     sub= "Note: Lines & colors represent different proportion of total read number",
     xlab = c("Read count of an ASV"),
     ylab= "Freq/Abu") +
rm(locations_per_seq_variant,freq_loc_table)

##############
#Frequency distribution 
ggplot(community, aes(x=as.numeric(total.reads), fill=type)) + #
  geom_histogram(  color="#e9ecef", alpha=0.8, position = 'identity') + # binwidth=0.1,
  labs(fill="", title= "Frequency distribution of Total Read no. across replicates", 
       x= "Total read no.", y="Freq")+
  plot_theme+
  scale_x_log10(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma)+
  scale_fill_manual(values=color_df$color, breaks = color_df$type)

ggplot(community, aes(x=total.asvs, fill=type)) + #
  geom_histogram(  color="#e9ecef", alpha=0.8, position = 'identity') + # binwidth=30,
  labs(fill="", title= "Frequency distribution of ASV no. across replicates", 
       x= "Total ASV no.", y="Freq")+
  plot_theme+
  scale_x_log10(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma)+
  scale_fill_manual(values=color_df$color, breaks = color_df$type)

#Frequency distribution  splitted after Sample Type
ggplot(community, aes(x=as.numeric(total.reads), fill=type)) + #
  geom_histogram(  color="#e9ecef", alpha=0.8, position = 'identity') + # binwidth=0.1,
  labs(fill="", title= "Frequency distribution of Read no. across replicates after Sample Type", 
       x= "Total read no.", y="Freq")+
  facet_wrap(~type, scales = "free_y")+theme_light()+ 
  theme(legend.position="none")+ 
  scale_x_log10(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma)+
  scale_fill_manual(values=color_df$color, breaks = color_df$type)

ggplot(community, aes(x=total.asvs, fill=type)) + #
  geom_histogram(  color="#e9ecef", alpha=0.8, position = 'identity') + # binwidth=0.1,
  labs(fill="", title= "Frequency distribution of ASV no. across replicates after Sample Type", 
       x= "Total read no.", y="Freq")+
  facet_wrap(~type, scales = "free_y")+
  plot_theme+
  theme(legend.position="none")+ 
  scale_y_continuous(labels = scales::comma)+
  scale_x_log10(labels = scales::comma) +
  scale_fill_manual(values=color_df$color, breaks = color_df$type)

# Richness ~ reads
ggplot(community, aes(x= total.reads, y= total.asvs, color= type)) + 
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

ggdraw(add_sub(nmds[[2]], "Low abundant replicates: replicates with less than 1% reads of max sample read number\n Colors represent samples", size=10))


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

##spider plot
nmds <- plot_simpleNDMS_spider(bray_data, community[index_samples,], "sample")

nmds[[2]] <- nmds[[2]] + 
  labs(title = "Community Similarity between the Sample replicates (seed no. 25) ",
       subtitle= "Controls were excluded, Bray Curtis Distance")+ 
  plot_theme+
  theme(legend.position="none")+
  scale_color_manual(values=color_samples)
print(nmds[[2]])

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

ggdraw(add_sub(nmds[[2]], "Low abundant replicates: replicates with less than 1% reads of max sample read number\n Colors represent samples", size=10))

##spider plot
nmds <- plot_simpleNDMS_spider(bray_data, community[index_samples,], "sample")

nmds[[2]] <- nmds[[2]] + 
  labs(title = "Community Similarity between the Sample replicates (seed no. 25) ",
       subtitle= "Controls were excluded & low abundant replicates, Bray Curtis Distance")+ 
  plot_theme+
  theme(legend.position="none")+
  scale_color_manual(values=color_samples)
ggdraw(add_sub(nmds[[2]], "Low abundant replicates: replicates with less than 1% reads of max sample read number\n Colors represent samples", size=10))

rm(bray_data, nmds)

dev.off()
