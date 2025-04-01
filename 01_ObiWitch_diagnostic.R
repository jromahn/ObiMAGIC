#!/usr/bin/env Rscript

########################################
# Diagnosis script fro ObiTools4
# Juliane Romahn
# January 7, 2024
# html_document
# execute via: 01_ObiWitch_diagnostic.R PROJECT_FOLDER/ FASTA-in.Obi-Format.fasta NGSFILE.tsv threads_number full/path/to/00_ObiScripts/00_OBIMAGIC_functions.R
########################################
  
rm(list=ls())  
#install.packages("jsonify")
#install.packages("devtools")
#devtools::install_git("https://git.metabarcoding.org/obitools/obitools4/robireadfasta.git")
library(ROBIFastread)
library(tidyverse)
library(vegan)
library(magrittr)
library(seqinr)
library(ggpmisc)
library(gapminder) # allow commas in numbers --> +scale_y_continuous(labels = scales::comma)
library(ggpubr) # ggplot stats
library(ghibli) ## ggplot colors
library(cowplot) ## add text under plot
library(plotly) ##ggtreemap
library(treemapify) ##ggtreemap
options(scipen = 999) # stop scientific notation & add comma to larger numbers


#############################
## Input
args = commandArgs(trailingOnly=TRUE)
path=args[1] 
input_file =args[2] 
ngs_file=args[3]
threads=args[4]
functions_path=args[5]
project=gsub("_results","", path)

#### read in all functions 
source(functions_path)
if (length(threads)==0) { threads=1}
output=gsub(".fasta", "", input_file)
#############################


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
##################


##############
## reformat data 
Large <- merge(Large, data.frame(id=old_asv_names, id_new = new_asv_names), by = "id")
Large$definition <-NULL
counttable <- as.data.frame(t(community))

total <- rowSums(counttable)
counttable <- cbind(total, counttable)
Representative_Sequence <- rownames(counttable)
counttable <- cbind(Representative_Sequence, counttable)
##################


################## create output files
write.fasta(sequences = as.list(Large$sequence), names= Large$id_new, file.path(path, paste(output, "renamed_sequences.fasta", sep="__")), 
            open = "w", nbchar = nchar(max(Large$sequence)), as.string = FALSE)
write.table(counttable, file=file.path(path, paste(output, "mothur_counttable.csv", sep="__")), quote=FALSE, sep = "\t", row.names = F)
write.table(Large, file=file.path(path, paste(output, "general_infos.tsv", sep="__")), quote=FALSE, sep = "\t", row.names = F)
save(community, file=file.path(path, paste(output, "community.RData", sep="__")))
rm(Large,counttable)     


################################## PLOTTING ################################## 
## prepare data and sample identification

#extract replicate and sample information from communtiy rownames
community <- extra_sample_names(community, type_identifier)

# get index from asv and samples
index <- grep(asv_start, colnames(community))
index_samples <-which(!grepl(grep_controls,community$type))

community$total.reads <- rowSums(community[,index])
community$total.asvs <- rowSums(1*(community[,index] >0))


################################## NGS Info ################################## 
## prepare data and sample identification if information exists

plate_info <- read.table(file.path(path, ngs_file), sep="\t")
index_position <- grep("position", plate_info[1,])

if (length(index_position) >0 ){
  plate_info <- extract_plate_information(plate_info, index_position)
  plate_plot <- community %>% select( total.reads,total.asvs, full_replicate,type )%>%
          right_join( plate_info, by = c("full_replicate"="replicate"))%>%
          arrange(plate_no,plate_position ) 

  #print(plate_plot)
  
}

############################################# PLOT #####################################################

pdf(file = file.path(path, paste(output, "diagnosis_plots_bioinformatics.pdf", sep="__")), paper = "a4r")

### first plots are a project overview

# secure that sample types are always named the same by sorting them before
type_List <-community %>% arrange(type)%>% pull(type) %>% unique()
## define colors
selected <- colorRampPalette(colors)(length(unique(community$type)))
color_df<- data.frame(type= type_List, color=selected)

###### project overview
print("categorize")
community_sum <- categorize_community(community,asv_start)
#rint(community_sum)
#quit()
print("Project plotting")
ggplot(community_sum, aes(area = replicates, fill = type, label= type, subgroup=type_long)) +
  geom_treemap()+
  geom_treemap_subgroup_border(color="white", size=4)+
  geom_treemap_subgroup_text(place ="center",color="white",reflow = T)+
  geom_treemap_text(color="black", alpha = 0.7)+
  scale_fill_manual(values=color_df$color, breaks = color_df$type)+
  labs(title= paste("Project overview:",output , "\nPrimer:", unique(community$primer)),
       subtitle = "Size represents replicate number\nColor coding is consistent within the document")+
  plot_theme  

#print it as table
ggplot() + theme_void() + annotate(geom="table",x=1,y=1,label=list(community_sum), size=3)+
  labs(title= paste("Project overview:",output ),
       subtitle = "Total includes multiplexed primers (if existing)")

##  plots about Metabarcoding success and quality
### print read numbers in plate format if infomration exists
print("Plate plotting")
if ( length(index_position) >0 ){
    for( no in unique(plate_plot$plate_no)){
        if(no >0){ # since we include a spaceholder if these information are missing within OBIMAGIC, loop should not executed in those cases
          plot <- plate_plot %>% filter(plate_no==no)%>%
            ggplot(aes(x=plate_column, y= forcats::fct_rev(factor(plate_row)), color= type, size= total.reads)) + 
            geom_point()+
            plot_theme+
            labs(title = paste("Metabarcoding PCR Plate setup - Plate:", no), 
                y= "Rows", x ="Columns", 
                color="Sample Type", size="Total read number per well")+
            #scale_size_manual(labels = scales::comma)+
            scale_color_manual( values=color_df$color, 
                              breaks = color_df$type)+
            guides(color = guide_legend(nrow = 2))+ theme(legend.box="vertical")
          
          print(plot)
    }
  }
  rm(plate_plot)
}


################################## print sequence length ##################################
ggplot(data, aes(x=nchar(sequence))) + 
  geom_histogram(fill="deeppink4", alpha=0.6) +
  plot_theme+
  labs(title= "Distribution of Amplicon Length (without primer pair sequence)", 
       subtitle = paste(" Average amplicon length",mean(nchar(data$sequence)) , sep=" - "),
       x= "Amplicon length (bp)",
       y= "Number of ASVs")+
  theme(axis.text.x = element_text(angle = 45,hjust=1))+
  scale_y_continuous(labels = scales::comma)

rm(data)

################################## ASV ##################################
plot <- ggplot(community, aes(x=as.factor(type), y=total.asvs, fill=type)) + 
  geom_boxplot( alpha=0.5) +
  labs(x = "Type", y="Total ASV no.", 
       title = "Control Overview - ASVs number distribution",
       subtitle = "ASVs per relicate splitted after Sample Type") +
  stat_compare_means(label = "p.signif",  method="t.test",ref.group = ".all.", geom="label")  +
  plot_theme+ theme(axis.text.x = element_text(angle = 45,hjust=1))+
  scale_y_continuous(labels = scales::comma)+
  scale_fill_manual(values=color_df$color, breaks = color_df$type)
ggdraw(add_sub(plot, "Asteriks represent mean multiple pairwise tests against all", size=10))


community%>%filter( mapply(grepl, grep_negative_controls, type, perl=T))%>%
  ggplot( aes(x=as.factor(type), y=total.asvs, fill=type)) + 
    geom_violin( alpha=0.5) +
    labs(x = "Type", y="Total ASV no.", 
         title = "Negative Control Overview - ASVs number distribution",
         subtitle = "ASVs per relicate splitted after Sample Type") +
    plot_theme+ 
    scale_y_continuous(labels = scales::comma)+
  scale_fill_manual(values=color_df$color, breaks = color_df$type)

ggplot(community, aes(x=as.factor(type), y=total.asvs, fill=type)) + 
  geom_violin( alpha=0.5) +
  labs(x = "Type", y="Total ASV no.", 
       title = "Control Overview - ASVs number distribution ",
       subtitle="ASVs per relicate splitted after Sample Type") +
  facet_wrap(~type, scales= "free") + plot_theme+ 
  scale_y_continuous(labels = scales::comma)+
  scale_fill_manual(values=color_df$color, breaks = color_df$type)


################################## Reads ##################################
plot <- ggplot(community, aes(x=as.factor(type), y=total.reads, fill=type)) + 
  geom_boxplot( alpha=0.5)+
  labs(x = "Type", y="Total Read no.", 
       title = "Control Overview - Read number distribution",
       subtitle = "Reads per relicate splitted after Sample Type") +
  stat_compare_means(label = "p.signif",  method="t.test",ref.group = ".all.", geom="label")  +
  plot_theme+ theme(axis.text.x = element_text(angle = 45,hjust=1))+ 
  scale_y_continuous(labels = scales::comma)+
  scale_fill_manual(values=color_df$color, breaks = color_df$type)
ggdraw(add_sub(plot, "Asteriks represent mean multiple pairwise tests against all", size=10))


community%>%filter( mapply(grepl, grep_negative_controls, type, perl=T))%>%
  ggplot(aes(x=as.factor(type), y=total.reads, fill=type)) + 
  geom_boxplot(alpha=0.5)+
  labs(x = "Type", y="Total Read no.", 
       title = "Negative Control Overview - Read number distribution",
       subtitle="Reads per relicate splitted after Sample Type",
       fill="Sample Type")+
  plot_theme+ theme(axis.text.x = element_text(angle = 45,hjust=1))+ 
  scale_y_continuous(labels = scales::comma)+
  scale_fill_manual(values=color_df$color, breaks = color_df$type)

ggplot(community, aes(x=as.factor(type), y=total.reads, fill=type)) + 
  geom_violin(alpha=0.5) +
  labs(x = "Type", y="Total Read no.", 
       title = "Control Overview - Reads",
       subtitle = "Reads per relicate splitted after Sample Type",
       fill="Sample Type") +
  facet_wrap(~type, scales= "free") + plot_theme+ 
  scale_y_continuous(labels = scales::comma)+
  scale_fill_manual(values=color_df$color, breaks = color_df$type)

## replication control
plot <- ggplot(community, aes(x=as.factor(repl_type), y=total.reads)) + 
  geom_boxplot(fill="deeppink4", alpha=0.5)+
  labs(x = "Type", y="Total Read no.", 
       title = "Control Overview - Reads",
       subtitle = "Total Read number distribution of samples for every replicate") +
  plot_theme+ theme(axis.text.x = element_text(angle = 45,hjust=1))+ 
  scale_y_continuous(labels = scales::comma)
ggdraw(add_sub(plot, "P represents Plates of PCR & Multiplex controls, R replicates of all other sample types", size=10))


plot <- community%>%filter( mapply(grepl, grep_negative_controls, type, perl=T))%>%
  ggplot(aes(x=as.factor(repl_type), y=total.reads, fill= type)) + 
  geom_boxplot( alpha=0.5)+
  labs(x = "Type", y="Total Read no.", 
       title = "Negative Control Overview ",
       subtitle = "Total Read number distribution for every replicate") +
  plot_theme+ theme(axis.text.x = element_text(angle = 45,hjust=1))+ 
  scale_y_continuous(labels = scales::comma)+
  scale_fill_manual(values=color_df$color, breaks = color_df$type)
ggdraw(add_sub(plot, "P represents Plates of PCR & Multiplex controls, R replicates of all other sample types", size=10))


plot <- ggplot(community, aes(x=as.factor(repl_type), y=total.asvs)) + 
  geom_boxplot(fill="deeppink4", alpha=0.5) +
  labs(x = "Type", y="Total ASVs. no.", 
       title = "Control Overview - Total ASV number distribution for every replicate",
       subtitle = "P represents Plates of PCR & Multiplex controls, R replciates of all other sample types") +
  plot_theme+ scale_y_continuous(labels = scales::comma)
ggdraw(add_sub(plot, "P represents Plates of PCR & Multiplex controls, R replciates of all other sample types", size=10))

##################################
#Frequency distribution of read numbers in **ASVs**
##################################

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
     xlab = c("Read count of an ASV"),
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

#dev.off()
#quit()


##################################
#Occupancy of ASVs in **replicates**
##################################

locations_per_seq_variant <- apply(community[,index]!= 0 , 2, sum)  #remove rows with any zero
freq_loc_table <- data.frame(table(locations_per_seq_variant))
freq_loc_table$locations_per_seq_variant <- as.vector(freq_loc_table$locations_per_seq_variant)
plot(freq_loc_table,log = c("xy"), 
     main="Occupancy/Frequency distribution of ASVs across replicates\n How often is each ASV represented?", 
     xlab = c("Abundance of an ASV in replicate no."), 
     ylab= "Frequency/Abundance") 
rm(locations_per_seq_variant,freq_loc_table)

##################################
#Frequency distribution of  **replicates**
##################################

ggplot(community, aes(x=as.numeric(total.reads), fill=type)) + #
  geom_histogram(  color="#e9ecef", alpha=0.8, position = 'identity') + # binwidth=0.1,
  labs(fill="", title= "Frequency distribution of Total Read No. across replicates", 
       x= "Total read no.", y="Freq")+
  scale_x_log10(labels = scales::comma) +plot_theme + 
  scale_y_continuous(labels = scales::comma)+
  scale_fill_manual(values=color_df$color, breaks = color_df$type)

ggplot(community, aes(x=total.asvs, fill=type)) + #
  geom_histogram(  color="#e9ecef", alpha=0.8, position = 'identity') + # binwidth=30,
  labs(fill="", title= "Frequency distribution of ASV no. in replicates", 
       x= "Total ASV no.", y="Freq")+
  scale_x_log10(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma)+
  plot_theme+
  scale_fill_manual(values=color_df$color, breaks = color_df$type)

#Frequency distribution  splitted after Sample Type
ggplot(community, aes(x=as.numeric(total.reads), fill=type)) + #
  geom_histogram(  color="#e9ecef", alpha=0.8, position = 'identity') + # binwidth=0.1,
  labs(fill="", title= "Frequency distribution of Read no. in replicates after Sample Type", 
       x= "Total read no.", y="Freq")+
  scale_x_log10(labels = scales::comma) +
  facet_wrap(~type, scales = "free_y")+plot_theme+
  theme(legend.position="none")+ 
  scale_y_continuous(labels = scales::comma)+
  scale_fill_manual(values=color_df$color, breaks = color_df$type)

ggplot(community, aes(x=total.asvs, fill=type)) + #
  geom_histogram(  color="#e9ecef", alpha=0.8, position = 'identity') + # binwidth=0.1,
  labs(fill="", title= "Frequency distribution of ASV no. in replicates after Sample Type", 
       x= "Total read no.", y="Freq")+
  scale_x_log10(labels = scales::comma) +
  facet_wrap(~type, scales = "free_y")+plot_theme+  
  theme(legend.position="none")+ 
  scale_y_continuous(labels = scales::comma)+
  scale_fill_manual(values=color_df$color, breaks = color_df$type)

################################## Richness ~ reads ##################################
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


################################## Ordination ##################################
set.seed(25)

# write full replciate name into rownames otherwise bray data can not be combined with original data
rownames(community) <- community$full_replicate

## do the ordination for each primer seperated
#pr="NA"
for (pr in unique(community$primer)){
  print(paste(pr, "- Start ordination"))
  
  #subset matrix and remove asv from other primer
  subset_community <-  community %>% filter( primer==pr)
  rownames(subset_community) <- subset_community$full_replicate
  
  zero_columns <- colSums(subset_community[, grep(asv_start, colnames(subset_community))]) == 0
  zero_column_names <- colnames(subset_community)[grep(asv_start, colnames(subset_community))][zero_columns]
  
  #remove from subset
  subset_community <- subset_community %>%select(-all_of(zero_column_names))
  
  #define indices new
  index <- grep(asv_start, colnames(subset_community))
  index_samples <-which(!grepl(grep_controls,subset_community$type))
  
  #caluclate
  bray_data <- metaMDS(subset_community[,index], distance="bray", parallel=threads)
  nmds <- plot_simpleNDMS_symbol(bray_data, subset_community, "type")
  
  #plot
  nmds[[2]] <- nmds[[2]] + 
    labs(title = paste(pr, "Similarity between all replicates ", sep=" - "),
         subtitle=" Bray Curtis Distance grouped after Sample Type (seed no. 25)",
         color="Sample Type")+ 
    scale_color_manual(values=color_df$color, breaks = color_df$type)
  print(nmds[[2]])
  
  ## remove low read replcicates
  print("low abundant replicates removed")
  sample_reads <- rowSums(subset_community[,index])
  sample_reads <- sample_reads > max(sample_reads)*0.01
  
  bray_data <- metaMDS(subset_community[sample_reads,index], distance ="bray", parallel=threads)
  nmds <- plot_simpleNDMS_symbol(bray_data, subset_community, "type")
  
  ##normal ordination plot
  nmds[[2]] <- nmds[[2]] + 
    labs(title = paste(pr, "Community Similarity between the Sample replicates (seed no. 25 )", sep=" - "),
         subtitle= "Low abundant replicates were excluded, Bray Curtis Distance",
         color="Sample Type")+ 
    plot_theme+
    scale_color_manual(values=color_df$color, breaks = color_df$type)
  
  nmds[[2]] <- ggdraw(add_sub(nmds[[2]], 
                 "Low abundant replicates: replicates with less than 1% reads of max sample read number\n Colors represent samples", size=10))
  print(nmds[[2]])
  

  ##### ordination with spiders for samples not controls
  #prepare samples for ghibli palette
  sample_no <- length(unique(subset_community[index_samples,"sample"]))
  color_samples <-as.vector(ghibli_palette(name = "PonyoMedium", n = sample_no, type = "continuous"))
  
  bray_data <- metaMDS(subset_community[index_samples,index], distance = "bray", parallel=threads)
  nmds <- plot_simpleNDMS_symbol(bray_data, subset_community[index_samples,], "sample")
  
  ##normal ordination plot
  nmds[[2]] <- nmds[[2]] + 
    labs(title = paste(pr, "Similarity between Sample replicates"),
         subtitle="Bray Curtis Distance grouped after Sample Type (seed no. 25), Controls were excluded")+ 
    plot_theme+
    theme(legend.position="none")+
    scale_color_manual(values=color_samples)
  nmds[[2]] <- ggdraw(add_sub(nmds[[2]], "Colors represent samples not sample types", size=10))
  print(nmds[[2]])
  
  ##spider plot
  nmds <- plot_simpleNDMS_spider(bray_data, subset_community[index_samples,], "sample")
  
  nmds[[2]] <- nmds[[2]] + 
    labs(title = paste(pr,"Community Similarity between the Sample replicates", sep=" - "),
         subtitle= "Controls were excluded, Bray Curtis Distance  (seed no. 25)")+ 
    plot_theme+
    theme(legend.position="none")+
    scale_color_manual(values=color_samples)
    nmds[[2]] <- ggdraw(add_sub(nmds[[2]], "Colors represent samples not sample types", size=10))
  print(nmds[[2]])
  
  ##############
  ## remove low read samples
  print("low abundant replicates removed")
  sample_reads <- rowSums(subset_community[index_samples,index])
  sample_reads <- sample_reads > max(sample_reads)*0.01
  
  bray_data <- metaMDS(subset_community[index_samples[sample_reads],index], distance ="bray", parallel=threads)
  nmds <- plot_simpleNDMS_symbol(bray_data, subset_community[index_samples,], "sample")
  
  ##normal ordination plot
  nmds[[2]] <- nmds[[2]] + 
    labs(title = paste(pr,"Community Similarity between the Sample replicates (seed no. 25 ) ", sep=" - "),
         subtitle= "Controls were excluded & low abundant replicates, Bray Curtis Distance")+ 
    plot_theme+
    theme(legend.position="none")+
    scale_color_manual(values=color_samples)
  
  nmds[[2]] <- ggdraw(add_sub(nmds[[2]], "Low abundant replicates: replicates with less than 1% reads of max sample read number\n Colors represent samples", size=10))
  print(nmds[[2]])
  
  ##spider plot
  nmds <- plot_simpleNDMS_spider(bray_data, subset_community[index_samples,], "sample")
  
  nmds[[2]] <- nmds[[2]] + 
    labs(title = "Community Similarity between the Sample replicates (seed no. 25) ",
         subtitle= "Controls were excluded & low abundant replicates, Bray Curtis Distance")+ 
    plot_theme+
    theme(legend.position="none")+
    scale_color_manual(values=color_samples)
  nmds[[2]] <- ggdraw(add_sub(nmds[[2]], "Low abundant replicates: replicates with less than 1% reads of max sample read number\n Colors represent samples", size=10))
  print(nmds[[2]])
  
  
  rm(bray_data, nmds,subset_community)
  
}

################################## MPrimer-multiplexing plots ##################################
## print proportion if different primers were multiplexed
if( length(unique(community$primer))>1){
  summary_primer <- community%>% group_by(primer)%>%
    dplyr::summarise( reads= sum(total.reads), asv= sum (total.asvs))
  
  # Basic piechart
  plot <- ggplot(summary_primer, aes(x="", y=reads, fill=primer)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    labs(title = "Proportion of reads of the different primers",
         x="",y="", fill="Primer")+ 
    plot_theme+ 
    scale_y_continuous(labels = scales::comma)+
    # ghibli stuff
    scale_fill_ghibli_d("LaputaMedium",direction = -1)
  
  print(plot)
  
  plot <- ggplot(summary_primer, aes(x="", y=asv, fill=primer)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    labs(title = "Proportion of ASVs of the different primers",
         x="",y="", fill="Primer")+ 
    plot_theme+ 
    scale_y_continuous(labels = scales::comma)+
    # ghibli stuff
    scale_fill_ghibli_d("LaputaMedium",direction = -1)
  
  print(plot)
}

dev.off()



