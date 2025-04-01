library(tidyverse)
library(gapminder) # allow commas in numbers --> +scale_y_continuous(labels = scales::comma)

options(scipen = 999) # stop scientific notation

rm(list=ls())

#setwd("/Users/juliane/Documents/00_Work_SGN/00_little_scripting/metabarcoding/2025_obimagic")
#inputfile<- "Johannes_MP4_dec24_results/3_unidentified__reasons.tsv"
#project <- "Johannes_MP4_dec24"

#### get input
args = commandArgs(trailingOnly=TRUE)
inputfile=args[1] #"Malawi_2023_eDNA_results/"
path=dirname(inputfile)
project=gsub("_results", "", path)
number=gsub("(\\d+)_.*", "\\1", basename(inputfile), perl =T)

data <- read.table(inputfile, header = T, sep="\t",comment.char = "")
colnames(data) <- gsub("X\\.", "", colnames(data))



pdf(file = file.path(path, paste(number, "unidentified__diagnosis_plots.pdf", sep="_")), paper = "a4r")

summary <- data %>% group_by(Error_message)%>%
  dplyr::summarise(count= n())



## summary of error messages
summary%>%
  mutate(name = fct_reorder(Error_message, desc(count))) %>%
  ggplot( aes(x=name, y=count)) + 
      geom_bar(stat = "identity",fill="#69b3a2", color="#e9ecef", alpha=0.8)+
      coord_flip()+
  scale_y_continuous(labels = scales::comma)+
  theme_light()+ theme(plot.title = element_text(face = "bold"))+
  labs(subtitle= "The abundance of the demulitplex error message", title=project, x="Error message", y="Abundance/Count")



## length distribution
ggplot(data, aes(x=Seq_length)) + 
  geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
  theme_light()+ theme(plot.title = element_text(face = "bold"))+
  labs(title= "Sequence Length distribution", x="Sequence Length (bp)", y="Abundance")


## length distribution
data %>% filter(Proportion_Ns > 0)%>%
    mutate(Proportion_Ns = Proportion_Ns *100)%>%
    ggplot( aes(x=Proportion_Ns)) + 
        geom_histogram( fill="#69b3a2", color="#e9ecef", alpha=0.8) +
        theme_light()+ theme(plot.title = element_text(face = "bold"))+
        coord_cartesian(xlim =c(0, max(data$Proportion_Ns)*1000))+ # makes 0 and 10% more of % as xlim
        labs(title= "Proportion of Ambiguous nucleotide distribution\n with tag & primer pair", x="Proportion (%)", y="Abundance")


# most common multiplexing problem

failure <- data %>% filter( grepl("assign|associate",Error_message ) & grepl("sample",Error_message ))%>%
  mutate(exact_problem = case_when( (Forward_Tag != "0" & Reverse_Tag != "0" ) ~ "Both tags missing",
                                    Forward_Tag != "0"                         ~ "Forward tag missing",
                                    Reverse_Tag != "0"                         ~ "Reverse tag missing",
                                    TRUE                                       ~ "Tag combination unknown"))%>%
  group_by(exact_problem) %>% summarise(count=n()) %>% ungroup()%>%
  mutate(exact_problem = fct_reorder(exact_problem, desc(count)))


ggplot(failure, aes(x=exact_problem, y=count)) + 
  geom_bar(stat = "identity",fill="#69b3a2", color="#e9ecef", alpha=0.8)+ 
  scale_y_continuous(labels = scales::comma)+
  theme_light()+ theme(plot.title = element_text(face = "bold"))+
  coord_flip()+
  labs(title= "Reason of unassociated tag combination", x="Reason", y="Abundance/Count")


# common unknown tag combination
failure <- data %>% filter( grepl("assign|associate",Error_message ) & grepl("sample",Error_message ))%>%
  filter(Forward_Tag != "0" & Reverse_Tag != "0")%>%
  mutate(Tag_combination = paste(Forward_Tag,":",Reverse_Tag, sep="")) %>%
  group_by(Tag_combination)%>%
  dplyr::summarise(count= n()) %>%
  ungroup()%>%
  slice_max(count, n = 10)%>%
  mutate(Tag_combination = fct_reorder(Tag_combination, desc(count)))


ggplot(failure, aes(x=Tag_combination, y=count)) + 
  geom_bar(stat = "identity",fill="#69b3a2", color="#e9ecef", alpha=0.8)+ 
  scale_y_continuous(labels = scales::comma)+
  theme_light()+ theme(plot.title = element_text(face = "bold"))+
  coord_flip()+
  labs(title= "Most common tag combination with multiplexing problems", x="PCR multiplexing tag", y="Abundance/Count")

## pattern of tag length
summary <- data %>% filter( grepl("assign|associate",Error_message ) & grepl("sample",Error_message ))%>%
    group_by(Forward_Tag_Length,Reverse_Tag_Length)%>%
    dplyr::summarise(count= n()) 

## define parameters for size legend
minimum <- round(min(summary$count), digits = -3) 
maximum <- round(max(summary$count), digits = -3) 
breaks_size <- seq(minimum, maximum,  length.out = 5)
ggplot(summary, aes(x=Forward_Tag_Length, y=Reverse_Tag_Length, size=count)) + 
          geom_point( )+
          theme_light()+ theme(plot.title = element_text(face = "bold"))+
          scale_size_continuous(limits= c(minimum, maximum), breaks = breaks_size, labels = scales::comma)+
          labs(title= "Tag length of unassignable sequences", x="Forward Tag Length (bp)", y="Reverse Tag Length (bp)", size= "Abundance/Count")


#### common forward tags
data %>% 
  filter(Forward_Tag !=  0)%>%
  filter(Forward_Tag_Length >=8)%>%
  group_by(Forward_Tag)%>%
  dplyr::summarise(count= n()) %>%
  arrange(desc(count))%>%
  dplyr::slice(1:20) %>%
  mutate(name = fct_reorder(Forward_Tag, desc(count))) %>%
  ggplot( aes(x=name, y=count)) + 
      geom_bar(stat = "identity", fill="#69b3a2", color="#e9ecef", alpha=0.8)+
      coord_flip()+
      scale_y_continuous(labels = scales::comma)+
      theme_light()+ theme(plot.title = element_text(face = "bold"))+
      labs(title= "Most common forward tags", x="Forward Tag", y="Abundance")


### stats for forward primer
summary <- data %>% filter(!is.na(Forward_Primer_Mimatches))%>%
  group_by(Forward_Primer_Mimatches, Forward_Tag_Length)%>%
  dplyr::summarise(count=n())

## define parameters for size legend
minimum <- round(min(summary$count), digits = -3) 
maximum <- round(max(summary$count), digits = -1) 
breaks_size <- seq(minimum, maximum,  length.out = 5)

summary%>%
  ggplot(aes(x=Forward_Tag_Length, y= Forward_Primer_Mimatches, size= count))+ 
  geom_point( )+
  theme_light()+ theme(plot.title = element_text(face = "bold"))+
  scale_size_continuous(limits= c(minimum, maximum), breaks = breaks_size, labels = scales::comma)+
  labs(title= "Forward primer & tag statistics", x="Forward Tag Length (bp)", y="Forward Primer Mismatches", size= "Abundance/Count")



#### common reverse tags
data %>% 
  filter(Reverse_Tag !=  0)%>%
  filter(Reverse_Tag_Length >=8)%>%
  group_by(Reverse_Tag)%>%
  dplyr::summarise(count= n()) %>%
  arrange(desc(count))%>%
  dplyr::slice(1:20) %>%
  mutate(name = fct_reorder(Reverse_Tag, desc(count))) %>%
  ggplot( aes(x=name, y=count)) + 
      geom_bar(stat = "identity",fill="#69b3a2", color="#e9ecef", alpha=0.8)+
      coord_flip()+
      scale_y_continuous(labels = scales::comma)+
      theme_light()+ theme(plot.title = element_text(face = "bold"))+
      labs(title= "Most common reverse tags", x="Reverse Tag", y="Abundance")




### stats for reverse primer
summary <- data %>% filter(!is.na(Reverse_Primer_Mimatches))%>%
  group_by(Reverse_Primer_Mimatches, Reverse_Tag_Length)%>%
  dplyr::summarise(count=n())

## define parameters for size legend
minimum <- round(min(summary$count), digits = -3) 
maximum <- round(max(summary$count), digits = -1) 
breaks_size <- seq(minimum, maximum,  length.out = 5)

summary%>%
  ggplot(aes(x=Reverse_Tag_Length, y= Reverse_Primer_Mimatches, size= count))+ 
    geom_point()+
    theme_light()+ theme(plot.title = element_text(face = "bold"))+
    scale_size_continuous(limits= c(minimum, maximum), breaks = breaks_size, labels = scales::comma)+
    labs(title= "Reverse primer & tag statistics", x="Reverse Tag Length (bp)", y="Reverse Primer Mismatches", size= "Abundance/Count")
  

dev.off()