
#############################
### pattern to extract from replicate name: sample (level above repl), type & ASV names
# matching pattern will be removed via gsub(so describe how to not find if)
repl_identifier="_(R|P).*" ## pattern to remove replicate & plate names
# new one if system changed to KAR23_0001_R01 --> first part includes information about project and sample type
type_identifier="_.*" ## patter to identify if it is a sample, a negative control or positive control
asv_start="ASV" ### prefix for rename and enumerating asv names
grep_negative_controls <- "FN$|SN$|EN$|MN$|PN$" # patter to find negative controls
grep_controls <- "FN$|SN$|SP$|EN$|FP$|EP$|MN$|PN$|PP$" # patter to find negative & positive controls
#############################

## categorize data by negative and positive control shortnings
categorize_community <- function(data,asv_start){
  #first calculcate total number of ASV present
  index <- grep(asv_start, colnames(community))
  asv <- data %>% select(type,all_of(index) )%>%
          pivot_longer(!type, names_to="ASV", values_to = "reads")%>%
          mutate(occurence=1*(reads>0))%>%
          select(-reads)%>%
          unique()%>%
          group_by(type)%>%
          summarise(total_asv=sum(occurence))%>%
          ungroup()

  data <- data %>% 
    group_by(type)%>%
    summarize(replicates_total=n(),
              replicates=length(unique(replicate)),
              samples_total=length(unique(full_sample)),
              samples= length(unique(sample)),
              total_reads=sum(total.reads))%>%
    ungroup()%>%
    left_join(asv, by ="type")%>%
    mutate(type_long=ifelse(grepl("FN$", type), "Field Neg. Control",
                     ifelse(grepl("SN$", type), "Sample Neg. Control",
                     ifelse(grepl("EN$", type), "Extraction Neg. Control",
                     ifelse(grepl("PN$", type), "PCR Neg. Control",
                     ifelse(grepl("MN$", type), "Multiplexing Control", 
                     ifelse(grepl("FP$", type), "Field Pos. Control",
                     ifelse(grepl("SP$", type), "Sample Pos. Control",
                     ifelse(grepl("EP$", type), "Extraction Pos. Control",
                     ifelse(grepl("PP$", type), "PCR Pos. Control","Sample"))))))))))%>%
  
  
  return(data)
}


###### extract replicate and sample information from communtiy rownames
extra_sample_names <- function(data, identifier){
  
  ## primer info
  data$primer <- gsub("__.*", "",  rownames(data))
  # if no primer name replace with NA
  if(data$primer[1] %in% rownames(data)[1]){
    data$primer= "NA"
  }
  
  ##sample info
  data <- data %>%
    rownames_to_column(var = "full_replicate")%>%
    mutate(replicate      = gsub(".*__", "",  full_replicate) )%>% # extract replicate number
    mutate(full_sample    = gsub(repl_identifier,"",full_replicate))%>% # extract sample name
    mutate(sample         = gsub(repl_identifier,"",replicate))%>% # extract sample name
    mutate(repl_type      = gsub(".*_([a-zA-Z0-9]+)$","\\1",replicate, perl = T)) %>%# identify replicates & plates
    mutate(type           = gsub(identifier,"",sample))%>% # extract sample type 
    mutate(type           = gsub("_","",type))
  
  return(data)
}

## extract plate number, rows and columns if those information exist

extract_plate_information <- function(plate_info, index_position){
    plate_info <- plate_info[,c(2,index_position)] # just get the neccessary columns
  colnames(plate_info) <- c("replicate", "plate_position")
  plate_info <- plate_info %>% 
    mutate(plate_position= gsub(".*position=(\\w+)", "\\1", plate_position, perl=T) ) %>%
    separate(plate_position, c("plate_no", "plate_position"))%>%
    mutate(plate_row    = gsub("\\d+","", plate_position),
          plate_column = gsub("(\\d+).*","\\1", plate_position, perl=T))
    
  return(plate_info)
}

################################## GGPLOT color & theme variables #############################
## define  plot theme
plot_theme <- theme_light() + 
  theme(plot.title= element_text(face="bold"),
        axis.text = element_text(size = 10),
        legend.position="bottom",
        strip.background =element_rect(fill="grey38"), #grey38
        panel.grid = element_line(color = "grey78", linewidth = 0.08))

## define colors
colors <- c("#104E8B", "#1F83B4FF", "#12A2A8FF", "#78A641FF", "#BCBD22FF", "#FFAA0EFF", "#FF7F0EFF", "#BA43B4FF",  "#6F63BBFF","#746455")


####################################       FUNCTIONS      #####################################
## taxonomy provided by ObiTools4
obirank_level= c("no rank","superkingdom", "kingdom", "clade","class","subclass","superorder","order" ,
                 "family","subfamily","tribe","subtribe", "genus" ,"species",  "subspecies")

## get ncbi taxonomy
ncbi_taxonomy <- function(unique_assigned){
  unique_assigned <- unique(unique_assigned)
  rank_information <- getTaxonomy(unique_assigned, sqlFile, 
                                  desiredTaxa = c("species","genus","family","order","class", "subphylum" ,"phylum", "kingdom","superkingdom"))
  rank_information <- data.frame(rank_information)
  rank_information$taxID <- unique_assigned
  return(rank_information)
}

####################################   ORDINATION   FUNCTIONS      #####################################
# ordination plot of communtiy data based on replicate data --> rownames
# column defines after which column the color is chosen
plot_simpleNDMS_symbol <- function(NMDS_result, data, column_c){
  replicate_coordinates <- data.frame(NMDS_result["points"])
  replicate_coordinates$full_replicate <- rownames(replicate_coordinates)
  replicate_coordinates <- merge(replicate_coordinates, data, by ="full_replicate", all.x = TRUE)
  
  #criteria color
  colnames(replicate_coordinates)[colnames(replicate_coordinates) == column_c] <- 'crit_color'
  replicate_coordinates$crit_color <- as.factor(replicate_coordinates$crit_color )
  
  # Now we can build the plot! 
  p <- ggplot() +
    geom_point(data = replicate_coordinates, aes(x = points.MDS1, y = points.MDS2, 
                                                 color = crit_color), size=2) + 
    annotate(geom = "label", 
             x =(sqrt(min(replicate_coordinates$points.MDS1)^2)*1.1), 
             y =min(replicate_coordinates$points.MDS2), 
             size =4,
             label = paste("Stress: ", round(NMDS_result$stress, digits = 3)))+ #+
    labs(color = column_c) + 
    plot_theme
  list_data <- list(replicate_coordinates, p)
  return(list_data)
}

#plot_simpleNDMS_spider(bray_data, community[index_samples,], "sample")
plot_simpleNDMS_spider <- function(NMDS_result, data, column_c){
  replicate_coordinates <- data.frame(NMDS_result["points"])
  replicate_coordinates$full_replicate <- rownames(replicate_coordinates)
  replicate_coordinates <- merge(replicate_coordinates, data, by ="full_replicate", all.x = TRUE)
  
  #criteria color
  colnames(replicate_coordinates)[colnames(replicate_coordinates) == column_c] <- 'crit_color'
  replicate_coordinates$crit_color <- as.factor(replicate_coordinates$crit_color )
  
  #data for spider
  cent <- aggregate(cbind(points.MDS1, points.MDS2) ~ crit_color, data = replicate_coordinates, FUN = mean)
  segs <- merge(replicate_coordinates, setNames(cent, c('crit_color','oNMDS1','oNMDS2')),
                by = 'crit_color', sort = FALSE)
  
  # Now we can build the plot! 
  
  p <- ggplot(replicate_coordinates, aes(x = points.MDS1, y = points.MDS2, 
                                         color = crit_color)) +
    geom_segment(data= segs, mapping = aes(xend = oNMDS1, yend = oNMDS2))+
    geom_point( size=2) + 
    annotate(geom = "label", 
             x =(sqrt(min(replicate_coordinates$points.MDS1)^2)*1.1),
             y =min(replicate_coordinates$points.MDS2), 
             size =4,
             label = paste("Stress: ", round(NMDS_result$stress, digits = 3)))+ #+
    labs(color = column_c) + 
    plot_theme+ 
    theme(legend.position = "None")
  list_data <- list(replicate_coordinates, p)
  return(list_data)
}
####################################################################################################
