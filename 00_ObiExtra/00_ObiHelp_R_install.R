#!/usr/bin/env Rscript

install.packages("jsonify", repos='http://cran.us.r-project.org')
install.packages("devtools", repos='http://cran.us.r-project.org')
devtools::install_git("https://git.metabarcoding.org/obitools/obitools4/robireadfasta.git")

install.packages(c("cowplot","dplyr","gapminder","ggplot2","ggpmisc","ggpubr","ghibli","magrittr","plotly","seqinr","taxize","taxonomizr","tidyverse","treemapify","vegan"),repos='http://cran.us.r-project.org')

install.packages('taxize', repos = c('https://ropensci.r-universe.dev', 'https://cloud.r-project.org'))
