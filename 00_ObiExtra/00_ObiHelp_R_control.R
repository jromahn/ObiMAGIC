#!/usr/bin/env Rscript

print("Start checking R packages:")
library(cowplot)
library(dplyr)
library(gapminder)
library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(ghibli)
library(magrittr)
library(seqinr)
library(taxize)
library(taxonomizr)
library(treemapify)
library(vegan)

print("Causes issues with newer R releases:")
library(tidyverse)
library(plotly)
print("Every R package is installed :)")

## not needed anymore
#library(ROBIFastread)
#library(jsonify)
