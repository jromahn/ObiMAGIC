#!/usr/bin/env bash

for i in $(find $(pwd) -type d -name "*_results" | sort); do cd $i/.. && tar -czvf $(echo $i | awk -F'/' '{print $NF}').final-results.tar.gz $(find $(basename $i) -type f -name "*RData" -or -name "*pdf" -or -name "*tsv" -or -name "*csv" -or -name "*txt" | sort) && tar -czvf $(echo $i | awk -F'/' '{print $NF}').raw-results.tar.gz $(find $(basename $i) -type f -name "*$(echo $i | awk -F'/' '{gsub("_results$","",$NF);print $NF}')*fasta" -or -name "*$(echo $i | awk -F'/' '{gsub("_results$","",$NF);print $NF}')*csv" | sort); done
