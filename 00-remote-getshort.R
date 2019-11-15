#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(readr))
setwd("./")

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} 

# print(args[1])
# stop("EXIT", call.=FALSE)

# fullpath <- "full.tsv.gz" 
file_name <- unlist(strsplit(args[1],"\\."))[1]
full_table <- read_tsv(args[1], col_names = FALSE)

## add names to columns
names(full_table) <- c("uni", "accesion", "title", "evalue", "score", "identity", "qseq", "hseq")

## separate uni -> uni,fullpep
system.time(
  SHORT <- full_table %>% 
    mutate(org = str_extract(title,"(?<=\\[)(([[:upper:]]){1}([[:lower:]]|[[:blank:]])+)(?=]$)")) %>%
    mutate(org = str_extract(org,"([[:alpha:]]+ [[:alpha:]]+)")) %>%
    filter(!grepl("^\\d|^\\(|,|-|\\&| x |\\.", org) & grepl("^[[:upper:]]",org), org != "Unknown", !is.na(org)) %>%
    separate(uni, sep = "_", into = c("uni", "fullpep")) %>%              ## for uni_fullpep
    # separate(uni, sep = "_", into = c("uni", "fullpep", "pname", "gi")) %>% ## for uni_fullpep_pname_gi
    group_by(uni, org, fullpep) %>% slice(which.min(evalue)) %>% ungroup()
)

write_csv(SHORT, path = paste(file_name,"_SHORT.csv.gz", sep = ""), append = FALSE)

## addlinages
# lin <- read_tsv("resp.tsv", col_names = FALSE)
# names(lin) <- c("org","taxid","linage")
# 
# sl <- SHORT %>%
#   left_join(y = lin) %>%
#   filter(!is.na(linage))

# extract unique orgs for taxonkit
print("Create uniq orgs table from short")
uniq_orgs <- SHORT %>%
  select(org) %>%
  distinct(org) %>%
  write_csv(. ,path = paste(file_name,"_TABLE_UNIQ_ORGS.csv", sep=""), append = FALSE, col_names = FALSE)

