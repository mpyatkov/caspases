# get 60AA qseq multiple alignment by group (uniprot, fullpep)

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")                   
BiocManager::install("DECIPHER")

source("./00-file_path.R") ## for pp - path function
library(DECIPHER)
library(tidyverse)

repr <- read_rds("./01-5-repr.rds")

repr_test <- repr %>%
  filter(grepl("CALU", pname),fullpep == "YSHDGNTD") %>%
  mutate(hseq=str_replace_all(hseq, '-',""))

## same for map_df functions
align_dataframe_map <- function(uni, fullpep, vorgs, vhseq) {
  names(vhseq) <- vorgs
  algn <- AlignSeqs(AAStringSet(vhseq))
  va <- as.vector(algn)
  names(va) <- c()
  tibble(
    uni = uni,
    fullpep = fullpep,
    orgs=vorgs,
    algn = va)
}

## Table with alignments
alignments_table <- repr_test %>% 
  select(uni, fullpep, org, hseq) %>% 
  group_by(uni, fullpep) %>% 
  summarise(org=list(org), hseq=list(hseq), cnt=n()) %>% 
  filter(cnt > 1) %>% 
  select(-cnt) %>% 
  ungroup() %>%
  rowwise() %>% 
  select(uni, fullpep, vorgs = org, vhseq = hseq) %>% 
  pmap_df(align_dataframe_map) %>% 
  write_csv(pp("021-TAB-align.csv.tar.bz2"))
