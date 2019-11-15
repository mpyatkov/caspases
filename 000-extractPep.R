#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# in: table with cleavage sites
# out: fasta file for Blast
# Rscript --vanilla 000-extractPep.R table_with_cleavage_sites.csv

# Additional functions for working with Uniprot
source("00-2019-uni-fasta.R")

if (length(args)==0) {
  stop("Table with sites required (input file)", call.=FALSE)
}

SITES_TABLE <- args[1] ## "Table_S1_Human_caspase_targets.csv"

# print(args[1])
# stop("EXIT", call.=FALSE)

# Main file which will contain all FASTA sequences
DBFILE <- "./FASTADB.csv"

### Create fname by current day
# createFname <- function(new = F){
#   prefix <-  format(Sys.time(), "%d%m%y")
#   #prefix <- "t1"
#   suffix <- 1
#   if (new) {
#     suffix <- ifelse(file.exists(paste0(prefix,"_",suffix)),
#                      as.integer(
#                        str_split(
#                          rev(dir(pattern = paste0(prefix,"_*")))[1], 
#                          pattern = "_", simplify = T)[2])+1,
#                      1)
#   }
#   paste0("run",prefix,"_",suffix)
# }
# createFname(new = F)

########### pipeline create db, or append to db
### all uni id
all_uni <- read_csv(SITES_TABLE) %>% 
  select(uni) %>% 
  distinct(uni)

### DB of downloaded FASTA sequences
# load file if created or create empty table
DB <- openDB(DBFILE) 

# test
# DB <- appendListToDB(tibble(uni = c("Q76657","Q76658")))

### fill DB by uni id
system.time(
  DB <- appendListToDB(all_uni, DB)
)

### write DB
write_csv(DB, path = DBFILE)    

######### create ....[-29]D[+30]... from full protein
uni_pep <- read_csv(SITES_TABLE) %>% 
  select(uni, fullpep) %>% 
  left_join(., y = DB) %>% 
  rowwise() %>% 
  mutate(header = paste0("> ",uni,"_",fullpep), 
         chunk = getNamino(fasta, fullpep, 30), 
         rec=paste0(header,"\n",chunk)) %>% 
  select(rec) %>% 
  pull() %>% 
  # write_lines(createFname(new = F))
  write_lines("SITES_60AA.txt")


