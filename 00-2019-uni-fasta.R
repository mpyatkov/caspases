library(httr)
library(dplyr)
library(stringr)
library(parallel)
library(purrr)
require(readr)

### Open DB file
openDB <- function(dbfile){
  if(file.exists(dbfile)) {
    DB <- read_csv(dbfile)
  }
  else {
    DB <- tibble(uni = character(), fasta = character())
  }
  DB
}

# Download the FASTA file from Uniprot and parse them, extracting seq, pname, gname
# input: Uniprot id
# output: c(uni, pname, gname, seq, unparsed_content)
getFasta <- function(uni) {
  r <- GET(paste0("https://www.uniprot.org/uniprot/",uni,".fasta"))
  p <- r %>% content(.,"text", encoding = "UTF-8") %>% str_split(., pattern = "\\n", simplify = TRUE)
  if (str_detect(p[1],"DOCTYPE")){
    warning(paste0(uni, " not found"))
    return(NA)
  }
  
  pname <- str_extract_all(p[1],"(?<=\\|)([[:alnum:]]+)", simplify = T)[2]
  gname <-  str_extract_all(p[1],"(?<=GN=)([[:alnum:]]+)", simplify = T)[1]
  seq <- str_c(p[2:length(p)], collapse = "")
  c(uni,pname,gname,seq, p[1])
}

### Append list of uni to db with sequences downloaded from uniprot
# input: data.frame with Uniprot id (colimn: uni)
# input: db file 
# output new db file
appendListToDB <- function(unitab, db) {
  
  uni_tibble <- anti_join(unitab, db)
  
  if (nrow(uni_tibble) < 1) {
    return(DB)
  }
  
  fasta_list <- mclapply(c(uni_tibble$uni), getFasta, mc.cores = 2)
  uni_seq <- tibble(uni = uni_tibble$uni, 
                    pname =map_chr(fasta_list, ~ .x[2]),
                    gname =map_chr(fasta_list, ~ .x[3]),
                    fasta = map_chr(fasta_list, ~ .x[4]),
                    header = map_chr(fasta_list, ~ .x[5]))
  
  tmp <- bind_rows(db, uni_seq)
  tmp
}

### get +/-30 amino acids from center Asp
# input: protein string
# input: 8AA peptide (site)
# input: size of left or right hand from central AA
getNamino <- function(st, pep, limit) {
  start_end <- str_locate(st,pep)
  start <- ifelse((start_end[1]-limit+4)<1,1,start_end[1]-limit+4)
  end <- start_end[2]+limit-4
  str_sub(st, start, end)
}

## tests
# getFasta("P10147")
# 
# appendListToDB(tibble(uni= c("P10147")), 
#                tibble(uni = character(), pname=character(), gname=character(),fasta=character()))
