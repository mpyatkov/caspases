# For processing we used mufold server:
# http://dslsrv2.eecs.missouri.edu/~zlht3/ss

# Comment: Long procedure ~ 2 days, because MUFOLD server allow only one 60AA
# sequence at the time (~1m:30s processing time). 
# Prepare sequences as FASTA strings, we can send no more then 10 sequence on request, 
# with length no more 1000

# TODO: future update: pack sequences by 16 (1000/60~16) in one string
# ex. site1AASEPsite2AASEPsite3AASEP...AASEPsite16 and pass
# like one protein, after that parse the result. It probably allows to process faster

library(tidyverse)
library(httr)
library(rvest)

source("00-file_path.R") ## for pp
source("00-2019-uni-fasta.R") ## for getNAmino

repr <- read_rds("./01-5-repr.rds")

# read table with sites
s1 <- read_csv("Table_S1_Human_caspase_targets.csv") %>% 
  select(uni,fullpep)

# uni, fasta
db <- read_csv("FASTADB.csv") %>% 
  select(uni,fasta)

# Prepare only human sequences for processing in remote server
seq_to_server <- s1 %>% 
  mutate(fullpep = str_replace_all(fullpep, "-","")) %>% 
  inner_join(db, y = .) %>% 
  rowwise() %>% 
  mutate(header = paste0("> ",uni,"_",fullpep), 
         chunk = getNamino(fasta, fullpep, 30), 
         rec=paste0(header,"\n",chunk)) %>% 
  ungroup() %>% 
  mutate(num = row_number()%%round(1+n()/10)) %>% 
  select(rec,num) %>% 
  group_by(num) %>% 
  summarise(fasta_line=str_c(c(rec), sep = "\n", collapse = "\n"))
  
# Email is just an identifier. At the time of November 13, 2019, 
# email is not required to be real. Just put any.
email = "testtest@test.com"

# send sequence to server
send_to_server <- function(seq, email){
  lbody = list(email=email, target_name = "", 'tasks[]'=c(1), sequence=seq)
  POST("http://dslsrv2.eecs.missouri.edu/~zlht3/ss/submit", 
       body = lbody,
       encode = c("form"))
  
  # This delay not required, but we need be polite
  Sys.sleep(3)
}

# send sequences to server
map2(seq_to_server$fasta_line, email, 
     function(fasta_line, email) send_to_server(fasta_line, email))

## ~900 sec (because delay in 3 sec)
send_to_server(testseq,email)

# Get page with results back using email identifier
retrive_results_status <- function(email){
  r <- POST("http://dslsrv2.eecs.missouri.edu/~zlht3/ss/retrieve",
            body = list(jobid="",email=email),
            encode = c("form"), verbose())
  r
}

## get results page with the table of finished sequence
results <- retrive_results_status(email)

## extract table from results using xpath
jobid_entries <- results %>% content() %>% 
  html_nodes(xpath='//html/body/div[2]/div[1]/div[2]/div/div/div[2]/table') %>%
  html_table()

# add right header to the table 
jobid_entries <- jobid_entries[[1]]
names(jobid_entries)<-c("num", "status", "target","time","len","sstype", "jobid")
write_csv(jobid_entries, pp("051-TMP-jobid_entries.csv"))

### get Q3 and Q8 from remote server
# http://dslsrv2.eecs.missouri.edu/~zlht3/ss/download_ss_results_only/ss_5d8936a4a0d86
get_prediction <- function(jobid){
  url <- "http://dslsrv2.eecs.missouri.edu/~zlht3/ss/download_ss_results_only/"
  ss <- GET(paste0(url, jobid))
  ss %>% 
    content(as="text", encoding = "UTF-8") %>% 
    str_trim(.) %>% 
    str_split(., "\n", simplify = T) %>% 
    list(Q3 = .[,2], Q8=.[,4])
}

# if nothing for download function returns null and terminates which 
# we use safely to analyse return type
safe_get_prediction <- possibly(get_prediction,list(Q3=NA,Q8=NA))

library(tictoc)
tic()
finished_prots <- jobid_entries %>% 
  select(target, jobid) %>% 
  pmap_df(., function(target, jobid) {
    Sys.sleep(0.3)
    print(target)
    lst <- safe_get_prediction(jobid)
    tibble(
      target = target,
      Q3 = lst$Q3,
      Q8 = lst$Q8
    )})
toc()

write_csv(finished_prots, pp("051-TMP-finished-prots.csv"))
finished_prots <- read_csv(pp("051-TMP-finished-prots.csv"))

## in our case one site wasn't calculated due to internal error, I calculated it in other service (YASPIN)
## but anyway we can use the main service to recalculate it if we desire
finished_prots <- finished_prots %>% 
  mutate(Q3 = ifelse(target == "P85037_VSGDSAVA_5",
                     "CCCCCCCCCCCCCCEECCCCCCCCCCCCCCHHHHCHHHHHHHHHHHHHHCCCCCHHHHCC",
                     Q3),
         Q8 = ifelse(target == "P85037_VSGDSAVA_5",Q3,Q8))

# final table which contain uni,fullpep, Q3,Q8
final_domains <- finished_prots %>% 
  separate(target,c("uni","fullpep","num")) %>%
  select(-num) %>% 
  write_csv("051-TAB-final-domains.csv") # must 051-TAB-full_domains.csv

# additional info in final domains
repr_domains <- read_csv(pp("051-TAB-final-domains.csv")) %>% 
  inner_join(repr %>% 
               select(uni, fullpep, qseq,org) %>% 
               filter(grepl("Homo", org)) %>% 
               group_by(uni,fullpep) %>% 
               distinct()) %>%
  select(-org) %>% 
  write_csv(pp("051-TAB-repr_domains.csv"))
