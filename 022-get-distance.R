library(tidyverse)
source("00-file_path.R")

# 1. For each pair of organims was obtained list (LP) of shared proteins 
# 2. For each protein (P) of this subset we compared coresponded sites 
#    using hamming distance.
# 3. All hamming distances related with one P from LP were averaged and divided on
#    maximum distance, for short sequences (sites) 8AA in length - 8, for 60AA (hseq) - 60
# 3.1 For 60AA sequences was provided multiple alignment for each group (Uniprot, Human site)
#     after that alignments were compared to each other using hamming distance.
# 4. As result distance matrix among all organisms
# 5. On this step used only highly represented organisms (~329)

repr <- read_rds("./01-5-repr.rds")

## numerate all (pname, fullpep. count) from 1 to ...
proteind_id <- repr %>% ## 
  select(gname, fullpep) %>%
  group_by(gname, fullpep) %>% 
  summarize(count = n()) %>% 
  ungroup() %>% 
  mutate(prot_id = row_number())

## org, vector_protid = c(pid), orgid
org_vpid_orgid <- repr %>% 
  inner_join(y = proteind_id, by = c("gname","fullpep")) %>%
  select(org, prot_id) %>%
  group_by(org) %>%
  summarise(vprot_id = list(prot_id)) %>%
  mutate(orgid = row_number()) %>% 
  ungroup()

## all combinations of orgid n(n-1)/2 (only lower diagonal)
expanded_orgs <- org_vpid_orgid %>% 
  select(orgid) %>%
  mutate(orgid1 = orgid) %>%
  expand(orgid,orgid1) %>%
  filter(orgid<orgid1)

## vectorized function of intersection between vectors of pid
prot_intersect <- function(n1,n2) {
  intersect(org_vpid_orgid[n1,]$vprot_id[[1]], org_vpid_orgid[n2,]$vprot_id[[1]])
}
prot_intersect <- Vectorize(prot_intersect) 

## table orgid1,orgid2, intersect_vector (9 seconds)
system.time(
  ex_orgs_dist <- expanded_orgs %>%
    mutate(orgdist = prot_intersect(orgid,orgid1))
)

# load hamming distance function
source("00-octet.R")   

## (org_index1, org_index2) -> (1,2,3,4) -> (pro1,prot2,prot3) -> compare octets of such proteins
counter <- 0
list_to_pid_tibble <- function(oi1,oi2) {
  if (oi1 != counter) {
    counter <<- oi1
    print(counter)
  }
  
  ## orgid to org
  line <- filter(ex_orgs_dist, orgid == oi1 & orgid1 == oi2)
  
  o1 <- org_vpid_orgid[oi1,]$org
  o2 <- org_vpid_orgid[oi2,]$org
  
  tmp <- tibble(prot_id = line$orgdist[[1]], fakecol=1L) %>%
    arrange(prot_id) %>%
    inner_join(proteind_id, ., by = "prot_id") %>% 
    select(gname, fullpep, fakecol)
}

# get distance between two organisms
get_distance <- function(oi1, oi2, col) { # col == algn or octet
  
  seq_len <- 8
  if (col == "algn") {seq_len <- 60 }

  tmp <- list_to_pid_tibble(oi1, oi2)
  
  oo1 <- org_vpid_orgid[oi1,]$org
  oo2 <- org_vpid_orgid[oi2,]$org
  
  repr %>% 
    filter(org == oo1 | org == oo2) %>%
    inner_join(x = ., y = tmp, by=c("gname", "fullpep")) %>%
    #select(gname, fullpep, octet, org) %>%
    select(gname, fullpep, octet=col, org) %>% 
    distinct() %>%
    spread(org, octet) %>%
    rename_at(3,~"o1") %>%
    rename_at(4,~"o2") %>%
    rowwise() %>%
    mutate(hamd = hammingStr(o1,o2)) %>%
    ungroup() %>%
    summarize(mn = mean(hamd,na.rm = T)/seq_len) %>% 
    pull(mn)
}

## long time operation , comparison all organisms to each other 8AA
system.time(
  qq <- ex_orgs_dist %>% 
    #sample_n(100) %>% 
    rowwise() %>%
    mutate(dst = get_distance(orgid,orgid1,"octet")) %>%
    arrange(dst) %>%
    mutate(org1 = org_vpid_orgid[orgid,]$org,
           org2 = org_vpid_orgid[orgid1,]$org) %>% 
    select(org1,org2, dst) %>% 
    write_csv(., pp("022-TAB-distance_bw_orgs_8AA.csv"))
)  
  
# Append alignment information about 60AA to repr table
algn <- read_csv(pp("021-TAB-align.csv.tar.bz2"))
repr <- repr %>% 
  inner_join(., algn, by = c("uni","fullpep","org"))

## long time operation, comparison all organisms to each other 60AA
system.time(
  qq <- ex_orgs_dist %>% 
    #sample_n(100) %>% 
    rowwise() %>%
    mutate(dst = get_distance(orgid,orgid1,"algn")) %>% 
    arrange(dst) %>%
    mutate(org1 = org_vpid_orgid[orgid,]$org,
           org2 = org_vpid_orgid[orgid1,]$org) %>% 
    select(org1,org2, dst) %>% 
    write_csv(., pp("022-TAB-distance_bw_orgs_60AA.csv"))
)  
  
  