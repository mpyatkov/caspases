library(tidyverse)
source("00-file_path.R")

repr <- read_rds("./01-5-repr.rds")

## Add information about how AA changed to each other at P1' position
## all types of changes (stab-> destab, destab->destab,...)
nrule_changes <- repr %>% 
  filter(found == "D") %>% 
  # filter(fullpep_ntype == octet_ntype, fullpep_ntype == 0, fullpep_nrule_amino != octet_nrule_amino) %>% 
  filter(fullpep_nrule_amino != octet_nrule_amino) %>% 
  rowwise() %>% 
  mutate(amino_status = paste0(fullpep_nrule_amino, " -> ",octet_nrule_amino, collapse = ""),
         status = paste0(fullpep_stype, " -> ",octet_stype, collapse = "")) %>% 
  select(uniprot = uni, gname,pname,org, Class, human_site = fullpep, org_site=octet, evalue, aa_from_to = amino_status, blast_description = title, human_site_long = qseq, org_site_long = hseq, status) %>% 
  arrange(uniprot, human_site, org) %>% 
  select(uniprot, gname, pname, org, Class, human_site, org_site, aa_from_to, status, evalue, blast_description) 

# in which sites more changes in P1'?
stats_by_nrule_changes <- nrule_changes %>% 
  group_by(uniprot, human_site, gname, pname) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count)) %>% 
  write_csv(pp("04-TAB-max_p1p_changes_full.csv"))

# what changes from AA to AA are more relevant?
stats_by_nrule_prevalent_changes <- nrule_changes %>% 
  group_by(aa_from_to) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count)) %>% 
  write_csv(pp("04-TAB-p1p_aa_change_count_full.csv"))

# what kind of changes are more releavant?
stats_by_nrule_stab_destab <- nrule_changes %>% 
  group_by(status) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count)) %>% 
  write_csv(pp("04-TAB-p1p_stab_destab_count_full.csv"))

### DESTAB
# only destab -> destab subset
nrule_changes_destab <- nrule_changes %>% 
  filter(status == "destab -> destab") %>% 
  arrange(evalue) %>% 
  write_csv(pp("04-TAB-destab_destab_only.csv"))

# in which sites more changes in P1' (destab -> destab only)?
# nrule_destab_most_prots <- nrule_changes_destab %>% 
#   group_by(uniprot, human_site, gname, pname) %>% 
#   summarise(count = n()) %>% 
#   arrange(desc(count)) %>% 
#   write_csv("04-TAB-uni_site_AAchanges_count.csv")

nrule_destab_most_prots <- nrule_changes_destab %>% 
  group_by(uniprot, human_site, gname, pname) %>% 
  summarise(count = n()) %>% 
  #inner_join(., y = nrule_changes %>% # was repr
  inner_join(., y = repr %>% 
               select(uniprot = uni,human_site = fullpep) %>% 
               #select(uniprot,human_site) %>% 
               group_by_all() %>% 
               summarise(total = n())) %>% 
  mutate(percent_of_changes = round(100*count/total, 2)) %>% 
  arrange(desc(percent_of_changes)) %>% 
  write_csv(pp("04-TAB-uni_site_AAchanges_count.csv"))

# what changes from AA to AA are more relevant (destab -> destab only)?
nrule_destab_fromto <- nrule_changes_destab %>% 
  ungroup() %>% 
  mutate(total = n()) %>% 
  group_by(aa_from_to) %>% 
  summarise(percent = round(100*n()/first(total), 2)) %>% 
  arrange(desc(percent)) %>% 
  write_csv(pp("04-TAB-fromto_AAchanges_count.csv"))
