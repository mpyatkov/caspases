## Aux. tables
source("00-file_path.R")

require(tidyverse)

### amino_count_procent

repr <- read_rds("./01-5-repr.rds")

#### TASK: Amino, count, procent
amino_count_procent <- repr %>% 
  filter(!grepl("Homo",org)) %>% 
  select(found) %>% 
  mutate(tot = n()) %>% 
  group_by(found) %>% 
  summarise(tot = first(tot), count = n(), procent= round(100*count/tot, digits = 5)) %>% 
  select(-tot) %>% 
  arrange(desc(count)) %>% 
  write_excel_csv(pp("03-TAB-amino_count_procent.csv"))
  
#### TASK: D estim, HammingDistance, count, procent
estim_hd_count_procent <- repr %>% 
  filter(!grepl("Homo",org)) %>% 
  select(found_type,hamdist) %>% 
  mutate(tot = n()) %>%
  group_by(found_type, hamdist) %>% 
  summarise(tot = first(tot), count = n(), procent= round(100*count/tot, digits = 2)) %>% 
  select(-tot) %>% 
  arrange(desc(found_type), hamdist,desc(count)) %>% 
  ungroup() %>% 
  mutate(found_type = ifelse(found_type == 2, "Asp_in_right_position", ifelse(found_type == 1, "Asp_in_near_sites", "No_asp")))%>% 
  write_excel_csv(pp("03-TAB-estim_hd_count_procent.csv"))
  
#### TASK: uni, gname, fullpep, count

uni_gname_fullpep_count <- repr %>% 
  filter(!grepl("Homo",org)) %>% 
  select(uni,gname, fullpep) %>% 
  group_by(uni,gname, fullpep) %>% 
  summarise(count = n()) %>% 
  arrange(desc(count)) %>% 
  write_excel_csv(pp("03-TAB-uni_gname_fullpep_count.csv"))


#### TASK: d_proc (06-dproc...)

# percent of Asp for each protein
# create a table with columns: pname, D%. SHORT, Evalue 1e-16, org !=  “Homo sapiens”

procD <- repr %>% 
  filter(!grepl("Homo",org)) %>% 
  select(uni, fullpep, found_type) %>% 
  group_by(uni, fullpep) %>% 
  mutate(found_type = if_else(found_type == 2,1,0)) %>% # changing D-estimate from 0,1,2 to 0,1 
  summarise(dproc = round(100*sum(found_type)/n(), digits = 2), 
            count = n()) %>% 
  arrange(dproc) %>% 
  write_excel_csv(pp("03-TAB-procD.csv"))

# take attention to znf598, only human site presented. In other organisms mutation in this position

### uni, fullpep, mean hamming distantce
meanHD <- repr %>% 
  filter(!grepl("Homo", org)) %>% 
  select(uni,fullpep,hamdist) %>% 
  group_by(uni,fullpep) %>% 
  summarise(meanhd = mean(hamdist), count = n()) %>% 
  ungroup() %>% 
  arrange(meanhd) %>% 
  write_excel_csv(pp("03-TAB-meanhd_for_uni_fullpep.csv"))


#### TASK: Agard data

# agard.csv contains filtered information from agard paper
agard <- read_csv("agard.csv", col_names = T) %>%
  select(uni, after_d4, elution, expect)

agard_cor <- repr %>% 
  filter(!grepl("Homo",org), found_type == 2) %>%
  group_by(uni,gname, fullpep) %>% 
  summarise(count = n()) %>%
  ungroup() %>% 
  mutate(after_d4 = str_sub(fullpep,5,8)) %>% 
  select(-fullpep) %>% 
  
  ## connect with agard data
  left_join(., agard, by=c("uni", "after_d4")) %>% 
  filter(!is.na(expect)) %>% 
  ungroup() %>% 
  arrange(count, elution) %>% 
  write_csv(pp("03-TAB-agard_cor.csv")) #,append = F, col_names = T)

# correlation between our data and agard data
acp <- ggplot(agard_cor, aes(x = count, y = elution))+
  geom_point(alpha=0.5, size = 2) +
  # scale_x_continuous(trans='log10') +
  # scale_y_continuous(trans='log10') +
  theme_bw() +
  labs(y="Elution time", 
       x="Conservation", 
       title="Correlation with Agard data. 1 point = 1 cleavage site") 

ggsave(filename = pp("03-FIG-agard.png"), plot = acp, width = 173, height = 160, units = "mm")

#### TASK: for each uni,fullpep count which max frequent octet_ntype
max_freq_octet_ntype <- repr %>% 
  select(uni,fullpep,octet_ntype) %>% 
  group_by_all() %>% 
  summarise(cnt = n()) %>% 
  filter(!is.na(octet_ntype)) %>% 
  ungroup() %>% 
  group_by(uni,fullpep) %>% 
  summarise(l = list(cnt), o = list(octet_ntype)) %>% 
  ungroup() %>% 
  rowwise() %>% 
  mutate(octet_ntype_max = o[which.max(l)]) %>% 
  select(uni, fullpep, octet_ntype_max) %>% 
  write_csv(pp("03-TAB-uni_fullpep_maxnrule.csv"))

#### TASK: Glutamat vs Aspartat in P1 position

# repr <- read_rds("./01-5-repr.rds")
strict_e <- repr %>%
  filter(org != "Homo sapiens", found == "E") # hamdist < 3

## Counting frequency of the E-cuts within the vertebrates
cons <- strict_e %>%
  select(uni, fullpep, pname) %>%
  group_by_all() %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>% 
  write_excel_csv(pp("03-TAB-Strict_conservation_range_E.csv"))

## Counting frequency of e-cuts within species to see trends in D->E substitutions
spece <- strict_e %>%
  select(uni, fullpep, pname, org) %>%
  group_by(org) %>%
  summarise(count = n()) %>%
  arrange(desc(count)) %>% 
  write_excel_csv(pp("03-TAB-Species_E.csv"))
