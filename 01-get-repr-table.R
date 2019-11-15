# ## Preprocessing
# All below tables obtained on preliminary stage of analysis. Such stage include next main steps:
#   
# * Obtaining BLAST NR database, extract Vertebrates subset
# * Converting table of sites to BLAST +/-30 amino acids format
# * BLAST processing and filtering some garbage from the result table (detailed information contained in CASPASE_PROJ.org)
# 
# ### Tables obtained after BLAST processing
# 
# * **SITES_60AA_SHORT.csv.gz** -- filter table without unknown proteins and other mess
# * **SITES_60AA_TABLE_ORG_PROT_COUNT.csv** -- org, count sites for each protein, how many proteins represented in BLAST
# * **SITES_60AA_TABLE_UNIQ_ORGS.csv** -- unique organisms, auxilary table which required for identifying lineages

source("00-file_path.R")
library(tidyverse)

fname <- "SITES_60AA"
path <- "./"
fp <- function(pth,fn,suffix) {
  paste0(pth ,fn, "_",suffix)
}

## Blast output filtered table
short <- read_csv(fp(path,fname,"SHORT.csv.gz"))

## Amount of proteins in Blast db for each organims
protein_representation <- read_csv(fp(path,fname,"TABLE_ORG_PROT_COUNT.csv"), 
                                   col_names = c("totprot", "org"))

## Attaching representation info to short table
orgCountFasta <- short %>%
  group_by(org) %>%
  summarise(count= n()) %>%
  left_join(y = protein_representation)

## temporary table with all organisms
saveRDS(orgCountFasta, "./01-orgCountFastaFull.rds")

##
## Scatter plot with representation of organisms
##

orgCountFastaFull <- read_rds("./01-orgCountFastaFull.rds")

## Matches vs coverage plot
g <- ggplot(orgCountFastaFull, aes(x = count, y = totprot))+
  geom_point(alpha=0.5, size = 2) +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') +
  theme_bw() +
  labs(y="Number of sequenced proteins from NR database (log10)", 
       x="Matches with human caspases targets (log10)", 
       title="Count of matches vs Total number of known proteins. 1 point = 1 organism", 
       caption = "Source: midwest") +
  geom_hline(yintercept = 8500, col="red", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = 500, col="red", linetype = "dashed", size = 0.5) +
  theme(axis.text.x = element_text(face="bold", color="black", 
                                   size=14),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=14))+
  annotation_logticks(sides = "lr")

ggsave(pp("01-FIG-count_vs_totcount.pdf"), plot=g, dpi = 600,width = 11, height = 8.5, units = "in")
ggsave(pp("01-FIG-count_vs_totcount.png"), plot=g, dpi = 600,width = 8, height = 5.5, units = "in")

# plot(g)

### REMOVE NON-REPRESENTATIVE ORGANISMS
orgCountFasta <- orgCountFastaFull %>%
  ## ORGANISMS
  filter(log10(count) > 2.5 & log10(totprot) > 3.9) %>%
  arrange(desc(count), desc(totprot))
rm(orgCountFastaFull)


DO_YOU_HAVE_NEW_SPECIES_IN_SHORT_TABLE <- FALSE

if (DO_YOU_HAVE_NEW_SPECIES_IN_SHORT_TABLE) {
  
  # Code below need only in case when we need to append new organisms which absent in **Basic_lineages**. 
  # For example in case when we have new vertebrates which very representative in FASTA file and have a lot of 
  # intersections with human proteins. In this step need to load lineages from wikipedia, ritis, ebi. 
  # The main function **taxon_adapter** takes one parameter long name of organism and return lineage 
  # vector c("Class", "Order", "Family", "Genus") (filename 00-2019-get-lineage.R)
  
  source("../00-2019-get-lineage.R")
  
  ## Lineages partially filled manually
  basic_lineages <- read_csv("Basic_lineages.csv")
  
  ## Append organisms which absent in Basic_lineage table
  absent_organisms <- anti_join(orgCountFasta, basic_lineages) %>% 
    select(org) %>% 
    rowwise() %>% 
    mutate(lin = taxon_adapter(org)) %>% 
    separate(lin, sep = ",", into = c("Class", "Order", "Family", "Genus")) %>% 
    select(Class,Order,Family,Genus,org)
  
  basic_lineages <- bind_rows(basic_lineages, absent_organisms) %>% 
    write_csv("./Basic_lineages.csv")
}


# Basic lineages for all organimsm
basic_lineages <- read_csv("Basic_lineages.csv")

# Create table repr by intersection **short** and **orgCountFasta**. In table only representative organisms. Also was added table **basic_lineages** 
repr1 <- repr %>% 
  inner_join(., y = orgCountFasta) %>% 
  inner_join(., y = basic_lineages)

### Extract orthologs cleavage sites (octet column) and add supplementary information
# Add to repr information about octets - orthologs cleavage sites
# Add difference between human site and found octet (**hamdist**)
# Add Asp in found position **found_type**, nrule amino acid after Asp.

source("00-octet.R")
system.time(
  repr <- repr %>%
    #select(-octet,-found,-nrule, -found_type, -hamdist) %>% 
    rowwise() %>%
    mutate(octet = getOctetByPositionNew(fullpep, qseq, hseq), 
           found = substr(octet, 4, 4), 
           hamdist = hammingStr(fullpep, octet),
           found_type = ifelse(found == "D", 2, ifelse(grepl("D", substr(octet, 2,6)), 1, 0)),
           fullpep_nrule_amino = substr(fullpep,5,5),
           octet_nrule_amino = substr(octet,5,5)) %>%
    ungroup()
)

# save intermediate result and load it again 
saveRDS(repr,"./01-1-repr.rds")
repr <- read_rds("./01-1-repr.rds")

#Load table with stabilazed and distabilazed amino acids. And append this information to table **repr**

## N-end_half-life.csv - stabilizing or destabilazing amino acids table
nrule_table <- read_csv("N-end_half-life.csv")  %>% 
  select(nrule_amino = nrule, nrule, stype = `2011_v_cyto`, ntype = stab_cyto)

## append to fullpep (human site) and octet (ortholog site) columns
repr1 <- left_join(repr, nrule_table, by = c("fullpep_nrule_amino" = "nrule_amino")) %>%
  rename(fullpep_stype = stype, fullpep_ntype = ntype) %>% 
  left_join(., nrule_table, by = c("octet_nrule_amino" = "nrule_amino")) %>% 
  rename(octet_stype = stype, octet_ntype = ntype)

# save intermediate results
saveRDS(repr1,"./01-2-repr.rds")
repr <- read_rds("./01-2-repr.rds")
rm(repr1)

## Add pname and gname information
source("00-2019-uni-fasta.R")
DB <- openDB("./FASTADB.csv")

repr <- inner_join(repr, DB) %>% 
  select(-header, -fasta)

saveRDS(repr,"./01-2-repr.rds")
repr <- read_rds("./01-3-repr.rds")

## Fixing some mistakes in basic lineages
repr1 <- repr %>% 
  mutate(Class = ifelse(Class == "Actinopteri", "Actinopterygii", Class)) %>% 
  mutate(Class = ifelse(Class == "Tetrapoda", "Amphibia", Class))

saveRDS(repr1,"./01-4-repr.rds")
repr <- read_rds("./01-4-repr.rds")
rm(repr1)

# remove sites which do not associated with human
# we have some artifacts in table because blast 
# marked some human sites with evalue more than 1e-16
all_human <- repr %>% 
  select(uni,fullpep,org) %>% 
  filter(grepl("Homo", org)) %>% 
  select(uni,fullpep) %>% 
  distinct()

repr5 <- inner_join(repr, all_human)
saveRDS(repr5,"./01-5-repr.rds")
