library(tidyverse)
source("00-file_path.R")
# Working with Hydrophobicity and domain structure of 60AA

#  Hydrophobicity index (kcal/mol) of amino acids in a distribution from non-polar to polar at pH = 7
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-S16-S8/tables/2

dfix <- "
AA,ix
I,4.92
L,4.92 	
V,4.04 	
P,4.04 	
F,2.98 	
M,2.35 	
W,2.33
A,1.81 	
C,1.28
G,0.94
Y,-0.14
R,-14.92
D,-8.72
E,-6.81
N,-6.64
K,5.55
T,-2.57
S,-3.40
H,-4.66
Q,-5.54
"
df<- read.csv(text = dfix, stringsAsFactors = F) %>% 
  as_tibble() %>% 
  write_csv(pp("052-TAB-AA_hidrophob_ix.csv"))


## remove chains of repetitive symbols which chain length less then `minseq` 
clean_from_garbage <- function(s, minseq = 3){
  
  removeNPrevElements <- function(st, cnt, from, elsym) {
    for (n in seq(1,cnt)) {
      st[from-n] <- elsym
    }
    st
  }
  
  cs <- str_split(s,"")[[1]]
  prev <- 1
  ep <- cs[prev]
  counter <- 1
  for (nxt in seq(2, length(cs))) {
    if (cs[nxt] == cs[prev]) {
      counter <- counter + 1
      prev = nxt
      next()
    }
    if (cs[nxt] != cs[prev]) {
      if (counter <= minseq) {
        # for (n in seq(1,counter)) {
        #   cs[nxt-n] <- "."
        # }
        # cs <- removeNPrevElements(cs, counter, nxt, cs[nxt])
        cs <- removeNPrevElements(cs, counter, nxt, "C")
      }
      counter <- 1
      prev = nxt
    }
    if ((nxt+1) == length(cs)) {
      # cs <- removeNPrevElements(cs, counter+1, length(cs)+1, cs[nxt-1])
      cs <- removeNPrevElements(cs, counter+1, length(cs)+1, "C")
    }
  }
  paste0(cs, collapse = "")
}

# get substring and count percent of symbol to whole length
proc <- function(st, sym, start = 0, end = 0) {
  
  # TODO: make more robust
  if (start == 0 || end == 0) {
    end <- str_length(st)
  }
  
  sub_s <- str_sub(st, start, end)
  sub_c <- str_count(sub_s, sym)
  round(100*sub_c/str_length(sub_s),2)
}

repr <- read_rds("./01-5-repr.rds")

hydrphobix <- read_csv(pp("052-TAB-AA_hidrophob_ix.csv"))

# represent sequence as sum of hydrophobity index (hix)
# where hix calculated for each AA in sequence
hydrophobix_sum <- function(s){
  cs <- str_split(s,"")[[1]]
  h <- tibble(AA =  cs) %>% 
    inner_join(hydrphobix, by = "AA") %>%
    select(ix) %>% 
    pull()
  sum(h)
}

# procent of Asp in P1 for each site
# procD <- read_csv("04-TAB-procD.csv") %>% 
procD <- read_csv("03-TAB-procD.csv") %>% 
  select(uni,fullpep,dproc)

# most common type of P1'
# maxnrule <- read_csv("04-TAB-uni_fullpep_maxnrule.csv")
maxnrule <- read_csv("04-TAB-uni_fullpep_maxnrule.csv")

# mean hamming distance for each site
# meanHD <- read_csv("04-TAB-meanhd_for_uni_fullpep.csv") %>% 
meanHD <- read_csv("03-TAB-meanhd_for_uni_fullpep.csv") %>% 
  select(uni, fullpep, meanhd) %>% 
  mutate(meanhd = round(meanhd,2))

## conservation rate
consRate <- repr %>% 
  select(uni, fullpep) %>% 
  group_by_all() %>% 
  summarise(conserv = n())

TS1 <- read_csv("Table_S1_Human_caspase_targets.csv") %>% 
  select(uni,fullpep, Q3)

# add qseq and domain, clean Q3
qseq_domain <- repr %>% 
  filter(grepl("Homo", org), !grepl("-",qseq)) %>%
  select(-org) %>%
  rowwise() %>% 
  mutate(slen = str_length(qseq)) %>% 
  filter(slen == 60) %>% 
  select(uni, fullpep, qseq) %>%
  inner_join(y = TS1, by = c("uni", "fullpep")) %>% 
  mutate(Q3clean = clean_from_garbage(Q3)) %>% 
  select(-Q3) %>% 
  rename(Q3 = Q3clean)

# join all tables
cortable <- qseq_domain %>% 
  inner_join(procD) %>% 
  inner_join(consRate) %>% 
  inner_join(maxnrule) %>% 
  inner_join(meanHD)

# calculate statistics for [left,center,right] for each 60AA
eps <- 0.0001
cortable_with_proc <- cortable %>% 
  rowwise() %>% 
  mutate( left_hydr20 = hydrophobix_sum(str_sub(qseq,1,20)),
          center_hydr20 = hydrophobix_sum(str_sub(qseq,21,40)),
          right_hydr20 = hydrophobix_sum(str_sub(qseq,41,60)),
          center_hydr8 = hydrophobix_sum(str_sub(qseq,27,34)),
          left_coil20 = proc(Q3,"C",1,20),
          center_coil20 = proc(Q3,"C",21,40),
          right_coil20 = proc(Q3,"C",41,60),
          center_coil8 = proc(Q3,"C",27,34),
          DN           = ifelse(str_count(str_sub(Q3,30,31), "C") == 2, 1,0),
          hydr_lr = (left_hydr20+right_hydr20)*0.5,
          coil_lr = (left_coil20+right_coil20)*0.5,
          more_coil = (eps+center_coil20)/(eps+coil_lr)) %>% 
  select(-qseq,-Q3)


# We moved our axis from {-N, M} to {1, M + N + 1} 
# for more convient calculation of result value. 
# And after that center_hydr20 and hydr_lr20 columns were scaled.

axes_shift_flangs <- cortable_with_proc$hydr_lr[which.min(cortable_with_proc$hydr_lr)]
axes_shift_center <- cortable_with_proc$center_hydr20[which.min(cortable_with_proc$center_hydr20)]
axes_shift <- abs(min(axes_shift_flangs, axes_shift_center))+1

cortable_scaled <- cortable_with_proc %>% 
  mutate(hydr_lr_scl = hydr_lr+axes_shift,
         hydr_center20_scl = center_hydr20+axes_shift,
         O_est = hydr_center20_scl/hydr_lr_scl) %>% 
  select(-hydr_lr_scl, -hydr_center20_scl, -hydr_lr, -coil_lr) %>% 
 write_csv(pp("052-TAB-cortable.csv"))

## DN must be 1 (P1 and P1' == C)
## coil percent >= 50 in center8 (HHCCCCHH <- good structure of domain)
## hydrophob < 1 (we measure ratio bw center20 and flangs20)

good_domain_str <- cortable_scaled %>% 
  select(uni,fullpep,DN,O_est,center_coil8) %>% 
  filter(DN == 1, O_est < 1, center_coil8 >=50) %>% 
  select(uni, fullpep) %>% 
  write_csv(pp("052-TAB-inner_cc8_hydr_cc2.csv"))
