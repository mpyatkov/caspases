library(ritis)
library(xml2)
library(stringr)
library(rjson)
library(wikitaxa)
library(htmltab)

## This is 5 functions which facilitate obtaining the lineage from wikipedia, ritis, ebi
## Download from www.ebi.ac.uk site
ebi_lineage <- function(org) {
  hd <- "http://www.ebi.ac.uk/ena/data/view/Taxon:"
  orgstr <- paste0(str_split(org," ",simplify = TRUE), collapse = "+")
  tl <- "&display=xml"
  url <- paste0(c(hd,orgstr,tl), collapse = "")
  
  request <- GET(url) %>%
    content("text", encoding = "UTF-8") %>%
    read_xml() %>% 
    xml_find_all("//lineage/taxon[@rank]")
  
  ranks <- request %>% 
    xml_find_all("@rank") %>% xml_contents() %>% xml_text(.)
  
  scientificName <- request %>% 
    xml_find_all("@scientificName") %>% xml_text(.) 
  
  lineage <- c("class", "order", "family", "genus")
  
  idx <- sapply(lineage, function(x) {which(ranks == x)})
  scientificName[idx]
  
  # # xml_attr("scientificName") %>%
  # paste0(collapse = ";") %>%unname(.)
}

## Using library(wikitaxa)
## Not always work on scientific names
## What the hell with output rank names ("Classis", "Ordo", "Familia", "Genus")?
wikitaxa_lineage <- function(org) {
  #common_name <- wt_wikipedia(name = "Bos mutus", wiki = "en")$common_names$name
  common_name <- wt_wikicommons_search(query = org)$query$search$title[1]
  r <- wt_wikicommons(common_name)$classification
  ranks <- r$rank
  names <- r$name
  lineage <- c("Classis", "Ordo", "Familia", "Genus")
  idx <- sapply(lineage, function(x) {which(ranks == x)})
  str_remove_all(names[idx], "[:][:space:]")
}
# wiki_lineage("Calypte anna")

## Parsing Template:Taxonomy page from organism
## problems: sometimes such pages don't exist
wiki_taxon <- function(org) {
  wurl <- "https://en.wikipedia.org/wiki/Template:Taxonomy/"
  #orgstr <- paste0(str_split(org," ",simplify = TRUE), collapse = "+")
  t <- htmltab(doc = paste0(wurl,org), which = "//table[1]", rm_nodata_cols = F)[,1:2]
  names(t) <- c("rank", "name")
  ranks <- t$rank
  names <- t$name
  lineage <- c("Class:", "Order:", "Family:", "Genus:")
  idx <- sapply(lineage, function(x) {which(ranks == x)})
  names[idx]
}

## Extracting lineage using ITIS site
## library(ritis)
## problems: not all species presented
ritis_lineage <- function(name) {
  tsn<- itis_search(q = paste0("nameWOInd:",name))$tsn[1]
  h <- hierarchy_full(tsn = tsn)
  lineage <- c("Superclass", "Order", "Family", "Genus")
  idx <- sapply(lineage, function(x) {which(h$rankname == x)})
  # paste0(h$taxonname[idx], collapse = ";")
  h$taxonname[idx]
}

## Extracting lineage from the main page about organism
main_wiki_lineage <- function(org) {
  wurl <- "https://en.wikipedia.org/w/api.php"
  r <- GET(wurl, query= list(action = "opensearch", search=org, format="json"))
  txt <- r %>% content("text", encoding = "UTF-8") 
  page_url <- fromJSON(txt)[4][[1]][[1]]
  t <- htmltab(doc = page_url, which = "//table[contains(@class, 'infobox biota')]", rm_nodata_cols = F)#[,1:2]
  names(t) <- c("rank", "name")
  ranks <- t$rank
  names <- t$name
  lineage <- c("Class:", "Order:", "Family:", "Genus:")
  idx <- sapply(lineage, function(x) {which(ranks == x)})
  #print(idx)
  # str_remove_all(names[idx], "[:][:space:]")
  #print(org)
  str_split(names[idx], pattern ="(?=[A-Z])", n = 3, simplify = T)[,2]
}
## Adapter which selected next function if previous one did not work
taxon_adapter <- function(fullorg) {
  shortorg <- str_split_fixed(fullorg, pattern = " ", n = 2)[1]
  r = tryCatch({wiki_taxon(shortorg)},
               error = function(e) {
                 tryCatch({wikitaxa_lineage(fullorg)},
                          error = function(e1) {
                            print("wikitaxa_lineage")
                            tryCatch({main_wiki_lineage(fullorg)},
                                     error = function(e2) {
                                       print("main_wiki_lineage")
                                       tryCatch({ebi_lineage(shortorg)},
                                                error = function (e3){
                                                  print("ebi_lineage")
                                                  ritis_lineage(shortorg)
                                                })
                                     })
                          })
               } , finally = {print(fullorg)})
  paste0(r,collapse = ",")
}

