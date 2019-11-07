#!/bin/bash
# processing.sh
# run: sh ./processing.sh ./FASTAFILE

FASTAFILE=$1
OUTFILE="SITES_60AA"
export BLASTDB=..`pwd`

# 1. Extract all potential orthologs for vertebrates
# ~ 3-4 hours depending on hardware

blastp -db nr.7742 -query $FASTAFILE \
       -outfmt "6 qaccver saccver stitle evalue score pident qseq sseq" \
       -out $OUTFILE.tsv -num_alignments 8000 \
       -num_threads 8 -evalue 1e-16

# 2. Filtering table FASTAFILE.tsv, removing duplicates
# output: $FASTAFILE_TABLE_UNIQ_ORGS.csv  -- list of unique organisms
# output: $FASTAFILE_1_SHORT.csv.gz       -- cleaned FASTAFILE.tsv table 
Rscript --vanilla 00-remote-getshort.R $OUTFILE.tsv

# 3. The number of sequences associated with each 
# organism in NR BLAST database (proteom representativeness)

## grep -Po '(?<=\[).*(?=\]$)'            -- extract  ex. [Homo sapiens]$
## grep -E -v "\.|\[|\]|\,|=|-|\(|\/"     -- parse garbage
## sed -e 's/^ *//;s/ /,/'                -- remove spaces which used uniq
## grep -f $OUTFILE_TABLE_UNIQ_ORGS.csv -- get organisms from file

cat nr.7742.fa | grep ">" | grep -Po '(?<=\[).*(?=\]$)' \
    | cut -d" " -f1,2 | sort | uniq -c | sort -n \
    | grep -E -v "\.|\[|\]|\,|=|-|\(|\/"  | sed -e 's/^ *//;s/ /,/' \
    | grep -f "$OUTFILE"_TABLE_UNIQ_ORGS.csv \
      > "$OUTFILE"_TABLE_ORG_PROT_COUNT.csv
