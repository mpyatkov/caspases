#s<-"MADKVLKEKRKLFIRSMGEGTINGLLDELLQTRVLNKEEMEKVKRENATVMDKTRALIDSVIPKGAQACQICITYICEEDSYLAGTLGLSADQTSGNYLNMQDSQG--VLSSFPAPQAVQDNPAMPTSSGSEGNVKLCSLEEAQRIWKQKSAEIYPIMDKSSRTRLALIICNEEFDSIPRRTGAEVDITGMTMLLQNLGYSVDVKKNLTASDMTTELEAFAHRPEHKTSDSTFLVFMSHGIREGICGKKHSEQVPDILQLNAIFNMLNTKNCPSLKDKPKVIIIQACRGDSPGVVWFKDSVGVSGNLSLPTTEEFEDDAIKKAHIEKDFIAFCSSTPDNVSWRHPTMGSVFIGRLIEHMQEYACSCDVEEIFRKVRFSFEQPDGRAQMPTTERVTLTRCFYLFPGH"
#p<-"NMQDSQGV"
# hammingStr(p, "NMQZSQG-")
#notExactHamming(p, s)

## For two exact by length patterns
hammingStr <- function(p1,p2) {
  p11 = strsplit(p1,"")[[1]]
  p22 = strsplit(p2,"")[[1]]
  r<-0
  for (i in 1:length(p11)){
    r<-r+ifelse(p11[i] == p22[i], 0, 1)  
  }
  r
}

# hammingStr <- Vectorize(hammingStr)

## For pattern and complete long line
notExactHamming <-function (p, s) {
  pl<-nchar(p)
  ns <-nchar(s)
  wc<-ns-pl+1
  res<-c()
  for (i in 1:wc) {
    res<-c(res, hammingStr(p, substr(s,i,i+pl-1)))
  }
  c(min(res), which.min(res))
}

## AAAAA-D-BBBBB --> AAADBBBB
getOctetWoSpaces <- function(x, pos) {
  ## remove "-"
  remIns <- function(s){return(paste0(strsplit(s,"-")[[1]], collapse = ""))}
  
  #left substring
  lsub <- function(x, n = 4) {
    l <- nchar(x)
    return(substr(x, l-n+1, l))
  }
  
  #right substring
  rsub <- function(x, n = 4) {
    return(substr(x, 1, n))
  }
  
  l <-  substr(x, 1,pos)
  r <-  substr(x, pos+1, 200)
  return(paste0(c(lsub(remIns(l)), rsub(remIns(r))), collapse = ""))
}

#getOctetWoSpaces <- Vectorize(getOctetWoSpaces)

getOctetByPosition <- function(pep, qseq, hseq) {
  
  qseq <- paste0(c("____",qseq,"____"),collapse = "")
  hseq <- paste0(c("____",hseq,"____"),collapse = "")
  
  ## find all D in qseq
  gr <- gregexpr('D', qseq)[[1]]
  
  ## if D not found
  if (gr[1]<0){
    return ("--------")
  }

  allAspPos <- as.vector(gr)
  
  ## compare pep and subseqs in all positions where found D
  hamres <- c()
  for (p in allAspPos) {
    hamres <- c(hamres, hammingStr(pep, getOctetWoSpaces(qseq, p)))
  }
  
  ## get position with minimal hamming dist
  
  minpos <- which.min(hamres)
  # При равных позициях отдавать предпочтение краевым (ham !=0)
  ## if same hamming dist found get one which close to end or start qseq
  if (hamres[minpos] != 0) {
    hamposits <- which(hamres == min(hamres))
    minposits <- allAspPos[hamposits]
    if (length(minposits) > 1) {
      lqseq2 <- nchar(qseq)/2
      minpos <- hamposits[which.max(abs(minposits-lqseq2))]
    }
  }
  
  #return(getOctetWoSpaces(hseq, allAspPos[minpos]))
  return(substr(hseq, allAspPos[minpos]-3,allAspPos[minpos]+4))
}

getOctetByPositionNew <- function(pep, qseq, hseq) {
  
  qseq <- paste0("____",qseq,"____")
  hseq <- paste0("____",hseq,"____")
  
  ## find all D in qseq
  allAspPos <- as.vector(gregexpr('D', qseq)[[1]])
  
  ## if D not found then BLAST alignment awful
  if (allAspPos[1]<0){
    return ("--------")
  }
  
  # return hamming distances for all Asp positions [1,13,16,29] --> [3,0,4,0]
  hamres <- vapply(allAspPos, function(p) {hammingStr(pep, getOctetWoSpaces(qseq, p))},FUN.VALUE = c(1))
  
  ## get position with minimal hamming dist. ex. [3,0,4,0] --> [2,4]
  minpos <-  which(hamres == min(hamres))
  
  ## if only one minpos was found then next two if condition does not need, finded_pos exact what we need
  finded_pos <- allAspPos[minpos]
  
  ## if we found two or many zero hamdist positions, we use hamming for hseq in this positions and take one with the best match
  if (length(minpos) > 1 && hamres[minpos[1]] == 0) {
    # if all minpos equal zero min(hamdist(substr))
    hseqhamdist <- vapply(minpos, function (pos) {hammingStr(pep, substr(hseq, allAspPos[pos]-3,allAspPos[pos]+4))}, FUN.VALUE = c(1))
    finded_pos <- allAspPos[minpos[which.min(hseqhamdist)]]

    # # Smart algo: search in hseq near finded_pos best match (+3/-3 shift)
    # if (min(hseqhamdist) != 0) {
    #   hseqsub <- substr(hseq, finded_pos-5,finded_pos+6)
    #   partAspPosHseq <- as.vector(gregexpr('D', hseqsub)[[1]])
    #   partHamdistHseq <- sapply(partAspPosHseq, function (pos) {hammingStr(pep, substr(hseqsub, pos-3,pos+4))})
    #   if(min(hseqhamdist) > min(partHamdistHseq)) {
    #     local_min_index <- partAspPosHseq[which.min(partHamdistHseq)]
    #     finded_pos <- finded_pos + local_min_index - 6 ## -----(-5)finded_pos(+6)--------
    #   }
    # }
  }
  
  ## if we found two or many non-zero positions then get one which closer to edges (begin or end qseq)
  if (length(minpos) > 1 && hamres[minpos[1]] != 0) {
    lqseq2 <- nchar(qseq)/2
    finded_pos <- allAspPos[minpos[which.max(abs(allAspPos[minpos]-lqseq2))]]
  }
  
  return(substr(hseq, finded_pos-3,finded_pos+4))
}

#getOctetByPosition <- Vectorize(getOctetByPosition)
# ###TESTS
# pep="SSSDAAQA"
# qseq="PKDQFCRRMGQKKQRPARAGQPHSSSDAAQAPAEQPHSSSDAAQAPCPRERCLGPPTTPGP"
# hseq="PDKQFCRRMGQKKQRPARAGQPH--------------SSSDAAQAPCPREHCLGPPTTPGP"
# getOctetByPosition(pep,qseq,hseq)
# getOctetByPositionNew(pep,qseq,hseq)
# 
# microbenchmark(getOctetByPosition(pep,qseq,hseq),
#                getOctetByPositionNew(pep,qseq,hseq))
# 
# pep="TSIDSVVI"
# qseq="VKLQPGQVLHMVNELYPYIVEFEEE"
# hseq="MKLQPGQVLYMVNELYPYIVEFEEE"
# getOctetByPosition("TSIDSVVI", "VKLQPGQVLHMVNELYPYIVEFEEE", "MKLQPGQVLYMVNELYPYIVEFEEE")
# getOctetByPositionNew(pep,qseq,hseq)
# 
# pep = "DVLDVLNE"
# qseq = "DELLGSDDSHDGESESNAKVADVLDVLNEVDEYSGSSEKIDLLASDPHEALICKSERVH"
# hseq = "DEILTSDVSPDGRSESN-----VVEVPIEVDGYSGASEKIALKANDPHGALMC--ERVH"
# getOctetByPosition(pep,qseq,hseq)
# getOctetByPositionNew(pep,qseq,hseq)
# 
# qseq = "DELLGSDDSHDGESESNAKVADVLDVLNEVDEYSGSSEKIDLLASDPHEALICKSERVH"
# hseq = "DKRLTSDVSPDGRSESNVEVPD------EVDGYSGASEKIALKASDPHGALIC--ERVH"
# getOctetByPosition(pep,qseq,hseq)
# getOctetByPositionNew(pep,qseq,hseq)
# 
# pep="DVLDVLNE"
# qseq="SDELLGSDD-SHDGESESNAKVADVLDVLNEVDEYSGSSEKIDLLASDPHEALICKSERVH"
# hseq="SDEMLTSDDDSHDGGSESKTEVAGAVEVPNEVHICSGSSEKIDLMASEPQDALIRERERVH"
# getOctetByPosition(pep,qseq,hseq)
# getOctetByPositionNew(pep,qseq,hseq)
# 
# pep="DLGDGLEN"
# qseq="RQLVKASSTGTESS-DDFEERDP---------DLGDGLENGLGSPFGKWTLSSAAQTHQLRRLRGPAKCR"
# hseq="RRLPRTPSTGTMSSADDLDEREPPSPSDPGLNDIVAEMPNSPG-PFRNTLMSKAAQTHKLRKLRAPSKCR"
# getOctetByPosition(pep,qseq,hseq)
# getOctetByPositionNew(pep,qseq,hseq)
# 
# pep="VGTDEEET"
# qseq="LRLSLNIDPDAKVEEEPEEEPEETAEDTTEDTEQDEDEEMDVGTD-EEETAKESTAEKDEL"
# hseq="LRLSLNIDPDAKVEEEPEEEPEETTEDTTEDTEQDEEDEMDAGADEEEETVKKSTAEKDEL"
# getOctetByPosition(pep,qseq,hseq)
# getOctetByPositionNew(pep,qseq,hseq)
# 
# pep="DEQDGASA"
# qseq = "TTAEEAM-SRGPPPAPEGGSRDE------------QDGASAETEPWAAAVPPEWVPIIQQDIQSQRKVKPQ"
# hseq = "TTAEEAMVSTGEGEAPPSGSSVDGSSQCSSGPQAGQDEVSPEMEPWAAAVPPEWVPIIRHDMLSQRKIKAQ"
# getOctetByPosition(pep,qseq,hseq)
# getOctetByPositionNew(pep,qseq,hseq)
# 
# pep="DETDSGAG"
# qseq="EKWAYNNS-DGGTAWISQES--FDVDETD---SGAGLKWKNVARQLKEDLSSIILLSEEDLQMLVD"
# hseq="EAWGPTNKIDGGTAWLSRESVALDGDEVDMMFSGTD-TWQSLAAQLKQDLASIILLSDSQLQTLVD"
# getOctetByPosition(pep,qseq,hseq)
# getOctetByPositionNew(pep,qseq,hseq)
# 
# pep="EERDGQYF"
# qseq="PAKLNQSGTSVGTDEESDVTQEEERD-GQYFEPVVPLPDLVEVSSGEENEQVVFSHRAKLY"
# hseq="PAKLNQSGTSVGTDEESDVTQEEERDSGQYFEPVVPLPDLVEVSSGEENEQVVFSHRAKLY"
# getOctetByPosition(pep,qseq,hseq)
# getOctetByPositionNew(pep,qseq,hseq)
#  
# pep="DEDDDDVD"
# qseq="MSDAAVDTSSEITTKDLKEKKEVVEEAENGRDAPANGNAENEENGEQEADNEV"
# hseq="MSDTAVDTSSEISTKDLKEKKEVVEETENGRDAPANGNA-NEENGEQEADNEV"
# getOctetByPosition(pep,qseq,hseq)
# getOctetByPositionNew(pep,qseq,hseq)
# 
# pep="STPDFGFG"
# qseq="IAAKIGGDAATTVNNSTP---DFG-FGGQKRQLEDGDQPESKKLASQGDSISSQLGP"
# hseq="IAAKIGGDAAPPVPNNNPPPPDFGGFGGQKRQLEDGDQPESKKLAAQGEALPAPMGP"
# getOctetByPosition(pep,qseq,hseq)
# getOctetByPositionNew(pep,qseq,hseq)
# 
# pep="QARDAQDV"
# qseq="ARAENSQLTERIRSIEALLEAGQARDAQDVQASQAEADQQQTRLKELESQVSGLEKEAIE"
# hseq="ARAENSQLTERIQSIEALLEAGQ---SQDAQARQAETDQQQDRVKELESQVRCMEKEATE"
# getOctetByPosition(pep,qseq,hseq)
# getOctetByPositionNew(pep,qseq,hseq)
#  
# pep="KESDDPMA"
# qseq="SDDPMAYIHFTAEGEVTFKSILFVPTSAPRGLFDEY"
# hseq="SDDPMAYIHFTAEGEVTFKSILFVPTSAPRGLFDEY"
# getOctetByPosition(pep,qseq,hseq)
# getOctetByPositionNew(pep,qseq,hseq)
# 
# pep="DAADAAAA"
# qseq="SGSSFEDMGELHQRLREEEVDADAADAAAAEEEDGEFLGMKGFKGQLSRQVADQMWQAGK"
# hseq="SGSSFEDMGELHQRLREEEVDADAA---AAEEEDGEFLGMKGFKGQLSRQVADQMWQAGK"
# getOctetByPosition(pep,qseq,hseq)
# getOctetByPositionNew(pep,qseq,hseq)
# 
# pep="TQGDEAEA"
# qseq="EKACSLAKTAFDEAIAELDTLSEESYKDSTLIMQLLRDNLTLWTSDTQGDE"
# hseq="DQAISLAKTTFDEAMGDLHTLSEDSYKDSTLIMQLLRDNLTLWTAECAGED"
# getOctetByPosition(pep,qseq,hseq)
# getOctetByPositionNew(pep,qseq,hseq)

