setwd("/projects/SARS-CoV-2/msa/")
library(msa)

seq <- Biostrings::readAAStringSet("txid2697049_AA_SGlycoprotein_sequence.txt")
seq.sg <- seq[grep("[protein=surface glycoprotein]", seq@ranges@NAMES, fixed = TRUE)]
Biostrings::writeXStringSet(seq.sg, "txid2697049_AA_SGlycoproteinONLY_sequence.txt")

seq.sg.align <- Biostrings::readAAMultipleAlignment("txid2697049_AA_SGlycoproteinONLY_sequence_align.txt", format = "clustal")

seq.sg.long <- seq.sg[seq.sg@ranges@width != 35]
Biostrings::writeXStringSet(seq.sg.long, "txid2697049_AA_SGlycoproteinONLY_long_sequence.txt")

# Genbank (17MAR2020): https://www.ncbi.nlm.nih.gov/nuccore/?term=(txid2697049)+AND+%22surface+glycoprotein%22%5BProtein+Name%5D++AND+%22complete+genome%22
## 81 items
seq <- Biostrings::readAAStringSet("txid2697049_AA_complete_genome_sequence.txt")
seq.sg <- seq[grep("[protein=surface glycoprotein]", seq@ranges@NAMES, fixed = TRUE)]
seq.sg.long <- seq.sg[seq.sg@ranges@width != 35]

Biostrings::writeXStringSet(seq.sg.long, "txid2697049_AA_SGlycoproteinONLY_long_sequence.txt")

# MSA with Clustal Omega
seq.sg.align <- Biostrings::readAAMultipleAlignment("txid2697049_AA_complete_genome_surface_glycoprotein_ONLY_sequence_align.txt", format = "clustal")


msaPrettyPrint(seq.sg.align, output="pdf", y=c(164, 213),subset=c(1:6), showNames="none", showLogo="none",consensusColor="ColdHot", showLegend=FALSE,askForOverwrite=FALSE)
msaPrettyPrint(seq.sg.align, output="pdf", showNames="left", showLogo="top",consensusColor="ColdHot", showLegend=TRUE,askForOverwrite=FALSE)

## Clustal Omega in RStudio
msaAlign <- msa(seq.sg.long, "ClustalOmega")
msaAlign2 <- msaConvert(msaAlign, type="seqinr::alignment")

# Pairwise Distances from Aligned Protein or DNA/RNA Sequences
d <- seqinr::dist.alignment(msaAlign2, "identity")
attr(d, "Labels") <- gsub(" \\[.+CDS\\]$", "", attr(d, "Labels"))

pheatmap::pheatmap(as.matrix(d), fontsize_row = 5, fontsize_col = 5)

#------------------------
# Genbank (18MAR2020): https://www.ncbi.nlm.nih.gov/nuccore/NC_045512,MN908947,MN970003,MN970004,MN938384,MN938385,MN938386,MN938387,MN938388,MN938389,MN938390,MN975262,MN975263,MN975264,MN975265,MN975266,MN975267,MN975268,MN985325,MN988713,MN994467,MN994468,MN997409,MN988668,MN988669,MN996527,MN996528,MN996529,MN996530,MN996531,MT007544,MT008022,MT008023,LR757995,LR757996,LR757997,LR757998,MT019529,MT019530,MT019531,MT019532,MT019533,MT020781,MT020880,MT020881,MT027062,MT027063,MT027064,LC522350,MT039873,MT039887,MT039888,MT039890,MT042773,MT042774,MT042775,MT042776,MT042777,MT042778,MT044257,MT044258,MT050414,MT050415,MT050416,MT050417,MT049951,LC523807,LC523808,LC523809,MT066157,MT066158,MT066159,MT066175,MT066176,MT072667,MT072668,MT072688,MT081059,MT081060,MT081061,MT081062,MT081063,MT081064,MT081065,MT081066,MT081067,MT081068,MT093571,MT093631,MT106052,MT106053,MT106054,MT111895,MT111896,MT118835,MT123290,MT123291,MT123292,MT123293,LC528232,LC528233,MT126808,MT127113,MT127114,MT127115,MT127116,MT135041,MT135042,MT135043,MT135044,MT152824,MT152900,MT012098,MT050493,MT121215,MT159705,MT159706,MT159707,MT159708,MT159709,MT159710,MT159711,MT159712,MT159713,MT159714,MT159715,MT159716,MT159717,MT159718,MT159719,MT159720,MT159721,MT159722,MT159778,MT066156,MT161607,MT163712,MT163714,MT163715,MT163716,MT163717,MT163718,MT163719,MT163720,MT163721,MT163737,MT163738,LC529905,MT184907,MT184908,MT184909,MT184910,MT184911,MT184912,MT184913,MT186676,MT186677,MT186678,MT186679,MT186680,MT186681,MT186682,MT187977,MT188339,MT188340,MT188341,MT192758,MT192759,MT192765,MT192772,MT192773,MT198651,MT198652,MT198653
## 174 items
seq <- Biostrings::readAAStringSet("20200318_EntrezNucleotide_SARSCoV2_AA_sequences.fasta")
seq.sg <- seq[grep("[protein=surface glycoprotein]", seq@ranges@NAMES, fixed = TRUE)]
length(seq.sg)
# [1] 98

require(msa)
## Clustal Omega in RStudio
msaAlign <- msa(seq.sg, "ClustalOmega")
msaAlign2 <- msaConvert(msaAlign, type="seqinr::alignment")

# Pairwise Distances from Aligned Protein or DNA/RNA Sequences
d <- seqinr::dist.alignment(msaAlign2, "identity")
attr(d, "Labels") <- gsub(" \\[.+CDS\\]$", "", attr(d, "Labels"))

pdf("20200318_EntrezNucleotide_SARSCoV2_AA_sequences_MSA_dist_mat.pdf")
pheatmap::pheatmap(as.matrix(d), fontsize_row = 4, fontsize_col = 4)
dev.off()

msaAlign@unmasked@ranges@NAMES <- gsub("lcl\\|([MNL][TNC]_?[0-9\\.]+)_(prot_[A-Z_0-9\\.]+) .+$", "\\1_\\2", msaAlign@unmasked@ranges@NAMES)
seq.sg.align.ids <- gsub("lcl\\|", "", seq.sg.align@unmasked@ranges@NAMES)
msaPrettyPrint(msaAlign, output="pdf", showNames="left", showLogo="top",consensusColor="ColdHot", showLegend=TRUE,askForOverwrite=FALSE)

saveWidth <- getOption("width")
options(width=5000)
sink("20200318_EntrezNucleotide_SARSCoV2_AA_sequences_align.txt")
print(msaAlign, show="complete", halfNrow=NA, nameWidth=50)
sink()
options(width=saveWidth)

system("sed -i 's/^ //g' 20200318_EntrezNucleotide_SARSCoV2_AA_sequences_align.txt")

# get mismatches in table format
aln.tab <- read.delim("20200318_EntrezNucleotide_SARSCoV2_AA_sequences_align.txt", skip=3, sep=" ", strip.white = TRUE, stringsAsFactors = FALSE, header=FALSE)
require(stringr)
cons <- unlist(str_split(aln.tab$V2[nrow(aln.tab)], ""))
mm.tab <- data.frame(sequence=character(), position = numeric(), mutation = character(), consensus = character())
for (i in 1:nrow(aln.tab)-1) {
  s <- unlist(str_split(aln.tab$V2[i], ""))
  mm <- which(s != cons)
  if(length(mm)<100){
    for (m in mm){
      tmp <- data.frame(sequence=aln.tab$V3[i], position=m, mutation=s[m], consensus=cons[m])
      mm.tab <- rbind(mm.tab, tmp)
    }
  } else {
    for (mi in 1:(length(mm)-1)){
      if(s[mm[mi]] == s[mm[mi+1]]) {next}
      tmp <- data.frame(sequence=aln.tab$V3[i], position=mm[mi], mutation=s[mm[mi]], consensus=cons[mm[mi]])
    }
  }
}

require(tidyverse)
mm.tab %>% filter(mutation != "X") %>% distinct(sequence)
write.csv2(mm.tab, "20200318_EntrezNucleotide_SARSCoV2_AA_sequences_align.csv", quote = FALSE, row.names = FALSE)

#------------------------
# GISAID (25MAR2020)

gisaid <- Biostrings::readDNAStringSet("/scratch/info/projects/SARS-CoV-2/gisaid/gisaid_cov2020_sequences.fasta")
refseq <- Biostrings::readDNAStringSet("../RefSeq_surface_glycoprotein.fasta")

# extract gene S / surface glycoprotein
gisaid.sg <- Biostrings::DNAStringSet()
for(g in 1:length(gisaid)){
  tmp.aln <- Biostrings::pairwiseAlignment(refseq, gisaid[g], type="local")
  tmp.sg <- subseq(gisaid[g], tmp.aln@subject@range)
  gisaid.sg <- c(gisaid.sg, tmp.sg)
}

gisaid.sg.amb <- replaceAmbiguities(gisaid.sg)

refseq.aa <- translate(refseq)
gisaid.sg.aa <- translate(gisaid.sg, if.fuzzy.codon = "solve")
writeXStringSet(gisaid.sg.aa, "../RefSeq_surface_glycoprotein_AA.fasta")

msaAlign.all <- msa(c(seq.sg, gisaid.sg.aa), "ClustalOmega")
msaAlign2.all <- msaConvert(msaAlign, type="seqinr::alignment")

# Pairwise Distances from Aligned Protein or DNA/RNA Sequences
d.all <- seqinr::dist.alignment(msaAlign2.all, "identity")
attr(d.all, "Labels") <- gsub(" \\[.+CDS\\]$", "", attr(d, "Labels"))

# ClustalW alignment
WmsaAlign.gisaid <- msa(gisaid.sg.aa, "ClustalW")

saveWidth <- getOption("width")
options(width=5000)
sink("20200330_GISAID_SARSCoV2_AA_sequences_align.txt")
print(WmsaAlign.gisaid, show="complete", halfNrow=NA, nameWidth=100)
sink()
options(width=saveWidth)

system("sed -i -r 's/^ +//g' 20200330_GISAID_SARSCoV2_AA_sequences_align.txt")
system("sed -i -r 's/\\/([A-Z][a-z]+) ([A-Z][a-z]+)\\//\\/\\1\\2\\//g' 20200330_GISAID_SARSCoV2_AA_sequences_align.txt")

# get mismatches in table format
aln.tab.gisaid <- read.delim("20200330_GISAID_SARSCoV2_AA_sequences_align.txt", skip=3, sep=" ", strip.white = TRUE, stringsAsFactors = FALSE, header=FALSE)
require(stringr)
cons.gisaid <- unlist(str_split(aln.tab.gisaid$V2[nrow(aln.tab.gisaid)], ""))
mm.tab.gisaid <- data.frame(sequence=character(), position = numeric(), mutation = character(), consensus = character())
for (i in 1:nrow(aln.tab.gisaid)-1) {
  s <- unlist(str_split(aln.tab.gisaid$V2[i], ""))
  mm <- which(s != cons.gisaid)
  if(length(mm)<100){
    for (m in mm){
      tmp <- data.frame(sequence=aln.tab.gisaid$V3[i], position=m, mutation=s[m], consensus=cons.gisaid[m])
      mm.tab.gisaid <- rbind(mm.tab.gisaid, tmp)
    }
  } else {
    for (mi in 1:(length(mm)-1)){
      if(s[mm[mi]] == s[mm[mi+1]]) {next}
      tmp <- data.frame(sequence=aln.tab.gisaid$V3[i], position=mm[mi], mutation=s[mm[mi]], consensus=cons.gisaid[mm[mi]])
    }
  }
}

write.csv2(mm.tab.gisaid, "20200330_GISAID_SARSCoV2_AA_sequences_align.csv", quote = FALSE, row.names = FALSE)

require(tidyverse)
mm.tab.gisaid %>% filter(mutation != "X",
                         !grepl("bat|pangolin", sequence)) %>% distinct(sequence)
# number of mutations in one sequence
mm.tab.gisaid %>% filter(mutation != "X") %>% group_by(sequence) %>% summarise(obs = n()) %>% arrange(desc(obs))
hist(mm.tab.gisaid %>% filter(mutation != "X", 
                              !grepl("bat|pangolin", sequence)) %>% group_by(sequence) %>% summarise(obs = n()) %>% arrange(desc(obs)) %>% pull(obs), xlab = "# mutations in one sequence", freq = F, main="Histogram")
# mutations in RBD (pos 322-533)
mm.tab.gisaid %>% filter(mutation != "X",
                         !grepl("bat|pangolin", sequence),
                         (position >= 322 & position <= 533)) %>% group_by(sequence) %>% summarise(obs = n()) %>% arrange(desc(obs)) %>% distinct(sequence)

# recurrent mutations 
mm.tab.gisaid %>% filter(mutation != "X",
                         !grepl("bat|pangolin", sequence)) %>% 
  group_by(position, consensus, mutation) %>% summarise(obs = n()) %>% arrange(desc(obs))


#------------------------
require(tidyverse)
require(Biostrings)
# GISAID (01APR2020)

gisaid.0104 <- Biostrings::readDNAStringSet("/scratch/info/projects/SARS-CoV-2/gisaid/gisaid_cov2020_sequences_01APR2020_human_host_low_cov_excl.fasta")
# warning invalid one-letter codes => refers to spaces in fasta seqs => removed in DNA string set
refseq <- Biostrings::readDNAStringSet("../RefSeq_surface_glycoprotein.fasta")

# extract gene S / surface glycoprotein
make_frame <- function(x) {
  mx <- ifelse(x %% 3 == 1, 0, 
               ifelse(x %% 2 == 1, 2, 
                      ifelse(x %% 3 == 0, 1)))
  return(mx)
}

# gisaid.0104.sg <- Biostrings::DNAStringSet()
# for(g in 1:length(gisaid.0104)){
#   print(paste0(g, ": ", gisaid.0104@ranges@NAMES[g]))
#   tmp.aln <- Biostrings::pairwiseAlignment(refseq, gisaid.0104[g], type="local")
#   # extract sequence in same reading frame as refseq seq:
#   tmp.sg <- subseq(gisaid.0104[g], start = tmp.aln@subject@range@start + make_frame(tmp.aln@pattern@range@start),
#                    width = tmp.aln@subject@range@width - (tmp.aln@pattern@range@width %% 3))
#   gisaid.0104.sg <- c(gisaid.0104.sg, tmp.sg)
# }

library(Biostrings)
library(doParallel)

# extract_sg <- function(st=1, en=2) {
#   cl <- makeCluster(6)
#   registerDoParallel(cl)
#   # ex <- Biostrings::DNAStringSet()
#   ex <-
#     foreach(g=st:en, .inorder = FALSE, .errorhandling = "remove", .combine = c) %dopar% {
#     
#       library(Biostrings)
#       tmp.aln <- Biostrings::pairwiseAlignment(refseq, gisaid.0104[g], type="local")
#       # extract sequence in same reading frame as refseq seq:
#       tmp.sg <- subseq(gisaid.0104[g], start = tmp.aln@subject@range@start + make_frame(tmp.aln@pattern@range@start),
#                        width = tmp.aln@subject@range@width - (tmp.aln@pattern@range@width %% 3))
#       # gisaid.0104.sg <- c(gisaid.0104.sg, tmp.sg)
#       tmp.sg
#     }
#   
#   stopCluster(cl)
# }

cl <- makeCluster(6)
registerDoParallel(cl)

st <- 1
en <- 589
gisaid.0104.sg.1 <- foreach(g=st:en, .inorder = FALSE, .errorhandling = "remove", .combine = c) %dopar% {
  
  library(Biostrings)
  tmp.aln <- Biostrings::pairwiseAlignment(refseq, gisaid.0104[g], type="local")
  # extract sequence in same reading frame as refseq seq:
  tmp.sg <- subseq(gisaid.0104[g], start = tmp.aln@subject@range@start + make_frame(tmp.aln@pattern@range@start),
                   width = tmp.aln@subject@range@width - (tmp.aln@pattern@range@width %% 3))
  # gisaid.0104.sg <- c(gisaid.0104.sg, tmp.sg)
  return(tmp.sg)
}

Biostrings::writeXStringSet(gisaid.0104.sg.1, "gisaid_sg_1.fasta")
rm(gisaid.0104.sg.1)

st <- 590
en <- 589*2
gisaid.0104.sg.2 <- foreach(g=st:en, .inorder = FALSE, .errorhandling = "remove", .combine = c) %dopar% {
  
  library(Biostrings)
  tmp.aln <- Biostrings::pairwiseAlignment(refseq, gisaid.0104[g], type="local")
  # extract sequence in same reading frame as refseq seq:
  tmp.sg <- subseq(gisaid.0104[g], start = tmp.aln@subject@range@start + make_frame(tmp.aln@pattern@range@start),
                   width = tmp.aln@subject@range@width - (tmp.aln@pattern@range@width %% 3))
  # gisaid.0104.sg <- c(gisaid.0104.sg, tmp.sg)
  return(tmp.sg)
}

Biostrings::writeXStringSet(gisaid.0104.sg.2, "gisaid_sg_2.fasta")
rm(gisaid.0104.sg.2)

st <- 1179
en <- 589*3
gisaid.0104.sg.3 <- foreach(g=st:en, .inorder = FALSE, .errorhandling = "remove", .combine = c) %dopar% {
  
  library(Biostrings)
  tmp.aln <- Biostrings::pairwiseAlignment(refseq, gisaid.0104[g], type="local")
  # extract sequence in same reading frame as refseq seq:
  tmp.sg <- subseq(gisaid.0104[g], start = tmp.aln@subject@range@start + make_frame(tmp.aln@pattern@range@start),
                   width = tmp.aln@subject@range@width - (tmp.aln@pattern@range@width %% 3))
  # gisaid.0104.sg <- c(gisaid.0104.sg, tmp.sg)
  return(tmp.sg)
}

Biostrings::writeXStringSet(gisaid.0104.sg.3, "gisaid_sg_3.fasta")
rm(gisaid.0104.sg.3)


st <- 1768
en <- 589*4
gisaid.0104.sg.4 <- foreach(g=st:en, .inorder = FALSE, .errorhandling = "remove", .combine = c) %dopar% {
  
  library(Biostrings)
  tmp.aln <- Biostrings::pairwiseAlignment(refseq, gisaid.0104[g], type="local")
  # extract sequence in same reading frame as refseq seq:
  tmp.sg <- subseq(gisaid.0104[g], start = tmp.aln@subject@range@start + make_frame(tmp.aln@pattern@range@start),
                   width = tmp.aln@subject@range@width - (tmp.aln@pattern@range@width %% 3))
  # gisaid.0104.sg <- c(gisaid.0104.sg, tmp.sg)
  return(tmp.sg)
}

Biostrings::writeXStringSet(gisaid.0104.sg.4, "gisaid_sg_4.fasta")
rm(gisaid.0104.sg.4)


st <- 2357
en <- 589*5
gisaid.0104.sg.5 <- foreach(g=st:en, .inorder = FALSE, .errorhandling = "remove", .combine = c) %dopar% {
  
  library(Biostrings)
  tmp.aln <- Biostrings::pairwiseAlignment(refseq, gisaid.0104[g], type="local")
  # extract sequence in same reading frame as refseq seq:
  tmp.sg <- subseq(gisaid.0104[g], start = tmp.aln@subject@range@start + make_frame(tmp.aln@pattern@range@start),
                   width = tmp.aln@subject@range@width - (tmp.aln@pattern@range@width %% 3))
  # gisaid.0104.sg <- c(gisaid.0104.sg, tmp.sg)
  return(tmp.sg)
}

Biostrings::writeXStringSet(gisaid.0104.sg.5, "gisaid_sg_5.fasta")
rm(gisaid.0104.sg.5)


stopCluster(cl)

gisaid.0104.sg <- c(readDNAStringSet("gisaid_sg_1.fasta"),
                    readDNAStringSet("gisaid_sg_2.fasta"),
                    readDNAStringSet("gisaid_sg_3.fasta"),
                    readDNAStringSet("gisaid_sg_4.fasta"),
                    readDNAStringSet("gisaid_sg_5.fasta"))


# cl <- makeCluster(6)
# registerDoParallel(cl)
# 
# gisaid.0104.sg <- foreach(g=1:length(gisaid.0104), .combine = c, .inorder = FALSE, .errorhandling = "remove") %dopar% {
#   
#   library(Biostrings)
#   print(paste0(g, ": ", gisaid.0104@ranges@NAMES[g]))
#   tmp.aln <- Biostrings::pairwiseAlignment(refseq, gisaid.0104[g], type="local")
#   # extract sequence in same reading frame as refseq seq:
#   tmp.sg <- subseq(gisaid.0104[g], start = tmp.aln@subject@range@start + make_frame(tmp.aln@pattern@range@start),
#                    width = tmp.aln@subject@range@width - (tmp.aln@pattern@range@width %% 3))
#   # gisaid.0104.sg <- c(gisaid.0104.sg, tmp.sg)
#   return(tmp.sg)
# }
# 
# stopCluster(cl)

refseq.aa <- translate(refseq)
gisaid.0104.sg.aa <- translate(gisaid.0104.sg, if.fuzzy.codon = "solve")
writeXStringSet(gisaid.0104.sg.aa, "../gisaid_cov2020_sequences_01APR2020_human_host_low_cov_excl.S_prot_AA.fasta")

# get mismatches in table format
# mm.tab.gisaid.0104 <- data.frame(sequence=character(), position = numeric(), 
#                                  mutation = character(), consensus = character())
# 
# for(g in 1:length(gisaid.0104.sg.aa)){
#   print(paste0(g, ": ", gisaid.0104.sg.aa@ranges@NAMES[g]))
#   tmp.aln.aa <- Biostrings::pairwiseAlignment(refseq.aa, gisaid.0104.sg.aa[g], type="local")
#   
#   for(m in 1:length(unlist(tmp.aln.aa@subject@mismatch))){
#     tmp <- data.frame(sequence=gisaid.0104.sg.aa@ranges@NAMES[g], 
#                       position=unlist(tmp.aln.aa@pattern@mismatch)[m], 
#                       mutation=as.character(subseq(as.character(tmp.aln.aa@subject@unaligned), unlist(tmp.aln.aa@subject@mismatch)[m], width=1)),
#                       consensus=as.character(subseq(as.character(tmp.aln.aa@pattern@unaligned), unlist(tmp.aln.aa@pattern@mismatch)[m], width=1)))
#     mm.tab.gisaid.0104 <- rbind(mm.tab.gisaid.0104, tmp)
#   }
# }

seq <- Biostrings::readAAStringSet("20200318_EntrezNucleotide_SARSCoV2_AA_sequences.fasta")
seq.sg <- seq[grep("[protein=surface glycoprotein]", seq@ranges@NAMES, fixed = TRUE)]

seq.0904 <- Biostrings::readAAStringSet("20200409_EntrezNucleotide_SARSCoV2_AA_sequences.fasta")
seq.0904.sg <- seq.0904[grep("[protein=surface glycoprotein]", seq.0904@ranges@NAMES, fixed = TRUE)]


cl <- makeCluster(10)
registerDoParallel(cl)
# mm.tab.gisaid.0104 <- data.frame(sequence=character(), position = numeric(), mutation = character(), consensus = character())

get_mms <- function(g=g, ref=refseq.aa, seqs=seq.sg){
  tmp.aln.aa <- pairwiseAlignment(ref, seqs[g], type="local")
  
  mms <- data.frame(sequence=character(), position = numeric(), mutation = character(), consensus = character())
  
  # number of mismatches
  nm <- length(unlist(tmp.aln.aa@subject@mismatch))
  if(nm > 0){
    for(m in 1:nm){
      tmp <- data.frame(sequence=seqs@ranges@NAMES[g], 
                        position=unlist(tmp.aln.aa@pattern@mismatch)[m], 
                        mutation=as.character(subseq(as.character(tmp.aln.aa@subject@unaligned), unlist(tmp.aln.aa@subject@mismatch)[m], width=1)),
                        consensus=as.character(subseq(as.character(tmp.aln.aa@pattern@unaligned), unlist(tmp.aln.aa@pattern@mismatch)[m], width=1)))
      mms <- rbind(mms, tmp)
    }
  }
  # number of indels:
  ni <- sum(nindel(tmp.aln.aa)@insertion[,1],nindel(tmp.aln.aa)@deletion[,1])
  if( ni > 0){
    indel.ranges <- sort(c(unlist(insertion(tmp.aln.aa)), unlist(deletion(tmp.aln.aa))))
    for(m in 1:ni){
      # if deletion in sample
      # if(as.character(subseq(unaligned(subject(tmp.aln.aa)), start(tmp.aln.aa@subject@indel[[m]]), width=5)) == 
         # as.character(subseq(unaligned(pattern(tmp.aln.aa)), start(indel.ranges[m])+1, width=5))) {
      if(unlist(start(indel.ranges)) %in% unlist(start(insertion(tmp.aln.aa)))){
           mut <- "-"
           co <- subseq(as.character(aligned(pattern(tmp.aln.aa))), start(indel.ranges[m]),end(indel.ranges[m]))
      # if insertion in sample
      # } else if(as.character(subseq(unaligned(subject(tmp.aln.aa)), start(tmp.aln.aa@subject@indel[[m]])+1, width=5)) == 
      #           as.character(subseq(unaligned(pattern(tmp.aln.aa)), start(indel.ranges[m]), width=5))) {
      } else if (unlist(start(indel.ranges)) %in% unlist(start(deletion(tmp.aln.aa)))){
        mut <- subseq(as.character(aligned(subject(tmp.aln.aa))), start(indel.ranges[m]),end(indel.ranges[m]))
        co <- "-"
      }
           
      tmp <- data.frame(sequence=seqs@ranges@NAMES[g], 
                        position=start(indel.ranges)[m], 
                        mutation=mut,
                        consensus=co)
      mms <- rbind(mms, tmp)
    }
  }
  
  return(mms)
}

mm.tab.genbank.1803 <- 
  foreach(g=1:length(seq.sg), .combine = rbind) %dopar% {
    library(Biostrings)
    get_mms(g, refseq.aa, seq.sg)
  }  

mm.tab.genbank.0904 <- 
  foreach(g=1:length(seq.0904.sg), .combine = rbind) %dopar% {
    library(Biostrings)
    get_mms(g, refseq.aa, seq.0904.sg)
  } 
    
mm.tab.gisaid.0104 <-
  foreach(g=1:length(gisaid.0104.sg.aa), .combine = rbind) %dopar% {
    library(Biostrings)
    get_mms(g, refseq.aa, gisaid.0104.sg.aa)
    # tmp.aln.aa <- pairwiseAlignment(refseq.aa, gisaid.0104.sg.aa[g], type="local")
    # 
    # if(length(unlist(tmp.aln.aa@subject@mismatch)) > 0){
    #   for(m in 1:length(unlist(tmp.aln.aa@subject@mismatch))){
    #     tmp <- data.frame(sequence=gisaid.0104.sg.aa@ranges@NAMES[g], 
    #                       position=unlist(tmp.aln.aa@pattern@mismatch)[m], 
    #                       mutation=as.character(subseq(as.character(tmp.aln.aa@subject@unaligned), unlist(tmp.aln.aa@subject@mismatch)[m], width=1)),
    #                       consensus=as.character(subseq(as.character(tmp.aln.aa@pattern@unaligned), unlist(tmp.aln.aa@pattern@mismatch)[m], width=1)))
    #     return(tmp)
    #   }
    # }
}

stopCluster(cl)

gisaid.0104.sg.aa.qual <- data.frame(sample = gisaid.0104.sg@ranges@NAMES,
                                  n_fraction=str_count(gisaid.0104.sg.aa, "X")/gisaid.0104.sg.aa@ranges@width,
                                  width=gisaid.0104.sg.aa@ranges@width)
seq.0904.sg.qual <- data.frame(sample = (seq.0904.sg@ranges@NAMES %>% str_split(., "_prot_", simplify = T))[,1] %>% gsub("lcl\\|", "", .),
                               n_fraction=str_count(seq.0904.sg, "X")/seq.0904.sg@ranges@width,
                               width=seq.0904.sg@ranges@width)

we <- rbind(gisaid.0104.sg.aa.qual, seq.0904.sg.qual) %>% filter(width < 100)


mm.tab.genbank.0904$sequence <- gsub("lcl\\|", "", str_split(mm.tab.genbank.0904$sequence, "_prot_", simplify = T)[,1])
mm.tab.genbank.0904$location <- NA
summary.file <- "20200409_EntrezNucleotide_SARSCoV2_summary.txt"
for(id in unique(mm.tab.genbank.0904$sequence)){
  s <- grep(paste0("^", id), readLines(summary.file))
  ann <- gsub("^[0-9]+\\. ", "", read.delim(summary.file, header = FALSE, sep="%", skip = s-3, nrows = 1)[1,1])
  mm.tab.genbank.0904$location[mm.tab.genbank.0904$sequence == id] <- ann
}
mm.tab.genbank.0904$location <- str_split(mm.tab.genbank.0904$location, " isolate ", simplify = T)[,2]

mm.tab.genbank.0904$location[grep("ESP", mm.tab.genbank.0904$location)] <- "Spain"
mm.tab.genbank.0904$location[grep("USA", mm.tab.genbank.0904$location)] <- "USA"
mm.tab.genbank.0904$location[grep("Wuhan", mm.tab.genbank.0904$location)] <- "Wuhan"
mm.tab.genbank.0904$location[grep("FIN", mm.tab.genbank.0904$location)] <- "Finland"
mm.tab.genbank.0904$location[grep("IND", mm.tab.genbank.0904$location)] <- "India"
mm.tab.genbank.0904$location[grep("PER", mm.tab.genbank.0904$location)] <- "Peru"
mm.tab.genbank.0904$location[grep("ISR", mm.tab.genbank.0904$location)] <- "Israel"
mm.tab.genbank.0904$location[grep("Australia", mm.tab.genbank.0904$location)] <- "Australia"
mm.tab.genbank.0904$location[grep("SWE", mm.tab.genbank.0904$location)] <- "Sweden"
mm.tab.genbank.0904$location[grep("CHN", mm.tab.genbank.0904$location)] <- "China"
mm.tab.genbank.0904$location[grep("COL", mm.tab.genbank.0904$location)] <- "Colombia"
mm.tab.genbank.0904$location[grep("SNU", mm.tab.genbank.0904$location)] <- "South Korea"



mm.tab.gisaid.0104$location <- unlist(lapply(str_split(mm.tab.gisaid.0104$sequence, "/"), "[[", 2))
mm.tab.gisaid.0104$location <- gsub("NetherlandsL", "Netherlands", mm.tab.gisaid.0104$location)


# require(stringr)
# ref <- unlist(str_split(refseq.aa, ""))
# mm.tab.gisaid.0104 <- data.frame(sequence=character(), position = numeric(), mutation = character(), consensus = character())
# for (i in 1:nrow(gisaid.0104.sg.aa)-1) {
#   s <- unlist(str_split(gisaid.0104.sg.aa[i], ""))
#   mm <- which(s != ref)
#   if(length(mm)<100){
#     for (m in mm){
#       tmp <- data.frame(sequence=gisaid.0104.sg.aa@ranges@NAMES[i], position=m, mutation=s[m], consensus=ref[m])
#       mm.tab.gisaid.0104 <- rbind(mm.tab.gisaid.0104, tmp)
#     }
#   } else {
#     for (mi in 1:(length(mm)-1)){
#       if(s[mm[mi]] == s[mm[mi+1]]) {next}
#       tmp <- data.frame(sequence=gisaid.0104.sg.aa@ranges@NAMES[i], position=mm[mi], mutation=s[mm[mi]], consensus=ref[mm[mi]])
#       mm.tab.gisaid.0104 <- rbind(mm.tab.gisaid.0104, tmp)
#     }
#   }
# }

write.csv2(mm.tab.gisaid.0104, "gisaid_cov2020_sequences_01APR2020_human_host_low_cov_excl.S_prot_AA.mm.csv", quote = FALSE, row.names = FALSE)
write.csv2(mm.tab.genbank.0904, "20200409_EntrezNucleotide_SARSCoV2_AA_sequences.mm.csv", quote = FALSE, row.names = FALSE)

all.qual <- rbind(seq.0904.sg.qual, gisaid.0104.sg.aa.qual)
mm.tab.all <- rbind(mm.tab.genbank.0904, mm.tab.gisaid.0104) %>% 
  filter(!(sequence %in% we)) %>% # exclude samples with width < 100 ==> no sample needs to be excluded
  mutate(qual=all.qual$n_fraction[match(.$sequence, all.qual$sample)])

mm.tab.all$location[mm.tab.all$location == "USA"] <- "United States"
mm.tab.all$location[mm.tab.all$location == "Congo"] <- "Congo [DRC]"

# longitude/latitude from https://developers.google.com/public-data/docs/canonical/countries_csv (16APR2020)
# cities and UK from www.latlong.net
geo.info <- read.table("geo_coordinates.tsv", stringsAsFactors = FALSE, row.names = 4, sep="\t", dec = ".", quote = "%", colClasses = c("character", "numeric", "numeric"), header=TRUE)
geo.info <- geo.info[(!is.na(geo.info$latitude) & !is.na(geo.info$longitude)),c("latitude", "longitude")]

geo.info <- rbind(geo.info,
                  read.table("geo_coordinates_UK.tsv", stringsAsFactors = FALSE, row.names = 3, sep="\t", dec = ".", quote = "%", header=TRUE),
                  read.table("geo_coordinates_cities.tsv", stringsAsFactors = FALSE, row.names = 3, sep="\t", dec = ".", quote = "%", header=TRUE))

geo.pca <- prcomp(geo.info)
geo.order <- sort(geo.pca$x[,1])


# number of sequences with mutations
mm.tab.gisaid.0104 %>% filter(mutation != "X") %>% distinct(sequence) %>% nrow(.)
# [1] 1233
mm.tab.genbank.1803 %>% filter(mutation != "X") %>% distinct(sequence) %>% nrow(.)
# [1] 14
mm.tab.genbank.0904 %>% filter(mutation != "X") %>% distinct(sequence) %>% nrow(.)
# [1] 186

mm.tab.all %>% filter(mutation != "X") %>% distinct(sequence) %>% nrow(.)
# [1] 1635

# number of mutations in one sequence
mm.tab.gisaid.0104 %>% filter(mutation != "X") %>% group_by(sequence) %>% summarise(obs = n()) %>% arrange(desc(obs)) %>%
  ggplot(aes(x=obs)) + geom_histogram(binwidth = 1) + theme_bw() + coord_trans(y="log1p") + 
  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-1.5) + ylim(c(0,2500)) + xlab("# mutations in one sequence") + 
  scale_x_continuous(breaks=seq(0, 8, 1))
hist(mm.tab.gisaid.0104 %>% filter(mutation != "X") %>% group_by(sequence) %>% summarise(obs = n()) %>% arrange(desc(obs)) %>% pull(obs), xlab = "# mutations in one sequence", main="Histogram")

mm.tab.genbank.1803 %>% filter(mutation != "X") %>% group_by(sequence) %>% summarise(obs = n()) %>% arrange(desc(obs))
# hist(mm.tab.genbank.1803 %>% filter(mutation != "X") %>% group_by(sequence) %>% summarise(obs = n()) %>% arrange(desc(obs)) %>% pull(obs), xlab = "# mutations in one sequence", freq = F, main="Histogram")
mm.tab.genbank.0904 %>% filter(mutation != "X") %>% group_by(sequence) %>% summarise(obs = n()) %>% arrange(desc(obs))



# mutations in RBD (pos 319-529)
mm.tab.gisaid.0104 %>% filter(mutation != "X",
                         (position >= 319 & position <= 529)) %>% 
                        group_by(position, mutation, consensus) %>% 
                        summarise(obs = n()) %>% arrange(desc(obs)) #%>% 
                        write.csv2("gisaid_cov2020_sequences_01APR2020_human_host_low_cov_excl.S_prot_AA.mm.RBD.csv", quote = FALSE, row.names = FALSE)

mm.tab.genbank.0904 %>% filter(mutation != "X",
                              (position >= 319 & position <= 529)) %>% 
  group_by(position, mutation, consensus) %>% 
  summarise(obs = n()) %>% arrange(desc(obs)) #%>% 
  write.csv2("20200409_EntrezNucleotide_SARSCoV2_AA_sequences.mm.RBD.csv", quote = FALSE, row.names = FALSE)
  
mm.tab.all %>% filter(mutation != "X", (position >= 319 & position <= 529)) %>% 
  group_by(position, mutation, consensus) %>% 
  summarise(obs = n()) %>% arrange(desc(obs))
# 25 variants in RBD

# recurrent mutations 
mm.tab.gisaid.0104 %>% filter(mutation != "X") %>% 
  group_by(position, mutation, consensus) %>% summarise(obs = n()) %>% arrange(desc(obs), position) %>% filter(obs > 1) %>%
  write.csv2("gisaid_cov2020_sequences_01APR2020_human_host_low_cov_excl.S_prot_AA.mm.recurrent.csv", quote = FALSE, row.names = FALSE)

mm.tab.genbank.1803 %>% filter(mutation != "X") %>% 
  group_by(position, mutation, consensus) %>% summarise(obs = n()) %>% arrange(desc(obs)) %>% filter(obs > 1) 

mm.tab.genbank.0904 %>% filter(mutation != "X") %>% 
  group_by(position, mutation, consensus) %>% summarise(obs = n()) %>% arrange(desc(obs)) %>% filter(obs > 1) 

mm.tab.all %>% filter(mutation != "X") %>% 
  group_by(position, mutation, consensus) %>% summarise(obs = n()) %>% arrange(desc(obs), position) %>% filter(obs > 1) 
# 70 recurrent events

# Wan et al., 2020, J. Virol., doi:10.1128/JVI.00127-20
critical <- c(455, 486, 493, 494, 501)
contact <- c(415, 439, 449, 453, 455, 486, 487, 489, 493, 498, 500, 501, 502, 505)

mm.tab.gisaid.0104 %>% filter(mutation != "X",
                              position %in% contact) 
#  hCoV-19/Malaysia/186197/2020|EPI_ISL_417919|2020-03-14      449        N         Y Malaysia
  
mm.tab.genbank.1803 %>% filter(mutation != "X",
                              position %in% critical)

mm.tab.genbank.0904 %>% filter(mutation != "X",
                               position %in% contact)


mm.tab.all %>% filter(mutation != "X") %>% 
  ggplot(aes(x=mutation, y=position)) + 
  geom_rect(xmin=-Inf, xmax=Inf, ymin=319, ymax=529, fill="red", alpha=0.8) +
  geom_point(shape=15) 

library(ggrepel)  
mm.tab.all %>% filter(mutation != "X") %>% 
  group_by(position, mutation, consensus) %>% summarise(obs = n()) %>%
  ggplot(aes(x=log10(obs), y=position)) + 
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=319, ymax=529, alpha=0.5, fill="red") +
    geom_point(shape=15, alpha=0.5) + theme_bw() + 
  xlab("# samples with mutation (log10)") + ylab("AA position S protein (red: RBD)") + 
  geom_text_repel( 
    data= mm.tab.all %>% 
      group_by(position, mutation, consensus) %>% summarise(obs = n()) %>% filter(mutation != "X",
                                                                                  obs > 5) , # Filter data first
    aes(label=str_c(consensus, position, mutation, sep="")),
    size=2.25, seed = 1234, min.segment.length = 0.1, 
    direction = "y", nudge_x = 0.05, hjust = 0,
    fontface = "bold"
  )

# sequences without mm
nmm <- (rbind(gisaid.0104.sg.aa.qual, seq.0904.sg.qual) %>% filter(width >= 100) %>% 
          pull(sample) %>% str_split(., "_prot_", simplify = T))[,1] %>% gsub("lcl\\|", "", .) 
nmm <- nmm[!(nmm %in% mm.tab.all$sequence)]

mm.tab.all %>% filter(mutation != "X") %>% group_by(sequence) %>% summarise(obs = n()) %>% arrange(desc(obs)) %>%
  add_row(sequence = nmm, obs=0) %>%
  ggplot(aes(x=obs)) + geom_histogram(binwidth = 1, color="black", fill="lightgrey") + theme_bw() + coord_trans(y="log1p") + 
  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-1.5, size=3) + ylim(c(0,2500)) + xlab("# mutations in one sequence") + 
  scale_x_continuous(breaks=seq(0, 8, 1)) 


library(RColorBrewer)
nb.cols <- length(unique(mm.tab.all$location)) # Define the number of colors needed
mycolors <- colorRampPalette(brewer.pal(8, "Spectral"))(nb.cols)
myorder <- mm.tab.all %>% filter(mutation != "X") %>% distinct(location) %>% select(location) %>%
  mutate(sorting = geo.order[match(location, names(geo.order))]) %>% arrange(sorting) %>% pull(location)

mm.tab.all %>% filter(mutation != "X") %>% 
  ggplot(aes(x=sequence, y = position, color = location)) + 
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=319, ymax=529, alpha=0.5, fill="grey") +
  geom_point(shape=15) + theme_bw() + theme(axis.text.x = element_blank()) +
    scale_x_discrete(limits=mm.tab.all %>% filter(mutation != "X") %>% mutate(sorting = geo.order[match(location, names(geo.order))]) %>%
                     arrange(position, mutation, sorting) %>% distinct(sequence) %>% pull(sequence),
                   breaks = NULL) + 
  scale_colour_manual(values = mycolors, limits = myorder)

ggsave("plot_all_variants_all_samples.png")

# RBD only
nb.cols <-  mm.tab.all %>% filter(mutation != "X",
                                  (position > 319 & position < 529)) %>% distinct(location) %>% nrow(.) # Define the number of colors needed
mycolors2 <- colorRampPalette(brewer.pal(8, "Spectral"))(nb.cols)
myorder2 <- mm.tab.all %>% filter(mutation != "X",
                                  (position > 319 & position < 529)) %>% distinct(location) %>% select(location) %>%
  mutate(sorting = geo.order[match(location, names(geo.order))]) %>% arrange(sorting) %>% pull(location)

mm.tab.all %>% filter(mutation != "X",
                      (position > 319 & position < 529)) %>%
  ggplot(aes(x=sequence, y = position, color = location)) + 
  # annotate("rect", xmin=-Inf, xmax=Inf, ymin=319, ymax=529, alpha=0.5, fill="grey") +
  geom_point(shape=15, size=4, show.legend = FALSE) + theme_bw() + theme(axis.text.x = element_blank()) +
  scale_x_discrete(limits=mm.tab.all %>% filter(mutation != "X", (position > 319 & position < 529)) %>% mutate(sorting = geo.order[match(location, names(geo.order))]) %>%
                     arrange(position, mutation, sorting) %>% distinct(sequence) %>% pull(sequence),
                   breaks = NULL) + 
  scale_colour_manual(values = mycolors, limits = myorder)
ggsave("plot_all_variants_all_samples.png")
                   