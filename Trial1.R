library(rentrez)
library(Biostrings)
library(tidyverse)
library(muscle)
library(DECIPHER)
library(ape)
#####C. lupus ctyb####
searchterm = "(canis lupus[ORGN] AND CytB[Gene]) NOT (canis lupus familiaris[ORGN]) NOT (genome[TITL])"
search <- entrez_search(db = "nuccore", term = searchterm)
max <- search$count
search <- entrez_search(db = "nuccore", term = searchterm, retmax = max, use_history = T)
fetch <- entrez_fetch(db = "nuccore", web_history=search$web_history, rettype = "fasta")
write(fetch, "Clupus_cytb.fasta", sep = "\n")
Clupus_cytb <- readDNAStringSet("Clupus_cytb.fasta")

#####C. familiaris ctyb ####
searchterm = "(canis familiaris[ORGN] AND CytB[Gene]) NOT (genome[TITL])"
search <- entrez_search(db = "nuccore", term = searchterm)
max <- search$count
search <- entrez_search(db = "nuccore", term = searchterm, retmax = max, use_history = T)
fetch <- entrez_fetch(db = "nuccore", web_history=search$web_history, rettype = "fasta")
write(fetch, "Cfamiliaris_cytb.fasta", sep = "\n")
Cfamiliaris_cytb <- readDNAStringSet("Cfamiliaris_cytb.fasta")

#####C. latrans ctyb ####
searchterm = "(canis latrans[ORGN] AND CytB[Gene]) NOT (genome[TITL])"
search <- entrez_search(db = "nuccore", term = searchterm)
max <- search$count
search <- entrez_search(db = "nuccore", term = searchterm, retmax = max, use_history = T)
fetch <- entrez_fetch(db = "nuccore", web_history=search$web_history, rettype = "fasta")
write(fetch, "Clatrans_cytb.fasta", sep = "\n")
Clatrans_cytb <- readDNAStringSet("Clatrans_cytb.fasta")

######C. simensis ctyb ####
####C. mesomelas ctyb ####
####C. lupaster ctyb ####
####C. aureus ctyb ####
####C. adustus ctyb ####

####C. lupus GHR ####
####C. familiaris GHR ####
####C. latrans GHR ####
####C. simensis GHR ####
####C. mesomelas GHR ####
####C. lupaster GHR ####
####C. aureus GHR ####
####C. adustus GHR ####




CFdf <- data.frame(Title = names(Cfamiliaris_cytb), Sequence = paste(Cfamiliaris_cytb))
CFdf$Species_Name <- paste(word(CFdf$Title, 2L, 4L), word(CFdf$Title, 1L))
CFdf <- CFdf %>%
  mutate(Sequence2 = str_remove_all(Sequence, "^N+|N+$|-")) %>%
  mutate(Undefined = str_count(Sequence2, "N")/nchar(Sequence2))
summary(CFdf)

summary(nchar(CFdf$Sequence2))
CFdf <- CFdf %>%
  filter(nchar(Sequence2) < 2000)
hist(nchar(CFdf$Sequence2))
length(nchar(CFdf$Sequence2))

CLdf <- data.frame(Title = names(Clupus_cytb), Sequence = paste(Clupus_cytb))
CLdf$Species_Name <- paste(word(CLdf$Title, 2L, 3L), word(CLdf$Title, 1L))
CLdf <- CLdf %>%
  mutate(Sequence2 = str_remove_all(Sequence, "^N+|N+$|-")) %>%
  mutate(Undefined = str_count(Sequence2, "N")/nchar(Sequence2))
summary(nchar(CLdf$Sequence2))
CLdf <- CLdf %>%
  filter(nchar(Sequence2) < 800)
  filter(nchar(Sequence2) > 600)


CLdf$Sequence2 <- DNAStringSet(CLdf$Sequence2)
names(CLdf$Sequence2) <- CLdf$Species_Name 
CLdf_alignment <- DNAStringSet(muscle::muscle(CLdf$Sequence2), use.names = TRUE)
BrowseSeqs(CLdf_alignment)
writeXStringSet(CLdf_alignment, file = "CLAlign.fasta")

CFdf$Sequence2 <- DNAStringSet(CFdf$Sequence2)
names(CFdf$Sequence2) <- CFdf$Species_Name 
CFdf_alignment <- DNAStringSet(muscle::muscle(CFdf$Sequence2), use.names = TRUE)
BrowseSeqs(CFdf_alignment)
writeXStringSet(CFdf_alignment, file = "CFAlign.fasta")


CanData <- apex::read.multiFASTA(c("CLAlign.fasta","CFAlign.fasta"))
plot(CanData, cex = 2)



