
library(Biostrings)
library(tidyverse)
library(muscle)
library(DECIPHER)
library(ape)
library(viridis)
library(cluster)
library(seqinr)
library(adegenet)
library(pegas)
library(apex)
library(mmod)
library(poppr)

Wolf <- read_tsv("http://www.boldsystems.org/index.php/API_Public/combined?taxon=Canis&format=tsv")
write_tsv(Wolf, "Wolf_BOLD_data.tsv")
Canis <- read_tsv("Wolf_BOLD_data.tsv")

#Reduced Canis 
Canis_rd <- Canis %>%
  select(processid, genus_name, species_name, subspecies_name, bin_uri, country, markercode, nucleotides, lon, lat, recordID) %>%
  filter(species_name == "Canis lupus" | species_name == "Canis familiaris") %>%
  filter(markercode == "COI-5P")

#Check
#count gets masked so do this
dplyr::count(Canis_rd, species_name)
dplyr::count(Canis_rd, subspecies_name)
dplyr::count(Canis_rd, markercode)
#Working on the null Canis familiaris is just a sub species of Canis and should be treated as such
#Relabel any Canis familiaris as such

for (n in 1:nrow(Canis_rd)){
  if (Canis_rd[n,"species_name"] == "Canis familiaris") {
    Canis_rd[n,"species_name"] <- "Canis lupus"
    Canis_rd[n,"subspecies_name"] <- "Canis lupus familiaris"
  } else if (is.na(Canis_rd[n,"subspecies_name"])){
    Canis_rd[n,"subspecies_name"] <- "Canis lupus"
  }
}

#Check
dplyr::count(Canis_rd, species_name)
dplyr::count(Canis_rd, subspecies_name)
#Increased by 4 as expected

#Quality checks and filtering
Canis_filtered <- Canis_rd %>%
  mutate(nucleotides2 = str_remove_all(nucleotides, "^N+|N+$|-")) %>%
  filter(!is.na(nucleotides2)) %>%
  mutate(species_name = str_replace(species_name, "Canis", "C.")) %>%
  mutate(subspecies_name = str_replace(subspecies_name, "Canis", "C.")) %>%
  mutate(ID = paste(subspecies_name, recordID, sep="_"))


summary(nchar(Canis_filtered$nucleotides2))

ggplot(Canis_filtered,aes(x=nchar(nucleotides2), y = subspecies_name, colour = subspecies_name), xmin = 0, xmax = max(nchar(nucleotides2))+100) +
  scale_x_continuous(breaks =seq(0, max(nchar(Canis_filtered$nucleotides2))+100 , by = 200)) +
  geom_jitter(show.legend = F) +
  geom_boxplot(outlier.shape = NA, colour="black", show.legend = F, fill = NA) +
  xlab("Sequence Length")

#Going to filter over 600 - 1600 
Canis_filtered <- Canis_filtered %>%
  filter(nchar(nucleotides2) >= 600) %>%
  filter(nchar(nucleotides2) <= 1600)

summary(nchar(Canis_filtered$nucleotides2))
dplyr::count(Canis_filtered, subspecies_name)

Canis_filtered <- Canis_filtered  %>%
  mutate(Undefined = str_count(nucleotides2, "N")/nchar(nucleotides2))

summary(Canis_filtered$Undefined)
sort(Canis_filtered$Undefined, decreasing = T)[1:20]

Canis_filtered <- Canis_filtered  %>%
  filter(Undefined < 0.01)

dplyr::count(Canis_filtered, subspecies_name)

#rm(Canis_rd)
# Align time
Canis_filtered <- as.data.frame(Canis_filtered)

class(Canis_filtered$nucleotides2)
Canis_filtered$nucleotides2 <- DNAStringSet(Canis_filtered$nucleotides2)

names(Canis_filtered$nucleotides2) <- Canis_filtered$ID

##Analysis with muscle. Dataset is small so will not adjust maxiters. Will run twice, with out setting a gap penalty and once with an intense gap penalty to see how that changes the alignment and the results
Canis_alignment <- DNAStringSet(muscle::muscle(Canis_filtered$nucleotides2, maxhours = 1), use.names = TRUE)
writeXStringSet(Canis_alignment, file = "CanisAlignment.fasta")
BrowseSeqs(Canis_alignment)


Canis_Bin <- as.DNAbin(Canis_alignment)
Canis_distanceMatrix <- dist.dna(Canis_Bin, model = "k80", as.matrix = TRUE, pairwise.deletion = TRUE)
Canis_clusters <- IdClusters(Canis_distanceMatrix, method = "NJ", cutoff = 0.02, showPlot = TRUE, type = "both", verbose = TRUE)

count(Canis_clusters[[1]], cluster)
Check <- filter(Canis_clusters[[1]], cluster == 1)
#Good might be a reverse 

# Lets start analysis

#Will help later 
Canis_Multi <- read.multiFASTA("CanisAlignment.fasta")

Canis_genind <- multidna2genind(Canis_Multi)

#Set pop
strata(Canis_genind) <- data.frame("Pop" = Canis_filtered$subspecies_name)
head(strata(Canis_genind))
setPop(Canis_genind) <- ~Pop
#Set of outputs to analysis 
diff_stats(Canis_genind, phi_st = T)

bootstrap <- chao_bootstrap(Canis_genind, nreps = 100)
#Could change statistic used here, probs should and justify 
summarise_bootstrap(bootstrap, Gst_Nei)

#amova
Canis_DistPair <- dist.multidna(Canis_Multi, pool = T)
#performing AMOVA
amova(Canis_DistPair ~ Pop, data = strata(Canis_genind), nperm = 100)

#PCA
X <- scaleGen(Canis_genind, NA.method="mean")
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))
s.label(pca1$li)

col <- funky(15)
CanFac <- as.factor(Canis_filtered$subspecies_name)
s.class(pca1$li, CanFac,xax=1,yax=3, col=transp(col,.6), axesell=FALSE,
        cstar=0, cpoint=3, grid=FALSE)

class(pop(Canis_genind))
