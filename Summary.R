#Lecture 19
#Requirred:
#install.packages("apex")
library("apex")
#install.packages("adegenet")
library("adegenet")
#install.packages("pegas")
library("pegas")
#install.packages("mmod")
library("mmod")
#install.packages("poppr")
library("poppr")

#Reading in the data from fasta 
#See tutorial for how the data was prepared to get to that state
#Function used:
beeData <- read.multiFASTA(c("Ebom_mt.fas", "Ebom_CAD.fas"))

#Editing the names, this i only have to do if i do multi genes which i am not!
getLocusNames(beeData)
(setLocusNames(beeData) <- gsub(".fas", "", getLocusNames(beeData)))

#Now extraxcting the SNP information for further anylsis
#I will set mlist to F because only one loci 
beeData.gid <- multidna2genind(beeData, mlst = TRUE)

#Check
class(beeData.gid)
summary(beeData.gid)
#Can take advantage of the type option and tab to maybe create a matrix that can be used for PCA
beeData.gid$tab[1:10,1:10]

#Now creating population structures
my_strata <- data.frame(regions = rep(c("West", "East"), each = 20), 
                        populations = rep(c("CA", "Ch", "Am", "AF"), each = 10))
#Link it
strata(beeData.gid) <- my_strata
beeData.gid@strata
#Offically set pop
setPop(beeData.gid) <- ~populations

#Visuialize results!
diff_stats(beeData.gid)
#If we only want PhiST (a variant for sequence data, which considers distances between sequences), we can use the following. Note we also could have added this as an argument to diff_stats().
Phi_st_Meirmans(beeData.gid)
#Could be good for part two questions

#Bootstrap to see if significantly different
bs <- chao_bootstrap(beeData.gid, nreps = 100)
summarise_bootstrap(bs, Gst_Nei)

#"computes pairwise genetic distances between individuals" (from documentation)
beeData_dist <- dist.multidna(beeData, pool = TRUE)
#performing AMOVA
amova(beeData_dist ~ populations, data = strata(beeData.gid), nperm = 100)

#PCA 
#using a genid object
X <- scaleGen(microbov, NA.method="mean")
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
#See tutorial for more ways to play with the image 
#An introduction to adegenet 2.0.0 play aroudn see what looks nice

