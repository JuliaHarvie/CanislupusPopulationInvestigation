#Identify all members of canis 

#Use a dendrogram to determine the closest to wolf and dog
#(Perhaps also try PCA here and compare)
#Perform F test and mavona on the 3 selected species

library(rentrez)
search <- entrez_search(db = "nuccore", term = "canis lupus[ORGN] cytb", retmax = 100)
search$count
search <- entrez_search(db = "nuccore", term = "canis familiaris[ORGN] cytb", retmax = 100)
search$count
search <- entrez_search(db = "nuccore", term = "canis latrans[ORGN] cytb", retmax = 100)
search 
search$count
search <- entrez_search(db = "nuccore", term = "canis simensis[ORGN] cytb", retmax = 100)
search 
search <- entrez_search(db = "nuccore", term = "canis mesomelas[ORGN] cytb", retmax = 100)
search
search <- entrez_search(db = "nuccore", term = "canis lupaster[ORGN] cytb", retmax = 100)
search 
search <- entrez_search(db = "nuccore", term = "canis aureus[ORGN] cytb", retmax = 100)
search 
search <- entrez_search(db = "nuccore", term = "canis adustus[ORGN] cytb", retmax = 100)
search 

search <- entrez_search(db = "nuccore", term = "Canis latrans[ORGN] GHR", retmax = 100)
search 
search <- entrez_search(db = "nuccore", term = "Canis lupus[ORGN] GHR", retmax = 100)
search 
search <- entrez_search(db = "nuccore", term = "Canis familiaris[ORGN] GHR", retmax = 100)
search 
search <- entrez_search(db = "nuccore", term = "canis simensis[ORGN] GHR", retmax = 100)
search



