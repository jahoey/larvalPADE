setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files")

library(ade4)
library(adegenet)
library(pegas)
library(hierfstat)
library(mmod)
library(poppr)

full <- read.structure("structure_input_Mar14_2017_528_1904.str",
                       n.ind = 528, n.loc = 1904, col.lab = 1, col.pop = 2, row.marknames = 1,
                       onerowperind = FALSE)

full.dataframe <- as.data.frame(full@tab)

#### Dealing with loci with multiple alleles ####
multiallelic <- names(which(full@loc.n.all > 2))

#### Calculate FST statistics ####
# Read in STRUCTURE file as a text file containing larvae & adults but multiallelic loci
all528_multiallelic <- read.table("structure_input_Mar14_2017_528_1904.txt", sep="\t", header = TRUE)
dim(all528_multiallelic) # 1056 x 1906

# Get rid of multiallelic loci
# Multiallelic loci list is in multiallelic
multiallelic <- names(which(full@loc.n.all > 2)) # 22 multiallelic loci
all528_biallelic <- all528_multiallelic[,!(colnames(all528_multiallelic) %in% multiallelic)]
dim(all528_biallelic) #1056 x 1884

# Split data into odd and even row dataframes
even_indexes<-seq(2,1056,2)
odd_indexes<-seq(1,1055,2)

odds <- data.frame(all528_biallelic[odd_indexes,]) # 528 x 1884
odds2 <- odds[,-c(1:2)] # 528 x 1882
evens <- data.frame(all528_biallelic[even_indexes,]) # 528 x 1884
evens2 <- evens[,-c(1:2)] # 528 x 1882

# Now paste each value in one dataframe with its corresponding in the other
s <- 1:length(colnames(evens2))
combo <- data.frame(matrix(nrow = 528, ncol = 1882))
for (i in s){
  combo[,i] <-paste(odds2[,i], evens2[,i], sep = '')
}

dim(combo) # 528 x 1882

combo[] <- lapply(combo, function(x) as.numeric(as.character(x)))# Convert to numeric, gives warning because replaces character 'NANA' with NA

pop.names <- as.numeric(as.character(evens$Pop))

# Combine the population numbers with the allele data
combo2 <- cbind(pop.names, combo)
dim(combo2) # 528 x 1883

pairwise.WCfst(combo2,diploid=TRUE) 
genet.dist(combo2, method = 'WC84')
