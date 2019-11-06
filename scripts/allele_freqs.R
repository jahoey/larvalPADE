#### Calculating MAF for 1904 loci (1st SNP on a contig) used for dispersal analysis vs. subsequent SNPs on a contig

#### This section is for the 1904 loci used in the dispersal paper ####
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files")

library(ade4)
library(adegenet)
library("hierfstat")
library(pegas)


# Reading in the allele frequency data for adults and larvae
full <- read.structure("structure_input_Mar14_2017_528_1904.str",
                       n.ind = 528, n.loc = 1904, col.lab = 1, col.pop = 2, row.marknames = 1,
                       onerowperind = FALSE)

full.dataframe <- as.data.frame(full@tab)

# Make vector of adult names
adults <- rownames(full.dataframe)[205:439]

# Remove adults from dataset
larvs <- full.dataframe[!(rownames(full.dataframe) %in% adults),]
larvs2 <- as.genind(larvs)

# Calculate allele frequencies
# Way 1 - Includes multiallelic loci
colSums(larvs,na.rm=TRUE)/(2*colSums(!is.na(larvs)))

# Way 2 - Only biallelic loci
# Multiallelic loci list is in multiallelic
multiallelic <- names(which(larvs2@loc.n.all > 2))
multiallelic <- paste(multiallelic, ".", sep = "")

multiallelic_data <- larvs2@tab[, grep(paste(multiallelic, "{3,4}", sep = "", collapse = "|"), colnames(larvs2@tab))] # mostly works, but still includes SNP_145*
multiallelic_data <- multiallelic_data[,-c(59:78)] # gets rid of additional SNP_145* that are actually biallelic but got picked up in the grep above
multiallelic_counts <- colSums(multiallelic_data, na.rm = TRUE)

multiallelic_names <- names(multiallelic_counts) #67
biallelic_only <- larvs[,!colnames(larvs) %in% multiallelic_names] # 293x3764 (before: 3831-67 = 3764)

# Split data into odd and even column dataframes
even_indexes<-seq(2,3764,2)
odd_indexes<-seq(1,3763,2)

# Biallelic
odds <- data.frame(biallelic_only[,odd_indexes]) # 293 x 1882
evens <- data.frame(biallelic_only[,even_indexes]) # 293 x 1882

odds.allelefreq <- colSums(odds,na.rm=TRUE)/(2*colSums(!is.na(odds)))
evens.allelefreq <- colSums(evens,na.rm=TRUE)/(2*colSums(!is.na(evens)))

# Pairwise minimums
maf <- pmin(odds.allelefreq, evens.allelefreq)

# Plot the MAF
hist(maf, breaks = 20, main = 'MAF: First SNPs on a contig')

#### This section is for SNPs that weren't the first one on the contig, and therefore got cut out ####
# Reading in the allele frequency data for adults and larvae
full2 <- read.structure("structure_input_Nov5_2019_528_1923.str",
                       n.ind = 528, n.loc = 1923, col.lab = 1, col.pop = 2, row.marknames = 1,
                       onerowperind = FALSE)

full2.dataframe <- as.data.frame(full2@tab)

# Make vector of adult names
adults <- rownames(full2.dataframe)[205:439]

# Remove adults from dataset
larvs <- full2.dataframe[!(rownames(full2.dataframe) %in% adults),]
larvs2 <- as.genind(larvs)

# Multiallelic loci list is in multiallelic
multiallelic <- names(which(larvs2@loc.n.all > 2))
multiallelic <- paste(multiallelic, ".", sep = "")

multiallelic_data <- larvs2@tab[, grep(paste(multiallelic, "{3,4}", sep = "", collapse = "|"), colnames(larvs2@tab))] # mostly works, but still includes SNP_110*, SNP_167*
multiallelic_data <- multiallelic_data[,-c(55:74, 84:103)] # gets rid of additional SNP_110*, SNP_167* that are actually biallelic but got picked up in the grep above
multiallelic_counts <- colSums(multiallelic_data, na.rm = TRUE)

multiallelic_names <- names(multiallelic_counts) #72
biallelic_only <- larvs[,!colnames(larvs) %in% multiallelic_names] # 293x3798 (before: 3870-72 = 3798)

# Split data into odd and even column dataframes
even_indexes<-seq(2,3798,2)
odd_indexes<-seq(1,3797,2)

# Biallelic
odds <- data.frame(biallelic_only[,odd_indexes]) # 293 x 1899
evens <- data.frame(biallelic_only[,even_indexes]) # 293 x 1899

odds.allelefreq <- colSums(odds,na.rm=TRUE)/(2*colSums(!is.na(odds)))
evens.allelefreq <- colSums(evens,na.rm=TRUE)/(2*colSums(!is.na(evens)))

# Pairwise minimums
maf <- pmin(odds.allelefreq, evens.allelefreq)

# Plot the MAF
hist(maf, breaks = 20, main = 'MAF: Excluding first SNPs on contig')

