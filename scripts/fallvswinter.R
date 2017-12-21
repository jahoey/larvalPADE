setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files")

library(ade4)
library(adegenet)
library(pegas)
library(hierfstat)
library(mmod)
library(poppr)
library(diveRsity)

full <- read.structure("structure_input_Mar14_2017_528_1904.str",
                       n.ind = 528, n.loc = 1904, col.lab = 1, col.pop = 2, row.marknames = 1,
                       onerowperind = FALSE)

full.dataframe <- as.data.frame(full@tab)

#### Dealing with loci with multiple alleles ####
multiallelic <- names(which(full@loc.n.all > 2))
multiallelic <- paste(multiallelic, ".", sep = "")

multiallelic_data <- full@tab[, grep(paste(multiallelic, "{3,4}", sep = "", collapse = "|"), colnames(full@tab))] # mostly works, but still includes SNP_145*
multiallelic_data <- multiallelic_data[,-c(59:78)] # gets rid of additional SNP_145* that are actually biallelic but got picked up in the grep above
multiallelic_counts <- colSums(multiallelic_data, na.rm = TRUE)

multiallelic_names <- names(multiallelic_counts) #67
biallelic_only <- full.dataframe[,!colnames(full.dataframe) %in% multiallelic_names] # 528x3764 (before: 3831-67 = 3764)

# Separate genind for adults and larvae
adults <- rownames(biallelic_only)[205:439] #528 fish total
larvs <- biallelic_only[!(rownames(biallelic_only) %in% adults),]
# write.table(larvs, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/larvs_allelecounts.txt", sep="\t")
adult_data <- biallelic_only[(rownames(biallelic_only) %in% adults),]
dim(larvs) #293 x 3764
larvs <- as.genind(larvs)
larvs@pop <- full@pop[which(full@pop != 4)]
adult_data <- as.genind(adult_data)
adult_data@pop <- full@pop[which(full@pop == 4)]

# Check the dimensions of the allele frequency table
dim(larvs@tab)

# Need to be able to categorize by region, so reading in larval database so that I can match larvae names to get location
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Larvae")
larvae_database <- read.csv("Larvae Sampling Database.csv", header=TRUE)

allelefreqs_larvs <- as.data.frame(larvs@tab) # Convert larval genind object to a dataframe
allelefreqs_larvs <- allelefreqs_larvs[-264,] # removed fish caught in June

# The names are all messed up because of the bioinformatics process...
names1 <- do.call(rbind, strsplit(as.character(larvae_database$ID..), '_'))
names2 <- do.call(rbind, strsplit(as.character(names1[,1]), 'E'))
correctnames <- paste(names2[,1], names2[,2], sep = 'E_')
correctnames2 <- paste(correctnames, names1[,2], sep = '')

# Now put back the fixed names in the larval database and the split names into the allele frequency database
larvae_database$ID.. <- correctnames2

# Need to modify names in allele count data (genind object) so that I can match them up with names in the larval database
freq_names <- as.vector(rownames(larvs@tab))
freq_names_split <- do.call(rbind, strsplit(as.character(freq_names), 'L'))
rownames(allelefreqs_larvs) <- freq_names_split[-264,1] # replaces rownames in dataframe allele frequency counts with modified format
rownames(larvs@tab) <- freq_names_split[,1] # replaces rownames in genind object with different formatting

# Subset the larval database to only larvae that are staying in the analysis, and order the names so that they match with the allele frequency matrix
larvae_database_sub <- larvae_database[larvae_database$ID.. %in% freq_names_split[,1],]
ordered_larvae_database_sub <- larvae_database_sub[order(larvae_database_sub$ID..),]
ordered_larvae_database_sub <- ordered_larvae_database_sub[-264,] # removed fish caught in June

# Make sure name order is good
ordered_larvae_database_sub$ID.. == rownames(allelefreqs_larvs)

# Creating strata to add to the genind object
time_period <- as.vector(larvs@pop) # Converts population numbers to time periods
time_period <- gsub('2', 'half', time_period)
time_period <- gsub('1', 'early', time_period)
time_period <- gsub('3', 'late', time_period)

geo <- ordered_larvae_database_sub[, "Place"] # Subset the ordered larval database and assign region based on place
geo <- gsub('Beaufort, NC', 'south', geo)
geo <- gsub('North Inlet, SC', 'south', geo)
geo <- gsub('Little Egg Inlet, NJ', 'north', geo)
geo <- gsub('York River, VA', 'north', geo)
geo <- gsub('Roosevelt Inlet, DE', 'north', geo)

season <- ordered_larvae_database_sub[,"Month"]
season <- gsub('10', 'fall', season)
season <- gsub('11', 'fall', season)
season <- gsub('12', 'fall', season)
season <- gsub('1', 'winter', season)
season <- gsub('2', 'winter', season)
season <- gsub('3', 'winter', season)
season <- gsub('4', 'winter', season)

# # Remove fish caught in June from allele frequency dataset and the metadata, and time_period, geo and season objects
# which(season == '6')
# allelefreqs_larvs <- allelefreqs_larvs[-which(season == '6'),]
# ordered_larvae_database_sub <- ordered_larvae_database_sub[-which(season == '6'),]
time_period <- time_period[-264]
# geo <- geo[-which(season == '6')]
# season <- season[-which(season == '6')]
# Need to remove this fish from the larvs genind object for this to work

# Make this into a genind object again and create strata
larvs2 <- as.genind(allelefreqs_larvs)
larvs2@pop <- larvs@pop[-264]

pop_strata <- data.frame(cbind(rownames(allelefreqs_larvs), time_period, geo, season))
strata(larvs2) <- pop_strata[,-1]

# Now that the genind object is squared away, how much missing data?
# Rows are individuals, columns are alleles at a locus
summary(colSums(is.na(larvs2@tab))) # Counts number of NAs for respective allele: Max is 22 (out of 292 fish)
summary(rowSums(is.na(larvs2@tab))) # Counts number of NAs for each individual fish: Max is 1252 (out of 3764 alleles)

#### Calculate mean allele frequecies between fall vs winter ####
allelefreqs <- as.data.frame(scaleGen(larvs2, center = FALSE, scale = FALSE, NA.method = "mean")) # this goes from counts to frequencies
rownames(allelefreqs)
dim(allelefreqs)
rownames(allelefreqs) == ordered_larvae_database_sub$ID..

allelefreq.larvs.mean2 <- aggregate(allelefreqs, list(pop_strata$season), mean, na.rm = TRUE)

# Calculate SE
se <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x))) # Function to calculate SE

allelefreq.larvs.se <- aggregate(allelefreqs, list(pop_strata$season), se) #take se by season
allelefreq.larvs.se.t <- t(allelefreq.larvs.se[,-1])

# Check that means make sense
mean(allelefreqs[which(pop_strata$season == 'fall'),3001], na.rm = T) # mean of snpX during fall
mean(allelefreqs[which(pop_strata$season == 'winter'),3001], na.rm = T) # mean of snpX during winter
allelefreq.larvs.mean2[,3002] # match!

# Calculate allele frequency differences between time periods
allelefreq.larvs.mean2.t <- allelefreq.larvs.mean2[,-1] # first need to get rid of the categories so there are fewer problems down the line (num turning to chr)
allelefreq.larvs.mean2.t <- t(allelefreq.larvs.mean2.t) # columns are "fall" and "winter"

diffs <- allelefreq.larvs.mean2.t[,2] - allelefreq.larvs.mean2.t[,1] # winter minus fall
hist(diffs)
summary(diffs)

diff_outliers <- which(as.vector(diffs) > mean(diffs) + (sd(diffs, na.rm = TRUE))*3)
diff_outliers_df <- diffs[diff_outliers] 

# Plot barplots using mean and SE for loci that are 3+ SD away from the mean
par(mfrow = c(3,5))
for (i in diff_outliers){
  barplot(allelefreq.larvs.mean2.t[i,], col = 'gray60', ylim = c(0,1.1), main = paste(colnames(allelefreq.larvs.mean2.t)[i]), yaxt = 'n')
  arrows(c(0.7,1.9), allelefreq.larvs.mean2.t[i,] - allelefreq.larvs.se.t[i,], c(0.7,1.9), allelefreq.larvs.mean2.t[i,] + allelefreq.larvs.se.t[i,], length = 0.05, angle = 90, code = 3)
  axis(1, at=c(0.7,1.9), labels = c('Fall \n(n = 56)', 'Winter \n(n = 236)'), mgp=c(3, 2, 0))
  axis(2, at=seq(0,1, by = 0.2), labels=seq(0,1, by= 0.2), las = 2)
  mtext('Allele frequency', side = 2, line = 2.5)
}

#### Calculating allele frequency differences between fall and winter in three time periods ####
allelefreq.larvs.tp.mean2 <- aggregate(allelefreqs, list(pop_strata$season, pop_strata$time_period), mean, na.rm = TRUE)
allelefreq.larvs.tp.mean2.t <- allelefreq.larvs.tp.mean2[,-c(1:2)] # first need to get rid of the categories so there are fewer problems down the line (num turning to chr)
allelefreq.larvs.tp.mean2.t <- t(allelefreq.larvs.tp.mean2.t) # 6 columns: fall/early, winter/early, fall/half, winter/half, fall/late, winter/late

allelefreq.larvs.tp.se <- aggregate(allelefreqs, list(pop_strata$season, pop_strata$time_period), se) #take se by season & time period
allelefreq.larvs.tp.se.t <- allelefreq.larvs.tp.se[,-c(1:2)] # first need to get rid of the categories so there are fewer problems down the line (num turning to chr)
allelefreq.larvs.tp.se.t <- t(allelefreq.larvs.tp.se.t) # 6 columns: fall/early, winter/early, fall/half, winter/half, fall/late, winter/late

# Plot the alleles where fall and winter were outliers when all time periods are lumped
for (i in diff_outliers){
  barplot(allelefreq.larvs.tp.mean2.t[i,], col = 'gray60', ylim = c(0,1.1), main = paste(rownames(allelefreq.larvs.mean2.t)[i]), yaxt = 'n')
  arrows(c(0.7,1.9,3.1,4.3,5.5,6.7), allelefreq.larvs.tp.mean2.t[i,] - allelefreq.larvs.tp.se.t[i,], c(0.7,1.9,3.1,4.3,5.5,6.7), allelefreq.larvs.tp.mean2.t[i,] + allelefreq.larvs.tp.se.t[i,], length = 0.05, angle = 90, code = 3)
  axis(1, at=c(0.7,1.9,3.1,4.3,5.5,6.7), labels = c('early+fall \n(n = 3)', 'early+winter \n(n = 56)', 'half+fall \n(n = 19)', 'half+winter \n(n = 68)', 'late+fall \n(n = 34)', 'late+winter \n(n = 112)'), mgp=c(3, 2, 0))
  axis(2, at=seq(0,1, by = 0.2), labels=seq(0,1, by= 0.2), las = 2)
  mtext('Allele frequency', side = 2, line = 2.5)
}

