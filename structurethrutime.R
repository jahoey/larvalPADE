setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files")

library(ade4)
library(adegenet)
library(pegas)
library(hierfstat)
library(mmod)
library(poppr)

# Reading in SNP data file containing only the first SNP at each locus
# full <- read.structure("structure_input_Sept7_2016.str",
#                          n.ind = 568, n.loc = 1651, col.lab = 1, col.pop = 2, row.marknames = 1,
#                          onerowperind = FALSE)

# full <- read.structure("structure_input_Feb10_2017_temppops.str",
#                        n.ind = 568, n.loc = 1651, col.lab = 1, col.pop = 2, row.marknames = 1,
#                        onerowperind = FALSE)

# full <- read.structure("structure_input_March12_2017_554_1109.str",
#                         n.ind = 554, n.loc = 1109, col.lab = 1, col.pop = 2, row.marknames = 1,
#                         onerowperind = FALSE) # too much missing data

# full <- read.structure("structure_input_Mar13_2017_528_3827.str",
#                        n.ind = 528, n.loc = 3827, col.lab = 1, col.pop = 2, row.marknames = 1,
#                        onerowperind = FALSE)

full <- read.structure("structure_input_Mar14_2017_528_1904.str",
                       n.ind = 528, n.loc = 1904, col.lab = 1, col.pop = 2, row.marknames = 1,
                       onerowperind = FALSE)

full.dataframe <- as.data.frame(full@tab)

# adults <- rownames(full.dataframe)[205:439] #528 fish total
# larvs <- full.dataframe[!(rownames(full.dataframe) %in% adults),]
# dim(larvs) #293 x 3831, 300 x 2229
# larvs <- as.genind(larvs)
# larvs@pop <- full@pop[which(full@pop != 4)]

###############################################################
#### How much missing data per locus in larvae and adults? ####
na.count <- sapply(full.dataframe, function(y) sum(length(which(is.na(y)))))
na.count <- data.frame(na.count)
head(na.count)

sum(is.na(full.dataframe[,2]))

# Histogram of SNPs with missing data
counts <- data.matrix(na.count)
hist(counts)

# Max number of missing individuals for a locus is 26. 26/528 = 5%

#### How much missing data per individual in larvae and adults? ####
loci.na.count <- rowSums(is.na(full.dataframe))
loci.na.count <- data.frame(loci.na.count)
head(loci.na.count)

sum(is.na(full.dataframe[6,]))

# Histogram of how much missing data within an individual
counts.byindiv <- data.matrix(loci.na.count)
hist(counts.byindiv)

# Max number of alleles missing within an individual is 1289. 1289/3831 = 34%
##############################
#### How much missing data per locus in larvae only? ####
na.count.larvs <- sapply(larvs, function(y) sum(length(which(is.na(y)))))
na.count.larvs <- data.frame(na.count.larvs)
head(na.count.larvs)

sum(is.na(larvs[,3]))

# Histogram of SNPs with missing data
counts.larvs <- data.matrix(na.count.larvs)
hist(counts.larvs)

# Max number of missing individuals for a locus is 22. 22/293 = 7.5%
# Max number of missing individuals for a locus is 26 using structure_input_March12_2017_554_1109.str. 26/300 = 8.6%

#### How much missing data per individual in larvae? ####
larvs.loci.na.count <- rowSums(is.na(larvs))
larvs.loci.na.count <- data.frame(larvs.loci.na.count)
head(larvs.loci.na.count)

sum(is.na(larvs[6,]))

# Histogram of how much missing data within an individual
counts.byindiv.larvs <- data.matrix(larvs.loci.na.count)
hist(counts.byindiv.larvs)
# Max number of alleles missing within an individual is 1289. 1289/3831 = 34%
# Max number of alleles missing within an individual is 1672 using structure_input_March12_2017_554_1109.str. 1672/2229 = 75%, this is too much

#######################################
#### Working with the full dataset ####
is.genind(full)
head(indNames(full),10)
locNames(full)
sum <- summary(full)
plot(sum$Hobs ~ sum$Hexp)

# # I want to plot a histogram of the number of fish in each year/How even are the years?
# ind_528 <- data.frame(matrix(nrow=528,ncol=2))
# names(ind_528) <- c("ID", "Year")
# ind_528$ID <- indNames(full)
# write.table(ind_528, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/ind_528.txt", sep="\t")
# ind_528 <- read.table("ind_528.txt")
# ind_528
# 
# # Plotting histogram of number of samples/year
# hist(ind_528$Year, breaks=337)
# table(ind_528$Year)

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

# # For adults & larvae. Turn it back into a genind
# full_bi <- as.genind(biallelic_only)
# full_bi@pop <- full@pop


########################################################################################################
# The rest of the script explores PCA (line 146), summary statistics (line 195), DAPC (line 203),      #
# combining larval allele frequencies with database data (line 225) for use in AMOVA (line 272) and    #
# larval allele frequency changes over time (line 346)                                                 #
########################################################################################################

#### PCA ####
sum(is.na(larvs$tab)) #9424 
X <- scaleGen(larvs, NA.method = "mean")
dim(X)
class (X)

# make PCA
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))

pca1

# Plotting PC1 and PC2
s.label(pca1$li)
title("PCA of summer flounder dataset\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

s.class(pca1$li, pop(larvs))
title("PCA of summer flounder dataset\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

col <- azur(3)
s.class(pca1$li, pop(larvs), xax=1,yax=2, col = transp(col,0.6), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = TRUE)

eig_percent <- round((pca1$eig/(sum(pca1$eig)))*100,2)
eig_percent [1:3]

#### Closer look at the monomorphic loci in larvae ####
biallelic_only[,2185:2186] #SNP_1111
biallelic_only[,3405:3406] #SNP_1723

biallelic_only[,c(2185,2186,3405,3406)]


#### Break PCA down by year/period, plotting PC1 vs PC2 528 fish and 1882 alleles####
plot(pca1$li[205:439,1], pca1$li[205:439,2], col = "blue", xlab = "PC1 (1.06%)", ylab = "PC2 (0.58%)", xlim = c(-50,10), ylim = c(-30,28)) # plots 2013-2014 adults
points(pca1$li[440:499,1], pca1$li[440:499,2], col = "tomato") # plot 1989-1993
points(pca1$li[500:528,1], pca1$li[500:528,2], col = "gold") # plot 1998-1999
points(pca1$li[1:58,1], pca1$li[1:58,2], col = "gold") # plot 2000-2002
points(pca1$li[59:204,1], pca1$li[59:204,2], col = "green") # plot 2008-2012

legend("topleft",
       legend=c("Larvae (1989-1993)", "Larvae (1998-2002)", "Larvae (2008-2012)", "Adults (2013-2014)"),
       pch=c(1, 1, 1, 1, 1),
       col=c("tomato", "gold", "green", "blue"))

# Plotting PC1 vs. PC3
plot(pca1$li[205:439,1], pca1$li[205:439,3], col = "blue", xlab = "PC1 (1.06%)", ylab = "PC3 (0.51%)", xlim = c(-50,10), ylim = c(-50,20)) # plots 2013-2014 adults
points(pca1$li[440:499,1], pca1$li[440:499,3], col = "tomato") # plot 1989-1993
points(pca1$li[500:528,1], pca1$li[500:528,3], col = "gold") # plot 1998-1999
points(pca1$li[1:58,1], pca1$li[1:58,3], col = "gold") # plot 2000-2002
points(pca1$li[59:204,1], pca1$li[59:204,3], col = "green") # plot 2008-2012

legend("bottomleft",
       legend=c("Larvae (1989-1993)", "Larvae (1998-2002)", "Larvae (2008-2012)", "Adults (2013-2014)"),
       pch=c(1, 1, 1, 1, 1),
       col=c("tomato", "gold", "green", "blue"))

#### PCA of larvae broken down by region ####
# See larvalstructure.R

#### Calculate statistics ####
stats <- diff_stats(larvs)
hist(stats$per.locus[,1])

pb <- as.data.frame(larvs@tab)
jelly <- as.data.frame(as.integer(larvs@pop))
fst <- wc(cbind(jelly, pb))


#### DAPC ####
# Don't think DAPC is actually that helpful because the PCs explain so little of the variation
# K-means
grp <- find.clusters(larvs, max.n.clust = 10) # no clear asymptote, so keep all PCs, 2
table.value(table(pop(larvs), grp$grp), col.lab=paste("grp", 1:2))
dapc1 <- dapc(larvs, grp$grp)
dapc1

scatter(dapc1)
test <- optim.a.score(dapc1) #80,1 --> 14

# Choosing the number of PCs to retain
dapc2 <- dapc(larvs, n.da = 1, n.pca = 14)
scatter(dapc2)


#### This bit builds a dataframe with hierarchical population strata for use in AMOVA and looking at changes in allele frequencies over time ####
# Differences between larvae?
dim(larvs@tab)

# Need to be able to categorize by region, so reading in larval database so that I can match larvae names to get location
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Larvae")
larvae_database <- read.csv("Larvae Sampling Database.csv", header=TRUE)

allelefreqs_larvs <- as.data.frame(larvs@tab) # Convert larval genind object to a dataframe

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
rownames(allelefreqs_larvs) <- freq_names_split[,1] # replaces rownames in dataframe allele frequency counts with modified format
rownames(larvs@tab) <- freq_names_split[,1] # replaces rownames in genind object with different formatting

# Subset the larval database to only larvae that are staying in the analysis, and order the names so that they match with the allele frequency matrix
larvae_database_sub <- larvae_database[larvae_database$ID.. %in% freq_names_split[,1],]
ordered_larvae_database_sub <- larvae_database_sub[order(larvae_database_sub$ID..),]

# Make sure name order is good
ordered_larvae_database_sub$ID.. == rownames(allelefreqs_larvs)

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

pop_strata <- data.frame(cbind(rownames(allelefreqs_larvs), time_period, geo))
strata(larvs) <- pop_strata[,-1]

#### AMOVA ####
larvs_dist <- dist(larvs)
larvs_strata <- strata(larvs)
larvs_amova <- pegas::amova(larvs_dist ~ geo, data = larvs_strata, nperm = 1000) # no difference between region or time period
larvs_amova

# AMOVA another way using the poppr package
table(strata(larvs, ~time_period))
table(strata(larvs, ~time_period/geo, combine = FALSE))
larvs_amova2 <- poppr.amova(larvs, ~time_period/geo, within = FALSE)
larvs_amova2
larvs_amova2_signif <- randtest(larvs_amova2, nrepet = 1000)
plot(larvs_amova2_signif)

fstat(larvs)

# AMOVA with just northern (or southern) fish broken out by time period
# First need to create different datasets for northern (and southern) fish
larvs_north <- pop_strata[which(pop_strata$geo == 'north'),]
larvs_south <- pop_strata[which(pop_strata$geo == 'south'),]

allelefreq.larvs.north <- allelefreqs_larvs[rownames(allelefreqs_larvs) %in% larvs_north$V1, ]
allelefreq.larvs.south <- allelefreqs_larvs[rownames(allelefreqs_larvs) %in% larvs_south$V1, ]

larvs_northonly <- as.genind(allelefreq.larvs.north)
strata(larvs_northonly) <- larvs_north[-1]
larvs_north_dist <- dist(larvs_northonly)
larvs_north_strata <- strata(larvs_northonly)
larvs_north_amova <- pegas::amova(larvs_north_dist ~ time_period, data = larvs_north_strata, nperm = 1000) # time period is not a significant contributor of explaining the variance in the north
larvs_north_amova

larvs_southonly <- as.genind(allelefreq.larvs.south)
strata(larvs_southonly) <- larvs_south[-1]
larvs_south_dist <- dist(larvs_southonly)
larvs_south_strata <- strata(larvs_southonly)
larvs_south_amova <- pegas::amova(larvs_south_dist ~ time_period, data = larvs_south_strata, nperm = 1000) # time period is not a significant contributor of explaining the variance in the south
larvs_south_amova



# This bit works but I'm not sure it is the AMOVA I want to be doing
full.freqs <- as.data.frame(full@tab)
larvs.dist <- dist(larvs.freqs)
strata(larvs.freqs) <- data.frame(pop(full)) #Setting the strata to populations
nameStrata(full) <- ~Population


full_dist <- dist(full) #Calculates a matrix of Euclidean distances between individuals using the genetic data, missing data is ignored
str(full_dist)
strata(full) <- data.frame(pop(full)) #Setting the strata to populations
nameStrata(full) <- ~Population
full_stra <- strata(full)

full_amova <- amova(full_dist ~ Population, data = full_stra, nperm = 0) #Just testing for population differentiation
full_amova #Populations can only explain a small amount of variance, most of the variance is within a population --> stability of population structure

# Permutation test of AMOVA
full_amova <- amova(full_dist ~ Population, data = full_stra) # 1000 permutations by default
full_amova # results are robust, populations are significantly different

# Trying another implementation of AMOVA
library("vegan")
amova.vegan <- adonis(full_dist ~ Population, data = full_stra, permutations = 1000)
amova.vegan #Populations are different

library("poppr")
amova.poppr <- poppr.amova(full, ~Population, within = FALSE)
amova.poppr.test <- randtest(amova.poppr, nrepet = 1000) #tests for significance
plot(amova.poppr.test) # variation between populations is significant, but this variation is small compared to variation within pops

# Might make sense to create another hierarchical level. Right now it's just "population", which refers to time period. Maybe time period and geographic location/region.
# Or grouping fish by larvae vs adults too


#### Allele frequecies over time ####
allelefreqs <- as.data.frame(scaleGen(larvs, center = FALSE, scale = FALSE, NA.method = "mean")) # this goes from counts to frequencies
rownames(allelefreqs)
dim(allelefreqs)
rownames(allelefreqs) == ordered_larvae_database_sub$ID..

allelefreq.larvs.mean2 <- aggregate(allelefreqs, list(pop_strata$time_period), mean, na.rm = TRUE)

# This way is not as good as above
# allelefreq.larvs.mean <- aggregate(allelefreqs.larvs.place[,3:3319], list(allelefreqs.larvs.place$period), mean) #way 1, this is bad if the column indexes change
# dim(allelefreq.larvs.mean) # 3 x 3318

# Check that means make sense
mean(allelefreqs[59:204,3]) # mean of snp3 during late period
allelefreq.larvs.mean2[,3] # match!
mean(allelefreqs[205:264,1300]) # early, snpn
allelefreq.larvs.mean2[,1301]

# Calculate allele frequency differences between time periods
allelefreq.larvs.mean2.t <- allelefreq.larvs.mean2[,-1] # first need to get rid of the categories so there are fewer problems down the line (num turning to chr)
allelefreq.larvs.mean2.t <- t(allelefreq.larvs.mean2.t) # columns are "early", "half" and "late"

half2early <- allelefreq.larvs.mean2.t[,2] - allelefreq.larvs.mean2.t[,1] # first half of the time period
late2half <- allelefreq.larvs.mean2.t[,3] - allelefreq.larvs.mean2.t[,2] # second half of the time period
late2early <- allelefreq.larvs.mean2.t[,3] - allelefreq.larvs.mean2.t[,1] # allele frequency difference over full time period

hist(half2early)
hist(late2half)
hist(late2early)

boxplot(half2early, xlab = "half2early")
boxplot(late2half, xlab = "late2half")
boxplot(late2early, xlab = "late2early")

summary(half2early)
summary(late2early) # when NA's are ignored, max difference between early and late time periods is 0.18
summary(late2half)

which(late2early > 0.10)
which(half2early > 0.13)
which(late2half > 0.12)

#### What about mean allele frequency changes between north & south and by time period?
# Merge the larval samping database and allele frequency dataframe
allelefreqs$ID <- rownames(allelefreqs) # adding an ID column so the IDs are not just rownames. The names are already sorted and match the place data, so another way to do this would be to cbind.
allelefreqs_database <- merge(pop_strata, allelefreqs, by.x = "V1", by.y = "ID")
dim(allelefreqs_database) # 293 x 3767

# Dividing allele frequency dataset into north and south based on sampling location
allelefreqs.larvs.north <- allelefreqs_database[which(allelefreqs_database$geo == 'north'), ]
dim(allelefreqs.larvs.north) # 146 x 3767

allelefreqs.larvs.south <- allelefreqs_database[which(allelefreqs_database$geo == 'south'), ]
dim(allelefreqs.larvs.south) # 147 x 3767

drops <- c("V1", "time_period", "geo")
allelefreq.larvs.mean.north <- allelefreqs.larvs.north[, !(colnames(allelefreqs.larvs.north) %in% drops)] # removing columns that aren't related to allele frequency
allelefreq.larvs.mean.north <- aggregate(allelefreq.larvs.mean.north, list(allelefreqs.larvs.north$time_period), mean, na.rm = TRUE) #take mean by period
allelefreq.larvs.mean.south <- allelefreqs.larvs.south[, !(colnames(allelefreqs.larvs.south) %in% drops)] # removing columns that aren't related to allele frequency
allelefreq.larvs.mean.south <- aggregate(allelefreq.larvs.mean.south, list(allelefreqs.larvs.south$time_period), mean, na.rm = TRUE) #take mean by period

# Calculate allele frequency differences between time periods & region
allelefreq.larvs.mean.north.t <- allelefreq.larvs.mean.north[,-1] # first need to get rid of the period categories so there are fewer problems down the line (num turning to chr)
allelefreq.larvs.mean.south.t <- allelefreq.larvs.mean.south[,-1]

allelefreq.larvs.mean.north.t <- t(allelefreq.larvs.mean.north.t) # Transpose: columns are "early", "half" and "late"
allelefreq.larvs.mean.south.t <- t(allelefreq.larvs.mean.south.t)

# Difference between time periods within a region
north.half2north.early <- allelefreq.larvs.mean.north.t[,2] - allelefreq.larvs.mean.north.t[,1]
north.late2north.early <- allelefreq.larvs.mean.north.t[,3] - allelefreq.larvs.mean.north.t[,1]
north.late2north.half <- allelefreq.larvs.mean.north.t[,3] - allelefreq.larvs.mean.north.t[,2]

south.half2south.early <- allelefreq.larvs.mean.south.t[,2] - allelefreq.larvs.mean.south.t[,1]
south.late2south.early <- allelefreq.larvs.mean.south.t[,3] - allelefreq.larvs.mean.south.t[,1]
south.late2south.half <- allelefreq.larvs.mean.south.t[,3] - allelefreq.larvs.mean.south.t[,2]

hist(north.half2north.early)
hist(north.late2north.early, breaks = 20)
hist(north.late2north.half) 

boxplot(north.half2north.early)
boxplot(north.late2north.early)
boxplot(north.late2north.half) 

hist(south.half2south.early)
hist(south.late2south.early)
hist(south.late2south.half)

boxplot(south.half2south.early)
boxplot(south.late2south.early)
boxplot(south.late2south.half)

summary(north.late2north.early) # when NA's ignored, max freq difference is 0.365
summary(south.late2south.early) # when NA's ignored, max freq difference is 0.294

# Plot loci whose allele frequencies have changed the most over time
# Allele frequencies that have changed the most over time in the south
# south_temporal_outliers <- which(abs(as.vector(south.late2south.early)) > (sd(south.late2south.early))*3) # Gives both alleles
# south_temporal_outliers_single <- south_temporal_outliers[seq(1, length(south_temporal_outliers), 2)]
south_temporal_outliers <- which(as.vector(south.late2south.early) > (sd(south.late2south.early))*3)

par(mfrow = c(3,4))
for (i in south_temporal_outliers){
  plot(allelefreq.larvs.mean.north.t[10,], ylim=c(0,1), col = NULL, xlab = 'Time Period', ylab = 'Allele Frequency', xaxt = 'n', main = paste(rownames(allelefreq.larvs.mean.south.t)[i]))
  axis(1, at=1:3, labels = c('1989-1993', '1998-2002', '2008-2012'))
  lines(allelefreq.larvs.mean.north.t[i,])
  lines(allelefreq.larvs.mean.south.t[i,], col = 'tomato')
}

# Allele frequencies that have changed the most over time in the north
north_temporal_outliers <- which(as.vector(north.late2north.early) > (sd(north.late2north.early))*3)

par(mfrow = c(4,4))
for (i in north_temporal_outliers){
  plot(allelefreq.larvs.mean.north.t[10,], ylim=c(0,1), col = NULL, xlab = 'Time Period', ylab = 'Allele Frequency', xaxt = 'n', main = paste(rownames(allelefreq.larvs.mean.north.t)[i]))
  axis(1, at=1:3, labels = c('1989-1993', '1998-2002', '2008-2012'))
  lines(allelefreq.larvs.mean.north.t[i,])
  lines(allelefreq.larvs.mean.south.t[i,], col = 'tomato')
}


# which(abs(north.late2north.early) > 0.15)
# png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/Larvae/regionalallelefreqthrutime.png", width=6, height=5, res=300, units="in")
# par(mar=c(4.1, 4.1, 2, 9), xpd=TRUE,
#     mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
#     tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
#     cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
#     ps=14, # point size, which is the font size
#     bg=NA)
# plot(allelefreq.larvs.mean.north.t[190,], ylim=c(0,1), col = NULL, xlab = 'Time Period', ylab = 'Allele Frequency', xaxt = 'n')
# axis(1, at=1:3, labels = c('1989-1993', '1998-2002', '2008-2012'))
# legend("topright",
#        c("North", "South"),
#        lty=c(1,1), 
#        lwd=c(2,2),col=c("black", "tomato"),
#        inset = c(-0.45,0))
# lines(allelefreq.larvs.mean.north.t[190,])
# lines(allelefreq.larvs.mean.north.t[289,])
# lines(allelefreq.larvs.mean.north.t[653,])
# lines(allelefreq.larvs.mean.north.t[840,])
# lines(allelefreq.larvs.mean.north.t[893,])
# lines(allelefreq.larvs.mean.north.t[929,])
# lines(allelefreq.larvs.mean.north.t[973,])
# lines(allelefreq.larvs.mean.north.t[1035,])
# lines(allelefreq.larvs.mean.north.t[1179,])
# lines(allelefreq.larvs.mean.north.t[1225,])
# lines(allelefreq.larvs.mean.north.t[1243,])
# lines(allelefreq.larvs.mean.north.t[2217,])
# lines(allelefreq.larvs.mean.north.t[2265,])
# lines(allelefreq.larvs.mean.north.t[2569,])
# lines(allelefreq.larvs.mean.north.t[2832,])
# lines(allelefreq.larvs.mean.north.t[2863,])
# lines(allelefreq.larvs.mean.north.t[3054,])
# lines(allelefreq.larvs.mean.north.t[3249,])
# 
# # Plotting the allele frequency of the same alleles in the south
# lines(allelefreq.larvs.mean.south.t[190,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[289,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[653,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[840,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[893,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[929,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[973,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[1035,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[1179,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[1225,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[1243,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[2217,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[2265,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[2569,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[2832,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[2863,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[3054,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[3249,], col='tomato')
# 
# dev.off()

# # Hard to see patterns when all plotted together, so here's just a few alleles over time so that you can see the convergence
# png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/Larvae/regionalallelefreqthrutime_sample.png", width=6, height=5, res=300, units="in")
# par(mar=c(4.1, 4.1, 2, 9), xpd=TRUE,
#     mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
#     tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
#     cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
#     ps=14, # point size, which is the font size
#     bg=NA)
# plot(allelefreq.larvs.mean.north.t[190,], ylim=c(0,1), col = NULL, xlab = 'Time Period', ylab = 'Allele Frequency', xaxt = 'n')
# axis(1, at=1:3, labels = c('1989-1993', '1998-2002', '2008-2012'))
# legend("topright",
#        c("North", "South"),
#        lty=c(1,1), 
#        lwd=c(2,2),col=c("black", "tomato"),
#        inset = c(-0.45,0))
# lines(allelefreq.larvs.mean.north.t[190,])
# lines(allelefreq.larvs.mean.north.t[289,])
# lines(allelefreq.larvs.mean.north.t[893,])
# lines(allelefreq.larvs.mean.north.t[929,])
# lines(allelefreq.larvs.mean.north.t[2863,])
# 
# lines(allelefreq.larvs.mean.south.t[190,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[289,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[893,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[929,], col='tomato')
# lines(allelefreq.larvs.mean.south.t[2863,], col='tomato')
# 
# dev.off()

# Same time period, different regions
plot(allelefreq.larvs.mean.north.t [,1]~ allelefreq.larvs.mean.south.t[,1])
plot(allelefreq.larvs.mean.north.t [,2]~ allelefreq.larvs.mean.south.t[,2])
plot(allelefreq.larvs.mean.north.t [,3]~ allelefreq.larvs.mean.south.t[,3])

# Same time period, different regions, single allele
plot(allelefreq.larvs.mean.north.t [,1][seq(1, length(allelefreq.larvs.mean.north.t [,1]), 2)]~ allelefreq.larvs.mean.south.t[,1][seq(1, length(allelefreq.larvs.mean.south.t [,1]), 2)])
abline(0,1)
plot(allelefreq.larvs.mean.north.t [,2][seq(1, length(allelefreq.larvs.mean.north.t [,1]), 2)]~ allelefreq.larvs.mean.south.t[,2][seq(1, length(allelefreq.larvs.mean.south.t [,1]), 2)])
abline(0,1)
plot(allelefreq.larvs.mean.north.t [,3][seq(1, length(allelefreq.larvs.mean.north.t [,1]), 2)]~ allelefreq.larvs.mean.south.t[,2][seq(1, length(allelefreq.larvs.mean.south.t [,1]), 2)])
abline(0,1)

# Spatial differences within a time period
north.late2south.late <- allelefreq.larvs.mean.north.t [,3] - allelefreq.larvs.mean.south.t[,3]
north.half2south.half <- allelefreq.larvs.mean.north.t [,2] - allelefreq.larvs.mean.south.t[,2]
north.early2south.early <- allelefreq.larvs.mean.north.t [,1] - allelefreq.larvs.mean.south.t[,1]

hist(north.early2south.early, breaks = 20)
hist(north.half2south.half, breaks = 20)
hist(north.late2south.late, breaks = 20) # the tails of the north/south differences are clearly smaller at the end of the time period

summary(north.early2south.early)
summary(north.half2south.half)
summary(north.late2south.late) 

which(abs(north.late2south.late) > 0.14)
which(abs(north.half2south.half) > 0.16)
which(abs(north.early2south.early) > 0.3) # spatial outliers between north vs south

# Visualizing allele frequencies over space/time
png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/Larvae/regional_allelefreq_dif_early.png", width=5, height=4.5, res=300, units="in")
hist(north.early2south.early, xlab = 'Regional allele frequency difference (1989-1993)', main = NULL)
dev.off()

# # Across time periods and regions
# plot(allelefreq.larvs.mean.north.t [,3]~ allelefreq.larvs.mean.south.t[,1])
# 
# # Difference in allele frequencies between north and south within a time period
# ns <- allelefreq.larvs.mean.north.t - allelefreq.larvs.mean.south.t
# which((ns[,3]- ns[,1]) > 0.20)
# 
# sn <- allelefreq.larvs.mean.south.t - allelefreq.larvs.mean.north.t
# which((sn[,3]- sn[,1]) > 0.20)

##### Plotting greatest allele frequency changes in the north for the 293 larval dataset #####
early_spatial_outliers <- which(abs(as.vector(north.early2south.early)) > (sd(north.early2south.early)*3))
early_spatial_outliers_single <- early_spatial_outliers[seq(1, length(early_spatial_outliers), 2)]

par(mfrow = c(4,5))
for (i in early_spatial_outliers_single){
  plot(allelefreq.larvs.mean.north.t[10,], ylim=c(0,1), col = NULL, xlab = 'Time Period', ylab = 'Allele Frequency', xaxt = 'n', main = paste(rownames(allelefreq.larvs.mean.south.t)[i]))
  axis(1, at=1:3, labels = c('1989-1993', '1998-2002', '2008-2012'))
  lines(allelefreq.larvs.mean.north.t[i,])
  lines(allelefreq.larvs.mean.south.t[i,], col = 'tomato')
  # legend("bottomright",
  #        c("North", "South"),
  #        lty=c(1,1), 
  #        lwd=c(2,2),col=c("black", "tomato"))
}


lines(allelefreq.larvs.mean.north.t[127,])
lines(allelefreq.larvs.mean.north.t[245,])
lines(allelefreq.larvs.mean.north.t[319,])
lines(allelefreq.larvs.mean.north.t[347,])
lines(allelefreq.larvs.mean.north.t[563,])
lines(allelefreq.larvs.mean.north.t[761,])
lines(allelefreq.larvs.mean.north.t[867,])
lines(allelefreq.larvs.mean.north.t[1979,])
lines(allelefreq.larvs.mean.north.t[2045,])
lines(allelefreq.larvs.mean.north.t[2111,])
lines(allelefreq.larvs.mean.north.t[2359,])

lines(allelefreq.larvs.mean.south.t[127,], col = 'tomato')
lines(allelefreq.larvs.mean.south.t[245,], col = 'tomato')
lines(allelefreq.larvs.mean.south.t[319,], col = 'tomato')
lines(allelefreq.larvs.mean.south.t[347,], col = 'tomato')
lines(allelefreq.larvs.mean.south.t[563,], col = 'tomato')
lines(allelefreq.larvs.mean.south.t[761,], col = 'tomato')
lines(allelefreq.larvs.mean.south.t[867,], col = 'tomato')
lines(allelefreq.larvs.mean.south.t[1979,], col = 'tomato')
lines(allelefreq.larvs.mean.south.t[2045,], col = 'tomato')
lines(allelefreq.larvs.mean.south.t[2111,], col = 'tomato')
lines(allelefreq.larvs.mean.south.t[2359,], col = 'tomato')

which(abs(south.late2south.early) > 0.15)

# Taking a look at spatial outliers from the latest time period
par(mfrow = c(3,3))
late_spatial_outliers <- which(abs(as.vector(north.late2south.late)) > 0.14)
late_spatial_outliers_single <- late_spatial_outliers[seq(1, length(late_spatial_outliers), 2)]
for (i in late_spatial_outliers_single){
  plot(allelefreq.larvs.mean.north.t[10,], ylim=c(0,1), col = NULL, xlab = 'Time Period', ylab = 'Allele Frequency', xaxt = 'n', main = paste('Locus', i))
  axis(1, at=1:3, labels = c('1989-1993', '1998-2002', '2008-2012'))
  lines(allelefreq.larvs.mean.north.t[i,])
  lines(allelefreq.larvs.mean.south.t[i,], col = 'tomato')
}

####################################################################
#### Plotting allele frequency by latitude for each time period ####
allelefreqs <- as.data.frame(scaleGen(larvs, center = FALSE, scale = FALSE, NA.method = "mean")) # this goes from counts to frequencies
rownames(allelefreqs)

early <- c('1989', '1990', '1991', '1992', '1993')
middle <- c('1998', '1999', '2000', '2001', '2002')
late <- c('2008', '2009', '2010', '2011', '2012')

early_allelefreqs <- allelefreqs[which(ordered_larvae_database_sub$Year %in% early),]
middle_allelefreqs <- allelefreqs[which(ordered_larvae_database_sub$Year %in% middle),]
late_allelefreqs <- allelefreqs[which(ordered_larvae_database_sub$Year %in% late),]

early_database <- ordered_larvae_database_sub[which(ordered_larvae_database_sub$Year %in% early),]
middle_database <- ordered_larvae_database_sub[which(ordered_larvae_database_sub$Year %in% middle),]
late_database <- ordered_larvae_database_sub[which(ordered_larvae_database_sub$Year %in% late),]

early_allelefreqs_mean <- aggregate(early_allelefreqs, list(early_database$Lat), mean, na.rm = TRUE)
dim(early_allelefreqs_mean)
middle_allelefreqs_mean <- aggregate(middle_allelefreqs, list(middle_database$Lat), mean, na.rm = TRUE)
dim(middle_allelefreqs_mean)
late_allelefreqs_mean <- aggregate(late_allelefreqs, list(late_database$Lat), mean, na.rm = TRUE)
dim(late_allelefreqs_mean)

plot(early_allelefreqs_mean$SNP_64.03 ~ early_allelefreqs_mean$Group.1, ylim = c(0, 1), xlab = "Latitude", col = NULL, xlim = c(33,40))
lines(early_allelefreqs_mean$SNP_64.03 ~ early_allelefreqs_mean$Group.1)
lines(middle_allelefreqs_mean$SNP_64.03 ~ middle_allelefreqs_mean$Group.1, col = "blue")
lines(late_allelefreqs_mean$SNP_64.03 ~ late_allelefreqs_mean$Group.1, col = "tomato")
legend("bottomleft",
       c("Early", "Middle", "Late"),
       lty=c(1,1,1),
       lwd=c(2,2, 2),col=c("black", "blue","tomato"))

plot(early_allelefreqs_mean$SNP_69.02 ~ early_allelefreqs_mean$Group.1, ylim = c(0, 1), xlab = "Latitude", col = NULL, xlim = c(33,40))
lines(early_allelefreqs_mean$SNP_69.02 ~ early_allelefreqs_mean$Group.1)
lines(middle_allelefreqs_mean$SNP_69.02 ~ middle_allelefreqs_mean$Group.1, col = "blue")
lines(late_allelefreqs_mean$SNP_69.02 ~ late_allelefreqs_mean$Group.1, col = "tomato")

