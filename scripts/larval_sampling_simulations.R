# This script randomly samples allele frequencies for the early, middle and late time periods for both north and south regions
# f2samp <- rbinom(1,c2,f2)/c2 # second sample allele frequency
# Then calculates regional allele frequency differences in the earliest and latest time periods
# And fits a linear regression for each locus, and then plots a summary figure

# c2 <- 293
# f2samp <- rbinom(n = 1882, size = c2, p = 0.5)/c2
# hist(f2samp)
# 
# north <- matrix(nrow = 1882, ncol = 3)
# north[,1] <- rbinom(n = 1882, size = c2, p = 0.5)/c2
# north[,2] <- rbinom(n = 1882, size = c2, p = 0.5)/c2
# north[,3] <- rbinom(n = 1882, size = c2, p = 0.5)/c2
# 
# south <- matrix(nrow = 1882, ncol = 3)
# south[,1] <- rbinom(n = 1882, size = c2, p = 0.5)/c2
# south[,2] <- rbinom(n = 1882, size = c2, p = 0.5)/c2
# south[,3] <- rbinom(n = 1882, size = c2, p = 0.5)/c2

#### Using a binomial distribution ####
n1 <- 8*2
n2 <- 38*2
n3 <- 100*2
n1samp <- rbinom(n = 1882, size = n1, p = 0.5)/n1
hist(n1samp)
n2samp <- rbinom(n = 1882, size = n2, p = 0.5)/n2
hist(n2samp)
n3samp <- rbinom(n = 1882, size = n3, p = 0.5)/n3
hist(n3samp)

north <- matrix(nrow = 1882, ncol = 3)
north[,1] <- rbinom(n = 1882, size = n1, p = 0.5)/n1
north[,2] <- rbinom(n = 1882, size = n2, p = 0.5)/n2
north[,3] <- rbinom(n = 1882, size = n3, p = 0.5)/n3

s1 <- 52*2
s2 <- 49*2
s3 <- 46*2

s1samp <- rbinom(n = 1882, size = s1, p = 0.5)/s1
hist(s1samp)
s2samp <- rbinom(n = 1882, size = s2, p = 0.5)/s2
hist(s2samp)
s3samp <- rbinom(n = 1882, size = s3, p = 0.5)/s3
hist(s3samp)

south <- matrix(nrow = 1882, ncol = 3)
south[,1] <- rbinom(n = 1882, size = s1, p = 0.5)/s1
south[,2] <- rbinom(n = 1882, size = s2, p = 0.5)/s2
south[,3] <- rbinom(n = 1882, size = s3, p = 0.5)/s3

# Calculate regional allele frequency differences for earliest and most recent time periods
northearly2southearly <- north[,1]-south[,1]
northlate2southlate <- north[,3]-south[,3]
hist(northearly2southearly)
hist(northlate2southlate)

# I need a for loop to fit a line to all simulated northern AND all southern allele frequencies
n <- 1882
north.lms <- lapply(1:n, function(x) lm(north[x,] ~ c(1,2,3)))
north.coeff <- sapply(north.lms, coef) # extract coefficients
north.slopes <- north.coeff[2,]

south.lms <- lapply(1:n, function(x) lm(south[x,] ~ c(1,2,3)))
south.coeff <- sapply(south.lms, coef) # extract coefficients
south.slopes <- south.coeff[2,]

axis1 <- abs(north.slopes)-abs(south.slopes)

axis2 <- abs(northearly2southearly)-abs(northlate2southlate)

plot(axis1 ~ axis2, ylab = '|Slope in north| - |Slope in south|', xlab = '|Difference between early NvS| - |Difference between late NvS|')
abline(h = 0, v = 0, col = 'blue')

#####################################################################################################
#### Simulating a null model for each locus ####
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
dim(larvs) #293 x 3764
larvs <- as.genind(larvs)
larvs@pop <- full@pop[which(full@pop != 4)]

#### Calculate mean allele frequencies for all fish ####
allelefreqs <- as.data.frame(scaleGen(larvs, center = FALSE, scale = FALSE, NA.method = "mean")) # this goes from counts to frequencies
dim(allelefreqs)

mean.allelefreqs <- as.vector(t(colMeans(allelefreqs, na.rm = TRUE)))
hist(mean.allelefreqs)

# Take only one allele per locus ####
odd_indexes<-seq(1,3763,2)

odds <- mean.allelefreqs[odd_indexes] # 1882

#### Simulate null distribution for locus 1 in each population ####
par(mfrow = c(2,3))
n1 <- 8*2
n2 <- 38*2
n3 <- 100*2
n1samp <- rbinom(n = 999, size = n1, p = odds[1])/n1
hist(n1samp)
n2samp <- rbinom(n = 999, size = n2, p = odds[1])/n2
hist(n2samp)
n3samp <- rbinom(n = 999, size = n3, p = odds[1])/n3
hist(n3samp)

s1 <- 52*2
s2 <- 49*2
s3 <- 46*2

s1samp <- rbinom(n = 999, size = s1, p = odds[1])/s1
hist(s1samp)
s2samp <- rbinom(n = 999, size = s2, p = odds[1])/s2
hist(s2samp)
s3samp <- rbinom(n = 999, size = s3, p = odds[1])/s3
hist(s3samp)

# Create empty matrices for north and south allele frequency means
north.means <- matrix(nrow = 1882, ncol = 3)
south.means <- matrix(nrow = 1882, ncol = 3)

# Populate means (or sample?)
north.means[1,1] <- mean(n1samp)
north.means[1,2] <- mean(n2samp)
north.means[1,3] <- mean(n3samp)
south.means[1,1] <- mean(s1samp)
south.means[1,2] <- mean(s2samp)
south.means[1,3] <- mean(s3samp)

##### For loop the simulations #####
north.sims <- array(dim = c(1882, 3, 999))
south.sims <- array(dim = c(1882, 3, 999))
n1 <- 8*2
n2 <- 38*2
n3 <- 100*2

s1 <- 52*2
s2 <- 49*2
s3 <- 46*2

for (i in 1:length(odds)){
north.sims[i,1,] <- rbinom(n = 999, size = n1, p = odds[i])/n1
north.sims[i,2,] <- rbinom(n = 999, size = n2, p = odds[i])/n2
north.sims[i,3,] <- rbinom(n = 999, size = n3, p = odds[i])/n3
south.sims[i,1,] <- rbinom(n = 999, size = s1, p = odds[i])/s1
south.sims[i,2,] <- rbinom(n = 999, size = s2, p = odds[i])/s2
south.sims[i,3,] <- rbinom(n = 999, size = s3, p = odds[i])/s3
}

# Calculate regional allele frequency differences for earliest and most recent time periods
northearly2southearly <- array(dim = c(1882,1,999))
for (j in 1:length(north.sims[1,1,])){
northearly2southearly[,1,j] <- north.sims[,1,j]-south.sims[,1,j]
}

northlate2southlate <- array(dim = c(1882,1,999))
for (k in 1:length(north.sims[1,3,])){
northlate2southlate[,1,k] <- north.sims[,3,k]-south.sims[,3,k]
}

hist(northearly2southearly[1,,]) # histrogram of differences for locus 1 from 999 simulations
hist(northlate2southlate[1,,]) # histrogram of differences for locus 1 from 999 simulations

# I need a for loop to fit a line to all simulated northern AND all southern allele frequencies
# Takes a long time because I'm fitting a line to 1882 x 2 x 999 times
north.slopes <- array(dim = c(1882,1,999)) # Create empty array to hold slope data
south.slopes <- array(dim = c(1882,1,999))

for (s in 1:length(north.sims[1,1,])){
n <- 1882
north.lms <- lapply(1:n, function(x) lm(north.sims[x,,s] ~ c(1,2,3)))
north.coeff <- sapply(north.lms, coef) # extract coefficients
north.slopes[,,s] <- north.coeff[2,]

south.lms <- lapply(1:n, function(x) lm(south.sims[x,,s] ~ c(1,2,3)))
south.coeff <- sapply(south.lms, coef) # extract coefficients
south.slopes[,,s] <- south.coeff[2,]
}

for (t in 1:length(north.sims[1,1,])){
axis1 <- abs(north.slopes[,,t])-abs(south.slopes[,,t]) #subtract north array1 from south array1, for loop through 999 arrays; end goal (1882 x 1 x 999)
}

for (u in 1:length(north.sims[1,,])){
axis2 <- abs(northearly2southearly[,,u])-abs(northlate2southlate[,,u]) # end goal: 1882 x 1 x 999
}

# Plot histograms (1882) for each locus made up of 999 simulations for each statistic
# Calculate 95% confidence intervals
# Plot observed statistic on histogram of simulated statistic

plot(axis1 ~ axis2, ylab = '|Slope in north| - |Slope in south|', xlab = '|Difference between early NvS| - |Difference between late NvS|')
abline(h = 0, v = 0, col = 'blue')



f2samp <- rbinom(n = 1882, size = 1, p = 0.5)/(sum(f2samp)/length(f2samp))

f2samp <- rbinom(n = 1900, size = 300, p = 0.5)/300

f2samp <- rbinom(n = 1, size = 300, p = 0.5)/300

for (i in 1:1882){
  f2samp[i] <- rbinom(n = 1, size = 293, p = 0.5)/293
}

