####################################################################################################################
#### This script does several things:                                                                           ####
#### • Reads in adult vcf, trims to contig/BP combos that are spatial outliers, writes new subset vcf           ####
#### • Reads in str file containing only spatial outliers, removes 9 of 10 GB fish                              ####
#### • Reads in location data, assigns 'north' or 'south' & calculates regional allele frequencies              ####
#### • PCA of 232 adult fish based on 15 spatial outliers                                                       ####
#### • Reads in larval&adult vcf, subsets to contig/BP combos that are spatial outliers, writes new subset vcf  ####
#### • Reads in str file containing only 10 spatial outliers, removed adult fish                                ####
#### • Calculates genotype probabilities based on these 10 loci                                                 ####
####################################################################################################################

library(ade4)
library(adegenet)
library(devtools)
library("hierfstat")
library(pegas)
library(fields)

# Doing conversions between arbitrary SNP # and contig # caused issues in the past, so can I start with adult vcf, subset to adult spatial outliers, and then use this to subset combined larval&adult vcf?
# Read in list of adult outliers -- need this to proceed with str file types
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files")
adult_contigs <- read.table('adultspatialoutliernames.txt', header = TRUE) # these are the 'SNP'names and locations of the adult spatial outliers
adult_contigs$join <- paste(adult_contigs$AdultContig, adult_contigs$AdultBP, sep = '_')

# This section only needs to be done once for str file creation
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/")

adults.vcf <- read.table("SNP.DP3g95maf05.FIL.FIL.recode.updated.vcf", skip = 58, sep = '\t') # Read in adult vcf
adults.vcf$V1 <- paste(adults.vcf$V3, adults.vcf$V4, sep = '_')# combine contig and bp# so that it's easy to subset
dim(adults.vcf) # 2120 x 252

adults.outliers.vcf <- adults.vcf[adults.vcf$V1%in% adult_contigs$join,] # keep only adult outliers
dim(adults.outliers.vcf) # 15 x 252

# write.table(adults.outliers.vcf[,-c(1:2)], 'structure_input_241_15outliers.vcf',row.names = FALSE, col.names = FALSE, sep = "\t") # get rid of first 2 columns so that vcf to str conversion will work

# Now use PGDspider to convert vcf to str file type

#### Reading in SNP data file containing all 241 PADE and the 15 environmentally-associated loci ####
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files")

outs <- read.structure("structure_input_241_15outliers.str",
                       n.ind = 241, n.loc = 15, col.lab = 1, col.pop = 2, row.marknames = 1,
                       onerowperind = FALSE)

# Fix column/SNP/contig names
outs.locusnames <- colnames(outs@tab)
outs.locusnames.split <- do.call(rbind, strsplit(as.character(outs.locusnames), '.', fixed = TRUE))

contig.vector <- vector()
for (i in 1:length(adult_contigs$join)){
  contig.vector <- append(contig.vector, rep(adult_contigs$join[i], 2))
}

colnames(outs@tab) <- paste(contig.vector, outs.locusnames.split[,2], sep = '.')

# Remove 9 GB fish
GB9 <- c("PADE_14230L1439", "PADE_14231L1440", "PADE_14232L1529", "PADE_14233L1588", "PADE_14234L1441", "PADE_14235L1442", "PADE_14236L1530", "PADE_14237L1531", "PADE_14238L1532")
PADE232_15loci <- as.data.frame(outs@tab[!rownames(outs@tab) %in% GB9,]) # remove 9 of 10 GB fish: 241 x 30
PADE232_15loci$names <- rownames(PADE232_15loci)

# Read in location data so that I can split fish into regions
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/Local Adaptation/")
locs <- read.csv('allpops_combo.csv', header = TRUE)

merge1 <- merge(locs, PADE232_15loci, by.x = "PinskyID", by.y = "names")
regions <- as.vector(merge1$bayenv_pop)

regions <- gsub('1', 'north', regions)
regions <- gsub('2', 'north', regions)
regions <- gsub('3', 'north', regions)
regions <- gsub('4', 'south', regions)
regions <- gsub('5', 'south', regions)

#### Calculate allele frequencies ####
outs.freqs <- colSums(as.data.frame(outs@tab), na.rm = TRUE)/(2*colSums(!is.na(outs@tab))) # for all loci

north.outs <- PADE232_15loci[which(regions == 'north'), -31] # remove column with names that was used for merge
south.outs <- PADE232_15loci[which(regions == 'south'), -31]

north.outs.freqs <- colSums(north.outs, na.rm = TRUE)/(2*colSums(!is.na(north.outs))) # just exact matches between larvae & adult datasets
south.outs.freqs <- colSums(south.outs, na.rm = TRUE)/(2*colSums(!is.na(south.outs)))

regional.outs.freqs <- rbind(north.outs.freqs, south.outs.freqs)

#### PCA ####
adults232_15loci <- as.genind(PADE232_15loci[,-31])
adults232_15loci@pop <- as.factor(regions)
sum(is.na(adults232_15loci$tab)) #24
X <- scaleGen(adults232_15loci, NA.method = "mean")
dim(X)
class (X)

# make PCA
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pca1$eig,main="PCA eigenvalues", col=heat.colors(50))

pca1

# Plotting PC1 and PC2
s.label(pca1$li)
title("PCA of summer flounder dataset\naxes 1-2")
add.scatter.eig(pca1$eig, 3,1,2)

s.class(pca1$li, pop(adults232_15loci))
title("PCA of summer flounder dataset\naxes 1-2")
add.scatter.eig(pca1$eig, 3,1,2)

col <- azur(15)
s.class(pca1$li, pop(adults232_15loci), xax=1,yax=2, col = transp(col,0.6), axesell=FALSE, cellipse=0, cstar=0,cpoint=3, grid=FALSE)

eig_percent <- round((pca1$eig/(sum(pca1$eig)))*100,2)
eig_percent [1:3]

#### Just plotting two regions ####
# Break PCA down by region, plotting PC1 vs PC2
plot(pca1$li[which(adults232_15loci@pop == 'north'),1], pca1$li[which(adults232_15loci@pop == 'north'),2], col = "blue", xlab = "PC1 (10.98%)", ylab = "PC2 (8.89%)", xlim = c(-6,6), ylim = c(-6,6)) # north
points(pca1$li[which(adults232_15loci@pop == 'south'),1], pca1$li[which(adults232_15loci@pop == 'south'),2], col = "red") # south

legend("topleft",
       legend=c("North", "South"),
       pch=c(1, 1),
       col=c("blue", "red"))

# Break PCA down by region, plotting PC1 vs P3
plot(pca1$li[which(adults232_15loci@pop == 'north'),1], pca1$li[which(adults232_15loci@pop == 'north'),3], col = "blue", xlab = "PC1 (10.98%)", ylab = "PC3 (8.77%)", xlim = c(-7,7), ylim = c(-4,6)) # north
points(pca1$li[which(adults232_15loci@pop == 'south'),1], pca1$li[which(adults232_15loci@pop == 'south'),3], col = "red") # south

legend("topleft",
       legend=c("North", "South"),
       pch=c(1, 1),
       col=c("blue", "red"))

##############################################################################################################################
#### I want to isolate the adult spatial outliers in the larval+adult dataset, and then assign larvae based on these loci ####
##############################################################################################################################
# This section only needs to be done once to create str file with only spatial outliers
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files")

# Read in larval/adult vcf
full.vcf <- read.table("SNP.DP3g95maf05lm75.FIL.recode.vcf", skip = 57, sep = '\t') # 3827 x 537
full.vcf$names <- paste(full.vcf$V1, full.vcf$V2, sep = '_')# combine contig and bp# so that it's easy to subset
dim(full.vcf) # 3827 x 538

# Remove all contigs that are not spatial outliers in adults
contig.vector # these are the contig/BP combos to keep
full.sub.vcf <- full.vcf[full.vcf$names %in% contig.vector,] # 10 x 538

# Compare the 'keeper' loci to adult_contigs. Expect 10 contig/BP combos
adult_contigs
adult_contigs2 <- adult_contigs[which(adult_contigs$AdultBP == adult_contigs$FullBP),] # subset adult contigs to those that exist in larvae
adult_contigs$action <- c("match", "match", "match", "retrieve from vcf", "snp not in vcf", "match", "match", "match", "match", "no such contig", "match", "match", "no such contig", "match, but 1 snp off", "no such contig")

# Write the subset vcf file to disk, then use PGDspider to convert to str file
# write.table(full.sub.vcf[,-c(538)], 'structure_input_528_10outliers.vcf',row.names = FALSE, col.names = FALSE, sep = "\t") # get rid of the last column containing names so that vcf to str conversion will work

# Now use PGDspider to convert vcf to str file type

#### Read in str file containing only 10 adult spatial outliers ####
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files")

fullsub <- read.structure("structure_input_528_10outliers.str",
                       n.ind = 528, n.loc = 10, col.lab = 1, col.pop = 2, row.marknames = 1,
                       onerowperind = FALSE)

full.freqs <- as.data.frame(fullsub@tab)
dim(full.freqs) #528 x 20

# Fix column/SNP/contig names
outs.full.locusnames <- colnames(full.freqs)
outs.full.locusnames.split <- do.call(rbind, strsplit(as.character(outs.full.locusnames), '.', fixed = TRUE))

full.contig.vector <- vector()
for (i in 1:length(full.sub.vcf$names)){
  full.contig.vector <- append(full.contig.vector, rep(full.sub.vcf$names[i], 2))
}

colnames(full.freqs) <- paste(full.contig.vector, outs.full.locusnames.split[,2], sep = '.')

# Remove adults
adults <- rownames(full.freqs)[205:439] #528 fish total

larvs.freqs <- full.freqs[!(rownames(full.freqs) %in% adults),]
dim(larvs.freqs) #293 x 20

# Subset adult spatial outliers to those only occuring in larvae
regional.outs.freqs
regional.outs.freqs10 <- regional.outs.freqs[,colnames(larvs.freqs)] # 2 x 20
colnames(regional.outs.freqs10) == colnames(larvs.freqs) # column names are in the same order
rbind(colnames(regional.outs.freqs10), colnames(larvs.freqs))

#### Allele frequency loci names between larvae & adults match ####
# First, let's look at histograms of allele frequencies. I'm guessing they should look approximately similar
dim(PADE232_15loci[,-31]) # minus names column: 232 x 30
dim(larvs.freqs) # 293 x 20

PADE232_10loci <- PADE232_15loci[,colnames(larvs.freqs)] # subsetting adult allele counts to those only in larvae&adult dataset
colnames(PADE232_10loci)==colnames(larvs.freqs) # column names are the same, so should be 1:1 for plotting

par(mfrow=c(2,2))
hist(PADE232_10loci[,1], main = paste(colnames(PADE232_10loci)[1]), xlab = "Adults")
hist(larvs.freqs[,1], main = paste(colnames(larvs.freqs)[1]), xlab = "Larvae")
hist(PADE232_10loci[,3], main = paste(colnames(PADE232_10loci)[3]), xlab = "Adults")
hist(larvs.freqs[,3], main = paste(colnames(larvs.freqs)[3]), xlab = "Larvae")
hist(PADE232_10loci[,5], main = paste(colnames(PADE232_10loci)[5]), xlab = "Adults")
hist(larvs.freqs[,5], main = paste(colnames(larvs.freqs)[5]), xlab = "Larvae")
hist(PADE232_10loci[,7], main = paste(colnames(PADE232_10loci)[7]), xlab = "Adults")
hist(larvs.freqs[,7], main = paste(colnames(larvs.freqs)[7]), xlab = "Larvae")
hist(PADE232_10loci[,9], main = paste(colnames(PADE232_10loci)[9]), xlab = "Adults")
hist(larvs.freqs[,9], main = paste(colnames(larvs.freqs)[9]), xlab = "Larvae")
hist(PADE232_10loci[,11], main = paste(colnames(PADE232_10loci)[11]), xlab = "Adults")
hist(larvs.freqs[,11], main = paste(colnames(larvs.freqs)[11]), xlab = "Larvae")
hist(PADE232_10loci[,13], main = paste(colnames(PADE232_10loci)[13]), xlab = "Adults")
hist(larvs.freqs[,13], main = paste(colnames(larvs.freqs)[13]), xlab = "Larvae")
hist(PADE232_10loci[,15], main = paste(colnames(PADE232_10loci)[15]), xlab = "Adults")
hist(larvs.freqs[,15], main = paste(colnames(larvs.freqs)[15]), xlab = "Larvae")
hist(PADE232_10loci[,17], main = paste(colnames(PADE232_10loci)[17]), xlab = "Adults")
hist(larvs.freqs[,17], main = paste(colnames(larvs.freqs)[17]), xlab = "Larvae")
hist(PADE232_10loci[,19], main = paste(colnames(PADE232_10loci)[19]), xlab = "Adults")
hist(larvs.freqs[,19], main = paste(colnames(larvs.freqs)[19]), xlab = "Larvae")

#### Calculating genotype likelihoods ####
regional.outs.freqs10

odds <- seq(1,20,2) # odd indicies to keep
regional.outs.freqs10.odds <- regional.outs.freqs10[,odds]

larvs.freqs.odds <- larvs.freqs[,odds] # only odds indicies in the larval dataset
larvs.freqs.odds[is.na(larvs.freqs.odds)] <- 9 # replace NA's with 9's to make the ifelse statements easier

colnames(larvs.freqs.odds) == colnames(regional.outs.freqs10.odds) # column names match?

# For loop to loop through each locus & multiply by the adult allele frequency, and then do this for all 293 larvae. NA's/9's get coded as 1's so they don't make a difference when each row's product is taken
north.likelihoods <- data.frame()
south.likelihoods <- data.frame()

for (j in 1:length(rownames(larvs.freqs.odds))){

for (i in 1:length(colnames(larvs.freqs.odds))){
  if(larvs.freqs.odds[j,i] == 2) {
    north.likelihoods[j,i] <- regional.outs.freqs10.odds[1,i]^2
  } else if (larvs.freqs.odds[j,i] == 1) {
      north.likelihoods[j,i] <- 2*(regional.outs.freqs10.odds[1,i] * (1-regional.outs.freqs10.odds[1,i]))
  } else if (larvs.freqs.odds[j,i] == 0) {
      north.likelihoods[j,i] <- ( 1-regional.outs.freqs10.odds[1,i])^2 
     } else {
        north.likelihoods[j,i] <- 1
    }
  }

for (i in 1:length(colnames(larvs.freqs.odds))){
  if(larvs.freqs.odds[j,i] == 2){
    south.likelihoods[j,i] <- regional.outs.freqs10.odds[2,i]^2
  } else if (larvs.freqs.odds[j,i] == 1) {
    south.likelihoods[j,i] <- 2*(regional.outs.freqs10.odds[2,i] * (1-regional.outs.freqs10.odds[2,i]))
  } else  if (larvs.freqs.odds[j,i] == 0) {
    south.likelihoods[j,i] <- (1-regional.outs.freqs10.odds[2,i])^2
  } else {
    south.likelihoods[j,i] <- 1
  }
}
}

# Multiply everything together
north.vector <- vector()
south.vector <- vector()

for (k in 1:length(north.likelihoods[,1])){
  north.vector[k] <- prod(north.likelihoods[k,])
  }

for (l in 1:length(south.likelihoods[,1])){
  south.vector[l] <- prod(south.likelihoods[l,])
  }

# Hand check a few likelihoods, including ones with NA's/9's
(0.7296296^2)*(2*0.7481481*(1-0.7481481))*(0.8777778^2)*(0.5814815^2)*(2*0.7925926*(1-0.7925926))*(0.962963^2)*(2*0.5962963*(1-0.5962963))*(0.8962963^2)*(0.9703704^2)*(0.7703704^2) #north fish1
(0.8402062^2)*(2*0.8578947*(1-0.8578947))*(0.9587629^2)*(0.4587629^2)*(2*0.6958763*(1-0.6958763))*(0.890625^2)*(2*0.6958763*(1-0.6958763))*(0.9639175^2)*(0.8762887^2)*(0.8842105^2)#south

(0.7296296^2)*(2*0.7481481*(1-0.7481481))*(0.8777778^2)*(0.5814815^2)*(0.7925926^2)*(0.962963^2)*(0.5962963^2)*(0.8962963^2)*(0.9703704^2)*(0.7703704^2) #north fish2
(0.8402062^2)*(2*0.8578947*(1-0.8578947))*(0.9587629^2)*(0.4587629^2)*(0.6958763^2)*(0.890625^2)*(0.6958763^2)*(0.9639175^2)*(0.8762887^2)*(0.8842105^2) #south

((1-0.7296296)^2)*(0.7481481^2)*(0.8777778^2)*(0.5814815^2)*(0.7925926^2)*(0.962963^2)*((1-0.5962963)^2)*(0.8962963^2)*(0.9703704^2)*(0.7703704^2) # north fish 4
((1-0.8402062)^2)*(0.8578947^2)*(0.9587629^2)*(0.4587629^2)*(0.6958763^2)*(0.890625^2)*((1-0.6958763)^2)*(0.9639175^2)*(0.8762887^2)*(0.8842105^2)

(2*0.7296296*(1-0.7296296))*(0.7481481^2)*(0.8777778^2)*(2*0.5814815*(1-0.5814815))*(2*0.7925926*(1-0.7925926))*(0.962963^2)*1*(0.8962963^2)*(0.9703704^2)*1 # north fish 29
(2*0.8402062*(1-0.8402062))*(0.8578947^2)*(0.9587629^2)*(2*0.4587629*(1-0.4587629))*(2*0.6958763*(1-0.6958763))*(0.890625^2)*1*(0.9639175^2)*(0.8762887^2)*1 #south

(0.7296296^2)*(0.7481481^2)*(0.8777778^2)*(2*0.5814815*(1-0.5814815))*(0.7925926^2)*(0.962963^2)*(0.5962963^2)*(0.8962963^2)*(0.9703704^2)*(0.7703704^2) # north fish 30
(0.8402062^2)*(0.8578947^2)*(0.9587629^2)*(2*0.4587629*(1-0.4587629))*(0.6958763^2)*(0.890625^2)*(0.6958763^2)*(0.9639175^2)*(0.8762887^2)*(0.8842105^2)

# Determine if each fish has a higher likelihood of coming from the north or the south
ratio <- north.vector/south.vector

assignment <- vector()
for (m in 1:length(ratio)){
  if (ratio[m] > 1){
    assignment[m] <- "north"
  } else {
    assignment[m] <- "south"
  }
}

table(assignment)

# Modify larvae names so that I can match them to their geographic location later
library(tidyr)

fishnames.split <- as.data.frame(do.call(rbind, strsplit(as.character(rownames(larvs.freqs.odds)), 'L', fixed = TRUE)))
fishnames.split2 <- separate(fishnames.split, V1, c("name1", "name2", "name3", "name4"),sep = c(4,5,7))
newnames <- paste(fishnames.split2$name1, fishnames.split2$name3, fishnames.split2$name2, fishnames.split2$name4, sep = '')

assignment.ids <- as.data.frame(cbind(newnames, north.vector, south.vector, assignment))

#### Compare genetic assignment vs. geographic location ####
# Read in file with larval IDs
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Larvae/")
larvs <- read.csv('Larvae Sampling Database.csv', header = TRUE)

larvs.merge <- merge(assignment.ids, larvs, by.x = "newnames", by.y = "ID..")
larvs.merge <- larvs.merge[, c("newnames", "north.vector", "south.vector", "assignment", "Sampling.Date", "Month", "Date", "Year", "Place")] # get rid of some columns to make data easier to look at
# write.table(larvs.merge, "~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/results/larval_assignment.txt", col.names = TRUE)

# Read in larval assignment
larvs.assign <- read.table("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/results/larval_assignment.txt", header = TRUE)

# Make "Place" eithr north or south of Hatteras: North = 1, South = 2
larvs.assign$Place <- gsub('Little Egg Inlet, NJ', 1, larvs.assign$Place)
larvs.assign$Place <- gsub('York River, VA', 1, larvs.assign$Place)
larvs.assign$Place <- gsub('Roosevelt Inlet, DE', 1, larvs.assign$Place)
larvs.assign$Place <- gsub('Chincoteague, VA', 1, larvs.assign$Place)
larvs.assign$Place <- gsub('Beaufort, NC', 2, larvs.assign$Place)
larvs.assign$Place <- gsub('North Inlet, SC', 2, larvs.assign$Place)

# Plot log likelihood for each population
north.log <- -log10(larvs.assign$north.vector)
south.log <- -log10(larvs.assign$south.vector)
geo.color <- as.numeric(larvs.assign$Place)

plot(south.log ~ north.log, xlab = '-log10(north likelihood)', ylab = '-log10(south likelihood)', col=ifelse(geo.color == 1, 'blue', 'tomato')) # blue is north, red is south
abline(a = 0,b=1)
legend("topleft",
       legend=c("Caught north of Hatteras", "Caught south of Hatteras"),
       pch=c(1, 1),
       col=c("blue", "tomato"))


a <- c(1,2,3,4)
b <- c(1,2,3,4)
c <- c(1,1,2,2)
plot(a~b, col= ifelse(c == 1, 'red', 'blue'))
plot(south.log[290:293] ~ north.log[290:293], xlab = '-log10(north likelihood)', ylab = '-log10(south likelihood)', col=ifelse(geo.color[290:293] == 1, 'blue', 'tomato'))