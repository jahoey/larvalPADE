setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files")

library(ade4)
library(adegenet)
library(devtools)
library("hierfstat")
library(pegas)
library(fields)

# Reading in SNP data file containing all 241 PADE and the 23 environmentally-associated loci
outs <- read.structure("structure_input_241_15outlier.str",
                       n.ind = 241, n.loc = 15, col.lab = 1, col.pop = 2, row.marknames = 1,
                       onerowperind = FALSE)

# Split str file into north vs south fish
north.outs <- as.genind(outs@tab[1:144,c(1,2,3,4,5,6,11,12,13,14,15,16,17,18,23,24)])
north.pop <- as.vector(outs@pop[1:144])
north.outs@pop <- as.factor(north.pop)

south.outs <- as.genind(outs@tab[145:241,c(1,2,3,4,5,6,11,12,13,14,15,16,17,18,23,24)])
south.pop <- as.vector(outs@pop[145:241])
south.outs@pop <- as.factor(south.pop)

# Calculate frequencies
outs.freqs <- colSums(as.data.frame(outs@tab), na.rm = TRUE)/(2*colSums(!is.na(outs@tab))) # for all loci

north.outs.freqs <- colSums(as.data.frame(north.outs@tab), na.rm = TRUE)/(2*colSums(!is.na(north.outs@tab))) # just exact matches between larvae & adult datasets
south.outs.freqs <- colSums(as.data.frame(south.outs@tab), na.rm = TRUE)/(2*colSums(!is.na(south.outs@tab)))

regional.outs.freqs <- rbind(north.outs.freqs, south.outs.freqs)

# Subset to only alleles found in adult and larval dataset (9 loci)
# regional.outs.freqs2 <- regional.outs.freqs[,c(1,2,3,4,5,6,11,12,13,14,15,16,17,18,23,24)] # SNP_919 = 27 & 28

#### PCA
sum(is.na(outs$tab)) #24
X <- scaleGen(outs, NA.method = "mean")
dim(X)
class (X)

# make PCA
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pca1$eig,main="PCA eigenvalues", col=heat.colors(50))

pca1

# Plotting PC1 and PC2
s.label(pca1$li)
title("PCA of summer flounder dataset\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

s.class(pca1$li, pop(outs))
title("PCA of summer flounder dataset\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

col <- azur(15)
s.class(pca1$li, pop(outs), xax=1,yax=2, col = transp(col,0.6), axesell=FALSE, cellipse=0, cstar=0,cpoint=3, grid=FALSE)

eig_percent <- round((pca1$eig/(sum(pca1$eig)))*100,2)
eig_percent [1:3]

# Break PCA down by region, plotting PC1 vs PC2
plot(pca1$li[1:49,1], pca1$li[1:49,2], col = "violet", xlab = "PC1 (10.61%)", ylab = "PC2 (8.79%)", xlim = c(-6,6), ylim = c(-5,5)) # pop1, 10 weird Georges Bank fish are outs@tab[25:34,]
points(pca1$li[50:103,1], pca1$li[50:103,2], col = "blue") # pop2
points(pca1$li[104:144,1], pca1$li[104:144,2], col = "green") # pop3
points(pca1$li[145:204,1], pca1$li[145:204,2], col = "gold") # pop4
points(pca1$li[205:241,1], pca1$li[205:241,2], col = "tomato") # pop5
points(pca1$li[25:34,1], pca1$li[25:34,2], col = "black") # 10 weird GB fish

legend("topleft",
       legend=c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5"),
       pch=c(1, 1, 1, 1, 1),
       col=c("violet", "blue", "green", "gold", "tomato"))

# Break PCA down by region, plotting PC1 vs P3
plot(pca1$li[1:49,1], pca1$li[1:49,3], col = "violet", xlab = "PC1 (10.61%)", ylab = "PC3 (8.73%)", xlim = c(-7,8), ylim = c(-6,5)) # pop1
points(pca1$li[50:103,1], pca1$li[50:103,3], col = "blue") # pop2
points(pca1$li[104:144,1], pca1$li[104:144,3], col = "green") # pop3
points(pca1$li[145:204,1], pca1$li[145:204,3], col = "gold") # pop4
points(pca1$li[205:241,1], pca1$li[205:241,3], col = "tomato") # pop5
points(pca1$li[25:34,1], pca1$li[25:34,3], col = "black") # 10 weird GB fish

legend("bottomright",
       legend=c("Pop1", "Pop2", "Pop3", "Pop4", "Pop5"),
       pch=c(1, 1, 1, 1, 1),
       col=c("violet", "blue", "green", "gold", "tomato"))

# 3D PCA
library(rgl)
color <- c(rep("violet", 49), rep("blue", 54), rep("green", 41), rep("gold", 60), rep("tomato", 37))
plot3d(pca1$li[,1:3], col = color, xlab = "PC1", ylab = "PC2", zlab = "PC3")

#### Just plotting two regions ####
# Break PCA down by region, plotting PC1 vs PC2
plot(pca1$li[1:49,1], pca1$li[1:49,2], col = "blue", xlab = "PC1 (10.61%)", ylab = "PC2 (8.79%)", xlim = c(-6,6), ylim = c(-5,5)) # pop1, 10 weird Georges Bank fish are outs@tab[25:34,]
points(pca1$li[50:103,1], pca1$li[50:103,2], col = "blue") # pop2
points(pca1$li[104:144,1], pca1$li[104:144,2], col = "blue") # pop3
points(pca1$li[145:204,1], pca1$li[145:204,2], col = "red") # pop4
points(pca1$li[205:241,1], pca1$li[205:241,2], col = "red") # pop5
points(pca1$li[25:34,1], pca1$li[25:34,2], col = "black") # 10 weird GB fish

legend("topleft",
       legend=c("North", "South"),
       pch=c(1, 1),
       col=c("blue", "red"))

# Break PCA down by region, plotting PC1 vs P3
plot(pca1$li[1:49,1], pca1$li[1:49,3], col = "blue", xlab = "PC1 (10.61%)", ylab = "PC3 (8.73%)", xlim = c(-7,8), ylim = c(-6,5)) # pop1
points(pca1$li[50:103,1], pca1$li[50:103,3], col = "blue") # pop2
points(pca1$li[104:144,1], pca1$li[104:144,3], col = "blue") # pop3
points(pca1$li[145:204,1], pca1$li[145:204,3], col = "red") # pop4
points(pca1$li[205:241,1], pca1$li[205:241,3], col = "red") # pop5
points(pca1$li[25:34,1], pca1$li[25:34,3], col = "black") # 10 weird GB fish

legend("bottomright",
       legend=c("North", "South"),
       pch=c(1, 1),
       col=c("blue", "red"))

#################################################################################################################
#### I want to isolate the adult spatial outliers in the larval+adult dataset, and then plot larvae in a PCA ####
#################################################################################################################
library(ade4)
library(adegenet)
library("hierfstat")
library(pegas)

setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files")

adult_contigs <- read.table('adultspatialoutliernames.txt', header = TRUE) # these are the 'SNP'names and locations of the adult spatial outliers
adult_contigs[] <- lapply(adult_contigs, as.character)
adult_contigs2 <- adult_contigs[which(adult_contigs$AdultBP == adult_contigs$FullBP),] # subset adult contigs to those that exist in larvae

# # Read in structure formatted file so that I can subset to adult spatial outliers
# full <- read.table('structure_528_adultspatialoutliers.txt', header = TRUE) # 1056 x 1906

#### Build the str file ####
# Subset the file
# snpnames <- c("ID", "Pop", "SNP_170", "SNP_458",  "SNP_685",  "SNP_833",  "SNP_999", "SNP_1062", "SNP_1129", "SNP_1171", "SNP_1194", "SNP_1275", "SNP_1281", "SNP_1441", "SNP_1639", "SNP_1720")
# full.sub <- full[,snpnames] # 1056 x 16
# 
# # Write the txt file and then change it to str file format and read back in
# write.table(full.sub, 'structure_528_adultspatialoutliers.txt',row.names = FALSE, col.names = TRUE, sep = "\t")

# And read it back in
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/data_files")

fullsub <- read.structure("structure_528_adultspatialoutliers.str",
                       n.ind = 528, n.loc = 14, col.lab = 1, col.pop = 2, row.marknames = 1,
                       onerowperind = FALSE)

full.freqs <- as.data.frame(fullsub@tab)
dim(full.freqs) #528 x 3764

adults <- rownames(full.freqs)[205:439] #528 fish total

larvs.freqs <- full.freqs[!(rownames(full.freqs) %in% adults),]
dim(larvs.freqs) #293 x 3764
larvs.freqs2 <- as.genind(larvs.freqs)

# Removing the adult populations and putting pop back into genin object
nope <- pop(fullsub)[pop(fullsub)==4]
pop <- pop(fullsub)[!as.vector(pop(fullsub)) %in% nope]
population <- as.vector(pop)
larvs.freqs2@pop <- as.factor(population)

#### PCA on for all fish using adult spatial outliers ####
sum(is.na(larvs.freqs2$tab)) #38
X <- scaleGen(larvs.freqs2, NA.method = "mean")
dim(X) #293 x 28
class (X)

# make PCA
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pca1$eig,main="PCA eigenvalues", col=heat.colors(50))

pca1

# Plotting PC1 and PC2
s.label(pca1$li)
title("PCA of summer flounder dataset\naxes 1-2")
add.scatter.eig(pca1$eig, 3,1,2)

s.class(pca1$li, pop(larvs.freqs2))
title("PCA of summer flounder dataset\naxes 1-2")
add.scatter.eig(pca1$eig, 3,1,2)

eig_percent <- round((pca1$eig/(sum(pca1$eig)))*100,2)
eig_percent [1:3]

# Break PCA down by region, plotting PC1 vs PC2
plot(pca1$li[c(1:6, 15:28, 40:48, 265:271, 290:291),1], pca1$li[c(1:6, 15:28, 40:48, 265:271, 290:291),2], col = "steelblue3", pch = 6, xlab = "PC1 (9.67%)", ylab = "PC2 (9.32%)", xlim = c(-8,6), ylim = c(-7,6)) # 38 fish (mid period/north)
points(pca1$li[c(7:14, 29:39, 49:58, 272:289, 292:293),1], pca1$li[c(7:14, 29:39, 49:58, 272:289, 292:293),2], col = "lightcoral", pch = 6) # 49 fish (mid period/south)
points(pca1$li[c(59:83, 95:130, 143:154, 178:204),1], pca1$li[c(59:83, 95:130, 143:154, 178:204),2], col = "dodgerblue4", pch = 0) # 100 fish (late period/north)
points(pca1$li[c(84:94, 131:142, 155:177),1], pca1$li[c(84:94, 131:142, 155:177),2], col = "red", pch = 0) # 46 fish (late period/south)
points(pca1$li[c(224:225, 240:244, 264),1], pca1$li[c(224:225, 240:244, 264),2], col = "blue") # 8 fish (early period/north)
points(pca1$li[c(205:223, 226:239, 245:263),1], pca1$li[c(205:223, 226:239, 245:263),2], col = "hotpink") # 52 (early period/south)

legend("topleft",
       legend=c("North", "South"),
       pch=c(1, 1),
       col=c("blue", "red"))

# Break PCA down by region, plotting PC1 vs PC3
plot(pca1$li[c(1:6, 15:28, 40:48, 265:271, 290:291),1], pca1$li[c(1:6, 15:28, 40:48, 265:271, 290:291),3], col = "steelblue3", pch = 6, xlab = "PC1 (9.67%)", ylab = "PC2 (9.32%)", xlim = c(-8,6), ylim = c(-7,6)) # 38 fish (mid period/north)
points(pca1$li[c(7:14, 29:39, 49:58, 272:289, 292:293),1], pca1$li[c(7:14, 29:39, 49:58, 272:289, 292:293),3], col = "lightcoral", pch = 6) # 49 fish (mid period/south)
points(pca1$li[c(59:83, 95:130, 143:154, 178:204),1], pca1$li[c(59:83, 95:130, 143:154, 178:204),3], col = "dodgerblue4", pch = 0) # 100 fish (late period/north)
points(pca1$li[c(84:94, 131:142, 155:177),1], pca1$li[c(84:94, 131:142, 155:177),3], col = "red", pch = 0) # 46 fish (late period/south)
points(pca1$li[c(224:225, 240:244, 264),1], pca1$li[c(224:225, 240:244, 264),3], col = "blue") # 8 fish (early period/north)
points(pca1$li[c(205:223, 226:239, 245:263),1], pca1$li[c(205:223, 226:239, 245:263),3], col = "hotpink") # 52 (early period/south)

legend("bottomleft",
       legend=c("North", "South"),
       pch=c(1, 1),
       col=c("blue", "red"))

#### Subset larvae allele frequencies to those that match with adults ####
larvs.freqs2.sub <- larvs.freqs2@tab[,c(1,2,3,4,5,6,11,12,13,14,15,16,17,18,24,23)] # SNP_1639 = index 25 & 26

# Do they match??
adult_contigs2
# colnames(regional.outs.freqs2)
# colnames(larvs.freqs2.sub)
rbind(colnames(regional.outs.freqs), colnames(larvs.freqs2.sub))

#### Allele frequency tables between larvae & adults match ####
# First, let's look at histograms of allele frequencies. I'm guessing they should look approximately similar
dim(as.data.frame(outs@tab[,c(1,2,3,4,5,6,11,12,13,14,15,16,17,18,23,24)])) # 241 x 16
dim(larvs.freqs2.sub) # 293 x 16

adult_contigs2
rbind(colnames(as.data.frame(outs@tab[,c(1,2,3,4,5,6,11,12,13,14,15,16,17,18,23,24)])), colnames(larvs.freqs2.sub)) # double check snp names, again

hist(as.data.frame(outs@tab)[,1])
hist(larvs.freqs2.sub[,1])
hist(as.data.frame(outs@tab)[,3])
hist(larvs.freqs2.sub[,3])
hist(as.data.frame(outs@tab)[,5])
hist(larvs.freqs2.sub[,5])
hist(as.data.frame(outs@tab)[,11])
hist(larvs.freqs2.sub[,7])
hist(as.data.frame(outs@tab)[,13])
hist(larvs.freqs2.sub[,9])
hist(as.data.frame(outs@tab)[,15])
hist(larvs.freqs2.sub[,11])
hist(as.data.frame(outs@tab)[,17])
hist(larvs.freqs2.sub[,13])
hist(as.data.frame(outs@tab)[,23])
hist(larvs.freqs2.sub[,15])

colnames(larvs.freqs2.sub) <- colnames(as.data.frame(outs@tab[,c(1,2,3,4,5,6,11,12,13,14,15,16,17,18,23,24)])) 

#### Calculating genotype likelihoods ####
regional.outs.freqs

odds <- seq(1,15,2) # odd indicies to keep
regional.outs.freqs.odds <- regional.outs.freqs[,odds]

larvs.freqs2.sub.odds <- larvs.freqs2.sub[,odds] # only odds indicies in the larval dataset
larvs.freqs2.sub.odds[is.na(larvs.freqs2.sub.odds)] <- 9 # replace NA's with 9's to make the ifelse statements easier

colnames(larvs.freqs2.sub.odds) == colnames(regional.outs.freqs.odds) # column names match?

# For loop to loop through each locus & multiply by the adult allele frequency, and then do this for all 293 larvae. NA's/9's get coded as 1's so they don't make a difference when each row's product is taken
north.likelihoods <- data.frame()
south.likelihoods <- data.frame()

for (j in 1:length(rownames(larvs.freqs2.sub.odds))){

for (i in 1:length(colnames(larvs.freqs2.sub.odds))){
  if(larvs.freqs2.sub.odds[j,i] == 2) {
    north.likelihoods[j,i] <- regional.outs.freqs.odds[1,i]^2
  } else if (larvs.freqs2.sub.odds[j,i] == 1) {
      north.likelihoods[j,i] <- 2*(regional.outs.freqs.odds[1,i] * (1-regional.outs.freqs.odds[1,i]))
  } else if (larvs.freqs2.sub.odds[j,i] == 0) {
      north.likelihoods[j,i] <- ( 1-regional.outs.freqs.odds[1,i])^2 
     } else {
        north.likelihoods[j,i] <- 1
    }
  }

for (i in 1:length(colnames(larvs.freqs2.sub.odds))){
  if(larvs.freqs2.sub.odds[j,i] == 2){
    south.likelihoods[j,i] <- regional.outs.freqs.odds[2,i]^2
  } else if (larvs.freqs2.sub.odds[j,i] == 1) {
    south.likelihoods[j,i] <- 2*(regional.outs.freqs.odds[2,i] * (1-regional.outs.freqs.odds[2,i]))
  } else  if (larvs.freqs2.sub.odds[j,i] == 0) {
    south.likelihoods[j,i] <- (1-regional.outs.freqs.odds[2,i])^2
  } else {
    south.likelihoods[j,i] <- 1
  }
}
}

for (i in 1:length(a)){
if(a[i] == 2){
  print('homo')
} else if (a[i] == 1){
  print('hetero')
} else if (a[i] == 0){
  print('homo')
} else {
  print('1')
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
(0.7187500^2)*(2*0.7534722*0.2465278)*(0.8854167^2)*(2*0.7881944*0.2118056)*(0.9513889^2)*(2*0.5937500*0.4062500)*(0.9027778^2)*(0.2256944^2) #north fish1
(0.8402062^2)*(2*0.8578947*(1-0.8578947))*(0.9587629^2)*(2*0.6958763*(1-0.6958763))*(0.8906250^2)*(2*0.6958763*(1-0.6958763))*(0.9639175^2)*((1-0.8842105)^2)#south

(0.7187500^2)*(2*0.7534722*(1-0.7534722))*(0.8854167^2)*(0.7881944^2)*(0.9513889^2)*(0.5937500^2)*(0.9027778^2)*((1-0.7743056)^2) #north fish2
(0.8402062^2)*(2*0.8578947*(1-0.8578947))*(0.9587629^2)*(0.6958763^2)*(0.8906250^2)*(0.6958763^2)*(0.9639175^2)*((1-0.8842105)^2) #south

(2*0.7187500*(1-0.7187500))*(0.7534722^2)*(0.8854167^2)*(2*0.7881944*(1-0.7881944))*(0.9513889^2)*(0.9027778^2)*((1-0.7743056)^2) # north fish 29
(2*0.8402062*(1-0.8402062))*(0.8578947^2)*(0.9587629^2)*(2*0.6958763*(1-0.6958763))*(0.8906250^2)*(0.9639175^2)*((1-0.8842105)^2) #south

(0.7187500^2)*(0.7534722^2)*(0.8854167^2)*(0.7881944^2)*(0.9513889^2)*(0.5937500^2)*(0.9027778^2) # north fish 30
(0.8402062^2)*(0.8578947^2)*(0.9587629^2)*(0.6958763^2)*(0.8906250^2)*(0.6958763^2)*(0.9639175^2)

# Determine if each fish has a higher likelihood of coming from the north or the south
regional.vector <- cbind(north.vector, south.vector)

assignment <- vector()
for (m in 1:length(regional.vector[,1])){
  if (regional.vector[m,1] > regional.vector[m,2]){
    assignment[m] <- "north"
  } else {
    assignment[m] <- "south"
  }
}
