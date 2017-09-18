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

#### PCA ####
sum(is.na(full$tab)) #41594 
X <- scaleGen(full, NA.method = "mean")
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

s.class(pca1$li, pop(full))
title("PCA of summer flounder dataset\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)

eig_percent <- round((pca1$eig/(sum(pca1$eig)))*100,2)
eig_percent [1:3]

# To make a nice PCA using all fish
png(file="~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/results/larvae_vs_adults_pca.png", width=8, height=7, res=300, units="in")

par(
  mar=c(5, 4, 4, 4), # panel magin size in "line number" units
  mgp=c(3, 1, 0), # default is c(3,1,0); line number for axis label, tick label, axis
  tcl=-0.5, # size of tick marks as distance INTO figure (negative means pointing outward)
  cex=1, # character expansion factor; keep as 1; if you have a many-panel figure, they start changing the default!
  ps=14, # point size, which is the font size
  bg=NA,
  bty = 'n'
)

col <- deepseasun(4)
# col <- c("#EE0000", "#F98F00", "#AA8F4F", rgb(1,1,1,0))
s.class(pca1$li, pop(full), xax=1,yax=2, col = transp(col,0.7), axesell=TRUE, cellipse=1.5, cstar=1,cpoint=1.75, grid=FALSE, addaxes = FALSE, xlim = c(-60,10), ylim = c(-40,40), clabel = 0)
axis(1, at=seq(-60,10, by=10), labels=seq(-60,10, by= 10), line = 1.5)
axis(2, at=seq(-30,30, by = 10), labels=seq(-30,30, by= 10), line = 0, las = 2)
mtext("PC1 (1.05%)", side = 1, line = 4)
mtext("PC2 (0.58%)", side = 2, line = 2.5)

legend(-60, 15,
       legend=c("1989-1993 Larvae (n = 60)", "1998-2002 Larvae (n = 87)", "2008-2012 Larvae (n = 146)", "2013-2014 Adults (n = 235)"),
       pch=c(19, 19, 19, 19),
       col=col,
       bty = "n",
       y.intersp = 1)

dev.off()

#### DAPC - Still seems kind of hokey to me ####
# K-means
grp <- find.clusters(full, max.n.clust = 10) # no clear asymptote, so keep all PCs, 2
table.value(table(pop(full), grp$grp), col.lab=paste("grp", 1:2))
dapc1 <- dapc(full, grp$grp)
dapc1

scatter(dapc1)
test <- optim.a.score(dapc1) #80,1 --> 14

# Choosing the number of PCs to retain
dapc2 <- dapc(full, n.da = 1, n.pca = 14)
scatter(dapc2)

#### AMOVA ####
full_pop <- as.data.frame(full@pop)
full_dist <- dist(full)
full_amova <- pegas::amova(full_dist ~ full@pop, data = full_pop, nperm = 1000) # different but which ones??
full_amova

#### Calculate mean allele frequencies for each 'population' ####
pops <- as.vector(full@pop)
full.df <- cbind(data.frame(scaleGen(full, center = FALSE, scale = FALSE, NA.method = "mean")), pops) # this goes from counts to frequencies
full.df.larvs <- full.df[which(full.df[,3832]== '3'),]
dim(full.df.larvs) # 146 x 3832
full.df.adults <- full.df[which(full.df[,3832]== '4'),]
dim(full.df.adults) # 235 x 3832
full.df2 <- rbind(full.df.larvs, full.df.adults)
dim(full.df2)

se <- function(x) sqrt(var(x,na.rm=TRUE)/length(na.omit(x))) # Function to calculate SE

full.df3 <- full.df2[, -3832] # removing the pop columns 
allelefreq.full.se <- aggregate(full.df3, list(full.df2$pops), se) #take se by period
allelefreq.full.mean <- aggregate(full.df3, list(full.df2$pops), mean, na.rm = TRUE) #take mean by period
allelefreq.full.mean.t <- t(allelefreq.full.mean[,-1]) # get rid of the categories so there are fewer problems down the line (num turning to chr), then transpose
dim(allelefreq.full.mean.t)

# Allele frequency differences between 1008-2012 larvae and 2013-2014 adults
diffs <- allelefreq.full.mean.t[,1] - allelefreq.full.mean.t[,2]
hist(diffs)
summary(diffs)

diff_outliers <- which(as.vector(diffs) > (sd(diffs, na.rm = TRUE))*3)
diff_outliers_df <- diffs[diff_outliers]

# Plotting boxplots of allele frequency of outliers in larvae vs adults
par(mfrow = c(3,4))
for (i in diff_outliers){
  lmts <- range(0,1)
  boxplot(full.df.larvs[,i], full.df.adults[,i], ylim = lmts, main = paste(colnames(full.df.adults)[i]))
  points(allelefreq.full.mean.t[i,], col = 'red', pch = 19)
  axis(1, at=1:2, labels = c('2008-2012 Larvae \n(n = 146)', '20013-2014 Adults \n(n = 235)'))
  # boxplot(full.df.adults[,i], ylim = lmts, add = TRUE)
}
