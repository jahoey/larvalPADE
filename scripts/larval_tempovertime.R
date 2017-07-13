# Reading in the environmental data and allele count data
setwd("/Users/jenniferhoey/Documents/Graduate School/Rutgers/Summer Flounder/Larvae")
larvae_database <- read.csv("Larvae Sampling Database.csv", header=TRUE)

setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis/")
larvae_allelecounts <- read.table("larvs_allelecounts.txt")
dim(larvae_allelecounts) # 293 x 3764

# Fixing the IDs so that I can match fish with allele counts and environmental data
# The names are all messed up because of the bioinformatics process...
names1 <- do.call(rbind, strsplit(as.character(larvae_database$ID..), '_'))
names2 <- do.call(rbind, strsplit(as.character(names1[,1]), 'E'))
correctnames <- paste(names2[,1], names2[,2], sep = 'E_')
correctnames2 <- paste(correctnames, names1[,2], sep = '')

# Now put back the fixed names in the larval database and the split names into the allele frequency database
larvae_database$ID.. <- correctnames2

# Need to modify names in allele count data (genind object) so that I can match them up with names in the larval database
freq_names <- rownames(larvae_allelecounts)
freq_names_split <- do.call(rbind, strsplit(as.character(freq_names), 'L'))
rownames(larvae_allelecounts) <- freq_names_split[,1] # replaces rownames in dataframe allele frequency counts with modified format

# Subset the larval database to only larvae that are staying in the analysis, and order the names so that they match with the allele frequency matrix
larvae_database_sub <- larvae_database[larvae_database$ID.. %in% freq_names_split[,1],]
ordered_larvae_database_sub <- larvae_database_sub[order(larvae_database_sub$ID..),]

# Do the names match?
rownames(larvae_allelecounts) == ordered_larvae_database_sub$ID..

# Breaking the dataset into NJ and NC and plotting
# New Jersey
nj <- ordered_larvae_database_sub[ordered_larvae_database_sub$Place %in% c("Little Egg Inlet, NJ"),]
nj$mdy_date <- paste(nj$Month, nj$Date, nj$Year, sep = '/')
nj$Sampling.Date <- as.Date(nj$mdy_date, "%m/%d/%Y")

plot(nj$Start.Temp ~ nj$Sampling.Date, xlab = "Date", ylab = 'Temperature', xaxt = 'n') # plot by date
axis(1, nj$Sampling.Date, format(nj$Sampling.Date, "%m/%d/%Y"), cex.axis = .7)
nj_lm1 <- lm(nj$Start.Temp ~ nj$Sampling.Date)
abline(nj_lm1) # decreasing slope

nj_tempmean <- aggregate(nj$Start.Temp, list(nj$Year), mean, na.rm = TRUE) #take mean by year
plot(nj_tempmean, xlab = "Year", ylab = "Mean Temperature") # plot by year
nj_lm2 <- lm(nj_tempmean$x ~ nj_tempmean$Group.1)
abline(nj_lm2) # decreasing slope

# North Carolina: about half the data is missing
nc <- ordered_larvae_database_sub[ordered_larvae_database_sub$Place %in% c("Beaufort, NC"),]

nc_tempmean <- aggregate(nc$Start.Temp, list(nc$Year), mean, na.rm = TRUE) #take mean by year
plot(nc_tempmean, xlab = "Year", ylab = "Mean Temperature")

