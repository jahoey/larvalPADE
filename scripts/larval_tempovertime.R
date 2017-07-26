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

########################################################################################
# Reading in NC temp data from 1994-2017
setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Larvae")

library(readxl)

read_excel_allsheets <- function(filename) {
  sheets <- readxl::excel_sheets(filename)
  x <-    lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  names(x) <- sheets
  x
}

nc_data <- read_excel_allsheets("NOAA BFT Temperature Data 19940715 to 20170720.xlsx") # creates a list of dataframes

# nc_df <- ldply(nc_data, data.frame)
nc_df <- as.data.frame(rbind(nc_data$`July 15 1994 to May 2006`, nc_data$`June 2006 to present`))
dim(nc_df) # 1947573 x 3

# The date and time are in a single column. Splitting them up.
nc_df$Date <- as.factor(format(nc_df$`DATE TIME`,"%m/%d/%Y",tz = "GMT"))

# Taking daily mean and ordering by date
nc_df_mean <- aggregate(nc_df$WATERTEMP, list(Date = nc_df[,"Date"]), mean, na.rm = TRUE) # length(unique(nc_df$Date)) = 8130
nc_df_mean_order <- nc_df_mean[order(as.Date(nc_df_mean$Date,format="%m/%d/%Y")),]
dim(nc_df_mean_order) # 8130 x 2

# Plot daily mean over time
nc_df_mean_order$Index <- 1:length(nc_df_mean_order$Date)
plot(nc_df_mean_order$x ~ nc_df_mean_order$Index)
lm_daily <- lm((nc_df_mean_order$x ~ nc_df_mean_order$Index))
abline(lm_daily)

########################################
# Getting years in a column by themselves
nc_df$Year <- as.factor(format(nc_df$`DATE TIME`,"%Y",tz = "GMT"))

# Taking yearly mean and ordering by year
nc_df_mean_year <- aggregate(nc_df$WATERTEMP, list(Year = nc_df[,"Year"]), mean, na.rm = TRUE) # length(unique(nc_df$Year)) = 24
dim(nc_df_mean_year) # 24 x 2

# Plot yearly mean over time
nc_df_mean_year$Index <- 1:length(nc_df_mean_year$Year)
plot(nc_df_mean_year$x ~ nc_df_mean_year$Year)
lm_year <- lm((nc_df_mean_year$x ~ nc_df_mean_year$Index))
abline(lm_year)

###################################
# Calculating temperature anomolies
mean(nc_df_mean_year$x)
mean(nc_df_mean$x, na.rm = TRUE)
mean(nc_df$WATERTEMP, na.rm = TRUE)

nc_df_mean_year$Anom <- mean(nc_df$WATERTEMP, na.rm = TRUE) - nc_df_mean_year$x

bar <- barplot(nc_df_mean_year$Anom, ylab = 'Temperature anomoly °C', axes = FALSE, ylim = c(-10, 6))
axis(1, at=bar, labels=nc_df_mean_year$Year)
axis(2, at=seq(-10,6, by = 2), labels=seq(-10,6, by= 2), las = 2)

