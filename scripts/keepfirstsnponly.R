setwd("~/Documents/Graduate School/Rutgers/Summer Flounder/Analysis/full_PADE_analysis")

data <- read.table("SNP.DP3g95maf05lm90.FIL.recode.vcf", skip = 61, sep = "\t", header = TRUE) # Need to delete the '#' in front of CHROM in order for the header to read in properly
dim(data) #3827 x 537

# Keeping only the first SNP at each contig
vcf_firstsnps <- data[!duplicated(data[,1]),] # these are the first snps at each contig
dim(vcf_firstsnps) # 1904 x 537

write.table(vcf_firstsnps, "SNP.DP3g95maf05lm90.FIL.recode.firstsnp.vcf", sep="\t", col.names = TRUE, row.names = FALSE)
# Manually removed all the " and replaced with nothing

# Can convert file types then read str file into R and analyze

# SNP.DP3g95maf05lm75.FIL.recode.vcf 3827 x 537 (firstsnp = 1904 x 537)
# SNP.DP3g95maf05lm90.FIL.recode.vcf 1967 x 563 (firstsnp = 1109 x 563)