# larvalPADE
Data files, code and analyses involving summer flounder larvae from 1989-2012.
#### scripts
- **15outlier_PCA.R** script examines 15 spatial outliers in adult PADE w/ PCA, calculates regional allele frequencies, identifies these SNPs in larvae, and calculates genotype likelihoods
- **fallvswinter.R** exploratory code examining allele frequency differences between larvae ingressing in the fall and the winter
- **keepfirstsnponly.R** script reads in a vcf, removes SNPs existing on the same contig & writes a new file
- **larvaevadults.R** exploratory code examining allele frequency differences between adult and larval PADE
- **larvaevadults_fst.R** calculates FST between adult and larval populations caught over time
- **larval_sampling_simulations.R** creates and plots a null model that examines allele frequency change over time that depends on sampling size
- **larval_tempovertime.R** preliminary analysis of temperature change over time using RUMFS bridgenetting & NOAA Beaufort datasets
- **larvalstructure.R** PCA broken down by region
- **structurethrutime.R** exploratory code for larval PADE, including PCA broken down by time period, summary statistics, DAPC, AMOVA and allele frequency change over time
