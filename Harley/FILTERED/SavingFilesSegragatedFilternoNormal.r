library(edgeR)
library(lumi)
# read in the raw data
rawdata <- read.csv("miRNA_raw_counts_Correct_Name_ReOrdered.csv", header= TRUE, row.names=1)
# create a DGEList object. In this case, raw data is in columns 2 to 55, and unique identifiers are in 1 (miRNA name)
# Filtering and Normalization
# compute the effective library size by using TMM normalization
keep <- rowSums(rawdata>10) >=2
y <- rawdata[keep,]
#write.csv(y, "miRNA_raw_counts_ALL_filtered.csv")

PC3data<- y[c(43:54)]
HEKdata<- y[c(16:30)]
write.csv(PC3data, "miRNA_raw_counts_PC3_filtered_exosome.csv")
write.csv(HEKdata, "miRNA_raw_counts_HEK_filtered_exosome.csv")

recombind<- read.csv("miRNA_raw_counts_PC3_filtered_exosome (2).csv", header=TRUE, row.names=1)
recomb<- read.csv("miRNA_raw_counts_HEK_filtered_exosome (2).csv", header= TRUE, row.names=1)
recombind<- cbind(recombind, recomb)
keep <- rowSums(recombind>10) >=2
y <- recombind[keep,]
PC3data<- y[c(1:12)]
HEKdata<- y[c(13:27)]
write.csv(PC3data, "miRNA_raw_counts_PC3_filtered_exosome.csv")
write.csv(HEKdata, "miRNA_raw_counts_HEK_filtered_exosome.csv")
