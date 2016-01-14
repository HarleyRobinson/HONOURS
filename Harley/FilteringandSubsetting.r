library(lumi)
library(edgeR)
rawdata <- read.delim("miRNA_raw_counts_PC3_pellet_G_1_2_3_FILTERED.txt", header= TRUE, row.names=1)
#rawdata<- read.csv("miRNA_raw_counts_Correct_Name_ReOrdered.csv", header= TRUE, row.names= 1)
#rawdata<- rawdata[c(31:42)]

#Below is only PC3 data
y <- DGEList(counts=rawdata[,1:12], genes=rownames(rawdata))
# Filtering and Normalization
# compute the effective library size by using TMM normalization
keep <- rowSums(cpm(y)>10) >=2
y <- y[keep,]
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y)
#export normalized data into tab delimited text for exploration in genespring
cpm<- cpm(y, lib.sizes=TRUE)
rownames(cpm)<-y$genes[,1]
write.table(cpm, "test.txt", sep="\t")

##To Filter and subset data 


