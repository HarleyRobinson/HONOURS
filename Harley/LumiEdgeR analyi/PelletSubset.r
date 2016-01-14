#To make life easier, load just a subset of the data. 
rawdata <- read.csv("miRNA_raw_counts_HEK_Pellet.csv")
# create a DGEList object. In this case, raw data is in columns 2 to 55, and unique identifiers are in 1 (miRNA name)
y <- DGEList(counts=rawdata[,2:13], genes=rawdata[,1])
# Filtering and Normalization
# compute the effective library size by using TMM normalization
keep <- rowSums(cpm(y)>10) >=2
y <- y[keep,]
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y)
#make QC plots
pdf(file="PelletOnlyDend.pdf", paper="a4r")
d <-cpm(y, normalized.lib.sizes=TRUE)
d<-t(d)
dist<-dist(d)
hc<-hclust(dist)
plot(hc)
dev.off()

#export normalized data into tab delimited text for exploration in genespring
cpm<- cpm(y, normalized.lib.sizes=TRUE)
rownames(cpm)<-y$genes[,1]
write.table(cpm, "cpm_miRNA_PC3_subset.txt", sep="\t")

#plot the MDS
tiff(file= "PelletOnlyMDS.tiff, width=480, height=480")
plotMDS(y, method="bcv")
dev.off()
#make a scatterplot
cpm<-replace(cpm, cpm==0, 0.01)
min<-ExpressionSet(assayData=cpm)
pairs(min, smoothScatter=TRUE)

# annotate groups
Groups <-factor(c(rep(c("GFP", "CAV1", "CAV2","CAV3"),3)))
design <- model.matrix(~Groups)
rownames(design) <- colnames(y)

#Perform the statistical analysis
# estimate the overall level of biological variability
y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
# estimate genewise and tagwise dispersions
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
# fit the linear model
fit <- glmFit(y, design)

# Identify differentially expressed genes first with individual comparisons
# define groups
y$samples$group <- c(rep(c("GFP", "CAV1", "CAV2","CAV3"),3))


# calculate differential expression
GFP_vs_CAV1 <- exactTest(y, pair=c("GFP","CAV1"))
# summarize results, -1 means down regulated, 1 means up regulated. 
# The following results mean that there was 1 down-regulated gene, 1 up-regulated gene, and 364 
# non-changing genes
summary(de <- decideTestsDGE(GFP_vs_CAV1, p=0.05, adjust="BH"))
topTags(GFP_vs_CAV1, n=10)

# calculate differential expression
GFP_vs_CAV2 <- exactTest(y, pair=c("GFP","CAV2"))
# summarize results, -1 means down regulated, 1 means up regulated. 
# The following results mean that there was 10 down-regulated genes, 5 up-regulated genes, and 351 
# non-changing genes
summary(de <- decideTestsDGE(GFP_vs_CAV2, p=0.05, adjust="BH"))
topTags(GFP_vs_CAV2, n=10)

# calculate differential expression
GFP_vs_CAV3 <- exactTest(y, pair=c("GFP","CAV3"))
# summarize results, -1 means down regulated, 1 means up regulated. 
# The following results mean that there was 1 down-regulated gene, 0 up-regulated genes, and 364 
# non-changing genes
summary(de <- decideTestsDGE(GFP_vs_CAV3, p=0.05, adjust="BH"))
topTags(GFP_vs_CAV3, n=10)
