# Install the needed libraries and load them
if(!requireNamespace("tweeDEseqCountData", quietly=TRUE)) BiocManager::install("tweeDEseqCountData")
if(!requireNamespace("DESeq2", quietly=TRUE)) BiocManager::install("DESeq2")
library(tweeDEseqCountData)
library(DESeq2)

# Load data from package tweeDEseqCountData
data(pickrell)
data(annotEnsembl63)
data(genderGenes)
ls()

# chrY genes should be up-regulated in male as compared to female samples
length(msYgenes)
msYgenes[1:4]

# chrX genes that escape inactivation, hence many are up-regulated in female as compared to male samples
length(XiEgenes)
XiEgenes[1:4]

# Ensembl gene annotations
head(annotEnsembl63)
annotEnsembl63[11000:11010,]

# The assay data in an expression set object
pickrell.eset
PickrellDataSet <- exprs(pickrell.eset)
dim(PickrellDataSet)
PickrellDataSet[1:6,1:4]
# gender labels available in the expression set object
gender <- pickrell.eset$gender
table(gender)

# Many genes are unexpressed with zero counts, and could be removed
# Pick genes with <50 zero counts under the 69 samples of the dataset
PickrellDataSet2 <- PickrellDataSet[which(rowSums(PickrellDataSet==0) < 50),]
dim(PickrellDataSet2)

# Build the DESeqDataSet object
coldata <- data.frame(gender, row.names=colnames(PickrellDataSet2))
des <- DESeqDataSetFromMatrix(countData=PickrellDataSet2, colData=coldata, design=~gender)
class(des)
des

# Add meta data columns to the DESeqDataSet object
featureData <- data.frame(annotEnsembl63[rownames(PickrellDataSet2), c("Symbol","Chr")])
mcols(des) <- DataFrame(mcols(des), featureData)
head(mcols(des))

# Perform Wald test
des <- DESeq(des, test="Wald")
res <- results(des, format="DataFrame")
head(res, n=10)

# Construct a data frame of results and sort it by adjusted p-value
df <- as.data.frame(cbind(mcols(des)[, c("Symbol","Chr")], as.data.frame(res)))
df <- df[sort(df$padj, decreasing=FALSE, index.return=TRUE)$ix,]
# Save results in a csv file that can be viewed in Excel
write.csv(df, file="DESeq_results.csv")

# Identify significant differentially expressed (DE) genes
DE.ind <- which(df$padj < 0.05 & abs(df$log2FoldChange) > 1)
df[DE.ind,]

# Show msY genes that are expressed in the dataset
intersect(rownames(PickrellDataSet2), msYgenes)

# Apply variance stabilization transform and generate PCA plot
vsd <- vst(des, blind=FALSE)
plotPCA(vsd, intgroup="gender")

# Generate PCA plot considering only the significant DE genes
gg <- rownames(df)[DE.ind]
plotPCA(vsd[gg,], intgroup="gender")
