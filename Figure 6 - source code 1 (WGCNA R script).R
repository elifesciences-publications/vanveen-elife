getwd()
source("http://bioconductor.org/biocLite.R") 
biocLite("BiocUpgrade")
biocLite(c("GO.db", "preprocessCore", "impute")) 
y
install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival"))
orgCodes = c("Hs", "Mm", "Rn", "Pf", "Sc", "Dm", "Bt", "Ce", "Cf", "Dr", "Gg"); 
orgExtensions = c(rep(".eg", 4), ".sgd", rep(".eg", 6)); 
packageNames = paste("org.", orgCodes, orgExtensions, ".db", sep=""); 

biocLite(c("GO.db", "KEGG.db", "topGO", packageNames, "hgu133a.db", "hgu95av2.db", "annotate", "hgu133plus2.db", "SNPlocs.Hsapiens.dbSNP.20100427", "minet", "OrderedList"))
biocLite(c("AnnotationDbi", "impute", "GO.db", "preprocessCore")) 
install.packages("WGCNA")
library(AnnotationDbi)
library(WGCNA)
getwd()





options(stringsAsFactors = FALSE)
lung3 = read.delim("~/Desktop/GSEA/Galaxy175-[Normalized_counts_file_on_data_145,_data_137,_and_others].tabular")
dim(lung3)
datExpr0 = as.data.frame(t(lung3[, -c(1)]))
names(datExpr0) = lung3$X
rownames(datExpr0) = names(lung3)[-c(1)]
means = colMeans(datExpr0, na.rm = FALSE)
keep <- means > 40
datExpr0 <- datExpr0[,keep]


gsg = goodSamplesGenes(datExpr0, verbose = 3, minNSamples = 13)
goodSamplesGenes(datExpr0, verbose = 3)

head(lung3)

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

dim(datExpr0)
head(datExpr0)
str(datExpr0)


sampleTree = hclust(dist(datExpr0), method = "average")
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dim(sampleTree)

str(sampleTree)

nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

net = blockwiseModules(datExpr0, power = 3, maxBlockSize = 14207,
                       TOMType = "signed", minModuleSize = 50,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE, TOMDenom = "min", 
                       saveTOMFileBase = "lung3", corType = "pearson",
                       verbose = 3)     



table(net$colors)

netgroups = as.data.frame(net$colors)
netgenes = as.data.frame(names(datExpr0))
write.csv(netgroups, file = "~/Desktop/temp/wgcnagroupsP3t.csv")
write.csv(netgenes, file = "~/Desktop/temp/wgcnagenesP3t.csv")
print(mergedColors)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)


#dendrogram
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)




moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]


datExpr0[moduleLabels=="0"]

cluster <- t(datExpr0[moduleLabels=="1"])
str(cluster)
str(net)
head(moduleLabels)


#modGenes = (moduleColors==module)

#str(net)
#netf = as.data.frame(net$colors)
#write.csv(netf, file = "netf.csv")
#exprf = as.data.frame(datExpr0$names)
#print(net$blockGenes)


str(datExpr0)


write.csv(datExpr0, file = "net.csv")

load(file = "~/Dropbox/JEV/R/lung3-block.1.RData")
dissTOM = 1-(as.matrix(TOM))


# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
# diag crashes computer (not after filtering)
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(15, 15)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot")
