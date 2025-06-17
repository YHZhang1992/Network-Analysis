# Install WGCNA from CRAN or Bioconductor
install.packages("WGCNA")  # Or:
# BiocManager::install("WGCNA")

# Load libraries
library(WGCNA)
options(stringsAsFactors = FALSE)

# Read expression matrix (genes x samples)
datExpr0 <- read.csv("expression_matrix.csv", row.names = 1)
datExpr0 <- as.data.frame(t(datExpr0))  # WGCNA requires samples as rows

# Check for good samples and genes
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree <- hclust(dist(datExpr0), method = "average")
plot(sampleTree, main = "Sample Clustering", sub = "", xlab = "")

powers <- c(1:20)
sft <- pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

# Plot scale-free fit
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit", type = "n")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, col = "red")
abline(h = 0.90, col = "blue")  # recommended threshold

softPower <- 6  # choose based on previous step
adjacency <- adjacency(datExpr0, power = softPower)

# Turn adjacency into topological overlap matrix (TOM)
TOM <- TOMsimilarity(adjacency)
dissTOM <- 1 - TOM

# Hierarchical clustering
geneTree <- hclust(as.dist(dissTOM), method = "average")
plot(geneTree, main = "Gene Clustering on TOM-based Dissimilarity")

# Module identification using dynamic tree cut
minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)

# Convert numeric labels into colors
moduleColors <- labels2colors(dynamicMods)
plotDendroAndColors(geneTree, moduleColors, "Dynamic Tree Cut")

MEList <- moduleEigengenes(datExpr0, colors = moduleColors)
MEs <- MEList$eigengenes

# Cluster module eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")
plot(METree, main = "Clustering of Module Eigengenes")

# Merge similar modules
mergeCutHeight <- 0.25
merge <- mergeCloseModules(datExpr0, moduleColors, cutHeight = mergeCutHeight, verbose = 3)

mergedColors <- merge$colors
mergedMEs <- merge$newMEs

# Update module assignments
moduleColors <- mergedColors
MEs <- mergedMEs

# Load trait data
traitData <- read.csv("trait_data.csv", row.names = 1)
traitData <- traitData[match(rownames(datExpr0), rownames(traitData)), ]

# Correlate module eigengenes with traits
moduleTraitCor <- cor(MEs, traitData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr0))

# Heatmap
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traitData),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = signif(moduleTraitCor, 2),
               setStdMargins = FALSE,
               cex.text = 0.5,
               main = "Module-Trait Relationships")

# Example: identify genes in the "blue" module
module <- "blue"
moduleGenes <- moduleColors == module

# Intramodular connectivity
geneModuleMembership <- cor(datExpr0, MEs, use = "p")
MMblue <- geneModuleMembership[moduleGenes, module]
topHubGenes <- names(sort(MMblue, decreasing = TRUE))[1:10]
print(topHubGenes)

# Export edges from TOM
modProbes <- names(datExpr0)[moduleGenes]
modTOM <- TOM[moduleGenes, moduleGenes]
dimnames(modTOM) <- list(modProbes, modProbes)

# Export for Cytoscape
cyt <- exportNetworkToCytoscape(modTOM,
                                edgeFile = "CytoscapeEdges.txt",
                                nodeFile = "CytoscapeNodes.txt",
                                weighted = TRUE,
                                threshold = 0.02,
                                nodeNames = modProbes,
                                nodeAttr = moduleColors[moduleGenes])

