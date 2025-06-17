# Install from Bioconductor (if not yet installed)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("GNET2")

# Load required libraries
library(GNET2)
library(data.table)
library(igraph)

# Load expression data
expr <- fread("expression_matrix.csv", data.table = FALSE)
rownames(expr) <- expr[,1]
expr <- expr[,-1]
expr <- as.matrix(expr)

# Optional: filter top variable genes
top_var_genes <- head(order(apply(expr, 1, var), decreasing = TRUE), 2000)
expr_filtered <- expr[top_var_genes, ]

# Read TRRUST TF-target pairs
trrust <- fread("trrust_rawdata.human.tsv", header = FALSE)  # for human
colnames(trrust) <- c("TF", "Target", "Type", "PMID")

# Filter for TFs in your expression data
trrust_filtered <- trrust[trrust$TF %in% rownames(expr_filtered) & trrust$Target %in% rownames(expr_filtered), ]

gnet_obj <- makeGNET(expr_filtered, tf_set = unique(trrust_filtered$TF))

# Supply prior network (TF-target edges from TRRUST)
tf_target_matrix <- as.matrix(trrust_filtered[, c("TF", "Target")])
gnet_obj <- set_prior(gnet_obj, tf_target_matrix)

# Run the GNET2 iterative inference
gnet_obj <- learn_network(gnet_obj, max_iter = 20, e_cutoff = 0.01)

# Extract modules and GRN
modules <- gnet_obj$module
network <- gnet_obj$network

# Convert to igraph object
g <- graph_from_data_frame(network, directed = TRUE)

# Plot
plot(g, vertex.label.cex = 0.6, edge.arrow.size = 0.2, layout = layout_with_fr)

# View module assignments
head(modules)

# Top regulators by number of targets
reg_counts <- table(network$reg)
head(sort(reg_counts, decreasing = TRUE), 10)

write.table(network, file = "GNET2_GRN.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(modules, file = "GNET2_modules.tsv", sep = "\t", quote = FALSE, row.names = FALSE)