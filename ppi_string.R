# Load necessary libraries
library(STRINGdb)
library(igraph)
library(dplyr)

# --- 1. Define Your List of Proteins/Genes of Interest ---
# This would typically be a list of differentially expressed genes,
# genes from a specific pathway, or any proteins central to your research.
# Use official gene symbols (e.g., HGNC symbols for human genes).
# For this example, we'll use a small set of well-known human genes.

my_proteins <- c(
  "TP53", "MYC", "EGFR", "AKT1", "STAT3", "JUN", "FOS", "RELA", "NFKB1",
  "TNF", "IL6", "VEGFA", "APP", "PSEN1", "MAPT"
)

cat("My list of proteins for PPI analysis:\n")
print(my_proteins)
cat("\n")

# --- 2. Initialize STRINGdb Object ---
# This step sets up the connection to the STRING database.
# 'version' refers to the STRING database version (e.g., "11.5", "12.0").
# 'species' refers to the NCBI Taxonomy ID (e.g., 9606 for Homo sapiens, 10090 for Mus musculus).
# 'score_threshold' is the minimum STRING combined score to consider an interaction.
# Higher scores indicate higher confidence interactions.

# Find species NCBI Taxonomy ID:
# You can search here: https://www.ncbi.nlm.nih.gov/taxonomy
# For Homo sapiens (human): 9606
# For Mus musculus (mouse): 10090

string_db <- STRINGdb$new(
  version = "12.0",            # Using STRINGdb version 12.0 (latest common release)
  species = 9606,              # Homo sapiens
  score_threshold = 400,       # Minimum combined score (0-999). 400 is medium confidence.
  input_directory = ""         # Optional: directory to save STRINGdb data
)

cat("STRINGdb object initialized for Homo sapiens (v12.0) with score threshold 400.\n\n")

# --- 3. Map Your Proteins to STRING Identifiers ---
# STRINGdb uses its own internal identifiers. This step maps your input gene
# symbols to STRING protein IDs. Genes that cannot be mapped will be excluded.

# The 'get_STRING_ids' function helps map gene symbols to STRING IDs.
# 'entrezGene' column is the default for input, but 'alias' also works well for gene symbols.
# Ensure your input column matches the 'via' argument.

mapped_proteins_df <- string_db$get_STRING_ids(
  my_proteins,
  data_type_id = "genes",  # Specify data type as genes
  alias_external_id = "my_genes", # Name for your input column in the resulting dataframe
  via = "alias"          # Use the 'alias' column in STRINGdb for mapping gene symbols
)

cat("Mapping status:\n")
print(head(mapped_proteins_df))
cat(paste("\nSuccessfully mapped", nrow(mapped_proteins_df), "out of", length(my_proteins), "proteins.\n\n"))

# Filter out unmapped proteins
mapped_proteins_only <- mapped_proteins_df %>%
  filter(!is.na(STRING_id)) %>% # Keep only rows where STRING_id is not NA
  select(STRING_id, my_genes)   # Select the STRING_id and your original gene symbol

# --- 4. Get PPI Interactions from STRINGdb ---
# This function retrieves all interactions for your mapped proteins above the
# specified score threshold.

ppi_interactions <- string_db$get_interactions(mapped_proteins_only$STRING_id)

cat("Retrieved", nrow(ppi_interactions), "PPI interactions above the threshold.\n\n")

# --- 5. Convert to an igraph Object for Network Analysis ---
# An igraph object is ideal for network visualization and further analysis.

# First, create an edge list (pairs of interacting proteins)
# We need to map STRING_id back to original gene symbols for better readability
edge_list_df <- ppi_interactions %>%
  left_join(mapped_proteins_only, by = c("node1" = "STRING_id")) %>%
  rename(protein1_symbol = my_genes) %>%
  left_join(mapped_proteins_only, by = c("node2" = "STRING_id")) %>%
  rename(protein2_symbol = my_genes) %>%
  select(protein1_symbol, protein2_symbol, combined_score)

# Create the igraph object
# The graph_from_data_frame function takes a dataframe with columns for edges
# and optionally, vertex (node) attributes.
# The 'directed=FALSE' means interactions are undirected (Protein A interacts with B, so B interacts with A).

ppi_graph <- graph_from_data_frame(
  d = edge_list_df,      # Edge list dataframe
  directed = FALSE,      # Undirected graph for typical PPIs
  vertices = data.frame(name = unique(c(edge_list_df$protein1_symbol, edge_list_df$protein2_symbol)))
)

cat("Network graph created with", vcount(ppi_graph), "nodes and", ecount(ppi_graph), "edges.\n\n")

# --- 6. Visualize the PPI Network ---
# Basic plot using igraph. For more advanced visualization, consider Cytoscape
# or other dedicated network visualization tools.

# Set plot parameters for better aesthetics (optional)
V(ppi_graph)$size <- 10 # Node size
V(ppi_graph)$label.cex <- 0.8 # Label font size
V(ppi_graph)$label.color <- "black" # Label color
E(ppi_graph)$width <- E(ppi_graph)$combined_score / 200 # Edge thickness based on score
E(ppi_graph)$color <- "gray" # Edge color

# Use a specific layout for better readability
# For example, 'layout_with_fr' for Fruchterman-Reingold layout
# or 'layout_nicely' for a good default layout.
set.seed(123) # for reproducibility of layout

# Plotting the network
plot(
  ppi_graph,
  layout = layout_with_fr(ppi_graph), # Force-directed layout
  vertex.label = V(ppi_graph)$name,   # Show gene symbols as labels
  vertex.frame.color = "darkgray",    # Node border color
  vertex.color = "lightblue",         # Node fill color
  main = "Protein-Protein Interaction Network" # Plot title
)

# --- 7. Basic Network Analysis (Optional) ---
# Calculate common network properties

cat("Basic Network Properties:\n")
cat(paste("Number of Nodes (Proteins):", vcount(ppi_graph), "\n"))
cat(paste("Number of Edges (Interactions):", ecount(ppi_graph), "\n"))
cat(paste("Network Density:", graph.density(ppi_graph), "\n")) # Proportion of actual edges to all possible edges

# Identify highly connected proteins (hubs)
node_degrees <- degree(ppi_graph)
top_hub_proteins <- sort(node_degrees, decreasing = TRUE)[1:min(5, length(node_degrees))]
cat("\nTop 5 Hub Proteins (by degree):\n")
print(top_hub_proteins)

# Find connected components (clusters)
components <- components(ppi_graph)
cat(paste("\nNumber of connected components (clusters):", components$no, "\n"))
if (components$no > 1) {
  cat("The network is disconnected, suggesting distinct functional modules.\n")
} else {
  cat("The network is fully connected.\n")
}