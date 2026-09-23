library(tidyverse)
library(tidygraph)
library(ggraph)
library(igraph)
library(biomaRt)
library(Matrix)
#alternative file name 
source('dev/utilities/dataprocessingHelperFunctions.R')
#If you have a set of complex data-cleaning steps that you use across ten different projects,
#use this source function 

#load single cell data of seedling_3d
a=load('data/GSE226097_seedling_3d_230221.RData')

#load the original AraBoXcis network that was trained on bulk RNA-seq in seedlings
araboxcis = read.csv("data/gboxNetwork22C.csv", header = T)


print(a)

dim(gbox)
# gives number of rows and columns 
length(clust)

rownames(gbox)[1:10]
#This is indexing. It tells R, "Don't show me all 10,000 names; just give me the ones from position 1 to position 10."

colnames(gbox)[1:10]

as.matrix(gbox[1:50, 1:3])

clust[1:8]

table(clust)

plot(table(sort(clust)), xlab = 'cluster name', ylab = 'number of cells', main = 'Seddling_3d')

dim(araboxcis)

#first 4 rows
araboxcis[1:4,]

hist(araboxcis[,3])

#transcription factor 
tfs= unique(araboxcis[,1])
#unique function is also useful it creates a new vector with all duplicates removed 

#lets filter this to only include trancription factors are also in the sinle cell RNA seq data set 
tfSubs=tfs[which(tfs %in% rownames(gbox))]

length(tfSubs)

#filter genes and cells with very low values 
dim(gbox)

#get rid of cells that have less than 1%of the genes expressed  
thresh = 0.01
numberGenesPerCell = apply(gbox, 2, function(i){length(which(i>0))})
includeCells = which(numberGenesPerCell>(thresh*dim(gbox)[1])) 
gbox_filtered = gbox[,includeCells]
dim(gbox_filtered)

#now lets get rid of genes expressed in less than 1% of the cells 
numberCellsPerGene=apply(gbox_filtered, 1, function(i){length(which(i>0))})
includeGenes=which(numberCellsPerGene>(thresh*dim(gbox_filtered)[2]))
gbox_filtered=gbox_filtered[includeGenes,]
dim(gbox_filtered)

#--------------------------------------UMAP Analysis-----------------------------------------------------
#lets visualise it using umap
library(umap)
gbox.umap <- umap(gbox_filtered)
#do the cell type clusters group together if we only look at G-box related genes 
colours=rainbow(length(unique(clust)))
plot(gbox.umap$layout[,1], gbox.umap$layout[,2], col=colours[clust[includeCells]],
     pch=20, main='UMAP SEEDLING_3D', xlab = 'UMAP Component 1', ylab='UMAP Component 2')     
     
#pca
pca=prcomp(gbox_filtered, scale. = T, rank. = 5)
gbox.pca.umap <- umap(pca$x)
colours=rainbow(length(unique(clust)))
plot(gbox.pca.umap$layout[,1], gbox.pca.umap$layout[,2], col=colours[clust[includeCells]],
     pch=20, main='PCA UMAP SEEDLING_3D', xlab='UMAP Component 1', ylab='UMAP Component 2')

#------------------------------------------------Network extraction -----------------------------------------------
library(GENIE3)

net2=GENIE3(as.matrix(gbox), regulators = tfSubs, nTrees = 20)
 save(net, file = "seedling3d_network_nTress_20.RData")

 # Convert the GENIE3 output into a list of edges
links <- getLinkList(net2)

# Look at the top 10 strongest connections
head(links, 1000)

# Save the first 1000 rows
write.csv(head(links, 1000), "links_subset.csv", row.names = F)

#------------------------------------------- coverts the format to matrix-------------------------------
ginieOutput= convertToAdjacency(net, 0.05)
dim(ginieOutput)
ginieOutput[1:10,]

#---------------------------------------------Making the data frame for the new genie list --------------------
#load the network 
a= load('data/seedling3d_network_nTress_20.RData')
newNet = GENIE3::getLinkList(net2) 

#-----------------------------
#loading this file to compare with the existing data about Arabidopsis
araboxcis = read.csv('data/gboxNetwork22C.csv', header = T)

#-----------------------making data frame to compare the old and new data set in a more organised way--------------------

#gets the set of unique genes in your new network 
genesInNet = unique(c(newNet[,1], newNet[,2]))

#filter the AraBOXcis network to only contain genes that are in your new network 
araboxcisFiltered = araboxcis[which(araboxcis[,1] %in% genesInNet & araboxcis[,2] %in% genesInNet),]

#extract the top edges in your new network, to make your network the same size as the araboxcis Filtered network
newNetTopEdges = newNet[1:length(araboxcisFiltered[,1]),]

#reformat edges so is it more straightforward to comapre them
edgesNew = paste(newNetTopEdges[,1], newNetTopEdges[,2], sep = '_')
edgesold = paste(araboxcisFiltered[,1], araboxcisFiltered[,2], sep = '_')

#----------------------------now, you can come up with all the different parts of venn diagram-------------

library(VennDiagram)

grid.newpage() # Clears the plotting area
draw.pairwise.venn(
  area1 = length(edgesold),
  area2 = length(edgesNew),
  cross.area = length(which(edgesNew %in% edgesold)),
  category = c("AraBOXcis", "Since AraBOXcis"),
  fill = c("skyblue", "pink"),
  alpha = c(0.5, 0.5),
  lty = "blank",
  cex = 1,
  cat.cex = 0.5
)

#the overlap
length(which(edgesNew %in% edgesold))

#in new network only 
length(which(! (edgesNew %in% edgesold)))

#in old network only 
length(which(! (edgesold %in% edgesNew)))

#find important genes in the network 
tfsNew = table(newNetTopEdges[,1])
tfsold = table(araboxcisFiltered[,1])[names(tfsNew)]

#histogram of degree should look like a exponential distribution because biological networks ...
hist(as.numeric(tfsNew), main = 'SinceAraBOXcis', xlab = 'degree of TFs')

#lets see if the same TFs have high degrees in AraBoxcis and our new network
plot(as.numeric(tfsNew), as.numeric(tfsold), xlab = 'degree in SinceAraBOXcis', ylab = 'degree in AraBOXcis')

#lets print out the 20 TFs with highest degrees 
sort(tfsNew, decreasing = TRUE)[1:20]

#installed igraph and network package
library(igraph)
library(network)
#installed packages"pheatmap"
library(pheatmap)
simple_network <- graph_from_edgelist(as.matrix(newNetTopEdges[,c(1,2)]))

#node betweeness
node_betweenness_all <- betweenness(simple_network)
node_betweenness = node_betweenness_all[which(node_betweenness_all>0)]
sort(node_betweenness, decreasing = TRUE)[1:10]
write.csv(sort(node_betweenness, decreasing = TRUE)[1:10],file = "betwenness.csv")
plot(sort(node_betweenness))

#node centrality 
node_centrality_all <- alpha_centrality(simple_network, alpha = 0.5)
node_centrality = node_centrality_all[which(node_centrality_all>0)]
sort(node_centrality, decreasing = TRUE)[1:20]

plot(sort(node_centrality))

#node hub
node_hub_all <- hub_score(simple_network)$vector
node_hub = node_hub_all[which(node_hub_all>0)]
sort(node_hub, decreasing = TRUE)[1:10]
write.csv(sort(node_hub, decreasing = TRUE)[1:10], file = "hub.csv")

plot(sort(node_hub))

#btw vs central
plot(node_betweenness_all, node_centrality_all)

#hub vs central
plot(node_hub_all, node_centrality_all)

#hub vs between 
plot(node_hub_all, node_betweenness_all)

#------------------------------PAFway analysis--------------------------------------------------

#find association between GO terms in the network 
a=load('data/functionalData.RData')
source('dev/utilities/dataprocessingHelperFunctions.R')

pafwayOut=pafway(GOconcat, newNetTopEdges, unique(goOfInterest))
rownames(pafwayOut)=colnames(pafwayOut)

# Filter to only include rows and columns with at least one significant factor
atLeastOneSigRow=which(apply(pafwayOut, 1, function(i){length(which(i<0.05))})>0)
atLeastOneSigCol=which(apply(pafwayOut, 2, function(i){length(which(i<0.05))})>0)

pafwayInterestingOnly=pafwayOut[atLeastOneSigRow, atLeastOneSigCol]


# ── Log-transform the p-value matrix ─────────────────────────────────────────
# Add a tiny offset to avoid log(0) = -Inf for any exact zero p-values
pafwayLog = log(pafwayInterestingOnly + 1e-300, 10)

# ── Plot: upstream on Y axis, downstream on X axis, log10 p-value colour bar ─
pheatmap(
  pafwayLog,
  
  # Axis labels
  # Columns = source (upstream regulators) → X axis
  # Rows    = target (downstream genes)    → Y axis
  # pheatmap: rows go on Y, columns go on X — this matches directly
  # so no transposition needed; just label them clearly
  
  # Axis title annotation via a custom draw — pheatmap does not support
  # axis titles natively, so we use annotation_col / annotation_row as proxies
  # and add titles via grid after plotting
  
  fontsize_row = 6,
  fontsize_col = 6,
  
  # Colour scale
  color = colorRampPalette(c("steelblue", "white", "firebrick"))(100),
  
  # Flip so most negative log10(p) (most significant) is darkest red
  # log10(p-value) is always negative (p < 1), so the range is negative
  # More negative = more significant = stronger regulatory coupling
  breaks = seq(
    min(pafwayLog[is.finite(pafwayLog)]),
    max(pafwayLog[is.finite(pafwayLog)]),
    length.out = 101
  ),
  
  # Legend label
  legend_breaks = pretty(pafwayLog[is.finite(pafwayLog)], n = 5),
  legend_labels = paste0(
    pretty(pafwayLog[is.finite(pafwayLog)], n = 5)
  ),
  
  # Clustering
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  
  # Explicit axis labels
  angle_col    = 45,
  
  # Main title
  main = "PAFway pathway co-regulation\n(log\u2081\u2080 p-value)"
)

# ── Add axis titles using grid (runs immediately after pheatmap) ──────────────
library(grid)

grid.text(
  "Downstream pathway (target)",
  x    = 0.55,         # horizontal centre of the heatmap body
  y    = 0.01,         # just below the bottom edge
  just = "centre",
  gp   = gpar(fontsize = 10, fontface = "bold")
)

grid.text(
  "Upstream pathway (regulator)",
  x    = 0.01,         # just left of the left edge
  y    = 0.45,         # vertical centre of the heatmap body
  just = "centre",
  rot  = 90,           # rotate 90° for Y axis
  gp   = gpar(fontsize = 10, fontface = "bold")
)
#----------------------------------------------------------------------------
#trying to create  gene network to look at the genes related closely 
# Load necessary libraries
install.packages("pacman")
pacman::p_load(igraph, ggraph, tidygraph, dplyr, tidyverse)

#-----------------------------General Network analysis -----------------------
library(tidyverse)
library(tidygraph)
library(ggraph)
library(igraph)
library(biomaRt)
select <- dplyr::select
filter <- dplyr::filter

# ── Step 1: Load data (your original code) ────────────────────────────────
edges <- read.csv("data/links_subset.csv")
btwn  <- read.csv("data/betwenness.csv")
hubs  <- read.csv("data/hub.csv")

colnames(btwn) <- c("gene_id", "betweenness")
colnames(hubs) <- c("gene_id", "hub_score")

# ── Step 2: Convert gene IDs → gene names (NEW) ───────────────────────────
all_gene_ids <- unique(c(edges$regulatoryGene, edges$targetGene))

mart <- useMart("plants_mart",
                dataset = "athaliana_eg_gene",
                host    = "https://plants.ensembl.org")

id_to_name <- getBM(
  attributes = c("tair_locus", "external_gene_name"),
  filters    = "tair_locus",
  values     = all_gene_ids,
  mart       = mart
) %>%
  distinct(tair_locus, .keep_all = TRUE) %>%
  mutate(display_name = ifelse(
    is.na(external_gene_name) | external_gene_name == "",
    tair_locus,        # fall back to locus ID if no name exists
    external_gene_name
  ))

# Replace IDs in edges with gene names
edges <- edges %>%
  left_join(id_to_name %>% select(tair_locus, display_name),
            by = c("regulatoryGene" = "tair_locus")) %>%
  rename(regulator_name = display_name) %>%
  left_join(id_to_name %>% select(tair_locus, display_name),
            by = c("targetGene" = "tair_locus")) %>%
  rename(target_name = display_name) %>%
  mutate(
    regulatoryGene = coalesce(regulator_name, regulatoryGene),
    targetGene     = coalesce(target_name,    targetGene)
  ) %>%
  select(-regulator_name, -target_name)

# Manual overrides for genes with no Ensembl common name
manual_names <- tibble::tibble(
  tair_locus   = c("AT5G48560", "AT1G10120"),
  display_name = c("CIB2", "CIB4")  # replace right side when you find names
)

# ── Step 3: Build node list (your original code + gene names) ─────────────
all_genes <- unique(c(edges$regulatoryGene, edges$targetGene))
nodes <- data.frame(name = all_genes)

nodes <- nodes %>%
  # Join centrality scores — match on display_name now
  left_join(
    btwn %>%
      left_join(id_to_name, by = c("gene_id" = "tair_locus")) %>%
      mutate(display_name = coalesce(display_name, gene_id)),
    by = c("name" = "display_name")
  ) %>%
  left_join(
    hubs %>%
      left_join(id_to_name, by = c("gene_id" = "tair_locus")) %>%
      mutate(display_name = coalesce(display_name, gene_id)),
    by = c("name" = "display_name")
  ) %>%
  select(name, betweenness, hub_score) %>%
  replace_na(list(betweenness = 0, hub_score = 0))

# ── Step 4: Build graph (your original code) ──────────────────────────────
graph <- tbl_graph(
  nodes    = nodes,
  edges    = edges,
  directed = TRUE,
  node_key = "name"
)

# ── Step 5: Plot — force-directed + gene names (your original structure) ──
network_plot <- ggraph(graph, layout = "fr") +  # CHANGED: stress → fr (force-directed)
  
  geom_edge_link(
    aes(alpha = weight),
    arrow       = arrow(length = unit(2, "mm")),
    end_cap     = circle(3, "mm"),
    color       = "grey70"
  ) +
  
  geom_node_point(
    aes(size = betweenness, color = hub_score)
  ) +
  
  scale_color_gradientn(
    colors = c("skyblue", "gold", "red")
  ) +
  
  # CHANGED: name → display_name (gene names not IDs)
  geom_node_text(
    aes(label = ifelse(betweenness > 4000 | hub_score > 0.9, name, "")),
    repel        = TRUE,
    size         = 4,
    fontface     = "bold",
    max.overlaps = Inf
  ) +
  
  labs(
    title    = "Gene Regulatory Network Analysis",
    subtitle = "Node size: Betweenness | Color: Hub Score | Labels: gene names",
    color    = "Hub Score",
    size     = "Betweenness"
  ) +
  theme_graph(base_family = "Arial")

print(network_plot)

ggsave("network_gene_names.png",
       plot = network_plot, width = 14, height = 10, dpi = 300, bg = "white")


#-----------------------TFs for photomorphonesis and visualise in Network--------------------

# ══════════════════════════════════════════════════════════════════════════
# BLOCK 1 — Libraries
# ══════════════════════════════════════════════════════════════════════════
library(Matrix)
library(tidyverse)
library(tidygraph)
library(ggraph)
library(igraph)
library(biomaRt)

# ══════════════════════════════════════════════════════════════════════════
# BLOCK 2 — Load all data
# ══════════════════════════════════════════════════════════════════════════
a <- load('data/functionalData.RData')
source('dev/utilities/dataprocessingHelperFunctions.R')
b <- load('data/seedling3d_network_nTress_20.RData')

edges <- read.csv("data/links_subset.csv")

# ══════════════════════════════════════════════════════════════════════════
# BLOCK 3 — Define GO lookup function
# ══════════════════════════════════════════════════════════════════════════
get_genes_by_GO <- function(term){
  names(GOconcat)[grep(term, GOconcat, ignore.case = TRUE)]
}

# ══════════════════════════════════════════════════════════════════════════
# BLOCK 4 — Recalculate betweenness for ALL genes
# ══════════════════════════════════════════════════════════════════════════
graph_full       <- graph_from_data_frame(
  edges[, c("regulatoryGene", "targetGene")],
  directed = TRUE)

full_betweenness <- betweenness(graph_full, directed = TRUE, normalized = FALSE)
full_hub         <- hub_score(graph_full)$vector

nodes_full <- data.frame(
  name        = names(full_betweenness),
  betweenness = as.numeric(full_betweenness),
  hub_score   = as.numeric(full_hub[names(full_betweenness)])
) %>% replace_na(list(betweenness = 0, hub_score = 0))

# ══════════════════════════════════════════════════════════════════════════
# BLOCK 4b — Fetch gene names automatically from biomaRt
# ══════════════════════════════════════════════════════════════════════════
all_gene_ids <- unique(c(edges$regulatoryGene, edges$targetGene))

mart <- tryCatch({
  useEnsemblGenomes(biomart = "plants_mart",
                    dataset = "athaliana_eg_gene")
}, error = function(e){
  useMart("plants_mart",
          dataset = "athaliana_eg_gene",
          host    = "https://plants.ensembl.org")
})

id_to_name <- getBM(
  attributes = c("tair_locus", "external_gene_name", "description"),
  filters    = "tair_locus",
  values     = all_gene_ids,
  mart       = mart
) %>%
  mutate(display_name = ifelse(
    is.na(external_gene_name) | external_gene_name == "",
    tair_locus,
    external_gene_name
  )) %>%
  distinct(tair_locus, .keep_all = TRUE)

cat("Gene names resolved  :", sum(id_to_name$display_name != id_to_name$tair_locus), "\n")
cat("Still using locus ID :", sum(id_to_name$display_name == id_to_name$tair_locus), "\n")

# Join names onto nodes
nodes_full <- nodes_full %>%
  left_join(id_to_name %>% select(tair_locus, display_name),
            by = c("name" = "tair_locus")) %>%
  mutate(display_name = ifelse(is.na(display_name), name, display_name))

# ══════════════════════════════════════════════════════════════════════════
# BLOCK 5 — GO lookups and filter
# ══════════════════════════════════════════════════════════════════════════
photo_genes <- get_genes_by_GO("photomorphogenesis")
tf_genes    <- get_genes_by_GO("transcription factor")
light_genes <- get_genes_by_GO("light")
photo_tfs   <- intersect(tf_genes, c(photo_genes, light_genes))

edges_photo         <- edges %>% filter(regulatoryGene %in% photo_tfs)
genes_in_subnetwork <- unique(c(edges_photo$regulatoryGene,
                                edges_photo$targetGene))
nodes_photo         <- nodes_full %>% filter(name %in% genes_in_subnetwork)

# ══════════════════════════════════════════════════════════════════════════
# BLOCK 6 — Plot with gene names
# ══════════════════════════════════════════════════════════════════════════
graph_photo <- tbl_graph(nodes    = nodes_photo,
                         edges    = edges_photo,
                         directed = TRUE,
                         node_key = "name")

network_plot <- ggraph(graph_photo, layout = "fr") +
  geom_edge_link(
    aes(alpha = weight),
    arrow   = arrow(length = unit(2, "mm")),
    end_cap = circle(3, "mm"),
    color   = "grey70"
  ) +
  geom_node_point(aes(size = betweenness, color = hub_score)) +
  scale_color_gradientn(colors = c("skyblue", "gold", "red"),
                        name   = "Hub Score") +
  scale_size_continuous(range = c(2, 12), name = "Betweenness") +
  geom_node_text(
    aes(label = ifelse(
      betweenness > quantile(betweenness, 0.75) |
        hub_score   > quantile(hub_score,   0.75),
      display_name, ""             # gene name shown here
    )),
    repel        = TRUE,
    size         = 4,
    fontface     = "bold",
    max.overlaps = Inf,
    family       = "sans"
  ) +
  
  labs(title    = "Photomorphogenesis TF regulatory network",
       subtitle = "Node size: Betweenness | Color: Hub Score | Labels: gene names") +
  theme_graph(base_family = "Arial")

print(network_plot)




