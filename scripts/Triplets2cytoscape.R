# ============================================================
# Convert triplet stability table to gene network
# Input  : triplet_stability_ranking.csv
# Output : triplet_network_edges.csv
#          triplet_network_nodes.csv
# ============================================================

cat("\nSelect triplet_stability_ranking.csv\n")
f <- file.choose()

df <- read.csv(f, stringsAsFactors = FALSE)

required <- c("A","B","C","stability_score")
if (!all(required %in% names(df))) {
  stop("CSV must contain columns: A, B, C, stability_score")
}

# ------------------------------------------------------------
# Convert each triplet into three edges
# ------------------------------------------------------------
edges <- rbind(
  data.frame(source=df$A, target=df$B, weight=df$stability_score),
  data.frame(source=df$A, target=df$C, weight=df$stability_score),
  data.frame(source=df$B, target=df$C, weight=df$stability_score)
)

# collapse duplicate edges
edge_key <- paste(pmin(edges$source,edges$target),
                  pmax(edges$source,edges$target),
                  sep="_")

edges$key <- edge_key

agg <- aggregate(weight ~ key, data=edges, mean)

split <- strsplit(agg$key,"_")

agg$source <- sapply(split,"[",1)
agg$target <- sapply(split,"[",2)

edges_final <- agg[,c("source","target","weight")]

# ------------------------------------------------------------
# Build node table (gene importance)
# ------------------------------------------------------------
genes <- unique(c(df$A,df$B,df$C))

node_degree <- table(c(edges_final$source,edges_final$target))

nodes <- data.frame(
  gene = genes,
  degree = as.numeric(node_degree[genes])
)

nodes$degree[is.na(nodes$degree)] <- 0

# ------------------------------------------------------------
# Write files
# ------------------------------------------------------------
outdir <- dirname(f)

write.csv(edges_final,
          file.path(outdir,"triplet_network_edges.csv"),
          row.names=FALSE)

write.csv(nodes,
          file.path(outdir,"triplet_network_nodes.csv"),
          row.names=FALSE)

cat("\nNetwork files written to:\n",outdir,"\n")