#!/usr/bin/env Rscript
# ======================================================================
# NLP-assisted triage for systematic review screening (PRISMA-safe)
# - Reads local references (CSV) exported from Zotero / PubMed tools
# - Uses OpenAI embeddings to RANK (not exclude) references for screening
# - Outputs ranked list + screening bins + optional clustering
#
# Requirements:
#   - Set OPENAI_API_KEY in your environment (recommended)
#       Sys.setenv(OPENAI_API_KEY="...")  # in .Renviron or before running
#
# Notes:
#   - This script does NOT automatically exclude papers.
#   - It produces a priority order and bins to speed human screening.
# ======================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("httr2", quietly = TRUE)) install.packages("httr2")
  if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")
  if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr")
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
  if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
  if (!requireNamespace("purrr", quietly = TRUE)) install.packages("purrr")
  if (!requireNamespace("tibble", quietly = TRUE)) install.packages("tibble")
  if (!requireNamespace("digest", quietly = TRUE)) install.packages("digest")
})

library(httr2)
library(jsonlite)
library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tibble)
library(digest)

# -----------------------------
# User parameters
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

# Usage:
#   Rscript nlp_triage.R /Users/schmidtc/Library/CloudStorage/Dropbox/RSTUDIO/2025_LiverMan_1/References/Systemic_Canalization.csv
in_path  <- ifelse(length(args) >= 1, args[[1]], "")
out_dir  <- ifelse(length(args) >= 2, args[[2]], file.path(getwd(), "outputs_nlp_triage"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!nzchar(in_path) || !file.exists(in_path)) {
  stop("Input file not found. Provide a local CSV path as the first argument.")
}

api_key <- Sys.getenv("OPENAI_API_KEY")
if (!nzchar(api_key)) {
  stop("OPENAI_API_KEY not found in environment. Set it via Sys.setenv() or .Renviron.")
}

# OpenAI settings (embeddings endpoint)
openai_base <- "https://api.openai.com/v1"
embedding_model <- "text-embedding-3-large"  # per OpenAI docs; high quality embeddings

# Cost/throughput controls
max_chars_per_doc <- 12000   # truncate long abstracts to control token usage
batch_size <- 64             # tune for your machine/network
sleep_between_batches_sec <- 0.2

# Clustering controls
do_clustering <- TRUE
k_clusters <- 12

# -----------------------------
# Broad-net query set (edit as desired)
# -----------------------------
queries <- c(
  # Core concept
  "developmental canalization definition operationalization quantitative metrics",
  "Waddington landscape canalization robustness developmental constraint",
  "developmental stability phenotypic robustness molecular data",
  
  # Variance / constraint
  "variance minimization developmental stability coefficient of variation CV IQR MAD",
  "expression variability across development variance constraints",
  
  # Dimensionality / manifolds
  "effective dimensionality intrinsic dimension developmental trajectory gene expression",
  "PCA eigenvalue spectrum low-dimensional structure developmental",
  "manifold learning developmental potency intrinsic dimension transcriptome",
  
  # Networks / structure
  "correlation structure stability network rewiring development module stability",
  "edge turnover pathway network development transcriptome metabolome",
  
  # Omics layers (broad)
  "transcriptomics RNA-seq developmental time course gene expression programs",
  "single-cell RNA-seq developmental trajectory differentiation",
  "metabolomics metabolome developmental time course LC-MS GC-MS",
  "proteomics proteome developmental time course LC-MS/MS",
  "epigenetics chromatin accessibility ATAC-seq ChIP-seq developmental",
  "DNA methylation developmental dynamics epigenome",
  "multi-omics integration developmental time course cross-layer alignment"
)

# -----------------------------
# Helpers: robust column picking for Zotero-ish exports
# -----------------------------
pick_first_existing <- function(df, candidates) {
  cand <- candidates[candidates %in% names(df)]
  if (length(cand) == 0) return(NA_character_)
  cand[[1]]
}

# Attempt to support common Zotero CSV column names
# Typical Zotero CSV columns include: "Title", "Abstract Note", "Author", "Date", "Publication Title", "DOI", "Url"
col_title   <- pick_first_existing(read_csv(in_path, show_col_types = FALSE, n_max = 1),
                                   c("Title", "title", "Article Title"))
col_abs     <- pick_first_existing(read_csv(in_path, show_col_types = FALSE, n_max = 1),
                                   c("Abstract Note", "Abstract", "abstract", "AbstractText"))
col_doi     <- pick_first_existing(read_csv(in_path, show_col_types = FALSE, n_max = 1),
                                   c("DOI", "doi"))
col_year    <- pick_first_existing(read_csv(in_path, show_col_types = FALSE, n_max = 1),
                                   c("Date", "Year", "year", "Publication Year"))
col_journal <- pick_first_existing(read_csv(in_path, show_col_types = FALSE, n_max = 1),
                                   c("Publication Title", "Journal", "journal"))
col_auth    <- pick_first_existing(read_csv(in_path, show_col_types = FALSE, n_max = 1),
                                   c("Author", "Authors", "author"))

if (is.na(col_title)) stop("Could not find a Title column. Please export a Zotero-style CSV with a Title field.")
if (is.na(col_abs)) message("Warning: No abstract column found. Ranking will rely on titles only (weaker).")

# -----------------------------
# Load and normalize input
# -----------------------------
raw <- read_csv(in_path, show_col_types = FALSE) %>%
  mutate(row_id = row_number()) %>%
  mutate(
    title = as.character(.data[[col_title]]),
    abstract = if (!is.na(col_abs)) as.character(.data[[col_abs]]) else "",
    doi = if (!is.na(col_doi)) as.character(.data[[col_doi]]) else NA_character_,
    year_raw = if (!is.na(col_year)) as.character(.data[[col_year]]) else NA_character_,
    journal = if (!is.na(col_journal)) as.character(.data[[col_journal]]) else NA_character_,
    authors = if (!is.na(col_auth)) as.character(.data[[col_auth]]) else NA_character_
  ) %>%
  mutate(
    year = suppressWarnings(as.integer(str_extract(year_raw, "\\b(19\\d{2}|20\\d{2})\\b"))),
    text_for_embed = str_squish(paste0("TITLE: ", title, "\nABSTRACT: ", abstract)),
    text_for_embed = str_sub(text_for_embed, 1, max_chars_per_doc),
    text_hash = vapply(text_for_embed, function(x) digest(x, algo = "sha1"), character(1))
  )

if (nrow(raw) == 0) stop("No rows found in input file.")
message(sprintf("Loaded %d references.", nrow(raw)))

# -----------------------------
# Embeddings API call
# -----------------------------
embed_texts <- function(text_vec, api_key, model = embedding_model) {
  req <- request(paste0(openai_base, "/embeddings")) %>%
    req_headers(
      "Authorization" = paste("Bearer", api_key),
      "Content-Type" = "application/json"
    ) %>%
    req_body_json(list(
      model = model,
      input = text_vec
    ))
  
  resp <- req_perform(req)
  if (resp_status(resp) >= 300) {
    stop("Embeddings request failed: ", resp_status(resp), "\n", resp_body_string(resp))
  }
  
  body <- resp_body_json(resp)
  # body$data is list of {embedding=...}
  embeddings <- lapply(body$data, function(d) d$embedding)
  embeddings
}

cosine_sim <- function(a, b) {
  # a,b numeric vectors
  denom <- sqrt(sum(a*a)) * sqrt(sum(b*b))
  if (denom == 0) return(0)
  sum(a*b) / denom
}

# -----------------------------
# Cache embeddings to disk (so you can re-run cheaply)
# -----------------------------
cache_path <- file.path(out_dir, "embeddings_cache.rds")
cache <- list()
if (file.exists(cache_path)) {
  cache <- readRDS(cache_path)
  if (!is.list(cache)) cache <- list()
}

get_cached <- function(hash) {
  if (!is.null(cache[[hash]])) return(cache[[hash]])
  NULL
}
set_cached <- function(hash, emb) {
  cache[[hash]] <<- emb
}

# Build doc embeddings, using cache
doc_embeddings <- vector("list", nrow(raw))
need_idx <- which(vapply(raw$text_hash, function(h) is.null(get_cached(h)), logical(1)))

if (length(need_idx) > 0) {
  message(sprintf("Embedding %d new documents (cached: %d).", length(need_idx), nrow(raw) - length(need_idx)))
  
  batches <- split(need_idx, ceiling(seq_along(need_idx) / batch_size))
  for (bi in seq_along(batches)) {
    idx <- batches[[bi]]
    texts <- raw$text_for_embed[idx]
    
    embs <- embed_texts(texts, api_key = api_key, model = embedding_model)
    
    for (j in seq_along(idx)) {
      h <- raw$text_hash[idx[[j]]]
      set_cached(h, embs[[j]])
    }
    
    saveRDS(cache, cache_path)
    Sys.sleep(sleep_between_batches_sec)
    message(sprintf("  Completed batch %d / %d", bi, length(batches)))
  }
} else {
  message("All document embeddings found in cache.")
}

for (i in seq_len(nrow(raw))) {
  doc_embeddings[[i]] <- get_cached(raw$text_hash[[i]])
}

# Embed query strings (small; no need to cache)
message("Embedding query set...")
query_embeddings <- embed_texts(queries, api_key = api_key, model = embedding_model)

# -----------------------------
# Score each document by best query similarity + mean similarity
# -----------------------------
message("Scoring documents...")
score_df <- tibble(row_id = raw$row_id) %>%
  mutate(
    best_sim = NA_real_,
    mean_sim = NA_real_,
    best_query = NA_character_
  )

for (i in seq_len(nrow(raw))) {
  d <- doc_embeddings[[i]]
  sims <- vapply(query_embeddings, function(q) cosine_sim(d, q), numeric(1))
  best_i <- which.max(sims)
  score_df$best_sim[i] <- sims[[best_i]]
  score_df$mean_sim[i] <- mean(sims)
  score_df$best_query[i] <- queries[[best_i]]
}

# -----------------------------
# Create screening bins (heuristic thresholds; adjust as you like)
# These bins are NOT exclusions; they are triage priority bands.
# -----------------------------
# Use quantiles so it adapts to corpus.
q95 <- quantile(score_df$best_sim, 0.95, na.rm = TRUE)
q80 <- quantile(score_df$best_sim, 0.80, na.rm = TRUE)
q50 <- quantile(score_df$best_sim, 0.50, na.rm = TRUE)

score_df <- score_df %>%
  mutate(
    triage_bin = case_when(
      best_sim >= q95 ~ "VeryLikely",
      best_sim >= q80 ~ "Likely",
      best_sim >= q50 ~ "Unclear",
      TRUE ~ "LowPriority"
    )
  )

# -----------------------------
# Optional clustering (coverage sanity-check)
# -----------------------------
cluster_df <- NULL
if (do_clustering) {
  message("Clustering embeddings (k-means)...")
  # Convert list of embeddings to matrix (may be large; workable for 2,500)
  emb_mat <- do.call(rbind, lapply(doc_embeddings, function(x) as.numeric(x)))
  
  set.seed(1)
  km <- kmeans(emb_mat, centers = k_clusters, iter.max = 50)
  cluster_df <- tibble(row_id = raw$row_id, cluster = as.integer(km$cluster))
}

# -----------------------------
# Assemble final ranked table
# -----------------------------
out <- raw %>%
  select(row_id, title, abstract, authors, year, journal, doi) %>%
  left_join(score_df, by = "row_id")

if (!is.null(cluster_df)) {
  out <- out %>% left_join(cluster_df, by = "row_id")
}

out_ranked <- out %>%
  arrange(desc(best_sim), desc(mean_sim)) %>%
  mutate(rank = row_number())

# -----------------------------
# Write outputs
# -----------------------------
ranked_path <- file.path(out_dir, "Ranked_References.csv")
bins_path   <- file.path(out_dir, "Screening_Bins.csv")

write_csv(out_ranked, ranked_path)
write_csv(out_ranked %>% select(rank, row_id, triage_bin, best_sim, best_query, title), bins_path)

message("Done.")
message("Ranked references: ", ranked_path)
message("Screening bins:    ", bins_path)
message("Embeddings cache:  ", cache_path)

# Suggested next step:
# - Start screening from triage_bin == "VeryLikely" then "Likely", etc.
# - Record counts per bin for transparency, but do NOT treat bins as exclusions.
