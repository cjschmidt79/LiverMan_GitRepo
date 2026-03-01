# --------------------------------------------
# 🔬 Refactored Triplet Analysis with Empirical FDR, Q-Value, and Bootstrapping
# --------------------------------------------

# Load Required Packages
required_packages <- c("ggplot2", "DT", "htmlwidgets", "tools", "parallel", "glue", "qvalue")
for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE)) install.packages(pkg)
    library(pkg, character.only = TRUE)
}
# Ensure qvalue is installed from Bioconductor if not found
if (!requireNamespace("qvalue", quietly = TRUE)) {
    message("📦 Installing qvalue from Bioconductor...")
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
    }
    BiocManager::install("qvalue")
}
library(qvalue)

# Input and Output Setup
file_path <- file.choose()
base_dir <- dirname(file.choose())
input_name <- tools::file_path_sans_ext(basename(file_path))
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(base_dir, paste0(input_name, "_", timestamp))
dir.create(output_dir, recursive = TRUE)
viz_dir <- file.path(output_dir, "visualizations")
dir.create(viz_dir, recursive = TRUE)

# Load Expression Data
cat("\n🧬 Loading data and assigning groups\n")
data_raw <- read.csv(file_path, row.names = 1, check.names = FALSE)
sample_names <- colnames(data_raw)
group_labels <- sub("^([A-Za-z0-9]+)[ _]?.*$", "\\1", sample_names)
unique_groups <- unique(group_labels)
print(unique_groups)

# Prompt for comma-separated early/late prefixes
early_input <- readline("Enter label prefix(es) for 'early' (e.g., 4,8): ")
late_input  <- readline("Enter label prefix(es) for 'late' (e.g., 20,24): ")

early_prefixes <- trimws(unlist(strsplit(early_input, ",")))
late_prefixes  <- trimws(unlist(strsplit(late_input, ",")))

# Helper to match any prefix
match_prefix <- function(label, prefixes) {
    any(sapply(prefixes, function(p) startsWith(label, p)))
}

# Assign groups using prefix matching
group_factor <- sapply(group_labels, function(lab) {
    if (match_prefix(lab, early_prefixes)) return("early")
    if (match_prefix(lab, late_prefixes)) return("late")
    return(NA)
})
group_factor <- factor(group_factor)

# Validate assignment
if (any(is.na(group_factor))) {
    cat("❌ The following samples could not be assigned to a group:\n")
    print(sample_names[is.na(group_factor)])
    stop("Fix group assignments and try again.")
}

cat("\n✅ Group assignment complete:\n")
print(table(group_factor))

# Set Parameters
min_expr <- as.numeric(readline("Minimum expression (default 10): "))
if (is.na(min_expr)) min_expr <- 10
min_var <- as.numeric(readline("Minimum variance (default 0.05): "))
if (is.na(min_var)) min_var <- 0.05
num_permutations <- as.numeric(readline("Permutations (default 1000): "))
if (is.na(num_permutations)) num_permutations <- 1000
num_bootstrap <- as.numeric(readline("Bootstraps (default 1000): "))
if (is.na(num_bootstrap)) num_bootstrap <- 1000
fdr_filter <- as.numeric(readline("Empirical FDR cutoff for further testing (default 0.2): "))
if (is.na(fdr_filter)) fdr_filter <- 0.2
fdr_qval <- as.numeric(readline("Q-value threshold (default 0.1): "))
if (is.na(fdr_qval)) fdr_qval <- 0.1

# Parallel
cores <- parallel::detectCores()
def_cores <- max(1, cores - 2)
user_cores <- as.numeric(readline(paste0("Cores to use (default ", def_cores, "): ")))
if (is.na(user_cores) || user_cores < 1 || user_cores > cores) user_cores <- def_cores

# Filter Expression
expr_filter <- rowSums(data_raw >= min_expr) >= 2
var_filter <- apply(data_raw, 1, var) >= min_var
data_filtered <- data_raw[expr_filter & var_filter, ]

group_factor <- droplevels(group_factor)
triplets <- combn(rownames(data_filtered), 3, simplify = FALSE)

# Main Triplet Function
analyze_triplet <- function(triplet) {
    A <- triplet[1]; B <- triplet[2]; C <- triplet[3]
    vals_A <- as.numeric(data_filtered[A, ])
    vals_B <- as.numeric(data_filtered[B, ])
    vals_C <- as.numeric(data_filtered[C, ])
    if (any(c(vals_A, vals_B, vals_C) <= 0)) return(NULL)
    
    df <- data.frame(A = log(vals_A), BDivC = log(vals_B / vals_C), group = group_factor)
    if (sd(df$BDivC) < 1e-4) return(NULL)
    
    model <- tryCatch(lm(A ~ BDivC * group, data = df), error = function(e) return(NULL))
    if (is.null(model)) return(NULL)
    
    s <- summary(model)
    p_int <- ifelse("BDivC:grouplate" %in% rownames(s$coefficients),
                    s$coefficients["BDivC:grouplate", 4], NA)
    main_slope <- coef(model)["BDivC"]
    int_slope <- ifelse(!is.na(p_int), coef(model)["BDivC:grouplate"], 0)
    
    perm_pvals_int <- replicate(num_permutations, {
        df$group <- sample(df$group)
        m_perm <- tryCatch(lm(A ~ BDivC * group, data = df), error = function(e) NULL)
        if (!is.null(m_perm) && "BDivC:grouplate" %in% rownames(coef(summary(m_perm)))) {
            coef(summary(m_perm))["BDivC:grouplate", 4]
        } else NA
    })
    
    emp_fdr_int <- mean(perm_pvals_int <= p_int, na.rm = TRUE)
    
    data.frame(A, B, C,
               slope_grp1 = round(main_slope, 3),
               slope_grp2 = round(main_slope + int_slope, 3),
               abs_slope = round(abs(main_slope), 3),
               r_squared = round(s$r.squared, 3),
               p_interaction = p_int,
               fdr_interaction = emp_fdr_int)
}

cat("\n⏳ Analyzing triplets...\n")
results <- parallel::mclapply(triplets, analyze_triplet, mc.cores = user_cores)
results_df <- do.call(rbind, Filter(Negate(is.null), results))

# Apply Empirical FDR Filter
filtered_results <- results_df[results_df$fdr_interaction < fdr_filter, ]

# ---- Validate and Clean P-values ----
valid_pvals <- !is.na(filtered_results$p_interaction) & 
    !is.nan(filtered_results$p_interaction) &
    is.finite(filtered_results$p_interaction) &
    filtered_results$p_interaction >= 0 & 
    filtered_results$p_interaction <= 1

num_invalid <- sum(!valid_pvals)
if (num_invalid > 0) {
    cat(paste0("⚠️ Removed ", num_invalid, " rows with invalid p-values before q-value adjustment.\n"))
    write.csv(filtered_results[!valid_pvals, ],
              file.path(output_dir, "dropped_invalid_pvals.csv"), row.names = FALSE)
}

# ---- Subset Valid P-values ----
valid_pvals_vec <- filtered_results$p_interaction[valid_pvals]

if (length(valid_pvals_vec) < 10) stop("❌ Too few valid p-values for q-value correction.")

# ---- Plot Distribution ----
png(file.path(output_dir, "pval_distribution.png"))
hist(valid_pvals_vec, breaks = 50, main = "Histogram of Valid p-values", xlab = "p-value")
dev.off()

#Apply qvalue() with full fault tolerance ----

# ---- Apply qvalue() with guaranteed fallback ----

filtered_results$qval_interaction <- NA

# Clean p-values
valid_idx <- !is.na(filtered_results$p_interaction) &
    !is.nan(filtered_results$p_interaction) &
    is.finite(filtered_results$p_interaction) &
    filtered_results$p_interaction >= 0 &
    filtered_results$p_interaction <= 1

valid_pvals_vec <- filtered_results$p_interaction[valid_idx]

# Safety check
if (length(valid_pvals_vec) < 10) stop("❌ Too few valid p-values to run qvalue().")

# Save for audit
write.csv(data.frame(p_interaction = valid_pvals_vec),
          file.path(output_dir, "valid_pvalues_used_for_qvalue.csv"),
          row.names = FALSE)

# Try qvalue with ultra-safe lambda range
qval_success <- FALSE
tryCatch({
    qobj <- qvalue::qvalue(valid_pvals_vec,
                           lambda = seq(0.1, 0.8, 0.1),  # safer lambda range
                           pi0.method = "smooth")
    if (any(is.na(qobj$qvalues))) stop("qvalue() output contains NA q-values.")
    filtered_results$qval_interaction[valid_idx] <- qobj$qvalues
    qval_success <- TRUE
}, error = function(e) {
    msg <- paste0("qvalue() failed: ", conditionMessage(e))
    cat("❌ ", msg, "\n")
    writeLines(msg, file.path(output_dir, "qvalue_error_log.txt"))
})

# ---- Optional Fallback: fdrtool ----
if (!qval_success) {
    if (!requireNamespace("fdrtool", quietly = TRUE)) {
        install.packages("fdrtool")
    }
    library(fdrtool)
    tryCatch({
        fdr <- fdrtool(valid_pvals_vec, statistic = "pvalue", plot = FALSE, verbose = FALSE)
        filtered_results$qval_interaction[valid_idx] <- fdr$qval
        qval_success <- TRUE
        cat("✅ fdrtool() used as fallback to compute q-values.\n")
    }, error = function(e2) {
        msg2 <- paste0("fdrtool() also failed: ", conditionMessage(e2))
        writeLines(msg2, file.path(output_dir, "qvalue_error_log.txt"))
        stop("🛑 Both qvalue() and fdrtool() failed. See qvalue_error_log.txt")
    })
}


# ---- Filter & Save ----
sig_results <- filtered_results[filtered_results$qval_interaction < fdr_qval, ]

write.csv(filtered_results, file.path(output_dir, "all_triplets.csv"), row.names = FALSE)
write.csv(sig_results, file.path(output_dir, "significant_triplets.csv"), row.names = FALSE)

# ---- Bootstrap & Save All to Single RDS ----
n_top <- as.integer(readline("How many top triplets to bootstrap and plot? "))
top_boot <- head(filtered_results[order(filtered_results$qval_interaction), ], n_top)

all_boots <- list()
bootstrap_rows <- list()  # summary rows

for (i in seq_len(nrow(top_boot))) {
    row <- top_boot[i, ]
    A <- row$A; B <- row$B; C <- row$C
    triplet_name <- paste(A, B, C, sep = "_")
    cat("🔁 Bootstrapping triplet", i, "of", nrow(top_boot), ":", triplet_name, "\n")
    
    df <- data.frame(
        A = log(as.numeric(data_filtered[A, ])),
        BDivC = log(as.numeric(data_filtered[B, ]) / as.numeric(data_filtered[C, ])),
        group = group_factor
    )
    
    slopes_early <- c(); slopes_late <- c()
    for (j in 1:num_bootstrap) {
        idx <- sample(seq_len(nrow(df)), replace = TRUE)
        d <- df[idx, ]
        m <- tryCatch(lm(A ~ BDivC * group, data = d), error = function(e) NULL)
        if (!is.null(m)) {
            coefs <- coef(m)
            if (all(c("BDivC", "BDivC:grouplate") %in% names(coefs))) {
                slopes_early <- c(slopes_early, coefs["BDivC"])
                slopes_late  <- c(slopes_late, coefs["BDivC"] + coefs["BDivC:grouplate"])
            }
        }
    }
    
    if (length(slopes_early) < 10 || length(slopes_late) < 10) {
        cat("⚠️ Skipping triplet:", triplet_name, "- too few valid slopes\n")
        next
    }
    
    slope_diff <- slopes_late - slopes_early
    ci <- quantile(slope_diff, c(0.025, 0.975), na.rm = TRUE)
    agree <- ifelse((row$qval_interaction < fdr_qval && ci[1] > 0) ||
                        (row$qval_interaction >= fdr_qval && ci[1] <= 0 && ci[2] >= 0), "yes", "no")
    
    all_boots[[triplet_name]] <- list(
        early = slopes_early,
        late = slopes_late,
        diff = slope_diff
    )
    
    bootstrap_rows[[length(bootstrap_rows) + 1]] <- data.frame(
        A = A, B = B, C = C,
        bootstrap_diff = mean(slope_diff, na.rm = TRUE),
        ci_lower = ci[1],
        ci_upper = ci[2],
        agree_with_qval = agree
    )
}

if (length(all_boots) == 0) {
    stop("❌ No bootstrap results were saved. All triplets failed or were skipped.")
}

# Save all bootstraps to single .rds
saveRDS(all_boots, file = file.path(output_dir, "bootstrap_results_all.rds"))

# Save summary CSV
bootstrap_summary_df <- do.call(rbind, bootstrap_rows)
write.csv(bootstrap_summary_df, file.path(output_dir, "bootstrap_summary_top.csv"), row.names = FALSE)

cat("\n✅ Bootstrapping complete. Saved", length(all_boots), "triplets.\n")
cat("📁 Output directory:", output_dir, "\n")

# ---- Generate Final Plots (Scatter + Bootstrap Density) ----

library(gridExtra)

plot_dir <- file.path(output_dir, "final_triplet_plots")
dir.create(plot_dir, showWarnings = FALSE)

cat("\n📊 Generating final plots with scatter + density for", length(all_boots), "triplets\n")

for (name in names(all_boots)) {
    A_B_C <- unlist(strsplit(name, "_"))
    A <- A_B_C[1]; B <- A_B_C[2]; C <- A_B_C[3]
    
    df <- data.frame(
        A = log(as.numeric(data_filtered[A, ])),
        BDivC = log(as.numeric(data_filtered[B, ]) / as.numeric(data_filtered[C, ])),
        group = group_factor
    )
    
    boot <- all_boots[[name]]
    
    # Left Panel: Scatter with regression lines
    p1 <- ggplot(df, aes(x = BDivC, y = A, color = group)) +
        geom_point(size = 2, alpha = 0.7) +
        geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
        labs(
            title = "Scatter Plot with Fitted Lines",
            x = paste0("log(", B, "/", C, ")"),
            y = paste0("log(", A, ")")
        ) +
        scale_color_manual(values = c("early" = "steelblue", "late" = "tomato")) +
        theme_minimal(base_size = 14)
    
    # Right Panel: Density of bootstrapped slopes
    df_density <- data.frame(
        Slope = c(boot$early, boot$late),
        Group = rep(c("Early", "Late"), times = c(length(boot$early), length(boot$late)))
    )
    
    p2 <- ggplot(df_density, aes(x = Slope, fill = Group)) +
        geom_density(alpha = 0.5) +
        labs(
            title = "Bootstrapped Slope Distributions",
            x = "Slope",
            y = "Density"
        ) +
        scale_fill_manual(values = c("Early" = "steelblue", "Late" = "tomato")) +
        theme_minimal(base_size = 14)
    
    combined <- grid.arrange(p1, p2, ncol = 2,
                             top = paste("Triplet:", A, "~ log(", B, "/", C, ")")
    )
    
    ggsave(
        filename = file.path(plot_dir, paste0("triplet_", name, ".png")),
        plot = combined,
        width = 12,
        height = 5
    )
}

cat("\n✅ Final scatter + density plots saved to:", plot_dir, "\n")

