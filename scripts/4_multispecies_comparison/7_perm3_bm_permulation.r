# ==============================================================================
# BM Permulation: Phylogenetic Control of Parallel Fold-Change Analysis
# (Gap 4.3 — Phylogenetic non-independence)
#
# Approach: Per-orthogroup Brownian Motion permulation (Schlötterer/Hsu framework)
#   For each replicate:
#     - For each OG: estimate BM rate from observed log2FC + tree,
#       simulate under BM, rank-match back to original values
#     - Recompute pairwise Spearman rho for shared vs background sets
#     - Test statistic: median rho(shared) - median rho(background)
#
# Two comparison groups:
#   (1) All background (n=1557)
#   (2) Magnitude-matched single-species DAPs (|FC| >= 0.53, n~155)
#
# Input:
#   - 3a_summary-by-OG_long_log2FC_DAP.csv
#   - pruned_pacmad_14tip_rooted.nwk
#
# Output:
#   - Perm3_BM_permulation.pdf       -- null distribution figure (2 panels)
#   - Perm3_BM_permulation_stats.csv -- observed delta, null mean/SD, P-value
#   - Perm3_BM_null_distributions.rds -- raw null vectors for reproducibility
#
# Author: Elad Oren
# Dependencies: tidyverse, ape, geiger, parallel, ggplot2, patchwork
# ==============================================================================

library(tidyverse)
library(ape)
library(geiger)
library(parallel)
library(ggplot2)
library(patchwork)

set.seed(42)
N_PERM  <- 10000
N_CORES <- detectCores() - 1

output_dir <- "/Users/eo235/Library/CloudStorage/OneDrive-Personal/Projects/BucklerLab/PACMAD-Rhizome-Proteomics/archive_revisions/manuscript_revised/evo_conservation"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==============================================================================
# 1. LOAD & PREPARE DATA
# ==============================================================================

df <- read_csv("3a_summary-by-OG_long_log2FC_DAP.csv", show_col_types = FALSE)

# Combine Trip19/Trip22 -> Td: average log2FC, DAP = TRUE if significant in either year
df_td <- df %>%
  mutate(species = if_else(species %in% c("Trip19", "Trip22"), "Td", species)) %>%
  group_by(orthogroup, species) %>%
  summarise(
    log2FC = mean(log2FC, na.rm = TRUE),
    is_dap  = max(if_else(DAP != "NS" & !is.na(DAP), 1L, 0L)),
    .groups = "drop"
  )

# Wide log2FC matrix
mat_fc <- df_td %>%
  select(orthogroup, species, log2FC) %>%
  pivot_wider(names_from = species, values_from = log2FC)

for (col in c("Ag", "Mg", "Pv", "Sn", "Td")) {
  if (!col %in% names(mat_fc)) mat_fc[[col]] <- NA_real_
}
mat_fc <- mat_fc %>% select(orthogroup, Ag, Mg, Pv, Sn, Td)

# 4-of-5 filter
n_present <- rowSums(!is.na(mat_fc %>% select(-orthogroup)))
mat_fc    <- mat_fc[n_present >= 4, ]
cat("OGs after 4-of-5 filter:", nrow(mat_fc), "\n")  # expect 1683

# DAP counts per OG (Td already collapsed)
dap_counts <- df_td %>%
  filter(orthogroup %in% mat_fc$orthogroup) %>%
  group_by(orthogroup) %>%
  summarise(n_dap = sum(is_dap), .groups = "drop")

ogs_shared <- dap_counts %>% filter(n_dap >= 2) %>% pull(orthogroup)
ogs_bg_all <- dap_counts %>% filter(n_dap <= 1) %>% pull(orthogroup)
cat("Shared DAPs:", length(ogs_shared), "\n")   # expect 126
cat("Background:", length(ogs_bg_all), "\n")    # expect 1557

# Magnitude-matched single-species DAPs:
# |FC| threshold = min mean |FC| of shared DAPs
mean_absfc <- mat_fc %>%
  rowwise() %>%
  mutate(mean_absFC = mean(abs(c_across(Ag:Td)), na.rm = TRUE)) %>%
  ungroup() %>%
  select(orthogroup, mean_absFC)

shared_min_fc <- mean_absfc %>%
  filter(orthogroup %in% ogs_shared) %>%
  pull(mean_absFC) %>% min(na.rm = TRUE)
cat("Shared DAP min mean |FC|:", round(shared_min_fc, 3), "\n")  # expect ~0.53

ogs_bg_matched <- dap_counts %>%
  filter(n_dap == 1) %>%           # confirmed single-species DAPs only
  inner_join(mean_absfc, by = "orthogroup") %>%
  filter(mean_absFC >= shared_min_fc) %>%
  pull(orthogroup)
cat("Magnitude-matched single-species DAPs:", length(ogs_bg_matched), "\n")  # expect ~155

# ==============================================================================
# 2. LOAD & PRUNE TREE
# ==============================================================================

# Strip [&R] prefix and internal ASTRAL annotations before parsing
nwk_raw <- readLines("../data/pruned_pacmad_14tip_rooted.nwk")
nwk_clean <- nwk_raw %>%
  paste(collapse = "") %>%
  gsub("^\\[&R\\]\\s*", "", .) %>%                     # remove [&R] prefix
  gsub("'\\[[^]]*\\]'", "", .) %>%                      # remove '[...]' internal labels
  gsub("\\[&[^]]*\\]", "", .)                            # remove any remaining [& ...] annotations

tree_full <- read.tree(text = nwk_clean)

# Tip name -> 2-letter code mapping
tip_map <- c(
  Andropogon_gerardii__TW           = "Ag",
  REF___Miscanthus_sinensis__genome  = "Mg",
  REF___Panicum_virgatum__genome     = "Pv",
  Sorghastrum_fuscescens__PANAND     = "Sn",
  Tripsacum_dactyloides__SRA         = "Td"
)

# Prune to 5 focal tips
tree_full <- ladderize(tree_full)
tips_to_keep <- names(tip_map)
tips_to_drop <- setdiff(tree_full$tip.label, tips_to_keep)
tree5        <- drop.tip(tree_full, tips_to_drop)

# Rename tips
tree5$tip.label <- tip_map[tree5$tip.label]
cat("5-tip tree tips:", paste(tree5$tip.label, collapse = ", "), "\n")
plot(tree5, main = "Pruned 5-tip PACMAD tree"); axisPhylo()

SPECIES <- c("Ag", "Mg", "Pv", "Sn", "Td")
PAIRS   <- combn(SPECIES, 2, simplify = FALSE)  # 10 pairs

# ==============================================================================
# 3. HELPER FUNCTIONS
# ==============================================================================

#' BM rank-match permulation for a single orthogroup
#' Returns a named numeric vector of permulated log2FC values (same names as input)
permulate_og <- function(fc_vec, tree) {
  # fc_vec: named numeric, names = species codes, may contain NAs
  obs <- fc_vec[!is.na(fc_vec)]
  sp  <- names(obs)
  
  if (length(sp) < 3) return(fc_vec)  # not enough tips for BM estimation
  
  # Prune tree to available species for this OG
  tr <- if (length(sp) < length(tree$tip.label)) {
    drop.tip(tree, setdiff(tree$tip.label, sp))
  } else {
    tree
  }
  
  # Estimate BM rate and simulate
  rate <- tryCatch(
    ratematrix(tr, obs[tr$tip.label]),
    error = function(e) NULL
  )
  if (is.null(rate)) return(fc_vec)
  
  sim <- tryCatch(
    sim.char(tr, rate, nsim = 1)[, , 1],
    error = function(e) NULL
  )
  if (is.null(sim)) return(fc_vec)
  
  # Rank-match: reassign observed values by ranks of simulated values
  perm_vals           <- sort(obs)[rank(sim[sp])]
  names(perm_vals)    <- sp
  result              <- fc_vec
  result[sp]          <- perm_vals[sp]
  return(result)
}

#' Compute median Spearman rho across species pairs for a given OG set
median_rho <- function(mat, og_set) {
  sub <- mat %>% filter(orthogroup %in% og_set) %>% select(-orthogroup)
  rhos <- sapply(PAIRS, function(p) {
    d <- sub[, p, drop = FALSE]
    d <- d[complete.cases(d), ]
    if (nrow(d) < 3) return(NA_real_)
    cor(d[[1]], d[[2]], method = "spearman")
  })
  median(rhos, na.rm = TRUE)
}

#' Test statistic: median rho(shared) - median rho(background)
delta_rho <- function(mat, shared, bg) {
  median_rho(mat, shared) - median_rho(mat, bg)
}

# ==============================================================================
# 4. OBSERVED TEST STATISTICS
# ==============================================================================

obs_delta_all     <- delta_rho(mat_fc, ogs_shared, ogs_bg_all)
obs_delta_matched <- delta_rho(mat_fc, ogs_shared, ogs_bg_matched)
cat("Observed delta rho (vs all background):      ", round(obs_delta_all,     4), "\n")
cat("Observed delta rho (vs magnitude-matched):   ", round(obs_delta_matched, 4), "\n")

# ==============================================================================
# 5. BM PERMULATION (10,000 replicates)
# ==============================================================================

# Remove any NA orthogroup rows before building list
mat_fc <- mat_fc %>% filter(!is.na(orthogroup))
cat("OGs after NA orthogroup removal:", nrow(mat_fc), "\n")  # expect 1683

# Recompute group membership on clean OG list
ogs_shared    <- dap_counts %>% filter(n_dap >= 2, orthogroup %in% mat_fc$orthogroup) %>% pull(orthogroup)
ogs_bg_all    <- dap_counts %>% filter(n_dap <= 1, orthogroup %in% mat_fc$orthogroup) %>% pull(orthogroup)
ogs_bg_matched <- dap_counts %>%
  filter(n_dap == 1, orthogroup %in% mat_fc$orthogroup) %>%
  inner_join(mean_absfc, by = "orthogroup") %>%
  filter(mean_absFC >= shared_min_fc) %>%
  pull(orthogroup)

cat("Shared DAPs:", length(ogs_shared), "\n")    # expect 126
cat("Background:", length(ogs_bg_all), "\n")     # expect 1557
cat("Magnitude-matched:", length(ogs_bg_matched), "\n")  # expect ~155

# Pre-extract OG fc vectors as a list for efficiency
og_list <- setNames(
  lapply(seq_len(nrow(mat_fc)), function(i) {
    unlist(mat_fc[i, c("Ag", "Mg", "Pv", "Sn", "Td")], use.names = FALSE)
  }),
  mat_fc$orthogroup
)

cat("Running BM permulation (", format(N_PERM, big.mark = ","), " replicates, ", N_CORES, " cores)...\n", sep = "")

null_results <- mclapply(seq_len(N_PERM), function(i) {
  # Permulate every OG independently
  perm_rows <- lapply(names(og_list), function(og) {
    fc_vec        <- og_list[[og]]
    names(fc_vec) <- SPECIES
    permulate_og(fc_vec, tree5)
  })
  
  # Rebuild matrix
  perm_mat <- as.data.frame(do.call(rbind, perm_rows))
  colnames(perm_mat) <- SPECIES
  perm_mat$orthogroup <- names(og_list)
  
  c(
    delta_all     = delta_rho(perm_mat, ogs_shared, ogs_bg_all),
    delta_matched = delta_rho(perm_mat, ogs_shared, ogs_bg_matched)
  )
}, mc.cores = N_CORES)

null_delta_all     <- sapply(null_results, `[[`, "delta_all")
null_delta_matched <- sapply(null_results, `[[`, "delta_matched")

# P-values
p_all     <- mean(null_delta_all     >= obs_delta_all,     na.rm = TRUE)
p_matched <- mean(null_delta_matched >= obs_delta_matched, na.rm = TRUE)

cat("Vs all background:    null mean =", round(mean(null_delta_all), 4),
    "| SD =", round(sd(null_delta_all), 4), "| P =", p_all, "\n")
cat("Vs magnitude-matched: null mean =", round(mean(null_delta_matched), 4),
    "| SD =", round(sd(null_delta_matched), 4), "| P =", p_matched, "\n")

# Save raw null distributions
saveRDS(
  list(null_all = null_delta_all, null_matched = null_delta_matched),
  file.path(output_dir, "Perm3_BM_null_distributions.rds")
)

# ==============================================================================
# 6. FIGURES
# ==============================================================================

X_RANGE <- c(-0.25, 0.42)   # shared x-axis across all 4 panels

make_panel <- function(null_vec, obs, p_val, n_perm, label, subtitle) {
  df_plot <- data.frame(delta_rho = null_vec)
  p_label <- if (p_val == 0) paste0("P < ", 1 / n_perm) else paste0("P = ", signif(p_val, 3))
  
  ggplot(df_plot, aes(x = delta_rho)) +
    geom_histogram(bins = 80, fill = "#95a5a6", color = "white", linewidth = 0.2) +
    geom_vline(xintercept = obs, color = "#e74c3c", linewidth = 1.2) +
    annotate("text", x = obs, y = Inf,
             label = paste0("Observed\nDelta.rho = ", round(obs, 3)),
             hjust = -0.1, vjust = 1.5, color = "#e74c3c", size = 3.5) +
    annotate("text", x = -Inf, y = Inf,
             label = paste0(label, "\n", p_label),
             hjust = -0.1, vjust = 1.5, size = 3.5, fontface = "bold") +
    coord_cartesian(xlim = X_RANGE) +
    labs(x = "Delta.rho (permulated)", y = "Count", subtitle = subtitle) +
    theme_classic(base_size = 11) +
    theme(plot.subtitle = element_text(size = 9, color = "grey40"))
}

p1 <- make_panel(
  null_delta_all, obs_delta_all, p_all, N_PERM,
  "A. Shared DAPs vs all background",
  paste0("Null mean = ", round(mean(null_delta_all), 3),
         " +/- ", round(sd(null_delta_all), 3),
         "  (n background = ", length(ogs_bg_all), ")")
)

p2 <- make_panel(
  null_delta_matched, obs_delta_matched, p_matched, N_PERM,
  "B. Shared DAPs vs magnitude-matched single-species DAPs",
  paste0("Null mean = ", round(mean(null_delta_matched), 3),
         " +/- ", round(sd(null_delta_matched), 3),
         "  (n matched = ", length(ogs_bg_matched), ")")
)

fig <- p1 + p2 +
  plot_annotation(
    title    = "BM Permulation: Phylogenetic control of parallel fold-change response",
    subtitle = paste0("Shared DAPs n = ", length(ogs_shared),
                      ", ", format(N_PERM, big.mark = ","), " BM permulation replicates",
                      " | Note: null rho > 0 expected under BM (conservative test)"),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(file.path(output_dir, "Perm3_BM_permulation.pdf"), fig, width = 12, height = 5)
cat("Saved: Perm3_BM_permulation.pdf\n")

# ==============================================================================
# 7. STATS TABLE
# ==============================================================================

stats_table <- data.frame(
  Comparison         = c("Shared vs all background", "Shared vs magnitude-matched"),
  N_shared           = length(ogs_shared),
  N_background       = c(length(ogs_bg_all), length(ogs_bg_matched)),
  Observed_delta     = round(c(obs_delta_all, obs_delta_matched), 4),
  Null_mean          = round(c(mean(null_delta_all), mean(null_delta_matched)), 4),
  Null_SD            = round(c(sd(null_delta_all),   sd(null_delta_matched)),   4),
  Residual_delta     = round(c(obs_delta_all - mean(null_delta_all),
                               obs_delta_matched - mean(null_delta_matched)), 4),
  P_value            = c(p_all, p_matched),
  N_perm             = N_PERM
)

print(stats_table)
write_csv(stats_table, file.path(output_dir, "Perm3_BM_permulation_stats.csv"))
cat("Saved: Perm3_BM_permulation_stats.csv\n")
