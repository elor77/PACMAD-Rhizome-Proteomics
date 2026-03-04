# ==============================================================================
# Permutation Test 2 (A + B): Effect-Size Stratified Shuffle
# 10,000 iterations, parallelized via parallel::mclapply
#
# Input:  3a_summary-by-OG_long_log2FC_DAP.csv
# Output: Perm2_AB_10k.pdf   — 2-panel null distribution figure
#         Perm2_AB_10k_stats.csv — observed Δρ, null mean/SD, P-value
#
# M3 note: uses mclapply with detectCores()-1 cores (fork-based, macOS-native)
# ==============================================================================

library(tidyverse)
library(parallel)
library(ggplot2)
library(patchwork)

set.seed(42)
N_PERM  <- 10000
N_CORES <- detectCores() - 1   # leave 1 core free; M3 Pro = 11, M3 Max = 15

output_dir <- "/Users/eo235/Library/CloudStorage/OneDrive-Personal/Projects/BucklerLab/PACMAD-Rhizome-Proteomics/archive_revisions/manuscript_revised/evo_conservation"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==============================================================================
# 1. LOAD & PREPARE DATA
# ==============================================================================

df <- read_csv("3a_summary-by-OG_long_log2FC_DAP.csv", show_col_types = FALSE)

# Combine Trip19/Trip22 → Td: average log2FC, DAP = TRUE if significant in either year
df_td <- df %>%
  mutate(species = if_else(species %in% c("Trip19", "Trip22"), "Td", species)) %>%
  group_by(orthogroup, species) %>%
  summarise(
    log2FC = mean(log2FC, na.rm = TRUE),
    is_dap  = max(if_else(DAP != "NS" & !is.na(DAP), 1L, 0L)),
    .groups = "drop"
  )

# Wide log2FC matrix (rows = OGs, cols = species)
mat_fc <- df_td %>%
  select(orthogroup, species, log2FC) %>%
  pivot_wider(names_from = species, values_from = log2FC)

# Ensure all 5 species columns exist
for (col in c("Ag", "Mg", "Pv", "Sn", "Td")) {
  if (!col %in% names(mat_fc)) mat_fc[[col]] <- NA_real_
}
mat_fc <- mat_fc %>% select(orthogroup, Ag, Mg, Pv, Sn, Td)

# 4-of-5 filter
keep_4of5 <- function(mat) {
  n_present <- rowSums(!is.na(mat %>% select(-orthogroup)))
  mat[n_present >= 4, ]
}
mat_fc <- keep_4of5(mat_fc)
mat_fc <- mat_fc %>% filter(!is.na(orthogroup))
cat("OGs after 4-of-5 filter:", nrow(mat_fc), "\n")   # expect 1683

# DAP count per OG (with Td already collapsed)
dap_counts <- df_td %>%
  filter(orthogroup %in% mat_fc$orthogroup) %>%
  group_by(orthogroup) %>%
  summarise(n_dap = sum(is_dap), .groups = "drop")

# Classify OGs
ogs_shared <- dap_counts %>% filter(n_dap >= 2) %>% pull(orthogroup)
ogs_bg     <- dap_counts %>% filter(n_dap <= 1) %>% pull(orthogroup)
cat("Shared DAPs:", length(ogs_shared), "\n")          # expect 126
cat("Background:", length(ogs_bg), "\n")               # expect 1557

# mean |log2FC| per OG (for stratification)
mean_absfc <- mat_fc %>%
  rowwise() %>%
  mutate(mean_absFC = mean(abs(c_across(Ag:Td)), na.rm = TRUE)) %>%
  ungroup() %>%
  select(orthogroup, mean_absFC)

# ==============================================================================
# 2. HELPER: compute Δρ for a given shared/background split
# ==============================================================================

SPECIES <- c("Ag", "Mg", "Pv", "Sn", "Td")
PAIRS   <- combn(SPECIES, 2, simplify = FALSE)   # 10 pairs

compute_delta_rho <- function(mat, shared_ogs, bg_ogs) {
  rho_shared <- sapply(PAIRS, function(p) {
    sub <- mat %>% filter(orthogroup %in% shared_ogs) %>% select(all_of(p))
    sub <- sub[complete.cases(sub), ]
    if (nrow(sub) < 3) return(NA_real_)
    cor(sub[[1]], sub[[2]], method = "spearman")
  })
  rho_bg <- sapply(PAIRS, function(p) {
    sub <- mat %>% filter(orthogroup %in% bg_ogs) %>% select(all_of(p))
    sub <- sub[complete.cases(sub), ]
    if (nrow(sub) < 3) return(NA_real_)
    cor(sub[[1]], sub[[2]], method = "spearman")
  })
  median(rho_shared, na.rm = TRUE) - median(rho_bg, na.rm = TRUE)
}

observed_delta <- compute_delta_rho(mat_fc, ogs_shared, ogs_bg)
cat("Observed Δρ:", round(observed_delta, 4), "\n")

# ==============================================================================
# 3. PERMUTATION A — naive label shuffle
# ==============================================================================

n_shared <- length(ogs_shared)
all_ogs  <- c(ogs_shared, ogs_bg)

cat("Running Permutation A (naive, 10k)...\n")
null_A <- unlist(mclapply(seq_len(N_PERM), function(i) {
  perm_shared <- sample(all_ogs, n_shared)
  perm_bg     <- setdiff(all_ogs, perm_shared)
  compute_delta_rho(mat_fc, perm_shared, perm_bg)
}, mc.cores = N_CORES))

p_A <- mean(null_A >= observed_delta, na.rm = TRUE)
cat("Perm A: null mean =", round(mean(null_A), 4),
    "| SD =", round(sd(null_A), 4),
    "| P =", p_A, "\n")

# ==============================================================================
# 4. PERMUTATION B — |log2FC|-stratified shuffle
# ==============================================================================

# Assign |FC| quintile bins
mat_fc_annot <- mat_fc %>%
  left_join(mean_absfc, by = "orthogroup") %>%
  mutate(fc_bin = ntile(mean_absFC, 5))

# Label each OG as shared (1) or background (0)
labels_df <- mat_fc_annot %>%
  select(orthogroup, fc_bin) %>%
  mutate(is_shared = orthogroup %in% ogs_shared)

cat("Running Permutation B (|FC|-stratified, 10k)...\n")
null_B <- unlist(mclapply(seq_len(N_PERM), function(i) {
  # shuffle shared/background labels within each |FC| quintile
  perm_labels <- labels_df %>%
    group_by(fc_bin) %>%
    mutate(is_shared = sample(is_shared)) %>%
    ungroup()
  perm_shared <- perm_labels %>% filter(is_shared) %>% pull(orthogroup)
  perm_bg     <- perm_labels %>% filter(!is_shared) %>% pull(orthogroup)
  compute_delta_rho(mat_fc, perm_shared, perm_bg)
}, mc.cores = N_CORES))

p_B <- mean(null_B >= observed_delta, na.rm = TRUE)
cat("Perm B: null mean =", round(mean(null_B), 4),
    "| SD =", round(sd(null_B), 4),
    "| P =", p_B, "\n")

# ==============================================================================
# 5. FIGURES
# ==============================================================================

X_RANGE <- c(-0.25, 0.42)   # shared x-axis across all 4 panels

make_panel <- function(null_vec, obs, p_val, label, subtitle) {
  df_plot <- data.frame(delta_rho = null_vec)
  p_label <- if (p_val == 0) paste0("P < ", 1 / N_PERM) else paste0("P = ", signif(p_val, 3))
  
  ggplot(df_plot, aes(x = delta_rho)) +
    geom_histogram(bins = 80, fill = "#95a5a6", color = "white", linewidth = 0.2) +
    geom_vline(xintercept = obs, color = "#e74c3c", linewidth = 1.2, linetype = "solid") +
    annotate("text", x = obs, y = Inf,
             label = paste0("Observed\nDelta.rho = ", round(obs, 3)),
             hjust = -0.1, vjust = 1.5, color = "#e74c3c", size = 3.5) +
    annotate("text", x = -Inf, y = Inf,
             label = paste0(label, "\n", p_label),
             hjust = -0.1, vjust = 1.5, size = 3.5, fontface = "bold") +
    coord_cartesian(xlim = X_RANGE) +
    labs(x = "Delta.rho (permuted)", y = "Count", subtitle = subtitle) +
    theme_classic(base_size = 11) +
    theme(plot.subtitle = element_text(size = 9, color = "grey40"))
}

p1 <- make_panel(null_A, observed_delta, p_A,
                 "A. Naive shuffle",
                 paste0("Null mean = ", round(mean(null_A), 3),
                        " +/- ", round(sd(null_A), 3)))

p2 <- make_panel(null_B, observed_delta, p_B,
                 "B. |log2FC|-stratified shuffle",
                 paste0("Null mean = ", round(mean(null_B), 3),
                        " +/- ", round(sd(null_B), 3)))

fig <- p1 + p2 +
  plot_annotation(
    title    = "Effect-size control: parallel fold-change response vs. magnitude-matched null",
    subtitle = paste0("Shared DAPs n = ", length(ogs_shared),
                      ", Background n = ", length(ogs_bg),
                      ", ", format(N_PERM, big.mark = ","), " permutations each"),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(file.path(output_dir, "Perm2_AB_10k.pdf"), fig, width = 10, height = 4.5)
cat("Saved: Perm2_AB_10k.pdf\n")

# ==============================================================================
# 6. STATS TABLE
# ==============================================================================

stats_table <- data.frame(
  Test            = c("A. Naive", "B. |FC|-stratified"),
  Observed_delta  = round(observed_delta, 4),
  Null_mean       = round(c(mean(null_A), mean(null_B)), 4),
  Null_SD         = round(c(sd(null_A),   sd(null_B)),   4),
  Residual_delta  = round(observed_delta - c(mean(null_A), mean(null_B)), 4),
  P_value         = c(p_A, p_B),
  N_perm          = N_PERM
)

print(stats_table)
write_csv(stats_table, file.path(output_dir, "Perm2_AB_10k_stats.csv"))
cat("Saved: Perm2_AB_10k_stats.csv\n")