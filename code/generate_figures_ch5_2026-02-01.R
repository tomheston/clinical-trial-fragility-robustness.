# =============================================================================
# CHAPTER 5 — FIGURE GENERATION (robust, self-contained)
# Produces Figures 5.1, 5.2, 5.3 and saves PNG + PDF to an /out_figures folder
# Input CSVs expected in the same folder (or FRAGILITY_OUT):
#   - roc_curve_mfq_sigReq.csv
#   - correlations_by_scenario.csv
# =============================================================================

options(warn = 1)

# ---- Package setup (loads or installs if missing) ---------------------------
ensure_pkg <- function(pkg) {
  if (!suppressWarnings(require(pkg, character.only = TRUE))) {
    message(sprintf("Package '%s' not found; attempting install…", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
    if (!suppressWarnings(require(pkg, character.only = TRUE))) {
      stop(sprintf("Could not load '%s'. Please install it and re-run.", pkg))
    }
  }
}

pkgs <- c("ggplot2","dplyr","tidyr","scales")
invisible(lapply(pkgs, ensure_pkg))

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# ---- Robust path discovery (drop-in replacement for setwd/output_dir) -------
needed_files <- c("roc_curve_mfq_sigReq.csv", "correlations_by_scenario.csv")

cand_roots <- c(
  Sys.getenv("FRAGILITY_OUT", unset = NA_character_),
  "~/out",
  "~/out_random_v2",
  "~/Documents/out",
  getwd()
) |> path.expand() |> unique()

exists_all <- function(root) all(file.exists(file.path(root, needed_files)))

input_dir <- NA_character_
for (r in cand_roots) {
  if (!is.na(r) && dir.exists(r) && exists_all(r)) { input_dir <- r; break }
}

if (is.na(input_dir)) {
  # Fallback: recursive search under Documents (Windows/macOS safe)
  docs <- path.expand("~/Documents")
  hits <- tryCatch(
    list.files(docs, pattern = paste(needed_files, collapse = "|"),
               recursive = TRUE, full.names = TRUE),
    error = function(e) character(0)
  )
  if (length(hits) > 0) {
    dirs <- dirname(hits)
    tab  <- table(dirs)
    best <- names(tab[tab >= length(needed_files)])[1]
    if (!is.na(best) && exists_all(best)) input_dir <- best
  }
}

if (is.na(input_dir) || !exists_all(input_dir)) {
  cat("\n*** Could not find required inputs. Tried these roots:\n")
  print(cand_roots)
  cat("\nMake sure these files live in the SAME folder (or set FRAGILITY_OUT):\n",
      paste0(" - ", needed_files, collapse = "\n"), "\n\n")
  stop("Missing expected input file(s).")
}

output_dir <- file.path(input_dir, "out_figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

infile  <- function(name) file.path(input_dir,  name)
outfile <- function(name) file.path(output_dir, name)

cat("Figure input dir : ", input_dir,  "\n")
cat("Figure output dir: ", output_dir, "\n\n")

# ---- Global theme -----------------------------------------------------------
theme_dissertation <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90"),
      plot.title       = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle    = element_text(size = 10, hjust = 0.5, color = "gray40"),
      axis.title       = element_text(face = "bold"),
      legend.position  = "right",
      legend.background= element_rect(color = "black", linewidth = 0.3),
      plot.margin      = margin(10, 10, 10, 10)
    )
}

## ===== Figure 5.1: MFQ ROC (Zoomed, no AUC annotation) =====

# --- paths (edit if needed) ---
csv_dir <- "~/out"          # folder that contains 'roc_curve_mfq_sigReq.csv'
fig_dir <- "~/out_figures"  # where to save the figure
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# --- minimal deps: install if missing, then load ---
need <- c("ggplot2","dplyr","scales")
to_install <- need[!need %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
lapply(need, library, character.only = TRUE)

# --- load data ---
roc_path <- file.path(csv_dir, "roc_curve_mfq_sigReq.csv")
roc <- read.csv(roc_path)

# ensure sorted & numeric
roc <- roc %>% mutate(FPR = as.numeric(FPR), TPR = as.numeric(TPR)) %>% arrange(FPR)

# find the MFQ = 0.05 point (fall back to nearest if exact not present)
i <- which.min(abs(roc$threshold - 0.05))
thr <- roc[i, , drop = FALSE]

# zoom region (keeps only the informative corner)
roc_zoom <- roc %>% filter(FPR <= 0.05, TPR <= 0.20)

# label text for the 0.05 point
lbl <- sprintf("MFQ = 0.05\nTPR = %.1f%%\nFPR = %.1f%%\nSpec = %.1f%%",
               100*thr$TPR, 100*thr$FPR, 100*(1 - thr$FPR))

# --- plot ---
p <- ggplot(roc_zoom, aes(x = FPR, y = TPR)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray60") +
  geom_line(linewidth = 1) +
  geom_point(data = thr, aes(FPR, TPR), size = 2.8) +
  geom_label(
    data = thr,
    aes(FPR, TPR, label = lbl),
    hjust = -0.05, vjust = 1.1, label.size = 0.25, label.r = unit(0.15,"lines")
  ) +
  coord_cartesian(xlim = c(0, 0.05), ylim = c(0, 0.20), expand = FALSE) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.5)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.5)) +
  labs(
    title = "Figure 5.1: MFQ ROC Curve for Fragility Detection (Zoomed)",
    subtitle = "Operating characteristics at varying MFQ thresholds (n=720,000 simulated trials)",
    x = "False Positive Rate",
    y = "True Positive Rate (Sensitivity)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40"),
    axis.title = element_text(face = "bold")
  )

# --- save ---
ggsave(file.path(fig_dir, "Figure_5_1_MFQ_ROC_Curve_Zoomed.png"), p,
       width = 8, height = 6, dpi = 600, bg = "white")
ggsave(file.path(fig_dir, "Figure_5_1_MFQ_ROC_Curve_Zoomed.pdf"), p,
       width = 8, height = 6, device = "pdf")

p  # prints in the console if you're in RStudio


# =============================================================================
# FIGURE 5.2 — SFI–RQ CORRELATION HEATMAP (RR=0.6 slice)
# =============================================================================
cat("Generating Figure 5.2: SFI–RQ Correlation Heatmap…\n")

cor_data <- read.csv(infile("correlations_by_scenario.csv"))

heatmap_data <- cor_data %>%
  dplyr::filter(RR_true == 0.6) %>%
  mutate(
    N_label   = paste0("N=", N),
    p2_label  = paste0("p\u2082=", p2),   # p₂
    rho_display = round(rho_SFI_RQ, 2)
  )

fig5_2 <- ggplot(heatmap_data, aes(x = p2_label, y = N_label, fill = rho_SFI_RQ)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", rho_display)),
            size = 3, fontface = "bold") +
  scale_fill_gradient2(
    low = "#0072B2", mid = "white", high = "#D55E00",
    midpoint = 0, limits = c(-1, 1),
    breaks = seq(-1, 1, 0.25),
    labels = sprintf("%.2f", seq(-1, 1, 0.25)),
    name = "Spearman \u03C1"
  ) +
  facet_wrap(~ alloc, ncol = 3,
             labeller = labeller(alloc = function(x) paste("Allocation:", x))) +
  labs(
    title    = "Figure 5.2: SFI–RQ Correlation Across Design Parameters",
    subtitle = "Heatmap shows context-dependent relationship (RR = 0.6 scenarios)",
    x = "Control Event Rate",
    y = "Sample Size"
  ) +
  theme_dissertation() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text       = element_text(face = "bold"),
    legend.position  = "bottom",
    legend.key.width = unit(2, "cm")
  )

ggsave(outfile("Figure_5_2_SFI_RQ_Correlation_Heatmap.png"),
       fig5_2, width = 10, height = 8, dpi = 600, bg = "white")
ggsave(outfile("Figure_5_2_SFI_RQ_Correlation_Heatmap.pdf"),
       fig5_2, width = 10, height = 8, device = "pdf")

cat("  Figure 5.2 saved.\n\n")

# =============================================================================
# FIGURE 5.3 — PATTERN (1,1,0) COMPONENT ANALYSIS (flow)
# =============================================================================
cat("Generating Figure 5.3: Pattern (1,1,0) Component Analysis…\n")

# Counts from null (RR=1.0) 120,000-trial subset
flow_data <- data.frame(
  stage = c("All Null Trials",
            "p ≤ 0.05",
            "p ≤ 0.05 AND\nMFQ ≤ 0.05",
            "Pattern (1,1,0):\np ≤ 0.05 AND\nMFQ ≤ 0.05 AND\nRQ < 0.075"),
  count = c(120000, 4958, 4459, 1603),
  x = c(1, 2, 3, 4),
  y = c(3, 3, 3, 3)
)

flow_data <- flow_data %>%
  mutate(
    pct_of_total = count / 120000 * 100,
    pct_of_prior = c(NA,
                     4958/120000*100,
                     4459/4958*100,
                     1603/4459*100),
    label = sprintf("%s\nn=%s\n%.2f%% of total",
                    stage,
                    format(count, big.mark = ","),
                    pct_of_total)
  )

fig5_3 <- ggplot(flow_data, aes(x = x, y = y)) +
  geom_rect(aes(xmin = x - 0.35, xmax = x + 0.35,
                ymin = y - 0.4, ymax = y + 0.4),
            fill = c("gray90", "#E69F00", "#CC79A7", "#D55E00"),
            color = "black", linewidth = 1) +
  geom_text(aes(label = label), size = 3.5, fontface = "bold", lineheight = 0.9) +
  geom_segment(aes(x = x + 0.35, xend = x + 0.65, y = y, yend = y),
               data = flow_data[1:3, ],
               arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
               linewidth = 1.1, color = "black") +
  annotate("text", x = 1.5, y = 3.5,
           label = sprintf("%.1f%%\nof trials", flow_data$pct_of_prior[2]),
           size = 3, fontface = "italic") +
  annotate("text", x = 2.5, y = 3.5,
           label = sprintf("%.1f%%\nof p≤0.05", flow_data$pct_of_prior[3]),
           size = 3, fontface = "italic") +
  annotate("text", x = 3.5, y = 3.5,
           label = sprintf("%.1f%%\nof fragile", flow_data$pct_of_prior[4]),
           size = 3, fontface = "italic") +
  annotate("rect", xmin = 0.5, xmax = 4.5, ymin = 1.8, ymax = 2.3,
           fill = "lightyellow", color = "black", linewidth = 0.8) +
  annotate("text", x = 2.5, y = 2.05,
           label = "Specificity: 98.66% (1 - 1.34%)\nBaseline Pattern (1,1,0) under null: 1.34%\nPharmaceutical rate: 18.0% (13.4×)",
           size = 3.6, fontface = "bold", lineheight = 1.1) +
  scale_x_continuous(limits = c(0.5, 4.5), expand = c(0, 0)) +
  scale_y_continuous(limits = c(1.5, 4), expand = c(0, 0)) +
  labs(
    title = "Figure 5.3: Pattern (1,1,0) Component Analysis",
    subtitle = "Flow showing progressive filtering to Pattern (1,1,0) in null simulations (120,000 trials)"
  ) +
  theme_void() +
  theme(
    plot.title    = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray40", margin = margin(b = 15)),
    plot.margin   = margin(20, 20, 20, 20)
  )

ggsave(outfile("Figure_5_3_Pattern_110_Component_Analysis.png"),
       fig5_3, width = 12, height = 6, dpi = 600, bg = "white")
ggsave(outfile("Figure_5_3_Pattern_110_Component_Analysis.pdf"),
       fig5_3, width = 12, height = 6, device = "pdf")

cat("  Figure 5.3 saved.\n\n")

# ---- Summary ----------------------------------------------------------------
cat("============================================================================\n")
cat("ALL FIGURES GENERATED\n")
cat("Output location:", output_dir, "\n\n")
cat("Files created:\n")
cat("  - Figure_5_1_MFQ_ROC_Curve_Zoomed.png (and .pdf)\n")
cat("  - Figure_5_2_SFI_RQ_Correlation_Heatmap.png (and .pdf)\n")
cat("  - Figure_5_3_Pattern_110_Component_Analysis.png (and .pdf)\n")
cat("============================================================================\n")

# ---- Optional: display in interactive sessions ------------------------------
if (interactive()) {
  print(fig5_1); print(fig5_2); print(fig5_3)
}
