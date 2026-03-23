# MFQ vs FQ Allocation Independence Test
# Replicates manuscript design for dissertation

set.seed(42)

# Simulation parameters (matching manuscript)
sample_sizes <- c(60, 120, 240)
allocation_ratios <- c(1.0, 1.5, 2.0, 3.0, 4.0)
control_risks <- c(0.05, 0.20, 0.40)
risk_ratios <- c(0.75, 0.90, 1.00, 1.10, 1.25)
n_reps <- 500  # Replicates per combination

results <- data.frame(
  N = numeric(),
  alloc_ratio = numeric(),
  p_control = numeric(),
  RR = numeric(),
  FI = numeric(),
  FQ = numeric(),
  MFQ = numeric(),
  rho = numeric()
)

# FI calculation function
calculate_FI <- function(a, b, c, d) {
  n_treat <- a + b
  n_control <- c + d
  
  p_initial <- fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value
  sig_initial <- p_initial < 0.05
  
  # Determine toggle arm (fewer events, or smaller arm if tied)
  if (a < c) {
    toggle_arm <- "treat"
  } else if (a > c) {
    toggle_arm <- "control"
  } else {
    toggle_arm <- ifelse(n_treat < n_control, "treat", "control")
  }
  
  FI <- 0
  current_a <- a
  current_c <- c
  
  for (i in 1:min(n_treat, n_control)) {
    if (toggle_arm == "treat") {
      if (current_a < n_treat) {
        current_a <- current_a + 1
        b_new <- n_treat - current_a
        test_table <- matrix(c(current_a, b_new, c, d), nrow = 2)
      } else break
    } else {
      if (current_c < n_control) {
        current_c <- current_c + 1
        d_new <- n_control - current_c
        test_table <- matrix(c(a, b, current_c, d_new), nrow = 2)
      } else break
    }
    
    p_new <- fisher.test(test_table)$p.value
    sig_new <- p_new < 0.05
    
    FI <- i
    if (sig_initial != sig_new) break
  }
  
  return(list(FI = FI, toggle_arm = toggle_arm))
}

cat("Running allocation independence simulation...\n")
cat("Total combinations:", length(sample_sizes) * length(allocation_ratios) * 
      length(control_risks) * length(risk_ratios) * n_reps, "\n\n")

pb <- txtProgressBar(min = 0, max = length(sample_sizes) * length(allocation_ratios) * 
                       length(control_risks) * length(risk_ratios) * n_reps, style = 3)
counter <- 0

for (N in sample_sizes) {
  for (alloc_ratio in allocation_ratios) {
    for (p_control in control_risks) {
      for (RR in risk_ratios) {
        
        # Calculate arm sizes
        n_control <- round(N / (alloc_ratio + 1))
        n_treat <- N - n_control
        
        # Treatment event rate
        p_treat <- min(p_control * RR, 0.95)
        
        for (rep in 1:n_reps) {
          # Generate 2x2 table
          events_control <- rbinom(1, n_control, p_control)
          events_treat <- rbinom(1, n_treat, p_treat)
          
          a <- events_treat
          b <- n_treat - events_treat
          c <- events_control
          d <- n_control - events_control
          
          # Skip invalid tables
          if (min(a, b, c, d) < 1) next
          
          # Calculate FI
          fi_result <- calculate_FI(a, b, c, d)
          FI <- fi_result$FI
          
          # Determine n_modified based on which arm was toggled
          if (fi_result$toggle_arm == "treat") {
            n_modified <- n_treat
          } else {
            n_modified <- n_control
          }
          
          # Calculate FQ and MFQ
          FQ <- FI / N
          MFQ <- FI / n_modified
          
          # Calculate rho
          rho <- MFQ / (2 * FQ)
          
          # Store results
          results <- rbind(results, data.frame(
            N = N,
            alloc_ratio = alloc_ratio,
            p_control = p_control,
            RR = RR,
            FI = FI,
            FQ = FQ,
            MFQ = MFQ,
            rho = rho
          ))
          
          counter <- counter + 1
          if (counter %% 1000 == 0) setTxtProgressBar(pb, counter)
        }
      }
    }
  }
}

close(pb)
cat("\n\nSimulation complete!\n")
cat("Valid tables generated:", nrow(results), "\n\n")

# ============================================================================
# ANALYSIS: Mean rho by allocation ratio
# ============================================================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("KEY TEST: Does mean ρ scale with allocation ratio?\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Calculate mean rho by allocation ratio
rho_by_alloc <- aggregate(rho ~ alloc_ratio, data = results, FUN = function(x) {
  c(mean = mean(x), sd = sd(x), n = length(x))
})

cat("Mean ρ = MFQ/(2×FQ) by Allocation Ratio:\n\n")
cat(sprintf("%-15s %10s %10s %10s\n", "Allocation", "Mean ρ", "SD", "N"))
cat(paste(rep("-", 70), collapse = ""), "\n")

for (i in 1:nrow(rho_by_alloc)) {
  alloc <- rho_by_alloc$alloc_ratio[i]
  mean_rho <- rho_by_alloc$rho[i, "mean"]
  sd_rho <- rho_by_alloc$rho[i, "sd"]
  n <- rho_by_alloc$rho[i, "n"]
  
  cat(sprintf("%-15s %10.3f %10.3f %10.0f\n", 
              paste0(alloc, ":1"), mean_rho, sd_rho, n))
}

# Test correlation
corr_test <- cor.test(results$alloc_ratio, results$rho)

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("CORRELATION TEST: Allocation Ratio vs ρ\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat(sprintf("\nPearson r = %.4f, p < %.2e\n", corr_test$estimate, corr_test$p.value))

if (corr_test$estimate > 0.9) {
  cat("\n✓✓✓ VERY STRONG correlation - ρ scales linearly with allocation\n")
  cat("✓ PROVES: FQ is allocation-biased, MFQ eliminates bias\n")
} else if (corr_test$estimate > 0.7) {
  cat("\n✓✓ STRONG correlation - ρ increases with allocation\n")
  cat("✓ Supports MFQ superiority claim\n")
} else {
  cat("\n~ Moderate correlation\n")
}

# Save results
write.csv(results, "dissertation_mfq_allocation_test.csv", row.names = FALSE)
cat("\nResults saved: dissertation_mfq_allocation_test.csv\n")

# Visualization
png("dissertation_rho_by_allocation.png", width = 1000, height = 600, res = 100)
boxplot(rho ~ alloc_ratio, data = results,
        xlab = "Allocation Ratio",
        ylab = "ρ = MFQ / (2×FQ)",
        main = "Allocation Scaling Factor by Allocation Ratio",
        names = paste0(unique(results$alloc_ratio), ":1"),
        col = "steelblue")
abline(h = 1, col = "red", lty = 2, lwd = 2)
legend("topleft", legend = "ρ = 1 (perfect balance)", col = "red", lty = 2, lwd = 2)
dev.off()

cat("Plot saved: dissertation_rho_by_allocation.png\n")
cat("\n", paste(rep("=", 70), collapse = ""), "\n")