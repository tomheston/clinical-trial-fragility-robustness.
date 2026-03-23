# RQ Simulation - Corrected Formula (Base R only, NA handling)
# Establishes empirical cutoffs for RQ based on realistic trial distributions

set.seed(42)

n_sims <- 50000 # 1 million simulations

results <- data.frame(
  p  = numeric(n_sims),
  RQ = numeric(n_sims)
)

valid <- 0
cat("Starting RQ simulation...\n")
pb <- txtProgressBar(min = 0, max = n_sims, style = 3)

while (valid < n_sims) {
  
  # Generate realistic trial parameters
  n_control <- sample(30:1200, 1)
  n_treat   <- sample(30:1200, 1)
  p_control <- runif(1, 0.03, 0.60)
  true_rr   <- runif(1, 0.20, 3.0)   # generates cutoffs for RR=1 only
  true_rr   <- 1.0 # generates cutoffs for RR=1 only
  
  # Cap probability at 1.0 to avoid NA in rbinom
  p_treat <- min(p_control * true_rr, 1.0)
  
  # Generate 2x2 table
  events_control <- rbinom(1, n_control, p_control)
  events_treat   <- rbinom(1, n_treat, p_treat)
  
  a <- events_treat
  b <- n_treat - events_treat
  c <- events_control
  d <- n_control - events_control
  
  # Validity checks - handle NAs explicitly
  if (any(is.na(c(a, b, c, d)))) next
  if (min(a, b, c, d) < 0) next
  if (a+b == 0 || c+d == 0 || a+c == 0 || b+d == 0) next
  
  # Calculate p-value
  ft <- tryCatch(
    fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "two.sided"),
    error = function(e) NULL
  )
  if (is.null(ft)) next
  
  p <- ft$p.value
  
  # CORRECT RQ CALCULATION (FRAGILITY_METRICS.md v10.2.2)
  total <- a + b + c + d
  
  # Expected values under independence
  expected_a <- (a + b) * (a + c) / total
  expected_b <- (a + b) * (b + d) / total
  expected_c <- (c + d) * (a + c) / total
  expected_d <- (c + d) * (b + d) / total
  
  # RRI = (1/k) * sum(|O - E|) where k=4 for 2x2 table
  rri <- (abs(a - expected_a) + abs(b - expected_b) + 
            abs(c - expected_c) + abs(d - expected_d)) / 4
  
  # RQ = RRI / (N/k)
  rq <- rri / (total / 4)
  
  # Validity check for RQ
  if (is.na(rq) || rq < 0 || rq > 1) next
  
  # Store results
  valid <- valid + 1
  results$p[valid] <- p
  results$RQ[valid] <- rq
  
  if (valid %% 10000 == 0) setTxtProgressBar(pb, valid)
}

close(pb)
cat("\nSimulation complete! Generated", valid, "valid tables.\n\n")

# Calculate reference percentiles
cat("RQ Percentiles:\n")
rq_percentiles <- quantile(results$RQ, 
                           probs = c(0.01, 0.05, 0.10, 0.25, 0.33, 0.50, 
                                     0.67, 0.75, 0.90, 0.95, 0.99))
print(round(rq_percentiles, 4))

# Summary statistics
cat("\nRQ Summary Statistics:\n")
cat("Mean:  ", round(mean(results$RQ), 4), "\n")
cat("Median:", round(median(results$RQ), 4), "\n")
cat("SD:    ", round(sd(results$RQ), 4), "\n")

# Proposed cutoffs based on simulation
cat("\nProposed Empirical Cutoffs:\n")
cat("Low robustness (near neutrality):  RQ < ", round(rq_percentiles["33%"], 3), "\n")
cat("Medium robustness:                 RQ ", round(rq_percentiles["33%"], 3), 
    " - ", round(rq_percentiles["67%"], 3), "\n")
cat("High robustness (far from null):   RQ > ", round(rq_percentiles["67%"], 3), "\n")

# Save results
write.csv(results, "rq_simulation_results.csv", row.names = FALSE)
cat("\nResults saved as 'rq_simulation_results.csv' in", getwd(), "\n")

# Create distribution plot (base R)
png("rq_distribution.png", width = 1000, height = 600, res = 100)
hist(results$RQ, breaks = 50, col = "steelblue", 
     main = "RQ Distribution from 1 Million Simulated 2x2 Tables",
     xlab = "RQ (Risk Quotient)", 
     ylab = "Frequency")
abline(v = rq_percentiles[c("33%", "67%")], col = "red", lty = 2, lwd = 2)
legend("topright", legend = "33rd & 67th percentiles", col = "red", lty = 2, lwd = 2)
dev.off()
cat("Distribution plot saved as 'rq_distribution.png'\n")