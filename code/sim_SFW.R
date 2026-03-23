# fragility_sims_v3.R
# Simulation pipeline for Chapter 5 (§5.5–§5.6) — SFI, PFI, MFQ, RQ + Pattern(110)=SFW
# Canonical version aligned with Python toolkit v.02-DEC-2025-GREEDY-VERIFICATION
#
# Adds: CONDITIONAL SFW rates among significant results by RR_true (and by scenario).
#
# Outputs (in OUT_DIR):
# - replicates_raw.csv
# - operating_characteristics.csv
# - table_5_9_property_checks.csv
# - correlations_by_scenario.csv
# - quartile_bin_probs.csv (if empirical_quartiles.csv exists)
# - thresholds_used.csv
# - fpr_overall_pattern110.csv
# - fpr_by_design_pattern110.csv
# - perf_mfq_*.csv, perf_rq_*.csv
# - roc_curve_mfq_sigReq.csv
# - sfw_unconditional_by_rr.csv
# - sfw_conditional_by_rr.csv
# - sfw_conditional_by_scenario.csv

# --------- User-configurable knobs ---------
SEED          <- as.integer(Sys.getenv("FRAGILITY_SEED",  "20251203"))
ALPHA         <- as.numeric(Sys.getenv("FRAGILITY_ALPHA", "0.05"))
REPS_PER_CELL <- as.integer(Sys.getenv("FRAGILITY_REPS",  "2000"))
OUT_DIR       <- path.expand(Sys.getenv("FRAGILITY_OUT",  "out"))
QFILE         <- file.path(OUT_DIR, "empirical_quartiles.csv")

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)
set.seed(SEED)

# Toolkit thresholds (match Python)
N_THRESHOLD        <- as.integer(Sys.getenv("FRAGILITY_N_THRESHOLD",        "50"))
MIN_CELL_THRESHOLD <- as.integer(Sys.getenv("FRAGILITY_MIN_CELL_THRESHOLD", "5"))

# --------- Design grid (edit if needed) ---------
N_TOTAL   <- c(60, 100, 200, 400, 800)
ALLOC     <- list(c(1,1), c(2,1), c(3,2))   # n1:n2 ratios
P_CTRL    <- c(0.05, 0.10, 0.20, 0.40)      # control risk
RR_TRUE <- c(0.60, 0.70, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.30, 1.40) # benefit (<1) and harm (>1) magnitudes



# --------- Core utilities ---------
n_total <- function(a,b,c,d) a + b + c + d

is_significant <- function(p, alpha = ALPHA){
  is.finite(p) && !is.na(p) && p <= alpha
}

test_p <- function(a,b,c,d, use_chi2 = FALSE){
  mat <- matrix(c(a,b,c,d), nrow = 2, byrow = TRUE)
  if (use_chi2){
    ct <- suppressWarnings(stats::chisq.test(mat, correct = FALSE))
    return(ct$p.value)
  } else {
    ft <- suppressWarnings(stats::fisher.test(mat))
    return(ft$p.value)
  }
}

# --------- FI / FQ / MFQ (Walsh-compliant) ---------
choose_arm <- function(a,b,c,d){
  # Arm with fewer events; if tie, fewer total; if still tie, A
  if (a < c) return("A")
  if (c < a) return("B")
  totA <- a + b; totB <- c + d
  if (totA < totB) return("A")
  if (totB < totA) return("B")
  "A"
}

toggle_once <- function(a,b,c,d, arm, direction){
  # direction: "up" (non->event), "down" (event->non) within arm
  if (arm == "A"){
    if (direction == "up"){
      if (b <= 0) return(NULL)
      return(list(a = a + 1L, b = b - 1L, c = c, d = d))
    } else {
      if (a <= 0) return(NULL)
      return(list(a = a - 1L, b = b + 1L, c = c, d = d))
    }
  } else { # arm B
    if (direction == "up"){
      if (d <= 0) return(NULL)
      return(list(a = a, b = b, c = c + 1L, d = d - 1L))
    } else {
      if (c <= 0) return(NULL)
      return(list(a = a, b = b, c = c - 1L, d = d + 1L))
    }
  }
}

steps_to_cross <- function(a,b,c,d, arm, direction,
                           use_chi2 = FALSE,
                           alpha = ALPHA,
                           max_iter = 1e6){
  base_p   <- test_p(a,b,c,d,use_chi2)
  base_sig <- is_significant(base_p, alpha)

  steps <- 0L
  A <- a; B <- b; C <- c; D <- d

  for (i in seq_len(max_iter)){
    nxt <- toggle_once(A,B,C,D, arm, direction)
    if (is.null(nxt)){
      last_p <- test_p(A,B,C,D,use_chi2)
      return(list(steps = NA_integer_, p = last_p))
    }
    A <- nxt$a; B <- nxt$b; C <- nxt$c; D <- nxt$d
    steps <- steps + 1L
    p_now <- test_p(A,B,C,D,use_chi2)
    if (is_significant(p_now, alpha) != base_sig){
      return(list(steps = steps, p = p_now))
    }
  }
  list(steps = NA_integer_, p = test_p(A,B,C,D,use_chi2))
}

compute_fi_fq_mfq <- function(a,b,c,d, alpha = ALPHA){
  N <- n_total(a,b,c,d)
  use_chi2 <- (N > N_THRESHOLD) && (min(a,b,c,d) >= MIN_CELL_THRESHOLD)

  base_p <- test_p(a,b,c,d,use_chi2)
  base_state <- if (is_significant(base_p, alpha)) "significant" else "non-significant"
  arm <- choose_arm(a,b,c,d)

  up_res   <- steps_to_cross(a,b,c,d, arm, "up",   use_chi2, alpha)
  down_res <- steps_to_cross(a,b,c,d, arm, "down", use_chi2, alpha)

  candidates <- list()
  if (!is.na(up_res$steps))   candidates[[length(candidates)+1]] <- list(direction = "non-events → events",
                                                                         steps = up_res$steps,
                                                                         p = up_res$p)
  if (!is.na(down_res$steps)) candidates[[length(candidates)+1]] <- list(direction = "events → non-events",
                                                                         steps = down_res$steps,
                                                                         p = down_res$p)

  n_mod <- if (arm == "A") a + b else c + d

  if (length(candidates) == 0L){
    return(list(
      FI = NA_integer_, FQ = NA_real_, MFQ = NA_real_,
      baseline_p = base_p, baseline_state = base_state,
      arm = arm, direction = NA_character_, n_mod = n_mod,
      final_p = base_p,
      test_used = if (use_chi2) "chi2" else "fisher_exact",
      note = "FI not attainable under Walsh constraints."
    ))
  }

  idx <- which.min(vapply(candidates, `[[`, integer(1), "steps"))
  best <- candidates[[idx]]

  FI  <- as.integer(best$steps)
  FQ  <- if (N > 0) FI / N else NA_real_
  MFQ <- if (n_mod > 0) FI / n_mod else NA_real_

  list(
    FI = FI, FQ = FQ, MFQ = MFQ,
    baseline_p = base_p, baseline_state = base_state,
    target_state = if (base_state == "significant") "non-significant" else "significant",
    arm = arm, direction = best$direction, n_mod = n_mod,
    final_p = best$p,
    test_used = if (use_chi2) "chi2" else "fisher_exact"
  )
}

# --------- SFI (larger arm, BFU = 1/n_large) ---------
compute_sfi <- function(a,b,c,d, alpha = ALPHA){
  N <- n_total(a,b,c,d)
  use_chi2 <- (N > N_THRESHOLD) && (min(a,b,c,d) >= MIN_CELL_THRESHOLD)

  n_A <- a + b
  n_B <- c + d

  if (n_A > n_B){
    larger_arm <- "A"; n_large <- n_A
  } else if (n_B > n_A){
    larger_arm <- "B"; n_large <- n_B
  } else {
    if (a < c){
      larger_arm <- "A"; n_large <- n_A
    } else {
      larger_arm <- "B"; n_large <- n_B
    }
  }

  if (n_large <= 0){
    return(list(
      SFI = NA_integer_, BFU = NA_real_,
      baseline_p = test_p(a,b,c,d,use_chi2),
      baseline_state = NA_character_,
      larger_arm = NA_character_, n_large = 0L,
      final_p = NA_real_,
      test_used = if (use_chi2) "chi2" else "fisher_exact",
      note = "Empty arm."
    ))
  }

  BFU <- 1.0 / n_large

  base_p <- test_p(a,b,c,d,use_chi2)
  base_state <- if (is_significant(base_p, alpha)) "significant" else "non-significant"

  up_res   <- steps_to_cross(a,b,c,d, larger_arm, "up",   use_chi2, alpha)
  down_res <- steps_to_cross(a,b,c,d, larger_arm, "down", use_chi2, alpha)

  candidates <- list()
  if (!is.na(up_res$steps))   candidates[[length(candidates)+1]] <- list(direction = "non-events → events",
                                                                         steps = up_res$steps,
                                                                         p = up_res$p)
  if (!is.na(down_res$steps)) candidates[[length(candidates)+1]] <- list(direction = "events → non-events",
                                                                         steps = down_res$steps,
                                                                         p = down_res$p)

  if (length(candidates) == 0L){
    return(list(
      SFI = NA_integer_, BFU = BFU,
      baseline_p = base_p, baseline_state = base_state,
      larger_arm = larger_arm, n_large = n_large,
      final_p = base_p,
      test_used = if (use_chi2) "chi2" else "fisher_exact",
      note = "SFI not attainable (larger arm cannot flip significance)."
    ))
  }

  idx <- which.min(vapply(candidates, `[[`, integer(1), "steps"))
  best <- candidates[[idx]]

  list(
    SFI = as.integer(best$steps),
    BFU = BFU,
    baseline_p = base_p, baseline_state = base_state,
    target_state = if (base_state == "significant") "non-significant" else "significant",
    larger_arm = larger_arm, n_large = n_large,
    direction = best$direction, final_p = best$p,
    test_used = if (use_chi2) "chi2" else "fisher_exact"
  )
}

# --------- Pearson χ² + PFI along fixed-margin path ---------
chisq_p_pearson <- function(a,b,c,d){
  tbl <- matrix(c(a,b,c,d), nrow = 2, byrow = TRUE)
  if (any(tbl < 0) || sum(tbl) == 0){
    return(c(chi2 = NA_real_, p = NA_real_))
  }
  ct <- suppressWarnings(stats::chisq.test(tbl, correct = FALSE))
  c(chi2 = unname(ct$statistic), p = ct$p.value)
}

p_value_pearson <- function(a,b,c,d){
  chisq_p_pearson(a,b,c,d)["p"]
}

p_along_x <- function(a,b,c,d,x){
  aa <- a + x
  bb <- b - x
  cc <- c - x
  dd <- d + x
  if (min(aa,bb,cc,dd) < 0) return(NA_real_)
  p_value_pearson(aa,bb,cc,dd)
}

compute_pfi_along_path <- function(a,b,c,d, p_func = p_along_x,
                                   alpha = ALPHA, tol = 1e-12, max_iter = 300){
  n <- as.numeric(a + b + c + d)
  if (n <= 0 || min(a,b,c,d) < 0) return(list(UFI_cont = NA_real_, PFI = NA_real_,
                                              p_flip = NA_real_, boundary_limited = FALSE))

  x_min <- -min(a,d)
  x_max <-  min(b,c)

  p0 <- p_func(a,b,c,d,0)
  if (!is.finite(p0)) return(list(UFI_cont = NA_real_, PFI = NA_real_,
                                  p_flip = NA_real_, boundary_limited = FALSE))
  sig0 <- (p0 <= alpha)

  f <- function(x) p_func(a,b,c,d,x) - alpha

  # independence / neutral point
  r1 <- a + b; r2 <- c + d
  c1 <- a + c; c2 <- b + d
  a_exp <- r1 * c1 / n
  x_neutral <- a_exp - a

  if (sig0){
    if (x_neutral >= 0){
      lo <- 0.0; hi <- x_neutral
    } else {
      lo <- x_neutral; hi <- 0.0
    }
  } else {
    if (x_neutral >= 0){
      lo <- 0.0; hi <- x_min
    } else {
      lo <- 0.0; hi <- x_max
    }
  }
  if (lo > hi){
    tmp <- lo; lo <- hi; hi <- tmp
  }
  lo <- max(lo, x_min)
  hi <- min(hi, x_max)
  boundary_x <- if (sig0) x_min else x_max

  # Degenerate -> boundary-limited
  if (abs(lo - hi) < tol){
    pfi_raw <- 4.0 * abs(boundary_x) / n
    pfi <- max(0.0, min(1.0, pfi_raw))
    return(list(UFI_cont = abs(boundary_x), PFI = pfi,
                p_flip = NA_real_, boundary_limited = TRUE))
  }

  flo <- f(lo); fhi <- f(hi)

  if (!is.finite(flo) || !is.finite(fhi) || flo * fhi > 0){
    pfi_raw <- 4.0 * abs(boundary_x) / n
    pfi <- max(0.0, min(1.0, pfi_raw))
    return(list(UFI_cont = abs(boundary_x), PFI = pfi,
                p_flip = NA_real_, boundary_limited = TRUE))
  }

  x_star <- NA_real_
  for (iter in seq_len(max_iter)){
    mid <- 0.5 * (lo + hi)
    fmid <- f(mid)
    if (!is.finite(fmid)){
      mid <- mid + .Machine$double.eps
      fmid <- f(mid)
      if (!is.finite(fmid)) break
    }
    if (abs(hi - lo) < tol || fmid == 0){
      x_star <- mid
      break
    }
    if (flo * fmid <= 0){
      hi <- mid; fhi <- fmid
    } else {
      lo <- mid; flo <- fmid
    }
  }
  if (is.na(x_star)){
    x_star <- 0.5 * (lo + hi)
  }

  p_flip <- p_func(a,b,c,d,x_star)
  pfi_raw <- 4.0 * abs(x_star) / n
  pfi <- max(0.0, min(1.0, pfi_raw))
  boundary_limited <- (abs(abs(x_star) - abs(boundary_x)) < tol)

  list(UFI_cont = abs(x_star), PFI = pfi,
       p_flip = p_flip, boundary_limited = boundary_limited)
}

# --------- RQ via RRI / expected counts ---------
RQ_metric <- function(a,b,c,d){
  N <- a + b + c + d
  if (N <= 0) return(NA_real_)

  n_A <- a + b
  n_B <- c + d
  events <- a + c
  non_events <- b + d

  E_a <- (n_A * events)     / N
  E_b <- (n_A * non_events) / N
  E_c <- (n_B * events)     / N
  E_d <- (n_B * non_events) / N

  k <- 4
  rri <- (1 / k) * (abs(a - E_a) + abs(b - E_b) + abs(c - E_c) + abs(d - E_d))
  rq  <- rri / (N / k)
  rq
}

# --------- Single replicate draw for a scenario ---------
one_rep <- function(n1,n2,p1,p2, alpha = ALPHA){
  a <- stats::rbinom(1, n1, p1)
  c <- stats::rbinom(1, n2, p2)
  b <- n1 - a
  d <- n2 - c

  N <- a + b + c + d
  use_chi2 <- (N > N_THRESHOLD) && (min(a,b,c,d) >= MIN_CELL_THRESHOLD)
  pval <- test_p(a,b,c,d,use_chi2)

  fi_res  <- compute_fi_fq_mfq(a,b,c,d, alpha)
  sfi_res <- compute_sfi(a,b,c,d, alpha)
  pfi_res <- compute_pfi_along_path(a,b,c,d, p_func = p_along_x, alpha = alpha)
  rq      <- RQ_metric(a,b,c,d)

  RR_obs <- if (c == 0 || n2 == 0 || n1 == 0) NA_real_ else (a/n1) / (c/n2)

  data.frame(
    a,b,c,d, n1, n2, p1, p2,
    RR_obs = RR_obs,
    pval = pval,
    SFI = sfi_res$SFI,
    MFQ = fi_res$MFQ,
    PFI = pfi_res$PFI,
    PFI_x = pfi_res$UFI_cont,
    PFI_boundary = pfi_res$boundary_limited,
    RQ = rq,
    stringsAsFactors = FALSE
  )
}

# --------- Scenario runner ---------
run_scenario <- function(Ntot, alloc_vec, p_ctrl, RR, reps = REPS_PER_CELL, alpha = ALPHA){
  r1 <- alloc_vec[1]; r2 <- alloc_vec[2]
  n1 <- round(Ntot * r1/(r1+r2))
  n2 <- Ntot - n1
  p2 <- p_ctrl
  p1 <- p_ctrl * RR
  p1 <- min(max(p1, 0), 1)  # clamp

  out <- do.call(rbind, lapply(seq_len(reps), function(i) one_rep(n1,n2,p1,p2, alpha)))
  out$N <- Ntot
  out$alloc <- sprintf("%d:%d", r1, r2)
  out$RR_true <- RR
  out
}

# --------- Aggregate OC summaries ---------
aggregate_OC <- function(df){
  agg <- function(x){
    if (all(is.na(x))) {
      return(c(mean = NA_real_, median = NA_real_, q25 = NA_real_, q75 = NA_real_))
    }
    c(
      mean   = mean(x, na.rm = TRUE),
      median = stats::median(x, na.rm = TRUE),
      q25    = stats::quantile(x, 0.25, na.rm = TRUE),
      q75    = stats::quantile(x, 0.75, na.rm = TRUE)
    )
  }
  res <- do.call(rbind, lapply(
    split(df, interaction(df$N,df$alloc,df$p2,df$RR_true, drop=TRUE)),
    function(d){
      sSFI <- agg(d$SFI); sPFI <- agg(d$PFI); sMFQ <- agg(d$MFQ); sRQ <- agg(d$RQ)
      data.frame(
        N = d$N[1], alloc = d$alloc[1], p2 = d$p2[1], RR_true = d$RR_true[1],
        SFI_mean = sSFI["mean"], SFI_med = sSFI["median"], SFI_q25 = sSFI["q25"], SFI_q75 = sSFI["q75"],
        PFI_mean = sPFI["mean"], PFI_med = sPFI["median"], PFI_q25 = sPFI["q25"], PFI_q75 = sPFI["q75"],
        MFQ_mean = sMFQ["mean"], MFQ_med = sMFQ["median"], MFQ_q25 = sMFQ["q25"], MFQ_q75 = sMFQ["q75"],
        RQ_mean  = sRQ["mean"],  RQ_med  = sRQ["median"],  RQ_q25  = sRQ["q25"],  RQ_q75  = sRQ["q75"],
        stringsAsFactors = FALSE)
    }))
  rownames(res) <- NULL
  res
}

# --------- Property checks for Table 5.9 ---------
property_checks <- function(){
  checks <- list()

  # SFI allocation invariance: swap arms but keep proportions identical
  a <- 30; b <- 70; c <- 20; d <- 80
  s1 <- compute_sfi(a,b,c,d)$SFI
  s2 <- compute_sfi(c,d,a,b)$SFI
  checks[[length(checks)+1]] <- data.frame(
    property="SFI allocation invariance",
    pass = isTRUE(s1 == s2),
    note=sprintf("SFI %s vs %s on swapped arms", as.character(s1), as.character(s2))
  )

  # MFQ allocation standardization: same proportions, different allocation ratio
  A <- list(a=60,b=40,c=40,d=60)   # 100/100
  B <- list(a=120,b=80,c=40,d=60)  # 200/100
  mfqA <- compute_fi_fq_mfq(A$a,A$b,A$c,A$d)$MFQ
  mfqB <- compute_fi_fq_mfq(B$a,B$b,B$c,B$d)$MFQ
  checks[[length(checks)+1]] <- data.frame(
    property="MFQ allocation standardization",
    pass = isTRUE(abs(mfqA - mfqB) < 1e-9),
    note=sprintf("MFQ %.6f vs %.6f", mfqA, mfqB)
  )

  # RQ bounds & neutrality
  rq0 <- RQ_metric(50,50,50,50)
  rq1 <- RQ_metric(100,0,0,100)
  checks[[length(checks)+1]] <- data.frame(
    property="RQ bounds & neutrality",
    pass = isTRUE((abs(rq0 - 0) < 1e-12) & (abs(rq1 - 1) < 1e-12)),
    note=sprintf("RQ neutrality=%.6f, max=%.6f", rq0, rq1)
  )

  # PFI boundary-limited behavior on perfect concordance (unflippable)
  pfi_bl <- compute_pfi_along_path(100,0,0,100, p_func = p_along_x, alpha = ALPHA)
  checks[[length(checks)+1]] <- data.frame(
    property="PFI boundary-limited case",
    pass = isTRUE(pfi_bl$boundary_limited),
    note=sprintf("PFI=%.6f, boundary=%s", pfi_bl$PFI, as.character(pfi_bl$boundary_limited))
  )

  do.call(rbind, checks)
}

# --------- MAIN: Run sims ---------
message("Running simulation grid…")
all_results <- list()
for (N in N_TOTAL){
  for (al in ALLOC){
    for (p2 in P_CTRL){
      for (RR in RR_TRUE){
        message(sprintf("Scenario N=%d, alloc %s:%s, p2=%.2f, RR_true=%.2f",
                        N, al[1], al[2], p2, RR))
        simdf <- run_scenario(N, al, p2, RR, reps = REPS_PER_CELL, alpha = ALPHA)
        all_results[[length(all_results)+1]] <- simdf
      }
    }
  }
}

sim_df <- do.call(rbind, all_results)
utils::write.csv(sim_df, file.path(OUT_DIR, "replicates_raw.csv"), row.names = FALSE)

# Aggregated operating characteristics
oc <- aggregate_OC(sim_df)
utils::write.csv(oc, file.path(OUT_DIR, "operating_characteristics.csv"), row.names = FALSE)

# Optional quartile bin probabilities
if (file.exists(QFILE)){
  qdf <- utils::read.csv(QFILE, stringsAsFactors = FALSE)
  get_bins <- function(x, m){
    qs <- qdf[qdf$metric==m, c("q1","q2","q3")]
    if (nrow(qs) == 0) return(rep(NA_real_, length(x)))
    cut(x, breaks = c(-Inf, qs$q1, qs$q2, qs$q3, Inf),
        labels = FALSE, right = TRUE)
  }
  sim_df$SFI_bin <- get_bins(sim_df$SFI, "SFI")
  sim_df$PFI_bin <- get_bins(sim_df$PFI, "PFI")
  sim_df$MFQ_bin <- get_bins(sim_df$MFQ, "MFQ")
  sim_df$RQ_bin  <- get_bins(sim_df$RQ,  "RQ")
  bin_probs <- aggregate(
    list(SFI=sim_df$SFI_bin,PFI=sim_df$PFI_bin,MFQ=sim_df$MFQ_bin,RQ=sim_df$RQ_bin),
    by=list(N=sim_df$N,alloc=sim_df$alloc,p2=sim_df$p2,RR_true=sim_df$RR_true),
    FUN=function(z){ mean(z==1, na.rm=TRUE) })
  utils::write.csv(bin_probs, file.path(OUT_DIR, "quartile_bin_probs.csv"), row.names = FALSE)
}

# Property checks table (Table 5.9 content starter)
props <- property_checks()
utils::write.csv(props, file.path(OUT_DIR, "table_5_9_property_checks.csv"), row.names = FALSE)

# Simple within-scenario correlations (diagnostics)
get_corrs <- function(d){
  M <- stats::cor(d[,c("SFI","PFI","MFQ","RQ")],
                  method = "spearman",
                  use = "pairwise.complete.obs")
  data.frame(
    N=d$N[1], alloc=d$alloc[1], p2=d$p2[1], RR_true=d$RR_true[1],
    rho_SFI_PFI=M["SFI","PFI"], rho_SFI_MFQ=M["SFI","MFQ"], rho_SFI_RQ=M["SFI","RQ"],
    rho_PFI_MFQ=M["PFI","MFQ"], rho_PFI_RQ=M["PFI","RQ"], rho_MFQ_RQ=M["MFQ","RQ"],
    stringsAsFactors = FALSE)
}
cr <- do.call(rbind, lapply(
  split(sim_df, interaction(sim_df$N,sim_df$alloc,sim_df$p2,sim_df$RR_true, drop=TRUE)),
  get_corrs))
utils::write.csv(cr, file.path(OUT_DIR, "correlations_by_scenario.csv"), row.names = FALSE)

# ---------- §5.6 add-ons (thresholds + FPR validation + sensitivity/specificity) ----------
null_df <- subset(sim_df, abs(RR_true - 1.0) < 1e-12)

RQ_t1 <- stats::quantile(null_df$RQ, probs = 1/3, na.rm = TRUE)
RQ_t2 <- stats::quantile(null_df$RQ, probs = 2/3, na.rm = TRUE)

thr <- data.frame(metric = c("RQ_t1","RQ_t2","MFQ_cut","alpha"),
                  value  = c(RQ_t1, RQ_t2, 0.05, ALPHA))
utils::write.csv(thr, file.path(OUT_DIR, "thresholds_used.csv"), row.names = FALSE)

sim_df$p_sig       <- sim_df$pval <= ALPHA
sim_df$mfq_fragile <- sim_df$MFQ  <= 0.05
sim_df$rq_low      <- sim_df$RQ   <= RQ_t1
sim_df$pattern110  <- with(sim_df, p_sig & mfq_fragile & rq_low)

null_df$p_sig       <- null_df$pval <= ALPHA
null_df$mfq_fragile <- null_df$MFQ  <= 0.05
null_df$rq_low      <- null_df$RQ   <= RQ_t1
null_df$pattern110  <- with(null_df, p_sig & mfq_fragile & rq_low)

fpr_overall <- mean(null_df$pattern110, na.rm = TRUE)
utils::write.csv(data.frame(FPR_overall_pattern110 = fpr_overall),
                 file.path(OUT_DIR, "fpr_overall_pattern110.csv"), row.names = FALSE)

fpr_by_design <- aggregate(list(FPR = null_df$pattern110),
                           by = list(N = null_df$N, alloc = null_df$alloc, p2 = null_df$p2),
                           FUN = function(z) mean(z, na.rm = TRUE))
utils::write.csv(fpr_by_design, file.path(OUT_DIR, "fpr_by_design_pattern110.csv"), row.names = FALSE)

perf_calc <- function(pos_flag, require_sig = TRUE){
  if (require_sig){
    pos <- pos_flag & (sim_df$pval <= ALPHA)
  } else {
    pos <- pos_flag
  }
  sim_df$POS <- pos
  TPR <- mean(subset(sim_df, abs(RR_true - 1.0) > 1e-12)$POS, na.rm = TRUE)
  FPR <- mean(subset(sim_df, abs(RR_true - 1.0) < 1e-12)$POS, na.rm = TRUE)
  overall <- data.frame(TPR = TPR, FPR = FPR)
  by_scen <- aggregate(list(POS = sim_df$POS),
                       by = list(N = sim_df$N, alloc = sim_df$alloc, p2 = sim_df$p2, RR_true = sim_df$RR_true),
                       FUN = function(z) mean(z, na.rm = TRUE))
  names(by_scen)[names(by_scen)=="POS"] <- if (require_sig) "Rate_sigReq" else "Rate_noSigReq"
  list(overall = overall, by_scen = by_scen)
}

res_mfq_sig   <- perf_calc(sim_df$mfq_fragile, require_sig = TRUE)
res_mfq_nosig <- perf_calc(sim_df$mfq_fragile, require_sig = FALSE)
utils::write.csv(res_mfq_sig$overall,   file.path(OUT_DIR, "perf_mfq_overall_sigReq.csv"), row.names = FALSE)
utils::write.csv(res_mfq_nosig$overall, file.path(OUT_DIR, "perf_mfq_overall_noSigReq.csv"), row.names = FALSE)
by_mfq <- merge(res_mfq_sig$by_scen, res_mfq_nosig$by_scen, by = c("N","alloc","p2","RR_true"))
utils::write.csv(by_mfq, file.path(OUT_DIR, "perf_mfq_by_scenario.csv"), row.names = FALSE)

res_rq_sig   <- perf_calc(sim_df$rq_low, require_sig = TRUE)
res_rq_nosig <- perf_calc(sim_df$rq_low, require_sig = FALSE)
utils::write.csv(res_rq_sig$overall,   file.path(OUT_DIR, "perf_rq_overall_sigReq.csv"), row.names = FALSE)
utils::write.csv(res_rq_nosig$overall, file.path(OUT_DIR, "perf_rq_overall_noSigReq.csv"), row.names = FALSE)
by_rq <- merge(res_rq_sig$by_scen, res_rq_nosig$by_scen, by = c("N","alloc","p2","RR_true"))
utils::write.csv(by_rq, file.path(OUT_DIR, "perf_rq_by_scenario.csv"), row.names = FALSE)

mfq_seq <- seq(0, 0.5, by = 0.005)
roc_rows <- lapply(mfq_seq, function(th){
  nonnull <- subset(sim_df, abs(RR_true - 1.0) > 1e-12)
  null    <- subset(sim_df, abs(RR_true - 1.0) < 1e-12)

  TPR <- mean(nonnull$pval <= ALPHA & nonnull$MFQ <= th, na.rm = TRUE)
  FPR <- mean(null$pval <= ALPHA & null$MFQ <= th, na.rm = TRUE)

  data.frame(threshold = th, TPR = TPR, FPR = FPR)
})
roc_mfq <- do.call(rbind, roc_rows)
utils::write.csv(roc_mfq, file.path(OUT_DIR, "roc_curve_mfq_sigReq.csv"), row.names = FALSE)

# ---------- NEW: SFW (pattern110) unconditional + conditional by RR_true ----------
sfw_uncond_by_rr <- aggregate(
  list(rate = sim_df$pattern110),
  by = list(RR_true = sim_df$RR_true),
  FUN = function(z) mean(z, na.rm = TRUE)
)
utils::write.csv(sfw_uncond_by_rr, file.path(OUT_DIR, "sfw_unconditional_by_rr.csv"), row.names = FALSE)

# Conditional: P(SFW | significant, RR_true)
sfw_cond_by_rr <- do.call(rbind, lapply(split(sim_df, sim_df$RR_true), function(d){
  sig <- d$p_sig
  rate <- if (!any(sig, na.rm = TRUE)) NA_real_ else mean(d$pattern110[sig], na.rm = TRUE)
  data.frame(RR_true = d$RR_true[1], rate = rate)
}))
rownames(sfw_cond_by_rr) <- NULL
utils::write.csv(sfw_cond_by_rr, file.path(OUT_DIR, "sfw_conditional_by_rr.csv"), row.names = FALSE)

# Conditional by full scenario: P(SFW | significant, N, alloc, p2, RR_true)
sfw_cond_by_scenario <- do.call(rbind, lapply(
  split(sim_df, interaction(sim_df$N, sim_df$alloc, sim_df$p2, sim_df$RR_true, drop=TRUE)),
  function(d){
    sig <- d$p_sig
    rate <- if (!any(sig, na.rm = TRUE)) NA_real_ else mean(d$pattern110[sig], na.rm = TRUE)
    data.frame(N=d$N[1], alloc=d$alloc[1], p2=d$p2[1], RR_true=d$RR_true[1], rate_sig = rate,
               n_sig = sum(sig, na.rm = TRUE), n_total = nrow(d))
  }
))
rownames(sfw_cond_by_scenario) <- NULL
utils::write.csv(sfw_cond_by_scenario, file.path(OUT_DIR, "sfw_conditional_by_scenario.csv"), row.names = FALSE)

message(sprintf(
  "Done. Artifacts written to %s :\n - replicates_raw.csv\n - operating_characteristics.csv\n - table_5_9_property_checks.csv\n - correlations_by_scenario.csv\n - quartile_bin_probs.csv (if empirical_quartiles.csv provided)\n - thresholds_used.csv\n - fpr_overall_pattern110.csv\n - fpr_by_design_pattern110.csv\n - perf_mfq_*.csv, perf_rq_*.csv\n - roc_curve_mfq_sigReq.csv\n - sfw_unconditional_by_rr.csv\n - sfw_conditional_by_rr.csv\n - sfw_conditional_by_scenario.csv",
  OUT_DIR
))
