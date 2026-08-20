# =============================================================================
# simulation.R
#
# Monte Carlo study for the Simplex Kibria-Lukman Estimator (SKLE).
#
# Full factorial design: p in {4, 8, 12, 16}, n in {25, 50, 100, 200},
# rho in {0.80, 0.90, 0.95, 0.99}, sigma^2 in {0.5, 1.0, 1.5, 2.0} and four
# link functions, giving 1,024 scenarios at 1,000 replications each.
#
# The dimension and the number of replications are read from the command
# line, one dimension per run:
#     Rscript simulation.R 8 1000
#
# Equation numbers refer to the manuscript.
# =============================================================================

library(MASS)
library(parallel)

# ==========================================
# 1. CONFIGURATION
# ==========================================
# p and replications read from the command line:
#   Rscript simulation.R 16 1000
.args <- commandArgs(trailingOnly = TRUE)
P_VAL        <- if (length(.args) >= 1) as.integer(.args[1]) else 8
N_VALUES     <- c(25, 50, 100, 200)
RHO_VALUES   <- c(0.80, 0.90, 0.95, 0.99)
SIGMA_VALUES <- c(0.5, 1.0, 1.5, 2.0)
LINKS        <- c("logit", "probit", "cloglog", "neglog")
REPLICATIONS <- if (length(.args) >= 2) as.integer(.args[2]) else 1000

MU_GUARD <- 1e-6
Y_GUARD  <- 1e-6

# Cluster workers are separate R processes and do not inherit the master RNG
# state, so each replication seeds itself from this base.
BASE_SEED <- 2025L

# ==========================================
# 2. FOLDER & FILE SETUP (SMART RESUME)
# ==========================================
# All four dimensions write one file each to the working directory.  Set
# SKLE_OUT to write them elsewhere.
output_dir <- Sys.getenv("SKLE_OUT", unset = ".")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(sprintf(">>> Created folder: %s\n", output_dir))
}

output_file <- file.path(output_dir, paste0("Simulation_Results_p", P_VAL, ".csv"))

if (!file.exists(output_file)) {
  headers <- data.frame(
    Link = character(), Sigma = numeric(), Rho = numeric(), n = numeric(),
    Successful_Reps = numeric(),
    MLE = numeric(),
    SRRE_k1 = numeric(), SKLE_k1 = numeric(),
    SRRE_k2 = numeric(), SKLE_k2 = numeric(),
    SRRE_k3 = numeric(), SKLE_k3 = numeric(),
    SRRE_k4 = numeric(), SKLE_k4 = numeric(),
    # median MSE, robust to the heavy tail at rho = 0.99
    MLE_Med = numeric(),
    SRRE_k1_Med = numeric(), SKLE_k1_Med = numeric(),
    SRRE_k2_Med = numeric(), SKLE_k2_Med = numeric(),
    SRRE_k3_Med = numeric(), SKLE_k3_Med = numeric(),
    SRRE_k4_Med = numeric(), SKLE_k4_Med = numeric(),
    # diagnostics
    Clamp_Rate = numeric(),           # fraction of mu hitting the guard
    Conv_Rate = numeric(),            # fraction of reps that converged
    Mean_k1 = numeric(), Mean_k2 = numeric(), Mean_k3 = numeric(),
    Mean_k4 = numeric(), K4_Binds = numeric(),
    Prov_Hold_k1 = numeric(),         # fraction with k1 < 2*lambda_min (regularity condition)
    Prov_Hold_k2 = numeric(),
    # Equation (33) itself
    # The Eq32_* columns hold the left-hand side of Equation (33), in the
    # scalar form of Equation (36); the column names do not track the
    # equation number.
    Eq32_Hold_k1 = numeric(),         # fraction with Eq (33) satisfied
    Eq32_Hold_k2 = numeric(),
    Eq32_Val_k1 = numeric(),          # median value of the Eq (33) LHS
    Eq32_Val_k2 = numeric(),
    # Eq (33) at the true alpha
    Eq32T_Hold_k1 = numeric(), Eq32T_Hold_k2 = numeric(),
    Eq32T_Val_k1 = numeric(),  Eq32T_Val_k2 = numeric(),
    SignFlip_k1 = numeric(),          # fraction with any SKLE sign reversal
    SignFlip_k2 = numeric(), SignFlip_k3 = numeric(), SignFlip_k4 = numeric()
  )
  write.csv(headers, output_file, row.names = FALSE)
  cat(sprintf(">>> Initialized fresh output file: %s\n", output_file))
  start_idx <- 1
} else {
  existing_data <- read.csv(output_file)
  start_idx <- nrow(existing_data) + 1
  cat(sprintf(">>> Found existing file. Resuming from scenario %d...\n", start_idx))
}

# ==========================================
# 3. HELPER FUNCTIONS
# ==========================================

# Inverse link
get_link_inverse <- function(eta, link_name) {
  eta <- pmax(-20, pmin(20, eta))
  if (link_name == "logit")        return(plogis(eta))
  else if (link_name == "probit")  return(pnorm(eta))
  else if (link_name == "cloglog") return(1 - exp(-exp(eta)))
  else if (link_name == "neglog")  return(exp(-exp(-eta)))
}

# d(mu)/d(eta) for each link, required by the weight in Equation (7).
get_mu_eta <- function(eta, link_name) {
  eta <- pmax(-20, pmin(20, eta))
  if (link_name == "logit") {
    mu <- plogis(eta); return(mu * (1 - mu))
  } else if (link_name == "probit") {
    return(dnorm(eta))
  } else if (link_name == "cloglog") {
    return(exp(eta) * exp(-exp(eta)))
  } else if (link_name == "neglog") {
    return(exp(-eta) * exp(-exp(-eta)))
  }
}
# Exact simplex sampler via Barndorff-Nielsen & Jorgensen (1991), Eq. (19):
#   d(y; mu) / sigma^2 ~ chi-squared(1).
# Always terminates, unlike rejection sampling, which stalls when mu
# approaches 0 or 1 under the cloglog and neglog links.
rsimplex_exact <- function(n, mu, dispersion) {
  cc <- rchisq(n, df = 1)
  t  <- dispersion * cc * mu^2 * (1 - mu)^2
  disc <- sqrt(t * (4 * mu * (1 - mu) + t))
  y1 <- ((2 * mu + t) - disc) / (2 * (1 + t))
  y2 <- ((2 * mu + t) + disc) / (2 * (1 + t))
  h1 <- mu + y1 * (1 - 2 * mu)
  h2 <- mu + y2 * (1 - 2 * mu)
  ifelse(runif(n) < h1 / (h1 + h2), y2, y1)
}
# Deviance-only objective. The remaining terms of Equation (4) do not
# involve beta, so the argmin is unchanged.
neg_log_lik_standard <- function(beta, X, y, link_name) {
  eta <- X %*% beta
  mu  <- get_link_inverse(eta, link_name)
  mu  <- pmax(1e-8, pmin(1 - 1e-8, mu))
  d   <- (y - mu)^2 / (y * (1 - y) * mu^2 * (1 - mu)^2)
  return(sum(d))
}

# ==========================================
# 4. WORKER FUNCTION (Executes 1 Replication)
# ==========================================
run_simulation_rep <- function(dummy_idx, n, rho, sigma, link_name, p_val,
                               scenario_idx) {

  # The seed depends only on the dimension, the scenario and the replication
  # index, so the results are independent of the number of cores, of how
  # parSapply distributes the work, and of where a resumed run restarts.  The
  # scenario multiplier exceeds the replication count, so the streams cannot
  # collide across scenarios.
  set.seed(BASE_SEED + p_val * 10000000L +
           (scenario_idx - 1L) * 100000L + dummy_idx)

  # Newhouse and Oman (1971): sum of squared coefficients equals 1
  beta_true <- rep(1, p_val) / sqrt(p_val)

  # --- A. GENERATE DATA (McDonald & Galarneau, 1975) ---
  Z <- matrix(rnorm(n * (p_val + 1)), nrow = n, ncol = p_val + 1)
  X_raw <- matrix(0, nrow = n, ncol = p_val)
  for (j in 1:p_val) {
    X_raw[, j] <- sqrt(1 - rho^2) * Z[, j] + rho * Z[, p_val + 1]
  }
  X <- scale(X_raw)

  eta_true <- X %*% beta_true

  # Scale the design rather than beta, so that sd(eta) = 1 and
  # ||beta||^2 = 1 hold simultaneously.
  eta_sd <- sd(as.vector(eta_true))
  if (!is.finite(eta_sd) || eta_sd < 1e-8) return(NULL)
  X <- X / eta_sd
  eta_true <- eta_true / eta_sd
  
  mu_true <- get_link_inverse(eta_true, link_name)
  
  # Numerical guard; the rate at which it binds is recorded.
  clamp_hits <- mean(mu_true < MU_GUARD | mu_true > 1 - MU_GUARD)
  mu_gen <- pmax(MU_GUARD, pmin(1 - MU_GUARD, mu_true))

  y <- rsimplex_exact(n, mu_gen, sigma)
  if (is.null(y) || any(!is.finite(y))) return(NULL)
  y <- pmax(Y_GUARD, pmin(1 - Y_GUARD, y))

  # --- B. COMPUTE SMLE ---
  beta_init <- rep(0, p_val)
  opt_res <- tryCatch({
    optim(par = beta_init, fn = neg_log_lik_standard, X = X, y = y,
          link_name = link_name, method = "BFGS",
          control = list(maxit = 200, reltol = 1e-8))
  }, error = function(e) NULL)

  # Record the failure rather than discarding the replication.
  if (is.null(opt_res) || opt_res$convergence != 0) {
    return(c(rep(NA_real_, 5), clamp_hits, 0, rep(NA_real_, 23)))
  }
  b_mle <- opt_res$par

  # --- C. WEIGHTS & INFORMATION MATRIX ---
  eta_hat <- X %*% b_mle
  mu_hat  <- get_link_inverse(eta_hat, link_name)
  mu_hat  <- pmax(1e-8, pmin(1 - 1e-8, mu_hat))

  d_i        <- (y - mu_hat)^2 / (y * (1 - y) * mu_hat^2 * (1 - mu_hat)^2)
  sigma2_est <- sum(d_i) / (n - p_val)

  P_val_vec <- mu_hat * (1 - mu_hat)
  d_mu      <- get_mu_eta(eta_hat, link_name)
  # IRLS weight, Equation (7):
  #   w_i = [ 3*sigma^2/(mu(1-mu)) + 1/(mu(1-mu))^3 ] * (dmu/deta)^2
  w <- (3 * sigma2_est / P_val_vec + 1 / P_val_vec^3) * d_mu^2

  if (any(!is.finite(w)) || any(w <= 0)) {
    return(c(rep(NA_real_, 5), clamp_hits, 0, rep(NA_real_, 23)))
  }

  X_weighted <- X * sqrt(w)
  S <- crossprod(X_weighted)          # Information matrix, Equation (10)

  # --- D. SHRINKAGE PARAMETERS ---
  eig  <- eigen(S, symmetric = TRUE)
  vals <- eig$values
  vecs <- eig$vectors
  if (min(vals) <= 0 || any(!is.finite(vals))) {
    return(c(rep(NA_real_, 5), clamp_hits, 0, rep(NA_real_, 23)))
  }

  alpha_sq <- as.vector(crossprod(vecs, b_mle))^2
  alpha_sq <- pmax(alpha_sq, 1e-12)   # guard against exact zeros

  # Equation (37), geometric mean estimator (Kibria, 2003).
  # Computed in logs to avoid underflow of the product at large p.
  k1 <- sigma2_est / exp(mean(log(alpha_sq)))
  # Equation (38), square-root median estimator (Muniz and Kibria, 2009).
  k2 <- median(sqrt(sigma2_est / alpha_sq))

  k1 <- max(1e-10, k1)
  k2 <- max(1e-10, k2)

  # Regularity condition k < 2*lambda_min
  thm2_bound <- 2 * min(vals)
  hold_k1 <- as.numeric(k1 < thm2_bound)
  hold_k2 <- as.numeric(k2 < thm2_bound)

  # Equation (33) of the manuscript, the actual Theorem 2 condition:
  #   alpha' [ (sigma^2/(3k^2)) diag((2*lambda_i*k - k^2)/lambda_i) ]^-1 alpha < 1
  # which simplifies to
  #   sum_i  3*k*lambda_i*alpha_i^2 / (sigma^2 * (2*lambda_i - k))  <  1
  # The regularity condition keeps every (2*lambda_i - k) positive; without it the
  # quadratic form is not defined as a positive-definite condition, so
  # Eq (33) is recorded as failing.
  eq32 <- function(k_val) {
    if (k_val >= thm2_bound) return(c(NA_real_, 0))
    val <- sum(3 * k_val * vals * alpha_sq / (sigma2_est * (2 * vals - k_val)))
    c(val, as.numeric(val < 1))
  }
  e1 <- eq32(k1); e2 <- eq32(k2)
  eq32_val_k1 <- e1[1]; eq32_hold_k1 <- e1[2]
  eq32_val_k2 <- e2[1]; eq32_hold_k2 <- e2[2]

  # The same condition evaluated at the TRUE canonical coefficients.
  # E[sum(alpha_hat^2)] = sum(alpha^2) + sigma^2 * sum(1/lambda), so alpha_hat
  # is inflated by roughly the SMLE MSE. Evaluating at the true alpha shows
  # whether the condition genuinely fails or is merely obscured by that noise.
  alpha_sq_true <- as.vector(crossprod(vecs, beta_true))^2
  eq32T <- function(k_val) {
    if (k_val >= thm2_bound) return(c(NA_real_, 0))
    val <- sum(3 * k_val * vals * alpha_sq_true / (sigma2_est * (2 * vals - k_val)))
    c(val, as.numeric(val < 1))
  }
  t1 <- eq32T(k1); t2 <- eq32T(k2)
  eq32T_val_k1 <- t1[1]; eq32T_hold_k1 <- t1[2]
  eq32T_val_k2 <- t2[1]; eq32T_hold_k2 <- t2[2]

  # k3 = the root of Equation (33): the largest k for which Theorem 2 holds.
  # The LHS is 0 at k = 0 and increasing, so a unique root exists below the
  # regularity condition bound whenever the LHS exceeds 1 there.
  f32 <- function(k_val) {
    sum(3 * k_val * vals * alpha_sq / (sigma2_est * (2 * vals - k_val))) - 1
  }
  hi <- thm2_bound * (1 - 1e-8)
  k3 <- tryCatch({
    if (f32(1e-12) * f32(hi) < 0) uniroot(f32, c(1e-12, hi), tol = 1e-12)$root else hi
  }, error = function(e) hi)
  k3 <- max(1e-10, k3)

  # k4: square-root median estimator truncated at the regularity condition bound.
  # Equals k2 when the regularity condition already holds, 2*lambda_min otherwise.
  k4 <- max(1e-10, min(k2, thm2_bound * (1 - 1e-8)))
  k4_binds <- as.numeric(k2 >= thm2_bound)

  # --- E. ESTIMATORS ---
  I_p <- diag(p_val)
  compute_estimators <- function(k_val) {
    inv_mat <- tryCatch(solve(S + k_val * I_p),
                        error = function(e) MASS::ginv(S + k_val * I_p))
    # SRRE, Equation (12)
    b_srre <- inv_mat %*% S %*% b_mle
    # SKLE, Equation (19): (S + kI)^-1 (S - kI) b = b - 2k (S + kI)^-1 b
    b_skle <- b_mle - 2 * k_val * (inv_mat %*% b_mle)
    list(srre = as.vector(b_srre), skle = as.vector(b_skle))
  }

  ests_k1 <- compute_estimators(k1)
  ests_k2 <- compute_estimators(k2)
  ests_k3 <- compute_estimators(k3)
  ests_k4 <- compute_estimators(k4)

  # Uncapped, k > lambda_min makes (S - kI) indefinite and the SKLE can
  # reverse coefficient signs. This is a real property of the KL family and
  # is worth reporting rather than preventing.
  flip_k1 <- mean(ests_k1$skle < 0)
  flip_k2 <- mean(ests_k2$skle < 0)
  flip_k3 <- mean(ests_k3$skle < 0)
  flip_k4 <- mean(ests_k4$skle < 0)

  # --- F. SQUARED ERRORS ---
  se_mle      <- sum((b_mle          - beta_true)^2)
  se_srre_k1  <- sum((ests_k1$srre   - beta_true)^2)
  se_skle_k1  <- sum((ests_k1$skle   - beta_true)^2)
  se_srre_k2  <- sum((ests_k2$srre   - beta_true)^2)
  se_skle_k2  <- sum((ests_k2$skle   - beta_true)^2)
  se_srre_k3  <- sum((ests_k3$srre   - beta_true)^2)
  se_skle_k3  <- sum((ests_k3$skle   - beta_true)^2)
  se_srre_k4  <- sum((ests_k4$srre   - beta_true)^2)
  se_skle_k4  <- sum((ests_k4$skle   - beta_true)^2)

  return(c(se_mle, se_srre_k1, se_skle_k1, se_srre_k2, se_skle_k2,
           clamp_hits, 1, k1, k2, hold_k1, hold_k2, flip_k1, flip_k2,
           eq32_hold_k1, eq32_hold_k2, eq32_val_k1, eq32_val_k2,
           se_srre_k3, se_skle_k3, k3, flip_k3,
           eq32T_hold_k1, eq32T_hold_k2, eq32T_val_k1, eq32T_val_k2,
           se_srre_k4, se_skle_k4, k4, flip_k4, k4_binds))
}

# ==========================================
# 5. MAIN PARALLEL EXECUTION
# ==========================================
grid <- expand.grid(n = N_VALUES, rho = RHO_VALUES,
                    sigma = SIGMA_VALUES, link = LINKS,
                    stringsAsFactors = FALSE)

n_cores <- max(1, detectCores() - 1)
cl <- makeCluster(n_cores)
clusterEvalQ(cl, { library(MASS) })
clusterExport(cl, c("run_simulation_rep", "get_link_inverse", "get_mu_eta",
                    "neg_log_lik_standard", "rsimplex_exact",
                    "MU_GUARD", "Y_GUARD", "BASE_SEED"))

cat(sprintf(">>> Running %d scenarios on %d cores.\n", nrow(grid), n_cores))

for (i in start_idx:nrow(grid)) {
  params <- grid[i, ]

  res_mat <- parSapply(cl, 1:REPLICATIONS, function(x, n, rho, sigma, link, p, scen) {
    out <- run_simulation_rep(x, n, rho, sigma, link, p, scen)
    if (is.null(out)) rep(NA_real_, 30) else out
  }, n = params$n, rho = params$rho, sigma = params$sigma,
     link = as.character(params$link), p = P_VAL, scen = i)

  res_mat <- matrix(as.numeric(res_mat), nrow = 30)

  conv_flag    <- res_mat[7, ]
  conv_flag[is.na(conv_flag)] <- 0
  n_attempted  <- sum(!is.na(res_mat[6, ]))
  n_successful <- sum(conv_flag == 1)

  if (n_successful < 2) {
    cat(sprintf("Scenario %d (%s, sigma=%.1f, rho=%.2f, n=%d): FAILED\n",
                i, params$link, params$sigma, params$rho, params$n))
    next
  }

  ok <- which(conv_flag == 1)
  idx   <- c(1:5, 18, 19, 26, 27)  ##### CHANGES 11 & 13 #####
  means <- rowMeans(res_mat[idx, ok, drop = FALSE], na.rm = TRUE)
  meds  <- apply(res_mat[idx, ok, drop = FALSE], 1, median, na.rm = TRUE)

  row_out <- data.frame(
    Link = as.character(params$link),
    Sigma = params$sigma,
    Rho = params$rho,
    n = params$n,
    Successful_Reps = n_successful,
    MLE     = means[1],
    SRRE_k1 = means[2], SKLE_k1 = means[3],
    SRRE_k2 = means[4], SKLE_k2 = means[5],
    SRRE_k3 = means[6], SKLE_k3 = means[7],
    SRRE_k4 = means[8], SKLE_k4 = means[9],
    MLE_Med     = meds[1],
    SRRE_k1_Med = meds[2], SKLE_k1_Med = meds[3],
    SRRE_k2_Med = meds[4], SKLE_k2_Med = meds[5],
    SRRE_k3_Med = meds[6], SKLE_k3_Med = meds[7],
    SRRE_k4_Med = meds[8], SKLE_k4_Med = meds[9],
    Clamp_Rate   = mean(res_mat[6, ], na.rm = TRUE),
    Conv_Rate    = n_successful / max(1, n_attempted),
    Mean_k1      = mean(res_mat[8,  ok], na.rm = TRUE),
    Mean_k2      = mean(res_mat[9,  ok], na.rm = TRUE),
    Mean_k3      = mean(res_mat[20, ok], na.rm = TRUE),
    Mean_k4      = mean(res_mat[28, ok], na.rm = TRUE),
    K4_Binds     = mean(res_mat[30, ok], na.rm = TRUE),
    Prov_Hold_k1 = mean(res_mat[10, ok], na.rm = TRUE),
    Prov_Hold_k2 = mean(res_mat[11, ok], na.rm = TRUE),
    Eq32_Hold_k1 = mean(res_mat[14, ok], na.rm = TRUE),
    Eq32_Hold_k2 = mean(res_mat[15, ok], na.rm = TRUE),
    Eq32_Val_k1  = median(res_mat[16, ok], na.rm = TRUE),
    Eq32_Val_k2  = median(res_mat[17, ok], na.rm = TRUE),
    Eq32T_Hold_k1 = mean(res_mat[22, ok], na.rm = TRUE),
    Eq32T_Hold_k2 = mean(res_mat[23, ok], na.rm = TRUE),
    Eq32T_Val_k1  = median(res_mat[24, ok], na.rm = TRUE),
    Eq32T_Val_k2  = median(res_mat[25, ok], na.rm = TRUE),
    SignFlip_k1  = mean(res_mat[12, ok], na.rm = TRUE),
    SignFlip_k2  = mean(res_mat[13, ok], na.rm = TRUE),
    SignFlip_k3  = mean(res_mat[21, ok], na.rm = TRUE),
    SignFlip_k4  = mean(res_mat[29, ok], na.rm = TRUE)
  )

  write.table(row_out, output_file, sep = ",", append = TRUE,
              row.names = FALSE, col.names = FALSE)

  cat(sprintf(paste0("Scenario %3d | %-7s sigma=%.1f rho=%.2f n=%3d | ",
                     "reps=%4d conv=%.2f clamp=%.3f | MLE=%.4f SKLE1=%.4f | ",
                     "prov=%.2f eq32T=%.2f | SKLE3=%.4f SKLE4=%.4f\n"),
              i, params$link, params$sigma, params$rho, params$n,
              n_successful, row_out$Conv_Rate, row_out$Clamp_Rate,
              means[1], means[3], row_out$Prov_Hold_k2,
              row_out$Eq32T_Hold_k2, means[7], means[9]))
  flush.console()
}

stopCluster(cl)
cat(sprintf("\n>>> Complete. Results written to:\n    %s\n", output_file))
