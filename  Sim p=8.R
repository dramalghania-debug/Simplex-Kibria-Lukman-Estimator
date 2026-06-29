library(VGAM)
library(MASS)
library(parallel)

# ==========================================
# 1. CONFIGURATION
# ==========================================
P_VAL <- 8 # Updated for the p=8 scenario
N_VALUES <- c(25, 50, 100, 200)
RHO_VALUES <- c(0.80, 0.90, 0.95, 0.99)  
SIGMA_VALUES <- c(0.5, 1.0, 1.5, 2.0)
LINKS <- c("logit", "probit", "cloglog", "neglog")
REPLICATIONS <- 1000

# ==========================================
# 2. FOLDER & FILE SETUP (SMART RESUME)
# ==========================================
# Robust Desktop Path Detection
if (.Platform$OS.type == "windows") {
  desktop_path <- file.path(Sys.getenv("USERPROFILE"), "Desktop")
} else {
  desktop_path <- file.path(Sys.getenv("HOME"), "Desktop")
}

# Static folder name with "New" added to prevent overwriting old runs
folder_name <- paste0("New_SKLE_Simulation_p", P_VAL)
output_dir <- file.path(desktop_path, folder_name)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(sprintf(">>> Created new folder on Desktop: %s\n", folder_name))
}

output_file <- file.path(output_dir, paste0("Simulation_Results_p", P_VAL, ".csv"))

# Check if file exists to determine if we are starting fresh or resuming
if (!file.exists(output_file)) {
  # Added the Successful_Reps column to the headers
  headers <- data.frame(Link=character(), Sigma=numeric(), Rho=numeric(), n=numeric(), 
                        Successful_Reps=numeric(), 
                        MLE=numeric(), 
                        SRE_k1=numeric(), SKLE_k1=numeric(), 
                        SRE_k2=numeric(), SKLE_k2=numeric())
  write.csv(headers, output_file, row.names = FALSE)
  cat(sprintf(">>> Initialized fresh output file: %s\n", output_file))
  start_idx <- 1
} else {
  # If file exists, count the rows to see where we left off
  existing_data <- read.csv(output_file)
  start_idx <- nrow(existing_data) + 1
  cat(sprintf(">>> Found existing file! Resuming simulation from scenario %d...\n", start_idx))
}

# Set seed for perfect reproducibility
set.seed(2025)

# ==========================================
# 3. HELPER FUNCTIONS
# ==========================================
# Link Functions Inverse
get_link_inverse <- function(eta, link_name) {
  eta <- pmax(-20, pmin(20, eta))
  if (link_name == "logit") return(plogis(eta))
  else if (link_name == "probit") return(pnorm(eta))
  else if (link_name == "cloglog") return(1 - exp(-exp(eta)))
  else if (link_name == "neglog") return(exp(-exp(-eta)))
}

# Standard Simplex Log-Likelihood
neg_log_lik_standard <- function(beta, X, y, link_name) {
  eta <- X %*% beta
  mu <- get_link_inverse(eta, link_name)
  mu <- pmax(1e-8, pmin(1 - 1e-8, mu))
  d <- (y - mu)^2 / (y * (1 - y) * mu^2 * (1 - mu)^2)
  return(sum(d))
}

# ==========================================
# 4. WORKER FUNCTION (Executes 1 Replication)
# ==========================================
run_simulation_rep <- function(dummy_idx, n, rho, sigma, link_name, p_val) {
  
  # Newhouse and Oman (1971) Condition: Sum of squared coefficients equals 1
  beta_true <- rep(1, p_val) / sqrt(p_val)
  
  # --- A. GENERATE DATA (McDonald & Galarneau, 1975) ---
  Z <- matrix(rnorm(n * (p_val + 1)), nrow = n, ncol = p_val + 1)
  X_raw <- matrix(0, nrow = n, ncol = p_val)
  
  # Apply the exact mathematical generation equation
  for(j in 1:p_val) {
    X_raw[, j] <- sqrt(1 - rho^2) * Z[, j] + rho * Z[, p_val + 1]
  }
  
  # Standardize Predictors (Required for Ridge/SKLE)
  X <- scale(X_raw)
  
  eta_true <- X %*% beta_true
  mu_true <- get_link_inverse(eta_true, link_name)
  mu_gen <- pmax(0.01, pmin(0.99, mu_true))
  
  y <- tryCatch({ VGAM::rsimplex(n, mu = mu_gen, dispersion = sigma) }, error = function(e) NULL)
  if (is.null(y)) return(NULL)
  y <- pmax(1e-6, pmin(1 - 1e-6, y))
  
  # --- B. COMPUTE STANDARD MLE ---
  beta_init <- rep(0, p_val)
  opt_res <- tryCatch({
    optim(par = beta_init, fn = neg_log_lik_standard, X = X, y = y, link_name = link_name,
          method = "BFGS", control = list(maxit = 200, reltol = 1e-8))
  }, error = function(e) NULL)
  
  if (is.null(opt_res) || opt_res$convergence != 0) return(NULL)
  b_mle <- opt_res$par
  
  # --- C. COMPUTE WEIGHTS & MATRICES ---
  eta_hat <- X %*% b_mle
  mu_hat <- get_link_inverse(eta_hat, link_name)
  mu_hat <- pmax(1e-8, pmin(1 - 1e-8, mu_hat))
  
  d_i <- (y - mu_hat)^2 / (y * (1 - y) * mu_hat^2 * (1 - mu_hat)^2)
  sigma2_est <- sum(d_i) / (n - p_val)
  
  P_val_vec <- mu_hat * (1 - mu_hat)
  
  # Using estimated dispersion (sigma2_est) instead of true population sigma
  w <- (3 * sigma2_est * P_val_vec) + (1 / P_val_vec)
  X_weighted <- X * sqrt(w)
  S <- t(X_weighted) %*% X_weighted # Information Matrix
  
  # --- D. CALCULATE k PARAMETERS & THEORETICAL BOUNDS ---
  eig <- eigen(S)
  vals <- eig$values
  vecs <- eig$vectors
  alpha_sq <- as.vector(crossprod(vecs, b_mle))^2
  
  # Theorem 2 Hard Ceiling (k < 2 * lambda_min)
  upper_limit <- (2 * min(vals)) - 1e-5
  
  # k1 (Geometric/HKB-type)
  k1 <- max(1e-10, min(sigma2_est / (prod(alpha_sq)^(1/p_val)), upper_limit))
  # k2 (Median-type)
  k2 <- max(1e-10, min(median(sqrt(sigma2_est / alpha_sq)), upper_limit))
  
  I_p <- diag(p_val)
  
  # --- E. COMPUTE ESTIMATORS (SRE & SKLE) ---
  compute_estimators <- function(k_val) {
    # Inverse of (S + kI)
    inv_mat <- tryCatch(solve(S + k_val * I_p), error = function(e) MASS::ginv(S + k_val * I_p))
    
    # SRE (Ridge) = (S + kI)^-1 * S * MLE
    b_sre <- inv_mat %*% S %*% b_mle
    
    # SKLE (Kibria-Lukman) = MLE - 2k(S + kI)^-1 * MLE
    b_skle <- b_mle - 2 * k_val * (inv_mat %*% b_mle)
    
    return(list(sre = b_sre, skle = b_skle))
  }
  
  ests_k1 <- compute_estimators(k1)
  ests_k2 <- compute_estimators(k2)
  
  # --- F. CALCULATE SQUARED ERRORS ---
  se_mle     <- sum((b_mle - beta_true)^2)
  se_sre_k1  <- sum((ests_k1$sre - beta_true)^2)
  se_skle_k1 <- sum((ests_k1$skle - beta_true)^2)
  se_sre_k2  <- sum((ests_k2$sre - beta_true)^2)
  se_skle_k2 <- sum((ests_k2$skle - beta_true)^2)
  
  return(c(se_mle, se_sre_k1, se_skle_k1, se_sre_k2, se_skle_k2))
}

# ==========================================
# 5. MAIN PARALLEL EXECUTION
# ==========================================
grid <- expand.grid(n = N_VALUES, rho = RHO_VALUES, sigma = SIGMA_VALUES, link = LINKS)
num_cores <- max(1, detectCores() - 1)

# Mac-optimized FORK cluster
cl <- makeCluster(num_cores, type = "FORK")

clusterEvalQ(cl, { library(VGAM); library(MASS) })
clusterExport(cl, c("get_link_inverse", "neg_log_lik_standard", "run_simulation_rep"))

cat(sprintf("\n>>> Starting Standard Algorithm Simulation (P=%d) [%d Total Scenarios]...\n", P_VAL, nrow(grid)))

if (start_idx > nrow(grid)) {
  cat("\n>>> Simulation is already 100% complete! No further action needed.\n")
} else {
  for(i in start_idx:nrow(grid)) {
    params <- grid[i, ]
    
    # Run Parallel Replications
    res_list <- parSapply(cl, 1:REPLICATIONS, function(x, n, rho, sigma, link, p) {
      run_simulation_rep(x, n, rho, sigma, link, p)
    }, n = params$n, rho = params$rho, sigma = params$sigma, link = as.character(params$link), p = P_VAL)
    
    # Handle failures gracefully
    if (is.list(res_list)) valid_res <- do.call(cbind, res_list[!sapply(res_list, is.null)])
    else valid_res <- res_list
    
    if (!is.null(valid_res) && ncol(valid_res) > 0) {
      means <- rowMeans(valid_res)
      actual_reps <- ncol(valid_res) # Count exactly how many converged
      
      # Create Data Frame for this row (Calculating SMSE)
      row_data <- data.frame(
        Link = params$link,
        Sigma = params$sigma,
        Rho = params$rho,
        n = params$n,
        Successful_Reps = actual_reps, # Add to CSV
        MLE = means[1],
        SRE_k1 = means[2],
        SKLE_k1 = means[3],
        SRE_k2 = means[4],
        SKLE_k2 = means[5]
      )
      
      # APPEND to CSV immediately
      write.table(row_data, output_file, sep = ",", row.names = FALSE, col.names = FALSE, append = TRUE)
      
      # Print progress to console including the successful reps
      cat(sprintf("Saved Scenario %-3d | %-7s | Sig=%.1f | Rho=%.2f | n=%-3d | Reps: %-4d -> MLE: %.4f | SKLE(k1): %.4f\n", 
                  i, params$link, params$sigma, params$rho, params$n, actual_reps, means[1], means[3]))
    }
  }
}

stopCluster(cl)
cat(sprintf("\n>>> DONE! All results saved to:\n%s\n", output_file))