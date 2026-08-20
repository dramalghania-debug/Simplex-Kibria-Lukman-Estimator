# =============================================================================
# body_fat.R
#
# Body Fat application for the Simplex Kibria-Lukman Estimator (SKLE).
# Fits the SMLE, SRRE and SKLE under four link functions and reports the
# estimated MSE of each, together with the shrinkage parameter diagnostics.
#
# The data are loaded from the mfp package; no data file is required.
#
# Equation numbers refer to the manuscript.
# =============================================================================

# Install required packages if missing
if(!require(mfp)) install.packages("mfp")

library(mfp)

# ==========================================
# 1. DATA PREPARATION & TRUNCATION
# ==========================================
data(bodyfat)

# Observation 182 records siri = 0, a physically impossible body fat
# percentage and a known recording error in the Penrose data. With the
# weight in Equation (7), which contains 1/P^3, a single boundary observation
# makes the weight matrix numerically singular and IRLS diverges.
bad <- which(bodyfat$siri <= 0 | bodyfat$siri >= 100)
if (length(bad)) {
  cat("Excluding", length(bad), "observation(s) with impossible body fat:", bad, "\n")
  bodyfat <- bodyfat[-bad, ]
}

n <- nrow(bodyfat)

# Convert Siri body fat percentage to a raw proportion
y_raw <- bodyfat$siri / 100

# Apply the strict 1e-5 truncation method 
y <- pmax(1e-5, pmin(1 - 1e-5, y_raw))

cat("=========================================================================\n")
cat("   DATASET: Body Fat Model (n =", n, ", p = 11)\n")
cat("   NOTE: 'density' excluded to prevent deterministic collinearity\n")
cat("=========================================================================\n\n")

# ==========================================
# 2. PREDICTOR SETUP & SCALING
# ==========================================
# 11 strictly physical/biological measurements
keep_vars <- c("age", "weight", "height", "neck", "chest", 
               "abdomen", "hip", "thigh", "knee", "forearm", "wrist")

X_raw <- as.matrix(bodyfat[, keep_vars])
p <- ncol(X_raw)

# Standardizing the matrix is mathematically required to separate 
# structural collinearity from simple unit-scale disparities.
X_std <- scale(X_raw)
X_aug <- cbind(1, X_std)
p_aug <- ncol(X_aug)

X_means <- colMeans(X_raw)
X_sds <- apply(X_raw, 2, sd)

# ==========================================
# 3. LINK FUNCTION DEFINITIONS
# ==========================================
get_link_funcs <- function(link_name) {
  if(link_name == "logit") {
    linkfun <- function(mu) qlogis(mu)
    linkinv <- function(eta) { eta <- pmax(-10, pmin(10, eta)); 1/(1+exp(-eta)) }
    mu.eta  <- function(eta) { mu <- linkinv(eta); mu*(1-mu) }
  } else if(link_name == "probit") {
    linkfun <- function(mu) qnorm(mu)
    linkinv <- function(eta) { eta <- pmax(-10, pmin(10, eta)); pnorm(eta) }
    mu.eta  <- function(eta) { eta <- pmax(-10, pmin(10, eta)); dnorm(eta) }
  } else if(link_name == "cloglog") {
    linkfun <- function(mu) log(-log(1-mu))
    linkinv <- function(eta) { eta <- pmax(-10, pmin(10, eta)); 1 - exp(-exp(eta)) }
    mu.eta  <- function(eta) { eta <- pmax(-10, pmin(10, eta)); exp(eta) * exp(-exp(eta)) }
  } else if(link_name == "neglog") { 
    linkfun <- function(mu) -log(-log(mu))
    linkinv <- function(eta) { eta <- pmax(-10, pmin(10, eta)); exp(-exp(-eta)) }
    mu.eta  <- function(eta) { eta <- pmax(-10, pmin(10, eta)); exp(-exp(-eta)) * exp(-eta) }
  }
  return(list(linkfun=linkfun, linkinv=linkinv, mu.eta=mu.eta))
}

# ==========================================
# 4. HELPER: BACK-TRANSFORM COEFFICIENTS
# ==========================================
back_transform <- function(beta_z) {
  beta_raw_preds <- beta_z[2:p_aug] / X_sds
  beta_raw_0 <- beta_z[1] - sum(beta_raw_preds * X_means)
  return(c(beta_raw_0, beta_raw_preds))
}

# ==========================================
# 5. MAIN IRLS LOOP & ESTIMATION
# ==========================================
links_to_test <- c("logit", "probit", "cloglog", "neglog")
all_tables <- list()
cn_values <- list()
ci_values <- list() # Added list to store condition indices
conv_values <- list()    # convergence status
k_diag <- data.frame()    # shrinkage diagnostics
beta_std_store <- list()    # standardised coefs for the ML check

mse_summary <- data.frame(
  Link = character(), Condition_Num = numeric(), SMLE = numeric(), 
  SRE_k1 = numeric(), SKL_k1 = numeric(), SRE_k2 = numeric(),
  SKL_k2 = numeric(), SRE_k3 = numeric(), SKL_k3 = numeric(),
  SRE_k4 = numeric(), SKL_k4 = numeric(), stringsAsFactors = FALSE
)

for(lnk in links_to_test) {
  
  funcs <- get_link_funcs(lnk)
  y_link_val <- funcs$linkfun(y)
  mod_lm <- lm(y_link_val ~ X_std)
  b_curr <- coef(mod_lm)
  sigma_curr <- sum(residuals(mod_lm)^2)/(n-p_aug)
  
  converged <- FALSE
  iter_used <- NA

  # Robust IRLS Algorithm for Simplex GLM
  for(iter in 1:100) {
    eta_it <- X_aug %*% b_curr
    mu_it  <- funcs$linkinv(eta_it)
    mu_it  <- pmax(1e-7, pmin(1 - 1e-7, mu_it)) 
    
    P_val <- mu_it * (1 - mu_it)         
    d_mu  <- funcs$mu.eta(eta_it)        
    
    y_var <- pmax(1e-9, y * (1 - y))
    u_dev <- (y - mu_it)^2 / (y_var * P_val^2)
    sigma_curr <- sum(u_dev) / (n - p_aug)
    # IRLS weight, Equation (7):
    #     w = [ 3*sigma^2/(mu(1-mu)) + 1/(mu(1-mu))^3 ] * (dmu/deta)^2
    w_ii <- (3 * sigma_curr / P_val + 1 / P_val^3) * d_mu^2
    W <- diag(as.vector(w_ii))
    
    u_adj <- ((y - mu_it) / P_val) * (u_dev + (1 / P_val^2))
    # Working response, Equation (8). Cov(beta) = sigma^2 (X'WX)^-1 makes
    # the Fisher information (X'WX)/sigma^2, so sigma^2 cancels here.
    y_star <- eta_it + (u_adj * d_mu) / w_ii
    
    XtWX <- t(X_aug) %*% W %*% X_aug
    XtWy_star <- t(X_aug) %*% W %*% y_star
    b_new <- solve(XtWX, XtWy_star)
    if(max(abs(b_new - b_curr)) < 1e-8) {
      b_curr    <- b_new
      converged <- TRUE
      iter_used <- iter
      break
    }
    b_curr <- b_new
  }

  if(!converged) {
    warning(sprintf("IRLS did NOT converge for the %s link. Results below are not valid.", lnk))
    iter_used <- 100
  }
  conv_values[[lnk]] <- list(converged = converged, iter = iter_used)

  
  beta_mle_std <- as.vector(b_curr)
  sigma_est <- sigma_curr 
  
  eta_final <- X_aug %*% beta_mle_std
  mu_final <- pmax(1e-7, pmin(1-1e-7, funcs$linkinv(eta_final)))
  
  P_final    <- mu_final * (1 - mu_final)
  d_mu_final <- funcs$mu.eta(eta_final)
  w_final <- (3 * sigma_est / P_final + 1 / P_final^3) * d_mu_final^2
  W_final <- diag(as.vector(w_final))
  
  # Fisher Information Matrix (S) & Diagnostics
  S <- t(X_aug) %*% W_final %*% X_aug
  eig <- eigen(S)
  vals <- pmax(1e-12, eig$values) # Safety net for minimum eigenvalue
  vecs <- eig$vectors
  
  cond_num <- sqrt(max(vals) / min(vals))
  cn_values[[lnk]] <- cond_num
  
  # Calculate Condition Indices (sqrt of max eigenvalue / each eigenvalue)
  cond_indices <- sqrt(max(vals) / vals)
  ci_values[[lnk]] <- cond_indices
  
  alpha <- t(vecs) %*% beta_mle_std
  alpha_sq <- as.vector(alpha^2)
  
  # ==========================================
  # 6. THEORETICAL BOUNDS & ESTIMATORS
  # ==========================================
  # Theorem 1 Bound
  k_thm1 <- sigma_est / sum(alpha_sq)
  
  # Theorem 2 Bound (Upper Limit)
  f_thm2 <- function(k) {
    sum((3 * k * vals * alpha_sq) / (sigma_est * (2 * vals - k))) - 1
  }
  
  upper_search <- 2 * min(vals) - 1e-10
  
  if(f_thm2(1e-10) * f_thm2(upper_search) < 0) {
    k_thm2 <- uniroot(f_thm2, interval = c(1e-10, upper_search))$root
  } else {
    k_thm2 <- upper_search 
  }
  
  k_max_theoretical <- min(k_thm1, k_thm2)
  
  k1_raw <- sigma_est / (prod(alpha_sq)^(1/p_aug))
  k2_raw <- median(sqrt(sigma_est / alpha_sq))
  
  # k1 and k2 are Equations (35) and (36). k_thm1 and k_thm2 are reported
  # as diagnostics, not imposed as constraints.
  k1 <- k1_raw
  k2 <- k2_raw

  # k3: the root of Equation (33), i.e. the largest k for which Theorem 2's
  #     condition is satisfied. Already computed above as k_thm2.
  # k4: the square-root median estimator truncated at the regularity
  #     condition k < 2*lambda_min. This is the rule recommended in the paper.
  lam_min2 <- 2 * min(vals)
  k3 <- max(1e-10, k_thm2)
  k4 <- max(1e-10, min(k2, lam_min2 * (1 - 1e-8)))
  k4_binds <- (k2 >= lam_min2)

  # record what k actually ended up being used
  # Reg_holds tests the REGULARITY CONDITION k < 2*lambda_min, which is what
  # Theorem 2 requires for D_Var to be positive definite. It is NOT Equation
  # (33) itself. Eq32_holds tests Equation (33).
  eq32_lhs <- function(kk) {
    if (kk >= lam_min2) return(NA_real_)
    sum(3 * kk * vals * alpha_sq / (sigma_est * (2 * vals - kk)))
  }
  k_diag <- rbind(k_diag, data.frame(
    Link = toupper(lnk), sigma2 = sigma_est,
    lam_min2 = lam_min2, k_thm1 = k_thm1, k_thm2 = k_thm2,
    k1 = k1, k2 = k2, k3 = k3, k4 = k4,
    Reg_holds_k1 = k1 < lam_min2,
    Reg_holds_k2 = k2 < lam_min2,
    ratio_k1 = k1 / lam_min2,
    ratio_k2 = k2 / lam_min2,
    Eq32_k1 = eq32_lhs(k1), Eq32_k2 = eq32_lhs(k2),
    k4_binds = k4_binds))
  
  b_mle_raw <- back_transform(beta_mle_std)
  mse_mle <- sigma_est * sum(1/vals)
  I_p <- diag(p_aug)
  
  # SRE Computations
  Inv_Sk1 <- solve(S + k1 * I_p)
  b_sre_k1_std <- as.vector(Inv_Sk1 %*% S %*% beta_mle_std)
  b_sre_k1_raw <- back_transform(b_sre_k1_std)
  mse_sre_k1 <- sigma_est * sum(vals/(vals+k1)^2) + k1^2 * sum(alpha_sq/(vals+k1)^2)
  
  # SKL Computations
  b_skl_k1_std <- as.vector(Inv_Sk1 %*% (S - k1 * I_p) %*% beta_mle_std)
  b_skl_k1_raw <- back_transform(b_skl_k1_std)
  mse_skl_k1 <- sigma_est * sum(((vals-k1)^2)/(vals*(vals+k1)^2)) + 4*k1^2 * sum(alpha_sq/(vals+k1)^2)
  
  Inv_Sk2 <- solve(S + k2 * I_p)
  b_sre_k2_std <- as.vector(Inv_Sk2 %*% S %*% beta_mle_std)
  b_sre_k2_raw <- back_transform(b_sre_k2_std)
  mse_sre_k2 <- sigma_est * sum(vals/(vals+k2)^2) + k2^2 * sum(alpha_sq/(vals+k2)^2)
  
  b_skl_k2_std <- as.vector(Inv_Sk2 %*% (S - k2 * I_p) %*% beta_mle_std)
  b_skl_k2_raw <- back_transform(b_skl_k2_std)
  mse_skl_k2 <- sigma_est * sum(((vals-k2)^2)/(vals*(vals+k2)^2)) + 4*k2^2 * sum(alpha_sq/(vals+k2)^2)

  # k3 and k4
  Inv_Sk3 <- solve(S + k3 * I_p)
  b_sre_k3_raw <- back_transform(as.vector(Inv_Sk3 %*% S %*% beta_mle_std))
  mse_sre_k3 <- sigma_est * sum(vals/(vals+k3)^2) + k3^2 * sum(alpha_sq/(vals+k3)^2)
  b_skl_k3_raw <- back_transform(as.vector(Inv_Sk3 %*% (S - k3 * I_p) %*% beta_mle_std))
  mse_skl_k3 <- sigma_est * sum(((vals-k3)^2)/(vals*(vals+k3)^2)) + 4*k3^2 * sum(alpha_sq/(vals+k3)^2)

  Inv_Sk4 <- solve(S + k4 * I_p)
  b_sre_k4_raw <- back_transform(as.vector(Inv_Sk4 %*% S %*% beta_mle_std))
  mse_sre_k4 <- sigma_est * sum(vals/(vals+k4)^2) + k4^2 * sum(alpha_sq/(vals+k4)^2)
  b_skl_k4_raw <- back_transform(as.vector(Inv_Sk4 %*% (S - k4 * I_p) %*% beta_mle_std))
  mse_skl_k4 <- sigma_est * sum(((vals-k4)^2)/(vals*(vals+k4)^2)) + 4*k4^2 * sum(alpha_sq/(vals+k4)^2)
  
  res_table <- data.frame(
    Variable = c("Intercept", keep_vars, "MSE"),
    SMLE = c(b_mle_raw, mse_mle),
    SRE_k1 = c(b_sre_k1_raw, mse_sre_k1),
    SKL_k1 = c(b_skl_k1_raw, mse_skl_k1),
    SRE_k2 = c(b_sre_k2_raw, mse_sre_k2),
    SKL_k2 = c(b_skl_k2_raw, mse_skl_k2),
    SRE_k3 = c(b_sre_k3_raw, mse_sre_k3),
    SKL_k3 = c(b_skl_k3_raw, mse_skl_k3),
    SRE_k4 = c(b_sre_k4_raw, mse_sre_k4),
    SKL_k4 = c(b_skl_k4_raw, mse_skl_k4)
  )
  all_tables[[lnk]] <- res_table
  beta_std_store[[lnk]] <- beta_mle_std
  
  mse_summary <- rbind(mse_summary, data.frame(
    Link = toupper(lnk), Condition_Num = cond_num, SMLE = mse_mle,
    SRE_k1 = mse_sre_k1, SKL_k1 = mse_skl_k1,
    SRE_k2 = mse_sre_k2, SKL_k2 = mse_skl_k2,
    SRE_k3 = mse_sre_k3, SKL_k3 = mse_skl_k3,
    SRE_k4 = mse_sre_k4, SKL_k4 = mse_skl_k4))
}

# ==========================================
# 7. OUTPUT RESULTS
# ==========================================
for(lnk in links_to_test) {
  cat(paste0("\n=========================================================================\n"))
  cat(paste0("   LINK FUNCTION: ", toupper(lnk), "\n"))
  cat(paste0("   CONVERGED: ", conv_values[[lnk]]$converged,
             "   (iterations: ", conv_values[[lnk]]$iter, ")\n"))
  cat(paste0("   WEIGHTED CONDITION NUMBER: ", round(cn_values[[lnk]], 4), "\n"))
  cat(paste0("   CONDITION INDICES:\n   ", paste(round(ci_values[[lnk]], 4), collapse = ", "), "\n"))
  cat(paste0("=========================================================================\n"))
  print(all_tables[[lnk]], row.names=FALSE, digits=4)
}

cat(paste0("\n=========================================================================\n"))
cat(paste0("   SUMMARY OF WEIGHTED DIAGNOSTICS & MEAN SQUARED ERRORS\n"))
cat(paste0("=========================================================================\n"))
print(mse_summary, row.names=FALSE, digits=4)

cat(paste0("\n=========================================================================\n"))
cat(paste0("   SHRINKAGE PARAMETER DIAGNOSTICS\n"))
cat(paste0("   Reg_holds = TRUE means k < 2*lambda_min (regularity condition).\n   Eq32 < 1 means Equation (33) itself is satisfied.\n"))
cat(paste0("=========================================================================\n"))
print(k_diag, row.names = FALSE, digits = 4)

# Independent check: refit each model by direct numerical maximum
# likelihood and compare. If the IRLS is right, these must agree.
cat(paste0("\n=========================================================================\n"))
cat(paste0("   INDEPENDENT CHECK: IRLS vs DIRECT MAXIMUM LIKELIHOOD\n"))
cat(paste0("=========================================================================\n"))
for(lnk in links_to_test) {
  funcs <- get_link_funcs(lnk)
  nll <- function(bb) {
    mu <- pmax(1e-8, pmin(1 - 1e-8, funcs$linkinv(X_aug %*% bb)))
    sum((y - mu)^2 / (y * (1 - y) * mu^2 * (1 - mu)^2))
  }
  op <- optim(beta_std_store[[lnk]], nll, method = "BFGS",
              control = list(maxit = 10000, reltol = 1e-15))
  op <- optim(op$par, nll, method = "BFGS",
              control = list(maxit = 10000, reltol = 1e-15))
  dif <- max(abs(beta_std_store[[lnk]] - op$par))
  cat(sprintf("   %-8s  max |IRLS - direct ML| = %.2e   %s\n",
              toupper(lnk), dif,
              ifelse(dif < 1e-3, "AGREE", "*** DISAGREE ***")))
}
cat("\n   (compared on the standardised scale, which is where both are fitted)\n\n")
