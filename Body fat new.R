# Install required packages if missing
if(!require(mfp)) install.packages("mfp")

library(mfp)

# ==========================================
# 1. DATA PREPARATION & TRUNCATION
# ==========================================
data(bodyfat)
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
  } else if(link_name == "nloglog") { 
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
links_to_test <- c("logit", "probit", "cloglog", "nloglog")
all_tables <- list()
cn_values <- list()
ci_values <- list() # Added list to store condition indices

mse_summary <- data.frame(
  Link = character(), Condition_Num = numeric(), SMLE = numeric(), 
  SRE_k1 = numeric(), SKL_k1 = numeric(), SRE_k2 = numeric(), 
  SKL_k2 = numeric(), stringsAsFactors = FALSE
)

for(lnk in links_to_test) {
  
  funcs <- get_link_funcs(lnk)
  y_link_val <- funcs$linkfun(y)
  mod_lm <- lm(y_link_val ~ X_std)
  b_curr <- coef(mod_lm)
  sigma_curr <- sum(residuals(mod_lm)^2)/(n-p_aug)
  
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
    
    w_ii <- 3 * sigma_curr * P_val + (1 / P_val)
    W <- diag(as.vector(w_ii))
    
    u_adj <- ((y - mu_it) / P_val) * (u_dev + (1 / P_val^2))
    y_star <- eta_it + (u_adj * d_mu) / (sigma_curr * w_ii)
    
    XtWX <- t(X_aug) %*% W %*% X_aug
    XtWy_star <- t(X_aug) %*% W %*% y_star
    b_new <- solve(XtWX, XtWy_star)
    
    if(max(abs(b_new - b_curr)) < 1e-8) break 
    b_curr <- b_new
  }
  
  beta_mle_std <- as.vector(b_curr)
  sigma_est <- sigma_curr 
  
  eta_final <- X_aug %*% beta_mle_std
  mu_final <- pmax(1e-7, pmin(1-1e-7, funcs$linkinv(eta_final)))
  
  P_final <- mu_final * (1 - mu_final)
  w_final <- 3 * sigma_est * P_final + (1 / P_final)
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
  
  # Restrict parameters to theoretical limits to guarantee positive-definiteness
  k1 <- min(k1_raw, k_max_theoretical)
  k2 <- min(k2_raw, k_max_theoretical)
  
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
  
  res_table <- data.frame(
    Variable = c("Intercept", keep_vars, "MSE"),
    SMLE = c(b_mle_raw, mse_mle),
    SRE_k1 = c(b_sre_k1_raw, mse_sre_k1),
    SKL_k1 = c(b_skl_k1_raw, mse_skl_k1),
    SRE_k2 = c(b_sre_k2_raw, mse_sre_k2),
    SKL_k2 = c(b_skl_k2_raw, mse_skl_k2)
  )
  all_tables[[lnk]] <- res_table
  
  mse_summary <- rbind(mse_summary, data.frame(
    Link = toupper(lnk), Condition_Num = cond_num, SMLE = mse_mle, 
    SRE_k1 = mse_sre_k1, SKL_k1 = mse_skl_k1, SRE_k2 = mse_sre_k2, SKL_k2 = mse_skl_k2
  ))
}

# ==========================================
# 7. OUTPUT RESULTS
# ==========================================
for(lnk in links_to_test) {
  cat(paste0("\n=========================================================================\n"))
  cat(paste0("   LINK FUNCTION: ", toupper(lnk), "\n"))
  cat(paste0("   WEIGHTED CONDITION NUMBER: ", round(cn_values[[lnk]], 4), "\n"))
  cat(paste0("   CONDITION INDICES:\n   ", paste(round(ci_values[[lnk]], 4), collapse = ", "), "\n"))
  cat(paste0("=========================================================================\n"))
  print(all_tables[[lnk]], row.names=FALSE, digits=4)
}

cat(paste0("\n=========================================================================\n"))
cat(paste0("   SUMMARY OF WEIGHTED DIAGNOSTICS & MEAN SQUARED ERRORS\n"))
cat(paste0("=========================================================================\n"))
print(mse_summary, row.names=FALSE, digits=4)