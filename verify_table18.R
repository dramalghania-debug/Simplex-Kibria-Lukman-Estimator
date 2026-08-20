## verify_table18.R -----------------------------------------------------
## Recomputes every number in Table 18 of the manuscript from the four
## Simulation_Results_p*.csv files, and prints the exact column arithmetic
## used for each one.
##
## HOW TO RUN
##   Put this file in the same folder as the four CSVs, then in RStudio:
##     Session -> Set Working Directory -> To Source File Location
##     Ctrl/Cmd + Shift + Enter
##
## Base R only. Nothing to install.
## -----------------------------------------------------------------------

files <- sprintf("Simulation_Results_p%d.csv", c(4, 8, 12, 16))
miss  <- files[!file.exists(files)]
if (length(miss))
  stop("Cannot find: ", paste(miss, collapse = ", "),
       "\nPut this script in the same folder as the four CSV files.")

d <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
cat("Read", nrow(d), "scenarios from", length(files), "files\n")
if (nrow(d) != 1024) cat("  NOTE: expected 1024\n")
cat("Replications per scenario:",
    paste(sort(unique(d$Successful_Reps)), collapse = ", "), "\n\n")

K  <- c("k1", "k2", "k3", "k4")
sk <- sapply(K, function(k) d[[paste0("SKLE_", k)]])
sr <- sapply(K, function(k) d[[paste0("SRRE_", k)]])
best <- apply(sk, 1, min)          # lowest of the four SKLE variants, per scenario

q90 <- function(x) as.numeric(quantile(x, 0.90, type = 7))

out <- data.frame(
  parameter   = c("k1", "k2", "k3", "k4"),
  cond_holds  = c(sprintf("%.0f%%", 100 * mean(d$Prov_Hold_k1)),
                  sprintf("%.0f%%", 100 * mean(d$Prov_Hold_k2)),
                  "100% (by construction)",
                  "100% (by construction)"),
  srre_ratio  = round(sapply(1:4, function(j) median(sr[, j] / d$MLE)), 3),
  skle_ratio  = round(sapply(1:4, function(j) median(sk[, j] / d$MLE)), 3),
  skle_lower  = sprintf("%.1f%%", 100 * sapply(1:4, function(j) mean(sk[, j] < sr[, j]))),
  regret_med  = round(sapply(1:4, function(j) median(sk[, j] / best)), 2),
  regret_90   = round(sapply(1:4, function(j) q90(sk[, j] / best)), 2),
  sign_rev    = sprintf("%.1f%%", 100 * sapply(K, function(k) mean(d[[paste0("SignFlip_", k)]]))),
  stringsAsFactors = FALSE
)
print(out, row.names = FALSE)

cat("\n---------------------------------------------------------------\n")
cat("Where each column comes from:\n")
cat("  cond_holds  mean of Prov_Hold_k1 / Prov_Hold_k2 across scenarios.\n")
cat("              k3 and k4 are not in the CSV: both are defined in\n")
cat("              simulation.R as strictly below 2*lambda_min, so the\n")
cat("              condition holds by construction. Verify by reading the\n")
cat("              definitions of k3 and k4 in the script.\n")
cat("  srre_ratio  median over scenarios of SRRE_k / MLE.\n")
cat("  skle_ratio  median over scenarios of SKLE_k / MLE.\n")
cat("  skle_lower  proportion of scenarios with SKLE_k < SRRE_k.\n")
cat("  regret_*    SKLE_k divided by the smallest of the four SKLE\n")
cat("              columns in the same scenario; median and 90th pct.\n")
cat("  sign_rev    mean of SignFlip_k across scenarios.\n")

cat("\n---------------------------------------------------------------\n")
cat("Other figures quoted in Sections 4.2 and 4.3:\n\n")
cat(sprintf("  SKLE beats SMLE: k1 %.1f%%, k2 %.1f%%, k3 %.1f%%, k4 %.1f%%\n",
            100 * mean(sk[, 1] < d$MLE), 100 * mean(sk[, 2] < d$MLE),
            100 * mean(sk[, 3] < d$MLE), 100 * mean(sk[, 4] < d$MLE)))

cat("\n  SKLE beats SRRE at k1, by dimension:\n")
for (p in c(4, 8, 12, 16)) {
  s <- read.csv(sprintf("Simulation_Results_p%d.csv", p))
  cat(sprintf("    p=%2d: k1 %.1f%%   k2 %.1f%%\n", p,
              100 * mean(s$SKLE_k1 < s$SRRE_k1),
              100 * mean(s$SKLE_k2 < s$SRRE_k2)))
}

cat("\n  Regularity condition hold rate, by dimension:\n")
for (p in c(4, 8, 12, 16)) {
  s <- read.csv(sprintf("Simulation_Results_p%d.csv", p))
  cat(sprintf("    p=%2d: at k1 %.1f%%   at k2 %.1f%%\n", p,
              100 * mean(s$Prov_Hold_k1), 100 * mean(s$Prov_Hold_k2)))
}

cat("\n  SKLE beats SRRE, by sample size (all dimensions pooled):\n")
for (n in c(25, 50, 100, 200)) {
  s <- d[d$n == n, ]
  cat(sprintf("    n=%3d: k1 %.0f%%   k2 %.0f%%\n", n,
              100 * mean(s$SKLE_k1 < s$SRRE_k1),
              100 * mean(s$SKLE_k2 < s$SRRE_k2)))
}

## ---- Table 17: the diagnostic banding ---------------------------------
cat("\n---------------------------------------------------------------\n")
cat("Table 17, recomputed:\n\n")
hold <- c(d$Prov_Hold_k1, d$Prov_Hold_k2)
win  <- c(sk[, 1] < sr[, 1], sk[, 2] < sr[, 2])
brk  <- c(0, .05, .10, .25, .50, .75, .95, 1.0000001)
lab  <- c("Below 0.05", "0.05 to 0.10", "0.10 to 0.25", "0.25 to 0.50",
          "0.50 to 0.75", "0.75 to 0.95", "0.95 and above")
g <- cut(hold, breaks = brk, right = FALSE, labels = lab, include.lowest = TRUE)
tab <- data.frame(band = lab,
                  cases = as.integer(table(g)),
                  skle_lower = sprintf("%.1f%%", 100 * tapply(win, g, mean)),
                  stringsAsFactors = FALSE)
print(tab, row.names = FALSE)
cat("\n  total entries:", length(hold), " (should be 2048)\n")
cat("  correlation between hold rate and SKLE winning:",
    round(cor(hold, as.numeric(win)), 3), "\n")
