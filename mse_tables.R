# --- mse_tables.R ---
#
# Writes the estimated MSE tables for p = 8 at sigma^2 = 0.5, 1 and 1.5 from
# Simulation_Results_p8.csv, one CSV per link function and dispersion level.
# The manuscript prints the four tables at sigma^2 = 2; these files give the
# same quantities at the three remaining dispersion levels.

# The CSV is read from the working directory by default, so the script runs
# from a clone with no editing.  Set SKLE_DATA to read it from elsewhere and
# SKLE_TABLES to write the output elsewhere:
#   SKLE_DATA=/path/to/csvs SKLE_TABLES=/path/to/tables Rscript mse_tables.R
data_dir  <- Sys.getenv("SKLE_DATA",   unset = ".")
table_dir <- Sys.getenv("SKLE_TABLES", unset = "./tables")

csv_path <- file.path(data_dir, "Simulation_Results_p8.csv")
if (!file.exists(csv_path)) {
  stop("Simulation_Results_p8.csv not found in ", normalizePath(data_dir))
}

df <- read.csv(csv_path, header = TRUE, stringsAsFactors = FALSE)

# Drop an exported row-number column if there is one.
if (colnames(df)[1] %in% c("X", "")) df <- df[, -1]

# The simulation writes the SMLE column as MLE.
if ("MLE" %in% colnames(df)) colnames(df)[colnames(df) == "MLE"] <- "SMLE"

df$Link <- tolower(trimws(df$Link))

est_cols <- c("SMLE",
              "SRRE_k1", "SKLE_k1", "SRRE_k2", "SKLE_k2",
              "SRRE_k3", "SKLE_k3", "SRRE_k4", "SKLE_k4")

required <- c("Link", "Sigma", "Rho", "n", est_cols)
missing  <- setdiff(required, colnames(df))
if (length(missing) > 0) {
  stop("missing columns in Simulation_Results_p8.csv: ",
       paste(missing, collapse = ", "))
}

if (!dir.exists(table_dir)) dir.create(table_dir, recursive = TRUE)

# Rows are ordered by sample size and then by correlation, as in the
# manuscript tables.  The dispersion level at sigma^2 = 2 is printed in the
# manuscript and is not written here.
links  <- c("logit", "probit", "cloglog", "neglog")
sigmas <- c(0.5, 1, 1.5)
tags   <- c("0p5", "1", "1p5")
ns     <- c(25, 50, 100, 200)
rhos   <- c(0.8, 0.9, 0.95, 0.99)

for (link in links) {
  for (j in seq_along(sigmas)) {

    out <- data.frame()

    for (n_val in ns) {
      for (rho_val in rhos) {
        row <- df[df$Link == link &
                  abs(df$Sigma - sigmas[j]) < 1e-9 &
                  df$n == n_val &
                  abs(df$Rho - rho_val) < 1e-9, ]

        if (nrow(row) != 1) {
          stop("expected one row for ", link, ", sigma^2 = ", sigmas[j],
               ", n = ", n_val, ", rho = ", rho_val,
               "; found ", nrow(row))
        }

        # Four decimal places, trailing zeros retained, as the tables are
        # printed in the manuscript.
        out <- rbind(out, data.frame(
          n   = n_val,
          rho = rho_val,
          setNames(
            lapply(est_cols, function(cc)
              formatC(as.numeric(row[[cc]]), format = "f", digits = 4)),
            est_cols
          ),
          stringsAsFactors = FALSE
        ))
      }
    }

    fname <- sprintf("MSE_p8_%s_sigma%s.csv", link, tags[j])
    write.csv(out, file.path(table_dir, fname),
              row.names = FALSE, quote = FALSE)
    cat("wrote", fname, "\n")
  }
}

cat("\n", length(links) * length(sigmas), "tables written to ",
    normalizePath(table_dir), "\n", sep = "")
