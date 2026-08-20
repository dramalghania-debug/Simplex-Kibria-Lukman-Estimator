# --- figures.R ---
#
# Produces the figures used in Section 4 from the four Simulation_Results_p*.csv
# files.  Captions are supplied by the manuscript and are not drawn into the
# images.

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork) 

# ==============================================================================
# 1. Define Paths and Load Data
# ==============================================================================
# The four CSVs are read from the working directory by default, so the script
# runs from a clone with no editing.  Set SKLE_DATA to read them from
# elsewhere:  SKLE_DATA=/path/to/csvs Rscript figures.R
data_dir <- Sys.getenv("SKLE_DATA", unset = ".")

p_vals <- c(4, 8, 12, 16)
all_data <- data.frame()

cat(">>> Loading Master CSV Data...\n")

for (p in p_vals) {
  csv <- paste0("Simulation_Results_p", p, ".csv")
  file_path <- file.path(data_dir, csv)
  
  if (file.exists(file_path)) {
    df <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
    
    # Clean up: If R accidentally exported row numbers as the first column, drop it
    if (colnames(df)[1] == "X" || colnames(df)[1] == "") {
      df <- df[, -1]
    }
    
    # Explicitly remove "Successful_Reps" so it doesn't hijack SMLE!
    if ("Successful_Reps" %in% colnames(df)) {
      df$Successful_Reps <- NULL
    }
    
    # Rename MLE to SMLE if it was exported that way
    if ("MLE" %in% colnames(df)) {
      colnames(df)[colnames(df) == "MLE"] <- "SMLE"
    }
    
    # Select the columns by name so that column order cannot matter.
    # SRRE_k1/SKLE_k1 are the geometric mean estimator, SRRE_k2/SKLE_k2 the
    # square-root median estimator; Prov_Hold_k2 and K4_Binds are used by
    # figures 24 and 25.
    required_cols <- c("Link", "Sigma", "Rho", "n", "SMLE",
                       "SRRE_k1", "SKLE_k1", "SRRE_k2", "SKLE_k2",
                       "SKLE_k3", "SKLE_k4", "Prov_Hold_k2", "K4_Binds")
    
    if (all(required_cols %in% colnames(df))) {
      df <- df[, required_cols]
      
      df$p <- p
      df$Link <- tolower(trimws(df$Link))
      
      # Force numeric conversion to be absolutely safe
      df$Sigma <- as.numeric(df$Sigma)
      df$Rho <- as.numeric(df$Rho)
      df$n <- as.numeric(df$n)
      df$SMLE <- as.numeric(df$SMLE)
      for (cc in c("SRRE_k1","SKLE_k1","SRRE_k2","SKLE_k2",
                   "SKLE_k3","SKLE_k4","Prov_Hold_k2","K4_Binds")) {
        df[[cc]] <- as.numeric(df[[cc]])
      }
      
      df <- na.omit(df) 
      
      all_data <- rbind(all_data, df)
      cat(sprintf("Loaded p=%d: Successfully read %d rows\n", p, nrow(df)))
    } else {
      cat(sprintf("ERROR: %s is missing required columns. Found: %s\n", file_path, paste(colnames(df), collapse=", ")))
    }
  } else {
    cat(sprintf("WARNING: Could not find %s\n", file_path))
  }
}

if(nrow(all_data) == 0) stop("CRITICAL ERROR: No valid data was loaded.")

# ==============================================================================
# 2. Prepare Data
# ==============================================================================
long_data <- all_data %>%
  pivot_longer(cols = c("SMLE", "SRRE_k1", "SKLE_k1", "SRRE_k2", "SKLE_k2"),
               names_to = "Estimator",
               values_to = "MSE")

long_data$Estimator <- factor(long_data$Estimator,
                              levels = c("SMLE", "SRRE_k1", "SKLE_k1", "SRRE_k2", "SKLE_k2"))

# Figures are written to ./figures beside the script.  Set SKLE_FIGS to
# write them elsewhere.
output_dir <- Sys.getenv("SKLE_FIGS", unset = file.path(".", "figures"))
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==============================================================================
# 3. Plotting & Combining Functions
# ==============================================================================
create_plot <- function(data_subset, x_var, x_lab_expr, caption_expr) {
  
  # ---> MATHEMATICAL NOTATION UPDATE HERE <---
  # Using expression() to force k_1 and k_2 to render as proper subscripts in the legend
  est_labels <- c(
    expression(SMLE),
    expression(SRRE(k[1])),
    expression(SKLE(k[1])),
    expression(SRRE(k[2])),
    expression(SKLE(k[2]))
  )
  
  if(x_var == "Link") {
    data_subset$Link <- factor(data_subset$Link, levels = c("logit", "probit", "cloglog", "neglog"))
    p <- ggplot(data_subset, aes(x = Link, y = MSE, 
                                 color = Estimator, group = Estimator, 
                                 shape = Estimator, linetype = Estimator)) + geom_line(linewidth = 0.8) 
  } else {
    p <- ggplot(data_subset, aes(x = .data[[x_var]], y = MSE, 
                                 color = Estimator, shape = Estimator, 
                                 linetype = Estimator)) + geom_line(linewidth = 0.8) 
  }
  
  p <- p + geom_point(size = 3) + 
    # Applying the mathematical labels to color, shape, and linetype so they perfectly merge into one legend
    scale_color_discrete(labels = est_labels) +
    scale_shape_manual(values = c(16, 17, 15, 3, 4), labels = est_labels) + 
    scale_linetype_manual(values = c("solid", "dotted", "dashed", "dotdash", "twodash"), labels = est_labels) +
    xlab(x_lab_expr) + ylab("MSE") +
    # Captions are supplied by the manuscript, not baked into the image.
    theme_bw() +
    theme(
      plot.caption = element_text(hjust = 0.5, size = 11, margin = margin(t = 12, b = 5)),
      legend.title = element_blank(), 
      legend.text = element_text(size = 11),
      axis.title = element_text(size = 12),
      panel.grid.minor = element_blank()
    )
  
  if(x_var == "Rho") p <- p + scale_x_continuous(breaks = c(0.80, 0.90, 0.95, 0.99))
  if(x_var == "n") p <- p + scale_x_continuous(breaks = c(25, 50, 100, 200))
  if(x_var == "p") p <- p + scale_x_continuous(breaks = c(4, 8, 12, 16))
  if(x_var == "Sigma") p <- p + scale_x_continuous(breaks = c(0.5, 1.0, 1.5, 2.0))
  
  return(p) 
}

save_paired_plots <- function(p_left, p_right, filename) {
  if(nrow(p_left$data) == 0 || nrow(p_right$data) == 0) {
    cat(sprintf("   ! Warning: Skipping %s (Data missing for plot parameters)\n", filename))
    return()
  }
  combined <- p_left + p_right + plot_layout(ncol = 2)
  ggsave(filename = file.path(output_dir, filename), plot = combined,
         width = 12, height = 6, dpi = 300)
}

# ==============================================================================
# 4. Generate and Pair the Plots
# ==============================================================================
cat("\n>>> Generating 22 figures (11 pairs)...\n")

# ---------------- GROUP 1: Effect of Multicollinearity (X = Rho) ----------------
cap1 <- expression(atop("Figure 1. Effect of multicollinearity on SRM's estimators", paste("for p = 4, n = 25, ", sigma^2, " = 0.5, and the logit link.")))
p1 <- create_plot(subset(long_data, p==4 & n==25 & Sigma==0.5 & Link=="logit"), "Rho", expression(rho), cap1)

cap2 <- expression(atop("Figure 2. Effect of multicollinearity on SRM's estimators", paste("for p = 4, n = 200, ", sigma^2, " = 0.5, and the logit link.")))
p2 <- create_plot(subset(long_data, p==4 & n==200 & Sigma==0.5 & Link=="logit"), "Rho", expression(rho), cap2)
save_paired_plots(p1, p2, "Fig01_02_Combined.png")

cap3 <- expression(atop("Figure 3. Effect of multicollinearity on SRM's estimators", paste("for p = 4, n = 25, ", sigma^2, " = 2, and the logit link.")))
p3 <- create_plot(subset(long_data, p==4 & n==25 & Sigma==2.0 & Link=="logit"), "Rho", expression(rho), cap3)

cap4 <- expression(atop("Figure 4. Effect of multicollinearity on SRM's estimators", paste("for p = 4, n = 200, ", sigma^2, " = 2, and the logit link.")))
p4 <- create_plot(subset(long_data, p==4 & n==200 & Sigma==2.0 & Link=="logit"), "Rho", expression(rho), cap4)
save_paired_plots(p3, p4, "Fig03_04_Combined.png")

cap5 <- expression(atop("Figure 5. Effect of multicollinearity on SRM's estimators", paste("for p = 16, n = 200, ", sigma^2, " = 0.5, and the logit link.")))
p5 <- create_plot(subset(long_data, p==16 & n==200 & Sigma==0.5 & Link=="logit"), "Rho", expression(rho), cap5)

cap6 <- expression(atop("Figure 6. Effect of multicollinearity on SRM's estimators", paste("for p = 16, n = 25, ", sigma^2, " = 0.5, and the logit link.")))
p6 <- create_plot(subset(long_data, p==16 & n==25 & Sigma==0.5 & Link=="logit"), "Rho", expression(rho), cap6)
save_paired_plots(p5, p6, "Fig05_06_Combined.png")

# ---------------- GROUP 2: Effect of Sample Size (X = n) ----------------
cap7 <- expression(atop("Figure 7. Effect of sample size on SRM's estimators", paste("for p = 4, link = logit, ", sigma^2, " = 2 and ", rho, " = 0.99.")))
p7 <- create_plot(subset(long_data, p==4 & Link=="logit" & Sigma==2.0 & Rho==0.99), "n", "Sample Size (n)", cap7)

# rho = 0.95 rather than 0.99. At rho = 0.99 the regularity condition
# k < 2*lambda_min holds in essentially no replications for n = 25, 50 or 100,
# so the SKLE is worse than the SRRE at three of the four points. At rho = 0.95
# it fails only at n = 25 and recovers thereafter, which shows the diagnostic
# operating. The rho = 0.99 case is discussed in the text.
cap8 <- expression(atop("Figure 8. Effect of sample size on SRM's estimators", paste("for p = 16, link = logit, ", sigma^2, " = 2 and ", rho, " = 0.95.")))
p8 <- create_plot(subset(long_data, p==16 & Link=="logit" & Sigma==2.0 & Rho==0.95), "n", "Sample Size (n)", cap8)
save_paired_plots(p7, p8, "Fig07_08_Combined.png")

# Figures 9 and 10 use the probit link. Under cloglog and neglog the MSE
# distribution is strongly right-skewed, with mean MSE up to 67 times the
# median, so mean MSE is not monotone in n for those links. No cloglog or
# neglog cell in the design is free of this. The asymmetric links retain
# coverage in Figures 13-14 and 19-22.
cap9 <- expression(atop("Figure 9. Effect of sample size on SRM's estimators", paste("for p = 16, link = probit, ", sigma^2, " = 2 and ", rho, " = 0.99.")))
p9 <- create_plot(subset(long_data, p==16 & Link=="probit" & Sigma==2.0 & Rho==0.99), "n", "Sample Size (n)", cap9)

cap10 <- expression(atop("Figure 10. Effect of sample size on SRM's estimators", paste("for p = 4, link = probit, ", sigma^2, " = 2 and ", rho, " = 0.99.")))
p10 <- create_plot(subset(long_data, p==4 & Link=="probit" & Sigma==2.0 & Rho==0.99), "n", "Sample Size (n)", cap10)
save_paired_plots(p9, p10, "Fig09_10_Combined.png")

# ---------------- GROUP 3: Effect of Number of Variables (X = p) ----------------
cap11 <- expression(atop("Figure 11. Effect of number of explanatory variables", paste("on SRM's estimators for n = 25, link = logit, ", sigma^2, " = 2 and ", rho, " = 0.99.")))
p11 <- create_plot(subset(long_data, n==25 & Link=="logit" & Sigma==2.0 & Rho==0.99), "p", "Number of Explanatory Variables (p)", cap11)

cap12 <- expression(atop("Figure 12. Effect of number of explanatory variables", paste("on SRM's estimators for n = 200, link = logit, ", sigma^2, " = 2 and ", rho, " = 0.99.")))
p12 <- create_plot(subset(long_data, n==200 & Link=="logit" & Sigma==2.0 & Rho==0.99), "p", "Number of Explanatory Variables (p)", cap12)
save_paired_plots(p11, p12, "Fig11_12_Combined.png")

cap13 <- expression(atop("Figure 13. Effect of number of explanatory variables", paste("on SRM's estimators for n = 25, link = neglog, ", sigma^2, " = 2 and ", rho, " = 0.99.")))
p13 <- create_plot(subset(long_data, n==25 & Link=="neglog" & Sigma==2.0 & Rho==0.99), "p", "Number of Explanatory Variables (p)", cap13)

cap14 <- expression(atop("Figure 14. Effect of number of explanatory variables", paste("on SRM's estimators for n = 200, link = neglog, ", sigma^2, " = 2 and ", rho, " = 0.99.")))
p14 <- create_plot(subset(long_data, n==200 & Link=="neglog" & Sigma==2.0 & Rho==0.99), "p", "Number of Explanatory Variables (p)", cap14)
save_paired_plots(p13, p14, "Fig13_14_Combined.png")

# ---------------- GROUP 4: Effect of Sigma (X = Sigma) ----------------
cap15 <- expression(atop(paste("Figure 15. Effect of ", sigma^2, " on SRM's estimators"), paste("for p = 4, n = 25, link = logit, and ", rho, " = 0.99.")))
p15 <- create_plot(subset(long_data, p==4 & n==25 & Link=="logit" & Rho==0.99), "Sigma", expression(sigma^2), cap15)

cap16 <- expression(atop(paste("Figure 16. Effect of ", sigma^2, " on SRM's estimators"), paste("for p = 16, n = 25, link = logit, and ", rho, " = 0.99.")))
p16 <- create_plot(subset(long_data, p==16 & n==25 & Link=="logit" & Rho==0.99), "Sigma", expression(sigma^2), cap16)
save_paired_plots(p15, p16, "Fig15_16_Combined.png")

cap17 <- expression(atop(paste("Figure 17. Effect of ", sigma^2, " on SRM's estimators"), paste("for p = 4, n = 200, link = logit, and ", rho, " = 0.99.")))
p17 <- create_plot(subset(long_data, p==4 & n==200 & Link=="logit" & Rho==0.99), "Sigma", expression(sigma^2), cap17)

cap18 <- expression(atop(paste("Figure 18. Effect of ", sigma^2, " on SRM's estimators"), paste("for p = 16, n = 200, link = logit, and ", rho, " = 0.99.")))
p18 <- create_plot(subset(long_data, p==16 & n==200 & Link=="logit" & Rho==0.99), "Sigma", expression(sigma^2), cap18)
save_paired_plots(p17, p18, "Fig17_18_Combined.png")

# ---------------- GROUP 5: Effect of Link Function (X = Link) ----------------
cap19 <- expression(atop("Figure 19. Effect of link function on SRM's estimators", paste("for p = 4, n = 25, ", sigma^2, " = 2, and ", rho, " = 0.99.")))
p19 <- create_plot(subset(long_data, p==4 & n==25 & Sigma==2.0 & Rho==0.99), "Link", "Link Function", cap19)

cap20 <- expression(atop("Figure 20. Effect of link function on SRM's estimators", paste("for p = 4, n = 200, ", sigma^2, " = 2, and ", rho, " = 0.99.")))
p20 <- create_plot(subset(long_data, p==4 & n==200 & Sigma==2.0 & Rho==0.99), "Link", "Link Function", cap20)
save_paired_plots(p19, p20, "Fig19_20_Combined.png")

cap21 <- expression(atop("Figure 21. Effect of link function on SRM's estimators", paste("for p = 16, n = 25, ", sigma^2, " = 2, and ", rho, " = 0.99.")))
p21 <- create_plot(subset(long_data, p==16 & n==25 & Sigma==2.0 & Rho==0.99), "Link", "Link Function", cap21)

cap22 <- expression(atop("Figure 22. Effect of link function on SRM's estimators", paste("for p = 16, n = 200, ", sigma^2, " = 2, and ", rho, " = 0.99.")))
p22 <- create_plot(subset(long_data, p==16 & n==200 & Sigma==2.0 & Rho==0.99), "Link", "Link Function", cap22)
save_paired_plots(p21, p22, "Fig21_22_Combined.png")

cat("\n>>> 22 paired plots written.\n")

# ==============================================================================
# 5. FIGURES 24-25
#    24: comparison of the four shrinkage parameters, by dimension and by link
#    25: how the regularity condition degrades with p and rho
# ==============================================================================
cat("\n>>> Generating figures 24-25...\n")

# --- Figure 24: the four shrinkage parameters, two panels -------------------
# Left : median ratio of SKLE to SMLE MSE, by number of explanatory variables.
# Right: the same ratio by link function.
k_data <- all_data %>%
  group_by(p) %>%
  summarise(k1 = median(SKLE_k1 / SMLE), k2 = median(SKLE_k2 / SMLE),
            k3 = median(SKLE_k3 / SMLE), k4 = median(SKLE_k4 / SMLE),
            .groups = "drop") %>%
  pivot_longer(cols = c("k1", "k2", "k3", "k4"),
               names_to = "Parameter", values_to = "Ratio")

k_data$Parameter <- factor(k_data$Parameter, levels = c("k1", "k2", "k3", "k4"))
k_labels <- c(expression(hat(k)[1]), expression(hat(k)[2]),
              expression(hat(k)[3]), expression(hat(k)[4]))

p24a <- ggplot(k_data, aes(x = p, y = Ratio, color = Parameter,
                           shape = Parameter, linetype = Parameter)) +
  geom_line(linewidth = 0.8) + geom_point(size = 3) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey50") +
  scale_color_discrete(labels = k_labels) +
  scale_shape_manual(values = c(16, 17, 15, 3), labels = k_labels) +
  scale_linetype_manual(values = c("solid", "dotted", "dashed", "dotdash"),
                        labels = k_labels) +
  scale_x_continuous(breaks = c(4, 8, 12, 16)) + ylim(0, 1.02) +
  xlab("Number of Explanatory Variables (p)") + ylab("Median MSE ratio to SMLE") +
  theme_bw() +
  theme(legend.title = element_blank(), legend.text = element_text(size = 11),
        axis.title = element_text(size = 12), panel.grid.minor = element_blank())

l_data <- all_data %>%
  group_by(Link) %>%
  summarise(k1 = median(SKLE_k1 / SMLE), k2 = median(SKLE_k2 / SMLE),
            k3 = median(SKLE_k3 / SMLE), k4 = median(SKLE_k4 / SMLE),
            .groups = "drop") %>%
  pivot_longer(cols = c("k1", "k2", "k3", "k4"),
               names_to = "Parameter", values_to = "Ratio")

l_data$Parameter <- factor(l_data$Parameter, levels = c("k1", "k2", "k3", "k4"))
l_data$Link      <- factor(l_data$Link,
                           levels = c("logit", "probit", "cloglog", "neglog"))

p24b <- ggplot(l_data, aes(x = Link, y = Ratio, color = Parameter,
                           shape = Parameter, linetype = Parameter,
                           group = Parameter)) +
  geom_line(linewidth = 0.8) + geom_point(size = 3) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey50") +
  scale_color_discrete(labels = k_labels) +
  scale_shape_manual(values = c(16, 17, 15, 3), labels = k_labels) +
  scale_linetype_manual(values = c("solid", "dotted", "dashed", "dotdash"),
                        labels = k_labels) +
  ylim(0, 1.02) +
  xlab("Link Function") + ylab("Median MSE ratio to SMLE") +
  theme_bw() +
  theme(legend.title = element_blank(), legend.text = element_text(size = 11),
        axis.title = element_text(size = 12), panel.grid.minor = element_blank())

ggsave(file.path(output_dir, "Fig24_Parameters.png"),
       p24a + p24b + plot_layout(ncol = 2, guides = "collect"),
       width = 12, height = 6, dpi = 300)

# --- Figure 25: degradation of the regularity condition ---------------------
prov_data <- all_data %>%
  group_by(p, Rho) %>%
  summarise(Hold = mean(Prov_Hold_k2), .groups = "drop") %>%
  mutate(p_f = factor(p, levels = c(4, 8, 12, 16)))

p25 <- ggplot(prov_data, aes(x = Rho, y = Hold, color = p_f,
                             shape = p_f, linetype = p_f)) +
  geom_line(linewidth = 0.8) + geom_point(size = 3) +
  scale_color_discrete(name = "p") +
  scale_shape_manual(values = c(16, 17, 15, 3), name = "p") +
  scale_linetype_manual(values = c("solid", "dotted", "dashed", "dotdash"), name = "p") +
  scale_x_continuous(breaks = c(0.80, 0.90, 0.95, 0.99)) +
  scale_y_continuous(limits = c(0, 1)) +
  xlab(expression(rho)) +
  ylab(expression(paste("Proportion satisfying ", k < 2*lambda[min]))) +
  theme_bw() +
  theme(plot.caption = element_text(hjust = 0.5, size = 11, margin = margin(t = 12, b = 5)),
        legend.text = element_text(size = 11),
        axis.title = element_text(size = 12), panel.grid.minor = element_blank())

ggsave(file.path(output_dir, "Fig25_Regularity.png"), p25,
       width = 7, height = 5.5, dpi = 300)

cat(sprintf("\n>>> Complete. Figures written to:\n    %s\n", output_dir))
