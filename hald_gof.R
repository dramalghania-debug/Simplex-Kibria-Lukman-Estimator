suppressMessages(library(goftest))
y_raw <- c(78.5,74.3,104.3,87.6,95.9,109.2,102.7,72.5,93.1,115.9,83.8,113.3,109.4)
n <- length(y_raw)
y <- pmax(1e-5, pmin(1-1e-5, y_raw/120))

dsimplex <- function(v, mu, sig) {
  d <- (v-mu)^2/(v*(1-v)*mu^2*(1-mu)^2)
  (2*pi*sig*(v*(1-v))^3)^(-1/2)*exp(-d/(2*sig))
}
psimplex <- function(q, mu, sig) sapply(q, function(v) {
  if (v <= 0) return(0); if (v >= 1) return(1)
  integrate(dsimplex, 1e-9, v, mu=mu, sig=sig, stop.on.error=FALSE)$value })
qsimplex <- function(p, mu, sig) sapply(p, function(pp)
  uniroot(function(v) psimplex(v,mu,sig)-pp, c(1e-6,1-1e-6), tol=1e-10)$root)
fit <- function(v) { m <- mean(v)
  list(mu=m, sig=mean((v-m)^2/(v*(1-v)*m^2*(1-m)^2))) }

f <- fit(y); mu <- f$mu; sig <- f$sig
cat(sprintf("n = %d   mu-hat = %.5f   sigma2-hat = %.5f\n\n", n, mu, sig))

u <- psimplex(y, mu, sig)
D  <- ks.test(y, psimplex, mu=mu, sig=sig)
cvm <- cvm.test(u, "punif")
ad  <- ad.test(u, "punif")
cat(sprintf("KS   D  = %.5f   asymptotic p = %.4f\n", D$statistic, D$p.value))
cat(sprintf("CvM  W2 = %.5f   asymptotic p = %.4f\n", cvm$statistic, cvm$p.value))
cat(sprintf("AD   A2 = %.5f   asymptotic p = %.4f\n\n", ad$statistic, ad$p.value))

cat("Chi-square, equiprobable bins, expected count = n/k:\n")
for (k in 2:4) {
  cut_pts <- qsimplex((1:(k-1))/k, mu, sig)
  obs <- as.vector(table(cut(y, c(-Inf, cut_pts, Inf))))
  exp_ct <- rep(n/k, k)
  X2 <- sum((obs-exp_ct)^2/exp_ct)
  cat(sprintf("  k=%d  expected=%.2f  obs=%s  X2=%.4f  df=%d  p=%.4f  (df=k-1-2: %s)\n",
      k, n/k, paste(obs, collapse=","), X2, k-1, pchisq(X2, k-1, lower.tail=FALSE),
      if (k-3 > 0) sprintf("p=%.4f", pchisq(X2, k-3, lower.tail=FALSE)) else "df<=0"))
}

cat("\nParametric bootstrap, parameters re-estimated on each sample (B = 4000):\n")
set.seed(20260815)
rsimplex <- function(m, mu, sig) qsimplex(runif(m), mu, sig)
B <- 4000
bs <- replicate(B, {
  ys <- pmax(1e-5, pmin(1-1e-5, rsimplex(n, mu, sig)))
  g <- fit(ys)
  us <- psimplex(ys, g$mu, g$sig)
  c(ks.test(ys, psimplex, mu=g$mu, sig=g$sig)$statistic,
    cvm.test(us, "punif")$statistic)
})
cat(sprintf("  KS   bootstrap p = %.4f   (asymptotic %.4f)\n", mean(bs[1,] >= D$statistic), D$p.value))
cat(sprintf("  CvM  bootstrap p = %.4f   (asymptotic %.4f)\n", mean(bs[2,] >= cvm$statistic), cvm$p.value))
