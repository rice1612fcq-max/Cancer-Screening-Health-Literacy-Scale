

suppressPackageStartupMessages({library(psych); library(GPArotation)})
args <- commandArgs(trailingOnly = TRUE); W <- args[1]; B_SIM <- 500L; SEED <- 20261103L
P <- readRDS(file.path(W, "prep.rds")); keep <- readRDS(file.path(W, "f02_kept_items.rds"))
X <- P$d[P$d$sample == "dev", keep]; N <- nrow(X); K <- ncol(X)
rho <- withCallingHandlers(polychoric(X, smooth = TRUE, correct = 0.5)$rho, warning = function(w) invokeRestart("muffleWarning"))
dimnames(rho) <- list(keep, keep)
R_obs <- cor(X); kmo <- KMO(R_obs)$MSA; bt <- cortest.bartlett(R_obs, n = N)
kmo_p <- KMO(rho)$MSA; bt_p <- cortest.bartlett(rho, n = N)
Rr <- rho; diag(Rr) <- smc(rho); ev_red <- sort(eigen(Rr, TRUE, TRUE)$values, TRUE)
thr <- lapply(X, function(v) qnorm(cumsum(prop.table(table(v)))[-length(table(v))]))
sim_one <- function(i) {
  Z <- matrix(rnorm(N * K), N); Y <- sapply(seq_len(K), function(j) 1L + rowSums(outer(Z[, j], thr[[j]], ">")))
  Rs <- suppressWarnings(psych::polychoric(Y, smooth = TRUE, correct = 0.5)$rho); diag(Rs) <- suppressWarnings(psych::smc(Rs))
  sort(eigen(Rs, TRUE, TRUE)$values, TRUE)
}
cl <- parallel::makeCluster(min(6L, parallel::detectCores() - 1L)); parallel::clusterSetRNGStream(cl, SEED)
parallel::clusterExport(cl, c("N", "K", "thr", "sim_one")); invisible(parallel::clusterEvalQ(cl, library(psych)))
sim <- do.call(rbind, parallel::parLapply(cl, seq_len(B_SIM), sim_one)); parallel::stopCluster(cl)
q95 <- apply(sim, 2, quantile, .95); nf <- which(ev_red <= q95)[1] - 1L
write.csv(data.frame(k = seq_len(K), reduced_eigen = ev_red, random_q95 = q95), file.path(W, "f03_parallel.csv"), row.names = FALSE)
cat("KMO (Pearson) =", round(kmo, 3), "; Bartlett chi2 =", round(bt$chisq, 1), "df =", bt$df, "\n")
cat("KMO (polychoric) =", round(kmo_p, 3), "; Bartlett chi2 (polychoric) =", round(bt_p$chisq, 1), "\n")
cat("PA factors:", nf, "| first 8 reduced eigen:", round(ev_red[1:8], 2), "| q95:", round(q95[1:8], 2), "\n")
efa_iter <- function(k) {
  it <- keep; log <- character()
  repeat {
    f <- fa(rho[it, it], k, n.obs = N, fm = "pa", rotate = "oblimin", max.iter = 500, warnings = FALSE)
    L <- abs(unclass(f$loadings)); s <- t(apply(L, 1, sort, decreasing = TRUE))
    bad <- s[, 1] <= .4 | (s[, 2] > .4 & (s[, 1] - s[, 2]) < .2)
    if (!any(bad)) break
    sev <- ifelse(s[, 1] <= .4, .4 - s[, 1] + 1, .2 - (s[, 1] - s[, 2])); w <- names(which.max(sev * bad))
    log <- c(log, w); it <- setdiff(it, w)
  }
  list(fit = f, items = it, removed = log)
}
e <- efa_iter(nf)
L <- unclass(e$fit$loadings); pf <- apply(abs(L), 1, which.max)
cat("\nEFA k =", nf, "| removed:", if (length(e$removed)) e$removed else "none", "\n")
for (j in seq_len(nf)) cat(" F", j, "(", sum(pf == j), "):", rownames(L)[pf == j], "\n")
cat("Variance explained (cum):", round(tail(e$fit$Vaccounted["Cumulative Var", ], 1), 3), "\n")
cat("Factor correlations:\n"); print(round(e$fit$Phi, 2))
write.csv(data.frame(item = rownames(L), round(L, 3), primary = paste0("F", pf), h2 = round(e$fit$communality, 3)),
          file.path(W, "f03_pattern.csv"), row.names = FALSE)
saveRDS(list(nf = nf, items = e$items, removed = e$removed, assign = split(rownames(L), pf), kmo = kmo, bt = bt),
        file.path(W, "f03_efa.rds"))
