
args <- commandArgs(trailingOnly = TRUE); W <- args[1]
P <- readRDS(file.path(W, "prep.rds")); d <- P$d[P$d$sample == "dev", ]; items <- P$items
X <- d[, items]; stopifnot(nrow(X) > 1L)
alpha <- function(m) { k <- ncol(m); (k/(k-1)) * (1 - sum(apply(m, 2, var)) / var(rowSums(m))) }
tot <- rowSums(X); q <- quantile(tot, c(.27, .73), type = 7)
lo <- tot <= q[1]; hi <- tot >= q[2]
cat("27% groups: low n =", sum(lo), " high n =", sum(hi), "\n")
res <- do.call(rbind, lapply(items, function(i) {
  D <- X[, startsWith(items, substr(i, 1, 1))]
  tt <- t.test(X[[i]][hi], X[[i]][lo])
  data.frame(item = i, dimension = substr(i, 1, 1),
             r_dim = cor(X[[i]], rowSums(D)), p_r = cor.test(X[[i]], rowSums(D))$p.value,
             r_total = cor(X[[i]], tot), CR_t = unname(tt$statistic), CR_p = tt$p.value,
             alpha_dim = alpha(D), alpha_if_deleted = alpha(D[, setdiff(names(D), i)]))
}))
res$fail_r <- abs(res$r_dim) < 0.30
res$fail_CR <- !(res$CR_t > 3 & res$CR_p < 0.05)
res$fail_alpha <- res$alpha_if_deleted > res$alpha_dim


keep_content <- if (length(args) >= 2 && startsWith(args[2], "keep=")) strsplit(sub("^keep=", "", args[2]), ",")[[1]] else character()
res$retained_content <- res$item %in% keep_content & res$fail_alpha & !res$fail_r & !res$fail_CR
res$deleted <- (res$fail_r | res$fail_CR | res$fail_alpha) & !res$retained_content
write.csv(res, file.path(W, "f02_ctt_items.csv"), row.names = FALSE)
cat("Dimension alpha:", tapply(res$alpha_dim, res$dimension, unique), "\n")
cat("r_dim range:", round(range(res$r_dim), 3), " CR t range:", round(range(res$CR_t), 2), "\n")
print(res[res$deleted, c("item", "r_dim", "CR_t", "CR_p", "alpha_dim", "alpha_if_deleted", "fail_r", "fail_CR", "fail_alpha")], row.names = FALSE, digits = 4)
saveRDS(res$item[!res$deleted], file.path(W, "f02_kept_items.rds"))
cat("Kept:", sum(!res$deleted), "items\n")
