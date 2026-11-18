
suppressPackageStartupMessages(library(mirt))
args <- commandArgs(trailingOnly = TRUE); W <- args[1]
P <- readRDS(file.path(W, "prep.rds")); E <- readRDS(file.path(W, "f03_efa.rds"))
d <- P$d[P$d$sample == "dev", ]; kept <- E$items
num <- function(x) as.integer(sub("^[ABC]", "", x))
theory <- list(PK = grep("^A", kept, value = TRUE), CC = grep("^B", kept, value = TRUE),
  PR = intersect(kept, c("C1", "C2")), PBC = intersect(kept, paste0("C", 3:9)),
  BAR = intersect(kept, paste0("C", 10:15)), SE = intersect(kept, paste0("C", 16:24)))
theory <- theory[lengths(theory) > 0]
if (any(lengths(theory) < 2)) stop("A content dimension has a single item: ", paste(names(theory)[lengths(theory) < 2], collapse = ","))
empirical <- unname(E$assign[lengths(E$assign) > 0])
facets <- list(
  PK = list(c("A1","A2","A4"), c("A6","A14","A17"), "A7", c("A8","A13","A15","A16"),
            c("A9","A10","A11","A12","A18","A20"), c("A19","A21","A22","A23")),
  CC = list(c("B1","B8","B9","B10"), c("B2","B5","B6"), c("B3","B4","B7","B13"),
            c("B19","B20","B21","B22","B23","B24"), c("B11","B12","B25"), c("B14","B15","B16","B17","B18","B26")),
  PR = list(c("C1","C2")), PBC = list(c("C3","C4","C5"), c("C6","C7","C8","C9")),
  BAR = list(c("C10","C11"), c("C12","C13","C14","C15")), SE = list(paste0("C", 16:24)))
subst <- list(unique(c(unlist(theory))))
Q3fun <- function(m) {
  set.seed(20261111); Th <- fscores(m, QMC = TRUE, full.scores = TRUE, verbose = FALSE)
  residuals(m, type = "Q3", Theta = Th, verbose = FALSE)
}
same_group <- function(i, j) any(vapply(subst, function(g) i %in% g && j %in% g, logical(1)))
facet_min <- list(PK = rep(1, 6), CC = c(1, 1, 1, 2, 1, 1), PR = 2, PBC = c(1, 1), BAR = c(1, 1), SE = 3)
fit <- function(G) {
  it <- unlist(G)
  spec <- paste(sapply(seq_along(G), function(j) paste0("F", j, " = ", paste(match(G[[j]], it), collapse = ","))), collapse = "\n")
  spec <- paste0(spec, "\nCOV = ", paste0("F", seq_along(G), collapse = "*"))
  set.seed(20261111)
  ans <- mirt(d[, it], mirt.model(spec), itemtype = "graded", method = "MHRM", verbose = FALSE,
              technical = list(NCYCLES = 4000, set.seed = 20261111))
  if (!isTRUE(extract.mirt(ans, "converged"))) {
    saveRDS(list(G = G, model = ans), file.path(W, "f04_failed_fit.rds"))
    stop("MHRM did not converge; do not use this fit for selection")
  }
  ans
}
stopifnot(!anyDuplicated(unlist(theory)), !anyDuplicated(unlist(empirical)))
if (!setequal(unlist(theory), kept) || !setequal(unlist(empirical), kept))
  stop("Structures must use identical retained items; unmapped: ", paste(setdiff(kept, unlist(theory)), collapse = ","))
norm_part <- function(L) { L <- lapply(unname(L), sort); L[order(sapply(L, `[`, 1))] }   
same <- identical(norm_part(theory), norm_part(empirical))
m_th <- fit(theory)
if (same) {
  
  ic <- data.frame(model = "content structure (identical to EFA assignment)", k = length(theory),
    AIC = extract.mirt(m_th, "AIC"), BIC = extract.mirt(m_th, "BIC"))
  print(ic); write.csv(ic, file.path(W, "f04_structure_IC.csv"), row.names = FALSE); G <- theory
  cat("Chosen: theory structure; identical to EFA 
")
} else {
  m_em <- fit(empirical)
  ic <- data.frame(model = c("theory (HBM-separated)", "empirical EFA"), k = c(length(theory), length(empirical)),
    AIC = c(extract.mirt(m_th, "AIC"), extract.mirt(m_em, "AIC")), BIC = c(extract.mirt(m_th, "BIC"), extract.mirt(m_em, "BIC")))
  print(ic); write.csv(ic, file.path(W, "f04_structure_IC.csv"), row.names = FALSE)
  delta_ic <- c(ic$AIC[1] - ic$AIC[2], ic$BIC[1] - ic$BIC[2]); stopifnot(all(is.finite(delta_ic)))
  G <- if (all(delta_ic <= 0) && any(delta_ic < 0)) theory else if (all(delta_ic >= 0) && any(delta_ic > 0)) empirical else
    stop("AIC and BIC disagree or are tied")
  if (!identical(G, theory)) stop("EFA-based structure selected: map content facets before item reduction.")
  cat("Chosen: theory structure; differs from EFA 
")
}
item_stats <- function(m, G) {
  it <- unlist(G); cf <- coef(m, simplify = TRUE)$items; fac <- rep(seq_along(G), lengths(G))
  a <- cf[cbind(seq_along(it), fac)]; dc <- cf[, grep("^d", colnames(cf)), drop = FALSE]
  b <- -dc / a; th <- seq(-2, 2, .05)
  info <- sapply(seq_along(it), function(i) {
    dd <- dc[i, !is.na(dc[i, ])]; Ps <- sapply(dd, function(dk) plogis(a[i] * th + dk))
    Ps <- cbind(1, Ps, 0); W_ <- Ps * (1 - Ps); Pc <- Ps[, -ncol(Ps)] - Ps[, -1]
    dP <- a[i] * (W_[, -ncol(W_)] - W_[, -1]); I <- rowSums(dP^2 / pmax(Pc, 1e-12))
    sum((I[-1] + I[-length(I)]) / 2 * diff(th)) })
  data.frame(item = it, factor = names(G)[fac], a = a, b_min = apply(b, 1, min, na.rm = TRUE),
             b_max = apply(b, 1, max, na.rm = TRUE), info = info, row.names = it)
}
protected <- function(x, G) {           
  f <- names(G)[sapply(G, function(g) x %in% g)]; rem <- setdiff(G[[f]], x)
  if (length(G[[f]]) <= 3) return(TRUE)
  any(mapply(function(fs, mn) x %in% fs && length(intersect(rem, fs)) < mn, facets[[f]], facet_min[[f]]))
}
m <- m_th; log <- list(); skipped <- character(); round <- 0
repeat {
  round <- round + 1; S <- item_stats(m, G); it <- rownames(S)
  Q <- Q3fun(m); Q[lower.tri(Q, TRUE)] <- NA
  L <- S$item[S$a <= 0.8]; L <- L[order(S[L, "a"])]
  q3cut <- 0.20
  pr <- which(Q > q3cut, arr.ind = TRUE)
  pr <- pr[vapply(seq_len(nrow(pr)), function(r) same_group(it[pr[r, 1]], it[pr[r, 2]]), logical(1)), , drop = FALSE]
  pr <- pr[order(-Q[pr]), , drop = FALSE]
  cand <- NULL; why <- NULL
  for (x in L) if (!protected(x, G)) { cand <- x; why <- sprintf("a=%.2f<=0.8", S[x, "a"]); break }
  if (is.null(cand) && nrow(pr)) for (r in seq_len(nrow(pr))) {
    i <- it[pr[r, 1]]; j <- it[pr[r, 2]]; q3 <- Q[pr[r, 1], pr[r, 2]]
    Ei <- abs(S[i, "b_min"]) > 3 | abs(S[i, "b_max"]) > 3; Ej <- abs(S[j, "b_min"]) > 3 | abs(S[j, "b_max"]) > 3
    first <- if (Ei != Ej) (if (Ei) i else j) else (if (S[i, "info"] <= S[j, "info"]) i else j)
    second <- setdiff(c(i, j), first)
    for (x in c(first, second)) if (!protected(x, G)) {
      cand <- x; why <- sprintf("Q3(%s,%s)=%.3f; %s; info %s=%.2f vs %s=%.2f", i, j, q3,
        if (x != first) paste(first, "retained for content coverage") else if (Ei != Ej) "extreme |b|>3 decides" else "lower information",
        i, S[i, "info"], j, S[j, "info"]); break }
    if (!is.null(cand)) break
    skipped <- union(skipped, sprintf("%s-%s(%.3f)", i, j, q3))
  }
  if (is.null(cand)) break
  cat(sprintf("round %d: delete %s | %s\n", round, cand, why))
  log[[round]] <- data.frame(round = round, deleted = cand, reason = why, n_Q3_pairs = nrow(pr), q3_cut = q3cut)
  G <- lapply(G, setdiff, cand); m <- fit(G)
}
logdf <- do.call(rbind, log); write.csv(logdf, file.path(W, "f04_deletion_log.csv"), row.names = FALSE)
S <- item_stats(m, G); write.csv(S, file.path(W, "f04_final_item_params.csv"), row.names = FALSE)
Q <- Q3fun(m); write.csv(round(Q, 3), file.path(W, "f04_final_Q3.csv"))
m2 <- M2(m, type = "C2", QMC = TRUE); sx <- itemfit(m, fit_stats = "S_X2", QMC = TRUE)
write.csv(sx, file.path(W, "f04_final_itemfit.csv"), row.names = FALSE)
cat("\nDeleted", nrow(logdf), "items; final", length(unlist(G)), "items:", lengths(G), "\n")
for (f in names(G)) cat(" ", f, ":", G[[f]], "\n")
cat("Final fit:", paste(names(m2), round(unlist(m2[1, ]), 3), sep = "=", collapse = " "), "\n")
cat("a range:", round(range(S$a), 2), " b range:", round(c(min(S$b_min), max(S$b_max)), 2), "\n")
Qn <- Q; Qn[lower.tri(Qn, TRUE)] <- NA
cat("Remaining Q3>0.20 pairs:", sum(Qn > .2, na.rm = TRUE), " skipped (protected):", skipped, "\n")
cat("S-X2 p<0.01:", sum(sx$p.S_X2 < .01, na.rm = TRUE), "of", nrow(sx), "\n")
cat("Latent correlations:\n"); print(round(summary(m, verbose = FALSE)$fcor, 3))
cat("Converged:", extract.mirt(m, "converged"), "\n")
saveRDS(list(G = G, model = m), file.path(W, "f04_final.rds"))
