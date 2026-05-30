## ================= Load the 6 chains (300k each) =================
load("./outputS30.RData") # Mod30
load("./outputS31.RData") # Mod31
load("./Revisor USPP\\outputS32.RData") # Mod32
load("./outputS33.RData") # Mod33
load("./Revisor USPP\\outputS34.RData") # Mod34
load("./Revisor USPP\\outputS35.RData") # Mod35

## ================= Convergence (coda): 6 chains, MBeta and MPsi =================
if (!requireNamespace("coda", quietly = TRUE)) install.packages("coda")
library(coda)

# ---------- Helpers ----------
align_columns_list <- function(mats, what = "MBeta") {
  stopifnot(length(mats) >= 2)
  cn1 <- colnames(mats[[1]])
  if (is.null(cn1)) {
    k <- ncol(mats[[1]])
    ok <- vapply(mats, function(M) ncol(M) == k, TRUE)
    if (!all(ok)) stop(sprintf("%s: número de colunas difere entre cadeias.", what))
    return(mats)
  }
  for (i in seq_along(mats)) {
    cni <- colnames(mats[[i]])
    if (is.null(cni) || !setequal(cni, cn1))
      stop(sprintf("%s: conjuntos de colunas diferentes entre cadeias.", what))
    mats[[i]] <- mats[[i]][, cn1, drop = FALSE]
  }
  mats
}

equalize_lengths <- function(mats) {
  min_n <- min(vapply(mats, nrow, 1L))
  lapply(mats, function(M) M[(nrow(M) - min_n + 1):nrow(M), , drop = FALSE])
}

analyze_component_all <- function(mods, comp = c("MBeta", "MPsi"), label = NULL) {
  comp  <- match.arg(comp)
  label <- if (is.null(label)) comp else label

  stopifnot(length(mods) >= 2)
  has_comp <- vapply(mods, function(m) comp %in% names(m), TRUE)
  if (!all(has_comp)) {
    cat(sprintf("[Info] '%s' ausente em alguma cadeia; %s não analisado.\n", comp, label))
    return(invisible(NULL))
  }

  mats <- lapply(mods, function(m) as.matrix(m[[comp]]))
  mats <- align_columns_list(mats, what = comp)
  mats <- equalize_lengths(mats)

  if (is.null(colnames(mats[[1]]))) {
    colnames(mats[[1]]) <- paste0(label, ".", seq_len(ncol(mats[[1]])))
    for (i in 2:length(mats)) colnames(mats[[i]]) <- colnames(mats[[1]])
  }

  mlist <- mcmc.list(lapply(mats, mcmc))
  gr    <- gelman.diag(mlist, autoburnin = FALSE)

  tab <- data.frame(
    Param        = colnames(mats[[1]]),
    Rhat         = as.numeric(gr$psrf[, 1]),
    Rhat_UpperCI = as.numeric(gr$psrf[, 2]),
    stringsAsFactors = FALSE
  )
  mpsrf <- unname(gr$mpsrf)

  cat("\n=================================================================\n")
  cat(sprintf("== %s: Gelman–Rubin com %d cadeias (todas as amostras) ==\n", label, length(mods)))
  cat(sprintf("Parâmetros: %d | Multivariate PSRF: %.4f\n\n", nrow(tab), mpsrf))
  print(tab, row.names = FALSE)

  invisible(list(table = tab, mpsrf = mpsrf, mlist = mlist))
}

# ---------- Run for 6 chains ----------
mods6 <- list(Mod30, Mod31, Mod32, Mod33, Mod34, Mod35)

res_beta <- analyze_component_all(mods6, "MBeta", "β (MBeta)")
res_psi  <- analyze_component_all(mods6, "MPsi",  "Ψ (MPsi)")

## ========= Prerequisites (same as already used) =========
if (!requireNamespace("coda", quietly = TRUE)) install.packages("coda")
library(coda)

## ---------- Helper functions (same) ----------
align_columns_list <- function(mats, what = "MBeta") {
  stopifnot(length(mats) >= 2)
  cn1 <- colnames(mats[[1]])
  if (is.null(cn1)) {
    k <- ncol(mats[[1]])
    ok <- vapply(mats, function(M) ncol(M) == k, TRUE)
    if (!all(ok)) stop(sprintf("%s: number of columns differ between chains.", what))
    return(mats)
  }
  for (i in seq_along(mats)) {
    cni <- colnames(mats[[i]])
    if (is.null(cni) || !setequal(cni, cn1))
      stop(sprintf("%s: different sets of columns between chains.", what))
    mats[[i]] <- mats[[i]][, cn1, drop = FALSE]
  }
  mats
}

equalize_lengths <- function(mats) {
  min_n <- min(vapply(mats, nrow, 1L))
  lapply(mats, function(M) M[(nrow(M) - min_n + 1):nrow(M), , drop = FALSE])
}

analyze_component_all <- function(mods, comp = c("MBeta","MPsi","Mbw","Mbm","Mvw","Mvm"), label = NULL) {
  comp  <- match.arg(comp)
  label <- if (is.null(label)) comp else label

  stopifnot(length(mods) >= 2)
  has_comp <- vapply(mods, function(m) comp %in% names(m), TRUE)
  if (!all(has_comp)) {
    cat(sprintf("[Info] '%s' missing in some chain; %s not analyzed.\n", comp, label))
    return(invisible(NULL))
  }

  mats <- lapply(mods, function(m) as.matrix(m[[comp]]))
  mats <- align_columns_list(mats, what = comp)
  mats <- equalize_lengths(mats)

  if (is.null(colnames(mats[[1]]))) {
    colnames(mats[[1]]) <- paste0(label, ".", seq_len(ncol(mats[[1]])))
    for (i in 2:length(mats)) colnames(mats[[i]]) <- colnames(mats[[1]])
  }

  mlist <- mcmc.list(lapply(mats, mcmc))
  gr    <- gelman.diag(mlist, autoburnin = FALSE)

  tab <- data.frame(
    Param        = colnames(mats[[1]]),
    Rhat         = as.numeric(gr$psrf[, 1]),
    Rhat_UpperCI = as.numeric(gr$psrf[, 2]),
    stringsAsFactors = FALSE
  )
  mpsrf <- unname(gr$mpsrf)

  cat("\n=================================================================\n")
  cat(sprintf("== %s: Gelman–Rubin with %d chains (all samples) ==\n", label, length(mods)))
  cat(sprintf("Parameters: %d | Multivariate PSRF: %.4f\n\n", nrow(tab), mpsrf))
  print(tab, row.names = FALSE)

  invisible(list(table = tab, mpsrf = mpsrf, mlist = mlist))
}

## ========= Run for the 6 chains =========
mods6 <- list(Mod30, Mod31, Mod32, Mod33, Mod34, Mod35)

## ----- Convergence (Gelman–Rubin) -----
res_b_w <- analyze_component_all(mods6, "Mbw", "φ_w (Mbw)")
res_b_m <- analyze_component_all(mods6, "Mbm", "φ_m (Mbm)")
res_v_w <- analyze_component_all(mods6, "Mvw", "σ²_w (Mvw)")
res_v_m <- analyze_component_all(mods6, "Mvm", "σ²_m (Mvm)")

## ---------- Pooled summaries (generic) ----------
summarize_pooled <- function(mods, comp, label,
                             param_names = NULL,
                             assume_3_deltas = TRUE) {
  has_comp <- vapply(mods, function(m) comp %in% names(m), TRUE)
  if (!all(has_comp)) {
    cat(sprintf("[Info] '%s' missing in some chain; %s not summarized.\n", comp, label))
    return(invisible(NULL))
  }
  mats <- lapply(mods, function(m) as.matrix(m[[comp]]))
  mats <- align_columns_list(mats, what = comp)
  mats <- equalize_lengths(mats)
  all_mat <- do.call(rbind, mats)

  k <- ncol(all_mat)
  cn <- colnames(all_mat)
  if (is.null(cn)) cn <- paste0(label, ".", seq_len(k))

  ## Common case: 3 deltas (k == 3) → one row with Δ1..Δ3
  if (assume_3_deltas && k == 3) {
    meds  <- apply(all_mat, 2, median)
    q025s <- apply(all_mat, 2, quantile, probs = 0.025)
    q975s <- apply(all_mat, 2, quantile, probs = 0.975)
    tab <- data.frame(
      Parameter   = label,
      `Δ1_median` = meds[1], `Δ1_2.5%` = q025s[1], `Δ1_97.5%` = q975s[1],
      `Δ2_median` = meds[2], `Δ2_2.5%` = q025s[2], `Δ2_97.5%` = q975s[2],
      `Δ3_median` = meds[3], `Δ3_2.5%` = q025s[3], `Δ3_97.5%` = q975s[3],
      check.names = FALSE
    )
    tab[, -1] <- round(tab[, -1], 3)
    print(tab, row.names = FALSE)
    return(invisible(tab))
  }

  ## MBeta/MPsi-like case: multiple parameters × 3 deltas (k multiple of 3)
  if (assume_3_deltas && k %% 3 == 0) {
    npar <- k / 3
    med_mat  <- t(matrix(apply(all_mat,  2, median),               3, npar))
    q025_mat <- t(matrix(apply(all_mat,  2, quantile, probs=0.025),3, npar))
    q975_mat <- t(matrix(apply(all_mat,  2, quantile, probs=0.975),3, npar))

    if (is.null(param_names)) param_names <- paste0(label, "[", seq_len(npar), "]")

    tab <- data.frame(
      Parameters   = param_names,
      `Δ1_median`  = med_mat[,1], `Δ1_2.5%` = q025_mat[,1], `Δ1_97.5%` = q975_mat[,1],
      `Δ2_median`  = med_mat[,2], `Δ2_2.5%` = q025_mat[,2], `Δ2_97.5%` = q975_mat[,2],
      `Δ3_median`  = med_mat[,3], `Δ3_2.5%` = q025_mat[,3], `Δ3_97.5%` = q975_mat[,3],
      check.names = FALSE
    )
    tab[, -1] <- round(tab[, -1], 3)
    print(tab, row.names = FALSE)
    return(invisible(tab))
  }

  ## General fallback: long table, one row per column
  med  <- apply(all_mat,  2, median)
  q025 <- apply(all_mat,  2, quantile, probs = 0.025)
  q975 <- apply(all_mat,  2, quantile, probs = 0.975)
  tab <- data.frame(Parameter = cn, Median = med, `2.5%` = q025, `97.5%` = q975, check.names = FALSE)
  tab[, -1] <- round(tab[, -1], 3)
  print(tab, row.names = FALSE)
  invisible(tab)
}

## ================= Table 1 (pooled across all CODA chains) =================
## Assumes you already have `mods6 <- list(Mod30, Mod31, Mod32, Mod33, Mod34, Mod35)`
## and the helper functions `align_columns_list()` and `equalize_lengths()` in memory.

## 1) Collect MBeta matrices from all chains and ensure consistent shape/order
mats_beta <- lapply(mods6, function(m) as.matrix(m$MBeta))   # one MBeta per chain
mats_beta <- align_columns_list(mats_beta, what = "MBeta")   # align columns across chains
mats_beta <- equalize_lengths(mats_beta)                     # truncate to common length

## 2) Pool samples by stacking rows (treating chains as independent draws)
MBeta_all <- do.call(rbind, mats_beta)

## 3) Compute posterior summaries using all pooled samples
med_mat  <- t(matrix(apply(MBeta_all,  2, median),              3, 9))
q025_mat <- t(matrix(apply(MBeta_all,  2, quantile, probs=0.025),3, 9))
q975_mat <- t(matrix(apply(MBeta_all,  2, quantile, probs=0.975),3, 9))

## 4) Assemble the table
param_nms <- c("intercept","Lon","Lat","East","West",
               "North","Demographic density","Military police","playground")

tab <- data.frame(
  Parameters   = param_nms,
  `Δ1_median`  = med_mat[,1], `Δ1_2.5%` = q025_mat[,1], `Δ1_97.5%` = q975_mat[,1],
  `Δ2_median`  = med_mat[,2], `Δ2_2.5%` = q025_mat[,2], `Δ2_97.5%` = q975_mat[,2],
  `Δ3_median`  = med_mat[,3], `Δ3_2.5%` = q025_mat[,3], `Δ3_97.5%` = q975_mat[,3],
  check.names = FALSE
)

## 5) Round and print
tab[, -1] <- round(tab[, -1], 3)
print(tab)

## ================= Table 2 (pooled across all CODA chains) =================
## Assumes you already have `mods6 <- list(Mod30, Mod31, Mod32, Mod33, Mod34, Mod35)`
## and the helper functions `align_columns_list()` and `equalize_lengths()` in memory.

## 1) Collect MPsi matrices from all chains and ensure consistent shape/order
mats_psi <- lapply(mods6, function(m) as.matrix(m$MPsi))     # one MPsi per chain
mats_psi <- align_columns_list(mats_psi, what = "MPsi")      # align columns across chains
mats_psi <- equalize_lengths(mats_psi)                       # truncate to common length

## 2) Pool samples by stacking rows (treating chains as independent draws)
MPsi_all <- do.call(rbind, mats_psi)

## 3) Compute posterior summaries using all pooled samples
med_mat  <- t(matrix(apply(MPsi_all,  2, median),               3, 9))
q025_mat <- t(matrix(apply(MPsi_all,  2, quantile, probs=0.025),3, 9))
q975_mat <- t(matrix(apply(MPsi_all,  2, quantile, probs=0.975),3, 9))

## 4) Assemble the table
param_nms <- c("intercept","Lon","Lat","East","West",
               "North","Demographic density","Military police","playground")

tab <- data.frame(
  Parameters   = param_nms,
  `Δ1_median`  = med_mat[,1], `Δ1_2.5%` = q025_mat[,1], `Δ1_97.5%` = q975_mat[,1],
  `Δ2_median`  = med_mat[,2], `Δ2_2.5%` = q025_mat[,2], `Δ2_97.5%` = q975_mat[,2],
  `Δ3_median`  = med_mat[,3], `Δ3_2.5%` = q025_mat[,3], `Δ3_97.5%` = q975_mat[,3],
  check.names = FALSE
)

## 5) Round and print
tab[, -1] <- round(tab[, -1], 3)
print(tab)

## ----- Summarized tables (all chains pooled) -----
## If each object has 3 columns (Δ1, Δ2, Δ3), one row per parameter:
tab_Mbw <- summarize_pooled(mods6, "Mbw", "φ_w (Mbw)", assume_3_deltas = TRUE)
tab_Mbm <- summarize_pooled(mods6, "Mbm", "φ_m (Mbm)", assume_3_deltas = TRUE)
tab_Mvw <- summarize_pooled(mods6, "Mvw", "σ²_w (Mvw)", assume_3_deltas = TRUE)
tab_Mvm <- summarize_pooled(mods6, "Mvm", "σ²_m (Mvm)", assume_3_deltas = TRUE)

# Figure 2

## ================= Build the figure using SIX chains with per-chain thinning =================
## English (R-style) comments; your plotting block remains EXACTLY the same.
## Target: from each model, keep 3,000 draws spaced by 100 (stride = 100) -> 18,000 total after pooling.

## ---------- Helpers (same behavior you used for CODA alignment) ----------
align_columns_list <- function(mats, what = "MBeta") {
  stopifnot(length(mats) >= 2)
  cn1 <- colnames(mats[[1]])
  if (is.null(cn1)) {
    k <- ncol(mats[[1]])
    ok <- vapply(mats, function(M) ncol(M) == k, TRUE)
    if (!all(ok)) stop(sprintf("%s: different number of columns across chains.", what))
    return(mats)
  }
  for (i in seq_along(mats)) {
    cni <- colnames(mats[[i]])
    if (is.null(cni) || !setequal(cni, cn1))
      stop(sprintf("%s: different column sets across chains.", what))
    mats[[i]] <- mats[[i]][, cn1, drop = FALSE]
  }
  mats
}

equalize_lengths <- function(mats) {
  min_n <- min(vapply(mats, nrow, 1L))
  lapply(mats, function(M) M[(nrow(M) - min_n + 1):nrow(M), , drop = FALSE])
}

## ---------- List of chains ----------
mods6 <- list(Mod30, Mod31, Mod32, Mod33, Mod34, Mod35)

## ---------- Extract MW/MM, align columns, equalize lengths ----------
mats_MW <- lapply(mods6, function(m) as.matrix(m$MW))  # MW per chain
mats_MM <- lapply(mods6, function(m) as.matrix(m$MM))  # MM per chain

mats_MW <- align_columns_list(mats_MW, what = "MW")
mats_MM <- align_columns_list(mats_MM, what = "MM")

mats_MW <- equalize_lengths(mats_MW)
mats_MM <- equalize_lengths(mats_MM)

## ---------- Per-chain thinning: 3,000 draws per chain with stride = 100 ----------
n_keep  <- 3000L   # draws per chain to retain
stride  <- 100L    # step between retained draws

thin_chain <- function(M, n_keep, stride) {
  N <- nrow(M)
  start_idx <- max(1L, N - (n_keep - 1L) * stride)   # back-fill from the end
  idx <- seq.int(from = start_idx, to = N, by = stride)
  ## Ensure exactly n_keep rows (trim or pad if off by one due to bounds)
  if (length(idx) > n_keep) idx <- tail(idx, n_keep)
  if (length(idx) < n_keep) idx <- c(rep(idx[1], n_keep - length(idx)), idx)  # very unlikely
  M[idx, , drop = FALSE]
}

mats_MW <- lapply(mats_MW, thin_chain, n_keep = n_keep, stride = stride)
mats_MM <- lapply(mats_MM, thin_chain, n_keep = n_keep, stride = stride)

## ---------- Pool thinned chains (final size: 6 * 3,000 = 18,000 draws) ----------
MW <- do.call(rbind, mats_MW)
MM <- do.call(rbind, mats_MM)

## ---------- Define time breaks and site/segment dimensions ----------
ini<-c(0,breaks[1:(length(breaks)-1)])
fim<-breaks

n_seg <- length(fim)
stopifnot(ncol(MW) %% n_seg == 0L)
n_sites <- ncol(MW) %/% n_seg

## Dummy shape anchors so nrow(Wm)/ncol(Wm) are numeric (prevents extent errors)
Wm <- matrix(NA_real_, nrow = n_seg, ncol = n_sites)
Mm <- matrix(NA_real_, nrow = n_seg, ncol = n_sites)

## ========================== YOUR ORIGINAL CODE (UNCHANGED) ==========================
mtnp <- function(dados) {
  a <- c(0, dados)
  m <- length(dados)
  k <- 1
  corr <- array(NA, dim = c(m, 2))
  for (i in 1:m) {
    corr[i, 1] <- (a[i] + a[i + 1]) / 2
    tt <- (a[i] + a[i + 1]) / 2
    rmt <- (1 / k) * (i - 1) + (tt - a[i]) / (k * (a[i + 1] - a[i]))
    corr[i, 2] <- rmt
  }
  corr
}
#
mf <- function(Wl, Ml, tau) {
  res <- exp(Wl) * tau^(exp(Ml))
  res
}
#
ml <- function(Wl, Ml, tini, tfim) {
  xx <- seq(tini, tfim, length = 100)
  res <- NULL

  for (i in 1:length(xx)) {
    aux <- exp(Wl) * xx[i]^(exp(Ml))
    res <- c(res, aux)
  }

  res1 <- list(res, xx)
  res1
}
#
mdil <- function(Wo, Mo, tauini, taufim) {
  L <- length(Wo)
  vecx <- NULL
  vecy <- NULL
  for (l in 1:L) {
    if (l == 1) {
      aux <- ml(Wo[l], Mo[l], tauini[l], taufim[l])
      vecx <- c(vecx, aux[[2]])
      vecy <- c(vecy, aux[[1]])
    } else {
      acum <- 0
      for (k in 1:(l - 1)) {
        aux1 <- mf(Wo[k], Mo[k], taufim[k]) - mf(Wo[k + 1], Mo[k + 1], taufim[k])
        acum <- acum + aux1
      }
      temp <- ml(Wo[l], Mo[l], tauini[l], taufim[l])
      aux2 <- acum + temp[[1]]
      vecx <- c(vecx, temp[[2]])
      vecy <- c(vecy, aux2)
    }
  }

  res <- list(vecx, vecy)
  res
}
#
# Sites l=1,...,34
Bairro<-c("Lagoa Nova","Nova Descoberta","Candelária","Capim Macio","Pitimbu","Neópolis","Ponta Negra","Santos Reis","Rocas","Ribeira","Praia do Meio","Petrópolis","Areia Preta","Cidade Alta","Mãe Luiza","Alecrim","Barro Vermelho","Tirol","Lagoa Seca","Cidade Esperança","Quintas","Dix_Sept Rosado","Bom Pastor","Nossa Sra. de Nazaré","Felipe Camarão","Cidade Nova","Guarapes","Planalto","Igapó","Potengi","Nossa Sra. da Apresentação","Lagoa Azul","Pajuçara","Redinha")

l <-3

y1 <- NULL
for (h in 1:nrow(MW)) {
  Wme <- matrix(MW[h, ], nrow(Wm), ncol(Wm))
  Mme <- matrix(MM[h, ], nrow(Mm), ncol(Mm))
  lll <- mdil(Wme[, l], Mme[, l], ini, fim)
  y1 <- rbind(y1, t(as.matrix(lll[[2]])))
}

medy1 <- apply(y1, 2, median)
per25 <- apply(y1, 2, quantile, probs = 0.025)
per975 <- apply(y1, 2, quantile, probs = 0.975)
infy <- min(per25)
supy <- max(per975)
retemp <- lll[[1]]
infx <- min(retemp)
supx <- max(retemp)
xx <- c(retemp, rev(retemp))
yy <- c(per25, rev(per975))
plot(xx, yy, type = "n", xlab = " ", ylab = " ", xlim = c(infx, supx), ylim = c(infy, supy), cex.lab = 1.5)
polygon(xx, yy, col = "gray", border = NA)
par(new = T)
plot(mtnp(natal_theft_times[, l]), type = "l", main = " ", xlim = c(infx, supx), ylim = c(infy, supy), xlab = " ", ylab = " ")
par(new = T)
plot(retemp, medy1,
  type = "l", xlim = c(infx, supx), ylim = c(infy, supy), xlab = " ",
  ylab = " ", main = Bairro[l], lwd = 2, cex.main = 2, col = 2
)

abline(v=fim[1],lty=2)
abline(v=fim[2],lty=2)

# Figure 3

## ================= alpha = exp(M) using ALL iterations from Mod30..Mod35 =================

## ---- Canonical neighborhood names (model order) ----
bairro_names <- c(
  "Lagoa Nova","Nova Descoberta","Candelária","Capim Macio","Pitimbu","Neópolis","Ponta Negra",
  "Santos Reis","Rocas","Ribeira","Praia do Meio","Petrópolis","Areia Preta","Cidade Alta",
  "Mãe Luiza","Alecrim","Barro Vermelho","Tirol","Lagoa Seca",
  "Cidade da Esperança","Quintas","Dix-Sept Rosado","Bom Pastor","Nossa Senhora de Nazaré",
  "Felipe Camarão","Cidade Nova","Guarapes","Planalto",
  "Igapó","Potengi","Nossa Senhora da Apresentação","Lagoa Azul","Pajuçara","Redinha"
)
stopifnot(length(bairro_names) == 34)

## ---- List your chains here (in any order you prefer) ----
mods <- list(Mod30, Mod31, Mod32, Mod33, Mod34, Mod35)

## ---- Dimensions: L deltas x N neighborhoods ----
L <- 3
N <- length(bairro_names)

## ---- Pull MM matrices and run basic checks ----
MM_list <- lapply(seq_along(mods), function(i) {
  x <- mods[[i]]$MM
  if (!is.matrix(x)) stop(sprintf("Mod%d$MM is not a matrix.", 29+i))
  if (ncol(x) < L * N) stop(sprintf("Mod%d$MM has fewer than %d columns.", 29+i, L*N))
  x
})

## ---- Helper: summarize exp(M) for a given column index across all chains ----
summ_expM_col <- function(MM_list, col_idx) {
  v <- unlist(lapply(MM_list, function(MM) MM[, col_idx]), use.names = FALSE)
  a <- exp(v)  # alpha = exp(M)
  qs <- quantile(a, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
  c(mean = mean(a, na.rm = TRUE), q2.5 = qs[1], q97.5 = qs[2])
}

## ---- Loop over neighborhoods (j) and deltas (k) using the correct column mapping ----
# Column mapping consistent with: MMiter <- matrix(MM[i, ], L, N)  (column-major fill)
# -> For neighborhood j in {1..N} and delta k in {1..L}, column = (j-1)*L + k
S_mean  <- matrix(NA_real_, nrow = N, ncol = L)
S_q025  <- matrix(NA_real_, nrow = N, ncol = L)
S_q975  <- matrix(NA_real_, nrow = N, ncol = L)

for (j in 1:N) {
  for (k in 1:L) {
    col_idx <- (j - 1) * L + k
    s <- summ_expM_col(MM_list, col_idx)
    S_mean[j, k] <- s[["mean"]]
    S_q025[j, k] <- s[["q2.5"]]
    S_q975[j, k] <- s[["q97.5"]]
  }
}

## ---- Build the wide table exactly in the requested format ----
alpha_table <- data.frame(
  Neighborhood = bairro_names,
  `Δ1_mean`   = S_mean[, 1], `Δ1_2.5%`  = S_q025[, 1], `Δ1_97.5%` = S_q975[, 1],
  `Δ2_mean`   = S_mean[, 2], `Δ2_2.5%`  = S_q025[, 2], `Δ2_97.5%` = S_q975[, 2],
  `Δ3_mean`   = S_mean[, 3], `Δ3_2.5%`  = S_q025[, 3], `Δ3_97.5%` = S_q975[, 3],
  check.names = FALSE
)

## ---- Nicely formatted version for console printing (3 decimals) ----
alpha_table_fmt <- alpha_table
alpha_table_fmt[-1] <- lapply(alpha_table_fmt[-1], function(v) sprintf("%.3f", v))
print(alpha_table_fmt, row.names = FALSE)

## ---- (Optional) save numeric CSV (not formatted) ----
# write.csv(alpha_table, "alpha_posterior_summary_all_chains_correct.csv", row.names = FALSE)

## ---- (Optional) long-format table handy for plotting ----
alpha_long <- rbind(
  data.frame(name = bairro_names, Interval = "Δ1", mean = S_mean[,1], q2.5 = S_q025[,1], q97.5 = S_q975[,1]),
  data.frame(name = bairro_names, Interval = "Δ2", mean = S_mean[,2], q2.5 = S_q025[,2], q97.5 = S_q975[,2]),
  data.frame(name = bairro_names, Interval = "Δ3", mean = S_mean[,3], q2.5 = S_q025[,3], q97.5 = S_q975[,3])
)
row.names(alpha_long) <- NULL

## ================== Maps of alpha = exp(M): Δ1, Δ2, Δ3 (shared scale, robust) ==================
pkgs <- c("sf","dplyr","ggplot2","osmdata","stringi","viridis","cowplot")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(sf); library(dplyr); library(ggplot2); library(osmdata); library(stringi)
library(viridis); library(cowplot)

## --- Checks: need alpha_long with columns name, Interval {"Δ1","Δ2","Δ3"}, mean
stopifnot(exists("alpha_long"))
stopifnot(all(c("name","Interval","mean") %in% names(alpha_long)))

## --- Canonical neighborhood names (model order) ---
bairro_names <- c(
  "Lagoa Nova","Nova Descoberta","Candelária","Capim Macio","Pitimbu","Neópolis","Ponta Negra",
  "Santos Reis","Rocas","Ribeira","Praia do Meio","Petrópolis","Areia Preta","Cidade Alta",
  "Mãe Luiza","Alecrim","Barro Vermelho","Tirol","Lagoa Seca",
  "Cidade da Esperança","Quintas","Dix-Sept Rosado","Bom Pastor","Nossa Senhora de Nazaré",
  "Felipe Camarão","Cidade Nova","Guarapes","Planalto",
  "Igapó","Potengi","Nossa Senhora da Apresentação","Lagoa Azul","Pajuçara","Redinha"
)

## --- Name normalizer ---
normalize_str <- function(x){
  x %>% stringi::stri_trans_general("Latin-ASCII") %>% tolower() %>%
    gsub("[^a-z0-9 ]"," ",.) %>% gsub("\\s+"," ",.) %>% trimws()
}

## --- Build polygons as 'natal_sf' if missing (avoid name clash with any 'plot_sf' function) ---
if (!exists("natal_sf")) {
  bb <- osmdata::getbb("Natal, Rio Grande do Norte, Brazil", format_out = "sf_polygon")
  stopifnot(!is.null(bb), nrow(bb) > 0); bbx <- sf::st_bbox(bb)

  q_admin_fun <- function(b) opq(b) |>
    add_osm_feature("boundary","administrative") |>
    add_osm_feature("admin_level", c("9","10"))
  q_place_fun <- function(b) opq(b) |>
    add_osm_feature("place", c("neighbourhood","suburb","quarter"))

  servers <- c("https://overpass-api.de/api/interpreter",
               "https://overpass.kumi.systems/api/interpreter")
  fetch_osm <- function(b, qfun){
    last <- NULL
    for (srv in servers){
      osmdata::set_overpass_url(srv)
      res <- try(osmdata::osmdata_sf(qfun(b)), silent = TRUE)
      if (!inherits(res, "try-error")) return(res)
      last <- res
    }
    stop("OSM download failed on all servers. Last error: ", as.character(last))
  }

  od_admin <- fetch_osm(bbx, q_admin_fun)
  od_place <- fetch_osm(bbx, q_place_fun)

  standardize_name_geom <- function(x){
    if (is.null(x) || nrow(x)==0) return(NULL)
    if (!inherits(x,"sf")) x <- st_as_sf(x)
    cand <- c("name","name:pt","name:en","official_name","short_name","alt_name")
    hit  <- intersect(cand, names(x))
    nm   <- if ("name" %in% hit) x[["name"]] else if (length(hit)>0) x[[hit[1]]] else rep(NA_character_, nrow(x))
    suppressWarnings(st_make_valid(st_sf(name = nm, geometry = st_geometry(x))))
  }

  poly_list <- list(od_admin$osm_polygons, od_admin$osm_multipolygons,
                    od_place$osm_polygons, od_place$osm_multipolygons) |>
    lapply(standardize_name_geom)
  poly_list <- poly_list[!vapply(poly_list, is.null, logical(1))]
  stopifnot(length(poly_list) > 0)

  osmp <- do.call(rbind, poly_list) |> suppressWarnings(st_make_valid())
  osmp <- st_intersection(osmp, st_union(bb)) |>
    mutate(name = ifelse(is.na(name), "", name),
           name_norm = normalize_str(name))

  bairro_norm <- normalize_str(bairro_names)
  match_idx <- match(osmp$name_norm, bairro_norm)
  unmatched <- which(is.na(match_idx) & nzchar(osmp$name_norm))
  if (length(unmatched) > 0){
    for (i in unmatched){
      cand <- agrep(osmp$name_norm[i], bairro_norm, max.distance = 0.25)
      if (length(cand) == 1) match_idx[i] <- cand
    }
  }
  matched <- !is.na(match_idx); stopifnot(any(matched))

  natal_sf <- osmp[matched, , drop = FALSE] |>
    mutate(name = bairro_names[match_idx[matched]]) |>
    select(name, geometry) |>
    st_make_valid() |>
    mutate(.area = as.numeric(st_area(geometry))) |>
    group_by(name) |>
    slice_max(.area, n = 1, with_ties = FALSE) |>
    ungroup() |>
    select(-.area)
}

## --- Extra safety: ensure natal_sf is an sf object, not a function ---
if (is.function(natal_sf)) stop("Object 'natal_sf' is a function. Choose another name for polygons.")

## --- Prepare alpha_long and join ---
alpha_long <- alpha_long |>
  mutate(Interval = factor(Interval, levels = c("Δ1","Δ2","Δ3")))

df_all <- dplyr::left_join(
  natal_sf,
  alpha_long |> dplyr::select(name, Interval, mean),
  by = c("name" = "name")
)

## --- Shared scale across Δ1, Δ2, Δ3 ---
vmin <- min(df_all$mean, na.rm = TRUE)
vmax <- max(df_all$mean, na.rm = TRUE)

## --- Map builder (no titles; legend title removed) ---
make_alpha_map <- function(interval_label){
  ggplot(df_all %>% dplyr::filter(Interval == interval_label)) +
    geom_sf(aes(fill = mean), color = "white", linewidth = 0.25, na.rm = TRUE) +
    scale_fill_viridis(
      name   = NULL,                # remove legend title
      option = "D",
      limits = c(vmin, vmax),       # SAME scale for all three maps
      oob    = scales::squish,
      na.value = "grey85"
    ) +
    coord_sf(datum = NA) +
    theme_void(base_size = 12) +
    theme(
      legend.position = "right",
      plot.title    = element_blank(),
      plot.subtitle = element_blank(),
      plot.caption  = element_blank()
    )
}

## --- Three maps (Δ1, Δ2, Δ3) ---
map_d1 <- make_alpha_map("Δ1")
map_d2 <- make_alpha_map("Δ2")
map_d3 <- make_alpha_map("Δ3")

## Show side-by-side (optional)
#cowplot::plot_grid(map_d1, map_d2, map_d3, ncol = 3, align = "hv")

## Save (optional)
# ggsave("alpha_delta1.png", map_d1, width=6, height=7, dpi=300)
# ggsave("alpha_delta2.png", map_d2, width=6, height=7, dpi=300)
# ggsave("alpha_delta3.png", map_d3, width=6, height=7, dpi=300)
# ggsave("alpha_3maps_shared_scale.png",
#        cowplot::plot_grid(map_d1, map_d2, map_d3, ncol=3, align="hv"),
#        width=18, height=7, dpi=300)

## ================= Map of alpha (posterior mean) — Δ1 only (legend title removed) =================
## Assumptions:
## - alpha_long: columns {name, Interval in {"Δ1","Δ2","Δ3"}, mean}
## - plot_sf: sf polygons for the 34 neighborhoods (built previously)
## - We compute vmin/vmax once to reuse the SAME scale across Δ1, Δ2, Δ3

library(dplyr)
library(ggplot2)
library(sf)
library(viridis)
library(cowplot)

stopifnot(exists("alpha_long"), exists("plot_sf"))
stopifnot(all(c("name","Interval","mean") %in% names(alpha_long)))

# Ensure factor levels order
alpha_long <- alpha_long %>%
  mutate(Interval = factor(Interval, levels = c("Δ1","Δ2","Δ3")))

# Join polygons with alpha summaries
df_all <- plot_sf %>%
  left_join(alpha_long %>% select(name, Interval, mean), by = c("name" = "name"))

# Global scale (shared across the three maps)
vmin <- min(df_all$mean, na.rm = TRUE)
vmax <- max(df_all$mean, na.rm = TRUE)

# Δ1 subset
df_d1 <- df_all %>% filter(Interval == "Δ1")

# Build Δ1 map: no title, legend on the right, bottom-centered "a)"
p_d1 <- ggplot(df_d1) +
  geom_sf(aes(fill = mean), color = "white", linewidth = 0.25, na.rm = TRUE) +
  scale_fill_viridis(
    name   = NULL,           # <- remove legend title
    option = "D",
    limits = c(vmin, vmax),  # <- SAME scale across Δ1, Δ2, Δ3
    oob    = scales::squish,
    na.value = "grey85"
  ) +
  coord_sf(datum = NA) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title    = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption  = element_blank(),
    plot.margin   = margin(t = 5.5, r = 5.5, b = 24, l = 5.5)  # space for "a)"
  )

# Add bottom-centered panel tag "a)"
p_d1_tag <- cowplot::ggdraw(p_d1) +
  cowplot::draw_label("a)", x = 0.50, y = 0.04, hjust = 0.5, vjust = 0,
                      fontface = "bold", size = 14)

print(p_d1_tag)

# To produce Δ2 / Δ3 with the SAME scale:
# df_d2 <- df_all %>% filter(Interval == "Δ2")
# df_d3 <- df_all %>% filter(Interval == "Δ3")
# ...repeat the same ggplot, keeping scale_fill_viridis(limits = c(vmin, vmax), name = NULL)


## ================= Map of alpha (posterior mean) — Δ2 only (shared scale, no legend title) =================
## Requirements:
## - alpha_long: data.frame with columns {name, Interval in {"Δ1","Δ2","Δ3"}, mean}
## - plot_sf: sf polygons for the 34 neighborhoods (already built)
## - We compute vmin/vmax from all intervals to keep the SAME scale across Δ1/Δ2/Δ3


stopifnot(exists("alpha_long"), exists("plot_sf"))
stopifnot(all(c("name","Interval","mean") %in% names(alpha_long)))

# Ensure factor level order
alpha_long <- alpha_long %>%
  mutate(Interval = factor(Interval, levels = c("Δ1","Δ2","Δ3")))

# Join polygons with alpha summaries (all intervals for global scale)
df_all <- plot_sf %>%
  left_join(alpha_long %>% select(name, Interval, mean), by = c("name" = "name"))

# Global continuous scale shared across Δ1, Δ2, Δ3
vmin <- min(df_all$mean, na.rm = TRUE)
vmax <- max(df_all$mean, na.rm = TRUE)

# Subset for Δ2
df_d2 <- df_all %>% filter(Interval == "Δ2")

# Build Δ2 map: no title, legend on the right, bottom-centered "b)"
p_d2 <- ggplot(df_d2) +
  geom_sf(aes(fill = mean), color = "white", linewidth = 0.25, na.rm = TRUE) +
  scale_fill_viridis(
    name   = NULL,            # remove legend title
    option = "D",
    limits = c(vmin, vmax),   # SAME scale across Δ1/Δ2/Δ3
    oob    = scales::squish,
    na.value = "grey85"
  ) +
  coord_sf(datum = NA) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title    = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption  = element_blank(),
    plot.margin   = margin(t = 5.5, r = 5.5, b = 24, l = 5.5)  # space for "b)"
  )

# Add bottom-centered panel tag "b)"
p_d2_tag <- cowplot::ggdraw(p_d2) +
  cowplot::draw_label("b)", x = 0.50, y = 0.04, hjust = 0.5, vjust = 0,
                      fontface = "bold", size = 14)

print(p_d2_tag)

# Optional: save to file
# ggsave("alpha_map_delta2.png", p_d2_tag, width = 6, height = 7, dpi = 300)

## ================= Map of alpha (posterior mean) — Δ3 only (shared scale, no legend title) =================
## Requirements:
## - alpha_long: data.frame with columns {name, Interval in {"Δ1","Δ2","Δ3"}, mean}
## - plot_sf: sf polygons for the 34 neighborhoods (already built)
## - We compute vmin/vmax from all intervals to keep the SAME scale across Δ1/Δ2/Δ3


stopifnot(exists("alpha_long"), exists("plot_sf"))
stopifnot(all(c("name","Interval","mean") %in% names(alpha_long)))

# Ensure factor level order
alpha_long <- alpha_long %>%
  mutate(Interval = factor(Interval, levels = c("Δ1","Δ2","Δ3")))

# Join polygons with alpha summaries (all intervals for global scale)
df_all <- plot_sf %>%
  left_join(alpha_long %>% select(name, Interval, mean), by = c("name" = "name"))

# Global continuous scale shared across Δ1, Δ2, Δ3
vmin <- min(df_all$mean, na.rm = TRUE)
vmax <- max(df_all$mean, na.rm = TRUE)

# Subset for Δ3
df_d3 <- df_all %>% filter(Interval == "Δ3")

# Build Δ3 map: no legend title, legend on the right, bottom-centered "c)"
p_d3 <- ggplot(df_d3) +
  geom_sf(aes(fill = mean), color = "white", linewidth = 0.25, na.rm = TRUE) +
  scale_fill_viridis(
    name   = NULL,            # remove legend title
    option = "D",
    limits = c(vmin, vmax),   # SAME scale across Δ1/Δ2/Δ3
    oob    = scales::squish,
    na.value = "grey85"
  ) +
  coord_sf(datum = NA) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title    = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption  = element_blank(),
    plot.margin   = margin(t = 5.5, r = 5.5, b = 24, l = 5.5)  # space for "c)"
  )

# Add bottom-centered panel tag "c)"
p_d3_tag <- cowplot::ggdraw(p_d3) +
  cowplot::draw_label("c)", x = 0.50, y = 0.04, hjust = 0.5, vjust = 0,
                      fontface = "bold", size = 14)

print(p_d3_tag)

# ============================ Figure 4(a): g_{2,1} (8-dec rounding) ============================
# Per-iteration reshaping:
#   MMiter <- matrix(MM[i, ], L, N); MWiter <- matrix(MW[i, ], L, N)
#   alpha = exp(t(MMiter))  # N x L
#   gamma = exp(t(MWiter))  # N x L
# Inputs expected: MW, MM (iters x N*L), fim = c(tau2_start, tau2_end)
# Output: tab_g21 (numeric) and tab_g21_fmt (strings)

# ---- Canonical neighborhood names (N = 34) ----
bairro_names <- c(
  "Lagoa Nova","Nova Descoberta","Candelária","Capim Macio","Pitimbu","Neópolis","Ponta Negra",
  "Santos Reis","Rocas","Ribeira","Praia do Meio","Petrópolis","Areia Preta","Cidade Alta",
  "Mãe Luiza","Alecrim","Barro Vermelho","Tirol","Lagoa Seca",
  "Cidade da Esperança","Quintas","Dix-Sept Rosado","Bom Pastor","Nossa Senhora de Nazaré",
  "Felipe Camarão","Cidade Nova","Guarapes","Planalto",
  "Igapó","Potengi","Nossa Senhora da Apresentação","Lagoa Azul","Pajuçara","Redinha"
)
stopifnot(length(bairro_names) == 34)

# ---- Basic dims ----
L <- 3
N <- length(bairro_names)
stopifnot(is.matrix(MM), is.matrix(MW))
stopifnot(ncol(MM) >= L*N, ncol(MW) >= L*N, nrow(MM) == nrow(MW))

# ---- Time points for g_{2,1}: 5% into Δ2 (start) and upper bound of Δ2 (end) ----
stopifnot(is.numeric(fim), length(fim) == 2, isTRUE(fim[2] > fim[1]))
Rdelta <- 0.05 * (fim[2] - fim[1])
TAini  <- fim[1] + Rdelta
TAfim  <- fim[2]

# ---- Accumulate g2ini and g2fim across iterations (rows = iters, cols = neighborhoods) ----
niter <- nrow(MW)
g2ini <- matrix(NA_real_, nrow = niter, ncol = N)
g2fim <- matrix(NA_real_, nrow = niter, ncol = N)

for (i in 1:niter) {
  MMiter <- matrix(MM[i, 1:(L*N)], nrow = L, ncol = N)  # column-major fill
  MWiter <- matrix(MW[i, 1:(L*N)], nrow = L, ncol = N)

  alp <- exp(t(MMiter))  # N x L (alpha)
  gam <- exp(t(MWiter))  # N x L (gamma)

  d12 <- alp[, 1] - alp[, 2]   # Δ1 - Δ2 per neighborhood
  g2ini[i, ] <- TAini^d12 - (gam[, 2] / gam[, 1])
  g2fim[i, ] <- TAfim^d12 - (gam[, 2] / gam[, 1])
}

# ---- 95% credible intervals (raw) per neighborhood ----
g2ini_ci_raw <- t(apply(g2ini, 2, quantile, c(0.025, 0.975), na.rm = TRUE))
g2fim_ci_raw <- t(apply(g2fim,  2, quantile, c(0.025, 0.975), na.rm = TRUE))

# ---- Round to EIGHT decimals BEFORE classification (closed intervals) ----
g2ini_ci <- round(g2ini_ci_raw, 8)
g2fim_ci <- round(g2fim_ci_raw,  8)
colnames(g2ini_ci) <- c("q2.5", "q97.5")
colnames(g2fim_ci) <- c("q2.5", "q97.5")

# ---- Closed-interval position relative to zero (inclusive at 0) ----
side_of_zero_closed <- function(lower, upper) {
  if (is.na(lower) || is.na(upper)) return(NA_character_)
  if (upper < 0) return("left")
  if (lower > 0) return("right")
  "cross"
}
s_ini <- apply(g2ini_ci, 1, function(z) side_of_zero_closed(z[1], z[2]))
s_fim <- apply(g2fim_ci,  1, function(z) side_of_zero_closed(z[1], z[2]))

# ---- Classification rules on the 8-dec rounded CLOSED intervals ----
symbol <- character(N)
for (k in seq_len(N)) {
  if (s_ini[k] == "left"  && s_fim[k] == "left") {
    symbol[k] <- "++"
  } else if ((s_ini[k] %in% c("right","cross")) && s_fim[k] == "left") {
    symbol[k] <- "+"
  } else if (s_ini[k] == "cross" && s_fim[k] == "cross") {
    symbol[k] <- "="
  } else if ((s_ini[k] %in% c("left","cross")) && s_fim[k] == "right") {
    symbol[k] <- "-"
  } else if (s_ini[k] == "right" && s_fim[k] == "right") {
    symbol[k] <- "--"
  } else {
    symbol[k] <- "="
  }
}

# ---- Final tables ----
tab_g21 <- data.frame(
  Neighborhood = bairro_names,
  g2ini_q2.5   = g2ini_ci[, "q2.5"],
  g2ini_q97.5  = g2ini_ci[, "q97.5"],
  g2fim_q2.5   = g2fim_ci[, "q2.5"],
  g2fim_q97.5  = g2fim_ci[, "q97.5"],
  symbol       = symbol,
  check.names  = FALSE
)

# Pretty print (fixed eight decimals)
tab_g21_fmt <- transform(
  tab_g21,
  g2ini_q2.5  = sprintf("%.8f", g2ini_q2.5),
  g2ini_q97.5 = sprintf("%.8f", g2ini_q97.5),
  g2fim_q2.5  = sprintf("%.8f", g2fim_q2.5),
  g2fim_q97.5 = sprintf("%.8f", g2fim_q97.5)
)
print(tab_g21_fmt, row.names = FALSE)


## =================== Build plot_sf (OSM) and map g_{2,1} classes ===================

# ---- Packages ----
pkgs <- c("sf","dplyr","ggplot2","osmdata","stringi","cowplot")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(sf); library(dplyr); library(ggplot2); library(osmdata); library(stringi); library(cowplot)

# ---- Require tab_g21 from the previous step ----
stopifnot(exists("tab_g21"), all(c("Neighborhood","symbol") %in% names(tab_g21)))

# ---- Canonical neighborhood names (N = 34) ----
bairro_names <- c(
  "Lagoa Nova","Nova Descoberta","Candelária","Capim Macio","Pitimbu","Neópolis","Ponta Negra",
  "Santos Reis","Rocas","Ribeira","Praia do Meio","Petrópolis","Areia Preta","Cidade Alta",
  "Mãe Luiza","Alecrim","Barro Vermelho","Tirol","Lagoa Seca",
  "Cidade da Esperança","Quintas","Dix-Sept Rosado","Bom Pastor","Nossa Senhora de Nazaré",
  "Felipe Camarão","Cidade Nova","Guarapes","Planalto",
  "Igapó","Potengi","Nossa Senhora da Apresentação","Lagoa Azul","Pajuçara","Redinha"
)
stopifnot(length(bairro_names) == 34)

# ---- String normalizer (ASCII, lower, strip punctuation) ----
normalize_str <- function(x) {
  x %>%
    stringi::stri_trans_general("Latin-ASCII") %>%
    tolower() %>%
    gsub("[^a-z0-9 ]", " ", .) %>%
    gsub("\\s+", " ", .) %>%
    trimws()
}

# ---- Get Natal bounding polygon and bbox ----
bb <- osmdata::getbb("Natal, Rio Grande do Norte, Brazil", format_out = "sf_polygon")
if (is.null(bb) || nrow(bb) == 0) stop("Could not resolve bounding box for Natal/RN via OSM.")
bbx <- sf::st_bbox(bb)

# ---- Rotate Overpass servers until one works ----
servers <- c(
  "https://overpass-api.de/api/interpreter",
  "https://overpass.kumi.systems/api/interpreter"
)

fetch_osm <- function(bbox, query_fun) {
  last_err <- NULL
  for (srv in servers) {
    osmdata::set_overpass_url(srv)
    res <- try({
      q <- query_fun(bbox)
      osmdata::osmdata_sf(q)
    }, silent = TRUE)
    if (!inherits(res, "try-error")) return(res)
    last_err <- res
  }
  stop("All Overpass servers failed. Last error: ", as.character(last_err))
}

q_admin_fun <- function(bbox) {
  opq(bbox = bbox) |>
    add_osm_feature(key = "boundary", value = "administrative") |>
    add_osm_feature(key = "admin_level", value = c("9","10"))
}

q_place_fun <- function(bbox) {
  opq(bbox = bbox) |>
    add_osm_feature(key = "place", value = c("neighbourhood","suburb","quarter"))
}

od_admin <- fetch_osm(bbx, q_admin_fun)
od_place <- fetch_osm(bbx, q_place_fun)

# ---- Keep ONLY name + geometry to avoid rbind mismatches ----
standardize_name_geom <- function(x) {
  if (is.null(x) || nrow(x) == 0) return(NULL)
  if (!inherits(x, "sf")) x <- st_as_sf(x)
  candidates <- c("name","name:pt","name:en","official_name","short_name","alt_name")
  hit <- intersect(candidates, names(x))
  nm <- if ("name" %in% hit) x[["name"]] else if (length(hit) > 0) x[[hit[1]]] else rep(NA_character_, nrow(x))
  out <- st_sf(name = nm, geometry = st_geometry(x))
  suppressWarnings(st_make_valid(out))
}

poly_list <- list(
  od_admin$osm_polygons, od_admin$osm_multipolygons,
  od_place$osm_polygons, od_place$osm_multipolygons
) |> lapply(standardize_name_geom)

poly_list <- poly_list[!vapply(poly_list, is.null, logical(1))]
if (length(poly_list) == 0) stop("No neighborhood-like polygons found on OSM for Natal/RN.")

osmp <- do.call(rbind, poly_list)
osmp <- suppressWarnings(st_make_valid(osmp))
bb_union <- sf::st_union(bb)
osmp <- sf::st_intersection(osmp, bb_union)

# ---- Normalize names and match to canonical list (exact then fuzzy) ----
osmp <- osmp |>
  mutate(name = ifelse(is.na(name), "", name),
         name_norm = normalize_str(name))

bairro_norm <- normalize_str(bairro_names)

# exact match
match_idx <- match(osmp$name_norm, bairro_norm)

# fuzzy match where still NA
unmatched <- which(is.na(match_idx) & nzchar(osmp$name_norm))
if (length(unmatched) > 0) {
  for (i in unmatched) {
    cand <- agrep(osmp$name_norm[i], bairro_norm, max.distance = 0.25)
    if (length(cand) == 1) match_idx[i] <- cand
  }
}

matched <- !is.na(match_idx)
if (!any(matched)) stop("Could not match any OSM polygon names to the 34-neighborhood list.")

natal_neighborhoods <- osmp[matched, , drop = FALSE] |>
  mutate(name = bairro_names[match_idx[matched]]) |>
  select(name, geometry) |>
  st_make_valid()

# If multiple polygons per neighborhood, keep the largest
natal_neighborhoods <- natal_neighborhoods |>
  mutate(.area = as.numeric(st_area(geometry))) |>
  group_by(name) |>
  slice_max(.area, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(-.area)

# ---- Build plot_sf and join with tab_g21 ----
plot_sf <- natal_neighborhoods

# Normalize symbols: treat "==" as "=" just in case
tab_g21 <- tab_g21 |>
  mutate(symbol = ifelse(symbol == "==", "=", symbol))

plot_g21 <- plot_sf |>
  left_join(tab_g21 |> select(Neighborhood, symbol),
            by = c("name" = "Neighborhood"))

# ---- Colors: ++ dark red, + light red, = very light purple (near white)
pal <- c(
  "++" = "#A50F15",  # dark red
  "+"  = "#FB6A4A",  # light red
  "="  = "#F5F0FF"   # very light purple
)

# Ensure factor levels (only the three of interest will show)
plot_g21 <- plot_g21 |>
  mutate(symbol = factor(symbol, levels = c("++","+","=")))

# ---- Draw categorical map (no title/subtitle/caption) ----
p_g21_map <- ggplot(plot_g21) +
  geom_sf(aes(fill = symbol), color = "white", linewidth = 0.25, na.rm = TRUE) +
  scale_fill_manual(
    name = expression(g[2,1]~"class"),
    values = pal,
    drop = TRUE,
    na.value = "grey85",
    guide = guide_legend(override.aes = list(color = NA))
  ) +
  coord_sf(datum = NA) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title    = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption  = element_blank(),
    plot.margin   = margin(t = 5.5, r = 5.5, b = 24, l = 5.5) # extra bottom space for "a)"
  )

# ---- Add bottom-centered panel tag "a)" ----
p_g21_map_tag <- cowplot::ggdraw(p_g21_map) +
  cowplot::draw_label("a)", x = 0.50, y = 0.04, hjust = 0.5, vjust = 0,
                      fontface = "bold", size = 14)

print(p_g21_map_tag)

# (optional) save
# ggsave("g21_map_panel_a.png", p_g21_map_tag, width = 6, height = 7, dpi = 300)


















# Figure 4 b)
# ============================ Figure 4(b): g_{3,2} (2-dec rounding) ============================
# Per-iteration reshaping:
#   MMiter <- matrix(MM[i, ], L, N);  MWiter <- matrix(MW[i, ], L, N)
#   alpha  = exp(t(MMiter))  # N x L (rows = neighborhoods, cols = deltas)
#   gamma  = exp(t(MWiter))  # N x L
# Inputs expected: 
#   - MW, MM: matrices with dims (iterations x (N*L))
#   - fim = c(tau1_end, tau2_end, tau3_end)
# Output: 
#   - tab_g32 (numeric) and tab_g32_fmt (strings): CI_95% for g3ini & g3fim and the symbol

# ---- Canonical neighborhood names (N = 34) ----
bairro_names <- c(
  "Lagoa Nova","Nova Descoberta","Candelária","Capim Macio","Pitimbu","Neópolis","Ponta Negra",
  "Santos Reis","Rocas","Ribeira","Praia do Meio","Petrópolis","Areia Preta","Cidade Alta",
  "Mãe Luiza","Alecrim","Barro Vermelho","Tirol","Lagoa Seca",
  "Cidade da Esperança","Quintas","Dix-Sept Rosado","Bom Pastor","Nossa Senhora de Nazaré",
  "Felipe Camarão","Cidade Nova","Guarapes","Planalto",
  "Igapó","Potengi","Nossa Senhora da Apresentação","Lagoa Azul","Pajuçara","Redinha"
)
stopifnot(length(bairro_names) == 34)

# ---- Basic dims and sanity checks ----
L <- 3
N <- length(bairro_names)
stopifnot(is.matrix(MM), is.matrix(MW))
stopifnot(ncol(MM) >= L*N, ncol(MW) >= L*N, nrow(MM) == nrow(MW))

# ---- Time anchors for the Δ2 -> Δ3 transition ----
# fim = c(tau1_end, tau2_end, tau3_end)
stopifnot(is.numeric(fim), length(fim) >= 3, isTRUE(fim[3] > fim[2]), isTRUE(fim[2] > fim[1]))
TBini <- fim[2] + 0.05 * (fim[3] - fim[2])  # 5% inside Δ3
TBfim <- fim[3]                              # end of Δ3

# ---- Accumulate g3ini and g3fim across iterations (rows = iters, cols = neighborhoods) ----
niter <- nrow(MW)
g3ini <- matrix(NA_real_, nrow = niter, ncol = N)
g3fim <- matrix(NA_real_, nrow = niter, ncol = N)

for (i in 1:niter) {
  # Reshape per-iteration, then transpose to N x L
  MMiter <- matrix(MM[i, 1:(L*N)], nrow = L, ncol = N)   # column-major fill
  MWiter <- matrix(MW[i, 1:(L*N)], nrow = L, ncol = N)

  alp <- exp(t(MMiter))  # N x L (alpha)
  gam <- exp(t(MWiter))  # N x L (gamma)

  # g_{3,2}(t) at TBini and TBfim, referenced to tau2_end = fim[2]
  # Formula: t^(α2) - τ2^(α2) - (γ3/γ2) * ( t^(α3) - τ2^(α3) )
  g3ini[i, ] <- TBini^(alp[, 2]) - fim[2]^(alp[, 2]) - (gam[, 3] / gam[, 2]) * (TBini^(alp[, 3]) - fim[2]^(alp[, 3]))
  g3fim[i, ] <- TBfim^(alp[, 2]) - fim[2]^(alp[, 2]) - (gam[, 3] / gam[, 2]) * (TBfim^(alp[, 3]) - fim[2]^(alp[, 3]))
}

# ---- 95% credible intervals (raw) per neighborhood ----
g3ini_ci_raw <- t(apply(g3ini, 2, quantile, c(0.025, 0.975), na.rm = TRUE))
g3fim_ci_raw <- t(apply(g3fim,  2, quantile, c(0.025, 0.975), na.rm = TRUE))

# ---- Round to TWO decimals BEFORE classification (closed intervals) ----
g3ini_ci <- round(g3ini_ci_raw, 2)
g3fim_ci <- round(g3fim_ci_raw,  2)
colnames(g3ini_ci) <- c("q2.5", "q97.5")
colnames(g3fim_ci) <- c("q2.5", "q97.5")

# ---- Closed-interval position relative to zero (inclusive at 0) ----
# "left"  if upper < 0;  "right" if lower > 0;  "cross" otherwise (includes touching 0)
side_of_zero_closed <- function(lower, upper) {
  if (is.na(lower) || is.na(upper)) return(NA_character_)
  if (upper < 0) return("left")
  if (lower > 0) return("right")
  "cross"
}
s_ini <- apply(g3ini_ci, 1, function(z) side_of_zero_closed(z[1], z[2]))
s_fim <- apply(g3fim_ci,  1, function(z) side_of_zero_closed(z[1], z[2]))

# ---- Classification rules (using the 2-dec rounded CLOSED intervals) ----
# 1) ini left  & fim left  -> "++"
# 2) (ini right OR ini cross) & fim left -> "+"
# 3) ini cross & fim cross -> "="
# 4) (ini left OR ini cross) & fim right -> "-"
# 5) ini right & fim right -> "--"
# Else -> "=" (conservative)
symbol <- character(N)
for (k in seq_len(N)) {
  if (s_ini[k] == "left"  && s_fim[k] == "left") {
    symbol[k] <- "++"
  } else if ((s_ini[k] %in% c("right","cross")) && s_fim[k] == "left") {
    symbol[k] <- "+"
  } else if (s_ini[k] == "cross" && s_fim[k] == "cross") {
    symbol[k] <- "="
  } else if ((s_ini[k] %in% c("left","cross")) && s_fim[k] == "right") {
    symbol[k] <- "-"
  } else if (s_ini[k] == "right" && s_fim[k] == "right") {
    symbol[k] <- "--"
  } else {
    symbol[k] <- "="
  }
}

# ---- Final table (no means; shows 2-dec rounded CIs) ----
tab_g32 <- data.frame(
  Neighborhood = bairro_names,
  g3ini_q2.5   = g3ini_ci[, "q2.5"],
  g3ini_q97.5  = g3ini_ci[, "q97.5"],
  g3fim_q2.5   = g3fim_ci[, "q2.5"],
  g3fim_q97.5  = g3fim_ci[, "q97.5"],
  symbol       = symbol,
  check.names  = FALSE
)

# ---- Pretty print (optional; fixed two decimals) ----
tab_g32_fmt <- transform(
  tab_g32,
  g3ini_q2.5  = sprintf("%.2f", g3ini_q2.5),
  g3ini_q97.5 = sprintf("%.2f", g3ini_q97.5),
  g3fim_q2.5  = sprintf("%.2f", g3fim_q2.5),
  g3fim_q97.5 = sprintf("%.2f", g3fim_q97.5)
)
print(tab_g32_fmt, row.names = FALSE)


## =================== Figura 4(b): mapa g_{3,2} com legenda completa e rótulo "b)" ===================

# --- Packages ---
pkgs <- c("sf","dplyr","ggplot2","osmdata","stringi","cowplot")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
library(sf); library(dplyr); library(ggplot2); library(osmdata); library(stringi); library(cowplot)

# --- Expected input: tab_g32 with columns Neighborhood and symbol ---
stopifnot(exists("tab_g32"), all(c("Neighborhood","symbol") %in% names(tab_g32)))

# --- Canonical neighborhood names (N = 34) ---
bairro_names <- c(
  "Lagoa Nova","Nova Descoberta","Candelária","Capim Macio","Pitimbu","Neópolis","Ponta Negra",
  "Santos Reis","Rocas","Ribeira","Praia do Meio","Petrópolis","Areia Preta","Cidade Alta",
  "Mãe Luiza","Alecrim","Barro Vermelho","Tirol","Lagoa Seca",
  "Cidade da Esperança","Quintas","Dix-Sept Rosado","Bom Pastor","Nossa Senhora de Nazaré",
  "Felipe Camarão","Cidade Nova","Guarapes","Planalto",
  "Igapó","Potengi","Nossa Senhora da Apresentação","Lagoa Azul","Pajuçara","Redinha"
)

# --- Name normalizer (ASCII, lower, strip punctuation) ---
normalize_str <- function(x) {
  x %>%
    stringi::stri_trans_general("Latin-ASCII") %>%
    tolower() %>%
    gsub("[^a-z0-9 ]", " ", .) %>%
    gsub("\\s+", " ", .) %>%
    trimws()
}

# --- Build 'plot_sf' from OSM if not already available ---
if (!exists("plot_sf")) {
  bb  <- osmdata::getbb("Natal, Rio Grande do Norte, Brazil", format_out = "sf_polygon")
  if (is.null(bb) || nrow(bb) == 0) stop("Could not resolve bbox for Natal/RN via OSM.")
  bbx <- sf::st_bbox(bb)

  servers <- c(
    "https://overpass-api.de/api/interpreter",
    "https://overpass.kumi.systems/api/interpreter"
  )
  fetch_osm <- function(bbox, query_fun) {
    last_err <- NULL
    for (srv in servers) {
      osmdata::set_overpass_url(srv)
      res <- try({ osmdata::osmdata_sf(query_fun(bbox)) }, silent = TRUE)
      if (!inherits(res, "try-error")) return(res)
      last_err <- res
    }
    stop("All Overpass servers failed. Last error: ", as.character(last_err))
  }
  q_admin_fun <- function(bbox) {
    opq(bbox = bbox) |>
      add_osm_feature("boundary", "administrative") |>
      add_osm_feature("admin_level", c("9","10"))
  }
  q_place_fun <- function(bbox) {
    opq(bbox = bbox) |>
      add_osm_feature("place", c("neighbourhood","suburb","quarter"))
  }

  od_admin <- fetch_osm(bbx, q_admin_fun)
  od_place <- fetch_osm(bbx, q_place_fun)

  standardize_name_geom <- function(x) {
    if (is.null(x) || nrow(x) == 0) return(NULL)
    if (!inherits(x, "sf")) x <- st_as_sf(x)
    candidates <- c("name","name:pt","name:en","official_name","short_name","alt_name")
    hit <- intersect(candidates, names(x))
    nm  <- if ("name" %in% hit) x[["name"]] else if (length(hit) > 0) x[[hit[1]]] else rep(NA_character_, nrow(x))
    st_sf(name = nm, geometry = st_geometry(x)) |> suppressWarnings(st_make_valid())
  }

  poly_list <- list(
    od_admin$osm_polygons, od_admin$osm_multipolygons,
    od_place$osm_polygons, od_place$osm_multipolygons
  ) |> lapply(standardize_name_geom)
  poly_list <- poly_list[!vapply(poly_list, is.null, logical(1))]
  if (length(poly_list) == 0) stop("No neighborhood-like polygons from OSM.")

  osmp <- do.call(rbind, poly_list) |> suppressWarnings(st_make_valid())
  osmp <- st_intersection(osmp, st_union(bb))

  osmp <- osmp |>
    mutate(name = ifelse(is.na(name), "", name),
           name_norm = normalize_str(name))

  bairro_norm <- normalize_str(bairro_names)
  match_idx   <- match(osmp$name_norm, bairro_norm)

  unmatched <- which(is.na(match_idx) & nzchar(osmp$name_norm))
  if (length(unmatched) > 0) {
    for (i in unmatched) {
      cand <- agrep(osmp$name_norm[i], bairro_norm, max.distance = 0.25)
      if (length(cand) == 1) match_idx[i] <- cand
    }
  }

  matched <- !is.na(match_idx)
  if (!any(matched)) stop("No OSM names matched the 34-neighborhood list.")

  plot_sf <- osmp[matched, , drop = FALSE] |>
    mutate(name = bairro_names[match_idx[matched]]) |>
    select(name, geometry) |>
    st_make_valid() |>
    mutate(.area = as.numeric(st_area(geometry))) |>
    group_by(name) |>
    slice_max(.area, n = 1, with_ties = FALSE) |>
    ungroup() |>
    select(-.area)
}

# --- Join classification; normalize potential "==" as "=" and set full level set ---
tab_g32  <- tab_g32 |> mutate(symbol = ifelse(symbol == "==", "=", symbol))
plot_g32 <- plot_sf |>
  left_join(tab_g32 |> select(Neighborhood, symbol),
            by = c("name" = "Neighborhood")) |>
  mutate(symbol = factor(symbol, levels = c("++","+","=","-","--")))

# --- Palette and legend setup (keep swatches for ++ and +, but no text labels) ---
pal <- c(
  "++" = "#A50F15",  # dark red
  "+"  = "#FB6A4A",  # light red
  "="  = "#F5F0FF",  # very light purple
  "-"  = "#74C476",  # light green
  "--" = "#238B45"   # dark green
)
legend_levels <- c("++","+","=","-","--")
legend_labels <- c("", "", "=", "-", "--")  # blank labels for ++ and +

# --- Map with legend on the right; only "b)" centered in the footer ---
p_g32_map <- ggplot(plot_g32) +
  geom_sf(aes(fill = symbol), color = "white", linewidth = 0.25, na.rm = TRUE) +
  scale_fill_manual(
    name   = expression(g[3,2]~"class"),
    values = pal,
    limits = legend_levels,
    labels = legend_labels,
    drop   = FALSE,
    na.value = "grey85",
    guide = guide_legend(
      direction = "vertical",
      byrow = TRUE,
      override.aes = list(color = NA),
      title.position = "top"
    )
  ) +
  coord_sf(datum = NA) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title    = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption  = element_blank(),
    plot.margin   = margin(t = 5.5, r = 5.5, b = 28, l = 5.5)  # space for the "b)" footer
  )

p_g32_b <- ggdraw(p_g32_map) +
  draw_label("b)", x = 0.50, y = 0.03, hjust = 0.5, vjust = 0,
             fontface = "bold", size = 14)

print(p_g32_b)

# (Optional) save
# ggsave("g32_map_panel_b.png", p_g32_b, width = 6.4, height = 7.0, dpi = 300)


# Figure 5 

## ======= Pool de TODAS as iterações das 6 cadeias =======
phiw_all <- c(Mod30$Mbw, Mod31$Mbw, Mod32$Mbw, Mod33$Mbw, Mod34$Mbw, Mod35$Mbw)  # phiw
phim_all <- c(Mod30$Mbm, Mod31$Mbm, Mod32$Mbm, Mod33$Mbm, Mod34$Mbm, Mod35$Mbm)  # phim

## (Opcional) checagens rápidas
#length(phiw_all); length(phim_all)


library(STPoissonSS)
data(natal_vehicle_theft_data)


natal_theft_times<-natal_vehicle_theft_data$natal_theft_times
sites<-natal_vehicle_theft_data$natal_sites_utm
Matdist<-as.matrix(dist(sites))

Matcorr<-round(exp(-median(phiw_all)*Matdist),2)
Matcorrshape<-round(exp(-median(phim_all)*Matdist),2)


## =================== Pré-requisitos ===================
library(ggplot2)
library(grid)  # para unit() na legenda

## =================== Matriz de correlação ===================
# Usa Matcorr existente; se não existir, recomputa (sem arredondar).
if (!exists("Matcorr")) {
  library(STPoissonSS)
  data(natal_vehicle_theft_data)
  sites   <- natal_vehicle_theft_data$natal_sites_utm
  Matdist <- as.matrix(dist(sites))
  stopifnot(exists("phiw_all"))
  Matcorr <- exp(-median(phiw_all) * Matdist)
}
stopifnot(is.matrix(Matcorr), nrow(Matcorr) == 34, ncol(Matcorr) == 34)

## =================== Layout (ordem canônica por zona) ===================
idx_blocks <- list(
  South = 1:7, East = 8:19, West = 20:28, North = 29:34
)

## =================== Dados longos com coordenadas "tabela" ===================
df <- data.frame(
  i = rep(1:34, each = 34),
  j = rep(1:34, 34),
  corr = as.vector(Matcorr)
)
df$y <- 35 - df$i   # topo = i=1
df$x <- df$j

## --------- Pintar somente o triângulo superior (inclui diagonal) ---------
df_upper <- subset(df, j >= i)

## Diagonal com rótulos 1..34
diag_df <- subset(df, i == j)
diag_df$label <- as.character(1:34)
diag_df$y <- 35 - diag_df$i
diag_df$x <- diag_df$j

## Linhas separadoras entre zonas
vlines <- c(7.5, 19.5, 28.5)
hlines <- c(27.5, 15.5, 6.5)

## Zonas no topo e na esquerda
zone_tops <- data.frame(
  zone = names(idx_blocks),
  x = sapply(idx_blocks, function(v) mean(v)),
  y = rep(34, 4)
)
zone_left <- data.frame(
  zone = names(idx_blocks),
  x = rep(1, 4),
  y = sapply(idx_blocks, function(v) 35 - mean(v))
)

## =================== Escala das cores (padrão ggplot2) ===================
min_corr <- min(Matcorr, na.rm = TRUE)

## =================== Figura ===================
p_tab <- ggplot(df_upper, aes(x = x, y = y, fill = corr)) +
  geom_tile(width = 1, height = 1) +
  geom_vline(xintercept = vlines, linewidth = 0.6, inherit.aes = FALSE) +
  geom_hline(yintercept = hlines, linewidth = 0.6, inherit.aes = FALSE) +
  geom_text(data = diag_df,
            aes(x = x, y = y, label = label),
            size = 3.5, inherit.aes = FALSE) +
  # Zonas no topo
  geom_text(data = zone_tops,
            aes(x = x, y = y, label = paste(zone, "Zone")),
            fontface = "bold", size = 4.0, vjust = -0.8,
            inherit.aes = FALSE) +
  # Zonas na esquerda
  geom_text(data = zone_left,
            aes(x = x, y = y, label = paste(zone, "Zone")),
            fontface = "bold", size = 4.0, hjust = 1.15,
            inherit.aes = FALSE) +
  # >>> escala padrão (branco → azul) do ggplot2
  scale_fill_gradient(limits = c(min_corr, 1), name = "Correlation") +
  coord_fixed(expand = FALSE, clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    axis.title      = element_blank(),
    axis.text       = element_blank(),
    panel.grid      = element_blank(),
    plot.margin     = margin(12, 24, 12, 60),
    legend.title    = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  guides(fill = guide_colorbar(
    direction      = "horizontal",
    title.position = "top",
    barwidth       = unit(10, "cm"),
    barheight      = unit(0.6, "cm"),
    ticks.colour   = "black"
  ))

p_tab


## =================== Pré-requisitos ===================
library(ggplot2)
library(grid)  # unit() na legenda

## =================== Distâncias e matriz de correlação (shape) ===================
# Garante Matdist; se não existir, recomputa a partir dos dados do pacote.
if (!exists("Matdist")) {
  library(STPoissonSS)
  data(natal_vehicle_theft_data)
  sites   <- natal_vehicle_theft_data$natal_sites_utm
  Matdist <- as.matrix(dist(sites))
}

# Precisamos de phim_all já criado previamente (todas as amostras concatenadas).
stopifnot(exists("phim_all"))

# Matriz de correlação para o parâmetro de forma M
Matcorrshape <- round(exp(-median(phim_all) * Matdist), 2)
stopifnot(is.matrix(Matcorrshape), nrow(Matcorrshape) == 34, ncol(Matcorrshape) == 34)

## =================== Layout (ordem canônica por zona) ===================
idx_blocks <- list(
  South = 1:7, East = 8:19, West = 20:28, North = 29:34
)

## =================== Dados longos com coordenadas "tabela" ===================
df <- data.frame(
  i = rep(1:34, each = 34),
  j = rep(1:34, 34),
  corr = as.vector(Matcorrshape)
)
df$y <- 35 - df$i   # topo = i=1
df$x <- df$j

## --------- Pintar somente o triângulo superior (inclui diagonal) ---------
df_upper <- subset(df, j >= i)

## Diagonal com rótulos 1..34
diag_df <- subset(df, i == j)
diag_df$label <- as.character(1:34)
diag_df$y <- 35 - diag_df$i
diag_df$x <- diag_df$j

## Linhas separadoras entre zonas
vlines <- c(7.5, 19.5, 28.5)
hlines <- c(27.5, 15.5, 6.5)

## Zonas no topo e na esquerda
zone_tops <- data.frame(
  zone = names(idx_blocks),
  x = sapply(idx_blocks, function(v) mean(v)),
  y = rep(34, 4)
)
zone_left <- data.frame(
  zone = names(idx_blocks),
  x = rep(1, 4),
  y = sapply(idx_blocks, function(v) 35 - mean(v))
)

## =================== Escala das cores (padrão ggplot2) ===================
min_corr <- min(Matcorrshape, na.rm = TRUE)

## =================== Figura ===================
p_shape <- ggplot(df_upper, aes(x = x, y = y, fill = corr)) +
  geom_tile(width = 1, height = 1) +
  geom_vline(xintercept = vlines, linewidth = 0.6, inherit.aes = FALSE) +
  geom_hline(yintercept = hlines, linewidth = 0.6, inherit.aes = FALSE) +
  geom_text(data = diag_df,
            aes(x = x, y = y, label = label),
            size = 3.5, inherit.aes = FALSE) +
  # Zonas no topo
  geom_text(data = zone_tops,
            aes(x = x, y = y, label = paste(zone, "Zone")),
            fontface = "bold", size = 4.0, vjust = -0.8,
            inherit.aes = FALSE) +
  # Zonas na esquerda
  geom_text(data = zone_left,
            aes(x = x, y = y, label = paste(zone, "Zone")),
            fontface = "bold", size = 4.0, hjust = 1.15,
            inherit.aes = FALSE) +
  scale_fill_gradient(limits = c(min_corr, 1), name = "Correlation") +
  coord_fixed(expand = FALSE, clip = "off") +
  theme_minimal(base_size = 12) +
  theme(
    axis.title      = element_blank(),
    axis.text       = element_blank(),
    panel.grid      = element_blank(),
    plot.margin     = margin(12, 24, 12, 60),
    legend.title    = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  guides(fill = guide_colorbar(
    direction      = "horizontal",
    title.position = "top",
    barwidth       = unit(10, "cm"),
    barheight      = unit(0.6, "cm"),
    ticks.colour   = "black"
  ))

p_shape
