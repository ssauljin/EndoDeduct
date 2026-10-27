################################################################################
##  Section 3 simulation study with a DISCRETE deductible
##
##  Design (everything except the deductible as in the continuous-deductible design):
##    x = (1, X1, X2),  X1 ~ N(0,1),  X2 ~ Bernoulli(0.6)
##    D  : proportional-odds logit on K = 5 levels (500, 1000, 2500, 5000, 10000)
##         P(D <= d_j | x) = logistic(kappa_j + x'beta_D),  kappa = (0.1, 1.0, 1.8, 2.6),  beta_D = (-0.5, -0.5)
##         (no intercept: it is absorbed in the thresholds; the convention is that of the
##          empirical model, so a positive coefficient means a SMALLER deductible)
##    N* : NB(r, exp(x'beta_N)/r),  beta_N = (-0.5, 0.5, 0.2),  r = 1
##    Y* : Pareto(exp(x'beta_Y), alpha) = GB2 with sigma = alpha1 = 1,  beta_Y = (8.5, 0.2, 0.5),  alpha = 3
##         (the intercept is moved from 3 to 8.5 so that the fraction of losses below the
##          deductible, about 42%, is the same as in the continuous design, about 45%)
##    dependence: Gaussian copulas between the LATENT uniform U_D of the deductible and N*,
##         and between U_D and each Y*_j (conditionally independent given U_D, x);
##         Cases 1-4 as before; Case 5 (optional) adds a partial correlation rho_{NY|D} between
##         the latent normals of N* and Y* given D, i.e. a violation of N* _||_ Y* | D.
##
##  Models: D (ordinal), SEM1-3, FEM1, FEM2, FEM2G, SDM1-3, FDM1-2 (Gaussian copula), fitted by
##  the sequential scheme of Section 2.3 with the per-unit likelihoods of the estimation pipeline.
################################################################################
suppressPackageStartupMessages({ library(VineCopula); library(numDeriv); library(MASS) })
has_kable <- FALSE
cfg <- list(kmax = 16, latent_mu = "groundup", hess_r = 2, hess_d = 0.01)

## ---- functions of the estimation pipeline (sections 1-3: distributions, likelihoods, fit_stage, sandwich) ----
if (!exists("pipeline_file")) pipeline_file <- "endogenous_deductible_pipeline.R"
.src <- readLines(pipeline_file)
.i1 <- grep("^## 1\\. distributions", .src); .i4 <- grep("^## 4\\. marginal models", .src)
eval(parse(text = .src[(.i1 - 1):(.i4 - 2)]))
rm(.src, .i1, .i4)

## ---- design ------------------------------------------------------------------
##  Cases 1-4: the model of Section 2.2 with a discrete deductible. Given the deductible level, N* and each Y*_j are
##  generated from their own uniform positions inside the atom (F_D(d-1), F_D(d)], independently of each other, which
##  is exactly eq. (5)-(6) with N* _||_ Y* | D.
##  Case 5: a single latent uniform U_D drives the deductible level, N* and all Y*_j (a common risk propensity); given
##  the observed level, N* and Y* are then dependent through their common position inside the atom (violation of
##  N* _||_ Y* | D, although N* _||_ Y* | U_D).
##  Case 6: as Case 5, plus a direct dependence between N* and Y*: the latent normal of each Y*_j has correlation
##  rho_NY = 0.3 with the residual of the latent normal of N* (frequency-severity dependence given the deductible).
design <- list(
  n = 5000,
  levels = c(500, 1000, 2500, 5000, 10000),
  kappa  = c(0.1, 1.0, 1.8, 2.6), beta_D = c(-0.5, -0.5),
  beta_N = c(-0.5, 0.5, 0.2), r = 1,
  beta_Y = c(8.5, 0.2, 0.5), alpha = 3,
  cases  = list(`1` = list(rho_N =  0.0, rho_Y =  0.0, rho_NY = 0,   shared = FALSE),
                `2` = list(rho_N = -0.5, rho_Y =  0.0, rho_NY = 0,   shared = FALSE),
                `3` = list(rho_N =  0.0, rho_Y = -0.5, rho_NY = 0,   shared = FALSE),
                `4` = list(rho_N = -0.6, rho_Y = -0.4, rho_NY = 0,   shared = FALSE),
                `5` = list(rho_N = -0.6, rho_Y = -0.4, rho_NY = 0,   shared = TRUE),
                `6` = list(rho_N = -0.6, rho_Y = -0.4, rho_NY = 0.3, shared = TRUE))
)
K_ded <- length(design$levels)
XD    <- function(X, dedtype) X[, -1, drop = FALSE]     # ordinal logit: no intercept
rho2tau <- function(rho) 2 * asin(rho) / pi
tau2rho <- function(tau) sin(pi * tau / 2)

## ---- data generation ---------------------------------------------------------
## Pareto (Lomax) with scale lambda and shape alpha: S(y) = (1 + y/lambda)^(-alpha);  quantile  lambda ((1-u)^(-1/alpha) - 1)
simulate_case <- function(case, seed, des = design) {
  set.seed(seed)
  n <- des$n; cs <- des$cases[[as.character(case)]]
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.6)
  X  <- cbind("(Intercept)" = 1, X1 = x1, X2 = x2)
  ## deductible level from the proportional-odds model (latent uniform uD)
  uD <- runif(n)
  Fd <- cbind(0, plogis(outer(drop(X[, -1] %*% des$beta_D), des$kappa, "+")), 1)   # n x (K+1): 0, F_D(d_1|x), ..., F_D(d_{K-1}|x), 1
  lev <- 1L + rowSums(uD > Fd[, 2:K_ded, drop = FALSE]); d <- des$levels[lev]
  um <- Fd[cbind(seq_len(n), lev)]; up <- Fd[cbind(seq_len(n), lev + 1L)]        # the atom (F_D(d-1), F_D(d)] of the realised level
  ## latent uniform position of N* inside the atom: its own draw (Cases 1-4) or the one of the deductible (Cases 5-6)
  uN <- if (cs$shared) uD else runif(n, um, up)
  eN <- rnorm(n); zN <- cs$rho_N * qnorm(uN) + sqrt(1 - cs$rho_N^2) * eN
  muN <- drop(exp(X %*% des$beta_N)); Nstar <- qnbinom(pnorm(zN), size = des$r, mu = muN)
  ## severities: own uniform position per claim (Cases 1-4) or the common one (Cases 5-6); optional direct N*-Y* dependence
  unit <- rep(seq_len(n), Nstar); M <- length(unit)
  uY <- if (cs$shared) uD[unit] else runif(M, um[unit], up[unit])
  eY <- cs$rho_NY * eN[unit] + sqrt(1 - cs$rho_NY^2) * rnorm(M)
  zY <- cs$rho_Y * qnorm(uY) + sqrt(1 - cs$rho_Y^2) * eY
  lam <- drop(exp(X[unit, , drop = FALSE] %*% des$beta_Y))
  Ystar <- lam * (pnorm(zY, lower.tail = FALSE)^(-1 / des$alpha) - 1)
  dcl <- d[unit]; pos <- Ystar > dcl
  Ntru <- tabulate(unit[pos], nbins = n)
  list(case = case, seed = seed, n = n, X = X, lev = lev, d = d, Nstar = Nstar, Ntru = Ntru, uD = uD,
       clm = list(unit = unit, Ystar = Ystar, d = dcl, lev = lev[unit], y = pmax(Ystar - dcl, 0), pos = pos))
}
## data lists in the format of the pipeline
sim_datasets <- function(sim) {
  X <- sim$X; cl <- sim$clm
  dat_pol <- list(X = X, n_cen = sim$Nstar, n_tru = sim$Ntru, d = sim$d, dratio = NA, dlev = sim$lev, unit = seq_len(sim$n))
  dat_clm <- list(X = X[cl$unit, , drop = FALSE], y = cl$y, ystar = cl$Ystar, d = cl$d, dratio = NA, dlev = cl$lev, unit = cl$unit)
  keep <- cl$pos
  dat_pos <- list(X = X[cl$unit[keep], , drop = FALSE], y = cl$y[keep], ystar = cl$Ystar[keep], d = cl$d[keep], dratio = NA,
                  dlev = cl$lev[keep], unit = cl$unit[keep])
  list(pol = dat_pol, clm = dat_clm, pos = dat_pos)
}

## ---- stages with explicit data and the Pareto restriction sigma = alpha1 = 1 ----------------
expand_pareto <- function(th, k) c(th[1:k], 1, 1, th[(k + 1):length(th)])   # (beta, alpha, [tau]) -> (beta, 1, 1, alpha, [tau])
sim_stage_D <- function(dat) {
  X <- XD(dat$X, "disc")
  make_stage("D", "D_disc", function(th, prev) log(pmax(D_disc_pmf(th, X, dat$dlev), 0)), dat$unit,   # crossed thresholds -> -Inf (rejected by minimize)
             lower = rep(-Inf, K_ded - 1 + ncol(X)), upper = rep(Inf, K_ded - 1 + ncol(X)),
             par_names = c(paste0("kappa", 1:(K_ded - 1)), colnames(X)), spec = list(dedtype = "disc"))
}
sim_stage_sev <- function(name, scenario, dedtype, fam, dat) {
  cb <- cop_bounds(fam); k <- ncol(dat$X)
  make_stage("S", name, function(th, prev) sev_loglik_i(expand_pareto(th, k), prev, scenario, dedtype, fam, dat), dat$unit,
             lower = c(rep(-Inf, k), 1e-3, cb$lo), upper = c(rep(Inf, k), Inf, cb$hi),
             par_names = c(colnames(dat$X), "alpha", cb$nm), spec = list(scenario = scenario, dedtype = dedtype, fam = fam))
}
sim_stage_freq <- function(name, scenario, dedtype, fam, sev_spec, dat) {
  cb <- cop_bounds(fam); k <- ncol(dat$X)
  make_stage("N", name, function(th, prev) {
    if (!is.null(prev$S)) prev$S <- expand_pareto(prev$S, k)
    freq_loglik_i(th, prev, scenario, dedtype, fam, sev_spec, dat) }, dat$unit,
    lower = c(rep(-Inf, k), 1e-3, cb$lo), upper = c(rep(Inf, k), Inf, cb$hi),
    par_names = c(colnames(dat$X), "r", cb$nm), spec = list(scenario = scenario, dedtype = dedtype, fam = fam, sev_spec = sev_spec))
}

## ---- starting values (data-driven; the true values are never used) ----------------------------
init_D <- function(dat) {
  X <- XD(dat$X, "disc"); df <- data.frame(lev = factor(dat$dlev, levels = 1:K_ded), X)
  po <- tryCatch(suppressWarnings(polr(lev ~ ., data = df, Hess = FALSE)), error = function(e) NULL)
  if (is.null(po)) c(qlogis(cumsum(prop.table(table(df$lev)))[-K_ded]), rep(0, ncol(X))) else c(po$zeta, -coef(po))
}
init_sev  <- function(dat, ystar) c(coef(lm(log(ystar) ~ dat$X - 1)), 2)          # beta from a log-linear regression, alpha = 2
init_freq <- function(dat, n, offset = NULL) {
  b <- coef(glm(n ~ dat$X - 1 + offset(if (is.null(offset)) rep(0, length(n)) else log(offset)), family = poisson))
  m <- mean(n); v <- var(n); r0 <- if (v > m) m^2 / (v - m) else 10
  c(b, min(max(r0, 0.1), 20))
}

## ---- fit all models of one replication ----------------------------------------------------
fit_models <- function(ds, models = c("D", "SEM1", "SEM2", "SEM3", "FEM1", "FEM2", "FEM2G", "SDM1", "SDM2", "SDM3", "FDM1", "FDM2"),
                       fam = 1, verbose = FALSE) {
  f <- list(); tm <- c()
  go <- function(nm, expr) { t0 <- proc.time()[3]; f[[nm]] <<- expr; tm[nm] <<- proc.time()[3] - t0
                             if (verbose) cat(sprintf("  %-6s loglik %12.3f  conv %d  (%.1f s)\n", nm, f[[nm]]$loglik, f[[nm]]$convergence, tm[nm])) }
  pol <- ds$pol; clm <- ds$clm; pos <- ds$pos
  if ("D" %in% models)    go("D",    fit_stage(sim_stage_D(pol), init_D(pol), naive_hessian = FALSE))
  if ("SEM1" %in% models) go("SEM1", fit_stage(sim_stage_sev("SEM1", "complete",  "none", 0, clm), init_sev(clm, clm$ystar), naive_hessian = FALSE))
  if ("SEM2" %in% models) go("SEM2", fit_stage(sim_stage_sev("SEM2", "censored",  "none", 0, clm), init_sev(pos, pos$ystar), naive_hessian = FALSE))
  if ("SEM3" %in% models) go("SEM3", fit_stage(sim_stage_sev("SEM3", "truncated", "none", 0, pos), init_sev(pos, pos$ystar), naive_hessian = FALSE))
  if ("FEM1" %in% models) go("FEM1", fit_stage(sim_stage_freq("FEM1", "censored", "none", 0, NULL, pol), init_freq(pol, pol$n_cen), naive_hessian = FALSE))
  sp_none <- list(dedtype = "none", fam = 0)
  if ("FEM2" %in% models) {                        # p(d;x) from SEM3 (truncated marginal severity model)
    p3 <- retention_prob(expand_pareto(f$SEM3$par, ncol(pol$X)), NULL, sp_none, pol)
    go("FEM2", fit_stage(sim_stage_freq("FEM2", "truncated", "none", 0, sp_none, pol), init_freq(pol, pol$n_tru, p3),
                         prev = list(S = f$SEM3$par), naive_hessian = FALSE)) }
  if ("FEM2G" %in% models) {                       # p(d;x) from SEM1 (complete-data marginal severity model)
    p1 <- retention_prob(expand_pareto(f$SEM1$par, ncol(pol$X)), NULL, sp_none, pol)
    go("FEM2G", fit_stage(sim_stage_freq("FEM2G", "truncated", "none", 0, sp_none, pol), init_freq(pol, pol$n_tru, p1),
                          prev = list(S = f$SEM1$par), naive_hessian = FALSE)) }
  prevD <- list(D = f$D$par)
  if ("SDM1" %in% models) go("SDM1", fit_stage(sim_stage_sev("SDM1", "complete",  "disc", fam, clm), c(f$SEM1$par, 0), prev = prevD, naive_hessian = FALSE))
  if ("SDM2" %in% models) go("SDM2", fit_stage(sim_stage_sev("SDM2", "censored",  "disc", fam, clm), c(f$SEM2$par, 0), prev = prevD, naive_hessian = FALSE))
  if ("SDM3" %in% models) go("SDM3", fit_stage(sim_stage_sev("SDM3", "truncated", "disc", fam, pos), c(f$SEM3$par, 0), prev = prevD, naive_hessian = FALSE))
  if ("FDM1" %in% models) go("FDM1", fit_stage(sim_stage_freq("FDM1", "censored", "disc", fam, NULL, pol), c(f$FEM1$par, 0), prev = prevD, naive_hessian = FALSE))
  if ("FDM2" %in% models)                          # p(d;x) from SDM3 (truncated joint severity model)
    go("FDM2", fit_stage(sim_stage_freq("FDM2", "truncated", "disc", fam, list(dedtype = "disc", fam = fam), pol), c(f$FEM2$par, 0),
                         prev = list(D = f$D$par, S = f$SDM3$par), naive_hessian = FALSE))
  attr(f, "time") <- tm
  f
}

## ---- conditional means given the deductible atom (Gaussian copula) ---------------------------
## Gauss-Legendre nodes on [-1, 1]
GL64 <- gl_nodes(64)
## E[(Y* - d)+ | D = d, x]  and  P(Y* > d | D = d, x)  for a Pareto(lambda, alpha) severity linked to U_D by a
## Gaussian copula with correlation rho; up/um = F_D(d | x), F_D(d- | x)  (rho = 0 or up - um = 1: marginal case)
cond_sev <- function(beta, alpha, rho, X, d, up, um) {
  lam <- drop(exp(X %*% beta)); vd <- 1 - (1 + d / lam)^(-alpha); zd <- qnorm(vd)
  L <- 14; Z <- zd + L * outer(rep(1, length(zd)), (GL64$x + 1) / 2); W <- GL64$w * L / 2
  Y <- lam * (pnorm(Z, lower.tail = FALSE)^(-1 / alpha) - 1)
  dens <- dnorm(Z)
  if (rho != 0) {
    s <- sqrt(1 - rho^2); zp <- qnorm(clamp01(up)); zm <- qnorm(clamp01(um))
    dens <- dens * (pnorm((zp - rho * Z) / s) - pnorm((zm - rho * Z) / s)) / (up - um)
  }
  Ec <- drop(((Y - d) * dens) %*% W); p <- drop(dens %*% W)
  list(Ecens = Ec, p = p, Etrunc = Ec / p)
}
## E[N* | D = d, x] for NB(r, mu) linked to U_D by a Gaussian copula with correlation rho
cond_freq <- function(beta, r, rho, X, up, um, kmax = 200) {
  mu <- drop(exp(X %*% beta)); if (rho == 0) return(mu)
  n <- length(mu); K <- kmax + 1
  FN <- pnbinom(matrix(0:kmax, K, n), size = r, mu = matrix(mu, K, n, byrow = TRUE))
  Up <- rep(up, each = K); Um <- rep(um, each = K)
  P  <- matrix(cop_rect(Um, Up, as.vector(FN), 1, rho) / (Up - Um), K, n)        # P(N* <= k | atom)
  colSums(1 - P)                                                                 # sum_k P(N* > k | atom)
}

## ---- prediction / ratemaking error of each model ----------------------------------------------
## For every policy (frequency) and every claim (severity) the conditional mean of the observed response given
## (d, x) is computed under the TRUE parameters and under each fitted model; two error measures are returned:
##   MAE_obs  = mean |observed response - predicted mean|           (the MAE of the paper)
##   MAE_prem = mean |true conditional mean - predicted mean|       (ratemaking error, no irreducible noise)
predict_all <- function(sim, ds, f, des = design) {
  case <- sim$case; rho <- des$cases[[as.character(case)]]
  pol <- ds$pol; clm <- ds$clm; pos <- ds$pos; k <- ncol(pol$X)
  atoms <- function(thD, dat) list(up = D_disc_cdf(thD, XD(dat$X, "disc"), dat$dlev), um = D_disc_cdf(thD, XD(dat$X, "disc"), dat$dlev - 1))
  thD0 <- c(des$kappa, des$beta_D)
  a_pol0 <- atoms(thD0, pol); a_clm0 <- atoms(thD0, clm); a_pos0 <- atoms(thD0, pos)
  a_pol <- atoms(f$D$par, pol); a_clm <- atoms(f$D$par, clm); a_pos <- atoms(f$D$par, pos)
  ## true conditional means
  sev0_clm <- cond_sev(des$beta_Y, des$alpha, rho[["rho_Y"]], clm$X, clm$d, a_clm0$up, a_clm0$um)
  sev0_pos <- cond_sev(des$beta_Y, des$alpha, rho[["rho_Y"]], pos$X, pos$d, a_pos0$up, a_pos0$um)
  p0_pol   <- cond_sev(des$beta_Y, des$alpha, rho[["rho_Y"]], pol$X, pol$d, a_pol0$up, a_pol0$um)$p
  EN0      <- cond_freq(des$beta_N, des$r, rho[["rho_N"]], pol$X, a_pol0$up, a_pol0$um)
  true <- list(N_cen = EN0, N_tru = EN0 * p0_pol, Y_cen = sev0_clm$Ecens, Y_tru = sev0_pos$Etrunc)
  obs  <- list(N_cen = pol$n_cen, N_tru = pol$n_tru, Y_cen = clm$y, Y_tru = pos$y)
  ## fitted severity models: conditional means and retention probabilities
  sev_fit <- function(nm, dat, a) { th <- f[[nm]]$par; rho_h <- if (length(th) > k + 1) tau2rho(th[k + 2]) else 0
                                    cond_sev(th[1:k], th[k + 1], rho_h, dat$X, dat$d, a$up, a$um) }
  pS <- list(SEM1 = sev_fit("SEM1", pol, a_pol)$p, SEM3 = sev_fit("SEM3", pol, a_pol)$p,
             SDM1 = sev_fit("SDM1", pol, a_pol)$p, SDM3 = sev_fit("SDM3", pol, a_pol)$p)
  frq_fit <- function(nm) { th <- f[[nm]]$par; rho_h <- if (length(th) > k + 1) tau2rho(th[k + 2]) else 0
                            cond_freq(th[1:k], th[k + 1], rho_h, pol$X, a_pol$up, a_pol$um) }
  pred <- list(
    N_cen = list(FEM1 = frq_fit("FEM1"), FDM1 = frq_fit("FDM1")),
    N_tru = list(FEM1 = frq_fit("FEM1") * pS$SEM1, FEM2G = frq_fit("FEM2G") * pS$SEM1, FEM2 = frq_fit("FEM2") * pS$SEM3,
                 FDM1 = frq_fit("FDM1") * pS$SDM1, FDM2 = frq_fit("FDM2") * pS$SDM3),
    Y_cen = list(SEM1 = sev_fit("SEM1", clm, a_clm)$Ecens, SEM2 = sev_fit("SEM2", clm, a_clm)$Ecens,
                 SDM1 = sev_fit("SDM1", clm, a_clm)$Ecens, SDM2 = sev_fit("SDM2", clm, a_clm)$Ecens),
    Y_tru = list(SEM1 = sev_fit("SEM1", pos, a_pos)$Etrunc, SEM3 = sev_fit("SEM3", pos, a_pos)$Etrunc,
                 SDM1 = sev_fit("SDM1", pos, a_pos)$Etrunc, SDM3 = sev_fit("SDM3", pos, a_pos)$Etrunc))
  out <- list()
  for (q in names(pred)) for (m in names(pred[[q]])) {
    out[[paste(q, m, "obs",  sep = ".")]] <- mean(abs(obs[[q]]  - pred[[q]][[m]]))
    out[[paste(q, m, "prem", sep = ".")]] <- mean(abs(true[[q]] - pred[[q]][[m]]))
    ## relative premium error by deductible level (frequency: per policy; severity: per claim)
    lv <- if (q %in% c("N_cen", "N_tru")) pol$dlev else if (q == "Y_cen") clm$dlev else pos$dlev
    lv <- factor(lv, levels = 1:K_ded)
    out[[paste(q, m, "bylev", sep = ".")]] <- tapply(pred[[q]][[m]] - true[[q]], lv, mean) / tapply(true[[q]], lv, mean)
  }
  out
}

## ---- aggregate insured loss and its components (notation of the revised manuscript) ----------------
##  In the revised manuscript N, Y denote the ground-up frequency / severity (N*, Y* in the comments above) and
##  N*, Y* the number and amounts of the insurer's payments, so that the aggregate insured loss of a policy-year is
##  S = Y*_1 + ... + Y*_{N*} and its pure premium is
##      pi(d, x) = E(S | D = d, x) = E(N | d, x) E((Y - d)+ | d, x) = E(N* | d, x) E(Y* | d, x)   (either scenario).
##  Ratemaking schemes (frequency model x severity model):
##      exogenous  full FEM1 x SEM1,  censored FEM1 x SEM2,  truncated FEM2 x SEM3 (p(d;x) from SEM3)
##      endogenous full FDM1 x SDM1,  censored FDM1 x SDM2,  truncated FDM2 x SDM3 (p(d;x) from SDM3)
##  Under truncated reporting E(N*|d,x) E(Y*|d,x) = p E(N|d,x) E(Y-d|Y>d,d,x) = E(N|d,x) E((Y-d)+|d,x), so every
##  scheme is evaluated as E(N|d,x) x E((Y-d)+|d,x) with its own frequency and severity model.
S_schemes <- list(EXO_full = c("FEM1", "SEM1"), EXO_cen = c("FEM1", "SEM2"), EXO_tru = c("FEM2", "SEM3"),
                  END_full = c("FDM1", "SDM1"), END_cen = c("FDM1", "SDM2"), END_tru = c("FDM2", "SDM3"))
## E[Y | D = d, x] of the ground-up loss: Gauss-Legendre quadrature of the Pareto quantile against the
## (copula-weighted) standard normal density on z in (-9, 11)
cond_gu_sev <- function(beta, alpha, rho, X, up, um) {
  lam <- drop(exp(X %*% beta)); n <- length(lam)
  L <- 20; z <- -9 + L * (GL64$x + 1) / 2; W <- GL64$w * L / 2
  Q <- pnorm(z, lower.tail = FALSE)^(-1 / alpha) - 1                     # Pareto quantile with unit scale
  dens <- matrix(dnorm(z), n, length(z), byrow = TRUE)
  if (rho != 0) {
    s <- sqrt(1 - rho^2); zp <- qnorm(clamp01(up)); zm <- qnorm(clamp01(um))
    dens <- dens * (pnorm((zp - rho * matrix(z, n, length(z), byrow = TRUE)) / s) -
                    pnorm((zm - rho * matrix(z, n, length(z), byrow = TRUE)) / s)) / (up - um)
  }
  lam * drop(dens %*% (Q * W))
}
## error measures of an estimated conditional mean 'est' against the true one 'tru' (both per policy-year):
##   MAE = mean |est - tru|,  MRE = mean (est - tru) / tru,  MSRE = mean ((est - tru) / tru)^2, overall and by level
err_measures <- function(est, tru, lev, prefix) {
  re <- (est - tru) / tru; lv <- factor(lev, levels = 1:K_ded); out <- list()
  out[[paste0(prefix, ".mae")]]  <- mean(abs(est - tru)); out[[paste0(prefix, ".mae.bylev")]]  <- tapply(abs(est - tru), lv, mean)
  out[[paste0(prefix, ".mre")]]  <- mean(re);             out[[paste0(prefix, ".mre.bylev")]]  <- tapply(re, lv, mean)
  out[[paste0(prefix, ".msre")]] <- mean(re^2);           out[[paste0(prefix, ".msre.bylev")]] <- tapply(re^2, lv, mean)
  out
}
##  f: fitted models as returned by fit_models (dependence parameters on the Kendall's tau scale)
##  Returns, with every conditional mean evaluated per policy-year at (d_i, x_i):
##   S.<scheme>.{mae, mre, msre}[.bylev]  pure premium E(S|d,x) against the true premium
##   S.<scheme>.mae_obs                    mean |S_i - estimated premium| (S_i: realised aggregate payment)
##   N.<m>, Y.<m>                          ground-up frequency E(N|d,x) and severity E(Y|d,x)
##   Nstar_tru.<m>                         number of payments under truncation E(N*|d,x) = p(d;x) E(N|d,x)
##                                         (under censoring N* = N)
##   Ystar_cen.<m>, Ystar_tru.<m>          payment per claim E((Y-d)+|d,x) and per payment E(Y-d|Y>d,d,x)
predict_aggregate <- function(sim, ds, f, des = design) {
  rho <- des$cases[[as.character(sim$case)]]
  pol <- ds$pol; k <- ncol(pol$X)
  atoms <- function(thD) list(up = D_disc_cdf(thD, XD(pol$X, "disc"), pol$dlev), um = D_disc_cdf(thD, XD(pol$X, "disc"), pol$dlev - 1))
  a0 <- atoms(c(des$kappa, des$beta_D)); a <- atoms(f$D$par)
  rho_of <- function(th) if (length(th) > k + 1) tau2rho(th[k + 2]) else 0
  ## true conditional means
  sev0 <- cond_sev(des$beta_Y, des$alpha, rho[["rho_Y"]], pol$X, pol$d, a0$up, a0$um)
  true <- list(N = cond_freq(des$beta_N, des$r, rho[["rho_N"]], pol$X, a0$up, a0$um),
               Y = cond_gu_sev(des$beta_Y, des$alpha, rho[["rho_Y"]], pol$X, a0$up, a0$um),
               p = sev0$p, Ystar_cen = sev0$Ecens, Ystar_tru = sev0$Etrunc)
  true$Nstar_tru <- true$N * true$p; true$S <- true$N * true$Ystar_cen
  ## fitted models
  sev <- lapply(setNames(nm = c("SEM1", "SEM2", "SEM3", "SDM1", "SDM2", "SDM3")), function(m) {
    th <- f[[m]]$par; r <- c(cond_sev(th[1:k], th[k + 1], rho_of(th), pol$X, pol$d, a$up, a$um),
                             list(EY = cond_gu_sev(th[1:k], th[k + 1], rho_of(th), pol$X, a$up, a$um))) })
  frq <- lapply(setNames(nm = c("FEM1", "FEM2G", "FEM2", "FDM1", "FDM2")), function(m) {
    th <- f[[m]]$par; cond_freq(th[1:k], th[k + 1], rho_of(th), pol$X, a$up, a$um) })
  p_of <- c(FEM1 = "SEM1", FEM2G = "SEM1", FEM2 = "SEM3", FDM1 = "SDM1", FDM2 = "SDM3")   # retention prob. used with each frequency model
  S_obs <- numeric(sim$n); agg <- rowsum(ds$clm$y, ds$clm$unit); S_obs[as.integer(rownames(agg))] <- agg[, 1]   # realised S_i
  lev <- pol$dlev; out <- list()
  for (s in names(S_schemes)) {
    est <- frq[[S_schemes[[s]][1]]] * sev[[S_schemes[[s]][2]]]$Ecens
    out <- c(out, err_measures(est, true$S, lev, paste0("S.", s)))
    out[[paste0("S.", s, ".mae_obs")]] <- mean(abs(S_obs - est))
  }
  for (m in names(frq)) out <- c(out, err_measures(frq[[m]], true$N, lev, paste0("N.", m)),
                                 err_measures(frq[[m]] * sev[[p_of[[m]]]]$p, true$Nstar_tru, lev, paste0("Nstar_tru.", m)))
  for (m in names(sev)) out <- c(out, err_measures(sev[[m]]$EY, true$Y, lev, paste0("Y.", m)),
                                 if (m != "SEM3" && m != "SDM3") err_measures(sev[[m]]$Ecens, true$Ystar_cen, lev, paste0("Ystar_cen.", m)),
                                 if (m != "SEM2" && m != "SDM2") err_measures(sev[[m]]$Etrunc, true$Ystar_tru, lev, paste0("Ystar_tru.", m)))
  out$S_true_mean <- tapply(true$S, factor(lev, levels = 1:K_ded), mean); out$S_obs_mean <- mean(S_obs)
  out
}

## ---- one replication ----------------------------------------------------------------------
##  sandwich: names of the models whose Godambe s.e. (naive / robust / sandwich) are stored, for coverage studies
##            ("FEM2" is cheap; "SDM3" and "FDM2" roughly double the cost of a replication)
run_replication <- function(case, rep, outdir = "results", sandwich = character(0), verbose = FALSE, overwrite = FALSE) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  file <- file.path(outdir, sprintf("case%s_rep%04d.rds", case, rep))
  if (file.exists(file) && !overwrite) return(invisible(readRDS(file)))
  seed <- 100000L * as.integer(case) + rep
  t0 <- proc.time()[3]
  sim <- simulate_case(case, seed); ds <- sim_datasets(sim)
  assign("m_units", sim$n, envir = globalenv())
  f <- fit_models(ds, verbose = verbose)
  par <- lapply(f, function(s) { p <- s$par; if ("tau" %in% names(p)) { p["tau"] <- tau2rho(p["tau"]); names(p)[names(p) == "tau"] <- "rho" }; p })
  res <- list(case = case, rep = rep, seed = seed, par = par,
              loglik = sapply(f, `[[`, "loglik"), convergence = sapply(f, `[[`, "convergence"),
              n_claims = length(ds$clm$y), n_pos = length(ds$pos$y), lev_prop = prop.table(table(factor(sim$lev, levels = 1:K_ded))),
              kendall = c(D_N = cor(sim$lev, sim$Nstar, method = "kendall"), D_Y = cor(sim$clm$lev, sim$clm$Ystar, method = "kendall")),
              metrics = predict_all(sim, ds, f), metrics_S = predict_aggregate(sim, ds, f), time_fit = attr(f, "time"))
  if (length(sandwich)) {
    chains <- list(FEM2 = list(f$SEM3, f$FEM2), SDM3 = list(f$D, f$SDM3), FDM2 = list(f$D, f$SDM3, f$FDM2))
    res$se <- lapply(chains[sandwich], function(ch) { t <- last_stage_table(sandwich_se(ch, verbose = FALSE))
                                                     t[, c("parameter", "estimate", "se_naive", "se_robust", "se_sandwich")] })
  }
  res$time_total <- proc.time()[3] - t0
  saveRDS(res, file)
  invisible(res)
}
