################################################################################
##  Summaries of the discrete-deductible simulation: Mean / s.e. / RMSE / Bias tables (Appendix layout of the
##  paper), MAE tables, coverage of the naive vs Godambe s.e., box plots.
##  usage:  Rscript summarize_sim.R [results_dir] [tex_dir]
################################################################################
args <- commandArgs(trailingOnly = TRUE)
resdir <- if (length(args) >= 1) args[1] else "results"
texdir <- if (length(args) >= 2) args[2] else "tex_sim"
dir.create(texdir, showWarnings = FALSE)
pipeline_file <- "endogenous_deductible_pipeline.R"; source("sim_functions.R")

files <- list.files(resdir, "^case[0-9]+_rep[0-9]+\\.rds$", full.names = TRUE)
res <- lapply(files, readRDS)
cases <- sort(unique(sapply(res, `[[`, "case")))
R <- sapply(cases, function(cs) sum(sapply(res, `[[`, "case") == cs)); names(R) <- cases
cat("replications per case:\n"); print(R)

## ---- true values ------------------------------------------------------------------------------
truth <- function(model, case) {
  rho <- design$cases[[as.character(case)]]
  switch(substr(model, 1, 1),
         D = c(design$kappa, design$beta_D),
         S = c(design$beta_Y, design$alpha, if (grepl("DM", model)) rho[["rho_Y"]]),
         F = c(design$beta_N, design$r,     if (grepl("DM", model)) rho[["rho_N"]]))
}
par_mat <- function(model, case) t(sapply(res[sapply(res, `[[`, "case") == case], function(r) r$par[[model]]))
stats <- function(model, case) {
  P <- par_mat(model, case); th <- truth(model, case); Rc <- nrow(P)
  m <- colMeans(P); v <- colMeans(sweep(P, 2, m)^2)
  data.frame(parameter = colnames(P), truth = th, mean = m, se = sqrt(v / Rc), rmse = sqrt(colMeans(sweep(P, 2, th)^2)),
             bias = m - th, var = v, row.names = NULL)
}
f4 <- function(x) formatC(x, format = "f", digits = 4)

## ---- Appendix-style tables -----------------------------------------------------------------------
## groups: list(label = <latex label>, idx = <parameter indices>)
tabular_block <- function(models, case, groups) {
  S <- lapply(models, stats, case = case); nm <- length(models)
  head <- c("\\toprule",
            paste0("& ", paste(sprintf("\\multicolumn{2}{c}{%s}", models), collapse = " & "), " \\\\"),
            paste(sprintf("\\cmidrule(lr){%d-%d}", 2 * seq_len(nm), 2 * seq_len(nm) + 1), collapse = " "),
            paste0("& ", paste(rep("Mean & RMSE", nm), collapse = " & "), " \\\\"),
            paste0("& ", paste(rep("(s.e.) & (Bias)", nm), collapse = " & "), " \\\\"), "\\midrule")
  body <- c()
  for (g in seq_along(groups)) {
    idx <- groups[[g]]$idx
    for (j in seq_along(idx)) {
      i <- idx[j]
      lab <- if (j == 1) sprintf("\\multirow{%d}{*}{%s}", 2 * length(idx), groups[[g]]$label) else ""
      l1 <- paste(sapply(S, function(s) sprintf("%s & %s", f4(s$mean[i]), f4(s$rmse[i]))), collapse = " & ")
      l2 <- paste(sapply(S, function(s) sprintf("(%s) & (%s)", f4(s$se[i]), f4(s$bias[i]))), collapse = " & ")
      body <- c(body, sprintf("%s & %s \\\\", lab, l1), sprintf(" & %s \\\\", l2))
      if (j < length(idx)) body <- c(body, sprintf("\\cmidrule{2-%d}", 2 * nm + 1))
    }
    body <- c(body, if (g < length(groups)) "\\midrule" else paste0("\\hhline{", strrep("=", 2 * nm + 1), "}"))
  }
  c(sprintf("\\begin{tabular}{l%s}", strrep("c", 2 * nm)), head, body, "\\end{tabular}")
}
two_panel_table <- function(left, right, caption, label, widths = c(0.65, 0.65)) {
  c("\\begin{table}[!ht]", "\\footnotesize", "\\centering", "{\\color{blue}", "\\resizebox{\\textwidth}{!}{%",
    sprintf("    \\begin{minipage}{%s\\textwidth}", widths[1]), "        \\centering", paste0("        ", left), "    \\end{minipage}",
    "    \\hspace{0.01\\textwidth}",
    sprintf("    \\begin{minipage}{%s\\textwidth}", widths[2]), "        \\centering", paste0("        ", right), "    \\end{minipage}",
    "}", "}", sprintf("\\caption{%s}\\label{%s}", caption, label), "\\end{table}")
}
g_sev  <- list(list(label = "$\\beta_Y$", idx = 1:3), list(label = "$\\alpha$", idx = 4))
g_sevJ <- c(g_sev, list(list(label = "$\\rho_Y$", idx = 5)))
g_frq  <- list(list(label = "$\\beta_N$", idx = 1:3), list(label = "$r$", idx = 4))
g_frqJ <- c(g_frq, list(list(label = "$\\rho_N$", idx = 5)))
g_D    <- list(list(label = "$\\kappa$", idx = 1:4), list(label = "$\\beta_D$", idx = 5:6))

for (cs in cases) {
  tab <- two_panel_table(tabular_block(c("SEM1", "SEM2", "SEM3"), cs, g_sev), tabular_block(c("FEM1", "FEM2G", "FEM2"), cs, g_frq),
                         sprintf("Simulation result of SEM and FEM in Case %d (discrete deductible)", cs), sprintf("tab:AppSimExoCase%d", cs))
  writeLines(tab, file.path(texdir, sprintf("tab_app_case%d_marginal.tex", cs)))
  tab <- two_panel_table(tabular_block(c("SDM1", "SDM2", "SDM3"), cs, g_sevJ), tabular_block(c("FDM1", "FDM2"), cs, g_frqJ),
                         sprintf("Simulation result of SDM and FDM in Case %d (discrete deductible)", cs), sprintf("tab:AppSimEndCase%d", cs),
                         widths = c(0.65, 0.5))
  writeLines(tab, file.path(texdir, sprintf("tab_app_case%d_joint.tex", cs)))
}
## deductible model: one table, one column pair per case
{
  S <- lapply(cases, function(cs) stats("D", cs)); nm <- length(cases)
  head <- c("\\toprule", paste0("& ", paste(sprintf("\\multicolumn{2}{c}{Case %d}", cases), collapse = " & "), " \\\\"),
            paste(sprintf("\\cmidrule(lr){%d-%d}", 2 * seq_len(nm), 2 * seq_len(nm) + 1), collapse = " "),
            paste0("& ", paste(rep("Mean & RMSE", nm), collapse = " & "), " \\\\"),
            paste0("& ", paste(rep("(s.e.) & (Bias)", nm), collapse = " & "), " \\\\"), "\\midrule")
  body <- c(); labs <- c(sprintf("$\\kappa_%d$", 1:4), "$\\beta_{D,1}$", "$\\beta_{D,2}$")
  for (i in 1:6) {
    l1 <- paste(sapply(S, function(s) sprintf("%s & %s", f4(s$mean[i]), f4(s$rmse[i]))), collapse = " & ")
    l2 <- paste(sapply(S, function(s) sprintf("(%s) & (%s)", f4(s$se[i]), f4(s$bias[i]))), collapse = " & ")
    body <- c(body, sprintf("%s & %s \\\\", labs[i], l1), sprintf(" & %s \\\\", l2), if (i < 6) sprintf("\\cmidrule{2-%d}", 2 * nm + 1))
  }
  tab <- c("\\begin{table}[!ht]", "\\footnotesize", "\\centering", "{\\color{blue}", "\\resizebox{\\textwidth}{!}{%", sprintf("\\begin{tabular}{l%s}", strrep("c", 2 * nm)), head, body,
           paste0("\\hhline{", strrep("=", 2 * nm + 1), "}"), "\\end{tabular}", "}", "}",
           "\\caption{Simulation result of the ordinal deductible model (thresholds $\\kappa_1<\\dots<\\kappa_4$ and $\\bm{\\beta}_D$) in Cases 1--4}\\label{tab:AppSimDed}", "\\end{table}")
  writeLines(tab, file.path(texdir, "tab_app_D.tex"))
}

## ---- MAE tables (main text) -----------------------------------------------------------------------
metric <- function(case, key) mean(sapply(res[sapply(res, `[[`, "case") == case], function(r) r$metrics[[key]]))
mae_table <- function(cs) {
  fq <- list(cen = c("FEM1", "FDM1"), tru = c("FEM1", "FEM2", "FDM1", "FDM2"))
  sv <- list(cen = c("SEM1", "SEM2", "SDM1", "SDM2"), tru = c("SEM1", "SEM3", "SDM1", "SDM3"))
  row <- function(q, mods, type, dig) paste(sapply(mods, function(m) formatC(metric(cs, paste(q, m, type, sep = ".")), format = "f", digits = dig)), collapse = " & ")
  c("\\begin{table}[!h]", "\\renewcommand{\\tabcolsep}{0.8mm}", "\\renewcommand{\\arraystretch}{1}", "\\footnotesize", "\\centering", "{\\color{blue}",
    "\\begin{tabular}[h]{ccccccc}", "\\toprule",
    "\\multirow{2}{*}{Frequency} & \\multicolumn{2}{c}{Censored} & \\multicolumn{4}{c}{Truncated}\\\\",
    "\\cmidrule(l{3pt}r{3pt}){2-3} \\cmidrule(l{3pt}r{3pt}){4-7}",
    " & FEM1 & FDM1 & FEM1 & FEM2 & FDM1 & FDM2\\\\", "\\midrule",
    sprintf("MAE & %s & %s \\\\", row("N_cen", fq$cen, "obs", 4), row("N_tru", fq$tru, "obs", 4)),
    if (cs <= 4) sprintf("MAE$_{\\text{prem}}$ & %s & %s \\\\", row("N_cen", fq$cen, "prem", 4), row("N_tru", fq$tru, "prem", 4)),
    "\\hhline{=======}", "\\end{tabular}", "\\\\", "\\vspace{0.3cm}",
    "\\begin{tabular}[h]{ccccccccc}", "\\toprule",
    "\\multirow{2}{*}{Severity} & \\multicolumn{4}{c}{Censored} & \\multicolumn{4}{c}{Truncated}\\\\",
    "\\cmidrule(l{3pt}r{3pt}){2-5} \\cmidrule(l{3pt}r{3pt}){6-9}",
    " & SEM1 & SEM2 & SDM1 & SDM2 & SEM1 & SEM3 & SDM1 & SDM3\\\\", "\\midrule",
    sprintf("MAE & %s & %s \\\\", row("Y_cen", sv$cen, "obs", 1), row("Y_tru", sv$tru, "obs", 1)),
    if (cs <= 4) sprintf("MAE$_{\\text{prem}}$ & %s & %s \\\\", row("Y_cen", sv$cen, "prem", 1), row("Y_tru", sv$tru, "prem", 1)),
    "\\hhline{=========}", "\\end{tabular}", "}",
    sprintf("\\caption{MAE of the frequency and severity models in Case %d}\\label{tab:maecase%d}", cs, cs), "\\end{table}")
}
for (cs in cases) writeLines(mae_table(cs), file.path(texdir, sprintf("tab_mae_case%d.tex", cs)))

## ---- coverage of 95% intervals: naive vs robust vs sandwich s.e. ----------------------------------------
cov_rows <- list()
for (cs in cases) {
  rr <- res[sapply(res, `[[`, "case") == cs]; rr <- rr[!sapply(rr, function(r) is.null(r$se))]
  if (!length(rr)) next
  for (mod in names(rr[[1]]$se)) {
    th <- truth(mod, cs); pn <- rr[[1]]$se[[mod]]$parameter
    Etau <- t(sapply(rr, function(r) r$se[[mod]]$estimate)); E <- Etau
    if ("tau" %in% pn) { j <- which(pn == "tau"); E[, j] <- tau2rho(Etau[, j]) }                        # report rho = sin(pi tau / 2)
    for (type in c("se_naive", "se_robust", "se_sandwich")) {
      SE <- t(sapply(rr, function(r) r$se[[mod]][[type]]))
      if ("tau" %in% pn) SE[, j] <- SE[, j] * (pi / 2) * cos(pi * Etau[, j] / 2)                       # delta method tau -> rho
      cover <- colMeans(abs(E - matrix(th, nrow(E), ncol(E), byrow = TRUE)) <= qnorm(0.975) * SE)
      cov_rows[[length(cov_rows) + 1]] <- data.frame(case = cs, model = mod, type = sub("se_", "", type), parameter = sub("tau", "rho", pn),
                                                     mean_se = colMeans(SE), emp_sd = apply(E, 2, sd), coverage = cover, R = nrow(E))
    }
  }
}
cov <- do.call(rbind, cov_rows); write.csv(cov, file.path(texdir, "coverage.csv"), row.names = FALSE)
## LaTeX: (case, model) rows: empirical sd, mean naive / robust / sandwich s.e., coverage
cov_tex <- function(rows) {
  lines <- c("\\begin{table}[!ht]", "\\footnotesize", "\\centering", "{\\color{blue}", "\\begin{tabular}{llccccccc}", "\\toprule",
             "Model & Parameter & Emp. SD & \\multicolumn{3}{c}{Mean s.e.} & \\multicolumn{3}{c}{Coverage of 95\\% interval} \\\\",
             "\\cmidrule(lr){4-6} \\cmidrule(lr){7-9}", " & & & Naive & Robust & Sandwich & Naive & Robust & Sandwich \\\\", "\\midrule")
  Rmin <- Inf
  for (k in seq_len(nrow(rows))) {
    cs <- rows$case[k]; mod <- rows$model[k]
    d <- cov[cov$case == cs & cov$model == mod, ]; if (!nrow(d)) next; pn <- unique(d$parameter); Rmin <- min(Rmin, d$R)
    for (i in seq_along(pn)) {
      dn <- d[d$type == "naive" & d$parameter == pn[i], ]; dr <- d[d$type == "robust" & d$parameter == pn[i], ]; ds <- d[d$type == "sandwich" & d$parameter == pn[i], ]
      sub_ <- if (substr(mod, 1, 1) == "S") "Y" else "N"
      lab <- switch(pn[i], "(Intercept)" = sprintf("$\\beta_{%s,0}$", sub_), X1 = sprintf("$\\beta_{%s,1}$", sub_), X2 = sprintf("$\\beta_{%s,2}$", sub_),
                    alpha = "$\\alpha$", r = "$r$", rho = sprintf("$\\rho_%s$", sub_))
      lines <- c(lines, sprintf("%s & %s & %s & %s & %s & %s & %.3f & %.3f & %.3f \\\\", if (i == 1) sprintf("%s (Case %d)", mod, cs) else "", lab,
                                f4(dn$emp_sd), f4(dn$mean_se), f4(dr$mean_se), f4(ds$mean_se), dn$coverage, dr$coverage, ds$coverage))
    }
    lines <- c(lines, "\\midrule")
  }
  lines[length(lines)] <- "\\hhline{=========}"
  c(lines, "\\end{tabular}", "}",
    sprintf("\\caption{Monte Carlo standard deviation of the sequential estimates, mean estimated standard errors (naive: stage-wise Hessian; robust: stage-wise sandwich $\\bm{H}_{kk}^{-1}\\bm{J}_{kk}\\bm{H}_{kk}^{-1}$; sandwich: Godambe information $\\bm{H}^{-1}\\bm{J}\\bm{H}^{-\\top}$ of Section 2.3) and empirical coverage of the nominal 95\\%% confidence intervals ($R=%d$ replications); the dependence parameters are on the correlation scale (delta method)}\\label{tab:simcoverage}", Rmin),
    "\\end{table}")
}
rows <- data.frame(case = c(1, 4, 4), model = c("FEM2", "SDM3", "FDM2"))
if (nrow(cov)) writeLines(cov_tex(rows), file.path(texdir, "tab_coverage.tex"))


## ---- relative error of the premium components by deductible level (Case 4) -------------------------------
bylev_table <- function(cs) {
  rr <- res[sapply(res, `[[`, "case") == cs]
  rel <- function(key) { M <- sapply(rr, function(r) { v <- r$metrics[[key]]; out <- rep(NA_real_, K_ded); out[as.integer(names(v))] <- v; out }); rowMeans(M, na.rm = TRUE) }
  row <- function(key) paste(sub("^-0\\.0$", "0.0", sprintf("%+.1f", 100 * rel(key))), collapse = " & ")
  c("\\begin{table}[!ht]", "\\footnotesize", "\\centering", "{\\color{blue}", "\\begin{tabular}{llrrrrr}", "\\toprule",
    " & & \\multicolumn{5}{c}{Deductible level} \\\\", "\\cmidrule(lr){3-7}",
    sprintf(" & Model & %s \\\\", paste(format(design$levels, big.mark = ","), collapse = " & ")), "\\midrule",
    sprintf("\\multirow{2}{*}{Frequency, censored} & FEM1 & %s \\\\", row("N_cen.FEM1.bylev")), sprintf(" & FDM1 & %s \\\\", row("N_cen.FDM1.bylev")), "\\midrule",
    sprintf("\\multirow{4}{*}{Frequency, truncated} & FEM1 & %s \\\\", row("N_tru.FEM1.bylev")), sprintf(" & FEM2 & %s \\\\", row("N_tru.FEM2.bylev")),
    sprintf(" & FDM1 & %s \\\\", row("N_tru.FDM1.bylev")), sprintf(" & FDM2 & %s \\\\", row("N_tru.FDM2.bylev")), "\\midrule",
    sprintf("\\multirow{4}{*}{Severity, censored} & SEM1 & %s \\\\", row("Y_cen.SEM1.bylev")), sprintf(" & SEM2 & %s \\\\", row("Y_cen.SEM2.bylev")),
    sprintf(" & SDM1 & %s \\\\", row("Y_cen.SDM1.bylev")), sprintf(" & SDM2 & %s \\\\", row("Y_cen.SDM2.bylev")), "\\midrule",
    sprintf("\\multirow{4}{*}{Severity, truncated} & SEM1 & %s \\\\", row("Y_tru.SEM1.bylev")), sprintf(" & SEM3 & %s \\\\", row("Y_tru.SEM3.bylev")),
    sprintf(" & SDM1 & %s \\\\", row("Y_tru.SDM1.bylev")), sprintf(" & SDM3 & %s \\\\", row("Y_tru.SDM3.bylev")),
    "\\hhline{=======}", "\\end{tabular}", "}",
    sprintf("\\caption{Relative error (in \\%%) of the estimated frequency and severity components of the pure premium by deductible level in Case %d: the average over the policy-years (frequency) or claims (severity) at each level of the difference between the estimated and the true component, divided by the average true component at that level, averaged over the replications}\\label{tab:bylevel%d}", cs, cs),
    "\\end{table}")
}
for (cs in intersect(cases, 1:4)) writeLines(bylev_table(cs), file.path(texdir, sprintf("tab_bylevel_case%d.tex", cs)))

## ---- figures: box plots of the estimates by model, one panel per parameter ---------------------------------
pdf_boxes <- function(cs) {
  pdf(file.path(texdir, sprintf("simcase%d.pdf", cs)), width = 11, height = 6.5)
  par(mfrow = c(2, 5), mar = c(4.2, 3.2, 2.2, 0.6), mgp = c(2, 0.7, 0), cex.axis = 0.8)
  sev <- c("SEM1", "SEM2", "SEM3", "SDM1", "SDM2", "SDM3"); frq <- c("FEM1", "FEM2G", "FEM2", "FDM1", "FDM2")
  box <- function(models, i, main, th) {
    L <- lapply(models, function(m) { P <- par_mat(m, cs); if (i <= ncol(P)) P[, i] else NA })
    boxplot(L, names = models, las = 2, main = main, outline = TRUE, pch = 20, cex = 0.5, col = ifelse(grepl("DM", models), "lightsteelblue", "grey90"))
    abline(h = th, col = "red", lty = 2)
  }
  rho <- design$cases[[as.character(cs)]]
  box(sev, 1, expression(beta[Y*","*0]), design$beta_Y[1]); box(sev, 2, expression(beta[Y*","*1]), design$beta_Y[2]); box(sev, 3, expression(beta[Y*","*2]), design$beta_Y[3])
  box(sev, 4, expression(alpha), design$alpha); box(sev[4:6], 5, expression(rho[Y]), rho[["rho_Y"]])
  box(frq, 1, expression(beta[N*","*0]), design$beta_N[1]); box(frq, 2, expression(beta[N*","*1]), design$beta_N[2]); box(frq, 3, expression(beta[N*","*2]), design$beta_N[3])
  box(frq, 4, expression(r), design$r); box(frq[4:5], 5, expression(rho[N]), rho[["rho_N"]])
  dev.off()
}
for (cs in cases) pdf_boxes(cs)

## ---- aggregate insured loss S = Y*_1 + ... + Y*_{N*} (notation of the revised manuscript) -------------------------
##  metrics_S of every replication: from the replication files (runs with predict_aggregate) or from
##  <results_dir>/aggregate_metrics.rds (written by add_aggregate_metrics.R for the earlier runs)
aggS <- lapply(Filter(function(r) !is.null(r$metrics_S), res), function(r) list(case = r$case, rep = r$rep, metrics_S = r$metrics_S))
agg_file <- file.path(resdir, "aggregate_metrics.rds")
if (file.exists(agg_file)) {
  old <- Filter(function(a) !any(sapply(aggS, function(b) b$case == a$case && b$rep == a$rep)), readRDS(agg_file))
  aggS <- c(aggS, old)
}
casesS <- intersect(cases, 1:4)                                             # Cases 1-4 are reported in the manuscript
if (length(aggS) && length(casesS)) {
  RS <- sapply(casesS, function(cs) sum(sapply(aggS, `[[`, "case") == cs)); names(RS) <- casesS
  mS  <- function(cs, key) mean(sapply(aggS[sapply(aggS, `[[`, "case") == cs], function(a) a$metrics_S[[key]]))
  mSl <- function(cs, key) { M <- sapply(aggS[sapply(aggS, `[[`, "case") == cs], function(a) { v <- a$metrics_S[[key]]
                               out <- rep(NA_real_, K_ded); out[as.integer(names(v))] <- v; out }); rowMeans(M, na.rm = TRUE) }
  fmt_mae  <- function(x, dig = 1) formatC(x, format = "f", digits = dig, big.mark = ",")
  fmt_mre  <- function(x) formatC(100 * x, format = "f", digits = 1, big.mark = ",", flag = "+")
  fmt_msre <- function(x) ifelse(x < 10, formatC(x, format = "f", digits = 3), formatC(x, format = "f", digits = 0, big.mark = ","))
  sch  <- names(S_schemes)
  sch_lab <- sapply(S_schemes, paste, collapse = "$\\times$")
  head6 <- c(" & & \\multicolumn{3}{c}{Exogenous deductible} & \\multicolumn{3}{c}{Endogenous deductible} \\\\",
             "\\cmidrule(lr){3-5} \\cmidrule(lr){6-8}",
             " & & Full & Censored & Truncated & Full & Censored & Truncated \\\\",
             sprintf(" & & %s \\\\", paste(sch_lab, collapse = " & ")), "\\midrule")
  ## (1) overall MAE / MRE / MSRE of the pure premium E(S|D=d;x), Cases 1-4
  body <- c()
  for (cs in casesS) {
    body <- c(body,
      sprintf("\\multirow{3}{*}{Case %d} & MAE & %s \\\\", cs, paste(fmt_mae(sapply(sch, function(s) mS(cs, paste0("S.", s, ".mae"))))), collapse = " & ")),
      sprintf(" & MRE (\\%%) & %s \\\\", paste(fmt_mre(sapply(sch, function(s) mS(cs, paste0("S.", s, ".mre")))), collapse = " & ")),
      sprintf(" & MSRE & %s \\\\", paste(fmt_msre(sapply(sch, function(s) mS(cs, paste0("S.", s, ".msre")))), collapse = " & ")),
      if (cs != tail(casesS, 1)) "\\midrule")
  }
  writeLines(c("\\begin{table}[!ht]", "\\footnotesize", "\\centering", "{\\color{blue}", "\\resizebox{\\textwidth}{!}{%",
               "\\begin{tabular}{llrrrrrr}", "\\toprule", head6, body, "\\hhline{========}", "\\end{tabular}", "}", "}",
               sprintf("\\caption{MAE, MRE and MSRE of the estimated pure premium ${\\rm E}(S|D=d;\\bm{x})$ of the aggregate insured loss $S=Y^*_1+\\cdots+Y^*_{N^*}$ under the exogenous (marginal) and endogenous (joint) ratemaking schemes in Cases 1--4, averaged over the $R=%d$ replications. Each scheme combines the frequency and severity models indicated in the header, fitted to the full, censored or truncated data.}\\label{tab:simS}", min(RS)),
               "\\end{table}"), file.path(texdir, "tab_S_overall.tex"))
  ## (2) MAE and MSRE of E(S|D=d;x) by deductible level, Cases 1-4
  lv <- format(design$levels, big.mark = ",", trim = TRUE)
  body <- c()
  for (cs in casesS) {
    tm <- mSl(cs, "S_true_mean")
    body <- c(body, sprintf("\\multirow{7}{*}{Case %d} & True premium & %s & \\multicolumn{5}{c}{} \\\\", cs, paste(fmt_mae(tm), collapse = " & ")), "\\cmidrule(lr){2-12}")
    for (s in sch) body <- c(body, sprintf(" & %s %s & %s & %s \\\\", if (grepl("^EXO", s)) "Exo." else "Endo.",
                                           switch(sub("^[A-Z]+_", "", s), full = "full", cen = "censored", tru = "truncated"),
                                           paste(fmt_mae(mSl(cs, paste0("S.", s, ".mae.bylev"))), collapse = " & "),
                                           paste(fmt_msre(mSl(cs, paste0("S.", s, ".msre.bylev"))), collapse = " & ")))
    if (cs != tail(casesS, 1)) body <- c(body, "\\midrule")
  }
  writeLines(c("\\begin{table}[!ht]", "\\footnotesize", "\\centering", "{\\color{blue}", "\\resizebox{\\textwidth}{!}{%",
               "\\begin{tabular}{llrrrrrrrrrr}", "\\toprule",
               " & & \\multicolumn{5}{c}{MAE by deductible level} & \\multicolumn{5}{c}{MSRE by deductible level} \\\\",
               "\\cmidrule(lr){3-7} \\cmidrule(lr){8-12}", sprintf(" & Scheme & %s & %s \\\\", paste(lv, collapse = " & "), paste(lv, collapse = " & ")),
               "\\midrule", body, "\\hhline{============}", "\\end{tabular}", "}", "}",
               "\\caption{MAE and MSRE of the estimated pure premium ${\\rm E}(S|D=d;\\bm{x})$ by deductible level in Cases 1--4: the averages over the policy-years at each level, averaged over the replications. ``True premium'' is the average true pure premium at each level. The frequency and severity models of each scheme are those of Table~\\ref{tab:simS}.}\\label{tab:simSbylevel}",
               "\\end{table}"), file.path(texdir, "tab_S_bylevel.tex"))
  ## (3) figure: mean ratio of the estimated to the true premium, 1 + MRE, by deductible level (log scale)
  pdf(file.path(texdir, "simS_bylevel.pdf"), width = 10, height = 3.2)
  par(mfrow = c(1, length(casesS)), mar = c(4, 3.6, 2, 0.5), mgp = c(2.2, 0.7, 0), cex.axis = 0.8)
  cols <- c(EXO_full = "grey55", EXO_cen = "darkorange", EXO_tru = "firebrick", END_full = "grey20", END_cen = "steelblue", END_tru = "navy")
  ltys <- ifelse(grepl("^EXO", sch), 2, 1); names(ltys) <- sch
  ratio <- lapply(casesS, function(cs) sapply(sch, function(s) 1 + mSl(cs, paste0("S.", s, ".mre.bylev"))))
  yl <- range(unlist(ratio), 1, na.rm = TRUE)
  for (i in seq_along(casesS)) {
    matplot(1:K_ded, ratio[[i]], type = "b", log = "y", pch = 19, cex = 0.6, col = cols[sch], lty = ltys, ylim = yl, xaxt = "n",
            xlab = "Deductible", ylab = if (i == 1) "Estimated / true premium" else "", main = sprintf("Case %d", casesS[i]))
    axis(1, at = 1:K_ded, labels = lv); abline(h = 1, col = "red", lty = 3)
    if (i == 1) legend("topleft", legend = c("Exo. full", "Exo. censored", "Exo. truncated", "Endo. full", "Endo. censored", "Endo. truncated"),
                       col = cols[sch], lty = ltys, pch = 19, cex = 0.62, bty = "n")
  }
  dev.off()
  ## (4) Appendix: MAE and MSRE of the components N, Y (ground-up) and N*, Y* (insurer's payments), Cases 1-4
  comp <- list(list(key = "N", lab = "$N$ (ground-up frequency)", models = c("FEM1", "FEM2G", "FEM2", "FDM1", "FDM2"), dig = 4),
               list(key = "Y", lab = "$Y$ (ground-up severity)", models = c("SEM1", "SEM2", "SEM3", "SDM1", "SDM2", "SDM3"), dig = 1),
               list(key = "Nstar_tru", lab = "$N^*$ (truncated)", models = c("FEM1", "FEM2G", "FEM2", "FDM1", "FDM2"), dig = 4),
               list(key = "Ystar_cen", lab = "$Y^*$ (censored)", models = c("SEM1", "SEM2", "SDM1", "SDM2"), dig = 1),
               list(key = "Ystar_tru", lab = "$Y^*$ (truncated)", models = c("SEM1", "SEM3", "SDM1", "SDM3"), dig = 1))
  nc <- length(casesS); body <- c()
  for (g in seq_along(comp)) {
    cg <- comp[[g]]
    for (j in seq_along(cg$models)) {
      m <- cg$models[j]; key <- paste(cg$key, m, sep = ".")
      vals <- unlist(lapply(casesS, function(cs) c(fmt_mae(mS(cs, paste0(key, ".mae")), cg$dig), fmt_msre(mS(cs, paste0(key, ".msre"))))))
      body <- c(body, sprintf("%s & %s & %s \\\\", if (j == 1) sprintf("\\multirow{%d}{*}{%s}", length(cg$models), cg$lab) else "", m, paste(vals, collapse = " & ")))
    }
    body <- c(body, if (g < length(comp)) "\\midrule" else paste0("\\hhline{", strrep("=", 2 + 2 * nc), "}"))
  }
  writeLines(c("\\begin{table}[!ht]", "\\footnotesize", "\\centering", "{\\color{blue}", "\\resizebox{\\textwidth}{!}{%",
               sprintf("\\begin{tabular}{ll%s}", strrep("rr", nc)), "\\toprule",
               sprintf("Component & Model & %s \\\\", paste(sprintf("\\multicolumn{2}{c}{Case %d}", casesS), collapse = " & ")),
               paste(sprintf("\\cmidrule(lr){%d-%d}", 1 + 2 * seq_len(nc), 2 + 2 * seq_len(nc)), collapse = " "),
               sprintf(" & & %s \\\\", paste(rep("MAE & MSRE", nc), collapse = " & ")), "\\midrule", body, "\\end{tabular}", "}", "}",
               "\\caption{MAE and MSRE of the estimated conditional means of the ground-up frequency ${\\rm E}(N|D=d;\\bm{x})$ and severity ${\\rm E}(Y|D=d;\\bm{x})$, of the number of payments under truncated reporting ${\\rm E}(N^*|D=d;\\bm{x})=p(d;\\bm{x}){\\rm E}(N|D=d;\\bm{x})$ (under censored reporting $N^*=N$), and of the payment per claim under censored reporting ${\\rm E}(Y^*|D=d;\\bm{x})={\\rm E}((Y-d)_+|D=d;\\bm{x})$ and truncated reporting ${\\rm E}(Y^*|D=d;\\bm{x})={\\rm E}(Y-d|Y>d,D=d;\\bm{x})$, evaluated for every policy-year and averaged over the replications. The truncated frequency models use $p(d;\\bm{x})$ from SEM1 (FEM1, FEM2G), SEM3 (FEM2), SDM1 (FDM1) and SDM3 (FDM2).}\\label{tab:AppSimComponents}",
               "\\end{table}"), file.path(texdir, "tab_app_components.tex"))
  ## (5) Appendix: realised-loss MAE, mean |S_i - estimated premium|
  writeLines(c("\\begin{table}[!ht]", "\\footnotesize", "\\centering", "{\\color{blue}", "\\begin{tabular}{lrrrrrr}", "\\toprule", " & \\multicolumn{3}{c}{Exogenous deductible} & \\multicolumn{3}{c}{Endogenous deductible} \\\\",
               "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7}", " & Full & Censored & Truncated & Full & Censored & Truncated \\\\", "\\midrule",
               sapply(casesS, function(cs) sprintf("Case %d & %s \\\\", cs, paste(fmt_mae(sapply(sch, function(s) mS(cs, paste0("S.", s, ".mae_obs"))))), collapse = " & "))),
               "\\hhline{=======}", "\\end{tabular}", "}",
               "\\caption{Mean absolute difference between the realised aggregate insured loss $S$ of a policy-year and its estimated pure premium ${\\rm E}(S|D=d;\\bm{x})$ (averaged over the replications); unlike Table~\\ref{tab:simS}, this measure contains the irreducible variability of $S$.}\\label{tab:AppSimSobs}",
               "\\end{table}"), file.path(texdir, "tab_app_S_obs.tex"))
}

## ---- key numbers for the text ----------------------------------------------------------------------------
sink(file.path(texdir, "key_numbers.txt"))
cat("replications per case:\n"); print(R)
cat("\nconvergence failures (any model):", sum(sapply(res, function(r) any(r$convergence != 0))), "\n")
cat("mean number of claims:", mean(sapply(res, `[[`, "n_claims")), " mean number of payments (Y*>d):", mean(sapply(res, `[[`, "n_pos")), "\n")
cat("fraction of losses below the deductible by case:\n"); print(tapply(sapply(res, function(r) 1 - r$n_pos / r$n_claims), sapply(res, `[[`, "case"), mean))
cat("mean level proportions:\n"); print(round(Reduce("+", lapply(res, `[[`, "lev_prop")) / length(res), 3))
cat("mean Kendall tau (level vs N*, level vs Y*) by case:\n"); print(t(sapply(cases, function(cs) colMeans(t(sapply(res[sapply(res, `[[`, "case") == cs], `[[`, "kendall"))))))
cat("mean time per replication (s) by case:\n"); print(tapply(sapply(res, `[[`, "time_total"), sapply(res, `[[`, "case"), mean))
for (cs in cases) { cat(sprintf("\n==== Case %d ====\n", cs))
  for (m in c("D", "SEM1", "SEM2", "SEM3", "FEM1", "FEM2G", "FEM2", "SDM1", "SDM2", "SDM3", "FDM1", "FDM2")) { cat(m, "\n"); print(stats(m, cs), digits = 4) }
  cat("MAE / premium MAE:\n"); keys <- grep("obs$|prem$", names(res[[1]]$metrics), value = TRUE)
  print(round(sapply(keys, function(k) metric(cs, k)), 4))
  cat("relative premium error by deductible level (mean over replications):\n")
  keys <- grep("bylev$", names(res[[1]]$metrics), value = TRUE)
  print(round(t(sapply(keys, function(k) rowMeans(sapply(res[sapply(res, `[[`, "case") == cs], function(r) { v <- r$metrics[[k]]; out <- rep(NA_real_, K_ded); out[as.integer(names(v))] <- v; out }), na.rm = TRUE))), 4))
}
if (nrow(cov)) { cat("\ncoverage:\n"); print(cov, digits = 4) }
if (length(aggS)) for (cs in casesS) {
  cat(sprintf("\n==== aggregate insured loss S, Case %d (R = %d) ====\n", cs, RS[as.character(cs)]))
  print(round(sapply(sch, function(s) c(mae = mS(cs, paste0("S.", s, ".mae")), mre = mS(cs, paste0("S.", s, ".mre")),
                                        msre = mS(cs, paste0("S.", s, ".msre")), mae_obs = mS(cs, paste0("S.", s, ".mae_obs")))), 4))
  cat("by level (true premium, then MAE / MRE / MSRE of each scheme):\n"); print(round(mSl(cs, "S_true_mean"), 2))
  for (s in sch) { cat(s, "\n"); print(round(rbind(mae = mSl(cs, paste0("S.", s, ".mae.bylev")), mre = mSl(cs, paste0("S.", s, ".mre.bylev")),
                                                   msre = mSl(cs, paste0("S.", s, ".msre.bylev"))), 4)) }
}
sink()
cat("written to", texdir, "\n")
