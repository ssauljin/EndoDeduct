## Aggregate insured loss S for replications that were run before predict_aggregate() was added to run_replication():
##   Rscript add_aggregate_metrics.R [results_dir] [n_cores]
## Each replication is regenerated from its seed (simulate_case is deterministic given the seed) and the stored
## estimates are plugged into predict_aggregate(); no model is refitted. The output, a list with one element per
## replication (case, rep, metrics_S), is saved as <results_dir>/aggregate_metrics.rds and read by summarize_sim.R.
args <- commandArgs(trailingOnly = TRUE)
resdir <- if (length(args) >= 1) args[1] else "results"
cores  <- if (length(args) >= 2) as.integer(args[2]) else 1L
source("sim_functions.R")
files <- list.files(resdir, "^case[0-9]+_rep[0-9]+\\.rds$", full.names = TRUE)
one <- function(file) {
  r <- readRDS(file)
  if (!is.null(r$metrics_S)) return(list(case = r$case, rep = r$rep, metrics_S = r$metrics_S))
  sim <- simulate_case(r$case, r$seed); ds <- sim_datasets(sim)
  ## stored estimates report the dependence parameters as rho; the prediction functions expect Kendall's tau
  f <- lapply(r$par, function(p) { if ("rho" %in% names(p)) p["rho"] <- rho2tau(p["rho"]); list(par = p) })
  stopifnot(length(ds$clm$y) == r$n_claims, length(ds$pos$y) == r$n_pos)       # the regenerated data are the original ones
  list(case = r$case, rep = r$rep, metrics_S = predict_aggregate(sim, ds, f))
}
agg <- parallel::mclapply(files, function(fl) tryCatch(one(fl), error = function(e) { message(fl, ": ", conditionMessage(e)); NULL }),
                          mc.cores = cores)
agg <- Filter(Negate(is.null), agg)
saveRDS(agg, file.path(resdir, "aggregate_metrics.rds"))
cat(length(agg), "of", length(files), "replications written to", file.path(resdir, "aggregate_metrics.rds"), "\n")
