# Precompute the operating-characteristics lookup table that the Shiny app reads
# instantly. cgr_operating(n_trials = 500) fits 2000 posteriors and takes ~9 s,
# so the app cannot call it on a slider - it reads this table instead.
#
# Sweep: n x p_cg x true_effect, at n_trials = 500 (the count below which Monte
# Carlo error on coverage is the same size as the effect being judged). Each cell
# stores cgr_operating()'s full 4-row output plus the empty-stratum diagnostics.
#
# Run once from the repo root:  Rscript data-raw/build_lookup.R
# Commits to inst/extdata/cgrc_lookup.rds (shipped as package data via
# cgrc_lookup()). Takes ~1 hour.

suppressMessages(library(cgrc.bayes))

N_GRID   <- c(60, 80, 120, 160, 200, 250, 300, 400, 500, 700, 1000)
PCG_GRID <- seq(0.50, 0.95, by = 0.05)
EFF_GRID <- c(0, 1.5, 3, 4.5)
N_TRIALS <- 500
SEED     <- 1

grid <- expand.grid(n = N_GRID, p_cg = PCG_GRID, true_effect = EFF_GRID,
                    KEEP.OUT.ATTRS = FALSE)
cat(sprintf("cells: %d  (n_trials = %d each)\n", nrow(grid), N_TRIALS))

t0 <- Sys.time(); out <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  op <- cgr_operating(n_trials = N_TRIALS, n = g$n, p_cg = g$p_cg,
                      mu_dte = g$true_effect, noise = "all", seed = SEED)
  op$n <- g$n; op$p_cg <- g$p_cg; op$true_effect <- g$true_effect
  out[[i]] <- op
  if (i %% 10 == 0 || i == nrow(grid)) {
    el <- as.numeric(Sys.time() - t0, units = "mins")
    cat(sprintf("  %3d/%d cells  %.1f min elapsed  (~%.0f min total)\n",
                i, nrow(grid), el, el / i * nrow(grid)))
  }
}

cgrc_lookup <- do.call(rbind, out)
rownames(cgrc_lookup) <- NULL
attr(cgrc_lookup, "meta") <- list(
  n_trials = N_TRIALS, seed = SEED, noise = "all",
  n_grid = N_GRID, p_cg_grid = PCG_GRID, eff_grid = EFF_GRID,
  built = Sys.time(), pkg_version = as.character(utils::packageVersion("cgrc.bayes")))

dir.create("inst/extdata", showWarnings = FALSE, recursive = TRUE)
saveRDS(cgrc_lookup, "inst/extdata/cgrc_lookup.rds")
cat(sprintf("wrote inst/extdata/cgrc_lookup.rds  (%d rows, %d cells)\n",
            nrow(cgrc_lookup), nrow(cgrc_lookup) / 4))
