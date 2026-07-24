# Precompute the UNKNOWN-preserving operating-characteristics lookup that Panel A
# reads for interactive UNKNOWN-design exploration. This is the six-stratum
# analogue of data-raw/build_lookup.R and uses cgr_unknown_operating() under the
# UNKNOWN-aware generative model (sim_aeb_unknown): UNKNOWN rate equal across arms,
# UNKNOWN responders carry no expectancy (see reports/UNRESOLVED.md U10).
#
# Sweep: n x p_cg x true_effect x u, at mu_aeb = 7.7 (the reference expectancy).
# mu_aeb is held fixed because, under assumption A2, the adjusted estimator's
# power and false-favourable rate are (to first order) independent of the
# expectancy magnitude - it mainly moves the UNADJUSTED comparator, which the
# design panel reports but does not adjust. The grid is deliberately coarser than
# the binary one to keep the build near an hour; the app snaps true_effect and u
# to the nearest grid level and interpolates bilinearly in (n, p_cg), exactly as
# for the binary lookup.
#
# RESUMABLE: loads any existing inst/extdata/cgrc_unknown_lookup.rds, skips cells
# already present, computes the rest, and saves incrementally.
#
# Run from the repo root:  Rscript data-raw/build_unknown_lookup.R   (~50 min).

suppressMessages(library(cgrc.bayes))

N_GRID   <- c(60, 120, 200, 300, 500, 1000)
PCG_GRID <- c(0.50, 0.60, 0.70, 0.80, 0.90)
EFF_GRID <- c(0, 1.5, 3)
U_GRID   <- c(0.10, 0.20, 0.30, 0.40)
MU_AEB   <- 7.7
N_TRIALS <- 500
SEED     <- 1
DEST     <- "inst/extdata/cgrc_unknown_lookup.rds"

target <- expand.grid(n = N_GRID, p_cg = PCG_GRID, true_effect = EFF_GRID,
                      u = U_GRID, KEEP.OUT.ATTRS = FALSE)
key <- function(d) paste(d$n, d$p_cg, d$true_effect, d$u, sep = "|")

existing <- if (file.exists(DEST)) readRDS(DEST) else NULL
have <- if (is.null(existing)) character(0) else unique(key(existing))
todo <- target[!key(target) %in% have, , drop = FALSE]
cat(sprintf("target cells: %d   already present: %d   to compute: %d\n",
            nrow(target), length(have), nrow(todo)))

acc <- if (is.null(existing)) list() else list(existing)
t0 <- Sys.time()
for (i in seq_len(nrow(todo))) {
  g  <- todo[i, ]
  op <- cgr_unknown_operating(n_trials = N_TRIALS, n = g$n, p_cg = g$p_cg,
                              u = g$u, mu_dte = g$true_effect, mu_aeb = MU_AEB,
                              noise = "all", seed = SEED)
  op$n <- g$n; op$p_cg <- g$p_cg; op$true_effect <- g$true_effect
  op$mu_aeb <- MU_AEB                          # u is already a column of op
  acc[[length(acc) + 1]] <- op
  if (i %% 10 == 0 || i == nrow(todo)) {       # incremental save
    dir.create("inst/extdata", showWarnings = FALSE, recursive = TRUE)
    saveRDS(do.call(rbind, acc), DEST)
    el <- as.numeric(Sys.time() - t0, units = "mins")
    cat(sprintf("  %4d/%d new cells  %.1f min  (~%.0f min total)\n",
                i, nrow(todo), el, el / i * nrow(todo)))
  }
}

lut <- do.call(rbind, acc); rownames(lut) <- NULL
attr(lut, "meta") <- list(
  n_trials = N_TRIALS, seed = SEED, noise = "all", mu_aeb = MU_AEB,
  n_grid = N_GRID, p_cg_grid = PCG_GRID, eff_grid = EFF_GRID, u_grid = U_GRID,
  model = "UNKNOWN-aware AEB (A1 equal-arm u, A2 no-expectancy UNKNOWN)",
  built = Sys.time(), pkg_version = as.character(utils::packageVersion("cgrc.bayes")))
saveRDS(lut, DEST)
cat(sprintf("wrote %s  (%d rows, %d cells)\n", DEST, nrow(lut), nrow(lut) / 4))
