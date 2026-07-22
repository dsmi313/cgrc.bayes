# Precompute the operating-characteristics lookup the Shiny app reads instantly.
# cgr_operating(n_trials = 500) fits 2000 posteriors (~9 s), so the app cannot
# call it on a slider - it reads this table and interpolates.
#
# Sweep: n x p_cg x true_effect x mu_aeb (expectancy magnitude), n_trials = 500.
# Expectancy is "a parameter of their data", not a constant, so three levels are
# swept: half / reference / double the microdose calibration.
#   - Reference level mu_aeb = 7.7: full grid (all n, all effects).
#   - New levels 3.85 and 15.4: reduced grid (drop n = 60 and the 4.5 effect,
#     the least informative cells) to keep the extra build near an hour.
#
# RESUMABLE: loads any existing inst/extdata/cgrc_lookup.rds, skips cells already
# present, computes only the rest, and saves incrementally. So this both builds
# the table from scratch and extends an existing one without recomputing it.
#
# Run from the repo root:  Rscript data-raw/build_lookup.R   (~1 h fresh; ~1.5 h
# to add the two new expectancy levels on top of the reference level).

suppressMessages(library(cgrc.bayes))

PCG_GRID <- seq(0.50, 0.95, by = 0.05)
N_FULL   <- c(60, 80, 120, 160, 200, 250, 300, 400, 500, 700, 1000)
N_RED    <- c(80, 120, 160, 200, 250, 300, 400, 500, 700, 1000)  # drop n = 60
EFF_FULL <- c(0, 1.5, 3, 4.5)
EFF_RED  <- c(0, 1.5, 3)                                          # drop 4.5
N_TRIALS <- 500
SEED     <- 1
DEST     <- "inst/extdata/cgrc_lookup.rds"

target <- rbind(
  expand.grid(n = N_FULL, p_cg = PCG_GRID, true_effect = EFF_FULL,
              mu_aeb = 7.7,          KEEP.OUT.ATTRS = FALSE),
  expand.grid(n = N_RED,  p_cg = PCG_GRID, true_effect = EFF_RED,
              mu_aeb = c(3.85, 15.4), KEEP.OUT.ATTRS = FALSE))

key <- function(d) paste(d$n, d$p_cg, d$true_effect, d$mu_aeb, sep = "|")

existing <- if (file.exists(DEST)) readRDS(DEST) else NULL
if (!is.null(existing) && !"mu_aeb" %in% names(existing)) existing$mu_aeb <- 7.7
have <- if (is.null(existing)) character(0) else unique(key(existing))
todo <- target[!key(target) %in% have, , drop = FALSE]
cat(sprintf("target cells: %d   already present: %d   to compute: %d\n",
            nrow(target), length(have), nrow(todo)))

acc <- if (is.null(existing)) list() else list(existing)
t0 <- Sys.time()
for (i in seq_len(nrow(todo))) {
  g  <- todo[i, ]
  op <- cgr_operating(n_trials = N_TRIALS, n = g$n, p_cg = g$p_cg,
                      mu_dte = g$true_effect, mu_aeb = g$mu_aeb,
                      noise = "all", seed = SEED)
  op$n <- g$n; op$p_cg <- g$p_cg; op$true_effect <- g$true_effect
  op$mu_aeb <- g$mu_aeb
  acc[[length(acc) + 1]] <- op
  if (i %% 20 == 0 || i == nrow(todo)) {              # incremental save
    dir.create("inst/extdata", showWarnings = FALSE, recursive = TRUE)
    saveRDS(do.call(rbind, acc), DEST)
    el <- as.numeric(Sys.time() - t0, units = "mins")
    cat(sprintf("  %4d/%d new cells  %.1f min  (~%.0f min for the new cells)\n",
                i, nrow(todo), el, el / i * nrow(todo)))
  }
}

lut <- do.call(rbind, acc); rownames(lut) <- NULL
attr(lut, "meta") <- list(
  n_trials = N_TRIALS, seed = SEED, noise = "all",
  p_cg_grid = PCG_GRID, mu_aeb_grid = c(3.85, 7.7, 15.4),
  built = Sys.time(), pkg_version = as.character(utils::packageVersion("cgrc.bayes")))
saveRDS(lut, DEST)
cat(sprintf("wrote %s  (%d rows, %d cells)\n", DEST, nrow(lut), nrow(lut) / 4))
