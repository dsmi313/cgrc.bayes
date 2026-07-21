# Pin the public microdose data locally so rendering never depends on a live
# network fetch. Run once; commit data/pacutes.csv and data/PROVENANCE.txt.

URL <- paste0("https://raw.githubusercontent.com/szb37/mcrds_public/",
              "master/data/pacutes.csv")
DEST <- "data/pacutes.csv"
EXPECTED_SHA256 <- "86aa784528ee045c61fadf3eacfd3e1897d16aae9839cee7cb4bfe839a7cc4e3"

dir.create("data", showWarnings = FALSE)
utils::download.file(URL, DEST, mode = "wb")

got <- if (requireNamespace("digest", quietly = TRUE)) {
  digest::digest(DEST, algo = "sha256", file = TRUE)
} else {
  warning("install.packages('digest') to verify the checksum")
  NA_character_
}

if (!is.na(got) && !identical(got, EXPECTED_SHA256)) {
  stop("checksum mismatch.\n  expected: ", EXPECTED_SHA256,
       "\n  got:      ", got,
       "\nThe upstream file changed. Do NOT silently accept this - the ",
       "published numbers were computed against the expected version.")
}

writeLines(c(
  "source_url:  https://raw.githubusercontent.com/szb37/mcrds_public/master/data/pacutes.csv",
  "accessed:    2026-07-21",
  paste0("sha256:      ", EXPECTED_SHA256),
  "bytes:       977449",
  "records:     15623",
  "note:        public data from the self-blinding psychedelic microdose trial,",
  "             as re-analysed by Szigeti et al. (2023) Sci Rep 13:12107."
), "data/PROVENANCE.txt")

message("done: ", DEST)
