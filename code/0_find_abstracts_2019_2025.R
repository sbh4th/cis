# 0_find_abstracts_2019_2025.R
#
# Fetch PubMed abstracts for the extended date range 2019-2025.
# This script is a drop-in replacement / extension of the original
# 0_find_abstracts.R from agbarnett/intervals.
#
# The original study covered abstracts up to end-2018 / early 2019.
# This script targets 2019-01-01 to 2025-12-31 so that new data
# can be appended to (or analysed alongside) the Georgescu/Wren dataset.
#
# Approach
# --------
# * Uses the rentrez package to query NCBI E-utilities.
# * Fetches abstracts in batches to stay within API rate limits.
# * Saves results as an .RData file with the same 8-column structure
#   as Georgescu.Wren.RData so the existing analysis code works unchanged.
#
# Output
# ------
#   data/abstracts_2019_2025.RData   — data frame `abstracts_new`
#   data/abstracts_2019_2025.txt     — first 30 rows, tab-delimited (preview)
#
# Requirements
# ------------
#   install.packages(c("rentrez", "XML", "dplyr", "stringr"))
#
# NCBI API key (recommended for higher rate limits, 10 req/s vs 3 req/s):
#   Register at https://www.ncbi.nlm.nih.gov/account/ then set:
#   Sys.setenv(ENTREZ_KEY = "your_key_here")
#   or add it to your ~/.Renviron file.
# ---------------------------------------------------------------------------

library(rentrez)
library(XML)
library(dplyr)
library(stringr)

# ── Configuration ────────────────────────────────────────────────────────────

START_YEAR  <- 2019
END_YEAR    <- 2025
BATCH_SIZE  <- 200      # records per efetch call (max 10000, keep modest)
SLEEP_SEC   <- 0.34     # pause between requests (~3 req/s without API key)
                        # set to 0.11 if you have an NCBI API key
OUT_DIR     <- "data"

# Set your NCBI API key here OR in your ~/.Renviron as ENTREZ_KEY=...
# ncbi_key <- "YOUR_KEY_HERE"
# set_entrez_key(ncbi_key)

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ── PubMed search query ──────────────────────────────────────────────────────
# Mirrors the original study:
#   - Journal articles with abstracts
#   - English language
#   - Date window: 2019-2025
#   - Must contain text suggesting a confidence interval
#
# The original paper used Medline (PubMed) abstracts broadly; we keep that
# scope rather than restricting to specific MeSH terms.

build_query <- function(year) {
  paste0(
    '("', year, '/01/01"[PDAT] : "', year, '/12/31"[PDAT]) ',
    'AND "journal article"[PT] ',
    'AND hasabstract ',
    'AND English[LA] ',
    'AND ("confidence interval"[TIAB] OR "confidence intervals"[TIAB] ',
    'OR "95% CI"[TIAB] OR "95% confidence"[TIAB])'
  )
}

# ── Helper: parse a single <PubmedArticle> XML node ─────────────────────────

parse_article <- function(node) {
  tryCatch({
    pmid    <- xmlValue(node[["MedlineCitation"]][["PMID"]])
    journal <- xmlValue(
      node[["MedlineCitation"]][["Article"]][["Journal"]][["Title"]]
    )
    year_node <- node[["MedlineCitation"]][["Article"]][["Journal"]][
      ["JournalIssue"]][["PubDate"]][["Year"]]
    pub_year <- if (!is.null(year_node)) as.integer(xmlValue(year_node)) else NA_integer_

    abstract_nodes <- getNodeSet(node, ".//AbstractText")
    abstract_text  <- paste(sapply(abstract_nodes, xmlValue), collapse = " ")

    list(
      pmid     = pmid,
      journal  = journal,
      year     = pub_year,
      abstract = abstract_text
    )
  }, error = function(e) NULL)
}

# ── Helper: extract confidence intervals from free text ─────────────────────
#
# Regex strategy follows the Barnett 1_find.intervals.R logic:
#   Match patterns like:
#     "1.23 (95% CI: 0.98 to 1.56)"
#     "1.23 (95% CI 0.98–1.56)"
#     "1.23 [95% CI: 0.98, 1.56]"
#     "OR 1.23, 95% CI 0.98-1.56"
#     etc.
#
# Returns a data.frame with columns: lower, mean, upper
# (multiple CIs may be extracted from one abstract)

extract_cis <- function(text) {

  # Normalise: replace en-dash, em-dash, minus variants with ASCII hyphen
  text <- str_replace_all(text, "\u2013|\u2014|\u2212", "-")
  # Normalise 'to' and commas as separators inside intervals
  text <- str_replace_all(text, "\\bto\\b", "-")

  # Number pattern (allow negative, decimals, scientific notation)
  num  <- "-?[0-9]+(?:\\.[0-9]+)?(?:[eE][+-]?[0-9]+)?"

  # Pattern A: mean (95% CI lower[-,]upper)  — brackets optional style
  patA <- paste0(
    "(", num, ")",                          # capture: mean / point estimate
    "\\s*[\\(\\[]",                         # opening bracket
    "\\s*(?:95%|90%|99%)?\\s*CI[:\\s]+",   # "CI" label
    "(", num, ")",                          # capture: lower
    "\\s*[-,]\\s*",                         # separator
    "(", num, ")",                          # capture: upper
    "\\s*[\\)\\]]"                          # closing bracket
  )

  # Pattern B: "95% CI mean (lower-upper)"  – less common but present
  patB <- paste0(
    "(?:95%|90%|99%)\\s*CI\\s+",
    "(", num, ")",                          # capture: mean (actually lower in this form)
    "\\s*[\\(\\[]",
    "(", num, ")",
    "\\s*[-,]\\s*",
    "(", num, ")",
    "\\s*[\\)\\]]"
  )

  # Pattern C: bare triplet "lower to upper" preceded by known keywords
  # e.g. "hazard ratio 1.23 (0.98 to 1.56)"
  patC <- paste0(
    "(", num, ")",                          # capture: mean
    "\\s*[\\(\\[]",
    "(", num, ")",                          # capture: lower
    "\\s+to\\s+",
    "(", num, ")",                          # capture: upper
    "\\s*[\\)\\]]"
  )

  results <- list()

  for (pat in c(patA, patC)) {
    m <- str_match_all(text, regex(pat, ignore_case = TRUE))[[1]]
    if (nrow(m) > 0) {
      df <- data.frame(
        mean  = as.numeric(m[, 2]),
        lower = as.numeric(m[, 3]),
        upper = as.numeric(m[, 4]),
        stringsAsFactors = FALSE
      )
      results <- c(results, list(df))
    }
  }

  # Pattern B has different column order
  m2 <- str_match_all(text, regex(patB, ignore_case = TRUE))[[1]]
  if (nrow(m2) > 0) {
    df2 <- data.frame(
      mean  = as.numeric(m2[, 2]),   # first number = lower in patB
      lower = as.numeric(m2[, 3]),
      upper = as.numeric(m2[, 4]),
      stringsAsFactors = FALSE
    )
    results <- c(results, list(df2))
  }

  if (length(results) == 0) return(NULL)
  bind_rows(results)
}

# ── Main extraction loop ─────────────────────────────────────────────────────

all_results <- list()

for (yr in START_YEAR:END_YEAR) {

  message("\n=== Year ", yr, " ===")
  query <- build_query(yr)

  # Count matching records
  search_res <- entrez_search(db = "pubmed", term = query, retmax = 0)
  total      <- search_res$count
  message("  Total records: ", total)

  if (total == 0) next

  # Use web_history for large result sets (avoids URL-length limits)
  search_hist <- entrez_search(
    db      = "pubmed",
    term    = query,
    retmax  = 0,
    use_history = TRUE
  )

  n_batches <- ceiling(total / BATCH_SIZE)
  year_rows <- list()

  for (b in seq_len(n_batches)) {
    retstart <- (b - 1) * BATCH_SIZE
    message("  Batch ", b, " / ", n_batches,
            "  (records ", retstart + 1, "-",
            min(retstart + BATCH_SIZE, total), ")")

    Sys.sleep(SLEEP_SEC)

    xml_raw <- tryCatch(
      entrez_fetch(
        db          = "pubmed",
        web_history = search_hist$web_history,
        rettype     = "xml",
        retmode     = "xml",
        retstart    = retstart,
        retmax      = BATCH_SIZE
      ),
      error = function(e) {
        message("    Fetch error: ", conditionMessage(e), " — retrying once")
        Sys.sleep(5)
        entrez_fetch(
          db          = "pubmed",
          web_history = search_hist$web_history,
          rettype     = "xml",
          retmode     = "xml",
          retstart    = retstart,
          retmax      = BATCH_SIZE
        )
      }
    )

    parsed  <- xmlParse(xml_raw, asText = TRUE)
    articles <- getNodeSet(parsed, "//PubmedArticle")

    for (art in articles) {
      meta <- parse_article(art)
      if (is.null(meta) || nchar(meta$abstract) < 10) next

      cis <- extract_cis(meta$abstract)
      if (is.null(cis) || nrow(cis) == 0) next

      cis$pubmed  <- meta$pmid
      cis$journal <- meta$journal
      cis$Year    <- meta$year
      cis$source  <- "abstract"
      cis$mistake <- (cis$mean < cis$lower) | (cis$mean > cis$upper)

      year_rows <- c(year_rows, list(cis))
    }

    free(parsed)   # release XML memory
  }

  if (length(year_rows) > 0) {
    yr_df <- bind_rows(year_rows)
    all_results <- c(all_results, list(yr_df))
    message("  CIs extracted this year: ", nrow(yr_df))
  }

  # Save a checkpoint after each year
  checkpoint <- bind_rows(all_results)
  save(checkpoint, file = file.path(OUT_DIR, "abstracts_2019_2025_checkpoint.RData"))
}

# ── Assemble final data frame ────────────────────────────────────────────────

abstracts_new <- bind_rows(all_results) %>%
  select(pubmed, journal, lower, mean, upper, Year, source, mistake) %>%
  distinct()   # remove any duplicates from overlapping batches

message("\nTotal CI rows extracted: ", nrow(abstracts_new))
message("Years represented: ", paste(sort(unique(abstracts_new$Year)), collapse = ", "))
message("Mistake rate: ", round(mean(abstracts_new$mistake, na.rm = TRUE) * 100, 2), "%")

# ── Save outputs ─────────────────────────────────────────────────────────────

save(abstracts_new,
     file = file.path(OUT_DIR, "abstracts_2019_2025.RData"))
message("Saved: ", file.path(OUT_DIR, "abstracts_2019_2025.RData"))

# Preview file (first 30 rows, tab-delimited) — mirrors Georgescu.Wren.30.txt
write.table(
  head(abstracts_new, 30),
  file      = file.path(OUT_DIR, "abstracts_2019_2025_30.txt"),
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)
message("Saved preview: ", file.path(OUT_DIR, "abstracts_2019_2025_30.txt"))
