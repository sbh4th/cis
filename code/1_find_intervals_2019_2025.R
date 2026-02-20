# 1_find_intervals_2019_2025.R
#
# Extract confidence intervals from PubMed abstracts, 2019–2025.
# Updated version of 1_find.intervals.R from agbarnett/intervals.
#
# This script takes the raw abstract data produced by
# 0_find_abstracts_2019_2025.R (or any data frame with a text column)
# and applies an improved regex pipeline to extract CI triplets.
#
# It can also be run as a STANDALONE script: set STANDALONE <- TRUE
# below and it will fetch a random sample of ~5,000 abstracts per year
# directly from PubMed (useful for quick validation before the full run).
#
# Output columns match the original Georgescu.Wren.RData exactly:
#   pubmed   — PubMed ID (character)
#   journal  — journal name
#   lower    — lower confidence limit (numeric)
#   mean     — point estimate (numeric)
#   upper    — upper confidence limit (numeric)
#   Year     — year published (integer)
#   source   — "abstract" (this script) or "full-text" (not in scope here)
#   mistake  — logical: TRUE when mean is outside [lower, upper]
#
# ---------------------------------------------------------------------------
# Usage
# ---------------------------------------------------------------------------
#   # Option 1 – use pre-fetched abstract data from 0_find_abstracts_2019_2025.R
#   source("0_find_abstracts_2019_2025.R")   # creates abstracts_new
#   source("1_find_intervals_2019_2025.R")
#
#   # Option 2 – standalone, fetches its own small sample
#   STANDALONE <- TRUE
#   source("1_find_intervals_2019_2025.R")
# ---------------------------------------------------------------------------

library(dplyr)
library(stringr)
library(rentrez)   # only needed in STANDALONE mode
library(XML)       # only needed in STANDALONE mode

STANDALONE <- FALSE   # set TRUE to fetch abstracts within this script
OUT_DIR    <- "data"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# ============================================================================
# PART 1 – Regex engine for CI extraction
# ============================================================================

# Number pattern (handles negatives, decimals, scientific notation)
.NUM <- "-?[0-9]+(?:\\.[0-9]+)?(?:[eE][+-]?[0-9]+)?"

# Confidence level variants: 95%, 90%, 99%, or unlabelled
.CONF_LEVEL <- "(?:9[0-9]%\\s*)?"

# Opening / closing bracket pairs (round or square)
.OPEN  <- "[\\(\\[]"
.CLOSE <- "[\\)\\]]"

# Build the master regex pattern list
# Each pattern is a named list: pattern (string), group order (mean/lower/upper)

CI_PATTERNS <- list(

  # ── Pattern 1 ──────────────────────────────────────────────────────────────
  # "1.23 (95% CI: 0.98 to 1.56)"  or  "1.23 (95% CI 0.98–1.56)"
  # Most common format in medical abstracts
  list(
    name    = "mean_CI_lower_sep_upper",
    pattern = paste0(
      "(", .NUM, ")",                                # [1] mean
      "\\s*", .OPEN,
      "\\s*", .CONF_LEVEL, "CI[:\\s]+",
      "(", .NUM, ")",                                # [2] lower
      "\\s*(?:to|-|–|,)\\s*",
      "(", .NUM, ")",                                # [3] upper
      "\\s*", .CLOSE
    ),
    groups = c(mean = 1, lower = 2, upper = 3)
  ),

  # ── Pattern 2 ──────────────────────────────────────────────────────────────
  # "95% CI: 0.98 to 1.56" — no explicit mean; mean set to NA
  # Still useful for counting/flagging intervals even without a point estimate
  list(
    name    = "CI_only_lower_upper",
    pattern = paste0(
      .CONF_LEVEL, "CI[:\\s]+",
      "(", .NUM, ")",                                # [1] lower
      "\\s*(?:to|-|–|,)\\s*",
      "(", .NUM, ")"                                 # [2] upper
    ),
    groups = c(mean = NA, lower = 1, upper = 2)
  ),

  # ── Pattern 3 ──────────────────────────────────────────────────────────────
  # "hazard ratio 1.23 (0.98 to 1.56)"  — no explicit CI label
  # Preceded by ratio/effect keywords
  list(
    name    = "keyword_mean_lower_to_upper",
    pattern = paste0(
      "(?:ratio|OR|HR|RR|IRR|coefficient|effect|estimate|SMD|MD|WMD|d)\\s+",
      "(", .NUM, ")",                                # [1] mean
      "\\s*", .OPEN,
      "(", .NUM, ")",                                # [2] lower
      "\\s+to\\s+",
      "(", .NUM, ")",                                # [3] upper
      "\\s*", .CLOSE
    ),
    groups = c(mean = 1, lower = 2, upper = 3)
  ),

  # ── Pattern 4 ──────────────────────────────────────────────────────────────
  # "mean 1.23, 95% CI 0.98-1.56" — comma-separated, no brackets
  list(
    name    = "mean_comma_CI_lower_upper",
    pattern = paste0(
      "(", .NUM, ")",                                # [1] mean
      "\\s*,\\s*",
      .CONF_LEVEL, "CI[:\\s]+",
      "(", .NUM, ")",                                # [2] lower
      "\\s*(?:to|-|–|,)\\s*",
      "(", .NUM, ")"                                 # [3] upper
    ),
    groups = c(mean = 1, lower = 2, upper = 3)
  ),

  # ── Pattern 5 ──────────────────────────────────────────────────────────────
  # "CI [0.98, 1.56]" — square-bracket style with comma
  list(
    name    = "CI_square_bracket",
    pattern = paste0(
      .CONF_LEVEL, "CI\\s*\\[",
      "(", .NUM, ")",                                # [1] lower
      "\\s*,\\s*",
      "(", .NUM, ")",                                # [2] upper
      "\\s*\\]"
    ),
    groups = c(mean = NA, lower = 1, upper = 2)
  )
)

#' Extract all CI triplets from a single abstract text string.
#'
#' @param text  Character string (one abstract).
#' @return      data.frame with columns: mean, lower, upper, pattern_name
#'              or NULL if nothing found.
extract_cis_from_text <- function(text) {
  if (is.na(text) || nchar(trimws(text)) == 0) return(NULL)

  # Normalise unicode dashes and the word "to" used as a range separator
  text <- str_replace_all(text, "\u2013|\u2014|\u2212|\uFE58|\uFE63|\uFF0D", "-")

  results <- list()

  for (pat in CI_PATTERNS) {
    m <- str_match_all(text, regex(pat$pattern, ignore_case = TRUE))[[1]]
    if (nrow(m) == 0) next

    g <- pat$groups
    df <- data.frame(
      mean  = if (is.na(g["mean"])) NA_real_
              else suppressWarnings(as.numeric(m[, g["mean"] + 1])),
      lower = suppressWarnings(as.numeric(m[, g["lower"] + 1])),
      upper = suppressWarnings(as.numeric(m[, g["upper"] + 1])),
      pattern_name = pat$name,
      stringsAsFactors = FALSE
    )
    results <- c(results, list(df))
  }

  if (length(results) == 0) return(NULL)

  out <- bind_rows(results)

  # Sanity filters
  out <- out %>%
    filter(
      !is.na(lower), !is.na(upper),
      lower < upper,                    # lower must be less than upper
      abs(lower) < 1e6,                 # exclude absurd values
      abs(upper) < 1e6
    )

  if (nrow(out) == 0) return(NULL)
  out
}

# ============================================================================
# PART 2 – Build the output data frame
# ============================================================================

#' Process a data frame of abstracts and return CI rows.
#'
#' @param df  data.frame with columns: pubmed, journal, Year, abstract
#' @return    data.frame matching the Georgescu.Wren.RData schema
process_abstracts <- function(df) {
  out_list <- vector("list", nrow(df))

  for (i in seq_len(nrow(df))) {
    if (i %% 10000 == 0) message("  Processing row ", i, " / ", nrow(df))

    cis <- extract_cis_from_text(df$abstract[i])
    if (is.null(cis)) next

    cis$pubmed   <- df$pubmed[i]
    cis$journal  <- df$journal[i]
    cis$Year     <- df$Year[i]
    cis$source   <- "abstract"
    cis$mistake  <- !is.na(cis$mean) &
                    ((cis$mean < cis$lower) | (cis$mean > cis$upper))

    out_list[[i]] <- cis
  }

  bind_rows(out_list) %>%
    select(pubmed, journal, lower, mean, upper, Year, source, mistake)
}

# ============================================================================
# PART 3 – Standalone fetch (optional)
# ============================================================================

if (STANDALONE) {
  message("Running in STANDALONE mode — fetching sample abstracts from PubMed")

  SAMPLE_PER_YEAR <- 5000
  START_YEAR      <- 2019
  END_YEAR        <- 2025
  BATCH_SIZE      <- 200
  SLEEP_SEC       <- 0.34

  build_query <- function(year) {
    paste0(
      '("', year, '/01/01"[PDAT] : "', year, '/12/31"[PDAT]) ',
      'AND "journal article"[PT] AND hasabstract AND English[LA] ',
      'AND ("confidence interval"[TIAB] OR "95% CI"[TIAB])'
    )
  }

  parse_article <- function(node) {
    tryCatch({
      pmid    <- xmlValue(node[["MedlineCitation"]][["PMID"]])
      journal <- xmlValue(
        node[["MedlineCitation"]][["Article"]][["Journal"]][["Title"]]
      )
      year_node <- node[["MedlineCitation"]][["Article"]][["Journal"]][
        ["JournalIssue"]][["PubDate"]][["Year"]]
      pub_year  <- if (!is.null(year_node)) as.integer(xmlValue(year_node))
                   else NA_integer_
      abs_nodes <- getNodeSet(node, ".//AbstractText")
      abstract  <- paste(sapply(abs_nodes, xmlValue), collapse = " ")
      list(pubmed = pmid, journal = journal, Year = pub_year, abstract = abstract)
    }, error = function(e) NULL)
  }

  raw_abstracts <- list()

  for (yr in START_YEAR:END_YEAR) {
    message("Fetching year ", yr)
    sh <- entrez_search(db = "pubmed", term = build_query(yr),
                        retmax = 0, use_history = TRUE)
    to_fetch <- min(sh$count, SAMPLE_PER_YEAR)
    n_batches <- ceiling(to_fetch / BATCH_SIZE)

    for (b in seq_len(n_batches)) {
      Sys.sleep(SLEEP_SEC)
      retstart <- (b - 1) * BATCH_SIZE
      retmax_b <- min(BATCH_SIZE, to_fetch - retstart)
      if (retmax_b <= 0) break

      xml_raw <- entrez_fetch(
        db = "pubmed", web_history = sh$web_history,
        rettype = "xml", retmode = "xml",
        retstart = retstart, retmax = retmax_b
      )
      parsed   <- xmlParse(xml_raw, asText = TRUE)
      articles <- getNodeSet(parsed, "//PubmedArticle")
      arts     <- lapply(articles, parse_article)
      raw_abstracts <- c(raw_abstracts, Filter(Negate(is.null), arts))
      free(parsed)
    }
  }

  abstracts_new <- bind_rows(lapply(raw_abstracts, as.data.frame,
                                    stringsAsFactors = FALSE))
}

# ============================================================================
# PART 4 – Run extraction and save
# ============================================================================

if (!exists("abstracts_new")) {
  stop(paste(
    "abstracts_new not found.\n",
    "Either run 0_find_abstracts_2019_2025.R first,",
    "or set STANDALONE <- TRUE at the top of this script."
  ))
}

message("Extracting CIs from ", nrow(abstracts_new), " abstracts...")
intervals_new <- process_abstracts(abstracts_new)

message("\n── Summary ──────────────────────────────────────────────")
message("Total CI rows:  ", nrow(intervals_new))
message("Unique PMIDs:   ", n_distinct(intervals_new$pubmed))
message("Years covered:  ", paste(sort(unique(intervals_new$Year)), collapse = ", "))
message("Mistake rate:   ",
        round(mean(intervals_new$mistake, na.rm = TRUE) * 100, 2), "%")
message("─────────────────────────────────────────────────────────")

# Save in the same format as the original Georgescu.Wren.RData
save(intervals_new,
     file = file.path(OUT_DIR, "intervals_2019_2025.RData"))
message("Saved: ", file.path(OUT_DIR, "intervals_2019_2025.RData"))

# Tab-delimited preview (first 30 rows)
write.table(
  head(intervals_new, 30),
  file      = file.path(OUT_DIR, "intervals_2019_2025_30.txt"),
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)
message("Saved preview: ", file.path(OUT_DIR, "intervals_2019_2025_30.txt"))

# ============================================================================
# PART 5 – (Optional) Combine with the original Georgescu.Wren dataset
# ============================================================================

# Uncomment the block below once you have the original .RData file to merge
# the two datasets for a unified 1976–2025 analysis.

# if (file.exists(file.path(OUT_DIR, "Georgescu.Wren.RData"))) {
#   load(file.path(OUT_DIR, "Georgescu.Wren.RData"))   # loads object `data`
#
#   # Standardise column names if needed
#   if ("source" %in% names(data) && !"source" %in% names(intervals_new)) {
#     intervals_new$source <- "abstract"
#   }
#
#   combined <- bind_rows(
#     data %>% mutate(dataset = "Georgescu_Wren"),
#     intervals_new %>% mutate(dataset = "Barnett_2019_2025")
#   ) %>% distinct(pubmed, lower, mean, upper, .keep_all = TRUE)
#
#   save(combined, file = file.path(OUT_DIR, "intervals_combined.RData"))
#   message("Combined dataset saved: ", nrow(combined), " rows")
# }
