## ---------------------------------------------------------------------------
## build_plos_submission.R
##
## Assembles the PLOS ONE submission package in manuscript/submission_plos/:
##   * renders manuscript.qmd (author details included, double-spaced,
##     continuous line and page numbers) and repairs it for Microsoft Word
##   * renders the supporting information: S1 File, S1 Checklist, S2 Checklist
##   * converts the cover letter
##   * leaves Fig1.tif ... Fig7.tif in place (written by make_plos_figures.R)
##   * writes plos_submission_summary.txt with the journal's own checks
##
## Run after run_all.R and make_plos_figures.R:
##   Rscript manuscript/build_plos_submission.R
## ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(xml2))

PROJ <- Sys.getenv("U5M_PROJ",
                   "D:/khoj_nso/NP-Under-Five-Mortality/u5m_ml_spatiotemporal")
Sys.setenv(U5M_PROJ = PROJ)
MS  <- file.path(PROJ, "manuscript")
OUT <- file.path(MS, "submission_plos")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

MANUSCRIPT <- "Manuscript.docx"      # PLOS ONE is single-blind: authors named

QUARTO <- Sys.which("quarto")
if (QUARTO == "") QUARTO <- "C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe"

run <- function(args) {
  out <- suppressWarnings(system2(QUARTO, args, stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  list(ok = is.null(status) || status == 0, log = out)
}

## ---- render ----------------------------------------------------------------
## U5M_PLOS=1 renders the upload version: captions in the text, figures uploaded
## separately, as PLOS requires.
Sys.setenv(U5M_PLOS = "1")
logs <- list()
for (q in c("manuscript.qmd", "additional_file_1_supplementary.qmd")) {
  message("rendering ", q)
  r <- run(c("render", shQuote(file.path(MS, q)), "--to", "docx"))
  logs[[q]] <- r$log
  if (!r$ok) { cat(r$log, sep = "\n"); stop("render failed: ", q) }
}
file.copy(file.path(MS, "manuscript.docx"), file.path(OUT, MANUSCRIPT), overwrite = TRUE)
file.copy(file.path(MS, "additional_file_1_supplementary.docx"),
          file.path(OUT, "S1_File.docx"), overwrite = TRUE)

## Quarto's table markup breaks a few Word schema rules, which makes Microsoft
## Word refuse to open the file even though pandoc and LibreOffice read it.
source(file.path(PROJ, "review", "repair_docx.R"))
for (f in c(MANUSCRIPT, "S1_File.docx"))
  message("repaired ", repair_docx(file.path(OUT, f)), " paragraphs in ", f)

## ---- reading copy: the same manuscript with the figures in the text --------
## manuscript.qmd places each figure under its caption unless U5M_PLOS is set,
## so the reading copy is the same file rendered again with that switch off.
## It stays outside submission_plos/ so it cannot be uploaded by mistake.
message("rendering manuscript.qmd (reading copy, figures in the text)")
Sys.unsetenv("U5M_PLOS")
r <- run(c("render", shQuote(file.path(MS, "manuscript.qmd")), "--to", "docx"))
if (!r$ok) { cat(r$log, sep = "\n"); stop("render failed: reading copy") }
reading_docx <- file.path(MS, "Manuscript_with_figures.docx")
file.copy(file.path(MS, "manuscript.docx"), reading_docx, overwrite = TRUE)
message("repaired ", repair_docx(reading_docx), " paragraphs in Manuscript_with_figures.docx")

md_docs <- c(additional_file_2_STROBE.md   = "S1_Checklist.docx",
             additional_file_3_TRIPOD_AI.md = "S2_Checklist.docx",
             cover_letter_plos.md           = "Cover_letter.docx")
for (md in names(md_docs)) {
  r <- run(c("pandoc", shQuote(file.path(MS, md)), "-o", shQuote(file.path(OUT, md_docs[[md]]))))
  if (!r$ok) { cat(r$log, sep = "\n"); stop("conversion failed: ", md) }
}

## ---- read the rendered manuscript ------------------------------------------
ns  <- c(w = "http://schemas.openxmlformats.org/wordprocessingml/2006/main")
doc <- read_xml(unz(file.path(OUT, MANUSCRIPT), "word/document.xml"))
paras <- xml_find_all(doc, "/w:document/w:body/w:p", ns)      # excludes table cells
style <- vapply(paras, function(p) {
  s <- xml_find_first(p, "./w:pPr/w:pStyle", ns)
  if (inherits(s, "xml_missing")) "" else xml_attr(s, "val")
}, character(1))
text <- vapply(paras, function(p) paste(xml_text(xml_find_all(p, ".//w:t", ns)), collapse = ""),
               character(1))
h1 <- which(style == "Heading1")
section <- rep(NA_character_, length(paras))
for (i in seq_along(h1)) {
  end <- if (i < length(h1)) h1[i + 1] - 1 else length(paras)
  section[h1[i]:end] <- text[h1[i]]
}
nwords <- function(x) sum(lengths(regmatches(x, gregexpr("\\S+", x, perl = TRUE))))
body_para <- !(seq_along(paras) %in% h1) & !grepl("Caption", style)
in_sec <- function(pattern) !is.na(section) & grepl(pattern, section)

title <- sub("^\\s*title:\\s*\"(.*)\"\\s*$", "\\1",
             grep("^title:", readLines(file.path(MS, "manuscript.qmd"), encoding = "UTF-8"),
                  value = TRUE)[1])
title <- gsub("--", "\u2013", title, fixed = TRUE)
short_title <- sub(".*Short title:\\*{0,2}\\s*", "", grep("Short title", text, value = TRUE)[1])

abstract_words <- nwords(text[in_sec("Abstract") & body_para])
main_words <- nwords(text[in_sec("Introduction|Materials and methods|Results|Discussion|Conclusions") & body_para])
n_tables <- length(xml_find_all(doc, "/w:document/w:body/w:tbl", ns))
n_refs   <- sum(style == "Bibliography")
all_text <- paste(text, collapse = " ")
has_lines <- length(xml_find_all(doc, "//w:sectPr/w:lnNumType", ns)) > 0

## PLOS cites figures as "Fig 1", in ascending order of first citation
fig_hits <- as.integer(unlist(regmatches(all_text, gregexpr("(?<=\\bFig )[0-9]+", all_text, perl = TRUE))))
first_cited <- fig_hits[!duplicated(fig_hits)]
figs_ok <- identical(first_cited, seq_len(7))
captions_ok <- vapply(1:7, function(i) any(grepl(sprintf("^Fig %d\\. ", i), text)), logical(1))

## supporting information: cited in the text and captioned at the end
si_items <- c("S1 File", "S1 Checklist", "S2 Checklist")
si_cited <- vapply(si_items, function(s)
  length(gregexpr(s, all_text, fixed = TRUE)[[1]][gregexpr(s, all_text, fixed = TRUE)[[1]] > 0]),
  integer(1))

## author details must be present (PLOS ONE is single-blind)
id_ok <- all(vapply(c("Acharya", "acharyajibesh", "0000-0001-7261-7729"),
                    function(s) grepl(s, all_text, fixed = TRUE), logical(1)))

## wording left over from the two previous journals
leftover <- unique(unlist(regmatches(
  paste(text[style != "Bibliography"], collapse = " "),
  gregexpr("Additional file|Supplementary material|Highlights|Declaration of generative|BMC Public Health|Elsevier|Fig\\.|Figure [0-9]",
           paste(text[style != "Bibliography"], collapse = " "), perl = TRUE))))

unresolved <- unique(unlist(lapply(logs, function(l) grep("not found|citeproc", l, value = TRUE))))
na_hits <- regmatches(all_text, gregexpr("\\bNA\\b|\\bNaN\\b|\\bInf\\b", all_text))[[1]]
empty_hits <- regmatches(all_text, gregexpr("\\(\\s*\\)|\\(\\s+\\w|,\\s+%[;,) ]|\\(\\s*per 1,000",
                                            all_text, perl = TRUE))[[1]]

## ---- figures ---------------------------------------------------------------
tifs <- sprintf("Fig%d.tif", 1:7)
missing_tif <- tifs[!file.exists(file.path(OUT, tifs))]
tif_mb <- if (!length(missing_tif)) round(file.size(file.path(OUT, tifs)) / 1e6, 2) else NA

flag <- function(bad) if (bad) "  <-- CHECK" else ""
summary <- c(
  "PLOS ONE submission package",
  paste("Built:", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "",
  sprintf("Title characters (<= 250):           %d%s", nchar(title), flag(nchar(title) > 250)),
  sprintf("Short title characters (<= 100):     %d%s", nchar(short_title), flag(nchar(short_title) > 100)),
  sprintf("Abstract words (<= 300):             %d%s", abstract_words, flag(abstract_words > 300)),
  sprintf("Main text words:                     %d", main_words),
  sprintf("Tables in manuscript:                %d", n_tables),
  sprintf("Figures cited as 'Fig n', in order:  %s%s", paste(first_cited, collapse = ", "), flag(!figs_ok)),
  sprintf("Figure captions Fig 1-7 present:     %d/7%s", sum(captions_ok), flag(!all(captions_ok))),
  sprintf("Figure files Fig1-7.tif:             %s%s",
          if (length(missing_tif)) paste("missing:", paste(missing_tif, collapse = ", "))
          else paste0("all present (", paste(tif_mb, collapse = ", "), " MB)"),
          flag(length(missing_tif) > 0)),
  sprintf("Supporting information mentions:     %s%s",
          paste(sprintf("%s x%d", names(si_cited), si_cited), collapse = "; "),
          flag(any(si_cited < 2))),
  sprintf("References:                          %d", n_refs),
  sprintf("Line numbering in manuscript:        %s%s", has_lines, flag(!has_lines)),
  sprintf("Author details present:              %s%s", id_ok, flag(!id_ok)),
  sprintf("Unresolved citation warnings:        %d", length(unresolved)),
  if (length(unresolved)) paste("   ", unresolved) else NULL,
  sprintf("NA/NaN/Inf in rendered text:         %d%s", length(na_hits), flag(length(na_hits) > 0)),
  sprintf("Empty inline values in text:         %d%s", length(empty_hits),
          if (length(empty_hits)) paste0("  <-- CHECK: ", paste(unique(empty_hits), collapse = " | ")) else ""),
  sprintf("Wording from earlier journals:       %d%s", length(leftover),
          if (length(leftover)) paste0("  <-- CHECK: ", paste(leftover, collapse = " | ")) else ""),
  "",
  "Upload as (PLOS ONE file types):",
  "  Cover Letter             Cover_letter.docx",
  "  Manuscript               Manuscript.docx        (with author details, line numbers)",
  "  Figure                   Fig1.tif ... Fig7.tif  (one per upload, in order)",
  "  Supporting Information   S1_File.docx, S1_Checklist.docx, S2_Checklist.docx",
  "",
  "Paste into the submission form (see plos_form_statements.md):",
  "  Financial Disclosure, Competing Interests, Data Availability, Ethics",
  "",
  "Files:",
  paste("  ", sort(list.files(OUT, pattern = "\\.(docx|tif)$")))
)
writeLines(summary, file.path(MS, "plos_submission_summary.txt"))
cat(summary, sep = "\n")
