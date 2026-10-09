## ---------------------------------------------------------------------------
## make_figure_file.R
##
## Writes submission_plos/All_figures.docx: the seven main figures, each under
## the caption exactly as it appears in the manuscript. Convenient for sending
## the figures to a co-author or reviewer in one file.
##
## PLOS wants the figures uploaded one TIFF at a time, so this file is NOT part
## of the submission.
##
##   Rscript manuscript/make_figure_file.R
## ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(xml2))

PROJ <- Sys.getenv("U5M_PROJ",
                   "D:/khoj_nso/NP-Under-Five-Mortality/u5m_ml_spatiotemporal")
MS  <- file.path(PROJ, "manuscript")
OUT <- file.path(MS, "submission_plos")
PNG <- file.path(PROJ, "figures", "plos")

## figure n -> the file behind it, in citation order
SRC <- c("fig1_national_trend", "fig3_clusters", "fig2_smoothed_u5mr_maps",
         "fig4_bayesian_maps", "fig5_ml_performance", "fig6_shap",
         "fig7_shap_by_province")

## captions are taken from the rendered manuscript, so the inline values in them
## are the ones that were actually published
ns  <- c(w = "http://schemas.openxmlformats.org/wordprocessingml/2006/main")
doc <- read_xml(unz(file.path(OUT, "Manuscript.docx"), "word/document.xml"))
txt <- vapply(xml_find_all(doc, "/w:document/w:body/w:p", ns),
              function(p) paste(xml_text(xml_find_all(p, ".//w:t", ns)), collapse = ""),
              character(1))
caps <- vapply(seq_along(SRC), function(i) {
  hit <- grep(sprintf("^Fig %d\\. ", i), txt, value = TRUE)
  if (!length(hit)) stop("caption not found for Fig ", i)
  hit[1]
}, character(1))

## bold the label and the title sentence, as in the manuscript
bold_head <- function(cap) {
  i <- regexpr("\\.\\s", substring(cap, 8))          # end of the title sentence
  head <- substr(cap, 1, 7 + i)
  sub(paste0("^", gsub("([][(){}.*+?^$|\\\\])", "\\\\\\1", head)),
      paste0("**", trimws(head), "** "), cap)
}

md <- c("# Figures",
        "",
        paste("Main-text figures of the manuscript, in citation order.",
              "Captions are reproduced from the manuscript.",
              "This file is for circulation only: PLOS ONE takes the figures as",
              "individual TIFF uploads (Fig1.tif ... Fig7.tif)."),
        "")
for (i in seq_along(SRC)) {
  img <- file.path(PNG, paste0(SRC[i], ".png"))
  if (!file.exists(img)) stop("missing figure: ", img, " (run make_plos_figures.R)")
  md <- c(md, bold_head(caps[i]), "",
          sprintf("![](../figures/plos/%s.png){width=6.3in}", SRC[i]), "")
}

md_file <- file.path(MS, "all_figures.md")
writeLines(md, md_file, useBytes = TRUE)

QUARTO <- Sys.which("quarto")
if (QUARTO == "") QUARTO <- "C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe"
out <- file.path(OUT, "All_figures.docx")
st <- system2(QUARTO, c("pandoc", shQuote(md_file), "--resource-path", shQuote(MS),
                        "--reference-doc", shQuote(file.path(MS, "reference.docx")),
                        "-o", shQuote(out)), stdout = TRUE, stderr = TRUE)
if (!is.null(attr(st, "status")) && attr(st, "status") != 0) { cat(st, sep = "\n"); stop("pandoc failed") }
unlink(md_file)

source(file.path(PROJ, "review", "repair_docx.R"))
repair_docx(out)
message("wrote ", out, "  (", round(file.size(out) / 1e6, 2), " MB, ", length(SRC), " figures)")
