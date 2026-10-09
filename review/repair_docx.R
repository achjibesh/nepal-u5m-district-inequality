## ---------------------------------------------------------------------------
## repair_docx.R
##
## repair_docx(path): fixes schema problems in Quarto's Word output that make
## Microsoft Word refuse to open the file ("the file appears to be corrupted"),
## although pandoc and LibreOffice read it:
##   * paragraphs with two <w:pPr> blocks (Quarto caption and callout markup):
##     merged into one, later settings winning, in schema order
##   * <w:tblPr> children out of schema order (callout boxes)
##   * list ids (<w:nsid>, <w:tmpl>) longer than 8 hex digits
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(xml2); library(zip) })

PPR_ORDER <- c("pStyle", "keepNext", "keepLines", "pageBreakBefore", "framePr", "widowControl", "numPr",
               "suppressLineNumbers", "pBdr", "shd", "tabs", "suppressAutoHyphens", "kinsoku", "wordWrap",
               "overflowPunct", "topLinePunct", "autoSpaceDE", "autoSpaceDN", "bidi", "adjustRightInd",
               "snapToGrid", "spacing", "ind", "contextualSpacing", "mirrorIndents", "suppressOverlap", "jc",
               "textDirection", "textAlignment", "textboxTightWrap", "outlineLvl", "divId", "cnfStyle",
               "rPr", "sectPr", "pPrChange")
TBLPR_ORDER <- c("tblStyle", "tblpPr", "tblOverlap", "bidiVisual", "tblStyleRowBandSize",
                 "tblStyleColBandSize", "tblW", "jc", "tblCellSpacing", "tblInd", "tblBorders", "shd",
                 "tblLayout", "tblCellMar", "tblLook", "tblCaption", "tblDescription", "tblPrChange")

reorder_kids <- function(node, ord) {
  kids <- xml_children(node)
  if (length(kids) < 2) return(invisible())
  pos <- match(xml_name(kids), ord); pos[is.na(pos)] <- length(ord) + 1
  if (!is.unsorted(pos)) return(invisible())
  keep <- kids[order(pos)]
  xml_remove(kids)
  for (k in keep) xml_add_child(node, k)
}

repair_docx <- function(path) {
  work <- file.path(tempdir(), paste0("docx_fix_", as.integer(Sys.time())))
  unlink(work, recursive = TRUE); dir.create(work)
  unzip(path, exdir = work)

  dp  <- file.path(work, "word", "document.xml")
  doc <- read_xml(dp)
  bad <- xml_find_all(doc, "//w:p[count(w:pPr) > 1]")
  for (p in bad) {
    pprs <- xml_find_all(p, "./w:pPr"); first <- pprs[[1]]
    for (extra in pprs[-1]) {
      for (k in xml_children(extra)) {
        same <- xml_find_first(first, paste0("./w:", xml_name(k)))
        if (!inherits(same, "xml_missing")) xml_remove(same)
        xml_add_child(first, k)
      }
      xml_remove(extra)
    }
    reorder_kids(first, PPR_ORDER)
  }
  for (tp in xml_find_all(doc, "//w:tblPr")) reorder_kids(tp, TBLPR_ORDER)
  write_xml(doc, dp)

  np_ <- file.path(work, "word", "numbering.xml")
  if (file.exists(np_)) {
    num <- read_xml(np_)
    for (n in xml_find_all(num, "//w:nsid | //w:tmpl")) {
      v <- xml_attr(n, "w:val")
      if (!is.na(v) && nchar(v) > 8) xml_attr(n, "w:val") <- substr(v, nchar(v) - 7, nchar(v))
    }
    write_xml(num, np_)
  }

  unlink(path)
  files <- list.files(work, recursive = TRUE, all.files = TRUE)
  files <- c("[Content_Types].xml", setdiff(files, "[Content_Types].xml"))
  zip::zip(path, files = files, root = work, mode = "mirror")
  unlink(work, recursive = TRUE)
  length(bad)
}
