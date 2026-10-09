## ---------------------------------------------------------------------------
## make_reference_docx.R
##
## Builds manuscript/reference.docx, the Word template Quarto applies when it
## renders manuscript.qmd to .docx. It gives the submission file the layout
## reviewers expect: 12 pt Times New Roman, double line spacing, continuous
## line numbers, centred page numbers, A4 with 2.5 cm margins. Tables keep
## single spacing so they stay legible.
##
## Run once (or after changing the layout):
##   Rscript manuscript/make_reference_docx.R
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(xml2); library(zip) })

MS     <- "D:/khoj_nso/NP-Under-Five-Mortality/u5m_ml_spatiotemporal/manuscript"
QUARTO <- Sys.which("quarto")
if (QUARTO == "") QUARTO <- "C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe"

work <- file.path(tempdir(), "refdoc")
unlink(work, recursive = TRUE); dir.create(work, recursive = TRUE)
base <- file.path(tempdir(), "pandoc_default_reference.docx")
unlink(base)
system2(QUARTO, c("pandoc", "-o", shQuote(base),
                  "--print-default-data-file", "reference.docx"))
stopifnot(file.exists(base))
unzip(base, exdir = work)

W  <- "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
RR <- "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
ns <- c(w = W, r = RR)

child <- function(parent, name) {
  node <- xml_find_first(parent, paste0("./", name), ns)
  if (inherits(node, "xml_missing")) node <- xml_add_child(parent, name)
  node
}

## ---- styles ---------------------------------------------------------------
sty_path <- file.path(work, "word", "styles.xml")
sty <- read_xml(sty_path)

set_font <- function(rpr, size_half_pts = NULL) {
  f <- child(rpr, "w:rFonts")
  for (a in c("w:asciiTheme", "w:hAnsiTheme", "w:eastAsiaTheme", "w:cstheme"))
    xml_attr(f, a) <- NULL
  for (a in c("w:ascii", "w:hAnsi", "w:eastAsia", "w:cs"))
    xml_attr(f, a) <- "Times New Roman"
  if (!is.null(size_half_pts)) {
    sz   <- child(rpr, "w:sz");   xml_attr(sz,   "w:val") <- size_half_pts
    szcs <- child(rpr, "w:szCs"); xml_attr(szcs, "w:val") <- size_half_pts
  }
}

dd  <- child(xml_find_first(sty, "//w:styles", ns), "w:docDefaults")
set_font(child(child(dd, "w:rPrDefault"), "w:rPr"), "24")
sp_def <- child(child(child(dd, "w:pPrDefault"), "w:pPr"), "w:spacing")
xml_attr(sp_def, "w:line") <- "480"
xml_attr(sp_def, "w:lineRule") <- "auto"

## headings and title: same face, black, no theme colour
for (st in xml_find_all(sty, "//w:style[starts-with(@w:styleId,'Heading') or @w:styleId='Title' or @w:styleId='Subtitle']", ns)) {
  rpr <- child(st, "w:rPr")
  set_font(rpr)
  col <- xml_find_first(rpr, "./w:color", ns)
  if (!inherits(col, "xml_missing")) xml_remove(col)
}

## tables (pandoc's "Compact" paragraph style) stay single-spaced, and are set
## two points smaller than the text: at 12 pt double-spaced a wide table does
## not fit the type area and Word wraps figures such as "15,094" mid-number
compact <- xml_find_first(sty, "//w:style[@w:styleId='Compact']", ns)
if (!inherits(compact, "xml_missing")) {
  cs <- child(child(compact, "w:pPr"), "w:spacing")
  xml_attr(cs, "w:line") <- "240"; xml_attr(cs, "w:lineRule") <- "auto"
  xml_attr(cs, "w:before") <- "20"; xml_attr(cs, "w:after") <- "20"
  set_font(child(compact, "w:rPr"), "20")          # 10 pt
}
write_xml(sty, sty_path)

## ---- footer with a PAGE field ----------------------------------------------
footer_xml <- sprintf(
'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:ftr xmlns:w="%s"><w:p><w:pPr><w:jc w:val="center"/></w:pPr><w:r><w:fldChar w:fldCharType="begin"/></w:r><w:r><w:instrText xml:space="preserve"> PAGE </w:instrText></w:r><w:r><w:fldChar w:fldCharType="separate"/></w:r><w:r><w:t>1</w:t></w:r><w:r><w:fldChar w:fldCharType="end"/></w:r></w:p></w:ftr>', W)
writeLines(footer_xml, file.path(work, "word", "footerU5M.xml"), useBytes = TRUE)

rels_path <- file.path(work, "word", "_rels", "document.xml.rels")
rels <- read_xml(rels_path)
if (length(xml_find_all(rels, "//*[@Id='rIdU5MFooter']")) == 0) {
  r <- xml_add_child(rels, "Relationship")
  xml_attr(r, "Id") <- "rIdU5MFooter"
  xml_attr(r, "Type") <- "http://schemas.openxmlformats.org/officeDocument/2006/relationships/footer"
  xml_attr(r, "Target") <- "footerU5M.xml"
}
write_xml(rels, rels_path)

ct_path <- file.path(work, "[Content_Types].xml")
ct <- read_xml(ct_path)
if (length(xml_find_all(ct, "//*[@PartName='/word/footerU5M.xml']")) == 0) {
  o <- xml_add_child(ct, "Override")
  xml_attr(o, "PartName") <- "/word/footerU5M.xml"
  xml_attr(o, "ContentType") <- "application/vnd.openxmlformats-officedocument.wordprocessingml.footer+xml"
}
## image types: pandoc's default template declares none, and Word refuses to
## open a file whose embedded figures have no declared content type
img_types <- c(png = "image/png", jpeg = "image/jpeg", jpg = "image/jpeg",
               gif = "image/gif", emf = "image/x-emf", svg = "image/svg+xml")
for (ext in names(img_types)) {
  if (length(xml_find_all(ct, sprintf("//*[local-name()='Default' and @Extension='%s']", ext))) == 0) {
    d <- xml_add_child(ct, "Default", .where = 0)
    xml_attr(d, "Extension") <- ext; xml_attr(d, "ContentType") <- img_types[[ext]]
  }
}
write_xml(ct, ct_path)

## ---- section properties, in the order the OOXML schema requires -----------
doc_path <- file.path(work, "word", "document.xml")
doc <- read_xml(doc_path)
new_sect <- read_xml(sprintf(
'<w:sectPr xmlns:w="%s" xmlns:r="%s">
  <w:footerReference w:type="default" r:id="rIdU5MFooter"/>
  <w:pgSz w:w="11906" w:h="16838"/>
  <w:pgMar w:top="1417" w:right="1417" w:bottom="1417" w:left="1417" w:header="708" w:footer="708" w:gutter="0"/>
  <w:lnNumType w:countBy="1" w:distance="360" w:restart="continuous"/>
  <w:pgNumType w:start="1"/>
  <w:cols w:space="708"/>
</w:sectPr>', W, RR))
old_sect <- xml_find_first(doc, "//w:body/w:sectPr", ns)
if (inherits(old_sect, "xml_missing")) {
  xml_add_child(xml_find_first(doc, "//w:body", ns), new_sect)
} else {
  xml_replace(old_sect, new_sect)
}
write_xml(doc, doc_path)

## ---- repackage ------------------------------------------------------------
out <- file.path(MS, "reference.docx")
unlink(out)
zip::zip(out, files = list.files(work, recursive = TRUE, all.files = TRUE),
         root = work)
message("wrote ", out)
