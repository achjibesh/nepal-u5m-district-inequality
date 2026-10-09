## ---------------------------------------------------------------------------
## make_plos_figures.R
##
## Writes the seven main figures as PLOS-compliant TIFFs in
## manuscript/submission_plos/, named Fig1.tif ... Fig7.tif in the order they
## are cited in the manuscript.
##
## PLOS figure rules applied here:
##   * TIFF only, RGB, LZW-compressed, 300 dpi
##   * 789-2250 px wide, at most 2625 px tall  (7.5 x 8.75 in at 300 dpi)
##   * 10 MB or less per file
##   * text in 8-12 pt at the printed size
##   * no figure number or title inside the image (those live in the caption)
##
## The figures in 06_figures.R were drawn for a 9.5-13.5 inch canvas. Simply
## saving them narrower clips the panel titles, and downsampling the wide
## version drops the type below 8 pt. So each figure is rebuilt here from its
## own panels: side-by-side panels are stacked where the width no longer allows
## them, titles are wrapped, and all text is reset to 8-10 pt. A PNG preview is
## written next to each TIFF so the rebuilt layout can be checked.
##
##   Rscript manuscript/make_plos_figures.R
## ---------------------------------------------------------------------------

PROJ <- Sys.getenv("U5M_PROJ", "D:/khoj_nso/NP-Under-Five-Mortality/u5m_ml_spatiotemporal")
Sys.setenv(U5M_PROJ = PROJ)

OUT <- file.path(PROJ, "manuscript", "submission_plos")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
PREV <- file.path(OUT, "previews")   # layout checks, not uploaded
dir.create(PREV, showWarnings = FALSE, recursive = TRUE)
PNG <- file.path(PROJ, "figures", "plos")   # readable PNGs for draft renders
dir.create(PNG, showWarnings = FALSE, recursive = TRUE)

## 06_figures.R builds every panel and composite figure; the objects are what
## this script rebuilds from.
source(file.path(PROJ, "R", "06_figures.R"))
suppressPackageStartupMessages(library(patchwork))

MAX_W_IN <- 7.5        # 2250 px at 300 dpi, the PLOS maximum
MAX_H_IN <- 8.75       # 2625 px at 300 dpi
DPI      <- 300
MAX_MB   <- 10

wrap_lab <- function(x, width) if (is.null(x) || !is.character(x)) x else
  paste(strwrap(x, width = width), collapse = "\n")

## 11 pt body text, 12.5 pt panel titles: the top of the PLOS 8-12 pt range, so
## the figures stay readable when Word shrinks them into a 6.3 inch column
fit_text <- function(p, base = 11, wrap = 52) {
  ## PLOS asks for no figure title inside the image: a whole-figure title
  ## repeats the caption. Panel headings ("A ...", "B ...") stay.
  t0 <- p$labels$title
  if (!is.null(t0) && is.character(t0) && !grepl("^[A-D]\\s", t0)) {
    p$labels$title <- NULL
    p$labels$subtitle <- NULL
  }
  p$labels$title    <- wrap_lab(p$labels$title, wrap)
  p$labels$subtitle <- wrap_lab(p$labels$subtitle, wrap)
  ## captions run under the full figure width, so they wrap wider than titles
  p$labels$caption  <- wrap_lab(p$labels$caption, 95)
  el <- list(
    text          = ggplot2::element_text(size = base),
    plot.title    = ggplot2::element_text(face = "bold", size = base + 1.5),
    plot.subtitle = ggplot2::element_text(colour = "grey30", size = base - 0.5),
    plot.caption  = ggplot2::element_text(colour = "grey30", size = base - 0.5, hjust = 0),
    plot.caption.position = "plot",
    plot.margin   = ggplot2::margin(8, 8, 8, 8),
    legend.title  = ggplot2::element_text(size = base),
    legend.text   = ggplot2::element_text(size = base - 0.5),
    legend.key.size = grid::unit(0.42, "cm"),
    strip.text    = ggplot2::element_text(face = "bold", size = base))
  ## the maps use theme_void: resetting the axis elements there would put the
  ## graticule labels back on
  if (!inherits(p$theme$axis.text, "element_blank"))
    el <- c(el, list(axis.title = ggplot2::element_text(size = base),
                     axis.text  = ggplot2::element_text(size = base - 0.5)))
  p + do.call(ggplot2::theme, el)
}

## Figures in the order they are cited in the manuscript. Panels that sat side
## by side on a 12-13.5 inch canvas are stacked here.
figs <- list(
  Fig1 = list(src = "fig1_national_trend", h = 7.2, plot = function()
    fit_text(p1a) / fit_text(p1b) + plot_layout(heights = c(1, 1.15))),

  Fig2 = list(src = "fig3_clusters", h = 7.5, plot = function()                    # scan and Gi* maps
    fit_text(p3a, wrap = 56) / fit_text(p3b, wrap = 56) + plot_layout(heights = c(1.25, 1))),

  Fig3 = list(src = "fig2_smoothed_u5mr_maps", h = 5.8, plot = function()                    # smoothed district maps
    fit_text(fig2, wrap = 56)),

  Fig4 = list(src = "fig4_bayesian_maps", h = 8.1, plot = function()                    # three maps, now stacked
    fit_text(p4a, wrap = 56) / fit_text(p4b, wrap = 56) / fit_text(p4c, wrap = 56)),

  ## Panels A-C share one colour scale, so a single legend is collected under
  ## the whole figure: four model names do not fit under a half-width panel.
  ## Panel B's duplicate copy (with PR-AUC values, which are in Table 7) is
  ## dropped, and the crowded axis ticks are thinned.
  Fig5 = list(src = "fig5_ml_performance", h = 6.9, plot = function() {
    thin <- ggplot2::scale_x_continuous(breaks = c(0, 0.5, 1))
    ## guides(colour = "none") drops a legend for good; a theme setting would be
    ## overridden by the composite theme applied with & below
    a <- fit_text(p5a, wrap = 34) + thin +
      ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2, byrow = TRUE))
    b <- fit_text(p5b, wrap = 34) + thin + ggplot2::guides(colour = "none")
    cc <- fit_text(p5c, wrap = 34) + ggplot2::guides(colour = "none")
    d <- fit_text(p5d, wrap = 34)
    ((a | b) / (cc | d)) + plot_layout(guides = "collect") &
      ggplot2::theme(legend.position = "bottom", legend.box = "vertical",
                     legend.margin = ggplot2::margin(2, 2, 2, 2))
  }),

  ## SHAP panels, now stacked; the five-group legend of panel A needs two rows
  ## at this text size
  Fig6 = list(src = "fig6_shap", h = 8.4, plot = function()
    (fit_text(p6a, wrap = 56) +
       ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE))) /
      fit_text(p6b, wrap = 56)),

  Fig7 = list(src = "fig7_shap_by_province", h = 5.4, plot = function() fit_text(p7, wrap = 56))
)

report <- data.frame()
for (nm in names(figs)) {
  spec <- figs[[nm]]
  p    <- spec$plot()
  tif  <- file.path(OUT, paste0(nm, ".tif"))
  w    <- MAX_W_IN
  h    <- min(spec$h, MAX_H_IN)

  ## shrink only if LZW cannot bring the file under the 10 MB limit
  for (shrink in c(1, 0.92, 0.85, 0.78)) {
    grDevices::tiff(tif, width = w * shrink, height = h * shrink, units = "in",
                    res = DPI, compression = "lzw", type = "cairo")
    print(p)
    dev.off()
    if (file.size(tif) / 1e6 <= MAX_MB) break
  }
  ggplot2::ggsave(file.path(PREV, paste0(nm, "_preview.png")), p,
                  width = w * shrink, height = h * shrink, dpi = 150, bg = "white")
  ## the same figure as a PNG, which the manuscript embeds in draft renders
  ggplot2::ggsave(file.path(PNG, paste0(spec$src, ".png")), p,
                  width = w * shrink, height = h * shrink, dpi = 300, bg = "white")

  report <- rbind(report, data.frame(
    figure = nm, width_px = round(w * shrink * DPI), height_px = round(h * shrink * DPI),
    MB = round(file.size(tif) / 1e6, 2)))
}

report$ok <- with(report, width_px >= 789 & width_px <= 2250 &
                    height_px <= 2625 & MB <= MAX_MB)
print(report, row.names = FALSE)
if (!all(report$ok)) stop("one or more figures fall outside the PLOS limits")
message("wrote ", nrow(report), " TIFFs to ", OUT)
