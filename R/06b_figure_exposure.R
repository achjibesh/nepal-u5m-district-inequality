## ---------------------------------------------------------------------------
## 06b_figure_exposure.R
##
## Methods figure (Additional file 1): how the discrete-time survival design
## turns a child into exposure and events, using real records. Panel A draws
## each child's observed lifetime to exit; panel B shows the same children
## after the person-period expansion.
## ---------------------------------------------------------------------------

source(file.path(Sys.getenv("U5M_PROJ",
                            "D:/khoj_nso/NP-Under-Five-Mortality/u5m_ml_spatiotemporal"),
                 "R", "00_setup.R"))

b  <- readRDS(file.path(DIR_DERIVED, "births.rds"))
pp <- readRDS(file.path(DIR_DERIVED, "person_period.rds"))

b <- b %>% mutate(status = case_when(
  u5_death == 1L   ~ "Died before 5",
  age_at_int >= 60 ~ "Reached age 5 alive",
  TRUE             ~ "Still under 5 at interview"))

## Real children spanning the range of exit times. `pick` selects which of the
## matching records to take, so the three survivors are three different children.
want <- tibble::tribble(
  ~status,                       ~target, ~pick,
  "Died before 5",                 0,      1,
  "Died before 5",                 2,      1,
  "Died before 5",                 4,      1,
  "Died before 5",                11,      1,
  "Died before 5",                31,      1,
  "Still under 5 at interview",    3,      1,
  "Still under 5 at interview",   14,      1,
  "Still under 5 at interview",   27,      1,
  "Still under 5 at interview",   41,      1,
  "Still under 5 at interview",   56,      1,
  "Reached age 5 alive",          60,      1,
  "Reached age 5 alive",          60,     40,
  "Reached age 5 alive",          60,    900
)

sel <- purrr::pmap_dfr(want, function(status, target, pick) {
  b %>% filter(status == !!status,
               abs(exit_month - target) == min(abs(exit_month - target))) %>%
    slice(pick)
}) %>%
  distinct(child_id, .keep_all = TRUE) %>%
  arrange(exit_month) %>%
  mutate(row = row_number(), lab = sprintf("Child %2d", row))

STCOL <- c("Died before 5"              = "#B2182B",
           "Reached age 5 alive"        = "#1B7837",
           "Still under 5 at interview" = "#2166AC")

pA <- ggplot(sel, aes(y = reorder(lab, row))) +
  geom_vline(xintercept = AGE_SEGMENTS$start[-1], colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 60, colour = "grey40", linetype = 2, linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = exit_month, yend = reorder(lab, row), colour = status),
               linewidth = 2.6, lineend = "butt") +
  geom_point(data = sel %>% filter(status == "Died before 5"),
             aes(x = exit_month), shape = 4, size = 3.4, stroke = 1.5, colour = "#B2182B") +
  geom_point(data = sel %>% filter(status == "Reached age 5 alive"),
             aes(x = exit_month), shape = 16, size = 3, colour = "#1B7837") +
  geom_segment(data = sel %>% filter(status == "Still under 5 at interview"),
               aes(x = exit_month, xend = exit_month + 3.2, yend = reorder(lab, row)),
               arrow = arrow(length = unit(0.16, "cm"), type = "open"),
               colour = "#2166AC", linewidth = 0.6) +
  geom_text(aes(x = exit_month, label = sprintf(" %.0f mo", exit_month)),
            hjust = 0, size = 2.9, colour = "grey30", nudge_x = 4.5) +
  scale_colour_manual(values = STCOL, name = NULL) +
  scale_x_continuous(breaks = c(0, AGE_SEGMENTS$end), limits = c(0, 74), expand = c(0, 0)) +
  labs(x = "Age of the child, in months", y = NULL,
       title = "A  Every child contributes exposure; only some contribute an event",
       subtitle = str_wrap(paste(
         "Each bar is one child, from birth until observation ends. Bar length is that",
         "child's exposure. Grey lines mark the eight DHS age segments; the dashed line",
         "is the fifth birthday."), 108)) +
  theme_u5m() +
  theme(axis.text.y = element_text(family = "mono"),
        panel.grid.major.y = element_blank(), legend.position = "bottom")

ppl <- pp %>%
  filter(child_id %in% sel$child_id) %>%
  left_join(sel %>% select(child_id, lab, row), by = "child_id") %>%
  mutate(cell = if_else(died == 1L, "Death occurs here", "At risk, survived"))

grid_all <- tidyr::crossing(sel %>% select(lab, row),
                            seg_f = factor(AGE_SEGMENTS$label, levels = AGE_SEGMENTS$label)) %>%
  left_join(ppl %>% select(lab, seg_f, expo, cell), by = c("lab", "seg_f")) %>%
  mutate(cell = tidyr::replace_na(cell, "Not at risk (not reached)"))

pB <- ggplot(grid_all, aes(seg_f, reorder(lab, row), fill = cell)) +
  geom_tile(colour = "white", linewidth = 1) +
  geom_text(aes(label = if_else(is.na(expo), "", sprintf("%.1f", expo))),
            size = 2.9, colour = "grey15") +
  scale_fill_manual(values = c("At risk, survived"         = "#D6E4EF",
                               "Death occurs here"         = "#E8A598",
                               "Not at risk (not reached)" = "grey96"), name = NULL) +
  labs(x = "DHS age segment (months)", y = NULL,
       title = "B  The same children after the person-period expansion",
       subtitle = str_wrap(paste(
         "One cell per child per segment it entered; numbers are months of exposure",
         "contributed. Blank cells are ages the child had not reached, about which",
         "the model is told nothing."), 108)) +
  theme_u5m() +
  theme(axis.text.y = element_text(family = "mono"),
        panel.grid = element_blank(), legend.position = "bottom")

save_fig(pA / pB + plot_layout(heights = c(1, 1.15)),
         "figS6_exposure_and_events", width = 11, height = 11)
message("[06b] done.")
