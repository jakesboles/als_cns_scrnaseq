library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggbeeswarm)
library(ggpubr)
library(rstatix)

setwd("/projects/b1169/boles/als_cns_scrnaseq")

in_dir <- "results/speckle/"

plots_dir <- "figures/"

# Stat figures ------------------------------------------------------------

stats <- list.files(in_dir,
                    pattern = "stats")
tissues <- str_split_i(stats, "_", i = 1)

stats <- map(paste0(in_dir, stats),
             read.csv)

for (i in seq_along(stats)) {
  stats[[i]]$tissue <- tissues[i]
}

stats <- list_rbind(stats)

stats <- stats %>% 
  mutate(sals_fc = log2(PropMean.sALS / PropMean.Control),
         c9_fc = log2(PropMean.C9orf72 / PropMean.Control),
         sig = if_else(FDR < 0.05, "sig", "ns") %>% 
           factor(levels = c("sig", "ns"),
                  labels = c("True", "False")),
         tissue = factor(tissue,
                         levels = c("brain", "sc", "muscle"),
                         labels = c("Motor cortex", "Cervical spinal cord", "Skeletal muscle")))
  
sig_cells <- stats %>% 
  filter(FDR < 0.05)  

stats  %>% 
  mutate(label = if_else(FDR < 0.05, X, "")) %>%
  ggplot(aes(x = sals_fc,
             y = c9_fc)) + 
  facet_wrap(. ~ tissue,
             nrow = 1) + 
  geom_point(aes(fill = sig,
                 size = BaselineProp),
             shape = 21) + 
  geom_label_repel(aes(label = label, color = sig),
                   show.legend = FALSE,
                   seed = 42,                    # reproducible layout
                   force = 5,
                   force_pull = 0.5,
                   max.overlaps = Inf,            # never drop a label
                   box.padding = 0.6,             # space around label
                   point.padding = 0.4,           # space around the point it's labeling
                   min.segment.length = 0,        # always draw the leader line
                   segment.color = "grey40",
                   # THIS is the part that avoids all points, not just labeled ones:
                   max.time = 5,
                   max.iter = 100000) +
  scale_color_manual(values = c("dodgerblue", "grey")) +
  scale_fill_manual(values = c("dodgerblue", "grey")) + 
  labs(fill = "Adjusted p < 0.05",
       size = "Average proportion of\ntotal tissue's cells",
       y = "log2FC in C9orf72-ALS vs Control",
       x = "log2FC in sALS vs Control") + 
  guides(fill = guide_legend(ncol = 1,
                             override.aes = list(size = 4))) + 
  theme_linedraw(base_size = 12) + 
  theme(strip.text = element_text(face = "bold", color = "black"),
        strip.background = element_rect(fill = "gray", color = "black"))
ggsave(filename = paste0(plots_dir, "speckle_volcanoes.png"),
       units = "in", dpi = 600,
       height = 3, width = 10)

# Dot plots ------------------------------------------------------------

props <- list.files(in_dir, pattern = "frequencies")
props <- map(paste0(in_dir, props),
             read.csv)
names(props) <- tissues

sig_cells_guide <- stats %>% 
  filter(FDR < 0.05) %>%
  dplyr::select(X, tissue) %>% 
  mutate(tissue2 = case_when(tissue == "Motor cortex" ~ "brain",
                             tissue == "Cervical spinal cord" ~ "sc",
                             tissue == "Skeletal muscle" ~ "muscle"))

# geom_pwc() stacks one bracket per significant pairwise comparison above
# the tallest point, but the panel's y-expansion didn't account for how
# many brackets end up stacked -- with up to 3 pairwise comparisons (3
# groups), 2-3 stacked brackets ran off the top of the saved file even
# though the plot looked fine interactively (where the device auto-sizes).
# The top expansion is made a function of n_sig (how many pairwise
# comparisons are actually significant, hide.ns = T drops the rest)
# instead of a fixed guess, so a cell with only 1 significant comparison
# isn't given unnecessary headroom and one with 3 gets enough.
# step.increase (bracket-to-bracket vertical spacing, as a fraction of the
# y-range) is set explicitly here so the expansion formula's assumptions
# about ggpubr's own spacing match what's actually drawn.
#
# Per the user, the saved file height stays fixed (ggsave() below is back
# to height = 3, not scaled by n_sig) -- expansion is the only lever now,
# so it has to reserve enough of that *fixed* panel for the brackets
# rather than getting help from a taller canvas. The real tradeoff: for a
# cell with 2-3 significant comparisons, the data points themselves get
# visually compressed into a smaller fraction of the fixed panel height,
# since more of that same fixed space is reserved for brackets at the
# top. The constants below (0.2 base, 0.22 per extra bracket) are a
# first-pass estimate, not a rendered/verified fit -- no R environment
# available here to check actual pixel output, so these may need
# further visual tuning once you can see real output, especially for the
# 3-significant-comparison case.

step_increase <- 0.15 # bracket-to-bracket vertical spacing fraction -- change as needed, and update top_expansion below to match if you do

for (i in 1:nrow(sig_cells_guide)) {
  tissue <- sig_cells_guide$tissue2[i]
  cell <- sig_cells_guide$X[i]
  cell <- str_replace_all(cell, " ", ".")

  df <- props[[tissue]] %>%
    mutate(group = factor(group,
                          levels = c("Control", "sALS", "C9orf72"),
                          labels = c("Control", "sALS", "C9orf72-ALS")))

  # Same test geom_pwc(method = "tukey_hsd") runs internally -- computed
  # explicitly here just to count how many brackets will actually be
  # drawn (hide.ns = T only shows p.adj < 0.05 ones).
  tukey_res <- df %>%
    rstatix::tukey_hsd(as.formula(paste0("`", cell, "` ~ group")))
  n_sig <- sum(tukey_res$p.adj < 0.05, na.rm = T)

  # Headroom above the tallest point: a base allowance for the first
  # bracket/label plus a per-bracket increment for each additional one --
  # larger than step_increase alone, since expansion is now the only
  # lever reserving space within a fixed-height panel (see note above).
  top_expansion <- 0.1 * max(n_sig - 1, 0)

  p <- df %>%
    ggplot(aes(x = group,
               y = !!sym(cell))) +
    geom_quasirandom(aes(fill = group),
                     show.legend = F,
                     shape = 21,
                     size = 4,
                     alpha = 0.7) +
    stat_summary(fun = mean,
                 geom = "crossbar") +
    stat_summary(fun.data = mean_se,
                 geom = "errorbar",
                 linewidth = 1.2,
                 width = 0.6) +
    geom_pwc(method = "tukey_hsd",
             hide.ns = F,
             step.increase = step_increase,
             label = "p = {p.adj.format}") +
    scale_y_continuous(expand = expansion(mult = c(0.05, top_expansion))) +
    scale_fill_manual(values = c("#b8b0a8", "#CC00FF", "#0CAA00")) +
    labs(y = "% of all cells") +
    ggtitle(paste0(str_replace_all(cell, "[.]", " "), " in\n", str_to_lower(sig_cells_guide$tissue[i]))) +
    theme_linedraw(base_size = 12) +
    theme(strip.text = element_text(face = "bold", color = "black"),
          strip.background = element_rect(fill = "gray", color = "black"),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5))
  ggsave(p,
         filename = paste0(plots_dir, "speckle_dots_", cell, "_", tissue, ".png"),
         units = "in", dpi = 600,
         height = 3, width = 3)
}
