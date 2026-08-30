library(tidyverse)
library(paletteer)
library(ggplot2)
library(scales)
library(janitor)
library(patchwork)

setwd("/projects/b1169/boles/als_cns_scrnaseq")

df <- read_csv("tab_data/target_als_demographics_compiled.csv")

df <- clean_names(df)

df <- df %>% 
  mutate(fm_hx_of_als_ftd = ifelse(is.na(fm_hx_of_als_ftd), "Unknown", fm_hx_of_als_ftd),
         c9orf72_mutation = ifelse(is.na(c9orf72_mutation), "N", c9orf72_mutation),
         race = ifelse(is.na(race), "Unknown", race),
         disease_duration = age_at_death - age_at_onset)

df <- df %>% 
  mutate(ffpe = case_when(case_number == "GBB-19-13" ~ "Y",
                             str_detect(case_number, "AU") ~ "N",
                             .default = "Y") %>% 
           factor(levels = c("Y", "N")),
         sc = case_when(case_number == "GBB-19-13" ~ "N",
                        .default = "Y") %>% 
           factor(levels = c("Y", "N")))

base_theme <- function(show_y_text = FALSE) {
  theme_minimal(base_size = 9) +
    theme(
      axis.text.x    = element_blank(),
      axis.ticks     = element_blank(),
      axis.title.x     = element_blank(),
      axis.title.y = element_text(angle = 0, vjust = 0.5),
      axis.text.y    = if (show_y_text) element_text(size = 7) else element_blank(),
      panel.grid     = element_blank(),
      # plot.title     = element_text(size = 8, angle = 90, hjust = 0, vjust = 0.5),
      # plot.title.position = "plot",
      legend.key.size = unit(0.35, "cm"),
      legend.text    = element_text(size = 6),
      legend.title   = element_text(size = 7, face = "bold")
    )
}

plot_continuous <- function(data, 
                            var, 
                            label, 
                            show_y_text = FALSE) {
  ggplot(data, aes(y = 1, 
                   x = case_number, 
                   fill = .data[[var]])) +
    geom_tile(color = "white", 
              linewidth = 0.4) +
    scale_fill_viridis_c(name = label, 
                         na.value = "grey85") +
    labs(y = label) +
    base_theme(show_y_text)
}

plot_categorical <- function(data, 
                             var, 
                             label, 
                             palette = NULL, 
                             show_y_text = FALSE) {
  p <- ggplot(data, aes(y = 1, 
                        x = case_number, 
                        fill = .data[[var]])) +
    geom_tile(color = "white", linewidth = 0.4) +
    labs(y = label) +
    base_theme(show_y_text)
  
  if (!is.null(palette)) {
    p + scale_fill_manual(values = palette, name = label, na.value = "grey85")
  } else {
    p + scale_fill_brewer(palette = "Set2", name = label, na.value = "grey85")
  }
}

df <- df %>% 
  arrange(clinical_diagnosis, c9orf72_mutation) %>%
  mutate(case_number = fct_inorder(case_number))

p_c9 <- plot_categorical(df,
                 "c9orf72_mutation",
                 label = "C9orf72 HRE")

p_age <- plot_continuous(df,
                "age_at_death",
                label = "Age at death (yr)")

p_group <- plot_categorical(df,
                 "clinical_diagnosis",
                 label = "Group")

p_duration <- plot_continuous(df,
                "disease_duration",
                label = "Disease duration")

p_group + p_c9 + p_age + p_duration +
  plot_layout(ncol = 1, widths = c(1.6, rep(1, 8))) +
  plot_annotation(title = "ALS Cohort Donor Demographics")

results_dir <- "results/demographics_figure/"
dir.create(results_dir, showWarnings = F, recursive = T)

ggsave(paste0(results_dir, "demographics_heatmap.png"),
       units = "in", dpi = 600,
       height = 5, width = 13)
