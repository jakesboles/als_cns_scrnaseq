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
         c9orf72_mutation = ifelse(is.na(c9orf72_mutation), "N", c9orf72_mutation) %>%
           factor(levels = c("Y", "N")),
         race = ifelse(is.na(race), "Unknown", race),
         disease_duration = age_at_death - age_at_onset)

df <- df %>% 
  mutate(ffpe = case_when(case_number == "GBB-19-13" ~ "Y",
                             str_detect(case_number, "AU") ~ NA,
                             .default = "Y"),
         sc = case_when(case_number == "GBB-19-13" ~ NA,
                        .default = "Y"))

unique(df$site_of_sx_onset)

df <- df %>% 
  mutate(onset_site = case_when(site_of_sx_onset == "bulbar/limb" ~ "Both",
                                str_detect(site_of_sx_onset, "(?i)leg|(?i)limb|(?i)hand") ~ "Limb",
                                str_detect(site_of_sx_onset, "(?i)bulbar") & 
                                  site_of_sx_onset != "bulbar/limb" ~ "Bulbar",
                                .default = NA) %>% 
           factor(levels = c("Bulbar", "Limb", "Both")))


base_theme <- function(show_y_text = FALSE) {
  theme_minimal(base_size = 12) +
    theme(
      axis.text.x    = element_blank(),
      axis.ticks     = element_blank(),
      axis.title.x     = element_blank(),
      axis.title.y = element_text(angle = 0, vjust = 0.5, hjust = 1),
      axis.text.y    = if (show_y_text) element_text(size = 7) else element_blank(),
      panel.grid     = element_blank(),
      plot.margin = margin(t = 1, r = 5, b = 1, l = 5, unit = "pt"),
      # plot.margin    = margin(t = 2, r = 5, b = 2, l = 5, unit = "pt"), 
      # plot.title     = element_text(size = 8, angle = 90, hjust = 0, vjust = 0.5),
      # plot.title.position = "plot",
      legend.key.size = unit(0.35, "cm"),
      legend.text    = element_text(size = 9),
      legend.title   = element_text(size = 12)
    )
}

add_na_strike <- function(plot, data, var, color = "red", linewidth = 0.5, inset = 0.42) {
  na_df <- data %>%
    mutate(y = 1, x = as.numeric(case_number)) %>%
    filter(is.na(.data[[var]]))
  
  if (nrow(na_df) == 0) return(plot)
  
  plot +
    geom_segment(
      data = na_df,
      aes(x = x - inset, xend = x + inset, y = y - inset, yend = y + inset),
      inherit.aes = FALSE, color = color, linewidth = linewidth, lineend = "round"
    )
}

plot_continuous <- function(data, 
                            var, 
                            label, 
                            palette = NULL,
                            show_y_text = FALSE) {
  p <- ggplot(data, aes(y = 1, 
                   x = case_number, 
                   fill = .data[[var]])) +
    geom_tile(color = "black", 
              linewidth = 0.4) +
    labs(y = label) +
    base_theme(show_y_text) + 
    theme(legend.title = element_blank())
  
  p <- add_na_strike(p, data, var)
  
  if (!is.null(palette)) { 
    p + scale_fill_paletteer_c(palette, na.value = "white",
                               guide = guide_colourbar(direction = "horizontal"))
  } else { 
    p + scale_fill_paletteer_c("viridis::viridis", na.value = "white",
                               guide = guide_colourbar(direction = "horizontal"))
  }
  
}

plot_categorical <- function(data, 
                             var, 
                             label, 
                             palette = NULL, 
                             show_y_text = FALSE) {
  p <- ggplot(data, aes(y = 1, 
                        x = case_number, 
                        fill = .data[[var]])) +
    geom_tile(color = "black", linewidth = 0.4) +
    labs(y = label) +
    base_theme(show_y_text) + 
    theme(legend.title = element_blank())
  
  p <- add_na_strike(p, data, var)
  
  if (!is.null(palette)) {
    p + scale_fill_manual(values = palette, name = label, na.value = "white",
                          na.translate = F)
  } else {
    p + scale_fill_brewer(palette = "Set2", name = label, na.value = "white",
                          na.translate = F)
  }
}

df <- df %>% 
  arrange(clinical_diagnosis, c9orf72_mutation) %>%
  mutate(case_number = fct_inorder(case_number))

p_c9 <- plot_categorical(df,
                 "c9orf72_mutation",
                 label = "C9orf72 HRE",
                 palette = c("#D6249F", "#E5E1DC"))

p_onset <- plot_categorical(df,
                            "onset_site",
                            label = "Site of onset",
                            palette = c("#FFD60A", "#8C7A1E", "#F0DE7D"))

p_age <- plot_continuous(df,
                "age_at_death",
                label = "Age at death",
                palette = "ggthemes::Red")

p_group <- plot_categorical(df,
                 "clinical_diagnosis",
                 label = "Group",
                 palette = c("#6a3d9a", "#c3a6e1", "#b8b0a8"))

p_duration <- plot_continuous(df,
                "disease_duration",
                label = "Disease duration",
                palette = "ggthemes::Orange")

p_sc <- plot_categorical(df,
                         "sc",
                         label = "scRNAseq",
                         palette = c("#2DC653")) + 
  theme(legend.position = "none")

p_ffpe <- plot_categorical(df,
                           "ffpe",
                           label = "Spatial biology",
                           palette = c("#00B4D8")) + 
  theme(legend.position = "none")

p <- p_group + p_c9 + p_age + p_duration + p_onset + p_sc + p_ffpe +
  plot_layout(ncol = 1, widths = c(1.6, rep(1, 8)))

results_dir <- "results/demographics_figure/"
dir.create(results_dir, showWarnings = F, recursive = T)

ggsave(p,
       filename = paste0(results_dir, "demographics_heatmap.png"),
       units = "in", dpi = 600,
       height = 2.5, width = 10)

write.csv(df,
          file = "tab_data/organized_metadata_for_plotting.csv",
          row.names = F)
