# Shared plot-styling helpers for consistent figures across this
# project: theme_als_project() is a theme() overlay (not a full theme)
# pinning plot/axis/legend/strip text and legend key size to one
# consistent look, usable on any ggplot-returning call (raw ggplot2,
# DimPlot(), dittoSeq, etc.). rotate_x_text() and no_key_outline() are
# separate opt-in helpers (axis-text rotation; stripping geom_tile()/
# geom_col()'s outline colour from just the legend key, not the plot)
# since not every plot needs them. save_figure() wraps ggsave() with
# cairo_pdf as the default device for more faithful PDF font embedding.
# Meant to be source()'d from every figure script, not copied into each
# one, as a deliberate exception to this project's per-script-helpers
# convention.

suppressMessages({
  library(ggplot2)
})

theme_als_project <- function(base_size = 6,
                              title_size = base_size * 1.5,
                              title_face = "plain",
                              strip_fill = "white",
                              text_color = "black",
                              legend_key_size = unit(8, "pt")){
  theme(
    plot.title = element_text(size = title_size, face = title_face,
                              hjust = 0.5, color = text_color),
    axis.title.x = element_text(size = base_size, color = text_color,
                                hjust = 0.5),
    axis.title.y = element_text(size = base_size, color = text_color),
    axis.text = element_text(size = base_size, color = text_color),
    legend.title = element_text(size = base_size, color = text_color),
    legend.text = element_text(size = base_size, color = text_color),
    legend.key.size = legend_key_size,
    legend.key = element_rect(fill = NA, color = NA),
    strip.text = element_text(size = base_size, face = "plain",
                              color = text_color),
    strip.background = element_rect(fill = strip_fill, color = "black")
  )
}

# Angled x-axis text with the matching hjust/vjust pair automatically --
# see header note above. angle = 0 (the default) returns an empty
# theme(), so it's always safe to add even when you don't need rotation.
rotate_x_text <- function(angle = 45){
  if (angle == 0){
    return(theme())
  }

  if (angle == 90){
    hjust <- 1
    vjust <- 0.5
  } else {
    hjust <- 1
    vjust <- 1
  }

  theme(axis.text.x = element_text(angle = angle, hjust = hjust, vjust = vjust))
}

# Strips the outline colour from a legend's key swatches without
# touching the outline drawn in the plot itself -- geom_tile()/
# geom_col()'s own colour= outline parameter bleeds into the legend key
# glyph, and that's not something theme() alone can override, since the
# key glyph is drawn from the layer's own colour/linewidth, not a theme
# element. aes_name is whichever aesthetic generates the legend --
# almost always "fill" for geom_tile()/geom_col() in this project.
no_key_outline <- function(aes_name = "fill"){
  do.call(guides, setNames(
    list(guide_legend(override.aes = list(colour = NA, linewidth = NA))),
    aes_name
  ))
}

umap_labs <- labs(y = "UMAP 2",
                  x = "UMAP 1")

umap_theme <- theme(axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    aspect.ratio = 1)

# Optional convenience wrapper around ggsave() for the PDF -> Illustrator
# workflow -- see header note above. All ggsave() arguments still work;
# device is the only one given a different default.
save_figure <- function(plot, filename, device = cairo_pdf, ...){
  ggsave(plot = plot, filename = filename, device = device, ...)
}

# Usage (source this file once near the top of a figure script, after
# your own library() calls, then add theme_als_project() to any plot --
# raw ggplot2 or a wrapper function's output -- the same way you'd add
# any other theme):
#
# source("r_scripts/plot_theme.R")
#
# ggplot(df, aes(x = group, y = value)) +
#   geom_jitter() +
#   theme_linedraw(base_size = 12) +
#   theme_als_project() +
#   rotate_x_text(45)
#
# DimPlot_scCustom(obj, group.by = "cell_type3") +
#   theme_als_project(title_size = 20)
#
# ggplot(df, aes(x = group, y = var, fill = value)) +
#   geom_tile(color = "black", linewidth = 0.4) +
#   theme_als_project() +
#   no_key_outline()
#
# save_figure(p, "figures/some_plot.pdf", height = 5, width = 6)
