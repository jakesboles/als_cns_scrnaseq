# Shared plot-styling helpers for consistent figures across this
# project: theme_als_project() is a theme() overlay (not a full theme --
# added on top of whatever base theme/panel look a plot already has)
# pinning plot/axis/legend/strip text to one consistent size/face/color/
# justification, usable on any ggplot-returning call (raw ggplot2,
# DimPlot(), dittoSeq, etc.). rotate_x_text() pairs an axis-text angle
# with its correct hjust/vjust automatically. save_figure() wraps
# ggsave() with cairo_pdf as the default device, for more faithful font
# embedding in PDFs headed into Illustrator. Meant to be source()'d from
# every figure script, not copied into each one -- deliberately
# centralized, unlike this project's usual per-script helpers, since
# copies of a consistency tool can drift exactly like this project's
# existing scattered theme() calls already have.

suppressMessages({
  library(ggplot2)
})

theme_als_project <- function(base_size = 12,
                              title_size = base_size * 1.5,
                              title_face = "plain",
                              strip_fill = "gray80",
                              text_color = "black"){
  theme(
    plot.title = element_text(size = title_size, face = title_face,
                              hjust = 0.5, color = text_color),
    axis.title.x = element_text(size = base_size, color = text_color,
                                hjust = 0.5),
    axis.title.y = element_text(size = base_size, color = text_color),
    axis.text = element_text(size = base_size, color = text_color),
    legend.title = element_text(size = base_size, color = text_color),
    legend.text = element_text(size = base_size, color = text_color),
    strip.text = element_text(size = base_size, face = "bold",
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
# save_figure(p, "figures/some_plot.pdf", height = 5, width = 6)
