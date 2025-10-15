base_theme <- theme_classic() +
  theme(
    plot.title = element_text(size = 14, color = 'black'),
    plot.subtitle = element_text(size = 12, color = 'black'),
    axis.text = element_text(size = 10, color = 'black'),
    axis.title = element_text(size = 12, color = 'black'),
    axis.ticks = element_line(color = 'black'),
    legend.title = element_text(size = 12, color = 'black'),
    legend.text = element_text(size = 10, color = 'black'),
    strip.background = element_blank(), # for faceted plots
    panel.grid = element_blank() # for faceted plots
  )
