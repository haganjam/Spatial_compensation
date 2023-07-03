#' @title: plot fig 3
#' 
#' @description: combine the two figures generated for the different types
#' of observational analyses
#' 

# load relevant libraries
library(ggplot2)
library(cowplot)

# load fig 3a
p1 <- readRDS(file = "output/fig_3a.rds")
plot(p1)

p2 <- readRDS(file = "output/fig_3b.rds")
p2 <- 
  p2 + 
  ggtitle("") +
  theme(legend.position = "none",
        plot.title = element_text(size = 36))
plot(p2)

# ggarrange the plots to keep the common legend
p12 <- plot_grid(p1, p2, nrow = 1, ncol = 2, align = "v",
                 labels = c("a", "b"), label_size = 11,
                 vjust = 3,
                 label_fontface = "plain"
                 )
plot(p12)

ggsave(filename = "figures-tables/fig_3.svg", p12,
       units = "cm", width = 20, height = 10)

### END
