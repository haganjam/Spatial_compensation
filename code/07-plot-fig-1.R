#' @title: plot fig 1
#' 
#' @description: this script combines three figures generated in the previous
#' scripts and plots fig 1 from the manuscript.
#' 

# load relevant libraries
library(ggplot2)
library(cowplot)

p1 <- readRDS(file = "output/fig_1a.rds")
p1 <- 
  p1 + 
  ggtitle("") +
  theme(legend.position = "none")

p2 <- readRDS(file = "output/fig_1b.rds")
p2 <- 
  p2 + 
  ggtitle("") +
  ylab(expression("Relative growth rate"~(g~g^{-1}~day^{-1}) )) +
  theme(legend.position = "none")

p3 <- readRDS(file = "output/fig_1c.rds")
p3 <- 
  p3 + 
  ggtitle("") +
  theme(legend.position = "none")

# ggarrange the plots to keep the common legend
p123 <- plot_grid(p1, p2, p3, nrow = 2, ncol = 2, align = "v",
                  labels = c("a", "b", "c"), label_size = 11,
                  label_fontface = "plain")
plot(p123)

ggsave(filename = "figures-tables/fig_1.svg", p123,
       units = "cm", width = 20, height = 18)

### END
