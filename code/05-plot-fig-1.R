
# Compensation analysis conceptual figure

# load relevant libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(here)

# load the plotting theme
source(here("scripts/function1_plotting_theme.R"))

p1 <- readRDS(file = here("figures/fig1a.rds"))
p1 <- 
  p1 + 
  ggtitle("") +
  theme(legend.position = "none")

p2 <- readRDS(file = here("figures/fig1b.rds"))
p2 <- 
  p2 + 
  ggtitle("") +
  theme(legend.position = "none")

p3 <- readRDS(file = here("figures/fig1c.rds"))
p3 <- 
  p3 + 
  ggtitle("") +
  theme(legend.position = "none")

# ggarrange the plots to keep the common legend
p123 <- plot_grid(p1, p2, p3, nrow = 2, ncol = 2, align = "v",
                  labels = c("a", "b", "c"), label_size = 11,
                  label_fontface = "plain")
plot(p123)

ggsave(filename = "figures/fig1.png", p123, dpi = 400,
       units = "cm", width = 20, height = 18)

### END
