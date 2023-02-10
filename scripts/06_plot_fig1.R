
# Compensation analysis conceptual figure

# load relevant libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(here)

# load the plotting theme
source(here("scripts/function_plotting_theme.R"))

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
p123 <- ggarrange(p1, p2, p3, common.legend = FALSE,
                 labels = c("a", "b", "c"), font.label = list(face = "plain", size = 11),
                 ncol = 2, nrow = 2,
                 widths = c(1, 1, 1))
plot(p123)

ggsave(filename = "figures/fig1.png", p123, dpi = 400,
       units = "cm", width = 20, height = 17)

### END
