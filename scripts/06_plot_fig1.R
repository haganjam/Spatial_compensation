
# Compensation analysis conceptual figure

# load relevant libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# set the working directory
setwd("C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_2_Fucus_landscape/compensation_analysis")
getwd()

# load the plotting theme
source("scripts/function_plotting_theme.R")

p1 <- readRDS(file = "figures/Fig_1a.rds")
p1 <- 
  p1 + 
  ggtitle("")
p2 <- readRDS(file = "figures/Fig_1b.rds")
p2 <- 
  p2 + 
  ggtitle("")
p3 <- readRDS(file = "figures/Fig_1c.rds")
p3 <- 
  p3 + 
  ggtitle("")

# ggarrange the plots to keep the common legend
p123 <- ggarrange(p1, p2, p3, common.legend = TRUE, legend = "top",
                 labels = c("a", "b", "c"), font.label = list(face = "plain", size = 11),
                 ncol = 3, nrow = 1,
                 widths = c(1, 1.1, 1))
plot(p123)

ggsave(filename = "figures/fig1.png", p123, dpi = 400,
       units = "cm", width = 21, height = 10)

### END

