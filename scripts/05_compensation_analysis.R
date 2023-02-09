
# Analyse the compensation data

# load relevant libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)
library(here)

# load the plotting theme
source(here("scripts/function_plotting_theme.R"))

# load the rds data
sp_dat <- readRDS(file = here("data/species_extinction_analysis.rds"))
head(sp_dat)

# convert to a tibble
sp_dat <- tibble(sp_dat)

# calculate the number how many zones maintain 50% of the productivity
# calculate the percent deviation from the initial productivity
sp_dat <- 
  sp_dat %>%
  mutate(thresh_50 = (prod_comp >= 0.50*prod_init),
         prod_change = ((prod_comp - prod_init)/prod_init)*100  ) %>%
  mutate(depth_zone = as.character(depth),
         comp_level = as.character(comp_level))

# remove the most extreme values
sp_dat <- 
  sp_dat %>%
  filter(prod_change > quantile(prod_change, 0.01),
         prod_change < quantile(prod_change, 0.99))

# check the summary data
summary(sp_dat)

# sample 100 points randomly from each category
sp_sub <- 
  sp_dat %>%
  group_by(species_ext, depth_zone, comp_level) %>%
  mutate(sd = sd(prod_change)) %>%
  sample_n(size = 100) %>%
  filter(comp_level %in% c(0.1, 0.5, 0.9)) %>%
  filter(sd != 0)
summary(sp_sub)

# summarise the data
sp_sum <- 
  sp_dat %>%
  group_by(species_ext, depth_zone, comp_level) %>%
  summarise(mu = mean(prod_change),
            PI_low = quantile(prod_change, 0.05),
            PI_high = quantile(prod_change, 0.95),
            n = n()) %>%
  filter(comp_level %in% c(0.1, 0.5, 0.9))

# plot these four graphs individually for each species
sp_code <- c("fu_sp", "fu_ve", "as_no", "fu_se")
sp_names <- c("F. spiralis Extinct", "F. vesiculosus Extinct", "A. nodosum Extinct", "F. serratus Extinct")

plots <- vector("list", length = length(sp_code))
for(i in 1:length(sp_code)) {
  
  px <- 
    ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 50, linetype = "dashed", colour = "red") +
    geom_hline(yintercept = -50, linetype = "dashed", colour = "red") +
    geom_quasirandom(data = sp_sub %>% filter(species_ext == sp_code[i]),
                     mapping = aes(x = depth_zone, y = prod_change, colour = comp_level),
                     dodge.width = 0.8, alpha = 0.2, shape = 1, size = 2.25,
                     show.legend=FALSE) +
    geom_point(data = sp_sum %>% filter(species_ext == sp_code[i]),
               mapping = aes(x = depth_zone, y = mu, colour = comp_level), 
               position = position_dodge(0.8), shape = 18, size = 3.5) +
    geom_errorbar(data = sp_sum %>% filter(species_ext == sp_code[i]),
                  mapping = aes(x = depth_zone, colour = comp_level,
                                ymin = PI_low, ymax = PI_high),
                  width = 0.05,
                  position = position_dodge(0.8),
                  show.legend=FALSE) +
    scale_colour_viridis_d(option = "E", end = 0.9) +
    ylab("Change in productivity (%)") +
    xlab("Depth zone") +
    ggtitle(label = sp_names[i]) +
    labs(colour = "Compensation (%)") +
    guides(colour = guide_legend(override.aes = list(shape = 18, size = 6))) +
    scale_y_continuous(limits = c(-480, 325)) +
    theme(legend.key = element_rect(fill = NA),
          legend.position = "top") +
    theme_meta() +
    theme(legend.title = element_text(size = 12),
          legend.text = element_text(size = 11.5),
          plot.title = element_text(size = 12, hjust = 0.5))
  
  # add to list
  plots[[i]] <- px
  
}

# join these graphs into a single graph
p1 <- 
  ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
            nrow = 2, ncol = 2,
            common.legend = TRUE, legend = "bottom",
            labels = c("a", "b", "c", "d"),
            font.label = list(face = "plain", size = 11))
plot(p1)

ggsave(filename = "figures/fig4.png", p1, dpi = 400,
       units = "cm", width = 18, height = 17)

### END
