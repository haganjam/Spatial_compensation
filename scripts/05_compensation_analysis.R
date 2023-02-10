
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

# calculate the change in productivity
sp_dat <- 
  sp_dat %>%
  mutate(prod_change = (prod_comp - prod_init)) %>%
  mutate(prod_change_perc = (prod_change/prod_init)*100) %>%
  mutate(depth = as.character(depth),
         comp_level = as.character(comp_level))

# check the summary
summary(sp_dat)

# check the maxima: These marginal cases
sp_dat %>%
  filter(prod_change_perc > 10000)

# change the factor levels of depth
sp_dat$depth <- factor(sp_dat$depth, levels = c("4", "3", "2", "1"))
levels(sp_dat$depth) <- c("-2 to -14", "-14 to -26", "-26 to -38", "-38 to -50")

# check the summary data
summary(sp_dat)

# sample 100 points randomly from each category
sp_sub <- 
  sp_dat %>%
  group_by(species_ext, depth, comp_level) %>%
  sample_n(size = 100) %>%
  filter(comp_level %in% c(0.1, 0.5, 0.9)) %>%
  filter(prod_change != 0)
summary(sp_sub)

# summarise the data
sp_sum <- 
  sp_dat %>%
  group_by(species_ext, depth, comp_level) %>%
  summarise(mu = mean(prod_change),
            PI_low = quantile(prod_change, 0.05),
            PI_high = quantile(prod_change, 0.95),
            n = n(), .groups = "drop") %>%
  filter(comp_level %in% c(0.1, 0.5, 0.9))
print(sp_sum)

# plot these four graphs individually for each species
sp_code <- c("fu_sp", "fu_ve", "as_no", "fu_se")
sp_names <- c("F. spiralis Extinct", "F. vesiculosus Extinct", "A. nodosum Extinct", "F. serratus Extinct")

# axis labels
ylabs <- c(expression("Dry biomass prod."~(g~day^{-1}) ), "", expression("Dry biomass prod."~(g~day^{-1}) ), "")
xlabs <- list(NULL, NULL, "Depth range (cm)", "Depth range (cm)")
x.text <- c("white", "white", "black", "black")
x.text.size <- c(1, 1, 10, 10)

plots <- vector("list", length = length(sp_code))
for(i in 1:length(sp_code)) {
  
  px <- 
    ggplot() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_errorbar(data = sp_sum %>% filter(species_ext == sp_code[i]),
                  mapping = aes(x = depth, colour = comp_level,
                                ymin = PI_low, ymax = PI_high),
                  width = 0,
                  position = position_dodge(0.8),
                  show.legend=FALSE) +
    geom_quasirandom(data = sp_sub %>% filter(species_ext == sp_code[i]),
                     mapping = aes(x = depth, y = prod_change, colour = comp_level),
                     dodge.width = 0.8, alpha = 0.2, shape = 1, size = 2.25,
                     show.legend = FALSE) +
    geom_point(data = sp_sum %>% filter(species_ext == sp_code[i]),
               mapping = aes(x = depth, y = mu, colour = comp_level), 
               position = position_dodge(0.8), shape = 18, size = 3.5) +
    scale_colour_viridis_d(option = "C", end = 0.9) +
    ylab(ylabs[i]) +
    xlab(xlabs[[i]]) +
    ggtitle(label = sp_names[i]) +
    labs(colour = "Compensation (%)") +
    guides(colour = guide_legend(override.aes = list(shape = 18, size = 6))) +
    scale_y_continuous(limits = c(-10, 10)) +
    theme(legend.key = element_rect(fill = NA),
          legend.position = "top") +
    theme_meta() +
    theme(legend.title = element_text(size = 12),
          legend.text = element_text(size = 11.5),
          plot.title = element_text(size = 12, hjust = 0.5),
          axis.text.x = element_text(colour = x.text[i], size = x.text.size[i]))
  
  # add to list
  plots[[i]] <- px
  
}

# join these graphs into a single graph
p1 <- 
  ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]],
            nrow = 2, ncol = 2,
            common.legend = TRUE, legend = "bottom",
            labels = c("a", "b", "c", "d"),
            font.label = list(face = "plain", size = 11),
            heights = c(1, 1.2))
plot(p1)

ggsave(filename = here("figures/fig4.png"), p1, dpi = 400,
       units = "cm", width = 18, height = 17)

### END
