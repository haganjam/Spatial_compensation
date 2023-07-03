#' @title: analyse the counterfactual compensation scenarios
#' 
#' @description: this script analyses the data generated from the counterfactual
#' simulations of species loss and species compensation
#' 

# load relevant libraries
library(dplyr)
library(readr)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)

# load the plotting theme
source("code/helper-plotting-theme.R")
source("code/helper-miscellaneous.R")

# load the rds data
sp_dat <- readRDS(file = "output/compensation_analysis_data.rds")
head(sp_dat)

# convert to a tibble
sp_dat <- tibble(sp_dat)
head(sp_dat)

# calculate the productivity across depth zones initially and under compensation
sp_all <- 
  sp_dat %>%
  mutate(depth = "All") %>%
  group_by(depth, species_ext, comp_level, sample_gr) %>%
  summarise(prod_init = sum(prod_init),
            prod_comp = sum(prod_comp), .groups = "drop") %>%
  select(species_ext, comp_level, sample_gr, depth, prod_init, prod_comp)

# bind the species all onto the sp_dat data.frame
sp_dat <- bind_rows(sp_dat, sp_all)

# arrange the data in sp_dat
sp_dat <- 
  sp_dat %>%
  arrange(species_ext, comp_level, sample_gr, depth)

# calculate the change in productivity
sp_dat <- 
  sp_dat %>%
  mutate(prod_change = (prod_comp - prod_init)) %>%
  mutate(depth = as.character(depth),
         comp_level = as.character(comp_level))

# check the summary
summary(sp_dat)

# change the factor levels of depth
sp_dat$depth <- factor(sp_dat$depth)
levels(sp_dat$depth) <- c(paste0("DZ", c("1", "2", "3", "4")), "All")

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

# check the maximum and minimum prod change for each species
sp_sub %>%
  group_by(species_ext) %>%
  summarise(range = range(prod_change))

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

# calculate percentage change from the mean
sp_PPE <- 
  sp_dat %>%
  group_by(species_ext, depth, comp_level) %>%
  mutate(perc_change = (prod_change/prod_init)*100 ) %>%
  summarise(mu_change = mean(perc_change),
            PI_low = quantile(perc_change, 0.05),
            PI_high = quantile(perc_change, 0.95), .groups = "drop") %>%
  filter(comp_level %in% c(0.1, 0.5, 0.9))
View(sp_PPE)

# generate summary tables: table 1 and table S5

# generate a supplementary table with the rest of the percentage changes
tab_sum <- 
  sp_PPE %>%
  mutate(mu_and_PI = paste0(round(mu_change, 1), " [",
                            round(PI_low, 1), " - (",
                            round(PI_high, 1), ")]" )) %>%
  select(species_ext, depth, comp_level, mu_and_PI) %>%
  pivot_wider(id_cols = c("species_ext", "depth"),
              names_from = "comp_level",
              values_from = "mu_and_PI")
print(tab_sum)

# rename the columns
names(tab_sum) <- c("Species extinct", "Depth zone", "10% compensation", "50% compensation", "90% compensation")

# reorder the rows
tab_sum$`Species extinct` <- factor(tab_sum$`Species extinct`, levels = c("fu_sp", "fu_ve", "as_no", "fu_se"))
levels(tab_sum$`Species extinct`) <- c("F. spiralis", "F. vesiculosus", "A. nodosum", "F. serratus")

# arrange and rename the species from their binomial codes
tab_sum <- arrange(tab_sum, `Species extinct`)

# table 1: main text
tab_1 <- 
  tab_sum %>%
  filter(`Depth zone` == "All") %>%
  select(-`Depth zone`)

# export this table to a .csv file
write_csv(x = tab_1, file = "figures-tables/table_1.csv")

# table S5: supplementary information
tab_s5 <- 
  tab_sum %>%
  filter(`Depth zone` != "All")

# export this table to a .csv file
write_csv(x = tab_s5, file = "figures-tables/table_S5.csv")


# plot these four graphs individually for each species
sp_code <- c("fu_sp", "fu_ve", "as_no", "fu_se")
sp_names <- c("F. spiralis Extinct", "F. vesiculosus Extinct", "A. nodosum Extinct", "F. serratus Extinct")

# axis labels
ylabs <- c(expression("Change in productivity"~(g~day^{-1}) ), "", expression("Change in productivity"~(g~day^{-1}) ), "")
xlabs <- list(NULL, NULL, NULL, NULL)
x.text <- c("white", "white", "black", "black")
x.text.size <- c(1, 1, 12, 12)

# set the ylims for the different species
ylims <- list(c(-3, 3),
              c(-3.5, 3),
              c(-3, 3),
              c(-6, 6))

# get the bar size for each plot
all_zone <- 
  lapply(ylims, function(x) {
  
  mp <- sum(x)/2
  h <- (max(x)-min(x))
  v <- c("midpoint" = mp, "height" = h)
  return(v)
  
})

# get a colour palette
col_pal <- wesanderson::wes_palette(name = "Royal1", n = 4, type = "discrete")
col_pal <- col_pal[c(1,2,4)]

plots <- vector("list", length = length(sp_code))
for(i in 1:length(sp_code)) {
  
  # set-up a data.frame for plotting tiles behind the ALL zone
  all_zone_bar <- 
    tibble(depth = factor(levels(sp_dat$depth), levels = levels(sp_dat$depth)),
           midpoint = all_zone[[i]][["midpoint"]],
           bar_height = c(rep(0, 4), all_zone[[i]][["height"]]))
  
  px <- 
    ggplot() +
    geom_tile(data = all_zone_bar,
              mapping = aes(x = depth, y = midpoint, 
                            width = 0.75, height = bar_height),
              alpha = 0.05) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 2.5, linetype = "dashed", colour = "red", alpha = 0.5) +
    geom_hline(yintercept = -2.5, linetype = "dashed", colour = "red", alpha = 0.5) +
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
    scale_colour_manual(values = col_pal) +
    ylab(ylabs[i]) +
    xlab(xlabs[[i]]) +
    ggtitle(label = sp_names[i]) +
    labs(colour = "Compensation (%)") +
    guides(colour = guide_legend(override.aes = list(shape = 18, size = 6))) +
    scale_y_continuous(limits = ylims[[i]], expand = c(0, 0)) +
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

ggsave(filename = "figures-tables/fig_4.svg", p1,
       units = "cm", width = 20, height = 18)

### END
