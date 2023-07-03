#' @title: perform the Schiettekatte et al. (2022) analysis on the macroalgae data
#' 
#' @description: this script analyses the macroalgae data using Schiettekatte et
#' al.'s (2022) method for calculating the the proportion of dominant species
#' at different sites

# load relevant libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

# load the winfree functions
source("code/helper-plotting-theme.R")
source("code/helper-miscellaneous.R")

# load the productivity data
prod_list <- readRDS("output/model_productivity.rds")

# plot the raw productivity across depth zones integrated across growth rates
prod <- 
  bind_rows(prod_list, .id = "sample_gr") %>%
  pivot_longer(cols = c("fu_se", "as_no", "fu_ve", "fu_sp"),
               names_to = "species",
               values_to = "prod")

# translate the productivity data so that it is above zero
prod <- 
  prod %>%
  group_by(sample_gr, depth_zone) %>%
  mutate(prod_trans = prod + abs(min(prod))) %>%
  ungroup()

# calculate species contributions
# calculate the species richness at each site for each growth rate sample
prod <- 
  prod %>%
  group_by(sample_gr, depth_zone) %>%
  mutate(prod_contr = prod_trans/sum(prod_trans),
         exp_contr = 1/sum(prod != 0)) %>%
  ungroup()

# classify species as dominant or not within a site
prod <- 
  prod %>%
  mutate(dom_01 = ifelse(prod_contr == 1, 1, 
                         ifelse(prod_contr > exp_contr, 1, 0)))

# how many species are dominant in at least one site
prop_sp_dom <- 
  prod %>%
  group_by(sample_gr) %>%
  filter(dom_01 == 1) %>%
  summarise(prop_sp_dom = length(unique(species))/4) %>%
  ungroup()
unique(prop_sp_dom$prop_sp_dom)

# check for missing values
any(is.na(prop_sp_dom$prop_sp_dom))
summary(prop_sp_dom$prop_sp_dom)

# get the mean and standard deviation
m_sd <- 
  prop_sp_dom %>%
  summarise(m = mean(prop_sp_dom),
            sd = sd(prop_sp_dom))

# get a decent red colour
col_pal <- wesanderson::wes_palette(name = "Darjeeling1", n = 1)

ggplot() +
  geom_histogram(data = prop_sp_dom,
                 mapping = aes(x = prop_sp_dom), 
                 bins = 20, fill = col_pal, colour = col_pal) +
  geom_point(data = m_sd,
             mapping = aes(x = m, y = 1100), size = 2,
             colour = col_pal) +
  geom_errorbarh(data = m_sd,
                mapping = aes(xmin = m-sd, xmax = m+sd, y = 1100),
                height = 0, colour = col_pal) +
  scale_x_continuous(limits = c(-0.1, 1.1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1200),
                     breaks = seq(0, 1200, 200)) +
  ylab("Number of samples") +
  xlab("Proportion of species being dominant") +
  theme_meta()
  







