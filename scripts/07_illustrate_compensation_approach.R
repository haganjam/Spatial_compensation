
# Figure illustrating compensation approach

# load relevant libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggpubr)

# set the working directory
setwd("C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_2_Fucus_landscape/compensation_analysis")
getwd()

# load the winfree functions
source("scripts/function_plotting_theme.R")

# set the number of days to simulate growth over
n_days <- 1

# load the summarised transect data
tra_dat <- read_csv(file = "data/transect_comp.csv")
head(tra_dat)

# copy tra_dat to modify it
df.tra <- tra_dat

# load the growth rate data
grow_dat_list <- readRDS(file = "data/growth_rate_data.rds")
grow_dat <- grow_dat_list[[500]]
head(grow_dat)

# choose a species to go extinct
sp_ext <- c("fu_se")

# choose a level of compensation
comp <- 0.30

# get the compensation level of standing biomass
comp.x <- tra_dat[[sp_ext]]*comp

# set the biomass of the extinct species to zero
df.tra[[sp_ext]] <- 0

# replace that biomass with one of the species that has positive growth in that zone
for(dz in 1:4) {
  
  # get the different names that have positive growth rates in that zone
  x.grow <- grow_dat[,names(grow_dat) != sp_ext]
  x.grow <- x.grow[x.grow$depth_zone == dz,-1]
  sp.grow <- names(x.grow)[x.grow > 0]
  
  # get the names of species with positive standing biomass in an adjacent zone
  
  # first, we get suitable zones
  zones <- c(dz-1, dz, dz+1)
  zones <- zones[(zones > 0) & (zones < 5)]
  
  x.ssb <- df.tra[, names(df.tra) != sp_ext]
  x.ssb <- x.ssb[x.ssb$depth_zone %in% zones,-1]
  sp.ssb <- names(x.ssb)[apply(x.ssb, 2, function(x) any(x >0))]

  # get intersecting species in this way
  sp_comp <- intersect(sp.grow, sp.ssb)
  
  if(length(sp_comp) > 0) {
    
    # sample one of those names
    sp_comp <- sample(sp_comp, 1)
    
    # the sp_comp is going to compensate for the dz depth zone
    df.tra[dz,][[sp_comp]] <- (df.tra[dz,][[sp_comp]] + comp.x[dz])
    
  }
  
}

fig_list <- 
  lapply(list(tra_dat, df.tra), function(x) {
  
  # plot a figure of the depth distribution
  tra_a_plot <- 
    x %>%
    rename(`F. serratus` = "fu_se",
           `A. nodosum` = "as_no",
           `F. vesiculosus` = "fu_ve",
           `F. spiralis` = "fu_sp") %>%
    pivot_longer(cols = c(`F. serratus`, `A. nodosum`, `F. vesiculosus`, `F. spiralis`),
                 names_to = "species", 
                 values_to = "biomass")
  
  # change levels of the factor
  tra_a_plot$species <- factor(tra_a_plot$species,
                               levels = c("F. serratus", "A. nodosum", "F. vesiculosus", "F. spiralis"))
  
  p1 <- 
    ggplot(data = tra_a_plot,
           mapping = aes(x = depth_zone, y = biomass, fill = species)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8,
             alpha = 0.75, colour = "black") +
    scale_fill_viridis_d(option = "A", end = 0.9) +
    ylab("Standing dry biomass (g)") +
    xlab("Depth zone") +
    scale_y_continuous(limits = c(0, 570)) +
    labs(fill = "Species") +
    theme_meta()
  
  return(p1)
  
} )

p1 <- 
  fig_list[[1]] +
  ggtitle("Intact community") +
  theme(plot.title = element_text(size = 11, hjust = 0.5))

p2 <- 
  fig_list[[2]] +
  ggtitle("F. serratus Extinction + (30% comp.)") +
  theme(plot.title = element_text(size = 11, hjust = 0.5))

# plot a single sample from the growth rate data
grow_plot <- 
  grow_dat %>%
  pivot_longer(cols = c("fu_se", "as_no", "fu_ve",  "fu_sp"),
               names_to = "species",
               values_to = "growth_rate")

# reorder the factor levels
grow_plot$species  <- factor(grow_plot$species, levels = c("fu_se", "as_no", "fu_ve",  "fu_sp"))

p3 <- 
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = grow_plot,
             mapping = aes(x = depth_zone, y = growth_rate, colour = species),
             position = position_dodge(0.25), size = 2) +
  geom_line(data = grow_plot,
            mapping = aes(x = depth_zone, y=growth_rate, colour = species),
            position = position_dodge(0.25)) +
  scale_colour_viridis_d(option = "A", end = 0.9) +
  xlab("Depth zone") +
  ylab(expression("Dry biomass change"~(g~g^{-1}~day^{-1}) )) +
  ggtitle("") +
  theme_meta() +
  theme(legend.position = "none",
        plot.title = element_text(size = 11))
plot(p3)

# calculate productivity in the intact and without compensation

# calculate raw productivity for each sample of growth rates
prod_raw <- 
  lapply(list(tra_dat, df.tra), function(ssb_dat) {
    
    # convert biomass to productivity units
    prodx <- mapply(function(x, y){x*y*n_days}, ssb_dat[,-1], grow_dat[,-1])
    prodx <- cbind(data.frame(depth_zone = as.character(1:4)), as.data.frame(prodx))
    return(prodx) }
  )

# plot the raw productivity across depth zones integrated across growth rates
prod <- 
  bind_rows(prod_raw, .id = "loss") %>%
  pivot_longer(cols = c("fu_se", "as_no", "fu_ve", "fu_sp"),
               names_to = "species",
               values_to = "prod") %>%
  group_by(depth_zone, species)

# change factor levels of species
prod$species <- factor(prod$species, levels = c("fu_se", "as_no", "fu_ve", "fu_sp"))

# change factor levels of loss
prod$loss <- factor(prod$loss)
levels(prod$loss) <- c("Intact", "F. serratus Extinct + (30% comp.)")

# summarise the data for each zone
prod_sum <- 
  prod %>%
  group_by(loss, depth_zone) %>%
  summarise(prod_sum = sum(prod))

p4 <- 
  ggplot(prod_sum,
         mapping = aes(x = depth_zone, y = prod_sum,
                       fill = loss)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_bar(width = 0.5, stat = "identity", position = "dodge",
           alpha = 0.65, colour = "black") +
  scale_fill_viridis_d(option = "D", begin = 0.3, end = 0.7) +
  xlab("Depth zone") +
  ylab(expression("Dry biomass prod."~(g~day^{-1}) )) +
  theme_meta() +
  ggtitle("") +
  theme(legend.position = c(0.60, 0.85),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        plot.title = element_text(size = 11))
plot(p4)

# arrange the figure into four panels

# arrange the first two figures with common legend
p12 <- 
  ggarrange(p1, p2,
            ncol = 2, nrow = 1, common.legend = TRUE,
            legend = "top",
            labels = c("a", "b"), font.label = list(face = "plain", size = 11))
plot(p12)

# arrange the second two figures without a common legend
p34 <- 
  ggarrange(p3, p4,
            ncol = 2, nrow = 1,
            labels = c("c", "d"), font.label = list(face = "plain", size = 11))
plot(p34)

# bind these two arrange plots
p1234 <- ggarrange(p12, p34, nrow = 2, ncol = 1)

ggsave(filename = "figures/fig3.png", p1234, dpi = 400,
       units = "cm", width = 20, height = 18)

### END
