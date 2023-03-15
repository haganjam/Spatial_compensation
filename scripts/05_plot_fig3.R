
# Figure illustrating compensation approach

# load relevant libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(here)

# load the winfree functions
source(here("scripts/function1_plotting_theme.R"))

# load the summarised transect data
tra_dat <- readRDS(file = here("data/transect_ssdb.rds"))
head(tra_dat)

# copy tra_dat to modify it
df.tra <- tra_dat

# load the growth rate data
grow_dat_list <- readRDS(file = "data/growth_rate_data.rds")
# n <- sample(1:length(grow_dat_list), 1)
# print(n)
grow_dat <- grow_dat_list[[1639]]
head(grow_dat)

# choose a species to go extinct
sp_ext <- c("fu_se")

# choose a level of compensation
comp <- 0.10

# get the compensation level of standing biomass
comp.x <- tra_dat[[sp_ext]]*comp

# set the biomass of the extinct species to zero
df.tra[[sp_ext]] <- 0

# replace that biomass with one of the species that has positive growth in that zone
for(dz in 1:4) {
  
  # get the different names that have positive growth rates in that zone
  x.grow <- grow_dat[,names(grow_dat) != sp_ext]
  x.grow <- x.grow[x.grow$depth == dz,-1]
  sp.grow <- names(x.grow)[x.grow > 0]
  
  # get the names of species with positive standing biomass in an adjacent zone
  
  # first, we get suitable zones
  zones <- c(dz-1, dz, dz+1)
  zones <- zones[(zones > 0) & (zones < 5)]
  
  x.ssb <- df.tra[, names(df.tra) != sp_ext]
  x.ssb <- x.ssb[x.ssb$depth %in% zones,-1]
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

# make a data list
d.list <- list(tra_dat, df.tra)
ylabs <- c("Standing dry biomass (g)", "")

# get an output list
fig_list <- vector("list", length = 2)

for(i in 1:length(fig_list)) {
  
  # plot a figure of the depth distribution
  tra_a_plot <- 
    d.list[[i]] %>%
    rename(`F. serratus` = "fu_se",
           `A. nodosum` = "as_no",
           `F. vesiculosus` = "fu_ve",
           `F. spiralis` = "fu_sp") %>%
    pivot_longer(cols = c(`F. serratus`, `A. nodosum`, `F. vesiculosus`, `F. spiralis`),
                 names_to = "species", 
                 values_to = "biomass")
  
  # change levels of the factor
  tra_a_plot$species <- factor(tra_a_plot$species,
                               c("F. spiralis", "F. vesiculosus", "A. nodosum", "F. serratus"))
  
  # change the levels of the depth factor
  tra_a_plot$depth <- factor(tra_a_plot$depth, levels = c("4", "3", "2", "1"))
  levels(tra_a_plot$depth) <- paste0("DZ", c("1", "2", "3", "4"))
  
  p1 <- 
    ggplot(data = tra_a_plot,
           mapping = aes(x = depth, y = biomass, fill = species)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8,
             alpha = 0.75, colour = "black") +
    scale_fill_viridis_d(option = "A", end = 0.9) +
    ylab(ylabs[i]) +
    xlab(NULL) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 570)) +
    labs(fill = "Species") +
    geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dashed", colour = "grey") +
    theme_meta() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 12))
  
  fig_list[[i]] <- p1
  
}

p1 <- 
  fig_list[[1]] +
  ggtitle("Intact community") +
  theme(plot.title = element_text(size = 11, hjust = 0.5))
plot(p1)

p2 <- 
  fig_list[[2]] +
  ggtitle("F. serratus Extinction + (10% comp.)") +
  theme(plot.title = element_text(size = 11, hjust = 0.5))
plot(p2)

# plot a single sample from the growth rate data
grow_plot <- 
  grow_dat %>%
  pivot_longer(cols = c("fu_se", "as_no", "fu_ve",  "fu_sp"),
               names_to = "species",
               values_to = "growth_rate")

# change levels of the species factor
grow_plot$species <- factor(grow_plot$species, 
                            levels = c("fu_sp", "fu_ve", "as_no", "fu_se"))
levels(grow_plot$species) <- c("F. spiralis", "F. vesiculosus", "A. nodosum", "F. serratus")

# change the levels of the depth factor
grow_plot$depth <- factor(grow_plot$depth, levels = c("4", "3", "2", "1"))
levels(grow_plot$depth) <-  paste0("DZ", c("1", "2", "3", "4"))

p3 <- 
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = grow_plot,
             mapping = aes(x = depth, y = growth_rate, colour = species),
             position = position_dodge(0.25), size = 2) +
  geom_line(data = grow_plot,
            mapping = aes(x = as.integer(depth), y=growth_rate, colour = species),
            position = position_dodge(0.25), show.legend = FALSE) +
  scale_colour_viridis_d(option = "A", end = 0.9) +
  xlab("") +
  ylab(expression("Dry biomass change"~(g~g^{-1}~day^{-1}) )) +
  ggtitle("") +
  theme_meta() +
  theme(plot.title = element_text(size = 11),
        legend.position = "none",
        axis.text.x = element_text(size = 12))
plot(p3)

# calculate productivity in the intact and without compensation

# calculate raw productivity for each sample of growth rates
prod_raw <- 
  lapply(list(tra_dat, df.tra), function(ssb_dat) {
    
    # convert biomass to productivity units
    prodx <- mapply(function(x, y){x*y}, ssb_dat[,-1], grow_dat[,-1])
    prodx <- cbind(data.frame(depth = as.character(1:4)), as.data.frame(prodx))
    return(prodx) }
  )

# plot the raw productivity across depth zones integrated across growth rates
prod <- 
  bind_rows(prod_raw, .id = "loss") %>%
  pivot_longer(cols = c("fu_se", "as_no", "fu_ve", "fu_sp"),
               names_to = "species",
               values_to = "prod") %>%
  group_by(depth, species)

# change factor levels of species
prod$species <- factor(prod$species, levels = c("fu_sp", "fu_ve", "as_no", "fu_se"))

# change factor levels of loss
prod$loss <- factor(prod$loss)
levels(prod$loss) <- c("Intact", "F. serratus Extinct + (10% comp.)")

# summarise the data for each zone
prod_sum <- 
  prod %>%
  group_by(loss, depth) %>%
  summarise(prod_sum = sum(prod))

# change the levels of the depth factor
prod_sum$depth <- factor(prod_sum$depth, levels = c("4", "3", "2", "1"))
levels(prod_sum$depth) <-  paste0("DZ", c("1", "2", "3", "4"))

p4 <- 
  ggplot(prod_sum,
         mapping = aes(x = depth, y = prod_sum,
                       fill = loss)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_bar(width = 0.5, stat = "identity", position = "dodge",
           alpha = 0.65, colour = "black") +
  scale_fill_manual(values = c("black", "grey")) +
  xlab("") +
  ylab(expression("Dry biomass prod."~(g~day^{-1}) )) +
  theme_meta() +
  ggtitle("") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 6)) +
  theme(legend.position = c(0.4, 0.85),
        legend.title = element_blank(),
        legend.text = element_text(size = 9),
        plot.title = element_text(size = 11),
        axis.text.x = element_text(size = 12),
        )
plot(p4)

# arrange the figure into four panels
p1234 <- plot_grid(p1, p2, p3, p4, nrow = 2, ncol = 2, align = "v",
                   labels = c("a", "b", "c", "d"), label_size = 11,
                   label_fontface = "plain")

ggsave(filename = "figures/fig3.png", p1234, dpi = 400,
       units = "cm", width = 20, height = 18)

### END
