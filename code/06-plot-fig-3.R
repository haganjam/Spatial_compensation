#' @title: illustrate approach to estimating compensation following species loss
#' 
#' @description: this script generates a set of four figures that is used
#' to illustrate how I estimate the change in dry biomass productivity following
#' species loss under different compensation scenarios.
#' 

# load relevant libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

# load the winfree functions
source("code/helper-plotting-theme.R")
source("code/helper-miscellaneous.R")

# load the summarised transect data
tra_dat <- readRDS(file = here("output/transect_ssdb.rds"))
head(tra_dat)

# copy the tra_dat that we will simulate species extinction on
tra_ext <- tra_dat

# load the growth rate data
grow_dat_list <- readRDS(file = "output/model_growth_rates.rds")

# get a random sample of growth rates
grow_dat <- grow_dat_list[[900]]
head(grow_dat)

# choose a species to go extinct
sp_ext <- c("fu_se")

# choose a level of compensation
comp <- 0.10

# get the compensation level of standing biomass
comp_ext <- tra_dat[[sp_ext]]*comp

# set the biomass of the extinct species to zero
tra_ext[[sp_ext]] <- 0

# replace that biomass with one of the species that has positive growth in that zone
for(dz in 1:nrow(grow_dat)) {

  # get the different names that have positive growth rates in that zone
  x_grow <- grow_dat[, names(grow_dat) != sp_ext]
  x_grow <- x_grow[x_grow$depth == dz, -1]
  sp_grow <- names(x_grow)[ x_grow > 0 ]
  
  # get the names of species with positive standing biomass in an adjacent zone
  
  # first, we get suitable zones
  zones <- c(dz-1, dz, dz+1)
  zones <- zones[(zones > 0) & (zones < 5)]
  
  # second, we get suitable zones with species present
  x_ssb <- tra_ext[, names(tra_ext) != sp_ext ]
  x_ssb <- x_ssb[x_ssb$depth %in% zones, -1]
  sp_ssb <- names(x_ssb)[apply(x_ssb, 2, function(x) any(x > 0))]

  # get intersecting species i.e. positive growth and adjacent zone
  sp_comp <- intersect(sp_grow, sp_ssb)
  
  # if there is more than one suitable species, we sample one randomly
  if( length(sp_comp) > 0 ) {
    
    # sample one of those names
    sp_comp <- sample(sp_comp, 1)
    
    # the sp_comp is going to compensate for the species lost in the dz depth zone
    tra_ext[dz, ][[sp_comp]] <- (tra_ext[dz, ][[sp_comp]] + comp_ext[dz])
    
  }
  
}

# make a data list
dat_list <- list(tra_dat, tra_ext)
ylabs <- c("Standing dry biomass (g)", "")

# get an output list
fig_list <- vector("list", length = length(dat_list))

# plot the standing stock with and without the extinct species
for(i in 1:length(fig_list)) {
  
  # plot a figure of the depth distribution
  tra_plot <- 
    dat_list[[i]] %>%
    rename(`F. serratus` = "fu_se",
           `A. nodosum` = "as_no",
           `F. vesiculosus` = "fu_ve",
           `F. spiralis` = "fu_sp") %>%
    pivot_longer(cols = c(`F. serratus`, `A. nodosum`, `F. vesiculosus`, `F. spiralis`),
                 names_to = "species", 
                 values_to = "biomass")
  
  # change levels of the factor
  tra_plot$species <- factor(tra_plot$species,
                             c("F. spiralis", "F. vesiculosus", "A. nodosum", "F. serratus"))
  
  # change the levels of the depth factor
  tra_plot$depth <- factor(tra_plot$depth, levels = c("4", "3", "2", "1"))
  levels(tra_plot$depth) <- paste0("DZ", c("1", "2", "3", "4"))
  
  p1 <- 
    ggplot(data = tra_plot,
           mapping = aes(x = depth, y = biomass, fill = species)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8,
             alpha = 0.75, colour = "black", size = 0.1) +
    # scale_fill_viridis_d(option = "A", end = 0.9) +
    scale_fill_manual(values = seaweed_pal()) +
    ylab(ylabs[i]) +
    xlab(NULL) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 570)) +
    labs(fill = "Species") +
    geom_vline(xintercept = c(1.5, 2.5, 3.5), 
               linetype = "longdash", colour = "black", 
               size = 0.25) +
    theme_meta() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 12))
  
  fig_list[[i]] <- p1
  
}

p1 <- 
  fig_list[[1]] +
  ggtitle("Intact community") +
  theme(plot.title = element_text(size = 12, hjust = 0.5))
plot(p1)

p2 <- 
  fig_list[[2]] +
  ggtitle("F. serratus Extinction + (10% comp.)") +
  theme(plot.title = element_text(size = 12, hjust = 0.5))
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
  # scale_colour_viridis_d(option = "A", end = 0.9) +
  scale_colour_manual(values = seaweed_pal()) +
  scale_y_continuous(limits = c(-0.01, 0.02)) +
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
  
  lapply(dat_list, function(ssdb_dat) {
    
    # convert biomass to productivity units
    prod_x <- 
      mapply(function(ssdb, growth_data){
        
        # multiply the ssdb by the growth rate data
        ssdb*growth_data
        
        }, 
        ssdb_dat[, -1], grow_dat[, -1])
    
    # add a depth variable to the data.frame
    prod_x <- cbind(data.frame(depth = as.character(1:4)), as.data.frame(prod_x))
    
    return(prod_x) 
    
    })

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

ggsave(filename = "figures-tables/fig_3.png", p1234, dpi = 400,
       units = "cm", width = 20, height = 18)

### END
