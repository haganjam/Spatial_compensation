
# Perform the Winfree et al. (2018) analysis on the macroalgae data

# load relevant libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(here)

# load the winfree functions
source("scripts/function2_winfree_2018.R")
source("scripts/function1_plotting_theme.R")

# load the summarised transect data
tra_dat <- readRDS(file = "data/transect_ssdb.rds")
head(tra_dat)

# load the growth rate data
grow_dat_list <- readRDS(file = "data/growth_rate_data.rds")

# calculate raw productivity for each sample of growth rates
prod_raw <- 
  lapply(grow_dat_list, function(grow_dat) {
  
  # convert biomass to productivity units
  prodx <- mapply(function(x, y){x*y}, tra_dat[,-1], grow_dat[,-1])
  prodx <- cbind(data.frame(depth_zone = as.character(1:4)), as.data.frame(prodx))
  return(prodx) }
  )

# plot the raw productivity across depth zones integrated across growth rates
prod <- 
  bind_rows(prod_raw, .id = "sample_gr") %>%
  pivot_longer(cols = c("fu_se", "as_no", "fu_ve", "fu_sp"),
               names_to = "species",
               values_to = "prod") %>%
  group_by(depth_zone, species) %>%
  summarise(prod_m = mean(prod, na.rm = TRUE),
            PI_low = quantile(prod, 0.05),
            PI_high = quantile(prod, 0.95))

# change factor levels of species
prod$species <- factor(prod$species, levels = c("fu_sp", "fu_ve", "as_no", "fu_se"))

# change the factor levels of depth
prod$depth <- factor(prod$depth_zone, levels = c("4", "3", "2", "1"))
levels(prod$depth) <- paste0("DZ", c("1", "2", "3", "4"))

prod_sum <- 
  bind_rows(prod_raw, .id = "sample_gr") %>%
  mutate(total_prod = (fu_se + as_no + fu_ve + fu_sp)) %>%
  group_by(depth_zone) %>%
  summarise(total_prod_m = mean(total_prod, na.rm = TRUE),
            PI_low = quantile(total_prod, 0.05),
            PI_high = quantile(total_prod, 0.95))

# change the factor levels of depth
prod_sum$depth <- factor(prod_sum$depth_zone, levels = c("4", "3", "2", "1"))
levels(prod_sum$depth) <- paste0("DZ", c("1", "2", "3", "4"))

p1 <- 
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_errorbar(data = prod, 
                mapping = aes(x = depth, ymin = PI_low, ymax = PI_high, colour = species),
                width = 0, position = position_dodge(0.5), alpha = 0.4) +
  geom_point(data = prod,
             mapping = aes(x = depth, y = prod_m, colour = species), 
             size = 2, position = position_dodge(0.5), alpha = 0.7,
             shape = 16) +
  geom_point(data = prod_sum,
             mapping = aes(x = depth, y = total_prod_m), 
             size = 3.5, colour = "red", shape = 18) +
  geom_errorbar(data = prod_sum, 
                mapping = aes(x = depth, ymin = PI_low, ymax = PI_high),
                width = 0, colour = "red") +
  scale_colour_viridis_d(option = "A", end = 0.9) +
  xlab("") +
  ylab(expression("Dry biomass prod."~(g~day^{-1}) )) +
  theme_meta() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12))
plot(p1)

saveRDS(object = p1, file = "figures/fig1c.rds")


# Winfree et al. (2018) approach"

# 1. for each sample of growth rates, calculate the productivity
# 2. use the Winfree functions to calculate number of species required

sp_pool_lev <- 
  
  lapply(c(0.25, 0.50, 0.75), function(x) {
  
  sp_pool_list <- 
    
    lapply(grow_dat_list, function(grow_dat) {
      
      # convert biomass to productivity units
      prodx <- mapply(function(x, y){x*y}, tra_dat[,-1], grow_dat[,-1])
      prodx <- cbind(data.frame(depth_zone = as.character(1:4)), as.data.frame(prodx))
      
      # calculate mean productivity
      mean_prod <- mean( rowSums(prodx[,-1]) )
      
      # generate different combinations of site numbers
      perms <- gtools::permutations(n = 4, r = 4, v = prodx$depth_zone)
      l.perms <- vector("list", length = nrow(perms))
      for(i in 1:nrow(perms)) {
        l.perms[[i]] <- perms[i,]
      }
      
      # split species by site matrix into list
      sp_site <- split(prodx[,-1], prodx$depth_zone)
      sp_site <- lapply(sp_site, function(x) {sapply(x, function(y) {y} )} )
      
      # calculate species required to reach 50% biomass
      sp_pool <- 
        
        lapply(l.perms, function(site_in) {
          
          df <- GA_optimiser(list_abundances = sp_site, 
                             func_threshold = x*mean_prod, 
                             site_vector = site_in
          )
          
          return(df)
          
        })
      
      # what does this mean?
      
      # if for a given permutation of 1:4, fu_se is alone first it means that
      # fu_se is enough to provide 75% of the biomass
      
      # if the second row is NA, it means that fu_se is still enough
      
      # if then the third is not NA, then those species are required as well to provide
      # 50% functioning
      
      # summarise this into the curves that Winfree et al. (2018) produce
      sp_pool <- 
        lapply(sp_pool, function(x) {
          
          y <- x[[1]]
          n_sp <- vector(length = length(x))
          n_sp[1] <- length(unique(y))
          for(i in 2:length(x)) {
            
            y <- c(y, x[[i]])
            n_sp[i] <- length(unique(y[!is.na(y)]))
            
          }
          
          df <- data.frame(n_sites = as.character(c(1:4)),
                           richness_required = n_sp)
          
          return(df)
          
        } )
      
      # bind the rows into a data.frame and summarise
      sp_pool <- bind_rows(sp_pool)
      
      return(sp_pool)
      
    } )
  
  # bind into a data.frame
  sp_pool_df <- bind_rows(sp_pool_list, .id = "sample")
  
  # add the threshold as a variable
  sp_pool_df$threshold <- x
  
  # return this data.frame
  return(sp_pool_df)
  
} )

# bind into a data.frame
sp_pool_lev <- bind_rows(sp_pool_lev)

# check the distribution of these different levels
ggplot(data = sp_pool_lev,
       mapping = aes(x = richness_required)) +
  geom_histogram() +
  facet_wrap(~n_sites, scales = "free")

# summarise into a mean and percentile interval
sp_pool_lev <- 
  sp_pool_lev %>%
  group_by(threshold, n_sites) %>%
  summarise(richness_required_m = mean(richness_required),
            richness_required_sd = sd(richness_required), .groups = "drop") %>%
  mutate(n_sites = as.numeric(n_sites),
         Thresh. = as.character(threshold))
print(sp_pool_lev)

# plot number of sites and richness required
p1 <- 
  ggplot(data = sp_pool_lev) +
  geom_line(mapping = aes(x = n_sites, y = richness_required_m, colour = Thresh.), 
            size = 0.5, show.legend = FALSE, position = position_dodge(0.25)) +
  geom_point(mapping = aes(x = n_sites, y = richness_required_m, colour = Thresh.),
             position = position_dodge(0.25), size = 2) +
  geom_errorbar(mapping = aes(x = n_sites, 
                              ymin = richness_required_m-richness_required_sd,
                              ymax = richness_required_m+richness_required_sd,
                              colour = Thresh.),
                width = 0, show.legend = FALSE, position = position_dodge(0.25)) +
  scale_colour_viridis_d(option = "E", end = 0.9) +
  guides(colour = guide_legend(override.aes = list(shape = 16, size = 4))) +
  ylab("Number of species required") +
  xlab("Number of depth zones") +
  scale_x_continuous(limits = c(0.8, 4.2)) +
  theme_meta() +
  theme(legend.key = element_rect(fill = NA))
plot(p1)

ggsave(filename = "figures/fig2.png", p1, dpi = 400,
       units = "cm", width = 13, height = 8)

### END
