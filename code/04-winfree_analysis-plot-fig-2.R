#' @title: perform the Winfree et al. (2018) analysis on the macroalgae data
#' 
#' @description: this script analyses the macroalgae data using Winfree et al.'s
#' (2018, Science) method to calculate the number of species required to maintain
#' ecosystem functioning across different sites
#' 

# load relevant libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

# load the winfree functions
source("code/helper-plotting-theme.R")
source("code/helper-miscellaneous.R")
source("code/helper-winfree-functions.R")

# load the summarised transect data
tra_dat <- readRDS(file = "output/transect_ssdb.rds")
head(tra_dat)

# load the growth rate data
grow_dat_list <- readRDS(file = "output/model_growth_rates.rds")

# calculate raw productivity for each sample of growth rates
prod_raw <- 
  
  lapply(grow_dat_list, function(grow_dat) {
  
  # convert biomass to productivity units
  prodx <- 
    mapply(function(ssdb, growth_rates){ 
    
    # multiply the standing stock dry biomass (ssdb) by the relative growth rates
    ssdb*growth_rates
    
    }, 
    tra_dat[,-1], grow_dat[,-1] )
  
  # add a depth zone variable back to the data
  prodx <- cbind(data.frame(depth_zone = as.character(1:4)), as.data.frame(prodx))
  
  return(prodx) 
  
  })

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

# calculate total productivity and summarise into the mean and percentile intervals
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
                width = 0, position = position_dodge(0.5), alpha = 0.7) +
  geom_point(data = prod,
             mapping = aes(x = depth, y = prod_m, colour = species), 
             size = 2.5, position = position_dodge(0.5), alpha = 0.9,
             shape = 16) +
  geom_point(data = prod_sum,
             mapping = aes(x = depth, y = total_prod_m), 
             size = 4.5, colour = "red", shape = 18) +
  geom_errorbar(data = prod_sum, 
                mapping = aes(x = depth, ymin = PI_low, ymax = PI_high),
                width = 0, colour = "red") +
  # scale_colour_viridis_d(option = "A", end = 0.9) +
  scale_colour_manual(values = seaweed_pal()) +
  xlab("") +
  ylab(expression("Dry biomass prod."~(g~day^{-1}) )) +
  theme_meta() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12))
plot(p1)

# save this as a .rds file
saveRDS(object = p1, file = "output/fig_1c.rds")


# Winfree et al. (2018) approach

# 1. for each sample of growth rates, calculate the productivity
# 2. use the Winfree functions to calculate number of species required

sp_pool_lev <- 
  
  lapply( c(0.25, 0.50, 0.75), function(func_lev) {
  
  sp_pool_list <- 
    
    lapply(grow_dat_list, function(grow_dat) {
      
      # convert biomass to productivity units
      prodx <- 
        mapply(function(ssdb, growth_rates){ 
        
          # multiply standing stock biomass by modelled growth rates
          ssdb*growth_rates 
          
        }, 
        tra_dat[,-1], grow_dat[,-1])
      
      # add a depth variable
      prodx <- cbind(data.frame(depth_zone = as.character(1:4)), as.data.frame(prodx))
      
      # calculate mean productivity across depth zones
      mean_prod <- mean( rowSums(prodx[,-1]) )
      
      # generate different combinations of site numbers
      perms <- gtools::permutations(n = 4, r = 4, v = prodx$depth_zone)
      l_perms <- vector("list", length = nrow(perms))
      for(i in 1:nrow(perms)) {
        l_perms[[i]] <- perms[i,]
      }
      
      # split species by site matrix into list
      sp_site <- split(prodx[,-1], prodx$depth_zone)
      
      # process the list so it has the correct format
      sp_site <- 
        
        lapply(sp_site, function(x) { 
        
        sapply(x, function(y) {
          
          y
          
          } ) 
          
        } )
      
      # calculate species required to reach 50% biomass
      sp_pool <- 
        
        lapply(l_perms, function(site_in) {
          
          df <- GA_optimiser(list_abundances = sp_site, 
                             func_threshold = func_lev*mean_prod, 
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
      
      # summarise this into the accumulation curves that Winfree et al. (2018) produce
      sp_pool <- 
        lapply(sp_pool, function(x) {
          
          # initialise first number of species
          y <- x[[1]]
          
          # set-up an empty vector equal to the number of site combinations
          n_sp <- vector(length = length(x))
          
          # initialise the empty vector with the number of 
          # unique species in first combination
          n_sp[1] <- length(unique(y))
          
          # run a loop to add the number of species required in the other site
          # and the first N sites
          for(i in 2:length(x)) {
            
            y <- c(y, x[[i]])
            n_sp[i] <- length(unique(y[!is.na(y)]))
            
          }
          
          # bring this output into a data.frame and add the number of sites
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
  sp_pool_df$threshold <- func_lev
  
  # return this data.frame
  return(sp_pool_df)
  
} )

# bind the Winfree et al. (2018) output into a data.frame
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
  ggtitle("Fucoid macroalgae data") +
  theme_meta() +
  theme(legend.key = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5, vjust = 2.5))
plot(p1)

ggsave(filename = "figures-tables/fig_2.png", p1, dpi = 400,
       units = "cm", width = 13, height = 8)

### END
