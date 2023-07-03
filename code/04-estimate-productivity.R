#' @title: calculate dry biomass productivity per depth zone
#' 
#' @description: this script uses the standing stock dry biomass data along with
#' the modelled growth rate data to estimate overall dry biomass productivity
#' per depth zone
#' 

# load relevant libraries
library(dplyr)
library(tidyr)
library(ggplot2)

# load the winfree functions
source("code/helper-plotting-theme.R")
source("code/helper-miscellaneous.R")

# load the summarised transect data
tra_dat <- readRDS(file = "output/transect_ssdb.rds")
head(tra_dat)

# load the growth rate data
grow_dat_list <- readRDS(file = "output/model_growth_rates.rds")
head(grow_dat_list[[1]])

# calculate raw productivity for each sample of growth rates
prod_list <- 
  
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

# save this as an output
saveRDS(object = prod_list, "output/model_productivity.rds")

# plot the raw productivity across depth zones integrated across growth rates
prod <- 
  bind_rows(prod_list, .id = "sample_gr") %>%
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
prod$depth <- factor(prod$depth_zone)
levels(prod$depth) <- paste0("DZ", c("1", "2", "3", "4"))

# calculate total productivity and summarise into the mean and percentile intervals
prod_sum <- 
  bind_rows(prod_list, .id = "sample_gr") %>%
  mutate(total_prod = (fu_se + as_no + fu_ve + fu_sp)) %>%
  group_by(depth_zone) %>%
  summarise(total_prod_m = mean(total_prod, na.rm = TRUE),
            PI_low = quantile(total_prod, 0.05),
            PI_high = quantile(total_prod, 0.95))

# change the factor levels of depth
prod_sum$depth <- factor(prod_sum$depth_zone)
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

### END
