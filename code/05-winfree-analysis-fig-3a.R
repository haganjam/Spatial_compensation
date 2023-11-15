#' @title: perform the Winfree et al. (2018) analysis on the macroalgae data
#' 
#' @description: this script analyses the macroalgae data using Winfree et al.'s
#' (2018, Science) method to calculate the number of species required to maintain
#' ecosystem functioning across different sites
#' 

# clear environment
rm(list = ls())

# load relevant libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

# load the winfree functions
source("code/helper-plotting-theme.R")
source("code/helper-miscellaneous.R")
source("code/helper-winfree-functions.R")

# load the productivity data
prod_list <- readRDS("output/model_productivity.rds")

# Winfree et al. (2018) approach

# 1. for each sample of growth rates, calculate the productivity
# 2. use the Winfree functions to calculate number of species required

sp_pool_lev <- 
  
  lapply( c(0.25, 0.50, 0.75), function(func_lev) {
  
  sp_pool_list <- 
    
    lapply(prod_list, function(prodx) {
      
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
      sp_pool <- dplyr::bind_rows(sp_pool)
      
      return(sp_pool)
      
    } )
  
  # bind into a data.frame
  sp_pool_df <- dplyr::bind_rows(sp_pool_list, .id = "sample")
  
  # add the threshold as a variable
  sp_pool_df$threshold <- func_lev
  
  # return this data.frame
  return(sp_pool_df)
  
} )

# bind the Winfree et al. (2018) output into a data.frame
sp_pool_lev <- dplyr::bind_rows(sp_pool_lev)

# check the distribution of these different levels
ggplot(data = sp_pool_lev,
       mapping = aes(x = richness_required)) +
  geom_histogram() +
  facet_wrap(~n_sites, scales = "free")

# summarise into a mean and percentile interval
sp_pool_lev <- 
  sp_pool_lev |>
  dplyr::group_by(threshold, n_sites) |>
  dplyr::summarise(richness_required_m = mean(richness_required),
                   richness_required_sd = sd(richness_required), .groups = "drop") |>
  dplyr::mutate(n_sites = as.numeric(n_sites),
                Thresh. = as.character(threshold))
print(sp_pool_lev)

# get a colour palette
col_pal <- gray.colors(3, start = 0.2, end = 0.8, rev = TRUE)

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
  scale_colour_manual(values = col_pal) +
  guides(colour = guide_legend(override.aes = list(shape = 16, size = 4))) +
  ylab("Number of species required") +
  xlab("Number of depth zones") +
  scale_x_continuous(limits = c(0.8, 4.2)) +
  theme_meta() +
  theme(legend.position = "top",
        legend.key = element_rect(fill = NA),
        plot.title = element_text(hjust = 0.5, vjust = 2.5))
plot(p1)

# export a figure for the defence
pd <- 
  p1 + 
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent', colour = 'transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent', colour = 'transparent') #transparent legend panel
  )
plot(pd)
ggsave(filename = "figures-tables/fig_3_def.pdf", pd,
       units = "cm", width = 10, height = 8.5, bg = "transparent")

# output as a .rds object
saveRDS(p1, "output/fig_3a.rds")

### END
