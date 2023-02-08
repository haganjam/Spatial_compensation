
# Compensation analysis based on the algae data

# load relevant libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

# set the working directory
setwd("C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_2_Fucus_landscape/compensation_analysis")
getwd()

# load the plotting theme
source("scripts/function_plotting_theme.R")

# set the number of days to simulate growth over
n_days <- 1

# load the summarised transect data
tra_dat <- read_csv(file = "data/transect_comp.csv")
head(tra_dat)

# load the growth rate data
grow_dat_list <- readRDS(file = "data/growth_rate_data.rds")

# set the number of species
sp_vec <- c("fu_se", "as_no", "fu_ve", "fu_sp")

sp_out <- vector("list", length = length(sp_vec))
for(i in 1:length(sp_vec)) {
  
  # set the extinct species
  sp_ext <- sp_vec[i]
  
  # assume that fu_se and then different levels of standing biomass compensation
  prod.out <- 
    
    lapply(seq(0.1, 0.9, 0.1), function(comp) {
      
      list.out <- 
        
        lapply(grow_dat_list, function(grow_dat) {
          
          # get a copy of the data that can be modified
          df.tra <- tra_dat
          
          # calculate the intact community productivity using the growth rate data
          prodx <- mapply(function(x, y){x*y*n_days}, df.tra[,-1], grow_dat[,-1])
          prodx <- cbind(data.frame(depth_zone = as.character(1:4)), as.data.frame(prodx))
          
          # get the compensation level of standing biomass
          comp.x <- df.tra[[sp_ext]]*comp
          
          # set the biomass of the extinct species to zero
          df.tra[[sp_ext]] <- 0
          
          # for each depth zone for the extinct species
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
          
          # calculate the productivity using the growth rate data
          # convert biomass to productivity units
          prodx2 <- mapply(function(x, y){x*y*n_days}, df.tra[,-1], grow_dat[,-1])
          prodx2 <- cbind(data.frame(depth_zone = as.character(1:4)), as.data.frame(prodx2))
          
          # do we reach 50% of the original biomass at all sites?
          thresh_50 <- rowSums(prodx2[,-1]) > 0.50*rowSums(prodx[,-1])
          perc_diff <- ((rowSums(prodx2[,-1]) - rowSums(prodx[,-1]))/rowSums(prodx[,-1]))*100
          
          # wrap this into a useful data.frame
          df.out <- data.frame(depth_zone = as.character(1:4),
                               prod_init = rowSums(prodx[,-1]),
                               prod_comp = rowSums(prodx2[,-1]))
          
          return(df.out)
          
        } )
      
      # bind into a data.frame
      list.out <- bind_rows(list.out, .id = "sample_gr")
      list.out$comp_level <- comp
      
      return(list.out)
      
    } )
  
  # bind the rows
  prod.out <- bind_rows(prod.out)
  prod.out$species_ext <- sp_ext
  
  # write into the list
  sp_out[[i]] <- prod.out
  
}

# bind into one big data.frame
sp_out <- bind_rows(sp_out)
head(sp_out)

# reorder the columns
sp_out <- 
  sp_out %>%
  select(species_ext, comp_level, sample_gr, depth_zone, 
         prod_init, prod_comp)

# save the data as a .rds file
saveRDS(sp_out, file = "data/species_extinction_analysis.rds")

### END
