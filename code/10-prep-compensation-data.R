#' @title: calculate dry biomass productivity with and without each species
#' 
#' @description: this script calculates the dry biomass productivity of
#' each depth zone assuming that each species has been lost individually.
#' 

# clear environment
rm(list = ls())

# set a seed
set.seed(490527)

# load relevant libraries
library(dplyr)

# load the plotting theme
source("code/helper-plotting-theme.R")
source("code/helper-miscellaneous.R")

# load the summarised transect data
tra_dat <- readRDS("output/transect_ssdb.rds" )
head(tra_dat)

# load the growth rate data
grow_dat_list <- readRDS(file = "output/model_growth_rates.rds")

# set the number of species
sp_vec <- c("fu_se", "as_no", "fu_ve", "fu_sp")

# iterate over each species
sp_out <- vector("list", length = length(sp_vec))
for(i in 1:length(sp_vec)) {
  
  # set the extinct species
  sp_ext <- sp_vec[i]
  
  # iterate over 9 levels of potential, assumed compensation
  prod <- 
    lapply(seq(0.1, 0.9, 0.1), function(comp_level) {
      
      # iterate over each set of growth rates from the posterior distribution
      comp <- 
        lapply(grow_dat_list, function(grow_dat) {
          
          # get a copy of the data that can be modified
          tra_ext <- tra_dat
          
          # calculate the intact community productivity using the growth rate data
          prod_x <- 
            mapply(function(ssdb, growth_data){ 
              
              # multiply the standing stock biomass and the growth rate sample
              ssdb*growth_data 
              
              }, 
              tra_ext[, -1], grow_dat[, -1])
          
          # add a depth zone column
          prod_x <- cbind(data.frame(depth = as.character(1:4)), as.data.frame(prod_x))
          
          # get the compensation level of standing biomass
          comp_ext <- tra_ext[[sp_ext]]*comp_level
          
          # set the biomass of the extinct species to zero
          tra_ext[[sp_ext]] <- 0
          
          # for each depth zone for the extinct species
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
            sp_comp <- dplyr::intersect(sp_grow, sp_ssb)
            
            # if there is more than one suitable species, we sample one randomly
            if( length(sp_comp) > 0 ) {
              
              # sample one of those names
              sp_comp <- sample(sp_comp, 1)
              
              # the sp_comp is going to compensate for the species lost in the dz depth zone
              tra_ext[dz, ][[sp_comp]] <- (tra_ext[dz, ][[sp_comp]] + comp_ext[dz])
              
            }
            
          }
          
          # calculate the productivity using the growth rate data
          # convert biomass to productivity units
          prod_x2 <- 
            mapply(function(ssdb, growth_data){ 
              
              # multiply the new ssdb by the growth rate data
              ssdb*growth_data 
              
              }, 
              tra_ext[,-1], grow_dat[,-1])
          
          # add a depth variable
          prod_x2 <- cbind(data.frame(depth = as.character(1:4)), as.data.frame(prod_x2))
          
          # wrap this into a useful data.frame
          df <- data.frame(depth = as.character(1:4),
                           prod_init = rowSums(prod_x[,-1]),
                           prod_comp = rowSums(prod_x2[,-1]))
          
          return(df)
          
        } )
      
      # bind into a data.frame
      comp <- dplyr::bind_rows(comp, .id = "sample_gr")
      comp$comp_level <- comp_level
      
      return(comp)
      
    } )
  
  # bind the rows
  prod <- dplyr::bind_rows(prod)
  prod$species_ext <- sp_ext
  
  # write into the list
  sp_out[[i]] <- prod
  
}

# bind into one big data.frame
sp_out <- dplyr::bind_rows(sp_out)
head(sp_out)

# reorder the columns
sp_out <- 
  sp_out |>
  dplyr::select(species_ext, comp_level, sample_gr, depth, prod_init, prod_comp)

# save the data as a .rds file
saveRDS(sp_out, file = "output/compensation_analysis_data.rds")

### END
