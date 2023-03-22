
# write functions to implement the Winfree et al. (2018) approach

# important note:
# in the Winfree et al. (2018) approach, even if the number of species
# does not meet the mean * 0.50 threshold when all species have been 
# added, they are considered required

# write a function to get a vector of required species
site_species_required <- function(sp_abundances, function_level) {
  
  if( length(sp_abundances) == 0 ) {
    return(c(NA))
  } else {
    
    sp_abundances_sort <- sort(sp_abundances, decreasing = TRUE)
    
    cum_sp_abundances <- cumsum(sp_abundances_sort)
    sp_needed <- vector()
    
    for(i in 1:length(cum_sp_abundances)){
      
      if (cum_sp_abundances[i] <= function_level) {
        sp_needed[i] <- names(cum_sp_abundances[i])
      }
      
      if(cum_sp_abundances[i] > function_level) {
        sp_needed[i] <- names(cum_sp_abundances[i])
        break
      }
    }
    
    return(sp_needed)
    
  }
  
}


# GA optimiser function
# list_abundances: named list (by site) of vectors of named species abundances
# func_threshold: threshold of minimum function
# site_vector: vector of site names

GA_optimiser <- function(list_abundances, func_threshold, site_vector) {
  
  # initialise the loop
  start_abun <- list_abundances[[site_vector[1]]]
  
  # set-up an output vector
  sp_required <- vector("list")
  
  # fill the output vector with the first set of species
  sp_required[[1]] <- site_species_required(sp_abundances = start_abun, function_level = func_threshold)
  
  # loop over all other sites and fill sequentially
  for(i in 2:length(site_vector)) {
    
    sp_abun <- list_abundances[[site_vector[i]]]
    
    sp_indices <- which( names(sp_abun) %in% unlist(sp_required[1:i] ) )
    
    if( sum(sp_abun[sp_indices]) <= func_threshold ) {
      
      sp_abun_subset <- sp_abun[-(sp_indices)]
      
      new_sp_required <- site_species_required(sp_abundances = sp_abun_subset, function_level = func_threshold)
      
      sp_required[[i]] <- new_sp_required
      
    } else { sp_required[[i]] <- NA }
    
  }
  
  return(sp_required)
  
}