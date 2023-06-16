
# how much compensation would there need to be to reverse this accumulation
# analysis?

# if we get one of these graphs and it says that 15 species are required to
# maintain 50% functioning in all of 10 sites, then what is the counterfactual?

# what I want to know is then given that set of species, how much compensation
# is required to maintain 50% functioning at all sites given the extinction of 
# a species

# load the winfree functions
source("code/helper-winfree-functions.R")

# simulate an example of three species and three sites
# calculate species required to reach 50% biomass

eg <- 
  list(p1 = c("sp1" = 10, "sp2" = 1, "sp3" = 1),
       p2 = c("sp1" = 1, "sp2" = 10, "sp3" = 1),
       p3 = c("sp1" = 1, "sp2" = 1, "sp3" = 10)
       )
  
# get the mean functioning across sites
mean_func <- mean(unlist(lapply(eg, sum)))
  
# get the species pool for each site
sp_pool <- GA_optimiser(list_abundances = eg, 
                        func_threshold = 0.5*mean_func, 
                        site_vector = c("p1", "p2", "p3"))

# summarise this into the accumulation curves that Winfree et al. (2018) produce
    
# initialise first number of species
y <- sp_pool[[1]]

# set-up an empty vector equal to the number of site combinations
n_sp <- vector(length = length(sp_pool))

# initialise the empty vector with the number of 
# unique species in first combination
n_sp[1] <- length(unique(y))

# run a loop to add the number of species required in the other site
# and the first N sites
for(i in 2:length(sp_pool)) {
  
  y <- c(y, sp_pool[[i]])
  n_sp[i] <- length(unique(y[!is.na(y)]))
  
}

# bring this output into a data.frame and add the number of sites
df <- data.frame(n_sites = as.character(c(1:3)),
                 richness_required = n_sp)
print(df)

# imagine that species 3 went extinct, how many communities would not reach the
# 50% threshold?

# in this case, let's set species three as extinct
eg_ext <- 
  lapply(eg, function(x) { 
  x[names(x) == "sp3"] <- 0
  return(x)
  } )

# how many 
sum(unlist(lapply(eg_ext, sum)) >= 6)/length(eg)

# to get to the average, these remaining species would need to 
# both triple their pollination
# we can ask, is this a credible? Could thjs feasibly happen?
eg_ext


