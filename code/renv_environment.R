#' @title - set-up a reproducible package environment using renv
#' 
#' @description - this scripts records the code used to create a local package dependency
#' structure using the renv package. The 'cmdstanr' and 'rethinking' packages were ignored
#' because they are not on CRAN and do not work well with renv
#' 

# view currently-ignored packaged
renv::settings$ignored.packages()

# ignore cmdstanr and the rethinking package because these didn't work with renv
renv::settings$ignored.packages(c("cmdstanr", "rethinking"), persist = FALSE)

# initialise the renv local library
renv::init()

# check the status
renv::status()

### END
