#' @title - Model the growth rates of the individuals across the tile experiment
#' 
#' @description - this script uses allometric equations to convert transects
#' of four fucoid macroalgae species to standing stock biomass in four
#' different depth zones.
#' 

# load relevant libraries
library(dplyr)
library(tidyr)
library(readr)
library(rstan)
library(ggplot2)

# load relevant helper functions
source("code/helper-plotting-theme.R")
source("code/helper-miscellaneous.R")

# load the growth rate data
gr_dat <- read_csv("data/growth_rate_data.csv")
head(gr_dat)

# calculate the number of missing individuals
gr_dat %>%
  filter(site_code %in% c("X", "Y") ) %>%
  group_by(binomial_code) %>%
  summarise(n = n(),
            not_NA = sum(!is.na(dry_weight_g_daily_change)),
            n_NA = sum(is.na(dry_weight_g_daily_change)),
            )

# remove growth rate NAs
gr_dat <- 
  gr_dat %>%
  filter(!is.na(dry_weight_g_daily_change))

# check the sites
unique(gr_dat$site_code)

# get site X and Y which is the same site as the transects
gr_dat <- 
  gr_dat %>%
  filter(site_code %in% c("X", "Y") )

# compare the growth rates overall among species
ggplot(data = gr_dat,
       mapping = aes(x = depth_treatment, y = dry_weight_g_daily_change)) +
  geom_point() +
  facet_wrap(~binomial_code) +
  theme_classic()

# compile the stan growth rate model
m1 <- rstan::stan_model("code/02-growth-rate-model.stan")
print(m1)

# create a list with the relevant data
dat_sp <- list(species = as.integer(factor(gr_dat$binomial_code)),
               depth = as.integer(factor(gr_dat$depth_treatment)),
               growth = gr_dat$dry_weight_g_daily_change,
               N = nrow(gr_dat),
               D_N = length(unique(gr_dat$depth_treatment)),
               S_N = length(unique(gr_dat$binomial_code)))

# sample the stan model: m1
m1_fit <- rstan::sampling(m1, data = dat_sp, 
                          iter = 1000, chains = 4, algorithm = c("NUTS"),
                          control = list(adapt_delta = 0.99),
                          seed = 54856)

# check the traceplots and rhat values to assess convergence
print(m1_fit)

# export the traceplots

# get a vector of parameters
pars <- m1_fit@model_pars 
pars <- pars[-length(pars)]
pars <- pars[pars != "vbar_mat"]

# get the parameter dimensions
par_dims <- m1_fit@par_dims
par_dims <- par_dims[-length(par_dims)]
par_dims <- par_dims[names(par_dims) != "vbar_mat"]
par_dims <- lapply(par_dims, function(x) if(sum(x) == 0) { 1 } else {x} )
par_dims <- lapply(par_dims, function(x) if(length(x) == 1 && x == 4) {c(2, 2)} else {x} )

# loop over each parameter
for(i in pars) {
  
  # generate the traceplot
  TP <- traceplot(m1_fit, pars = i )
  
  # get the correct plot dimensions
  w <- par_dims[[i]][1]
  h <- if( length(par_dims[[i]]) > 1) {
    par_dims[[i]][2]
  } else {
    1
  }
  
  # sigma plot needs to be wider
  if(i == "sigma") {
    w <- 1.5
  }
  
  # export the traceplot
  ggsave(filename = paste0("figures-tables/", paste0("TP_",i, ".svg")), 
         TP, 
         units = "cm", width = (5*w), height = (4.5*h))
  
}

# generate diagnostic information output
m1_sum <- summary(m1_fit)

# convert to tibble and add parameter names
tab_s4 <- as_tibble(m1_sum$summary)
tab_s4 <- bind_cols(tibble(parameter = row.names(m1_sum$summary)),
                    tab_s4)

# select relevant columns
tab_s4 <- select(tab_s4, parameter, mean, sd, n_eff, Rhat)

# rename the columns
names(tab_s4) <- c("Parameter", "Mean", "SD", "Effective number of samples", "R-hat")

# remove the final row
tab_s4 <- tab_s4[-nrow(tab_s4),]

# remove the vbar_mat parameter as these are only for calculation purposes
tab_s4 <- filter(tab_s4, !grepl(pattern = "vbar_mat", Parameter) )

# export this table as a .csv file
write_csv(x = tab_s4, file = "figures-tables/table_S4.csv")

# extract the samples from the estimated posterior distribution
m1_post <- rstan::extract(m1_fit)

# check the fit of the model to the observed data

# extract the relevant parameters
alpha <- m1_post$a

# get predictions from the model for the data
pred_list <- vector("list", length = length(m1_post$sigma))
for (i in 1:length(m1_post$sigma)) {
  
  mu <- 
    
    sapply(1:length(dat_sp$species), function(x) {
    
    x <- alpha[, , dat_sp$depth[x]][i , dat_sp$species[x]]
    return(x)
    
  })
  
  pred_list[[i]] <- mu
  
}

# bind into a matrix
mu <- do.call("rbind", pred_list)

# get the mu value 
dat_sp$mu <- apply(mu, 2, mean)

# get the minimum and maximum values
dat_sp$PI_low <- apply(mu, 2, min)
dat_sp$PI_high <- apply(mu, 2, max)

# bind the data list with the predictions into a data.frame
obs <- bind_cols(dat_sp)
pred <- 
  obs %>%
  group_by(species, depth) %>%
  summarise(mu = first(mu), 
            PI_low = first(PI_low),
            PI_high = first(PI_high))

# plot the observed data versus the model predictions
ggplot() +
  geom_jitter(data = obs,
              mapping = aes(x = depth, y = growth),
              alpha = 0.25, shape = 16, size = 2, width = 0.1) +
  geom_point(data = pred,
             mapping = aes(x = depth, y = mu), colour = "red") +
  geom_line(data = pred,
            mapping = aes(x = depth, y = mu), colour = "red") +
  geom_errorbar(data = pred,
                mapping = aes(x = depth, ymin = PI_low, ymax = PI_high),
                width = 0, colour = "red") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~species, scales = "free") +
  theme_meta()

# calculate the r2 of this model
r2 <- 
  
  apply(mu, 1, function(x) {
    
    r <- x - obs$growth
    r <- 1 - (var2(r)/var2(obs$growth))
    
    return(r)
    
  } )

# calculate the summary statistics of the r2 value
mean(r2)
range(r2)

# plot the relationship between observed and predicted values with a one-to-one line
plot(obs$mu, obs$growth)
abline(a = 0, b = 1)


# simulate a posterior distribution for each species at each transect depth

# set-up a data.frame with the four depth levels crossed with the four species
df_pred <- expand.grid(depth = c(1:4), species = c(1:4))

# simulate the expected growth rates growth rates
pred_list <- vector("list", length = length(m1_post$sigma))
for (i in 1:length(m1_post$sigma)) {
  
  mu <- 
    
    sapply(1:length(df_pred$species), function(x) {
      
      x <- alpha[, , df_pred$depth[x]][i , df_pred$species[x]]
      return(x)
      
    })
  
  pred_list[[i]] <- mu
  
}

# bind into a matrix
sim_gr <- do.call("rbind", pred_list)

# add mu and percentile intervals onto the df_pred data.frame
# we divide by 100 to get back to the original growth-rate scale
df_pred$mu <- apply(sim_gr, 2, function(x) mean(x/100))
df_pred$PI_low <- apply(sim_gr, 2, function(x) PI(x/100, prob = 0.90)[1] )
df_pred$PI_high <- apply(sim_gr, 2, function(x) PI(x/100, prob = 0.90)[2] )

# convert the species variable into a factor to make sure the colours are equalised
df_pred$species <- factor(df_pred$species, levels = c("3", "4", "1", "2"))
levels(df_pred$species) <- c("fu_sp", "fu_ve", "as_no", "fu_se")


# plot the modeled growth rates

# get the raw data
gr_dat$depth <- as.integer(as.factor(gr_dat$depth_treatment))
gr_dat$growth <- (gr_dat$dry_weight_g_daily_change/100)
gr_dat$species <- factor(gr_dat$binomial_code, levels = c("fu_sp", "fu_ve", "as_no", "fu_se"))

# modify the depth factors

# df_pred
df_pred_plot <- df_pred
df_pred_plot$depth <- factor(df_pred_plot$depth, levels = c("4", "3", "2", "1"))
levels(df_pred_plot$depth) <- paste0("DZ", c("1", "2", "3", "4"))

# raw data gr_dat
gr_dat$depth <- factor(gr_dat$depth, levels = c("4", "3", "2", "1"))
levels(gr_dat$depth) <- paste0("DZ", c("1", "2", "3", "4"))

p1 <- 
  ggplot() +
  geom_point(data = gr_dat,
             mapping = aes(x = depth, y = growth, colour = species),
             shape = 1, alpha = 0.5, size = 1.2, 
             position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.25)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = df_pred_plot,
             mapping = aes(x = depth, y = mu, colour = species),
             position = position_dodge(0.25), size = 3, shape = 18) +
  geom_errorbar(data = df_pred_plot,
                mapping = aes(x = depth, ymin = PI_low, ymax = PI_high, colour = species),
                width = 0,
                position = position_dodge(0.25)) +
  geom_line(data = df_pred_plot,
            mapping = aes(x = as.integer(depth), y=mu, colour = species),
            position = position_dodge(0.25)) +
  # scale_colour_viridis_d(option = "A", end = 0.9) +
  scale_colour_manual(values = seaweed_pal()) +
  xlab("") +
  ylab(expression("Dry biomass change"~(g~g^{-1}~day^{-1}) )) +
  theme_meta() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12))
plot(p1)

# save the plot as a .rds file
saveRDS(object = p1, file = "output/fig_1b.rds")

# convert these growth rates into a list
sim_gr <- 
  
  apply(sim_gr, 1, function(x) {
  
  # bind the meta data and relevel the species
  x <- cbind(df_pred[,c(1, 2)], data.frame(growth = x/100))
  
  # manipulate the data.frame to the appropriate long format
  x <- 
    x %>%
    pivot_wider(id_cols = "depth",
                names_from = "species",
                values_from = "growth") %>%
    select(depth, fu_se, as_no, fu_ve, fu_sp)
  
  return(x)
  
} )

# check this list
sim_gr[[sample(1:length(sim_gr), 1)]]

# save this growth rate list as a .rds file
saveRDS(object = sim_gr, file = "output/model_growth_rates.rds")

### END
