
# Model the growth rates of the individuals across the tile experiment

# load relevant libraries
library(dplyr)
library(tidyr)
library(readr)
library(readr)
library(rethinking)
library(ggplot2)

# load the plotting theme
source("scripts/function_plotting_theme.R")

# load the growth rate data
gr_dat <- read_csv("data/compensation_data.csv")
head(gr_dat)

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
       mapping = aes(x = binomial_code, dry_weight_g_daily_change)) +
  geom_boxplot() +
  theme_classic()

# relationship between growth and depth among species
ggplot(data = gr_dat,
       mapping = aes(x = depth_treatment, y = dry_weight_g_daily_change)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~binomial_code) +
  theme_classic()

# create a list with the relevant data

# subset fucus vesiculosus
dat_sp <- list(species = as.integer(factor(gr_dat$binomial_code)),
               depth = as.integer(factor(gr_dat$depth_treatment)),
               growth = gr_dat$dry_weight_g_daily_change)

# fit a model with partial pooling by across depths within each species
m1 <- ulam(
  
  alist(
    
    growth ~ dnorm(mu, sigma),
    
    mu <- alpha[species, depth],
    
    # adaptive priors
    vector[4]:alpha[species] ~ multi_normal(0, Rho_species, sigma_species),
    
    # fixed priors
    sigma_species ~ dexp(1),
    Rho_species ~ dlkjcorr(4),
    sigma ~ dexp(1)
    
  ) , data = dat_sp, chains = 4 , cores = 4 )

# check the precis output and the trace-plots

# note the diagonals in the correlation matrix should be NaN because they are constants
precis(m1, depth = 3)
traceplot(m1)

# check the fit of the model to the data

# extract the predictions from this model
mu <- link(m1)

# get the mu value 
dat_sp$mu <- apply(mu, 2, mean)

# get the PI value
dat_sp$PI_low <- apply(mu, 2, min)
dat_sp$PI_high <- apply(mu, 2, max)

x <- bind_cols(dat_sp)
y <- 
  x %>%
  group_by(species, depth) %>%
  summarise(mu = first(mu), 
            PI_low = first(PI_low),
            PI_high = first(PI_high))

ggplot() +
  geom_jitter(data = x,
              mapping = aes(x = depth, y = growth),
              alpha = 0.25, shape = 16, size = 2, width = 0.1) +
  geom_point(data = y,
             mapping = aes(x = depth, y = mu), colour = "red") +
  geom_line(data = y,
            mapping = aes(x = depth, y = mu), colour = "red") +
  geom_errorbar(data = y,
                mapping = aes(x = depth, ymin = PI_low, ymax = PI_high),
                width = 0, colour = "red") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~species, scales = "free") +
  theme_meta()

# we just treat the depth zones as discrete variables instead of predicting to this new data

# simulate a posterior distribution for each species at each transect depth
df.pred <- expand.grid(depth = c(1:4), species = c(1:4))

# simulate the growth rates
sim.gr <- link(fit = m1, data = df.pred)

# add mu and percentile intervals onto the df.pred data.frame
df.pred$mu <- apply(sim.gr, 2, function(x) mean(x/100))
df.pred$PI_low <- apply(sim.gr, 2, function(x) PI(x/100, prob = 0.90)[1] )
df.pred$PI_high <- apply(sim.gr, 2, function(x) PI(x/100, prob = 0.90)[2] )

# convert the species variable into a factor to make the colours are equalised
df.pred$species <- factor(df.pred$species, levels = c("2", "1", "4", "3"))
levels(df.pred$species) <- c("fu_se", "as_no", "fu_ve", "fu_sp")

# plot the modeled growth rates

# get the raw data
gr_dat$depth <- as.integer(as.factor(gr_dat$depth_treatment))
gr_dat$growth <- (gr_dat$dry_weight_g_daily_change/100)
gr_dat$species <- factor(gr_dat$binomial_code, levels = c("fu_se", "as_no", "fu_ve", "fu_sp"))

p1 <- 
  ggplot() +
  geom_point(data = gr_dat,
             mapping = aes(x = depth, y = growth, colour = species),
             shape = 1, alpha = 0.3, size = 2, 
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.25)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = df.pred,
             mapping = aes(x = depth, y = mu, colour = species),
             position = position_dodge(0.25), size = 2) +
  geom_errorbar(data = df.pred,
                mapping = aes(x = depth, ymin = PI_low, ymax = PI_high, colour = species),
                width = 0,
                position = position_dodge(0.25)) +
  geom_line(data = df.pred,
            mapping = aes(x = depth, y=mu, colour = species),
            position = position_dodge(0.25)) +
  scale_colour_viridis_d(option = "A", end = 0.9) +
  xlab("Depth zone") +
  ylab(expression("Dry biomass change"~(g~g^{-1}~day^{-1}) )) +
  theme_meta() +
  theme(legend.position = "none")
plot(p1)

saveRDS(object = p1, file = "figures/fig1b.rds")

# convert these growth rates into a list
sim.gr <- 
  
  apply(sim.gr, 1, function(x) {
  
  # bind the meta data and relevel the species
  x <- cbind(df.pred[,c(1, 2)], data.frame(growth = x/100))
  
  # manipulate the data.frame to the appropriate format
  x <- 
    x %>%
    pivot_wider(id_cols = "depth",
                names_from = "species",
                values_from = "growth") %>%
    select(depth, fu_se, as_no, fu_ve, fu_sp)
  
  return(x)
  
} )

# check this list
sim.gr[[sample(1:length(sim.gr), 1)]]

# save this growth rate list as a .rds file
saveRDS(object = sim.gr, file = "data/growth_rate_data.rds")

### END
