
# Model the growth rates of the individuals across the tile experiment

# load relevant libraries
library(dplyr)
library(tidyr)
library(readr)
library(readr)
library(rethinking)
library(ggplot2)

# set the working directory
setwd("C:/Users/james/OneDrive/PhD_Gothenburg/Chapter_2_Fucus_landscape/compensation_analysis")
getwd()

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

# get site Y which is the same site as the transects
gr_dat <- 
  gr_dat %>%
  filter(site_code == "Y")

# compare the growth rates overall among species
ggplot(data = gr_dat,
       mapping = aes(x = binomial_code, dry_weight_g_daily_change)) +
  geom_boxplot()

# model the growth rates as a function of depth and species
gr_list <- list(depth = gr_dat$depth_treatment,
                species = as.integer(as.factor(gr_dat$binomial_code)),
                growth = gr_dat$dry_weight_g_daily_change)

# fit a model to these data
m1 <- 
  ulam(
    
    alist(growth ~ dnorm(mu, sigma),
          mu <- a[species] + b1[species]*depth,
          
          a[species] ~ dnorm(0, 2),
          b1[species] ~ dnorm(0, 2),
          sigma ~ dexp(1) ),
       
       data = gr_list, chains = 4)
precis(m1, depth = 2)
traceplot(m1)
dev.off()

# check the predictions
pred.m1 <- link(fit = m1)

# calculate the mean and percentile interval
mu <- apply(pred.m1, 2, mean)
PI <- apply(pred.m1, 2, PI)

# check fit to sample
plot(gr_list$growth, mu)

# add the predictions to the list
df <- bind_cols(gr_list)
df$growth_pred <- mu

# check relationship within species
ggplot(data = df,
       mapping = aes(x = growth, y = growth_pred, colour = depth)) +
  geom_point() +
  facet_wrap(~species, scales = "free") +
  geom_abline(intercept = 0, slope = 1) +
  theme_classic()

# simulate a posterior distribution for each species at each transect depth
df.pred <- expand.grid(depth = c(-44, -32, -20, -8),
                       species = c(1:4))

# simulate the growth rates
sim.gr <- link(fit = m1, data = df.pred)

# add mu and percentile intervals onto the df.pred data.frame
df.pred$mu <- apply(sim.gr, 2, function(x) mean(x/100))
df.pred$PI_low <- apply(sim.gr, 2, function(x) PI(x/100, prob = 0.90)[1] )
df.pred$PI_high <- apply(sim.gr, 2, function(x) PI(x/100, prob = 0.90)[2] )

# convert the species variable into a factor to make the colours are equalised
df.pred$species <- factor(df.pred$species)
levels(df.pred$species) <-  list("as_no" = "1", "fu_se" = "2", "fu_sp" = "3", "fu_ve" = "4")

# reorder the levels
df.pred$species  <- factor(df.pred$species, levels = c("fu_se", "as_no", "fu_ve",  "fu_sp"))

# add a depth zone variable
df.pred$depth_zone <- (rep(1:4, 4))

# plot the modeled growth rates
fig_1b <- 
  ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = df.pred,
             mapping = aes(x = depth_zone, y = mu, colour = species),
             position = position_dodge(0.25), size = 2) +
  geom_errorbar(data = df.pred,
                mapping = aes(x = depth_zone, ymin = PI_low, ymax = PI_high, colour = species),
                width = 0,
                position = position_dodge(0.25)) +
  geom_line(data = df.pred,
            mapping = aes(x = depth_zone, y=mu, colour = species),
            position = position_dodge(0.25)) +
  scale_colour_viridis_d(option = "A", end = 0.9) +
  xlab("Depth zone") +
  ylab(expression("Dry biomass change"~(g~g^{-1}~day^{-1}) )) +
  theme_meta() +
  theme(legend.position = "none")
plot(fig_1b)

saveRDS(object = fig_1b, file = "figures/Fig_1b.rds")


# convert these growth rates into a list
sim.gr <- 
  
  apply(sim.gr, 1, function(x) {
  
  # bind the meta data and relevel the species
  x <- cbind(df.pred[,c(1, 2)], data.frame(growth = x/100))
  x$species <- factor(x$species)
  levels(x$species) <-  list(as_no = "1", fu_se = "2", fu_sp = "3", fu_ve = "4")
  
  # change the depth factors
  x$depth <- as.integer(as.factor(x$depth))
  
  # manipulate the data.frame to the appropriate format
  x <- 
    x %>%
    pivot_wider(id_cols = "depth",
                names_from = "species",
                values_from = "growth") %>%
    select(depth, fu_se, as_no, fu_ve, fu_sp) %>%
    rename(depth_zone = depth)
  
  return(x)
  
} )

# save this growth rate list as a .rds file
saveRDS(object = sim.gr, file = "data/growth_rate_data.rds")

### END
