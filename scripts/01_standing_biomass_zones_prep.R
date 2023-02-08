
# compensation analysis test

# load relevant libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(here)

# load the plotting theme
source("scripts/function_plotting_theme.R")

# import the raw transect data
tra_dat <- read_csv(file = "data/transect_data.csv",
                    col_types = c("ccccccccnnncnncc"),
                    na = c("NA"))

# check the names of the transect data
names(tra_dat)

# 1. subset the data needed to interpolate the depth to the different points
# 2. correct the depth by the water level from the nearby station
depth_data <- 
  tra_dat %>% 
  mutate(depth_correct = (water_level_cm + depth_cm) ) %>%
  select(date, transect_id, position, depth_correct) %>%
  distinct()

# problems with the start of transect 2, remove positions 0 to 4 because we do not have a starting depth
depth_data <- 
  depth_data %>%
  filter( !(transect_id == 2 & position %in% 0:4) )

# split into a list
depth_list <- split(depth_data, depth_data$transect_id)

# loop over all transects
depth_out <- vector("list", length = length(depth_list))
for(i in 1:length(depth_list)) {
  
  # initialise a data.frame to work with
  df <- depth_list[[i]]
  
  # get the dividers
  dividers <- which(!is.na(df$depth_correct) )
  
  # duplicate middles for which the data are needed for multiple calculations
  if (length(dividers) > 2) {
    
    dups <- dividers[-c(1, length(dividers))]
    
    df <- 
      df[c(1:nrow(df), dups), ] %>%
      arrange(date, transect_id, position)
    
  }
  
  # add an ID column
  x <- vector("list", length = (length(dividers)-1))
  for(j in 1:(length(dividers)-1) ) {
    
    x[[j]] <- rep(j, ( (dividers[j+1] - dividers[j])+1 ) )
    
  }
  
  # write the ID column into the data.frame
  df$position_id <- unlist(x)
  
  # split data by the position ID
  df_list <- split(df, df$position_id)
  
  # for each block, interpolate the depth between the points
  int_depth <- 
    lapply(df_list, function(y){
      
      y %>%
        group_by(position_id) %>%
        mutate(y_int = first(depth_correct),
               d_depth = (last(depth_correct) - first(depth_correct)),
               d_position = last(10*position) - first(10*position)) %>%
        ungroup() %>%
        mutate(slope = d_depth/d_position) %>%
        select(-d_depth, -d_position) %>%
        mutate(distance = 0:(length(position)-1)*10 ) %>%
        mutate(depth_interpolated = y_int + ((distance)*slope) )  %>%
        select(-distance) %>%
        mutate(depth_interpolated = if_else(!is.na(depth_correct), depth_correct, depth_interpolated)) %>%
        select(-y_int, -slope, -position_id)
      
    }) %>%
    bind_rows(.,) %>%
    distinct()
  
  depth_out[[i]] <- int_depth
  
} 

# bind the list into a data.frame
depth_out <- bind_rows(depth_out)

# join this back to the full dataset and reorder the columns
tra_a <- 
  full_join(tra_dat, depth_out, by = c("date", "transect_id", "position")) %>%
  select(date, transect_id, site_code, time, position, water_level_cm, 
         depth_cm, depth_correct, depth_interpolated, binomial_code, length_cm, circum_cm, 
         field_observer, notes)

# check the summary statistics
summary(tra_a)

# check for unique values, especially for the binomial codes
lapply(tra_a, function(x) unique(x))

# remove missing binomial codes and the missing values "", NA
tra_a <- 
  tra_a %>%
  filter( !(binomial_code %in% c("-9999", "") | is.na(binomial_code) | is.na(depth_interpolated) )  )
names(tra_a)

# reorganise the column names
tra_a <- 
  tra_a %>%
  select(site_code, date, transect_id, time, position,
         depth_interpolated, binomial_code, length_cm, circum_cm) %>%
  rename(depth_RH2000_cm = depth_interpolated)

# check the depth distribution
hist(tra_a$depth_RH2000_cm)

# add a depth zone variable

# check the depths of the cut variable
cut(tra_a$depth_RH2000_cm, breaks = 4)

# midpoints of these depth zones are: -44, -32, -20, -8
tra_a$depth_zone <- cut(tra_a$depth_RH2000_cm, breaks = 4, labels = 1:4)

# load the allometry data
allo <- read_csv(file = "data/sample_data_biomass_allometry.csv")
head(allo)

lm.allo <- function(data, slope, e.vars) {
  
  # set an output list for the model coefficients
  est.lm <- vector("list", length(e.vars))
  names(est.lm) <- seq(1:length(e.vars))
  
  # set an output list for the model fit statistics
  fit.lm <- vector("list", length(e.vars))
  names(fit.lm) <- seq_along(1:length(e.vars))
  
  for (i in 1:length(e.vars) ) {
    
    # fit model using chosen predictors
    lm.e.vars <- lm(reformulate(e.vars[[i]], slope), data = data)
    
    # write coefficients to the est.lm list
    est.lm[[i]] <- broom::tidy(lm.e.vars)
    
    # write fit statistics to the fit.lm list
    fit.lm[[i]] <- broom::glance(lm.e.vars)
  }
  
  # convert lists to data.frames and join
  full_join(bind_rows(est.lm, .id = "model"), 
            bind_rows(fit.lm, .id = "model"),
            by = "model")
  
}

# set up the explanatory variables for the different models
exp.vars <- list(c("log_length_cm*log_circum_cm"),
                 c("log_length_cm", "log_circum_cm", "log_circum_cm2"),
                 c("log_length_cm", "log_length_cm", "log_length_cm2"),
                 c("log_length_cm","log_circum_cm"),
                 c("log_length_cm", "log_length_cm2"),
                 c("log_circum_cm", "log_circum_cm2"),
                 c("log_length_cm"),
                 c("log_circum_cm2"),
                 "1")

# add model variables to the data
allo <- 
  allo %>%
  mutate(log_dry_weight_g = log(dry_weight_g),
         log_length_cm = log(length_cm),
         log_length_cm2 = (log(length_cm))^2,
         log_circum_cm = log(circum_cm),
         log_circum_cm2 = (log(circum_cm)^2) )

# run the different models for the different species
lm.x <- 
  lapply(c("fu_sp", "fu_ve", "as_no", "fu_se"), function(x) {
  
  lm.x <- lm.allo(data = filter(allo, binomial_code == x), slope = "log_dry_weight_g", e.vars = exp.vars)
  lm.x <- lm.x[,c("model", "term", "r.squared", "AIC", "nobs")]
  return(lm.x)
  
} )

# check the best models
lapply(lm.x, function(x) print(x))

# run the best models for each species

# F. spiralis
fu_sp <- lm(log_dry_weight_g ~ log_length_cm*log_circum_cm, 
            data = filter(allo, binomial_code == "fu_sp"))
summary(fu_sp)

# check the model fit
df1 <- filter(allo, binomial_code == "fu_sp")
df1$dry_weight_pred <- exp(predict(fu_sp, data = df1))
plot(df1$dry_weight_g, df1$dry_weight_pred)
abline(a = 0, b = 1, col = "red")

# F. vesiculosus
fu_ve <- lm(log_dry_weight_g ~ log_length_cm*log_circum_cm, 
            data = filter(allo, binomial_code == "fu_ve"))
summary(fu_ve)

# check the model fit
df2 <- filter(allo, binomial_code == "fu_ve")
df2$dry_weight_pred <- exp(predict(fu_ve, data = df2))
plot(df2$dry_weight_g, df2$dry_weight_pred)
abline(a = 0, b = 1, col = "red")

# A. nodosum
as_no <- lm(log_dry_weight_g ~ log_length_cm*log_circum_cm, 
            data = filter(allo, binomial_code == "as_no"))
summary(as_no)

# check the model fit
df3 <- filter(allo, binomial_code == "as_no")
df3$dry_weight_pred <- exp(predict(as_no, data = df3))
plot(df3$dry_weight_g, df3$dry_weight_pred)
abline(a = 0, b = 1, col = "red")

# F. serratus
fu_se <- lm(log_dry_weight_g ~ log_length_cm*log_circum_cm, 
            data = filter(allo, binomial_code == "fu_se"))
summary(fu_se)

# check the model fit
df4 <- filter(allo, binomial_code == "fu_se")
df4$dry_weight_pred <- exp(predict(fu_se, data = df4))
plot(df4$dry_weight_g, df4$dry_weight_pred)
abline(a = 0, b = 1, col = "red")

# plot these model fits as a supplementary figure
df.fit <- do.call("rbind", list(df1, df2, df3, df4))
df.fit <- df.fit[,c("binomial_code", "dry_weight_g", "dry_weight_pred") ]
df.fit$species <- factor(df.fit$binomial_code, levels = c("fu_sp", "fu_ve", "as_no", "fu_se"))
levels(df.fit$species) <- c("F. spiralis", "F. vesiculosus", "A. nodosum", "F. serratus")

p1 <- 
  ggplot(data = df.fit,
       mapping = aes(x = dry_weight_g, y = dry_weight_pred, colour = species)) +
  geom_point(size = 2.5, shape = 16, alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, colour = "black", linetype = "dashed") +
  scale_colour_viridis_d(option = "A", end = 0.9) +
  facet_wrap(~species, scales = "free") +
  ylab("Predicted dry weight (g)") +
  xlab("Observed dry weight (g)") +
  guides(colour = guide_legend(override.aes = list(size = 3.5))) +
  theme_meta() +
  theme(legend.title = element_blank(),
        legend.key = element_blank())
plot(p1)

ggsave(filename = here("figures/figS1.png"), p1, dpi = 400,
       units = "cm", width = 18, height = 14)

# calculate the average and standard deviation of the error
df.fit %>%
  mutate(error = (abs(dry_weight_g-dry_weight_pred)/dry_weight_g)*100 ) %>%
  group_by(binomial_code) %>%
  summarise(error_med = median(error),
            error_m = mean(error),
            error_sd = sd(error))


lm.coef <- 
  lapply(list(fu_sp, fu_ve, as_no, fu_se), function(x) {
  
  y <- coef(x)
  names(y) <- NULL
  df <- data.frame(int = y[1],
                   s1 = y[2],
                   s2 = y[3],
                   s3 = y[4])
  return(df)
  
} )

# bind these coefficients into a data.frame
lm.coef <- bind_rows(lm.coef)
lm.coef$binomial_code <- c("fu_sp", "fu_ve", "as_no", "fu_se")

# bind these coefficients to the transect data
tra_a <- full_join(tra_a, lm.coef, by = "binomial_code")

# use these models to get the dry mass of each individual sampled
tra_a <- 
  tra_a %>%
  mutate(dry_mass_g = exp(int + (s1*log(length_cm)) + (s2*log(circum_cm)) + (s3*log(length_cm)*log(circum_cm)) )) %>%
  select(-int, -s1, -s2, -s3) %>%
  filter(!is.na(dry_mass_g))

# check distribution of biomass per transect
tra_a %>%
  group_by(transect_id, depth_zone, binomial_code) %>%
  summarise(dry_mass_g = sum(dry_mass_g, na.rm = TRUE), .groups = "drop") %>%
  ggplot(data = .,
         mapping = aes(x = depth_zone, y = dry_mass_g, fill = binomial_code)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  facet_wrap(~transect_id)

# sum up biomass per zone per species
tra_a <- 
  tra_a %>%
  group_by(depth_zone, binomial_code) %>%
  summarise(dry_mass_g = sum(dry_mass_g, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(id_cols = "depth_zone",
              names_from = "binomial_code",
              values_from = "dry_mass_g",
              values_fill = 0)

# convert the depth zone to a character
tra_a <- 
  tra_a %>%
  mutate(depth_zone = as.character(depth_zone)) %>%
  select(depth_zone, fu_se, as_no, fu_ve, fu_sp)

# write this into a .rdsfile
saveRDS(object = tra_a, file = "data/transect_ssdb.csv")

# plot a figure of the depth distribution
tra_a_plot <- 
  tra_a %>%
  rename(`F. serratus` = "fu_se",
         `A. nodosum` = "as_no",
         `F. vesiculosus` = "fu_ve",
         `F. spiralis` = "fu_sp") %>%
  pivot_longer(cols = c(`F. serratus`, `A. nodosum`, `F. vesiculosus`, `F. spiralis`),
               names_to = "species", 
               values_to = "biomass")

# change levels of the factor
tra_a_plot$species <- factor(tra_a_plot$species,
                             levels = c("F. serratus", "A. nodosum", "F. vesiculosus", "F. spiralis"))

p2 <- 
  ggplot(data = tra_a_plot,
       mapping = aes(x = depth_zone, y = biomass, fill = species)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8,
           alpha = 0.75, colour = "black") +
  scale_fill_viridis_d(option = "A", end = 0.9) +
  ylab("Standing dry biomass (g)") +
  xlab("Depth zone") +
  labs(fill = "Species") +
  theme_meta()
plot(p2)

saveRDS(object = p2, file = "figures/Fig_1a.rds")

### END
