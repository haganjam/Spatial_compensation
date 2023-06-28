#' @title: calculate standing biomass of macroalgae in different depth zones
#' 
#' @description: this script uses allometric equations to convert transects
#' of four fucoid macroalgae species to standing stock biomass in four
#' different depth zones.
#' 

# load relevant libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

# load the plotting theme
source("code/helper-plotting-theme.R")
source("code/helper-miscellaneous.R")

# import the raw transect data
tra_dat <- read_csv(file = "data/transect_data.csv",
                    col_types = c("ccccccccnnncnncc"),
                    na = c("NA"))

# check the names of the transect data
names(tra_dat)

# 1. subset the data needed to interpolate the depth to the different points
# 2. correct the depth by the water level from the nearby station
dep_dat <- 
  tra_dat %>% 
  mutate(depth_correct = (water_level_cm + depth_cm) ) %>%
  select(date, transect_id, position, depth_correct) %>%
  distinct()

# problems with the start of transect 2, remove positions 0 to 4 because we do not have a starting depth
dep_dat <- 
  dep_dat %>%
  filter( !(transect_id == 2 & position %in% 0:4) )

# split into a list
dep_list <- split(dep_dat, dep_dat$transect_id)

# loop over all transects

# intialise an output list to store interpolated depth data
dep_int <- vector("list", length = length(dep_list))
for(i in 1:length(dep_list)) {
  
  # initialise a data.frame to work with
  df <- dep_list[[i]]
  
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
  int <- 
    
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
  
  dep_int[[i]] <- int
  
} 

# bind the list into a data.frame
dep_int <- bind_rows(dep_int)

# join this back to the full dataset and reorder the columns
tra_dat <- full_join(tra_dat, dep_int, by = c("date", "transect_id", "position")) 

# reorder and rename the columns
tra_dat <- 
  tra_dat %>%
  select(site_code, date, transect_id, time, position,
         depth_interpolated, binomial_code, length_cm, circum_cm) %>%
  rename(depth_RH2000_cm = depth_interpolated)

# check the summary statistics
summary(tra_dat)

# check for unique values, especially for the binomial codes
lapply(tra_dat, function(x) unique(x))

# remove missing binomial codes and the missing values "", NA
tra_dat <- 
  tra_dat %>%
  filter( !(binomial_code %in% c("-9999", "") | 
              is.na(binomial_code) | 
              is.na(depth_RH2000_cm) )  )

# check the names of the columns
names(tra_dat)

# check the summary
summary(tra_dat)

# depth range of each transect
min_max_depths <- 
  tra_dat %>%
  group_by(transect_id) %>% 
  summarise(min_depth_RH2000_cm = min(depth_RH2000_cm, na.rm = TRUE),
            max_depth_RH2000_cm = max(depth_RH2000_cm, na.rm = TRUE))
print(min_max_depths)

for(i in 1:nrow(min_max_depths)) {
  
  x <- 
    tra_dat %>%
    filter(transect_id == as.character(i)) %>%
    filter(depth_RH2000_cm == min_max_depths[i,]$min_depth_RH2000_cm |
           depth_RH2000_cm == min_max_depths[i,]$max_depth_RH2000_cm)
  print(x)
  
}

  

# check the depth distribution
hist(tra_dat$depth_RH2000_cm)

# check the minimum and maximum algal lengths
summary(tra_dat$length_cm)


# add a depth zone variable

# check the depths of the cut variable
cut(tra_dat$depth_RH2000_cm, breaks = 4)

# midpoints of these depth zones are: -44, -32, -20, -8
tra_dat$depth_zone <- cut(tra_dat$depth_RH2000_cm, breaks = 4, labels = 1:4)


# load the allometry data
allo_dat <- read_csv(file = "data/allometry_data.csv")
head(allo_dat)

# check the sites from which these individuals were drawn
tab_s1 <- 
  allo_dat %>%
  group_by(site_code, binomial_code) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(id_cols = "site_code",
              names_from = "binomial_code",
              values_from = "n",
              values_fill = 0)

# get the site lat-lon coordinates
sites <- read_csv("data/site_code_coordinates.csv")
head(sites)

# get the relevant sites
sites <- 
  sites %>%
  filter(site_code %in% tab_s1$site_code)

# bind the latitude-longitude information to the number of individuals data
tab_s1 <- full_join(sites, tab_s1, by = "site_code")

# reorder the columns
tab_s1 <- 
  tab_s1 %>%
  select(site_code, lat_dd, lon_dd, fu_sp, fu_ve, as_no, fu_se)

# check whether the number of individuals is correct
colSums(tab_s1[,c( 4, 5, 6, 7)])

# rename the columns
names(tab_s1) <- c("Site code", "Lat (DD)", "Lon (DD)", 
                   "F. spiralis", "F. vesiculosus", "A. nodosum", "F. serratus")

# export this table as a .csv file
write_csv(x = tab_s1, file = "figures-tables/table_S1.csv")

# set up the list of explanatory variables for the different models
exp_vars <- list(c("log_length_cm", "log_circum_cm", "log_length_cm:log_circum_cm"),
                 c("log_length_cm", "log_circum_cm", "log_circum_cm2"),
                 c("log_length_cm", "log_length_cm", "log_length_cm2"),
                 c("log_length_cm","log_circum_cm"),
                 c("log_length_cm", "log_length_cm2"),
                 c("log_circum_cm", "log_circum_cm2"),
                 c("log_length_cm"),
                 c("log_circum_cm2"),
                 "1")

# add model variables to the data
allo_dat <- 
  allo_dat %>%
  mutate(log_dry_weight_g = log(dry_weight_g),
         log_length_cm = log(length_cm),
         log_length_cm2 = (log(length_cm))^2,
         log_circum_cm = log(circum_cm),
         log_circum_cm2 = (log(circum_cm)^2) )

# run the different models for the different species
lm_sp <- 
  
  lapply(c("fu_sp", "fu_ve", "as_no", "fu_se"), function(x) {
  
    # prepare the input data
    input <- 
      allo_dat %>%
      filter(binomial_code == x)
    
    # fit the different linear models using lm_compare
    lm_sp <- lm_compare(data = input, resp = "log_dry_weight_g", pred = exp_vars)
    lm_sp <- lm_sp[, c("model", "term", "r.squared", "AIC", "nobs")]
    
    # return the linear model table
    return(lm_sp)
  
} )

# generate a table for these models
names(lm_sp) <- c("fu_sp", "fu_ve", "as_no", "fu_se")
tab_s2 <- bind_rows(lm_sp, .id = "species")

# reorder the factors
tab_s2$species <- factor(tab_s2$species, levels = c("fu_sp", "fu_ve", "as_no", "fu_se"))
levels(tab_s2$species) <- c("F. spiralis", "F. vesiculosus", "A. nodosum", "F. serratus")

# generate a summary table
tab_s2 <- 
  tab_s2 %>%
  mutate(term = ifelse(term == "(Intercept)", "(Int.)", term)) %>%
  group_by(species, model) %>%
  summarise(terms = paste(term, collapse = " + "), 
            N = first(nobs),
            r2 = round(first(r.squared), 2), 
            AIC = round(first(AIC), 1), .groups = "drop") %>%
  select(-model) %>%
  arrange(species, AIC)

# calculate the delta AIC values
tab_s2 <- 
  tab_s2 %>%
  group_by(species) %>%
  mutate(dAIC = AIC - first(AIC)) %>%
  ungroup()

# export this table as a .csv file
write_csv(x = tab_s2, file = "figures-tables/table_S2.csv")

# run the best models for each species

# which models are best for each species?
tab_s2 %>%
  group_by(species) %>%
  filter(AIC == min(AIC)) %>%
  View()

# F. spiralis

# subset the data
df1 <- 
  allo_dat %>%
  filter(binomial_code == "fu_sp")

# fit the model
fu_sp <- lm(log_dry_weight_g ~ log_length_cm + log_circum_cm + log_circum_cm2, 
            data = df1)

# graphical analyses of residuals
# plot(fu_sp)

# check the model fit
df1$dry_weight_pred <- exp(predict(fu_sp, data = df1)) * df1_dsf

# calculate Duan's smearing factor
df1_dsf <- mean(exp(residuals(fu_sp)))

plot(df1$dry_weight_g, df1$dry_weight_pred )
abline(a = 0, b = 1, col = "red")

# F. vesiculosus

# subset the data
df2 <- 
  allo_dat %>%
  filter(binomial_code == "fu_ve")

# fit the model
fu_ve <- lm(log_dry_weight_g ~ log_length_cm + log_circum_cm + log_circum_cm2, 
            data = df2)

# graphical analyses of residuals
# plot(fu_ve)

# calculate Duan's smearing factor
df2_dsf <- mean(exp(residuals(fu_ve)))

# check the model fit
df2$dry_weight_pred <- exp(predict(fu_ve, data = df2)) * df2_dsf
plot(df2$dry_weight_g, df2$dry_weight_pred)
abline(a = 0, b = 1, col = "red")

# A. nodosum

# subset the data
df3 <- 
  allo_dat %>%
  filter(binomial_code == "as_no")

# fit the model
as_no <- lm(log_dry_weight_g ~ log_length_cm + log_circum_cm, 
            data = df3)

# graphical analyses of residuals
# plot(as_no)

# calculate Duan's smearing factor
df3_dsf <- mean(exp(residuals(as_no)))

# check the model fit
df3$dry_weight_pred <- exp(predict(as_no, data = df3))* df3_dsf
plot(df3$dry_weight_g, df3$dry_weight_pred )
abline(a = 0, b = 1, col = "red")

# F. serratus

# subset the data
df4 <- 
  allo_dat %>%
  filter(binomial_code == "fu_se")

# fit the model
fu_se <- lm(log_dry_weight_g ~ log_length_cm + log_circum_cm, 
            data = df4)

# graphical analyses of residuals
# plot(fu_se)

# calculate Duan's smearing factor
df4_dsf <- mean(exp(residuals(fu_se)))

# check the model fit
df4$dry_weight_pred <- exp(predict(fu_se, data = df4)) * df4_dsf 
plot(df4$dry_weight_g, df4$dry_weight_pred)
abline(a = 0, b = 1, col = "red")

# make a table with the model coefficients: Table S3
mod_list <- c("fu_sp", "fu_ve", "as_no", "fu_se") 
dsf <- list(df1_dsf, df2_dsf, df3_dsf, df4_dsf)

ts3_list <- vector("list", length = length(mod_list))
for(i in 1:length(mod_list)) {
  
  x <- eval(parse(text = mod_list[i]))
  y <- broom::tidy(x)
  y[["species"]] <- mod_list[i]
  y[["duan_smearing_factor"]] <- dsf[[i]]
  z <- dplyr::select(y, species, term, estimate, std.error, statistic, p.value, duan_smearing_factor)
  ts3_list[[i]] <- z
  
}
tab_s3 <- bind_rows(ts3_list)

# export this table as a .csv file
write_csv(x = tab_s3, file = "figures-tables/table_S3.csv")

# plot these model fits as a supplementary figure
df_fit <- do.call("rbind", list(df1, df2, df3, df4))
df_fit <- df_fit[,c("binomial_code", "dry_weight_g", "dry_weight_pred") ]
df_fit$species <- factor(df_fit$binomial_code, levels = c("fu_sp", "fu_ve", "as_no", "fu_se"))
levels(df_fit$species) <- c("F. spiralis", "F. vesiculosus", "A. nodosum", "F. serratus")

p1 <- 
  ggplot(data = df_fit,
       mapping = aes(x = dry_weight_g, y = dry_weight_pred, colour = species)) +
  geom_point(size = 2.5, shape = 16, alpha = 0.75) +
  geom_abline(intercept = 0, slope = 1, colour = "black", linetype = "dashed") +
  # scale_colour_viridis_d(option = "A", end = 0.9, direction = 1) +
  scale_colour_manual(values = seaweed_pal()) +
  facet_wrap(~species, scales = "free") +
  ylab("Predicted dry weight (g)") +
  xlab("Observed dry weight (g)") +
  guides(colour = guide_legend(override.aes = list(size = 3.5))) +
  theme_meta() +
  theme(legend.position = "none",
        strip.background = element_blank())
plot(p1)

ggsave(filename = "figures-tables/fig_S1.svg", p1,
       units = "cm", width = 15, height = 14)

# calculate the average and standard deviation of the error
df_fit %>%
  mutate(error = (abs(dry_weight_g-dry_weight_pred)/dry_weight_g)*100 ) %>%
  group_by(binomial_code) %>%
  summarise(n = n(),
            error_med = median(error),
            error_m = mean(error),
            error_sd = sd(error))


# use these models to predict the dry weight for the individuals in the transect

# first, get the relevant variables
tra_dat <- 
  tra_dat %>%
  mutate(log_length_cm = log(length_cm),
         log_length_cm2 = (log(length_cm))^2,
         log_circum_cm = log(circum_cm),
         log_circum_cm2 = (log(circum_cm)^2))

# split by the binomial code
tra_dat <- split(tra_dat, tra_dat$binomial_code)
lm_list <- list(as_no, fu_se, fu_sp, fu_ve)
dfs_list <- list(df1_dsf, df2_dsf, df3_dsf, df4_dsf)

# predict the dry weights
tra_dat <- 
  
  mapply(function(x, y, z) {
  
  x[["dry_mass_g"]] <- exp(predict(object = y, x)) * z
  return(x)
  
}, tra_dat, lm_list, dfs_list, SIMPLIFY = FALSE)

# bind back into a data.frame
tra_dat <- bind_rows(tra_dat)

# 1. remove log-transformed columns 
# 2. remove individuals for which predictions were unavailable
tra_dat <- 
  tra_dat %>%
  select(-starts_with("log")) %>%
  filter(!is.na(dry_mass_g))

# check distribution of biomass per transect
tra_dat %>%
  group_by(transect_id, depth_zone, binomial_code) %>%
  summarise(dry_mass_g = sum(dry_mass_g, na.rm = TRUE), .groups = "drop") %>%
  ggplot(data = .,
         mapping = aes(x = depth_zone, y = dry_mass_g, fill = binomial_code)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8) +
  facet_wrap(~transect_id)

# sum up biomass per zone per species
tra_dat <- 
  tra_dat %>%
  group_by(depth_zone, binomial_code) %>%
  summarise(dry_mass_g = sum(dry_mass_g, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(id_cols = "depth_zone",
              names_from = "binomial_code",
              values_from = "dry_mass_g",
              values_fill = 0)

# convert the depth zone to a character
tra_dat <- 
  tra_dat %>%
  mutate(depth_zone = as.character(depth_zone)) %>%
  select(depth_zone, fu_se, as_no, fu_ve, fu_sp) %>%
  rename(depth = depth_zone)

# write this into a .rdsfile
saveRDS(object = tra_dat, file = "output/transect_ssdb.rds")

# plot a figure of the depth distribution
tra_plot <- 
  tra_dat %>%
  rename(`F. serratus` = "fu_se",
         `A. nodosum` = "as_no",
         `F. vesiculosus` = "fu_ve",
         `F. spiralis` = "fu_sp") %>%
  pivot_longer(cols = c(`F. serratus`, `A. nodosum`, `F. vesiculosus`, `F. spiralis`),
               names_to = "species", 
               values_to = "biomass")

# change levels of the factor
tra_plot$species <- factor(tra_plot$species,
                           levels = c("F. spiralis", "F. vesiculosus", "A. nodosum", "F. serratus"))

# modify the depth zone factor
tra_plot$depth <- factor(tra_plot$depth)
levels(tra_plot$depth) <- paste0("DZ", c("4", "3", "2", "1"))
tra_plot$depth <- factor(tra_plot$depth, levels = paste0("DZ", c("1", "2", "3", "4")) )

p2 <- 
  ggplot(data = tra_plot,
       mapping = aes(x = depth, y = biomass, fill = species)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8,
           alpha = 0.75, colour = "black", size = 0.1) +
  # scale_fill_viridis_d(option = "A", end = 0.9) +
  scale_fill_manual(values = seaweed_pal()) +
  ylab("Standing dry biomass (g)") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 575)) +
  xlab("") +
  labs(fill = "Species") +
  theme_meta() +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), 
             linetype = "longdash", colour = "black", 
             size = 0.25) +
  theme(axis.text.x = element_text(size = 12),
        axis.ticks.x = element_blank())
plot(p2)

saveRDS(object = p2, file = "output/fig_1a.rds")

### END
