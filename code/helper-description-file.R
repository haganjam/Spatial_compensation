
# description file code

# load the holepunch library
library(holepunch)

# set the author name etc.
options(
  usethis.full_name = "James Hagan",
  usethis.description = list(
    `Authors@R` = 'person("James", "Hagan", email = "james_hagan@outlook.com", role = c("aut", "cre"), 
    comment = c(ORCID = "0000-0002-7504-3393"))',
    License = "MIT",
    Version = "1.0.0"
  ))

# write the description file
write_compendium_description(
  package = "Spatial compensation",
  description = "Spatial compensation of biomass production in Fucoid macroalgae",
  version = "1.0")
