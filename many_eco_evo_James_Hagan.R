################################################################################################
### Many Eco Evo analysis: James G. Hagan ######################################################
################################################################################################

# Load relevant libraries
library(tidyverse)

# Load the data
euc_dat <- readr::read_delim("data/Euc_data.csv", delim = ",")

# Check the variables
euc_dat %>% names()

# Rename problematic variables names
euc_dat <- euc_dat %>%
  rename(Easting = "Easting ",
         quadrat_no = "Quadrat no",
         landscape_position = "Landscape position",
         distance_to_Eucalypt_canopy_m = "Distance_to_Eucalypt_canopy(m)",
         euc_sdlgs_50cm_2m = "euc_sdlgs50cm-2m",
         euc_sdlgs_2m = "euc_sdlgs>2m"
         )

# Perform various checks and balances on the data

# Check the structure of the different variables

str(euc_dat) # variables are structured as expected

# Check various individual variables
euc_dat$Season %>% 
  unique()

euc_dat$Property %>% 
  unique() %>%
  length()

euc_dat$Aspect %>%
  unique()

# Check the dimensions of the data
nrow(euc_dat)

ncol(euc_dat)

# Check for missing data
lapply(euc_dat, function(x) { if_else(is.na(x), 1, 0) %>% sum() } )

# Which rows have missing data for different variables
lapply(euc_dat, function(x) { 
  bind_cols(SurveyID = euc_dat$SurveyID, na_row = is.na(x) ) %>%
    filter(na_row == TRUE) %>%
    pull(SurveyID) 
  } )

# Therefore, different rows have different types of missing data


# How is Euc_canopy_cover and distance_to_Eucalypt_canopy_m correlated?

















