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
str(euc_dat)

nrow()

ncol()

unique()

# Check for NAs


















