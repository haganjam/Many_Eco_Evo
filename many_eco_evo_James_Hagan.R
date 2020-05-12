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

# Remove rows without complete data
euc_dat <- euc_dat %>%
  filter_all( all_vars( (!is.na(.)) ) )

# Check basic summary statistics
summary(euc_dat)

# Site variables
site_vars <- c("SurveyID", "Date", "Season", "Property", "quadrat_no", 
               "Easting", "Northing" )

# Cover variables
cov_vars <- c("ExoticAnnualGrass_cover", "ExoticAnnualHerb_cover", "ExoticPerennialHerb_cover", 
              "ExoticPerennialGrass_cover", "ExoticShrub_cover", "NativePerennialFern_cover", 
              "NativePerennialGrass_cover", "NativePerennialHerb_cover", "NativePerennialGraminoid_cover", 
              "NativeShrub_cover", "BareGround_cover", "Litter_cover", "MossLichen_cover", "Rock_cover"
              )

# Eucalyptus variables
euc_vars <- c("Euc_canopy_cover", "distance_to_Eucalypt_canopy_m")

# Soil variables
soil_vars <- c("annual_precipitation", "precipitation_warmest_quarter", 
              "precipitation_coldest_quarter", "PET", "MrVBF", "K_perc", "Th_ppm", 
              "U_ppm", "SRad_Jan", "SRad_Jul")

# Precipitation variables
prec_vars <- c("annual_precipitation", "precipitation_warmest_quarter", 
               "precipitation_coldest_quarter", "PET")

# Landscape position
land_vars <- c("landscape_position", "MrVBF")

# Sun variables and aspect
sun_vars <- c("Aspect", "SRad_Jan", "SRad_Jul")


# Do the remotely sensed sun variables represent aspect?
ggplot(data = euc_dat %>%
         filter(landscape_position == "slope"),
       mapping = aes(x = SRad_Jul, y = SRad_Jan, colour = Aspect)) +
  geom_point()

sun_clust <- dist(select(euc_dat, sun_vars, -Aspect) %>% 
       scale(scale = TRUE, center = TRUE),
     method = "euclidean")

sun_clust <- hclust(sun_clust)

plot(sun_clust)

# Do the remotely sensed landscape variables represent landscape position?
ggplot(data = euc_dat,
       mapping = aes(x = landscape_position, y = MrVBF)) +
  geom_point() +
  theme_classic()








# How is Euc_canopy_cover and distance_to_Eucalypt_canopy_m correlated?

















