
# Title: ManyEcoEvo analysis (James G. Hagan)


# to do:

# check weird outlier with canopy cover and distance to eucalypt

# choose the variables

# perform the paired analysis and then run some models


# load relevant libraries
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(broom)
library(RColorBrewer)
library(viridis)
library(here)
library(corrplot)
library(lme4)
library(piecewiseSEM)

# load the data
euc_dat <- readr::read_delim("data/Euc_data.csv", delim = ",")

# check the variables
names(euc_dat)

# rename problematic variables names
euc_dat <- 
  euc_dat %>%
  rename(Easting = "Easting ",
         quadrat_no = "Quadrat no",
         landscape_position = "Landscape position",
         distance_to_Eucalypt_canopy_m = "Distance_to_Eucalypt_canopy(m)",
         euc_sdlgs_50cm_2m = "euc_sdlgs50cm-2m",
         euc_sdlgs_2m = "euc_sdlgs>2m"
         )

### perform various checks and balances on the data

# check the structure of the different variables

str(euc_dat) # variables are structured as expected

# check the dimensions of the data
nrow(euc_dat)

ncol(euc_dat)

# check for missing data
lapply(euc_dat, function(x) { if_else(is.na(x), 1, 0) %>% sum() } )

# which rows have missing data for different variables
lapply(euc_dat, function(x) { 
  bind_cols(SurveyID = euc_dat$SurveyID, na_row = is.na(x) ) %>%
    filter(na_row == TRUE) %>%
    pull(SurveyID) 
  } )

# different rows have different types of missing data

# remove rows without complete data
euc_dat <- 
  euc_dat %>%
  filter_all( all_vars( (!is.na(.)) ) )

# check basic summary statistics
summary(euc_dat)


### create vectors with different variable combinations

# site variables
site_vars <- c("SurveyID", "Date", "Season", "Property", "quadrat_no", 
               "Easting", "Northing" )

# cover variables
cov_vars <- c("ExoticAnnualGrass_cover", "ExoticAnnualHerb_cover", "ExoticPerennialHerb_cover", 
              "ExoticPerennialGrass_cover", "ExoticShrub_cover", "NativePerennialFern_cover", 
              "NativePerennialGrass_cover", "NativePerennialHerb_cover", "NativePerennialGraminoid_cover", 
              "NativeShrub_cover", "BareGround_cover", "Litter_cover", "MossLichen_cover", "Rock_cover"
              )

# eucalyptus variables
euc_vars <- c("Euc_canopy_cover", "distance_to_Eucalypt_canopy_m")

# soil variables
soil_vars <- c("K_perc", "Th_ppm", "U_ppm")

# precipitation variables
prec_vars <- c("annual_precipitation", "precipitation_warmest_quarter", 
               "precipitation_coldest_quarter", "PET")

# landscape position
land_vars <- c("landscape_position", "MrVBF")

# sun variables and aspect
sun_vars <- c("Aspect", "SRad_Jan", "SRad_Jul")

# seedling variables
seed_vars <- c("euc_sdlgs0_50cm", "euc_sdlgs_50cm_2m", "euc_sdlgs_2m")



### data exploration

### do remotely sensed sun variables represent aspect?

euc_dat %>% 
  group_by(Aspect) %>%
  summarise(n = n())

euc_dat %>% 
  select(site_vars, sun_vars) %>%
  ggplot(data = .,
         mapping = aes(x = Aspect, y = SRad_Jul, colour = SRad_Jan)) +
  geom_point() +
  facet_wrap(~Season, scales = "free")

euc_dat %>% 
  select(site_vars, sun_vars) %>%
  ggplot(data = .,
         mapping = aes(x = SRad_Jan, y = SRad_Jul)) +
  geom_point() +
  facet_wrap(~Season, scales = "free")



### do the remotely sensed variables represent landscape position?
ggplot(data = euc_dat,
       mapping = aes(x = landscape_position, y = MrVBF)) +
  geom_point() +
  theme_classic()


### how is Euc_canopy_cover and distance_to_Eucalypt_canopy_m correlated?
euc_vars

ggplot(data = euc_dat,
       mapping = aes(x = Euc_canopy_cover, y = distance_to_Eucalypt_canopy_m)) +
  geom_point() +
  theme_classic()

euc_dat %>%
  filter(Euc_canopy_cover > 0, distance_to_Eucalypt_canopy_m > 0) %>%
  View()

# strange value where canopy cover is high despite nearest Eucalyptus 38 m away?

# this doesn't seem possible so I will remove this data point

# survey 239


### check the precipitation variables
prec_vars

euc_dat %>%
  split(., .$Season) %>%
  lapply(function (x) {select(x, prec_vars) %>% pairs()} )

euc_dat %>%
  split(., .$Season) %>%
  lapply(function (x) {select(x, prec_vars) %>% cor() } )

# all precipitation variables are highly correlated


### check the soil variables
soil_vars

euc_dat %>%
  split(., .$Season) %>%
  lapply(function (x) {select(x, soil_vars) %>% pairs()} )

euc_dat %>%
  split(., .$Season) %>%
  lapply(function (x) {select(x, soil_vars) %>% cor() } )


### cover variables
cov_vars

euc_dat %>%
  select(cov_vars) %>%
  rowSums()

# check correlations among the cover variables
euc_dat %>%
  split(., .$Season) %>%
  lapply(function (x) {
    select(x, cov_vars[grepl("Herb", cov_vars)]) %>% pairs()
    } 
    )

# check the bare ground ish variables
euc_dat %>%
  split(., .$Season) %>%
  lapply(function (x) {
    select(x, BareGround_cover, Litter_cover, MossLichen_cover, Rock_cover) %>% 
      pairs()
  } 
  )

# check proportion of different cover variables
euc_dat$total_cover <- 
  euc_dat %>%
  select(cov_vars) %>%
  rowSums()

euc_dat %>%
  mutate_at(vars(cov_vars), ~(./total_cover)*100 ) %>%
  select(site_vars, cov_vars) %>%
  group_by(Property) %>%
  summarise_at(vars(cov_vars), ~mean(., na.rm = TRUE)) %>%
  gather(cov_vars, key = "key", value = "value") %>%
  ggplot(data = .,
         mapping = aes(x = value, fill = key)) +
  geom_histogram() +
  geom_vline(xintercept = 1) +
  facet_wrap(~key, scales = "free") +
  theme_classic() +
  theme(legend.position = "none")

euc_dat %>%
  mutate_at(vars(cov_vars), ~(./total_cover)*100 ) %>%
  select(cov_vars) %>%
  summary()

euc_dat %>%
  mutate_at(vars(cov_vars), ~(./total_cover)*100 ) %>%
  select(site_vars, cov_vars) %>%
  group_by(Property) %>%
  summarise_at(vars(cov_vars), ~mean(., na.rm = TRUE)) %>%
  gather(cov_vars, key = "key", value = "value") %>%
  ggplot(data = .,
         mapping = aes(y = value, x = key)) +
  geom_point() +
  geom_hline(yintercept = 1) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))


### seedling variables
seed_vars

euc_dat %>%
  split(., .$Season) %>%
  lapply(function (x) {
    select(x, seed_vars) %>% pairs()
  } 
  )

euc_dat %>%
  filter(euc_sdlgs_2m > 0) %>%
  nrow()

euc_dat %>%
  select(site_vars, seed_vars) %>%
  mutate(sdl_0_2m = (euc_sdlgs0_50cm + euc_sdlgs_50cm_2m),
         sdl_0_more_2m = (euc_sdlgs0_50cm + euc_sdlgs_50cm_2m + euc_sdlgs_2m),
         sdl_2m_more = (euc_sdlgs_50cm_2m + euc_sdlgs_2m)) %>% 
  split(., .$Season) %>%
  lapply(function (x) {
    select(x, contains("sdl")) %>% pairs()
  } 
  )


### analysis data

euc_ana <- 
  euc_dat %>%
  filter(SurveyID != 239) # removes data point with strange canopy cover value


### create new analysis variables

# (1) sum-up the three seedling size classes

euc_ana$seedlings_all <- 
  euc_ana %>%
  select(seed_vars) %>%
  rowSums()


# (2) sum-up the'bare ground' type variables:
# BareGround_cover
# MossLichen_cover

euc_ana$bare_all <- 
  euc_ana %>%
  select(BareGround_cover, MossLichen_cover) %>%
  rowSums()


# (3) create an exotic annual plant cover variable (exotic grass and herbs)
euc_ana$exotic_annual <- 
  euc_ana %>%
  select(contains("ExoticAnnual")) %>%
  rowSums()


# (4) create perennial native grass/graminoid variable
euc_ana$native_perennial_grass <- 
  euc_ana %>%
  select(contains("NativePerennialGr")) %>%
  rowSums()

# proportion of graminoids is very low
euc_ana$NativePerennialGraminoid_cover %>%
  summary()

# proportion of exotic perennial to native perennial
euc_ana$ExoticPerennialGrass_cover %>%
  summary()

(euc_ana$ExoticPerennialGrass_cover/euc_ana$native_perennial_grass) %>%
  summary()


# (5) create a total perennial grass cover variable
euc_ana$total_perennial_grass <- 
  euc_ana %>%
  select(contains("PerennialGr")) %>%
  rowSums()


# (6) create a perennial herb cover variable
euc_ana$perennial_herb <- 
  euc_ana %>%
  select(contains("PerennialHerb")) %>%
  rowSums()


# (7) create a shrub cover variable
euc_ana$shrub <- 
  euc_ana %>%
  select(contains("Shrub")) %>%
  rowSums()

# proportion of exotic shrubs is very low
euc_ana$ExoticShrub_cover %>%
  max()


# (8) create an exotic grass cover variable
euc_ana$exotic_grass <- 
  euc_ana %>%
  select(intersect(contains("Exotic"), contains("Grass"))) %>%
  rowSums()
  

# (9) create a property-season combination variables
euc_ana <- 
  euc_ana %>%
  mutate(Property_season = paste(Property, Season, sep = "_"))


# (10) create a seedling_y_n variables for all seedlings
euc_ana <- 
  euc_ana %>%
  mutate(seedling_y_n = if_else(seedlings_all > 0, 1, 0))


# (11) create a young seedling variable
euc_ana <- 
  euc_ana %>%
  mutate(young_seedling_y_n = if_else(euc_sdlgs0_50cm > 0, 1, 0))


# (12) create a total cover variable
euc_ana$total_cover <- 
  euc_ana %>%
  select(cov_vars) %>%
  rowSums()



### examine the seedlings variables because some plots do not have seedlings

# create a list of property season combinations where recruitment occurred
recruit_sites <- 
  euc_ana %>%
  group_by(Property_season) %>%
  summarise(seedling_y_n = sum(seedling_y_n),
            n = n()) %>%
  ungroup() %>%
  filter(seedling_y_n > 0) %>%
  pull(Property_season)

recruit_sites


### which cover variables are important?

cov_sub <- 
  c("bare_all", 
    "total_perennial_grass", "native_perennial_grass", "ExoticAnnualGrass_cover", "ExoticPerennialGrass_cover",
    "exotic_annual", "ExoticAnnualHerb_cover", 
    "perennial_herb", 
    "Litter_cover",
    "shrub")

euc_ana %>%
  select(cov_sub, total_cover) %>%
  mutate_at(vars(cov_sub), ~(./total_cover)*100) %>%
  summary()

# exotic perennial grass cover is generally low and therefore can use total grass cover

euc_ana %>%
  filter(shrub > 0) %>%
  View()

# keep the shrubs in mind but they are only in a few plots


### reduce the euc_ana variables
euc_ana %>%
  names()

# site variables
site_vars <- 
  c(site_vars[!(site_vars %in% c("Easting", "Northing"))], "Property_season")

# soil variables
soil_vars

# precipitation variables
prec_vars <- 
  prec_vars[prec_vars == "PET"]

# eucalyptus variables
euc_vars

# cover variables
cov_sub <- 
  c("bare_all", 
    "total_perennial_grass", "native_perennial_grass", "ExoticAnnualGrass_cover",
    "exotic_annual", "ExoticAnnualHerb_cover", 
    "perennial_herb", 
    "Litter_cover",
    "shrub",
    "total_cover")

# seedling variables
seed_vars <-
  c(seed_vars, "seedlings_all", "seedling_y_n", "young_seedling_y_n")

# reduce the variables in euc_ana
euc_ana <- 
  euc_ana %>%
  select(site_vars, soil_vars, prec_vars, euc_vars, cov_sub, seed_vars)

# create a total non-woody plant cover variable (excluding exotic annual grass and total perennial grass)
euc_ana$non_wood_plant <- 
  euc_ana %>%
  select(ExoticAnnualHerb_cover, perennial_herb) %>%
  rowSums()



### which variables are the most important?

# random effects:

# property
# season

# or just property because they were measured on different plots over time


# fixed effects:

# eucalyptus variables
# distance_to_Eucalypt_canopy_m

# precipitation variables
# PET

# cover variables
# bare_all
# non_wood_plant
# total_perennial_grass
# ExoticAnnualGrass_cover
# Litter_cover
# total_cover

((euc_ana %>%
  select(bare_all, non_wood_plant, total_perennial_grass,
         ExoticAnnualGrass_cover, Litter_cover) %>%
  rowSums()/euc_ana$total_cover)*100) %>%
  summary()

# important summary statistics to report


# explore these variables for distribution etc. within each site

vars <- 
  c("distance_to_Eucalypt_canopy_m",
    "PET",
    "bare_all",
    "non_wood_plant",
    "total_perennial_grass",
    "native_perennial_grass",
    "ExoticAnnualGrass_cover",
    "Litter_cover",
    "total_cover",
    seed_vars)

ggplot(data = euc_ana %>%
         gather(vars, key = "variable", value = "val"),
       mapping = aes(x = val)) +
  geom_histogram() +
  facet_wrap(~variable, scales = "free")


# some variables require some transformations to improve their properties
# distance_to_Eucalypt_canopy_m
# bare_all",
# non_wood_plant",
# total_perennial_grass",
# native_perennial_grass",
# ExoticAnnualGrass_cover

ggplot(data = euc_ana %>%
         mutate_at(vars(c("distance_to_Eucalypt_canopy_m", "bare_all",
                             "non_wood_plant", "total_perennial_grass",
                             "native_perennial_grass")), 
                      ~ log10(1 + .)) %>%
         gather(vars, key = "variable", value = "val"),
       mapping = aes(x = val)) +
  geom_histogram() +
  facet_wrap(~variable, scales = "free")



vars %>%
  length()

euc_ana %>%
  select(vars[c(3, 4, 5, 6, 7, 8)]) %>%
  pairs()

names(euc_ana)

ggplot(data = euc_ana, 
       mapping = aes(x = ExoticAnnualGrass_cover,
                     y = euc_sdlgs0_50cm,
                     colour = Property)) +
  geom_jitter(width = 0.05) +
  facet_wrap(~Season, scales = "free") +
  theme(legend.position = "none")


# get site-year combinations where seedlings were observed
p_sites <- 
  euc_ana %>%
  group_by(Property_season) %>%
  summarise(seedlings = sum(seedlings_all)) %>%
  ungroup() %>%
  filter(seedlings > 0) %>%
  pull(Property_season)

ggplot(data = euc_ana %>%
         filter(Property_season %in% p_sites), 
       mapping = aes(x = log(1+distance_to_Eucalypt_canopy_m),
                     y = young_seedling_y_n,
                     colour = PET)) +
  geom_jitter(width = 0.05) +
  geom_smooth(method = "lm") +
  facet_wrap(~Property_season, scales = "free") +
  theme(legend.position = "none")


### fit a mean-level model (i.e. each property-season combination)

mean_dat <- 
  euc_ana %>%
  group_by(Property) %>%
  summarise_at(vars(c("distance_to_Eucalypt_canopy_m",
                      "PET",
                      "bare_all",
                      "non_wood_plant",
                      "total_perennial_grass",
                      "ExoticAnnualGrass_cover",
                      "Litter_cover",
                      seed_vars)),
               ~ mean(., na.rm = TRUE)) %>%
  ungroup()


### use a paired test at each site-time combination

vars <- 
  c("distance_to_Eucalypt_canopy_m",
  "PET",
  "bare_all",
  "non_wood_plant",
  "total_perennial_grass",
  "native_perennial_grass",
  "ExoticAnnualGrass_cover",
  "Litter_cover",
  "total_cover",
  seed_vars[-5])

pair_dat <- 
  euc_ana %>%
  group_by(Property_season, seedling_y_n) %>%
  summarise_at(vars(vars),
               ~ mean(., na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Property_season) %>%
  mutate(n = n()) %>%
  filter(n > 1) %>%
  summarise_at(vars(vars),
               ~ diff(.) ) %>%
  ungroup()

ggplot(data = pair_dat %>%
         gather(vars, key = "variable", value = "val"),
       mapping = aes(x = val) ) +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  facet_wrap(~ variable, scales = "free") +
  theme_classic()

pair_dat %>%
  gather(vars, key = "variable", value = "val") %>%
  group_by(variable) %>%
  summarise(val_mean = mean(val))

