
# Title: ManyEcoEvo analysis (James G. Hagan)

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

# load the data
euc_dat <- readr::read_delim("data/Euc_data.csv", delim = ",")

# check the variables
euc_dat %>% names()

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

# check various individual variables
euc_dat$Season %>% 
  unique()

euc_dat$Property %>% 
  unique() %>%
  length()

euc_dat$Aspect %>%
  unique()

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
euc_ana <- euc_dat

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


# (4) create perennial grass/graminoid variable
euc_ana$perennial_grass <- 
  euc_ana %>%
  select(contains("PerennialGr")) %>%
  rowSums()

# very few perennial graminoids


# (5) create perennial native grass/graminoid variable
euc_ana$native_perennial_grass <- 
  euc_ana %>%
  select(contains("NativePerennialGr")) %>%
  rowSums()

# proportion of graminoids is very low
euc_ana$NativePerennialGraminoid_cover %>%
  summary()


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



### examine the seedlings variables because some plots do not have seedlings

# create a list of property season combinations were recruitment was 
recruit_sites <- 
  euc_ana %>%
  group_by(Property_season) %>%
  summarise(seedling_y_n = sum(seedling_y_n),
            n = n()) %>%
  ungroup() %>%
  filter(seedling_y_n > 0) %>%
  pull(Property_season)


### explore the grass variables

# What grass variables are there?
euc_ana %>%
  select(contains("grass")) %>%
  names()

# explore the relationship between grass cover variables
euc_ana %>%
  select(contains("grass")) %>%
  pairs()

euc_ana %>%
  select(contains("grass")) %>%
  cor() %>%
  corrplot(method = "number")

# these variables are uncorrelated and capture all grass proportion (see below)
# native_perennial_grass
# ExoticAnnualGrass_cover
# ExoticPerennialGrass_cover

# see what proportion of total grass cover they occupy:
tot_grass <- 
  euc_ana %>%
  select(cov_vars) %>%
  select(contains("gra")) %>%
  rowSums()

grass_3 <- 
  euc_ana %>%
  select(native_perennial_grass, ExoticAnnualGrass_cover, ExoticPerennialGrass_cover) %>%
  rowSums()

((grass_3/tot_grass)*100) %>%
  summary()

# NAs are cases where grass cover overall is zero


### explore the herb variables
euc_ana %>%
  select(contains("herb")) %>%
  pairs()

# generally more exotic perennial herbs
euc_ana %>%
  select(contains("PerennialHerb")) %>%
  summary()

# boundary correlations
euc_ana %>%
  select(contains("PerennialHerb")) %>%
  filter_all(all_vars(.>0)) %>%
  pairs()

euc_ana %>%
  select(contains("herb")) %>%
  names()

# use two herb variables:
# ExoticAnnualHerb_cover
# perennial_herb


### what subset of cover variables capture most important information (IMO)?

euc_ana %>%
  names()

# bare_all
# shrub
# ExoticAnnualHerb_cover
# perennial_herb
# Litter_cover
# native_perennial_grass
# ExoticAnnualGrass_cover
# ExoticPerennialGrass_cover
# perennial_grass

cov_sub <- c("bare_all", "shrub", "ExoticAnnualHerb_cover", "perennial_herb",
             "Litter_cover", "native_perennial_grass", "ExoticAnnualGrass_cover", 
             "ExoticPerennialGrass_cover", "perennial_grass")

# these variables account for almost all cover in almost all data points
((euc_ana %>%
  select(cov_sub) %>%
  rowSums() )/(euc_ana$total_cover)) %>%
  summary()

# how are these cover variables correlated?
euc_ana %>%
  select(cov_sub) %>%
  pairs()

euc_ana %>%
  select(cov_sub) %>%
  cor() %>%
  corrplot(method = "number")


# when do Eucalypt seedlings recruit?
euc_ana %>%
  group_by(Property, Season) %>%
  summarise(young_seedling_y_n = sum(young_seedling_y_n, na.rm = TRUE)/n(),
            seedling_y_n = sum(seedling_y_n, na.rm = TRUE)/n()) %>%
  View()



### explore how grass affects seedlings
euc_ana %>%
  names()

ggplot(data = euc_ana,
       mapping = aes(x = exotic_annual, y = seedlings_all, colour = Season)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~Property, scales = "free") +
  theme_classic()

euc_ana %>%
  group_by(Property) %>%
  mutate(PET = mean(PET, na.rm = TRUE)) %>%
  ggplot(data = .,
         mapping = aes(x = perennial_grass, y = seedling_y_n, colour = Season)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~as.character(PET), scales = "free") +
  theme_classic()


# within each property season combination, compare characteristics of plots with and without seedlings
euc_ana %>%
  group_by(Property) %>%
  mutate(PET = mean(PET, na.rm = TRUE)) %>%
  ungroup() %>%
  gather(cov_sub, key = "var", value = "cover") %>%
  group_by(var, PET, Season, young_seedling_y_n) %>%
  summarise(mean_cover = mean(cover, na.rm = TRUE)) %>%
  ungroup() %>%
  spread(young_seedling_y_n, mean_cover) %>%
  mutate(diff = `1`-`0`) %>%
  filter(grepl("rass", var)) %>%
  ggplot(data = .,
         mapping = aes(x = PET, y = diff)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  geom_hline(yintercept = 0) +
  facet_wrap(~var, scales = "free")


euc_ana %>%
  gather(cov_sub, euc_vars, PET, SRad_Jul, soil_vars, MrVBF,
         key = "cover_var", value = "cover") %>%
  group_by(cover_var, Property, Season, young_seedling_y_n) %>%
  summarise(mean_cover = mean(cover, na.rm = TRUE)) %>%
  ungroup() %>%
  spread(key = "young_seedling_y_n", value = "mean_cover") %>%
  mutate(cover_diff = `1`-`0`) %>%
  group_by(cover_var, Property) %>%
  summarise(mean_cover_diff = mean(cover_diff, na.rm = TRUE)) %>%
  ggplot(data = .,
         mapping = aes(x = Property, y = mean_cover_diff)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  facet_wrap(~cover_var, scales = "free") +
  theme_classic()

# which variables?
names(euc_ana)

# random effects:

# property
# season

# fixed effects:

# distance_to_Eucalypt_canopy_m
# PET
# perennial_grass
# bare_all
# ExoticAnnualGrass_cover
# ExoticAnnualHerb_cover
# perennial_herb

# how much total cover do these variables account for?






