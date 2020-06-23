
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
library(nlme)

# load the overdispersion function (GLMM FAQ)

overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}



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

# check the seedlings all variable
euc_ana$seedlings_all %>%
  range()

euc_ana %>%
  filter(seedlings_all > 80) %>%
  View()

# surveyID = 44 has a very high number of small seedlings


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
  c(site_vars[!(site_vars %in% c("Easting", "Northing"))], "Property_season", "MrVBF")

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

# how many plots in total?
nrow(euc_ana)

# how many plots in total have any seedlings?
euc_ana %>%
  filter(seedling_y_n > 0) %>%
  nrow()

# 25% of plots have any recruitment

euc_ana %>%
  filter(Euc_canopy_cover > 0) %>%
  nrow()

ggplot(data = euc_ana %>%
         filter(seedling_y_n > 0),
       mapping = aes(x = (distance_to_Eucalypt_canopy_m) )) +
  geom_histogram()

euc_ana %>%
  filter(seedling_y_n > 0, distance_to_Eucalypt_canopy_m < 10) %>%
  nrow()

# 72 of 89 (i.e. 81%) are within 10 m of a Eucalypt canopy

euc_ana %>%
  filter(seedling_y_n > 0, distance_to_Eucalypt_canopy_m < 20) %>%
  nrow()

# 80 of 89 (i.e. 90%) are within 20 m of a Eucalypt canopy

# how are grass cover and distance to Eucalypt canopy cover related?
ggplot(data = euc_ana,
       mapping = aes(x = (distance_to_Eucalypt_canopy_m), y = total_perennial_grass,
                     colour = Property)) +
  geom_point()

# what is the relationship between grass and seedlings within a 20 m eucalypt radius?
ggplot(data = euc_ana %>%
         filter(distance_to_Eucalypt_canopy_m < 20),
       mapping = aes(x = log10(total_perennial_grass), y = seedling_y_n,
                     colour = Season)) +
  geom_point() +
  facet_wrap(~Property, scales = "free")

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


# get site-year combinations where seedlings were observed
recruit_sites


### do a paired analysis at the site-season scale

g_vars <- c("total_perennial_grass",
            "ExoticAnnualGrass_cover",
            "native_perennial_grass")

pair_dat <- 
  euc_ana %>%
  filter(Property_season %in% recruit_sites) %>%
  group_by(Property_season, seedling_y_n) %>%
  summarise_at(vars(g_vars),
               list(m = ~ mean(., na.rm = TRUE),
                    n = ~ n()) ) %>%
  ungroup() %>%
  group_by(Property_season) %>%
  mutate(n = n()) %>%
  filter(n > 1) %>%
  summarise_at(vars(paste(g_vars, c("m"), sep = "_")),
               ~ diff(., na.rm = TRUE)) %>%
  ungroup() %>%
  gather(paste(g_vars, c("m"), sep = "_"),
         key = "grass_variable", value = "cover")

euc_ana$Property_season %>%
  unique() %>%
  length()

pair_dat$Property_season %>%
  unique() %>%
  length()

# plot out these differences

ggplot(data = pair_dat, 
       mapping = aes(x = cover, fill = grass_variable)) +
  geom_histogram() +
  geom_vline(xintercept = 0) +
  theme_classic()

split(pair_dat, pair_dat$grass_variable) %>%
  lapply(., function(x) { t.test(x = x$cover, alternative = c("two.sided"), mu = 0) })


### fit a mean-level model (i.e. each property-season combination)

fit_vars <- 
  c("distance_to_Eucalypt_canopy_m",
    "PET",
    "bare_all",
    "non_wood_plant",
    "total_perennial_grass",
    "ExoticAnnualGrass_cover",
    "Litter_cover",
    seed_vars)

mean_dat <- 
  euc_ana %>%
  group_by(Property, Season) %>%
  summarise_at(vars(fit_vars),
               ~ mean(., na.rm = TRUE)) %>%
  ungroup()

mean_dat %>%
  select(fit_vars[1:round((length(fit_vars)/2), 0)]) %>%
  pairs()

mean_dat %>%
  select(fit_vars[1:round((length(fit_vars)/2), 0)]) %>%
  cor(method = "spearman") %>%
  corrplot(method = "number")

mean_dat %>%
  select(fit_vars[(round((length(fit_vars)/2), 0)+1):length(fit_vars)]) %>%
  pairs()

mean_dat %>%
  select(fit_vars[(round((length(fit_vars)/2), 0)+1):length(fit_vars)]) %>%
  cor(method = "spearman") %>%
  corrplot(method = "number")

# check the variable distributions

mean_dat %>%
  gather(fit_vars[1:10],
         key = "variable", value = "value") %>%
  ggplot(data = .,
         mapping = aes(x = value)) +
  geom_histogram() +
  facet_wrap(~ variable, scales = "free") +
  theme_classic()

# distance to Eucalypt canopy is a bit skewed
# Exotic annual grass cover is also a bit skewed
# total perennial grass is also a bit skewed

# plot some of these relationships
ggplot(data = mean_dat,
       mapping = aes(x = PET, y = log(1+seedlings_all), colour = Season)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  theme_classic()

ggplot(data = mean_dat,
       mapping = aes(x = sqrt(distance_to_Eucalypt_canopy_m), y = log(1+seedlings_all), 
                     colour = PET)) +
  geom_point() +
  geom_smooth(se = FALSE, method = "lm") +
  facet_wrap(~ Season) +
  theme_classic()

names(mean_dat)

# fit a model without random effects using restricted maximum likelihood (function gls)
lm_mean1 <- nlme::gls(log(1+seedlings_all) ~ 
                   sqrt(distance_to_Eucalypt_canopy_m) + 
                   sqrt(ExoticAnnualGrass_cover) +
                   sqrt(total_perennial_grass) +
                   Litter_cover +
                   bare_all +
                   PET +
                   Season,
                   data = mean_dat, method = "REML")
AIC(lm_mean1)

# fit a model with the same fixed effects and a random effect using restricted maximum likelihood (function lmer)
lm_mean2 <- lme4::lmer(log(1+seedlings_all) ~ 
                    sqrt(distance_to_Eucalypt_canopy_m) + 
                    sqrt(ExoticAnnualGrass_cover) +
                    sqrt(total_perennial_grass) +
                    Litter_cover +
                    bare_all +
                    PET +
                    Season +
                    (1|Property),
                  data = mean_dat, REML = TRUE)

AIC(lm_mean2)

# we do not accept the model with the random effect

# check the model assumptions
plot(lm_mean1)
residuals(lm_mean1) %>% 
  hist()

car::vif(lm_mean1)

rsquared(lm_mean1)

# check model coefficients
summary(lm_mean1)

# plot the predictions

plot(log(1+mean_dat$seedlings_all), predict(lm_mean1))

ggplot(data = mean_dat %>%
         mutate(pred = predict(lm_mean1)),
       mapping = aes(x = sqrt(distance_to_Eucalypt_canopy_m), 
                     y = log(1+mean_dat$seedlings_all))) +
  geom_point() +
  geom_smooth(mapping = aes(y = pred), method = "lm") +
  theme_classic()



### fit a model to the full dataset

# negative binomial is not a bad fit but first should fit a normal Poisson

# however, there is a lot of unexplained variance when canopy cover is very low
# I want to see which variables explain this variation before modelling further

fit_vars2 <- 
  c("distance_to_Eucalypt_canopy_m",
    "PET",
    "bare_all",
    "non_wood_plant",
    "total_perennial_grass",
    "ExoticAnnualGrass_cover",
    "Litter_cover", 
    "Euc_canopy_cover",
    "K_perc",
    "Th_ppm",
    "MrVBF")

scale_this <- function(x){
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm=TRUE)
}

full_dat <- 
  euc_ana %>%
  mutate_at(vars(fit_vars2),
               ~ scale_this(.))
  
lm_full1 <- 
  lme4::glmer.nb(seedlings_all ~ 
             (ExoticAnnualGrass_cover) +
             (total_perennial_grass) +
             Litter_cover +
             PET +
             Season +
             Euc_canopy_cover:(distance_to_Eucalypt_canopy_m) +
             (1|Property),
           data = full_dat %>%
             filter(SurveyID != 44), verbose = TRUE)

overdisp_fun(lm_full1)

ggplot(data = full_dat %>%
         filter(SurveyID != 44) %>%
         mutate(pred = predict(lm_full1, type = "response")),
       mapping = aes(x = distance_to_Eucalypt_canopy_m, 
                     y = seedlings_all)) +
  geom_point() +
  geom_point(mapping = aes(y = pred), colour = "red") +
  facet_wrap(~Property, scales = "free") +
  theme_classic()

piecewiseSEM::rsquared(lm_full1)


### explore dataset with low distance to eucalypt canopy

can_dat <- 
  euc_ana %>%
  filter(distance_to_Eucalypt_canopy_m == 0) %>%
  select(SurveyID, Season, Property, Property_season, quadrat_no,
         fit_vars2, seedlings_all)

can_dat <- 
  can_dat %>%
  gather(fit_vars2, 
         key = "variable", value = "value") %>%
  split(., .$variable)

names(can_dat)

ggplot(data = can_dat[[11]],
       mapping = aes(x = value, 
                     y = seedlings_all)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~Property, scales = "free") +
  theme_classic()


### model this for properties where there was some recruitment

fit_vars2 <- 
  c("distance_to_Eucalypt_canopy_m",
    "PET",
    "bare_all",
    "non_wood_plant",
    "total_perennial_grass",
    "ExoticAnnualGrass_cover",
    "Litter_cover", 
    "Euc_canopy_cover",
    "K_perc",
    "Th_ppm",
    "MrVBF")

recruits_prop <- 
  euc_ana %>%
  group_by(Property) %>%
  summarise(sum = sum(seedlings_all) ) %>%
  ungroup() %>%
  filter(sum > 0) %>%
  pull(Property)

scale_this <- function(x){
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm=TRUE)
}

rec_dat <- 
  euc_ana %>%
  filter(Property %in% recruits_prop) %>%
  filter(distance_to_Eucalypt_canopy_m == 0) %>%
  select(SurveyID, Season, Property, Property_season, quadrat_no,
         fit_vars2, seedlings_all, seedling_y_n) %>%
  mutate_at(vars(fit_vars2),
            ~ scale_this(.))

lm_rec1 <- 
  lme4::glmer.nb(seedlings_all ~ 
                   total_perennial_grass +
                   Litter_cover +
                   ExoticAnnualGrass_cover +
                   (1|Property) +
                   (total_perennial_grass+0|Property),
                 data = rec_dat, verbose = TRUE)

overdisp_fun(lm_rec1)

ggplot(data = rec_dat %>%
         mutate(pred = predict(lm_rec1, type = "response")),
       mapping = aes(x = total_perennial_grass, 
                     y = seedling_y_n)) +
  geom_point() +
  geom_point(mapping = aes(y = pred), colour = "red") +
  facet_wrap(~Property, scale = "free") +
  theme_classic()

piecewiseSEM::rsquared(lm_rec1)

