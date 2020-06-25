
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
library(lme4)
library(piecewiseSEM)
library(nlme)

# make a folder to export figures and tables
if(! dir.exists(here("figures_tables"))){
  dir.create(here("figures_tables"))
}


# define some useful functions

# load the overdispersion function (GLMM FAQ)
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

# AICc function from qPCR package
AICc <- function(object)
{
  aic <- AIC(object)
  if (!is.numeric(aic)) stop("Cannot calculate AIC!")
  k <- length(coef(object))
  n <- length(residuals(object))
  aic + ((2 * k * (k + 1))/(n - k - 1))
}

# akaike.weights function from qPCR package
akaike.weights <- function(x)
{
  x <- x[!is.na(x)]
  delta.aic <- x - min(x, na.rm = TRUE)
  rel.LL <- exp(-0.5 * delta.aic)
  sum.LL <- sum(rel.LL, na.rm = TRUE)
  weights.aic <- rel.LL/sum.LL
  return(list(deltaAIC = delta.aic, rel.LL = rel.LL, weights = weights.aic))
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

# sum-up the three seedling size classes

euc_ana$seedlings_all <- 
  euc_ana %>%
  select(all_of(seed_vars) ) %>%
  rowSums()

# sum-up the'bare ground' type variables:
# BareGround_cover
# MossLichen_cover

euc_ana$bare_all <- 
  euc_ana %>%
  select(BareGround_cover, MossLichen_cover) %>%
  rowSums()


# exotic plant variables

# create an exotic annual plant cover variable (exotic grass and herbs)
euc_ana$exotic_annual <- 
  euc_ana %>%
  select("ExoticAnnualGrass_cover", "ExoticAnnualHerb_cover") %>%
  rowSums()

# create an exotic grass cover variable
euc_ana$exotic_grass <- 
  euc_ana %>%
  select("ExoticAnnualGrass_cover", "ExoticPerennialGrass_cover") %>%
  rowSums()

# create an exotic herb cover variable
euc_ana$exotic_herbs <- 
  euc_ana %>%
  select("ExoticAnnualHerb_cover", "ExoticPerennialHerb_cover") %>%
  rowSums()

# create an exotic non-woody plant cover variable
euc_ana$exotic_non_woody <- 
  euc_ana %>%
  select("ExoticAnnualHerb_cover", "ExoticPerennialHerb_cover",
         "ExoticAnnualGrass_cover", "ExoticPerennialGrass_cover") %>%
  rowSums()


# native plant variables

# create perennial native grass/graminoid variable
euc_ana$native_perennial_grass <- 
  euc_ana %>%
  select("NativePerennialGrass_cover", "NativePerennialGraminoid_cover") %>%
  rowSums()

# create a total native non-woody plant variable
euc_ana$native_non_woody <- 
  euc_ana %>%
  select("NativePerennialGrass_cover", "NativePerennialGraminoid_cover",
         "NativePerennialFern_cover", "NativePerennialHerb_cover") %>%
  rowSums()


# composite plant variables

# create a total grass variable

euc_ana$total_grass <- 
  euc_ana %>%
  select("NativePerennialGrass_cover", "NativePerennialGraminoid_cover",
         "ExoticAnnualGrass_cover", "ExoticPerennialGrass_cover"
  ) %>%
  rowSums()

# create a total non-woody plant variable
euc_ana$total_non_woody <- 
  euc_ana %>%
  select("NativePerennialGrass_cover", "NativePerennialGraminoid_cover",
         "NativePerennialFern_cover", "NativePerennialHerb_cover",
         "ExoticAnnualHerb_cover", "ExoticPerennialHerb_cover",
         "ExoticAnnualGrass_cover", "ExoticPerennialGrass_cover"
         ) %>%
  rowSums()


# create a shrub cover variable
euc_ana$shrub <- 
  euc_ana %>%
  select(contains("Shrub")) %>%
  rowSums()

# proportion of exotic shrubs is very low
euc_ana$ExoticShrub_cover %>%
  max()

# create a property-season combination variables
euc_ana <- 
  euc_ana %>%
  mutate(Property_season = paste(Property, Season, sep = "_"))


# create a seedling_y_n variables for all seedlings
euc_ana <- 
  euc_ana %>%
  mutate(seedling_y_n = if_else(seedlings_all > 0, 1, 0))


# create a young seedling variable
euc_ana <- 
  euc_ana %>%
  mutate(young_seedling_y_n = if_else(euc_sdlgs0_50cm > 0, 1, 0))

# create a total plant cover variable
euc_ana$total_cover <- 
  euc_ana %>%
  select("ExoticAnnualGrass_cover", "ExoticAnnualHerb_cover",      
         "ExoticPerennialHerb_cover", "ExoticPerennialGrass_cover", 
         "ExoticShrub_cover", "NativePerennialFern_cover",
         "NativePerennialGrass_cover", "NativePerennialHerb_cover",    
         "NativePerennialGraminoid_cover","NativeShrub_cover") %>%
  rowSums()

# create a proportion exotic variable (of non woody plants)
euc_ana <- 
  euc_ana %>%
  mutate(exotic_proportion = (exotic_non_woody/total_non_woody)*100)
  

# check the seedling outlier
euc_ana %>%
  filter(seedlings_all > 80)

# surveyID = 44 has a very high number of small seedlings


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

### reduce the euc_ana variables
names(euc_ana)

euc_ana <- 
  euc_ana %>%
  select("SurveyID", "Date", "Season", "Property", "quadrat_no",
         "Property_season", 
         "MrVBF",
         "K_perc", "Th_ppm", "U_ppm",
         "annual_precipitation",
         "bare_all", "Rock_cover",
         "native_perennial_grass", "NativePerennialHerb_cover", "native_non_woody",
         "Litter_cover",
         "ExoticAnnualGrass_cover", "ExoticAnnualHerb_cover", "exotic_annual",
         "ExoticPerennialHerb_cover", "ExoticPerennialGrass_cover",
         "exotic_herbs", "exotic_grass", "exotic_non_woody",
         "total_grass", "total_non_woody",
         "Litter_cover",
         "shrub",
         "total_cover", 
         "exotic_proportion",
         "Euc_canopy_cover", "distance_to_Eucalypt_canopy_m",
         "euc_sdlgs0_50cm", "euc_sdlgs_50cm_2m", "euc_sdlgs_2m",
         "seedlings_all", "seedling_y_n", "young_seedling_y_n")


### summary statistics

# how much do different cover variables account for?

covs <- c("exotic_grass",
          "native_perennial_grass",
          "exotic_herbs")

prop_cov <- 
  euc_ana %>%
  gather(all_of(covs), key = "var", value = "cov") %>%
  group_by(SurveyID) %>%
  summarise(total_cover = mean(total_cover, na.rm = TRUE),
            cov = sum(cov, na.rm = TRUE) ) %>%
  mutate(prop = (cov/total_cover)*100)

summary(prop_cov)
hist(prop_cov$prop)

prop_cov %>% 
  filter(prop < 80) %>%
  View()

# important summary statistics to report

# how many plots in total?
nrow(euc_ana)

# how many plots in total have any seedlings?
euc_ana %>%
  filter(seedling_y_n > 0) %>%
  nrow()

# 25% of plots have any recruitment


# check the Sharrock property because it has very high recruitment
euc_ana %>%
  gather("MrVBF",
         "K_perc", "Th_ppm", "U_ppm",
         "annual_precipitation",
         "bare_all", "Rock_cover",
         "Euc_canopy_cover", "distance_to_Eucalypt_canopy_m",
         "Litter_cover", "total_cover",
         key = "variable", value = "val") %>%
  ggplot(data = .,
         mapping = aes(x = Property, y = val, colour = Season)) +
  geom_boxplot() +
  facet_wrap(~variable, scales = "free") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90))

euc_ana %>%
  group_by(Property) %>%
  summarise(n = n(),
            distance_to_Eucalypt_canopy_m = mean(distance_to_Eucalypt_canopy_m),
            total_cover = mean(total_cover),
            Litter_cover = mean(Litter_cover))



### Analysis 1 (# redo with the new variables!)

### Paired analysis: recruitment vs. no recruitment at the site scale (# needs redoing)

# get site-year combinations where seedlings were observed
recruit_sites

g_vars <- c("total_perennial_grass",
            "ExoticAnnualGrass_cover",
            "native_perennial_grass")

pair_dat <- 
  euc_ana %>%
  filter(Property_season %in% recruit_sites) %>%
  group_by(Property, Season, Property_season, seedling_y_n) %>%
  summarise_at(vars(g_vars),
               list(m = ~ mean(., na.rm = TRUE),
                    n = ~ n()) ) %>%
  ungroup() %>%
  group_by(Property, Season, Property_season) %>%
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
  facet_wrap(~Season, scales = "free") +
  theme_classic()

split(pair_dat, pair_dat$grass_variable) %>%
  lapply(., function(x) { t.test(x = x$cover, alternative = c("two.sided"), mu = 0) })


# get moderator variables

mod_vars <- 
  c("distance_to_Eucalypt_canopy_m",
    "PET",
    "bare_all",
    "non_wood_plant",
    "MrVBF",
    "K_perc",
    "Th_ppm",
    "U_ppm",
    "shrub",
    "Litter_cover")

mod_dat <- 
  euc_ana %>%
  group_by(Property, Season, Property_season) %>%
  summarise_at(vars(mod_vars),
               ~ mean(., na.rm = TRUE)) %>%
  ungroup()


# join the mod_dat data to the paired data
pair_dat_mod <- 
  full_join(pair_dat, mod_dat, by = c("Property", "Season", "Property_season")) %>%
  split(., .$grass_variable)

pair_dat_mod %>% names()

pair_dat_mod[[3]] %>%
  gather(mod_vars,
         key = "mod", value = "val") %>%
  ggplot(data = .,
         mapping = aes(x = val, y = cover)) +
  geom_point() +
  geom_hline(yintercept = 0) +
  geom_smooth(se = FALSE, method = "lm") +
  facet_wrap(~ mod, scales = "free") +
  theme_classic()




### examine mean-level trends at the property scale (i.e. each Property)

names(euc_ana)

fit_vars <- 
  c("MrVBF",
    "K_perc", "Th_ppm", "U_ppm",
    "annual_precipitation",
    "bare_all", "Rock_cover",
    "native_perennial_grass", "NativePerennialHerb_cover", "native_non_woody",
    "Litter_cover",
    "ExoticAnnualGrass_cover", "ExoticAnnualHerb_cover", "exotic_annual",
    "ExoticPerennialHerb_cover", "ExoticPerennialGrass_cover",
    "exotic_herbs", "exotic_grass", "exotic_non_woody",
    "total_grass", "total_non_woody",
    "Litter_cover",
    "shrub",
    "total_cover", 
    "exotic_proportion",
    "Euc_canopy_cover", "distance_to_Eucalypt_canopy_m",
    "euc_sdlgs0_50cm", "euc_sdlgs_50cm_2m", "euc_sdlgs_2m",
    "seedlings_all", "seedling_y_n", "young_seedling_y_n")

# three cover variables account for the vast majority of plant cover
covs

# in addition, bare_all and Litter_cover are probably also important

# seedling variables: 
# - mean seedlings across plots
# - proportion of plots with seedlings

mean_dat <- 
  euc_ana %>%
  group_by(Property, Season) %>%
  mutate(seedling_y_n = sum(seedling_y_n, na.rm = TRUE)/n()) %>%
  summarise_at(vars(fit_vars),
               ~ mean(., na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Property) %>%
  summarise_at(vars(all_of(fit_vars)), list( ~ mean(., na.rm = TRUE),
                       ~ sd(., na.rm = TRUE) ) )


# check the range of mean exotic cover
range(mean_dat$exotic_proportion_mean)

mean_dat %>%
  filter(exotic_proportion_mean == max(exotic_proportion_mean) | exotic_proportion_mean == min(exotic_proportion_mean)) %>%
  select(contains("exotic_proportion"))

mean_dat %>%
  mutate(egrass_prop = ExoticAnnualGrass_cover_mean/exotic_non_woody_mean) %>%
  summarise(m = mean(egrass_prop), sd = sd(egrass_prop))

mean_dat %>%
  filter(annual_precipitation_mean == max(annual_precipitation_mean) | annual_precipitation_mean == min(annual_precipitation_mean)) %>%
  select(contains("annual_precipitation"))



# characterise the gradient:

# litter cover
ggplot(data = mean_dat,
       mapping = aes(x = annual_precipitation_mean, y = Litter_cover_mean)) +
  geom_point()

cor.test(mean_dat$annual_precipitation_mean, mean_dat$Litter_cover_mean,
         method = "spearman")

# clear increase in litter cover along the rainfall gradient

# total non-woody plant cover
ggplot(data = mean_dat,
       mapping = aes(x = annual_precipitation_mean, y = total_non_woody_mean)) +
  geom_point()

cor.test(mean_dat$annual_precipitation_mean, mean_dat$total_non_woody_mean,
         method = "spearman")


# exotic plant cover
ggplot(data = mean_dat,
       mapping = aes(x = annual_precipitation_mean, y = exotic_non_woody_mean)) +
  geom_point()

cor.test(mean_dat$annual_precipitation_mean, mean_dat$exotic_non_woody_mean,
         method = "spearman")

# non-woody exotic plant cover increases with precipitation

# these increases are due to increases in both exotic grass cover and exotic herb cover
ggplot(data = mean_dat,
       mapping = aes(x = annual_precipitation_mean, y = exotic_grass_mean)) +
  geom_point()

ggplot(data = mean_dat,
       mapping = aes(x = annual_precipitation_mean, y = exotic_herbs_mean)) +
  geom_point()

# native plant cover
ggplot(data = mean_dat,
       mapping = aes(x = annual_precipitation_mean, y = native_non_woody_mean)) +
  geom_point()

ggplot(data = mean_dat,
       mapping = aes(x = annual_precipitation_mean, y = native_perennial_grass_mean)) +
  geom_point()

# native non-woody plant cover also increases markedly along the rainfall gradient

# total grass cover
ggplot(data = mean_dat,
       mapping = aes(x = annual_precipitation_mean, y = total_grass_mean)) +
  geom_point()

cor.test(mean_dat$annual_precipitation_mean, mean_dat$total_grass_mean,
         method = "spearman")

# do exotics dominate at wetter or drier sites?
ggplot(data = mean_dat,
       mapping = aes(x = annual_precipitation_mean, y = exotic_proportion_mean)) +
  geom_point()

cor.test(mean_dat$annual_precipitation_mean, mean_dat$exotic_proportion_mean,
         method = "spearman")

# sites vary in how dominated they are by exotics but this isn't linked to rainfall

# how do these cover variables relate to seedling recruitment
ggplot(data = mean_dat,
       mapping = aes(x = exotic_proportion_mean, y = seedling_y_n_mean)) +
  geom_point()

cor.test(mean_dat$exotic_proportion_mean, mean_dat$seedling_y_n_mean,
         method = "spearman")

ggplot(data = mean_dat,
       mapping = aes(x = total_grass_mean, y = seedling_y_n_mean)) +
  geom_point()

ggplot(data = mean_dat,
       mapping = aes(x = exotic_grass_mean, y = seedling_y_n_mean)) +
  geom_point()

ggplot(data = mean_dat,
       mapping = aes(x = native_perennial_grass_mean, y = seedling_y_n_mean)) +
  geom_point()

ggplot(data = mean_dat,
       mapping = aes(x = annual_precipitation_mean, y = seedling_y_n_mean)) +
  geom_point()

ggplot(data = mean_dat,
       mapping = aes(x = distance_to_Eucalypt_canopy_m_mean, y = seedling_y_n_mean)) +
  geom_point()

ggplot(data = mean_dat,
       mapping = aes(x = annual_precipitation_mean, y = distance_to_Eucalypt_canopy_m_mean)) +
  geom_point()

cor.test(mean_dat$annual_precipitation_mean, mean_dat$distance_to_Eucalypt_canopy_m_mean,
         method = "spearman")


# which variables are correlated with exotic proportion?
ggplot(data = mean_dat,
       mapping = aes(x = exotic_proportion_mean, y = total_grass_mean)) +
  geom_point()

ggplot(data = mean_dat,
       mapping = aes(x = exotic_proportion_mean, y = native_perennial_grass_mean)) +
  geom_point()

ggplot(data = mean_dat,
       mapping = aes(x = exotic_proportion_mean, y = exotic_grass_mean)) +
  geom_point()

ggplot(data = mean_dat,
       mapping = aes(x = exotic_proportion_mean, y = ExoticAnnualGrass_cover_mean)) +
  geom_point()

# check if I'm missing any important variables by correlating exotic proportion with other explanatory variables
mean_dat %>%
  gather("MrVBF_mean", "K_perc_mean", "Th_ppm_mean", "U_ppm_mean", "shrub_mean",
         "Euc_canopy_cover_mean",
         key = "pred", val = "value") %>%
  ggplot(data = .,
         mapping = aes(x = value, y = seedling_y_n_mean)) +
  geom_point() +
  facet_wrap(~pred, scales = "free")

# other variables don't seem very important


### fit a generalised linear model to these data

# how are the other three variables distributed?
mean_dat %>%
  select(exotic_proportion_mean, Litter_cover_mean, annual_precipitation_mean,
         distance_to_Eucalypt_canopy_m_mean, total_grass_mean, exotic_grass_mean,
         native_perennial_grass_mean, exotic_herbs_mean) %>%
  gather(key = "var", value = "val") %>%
  ggplot(data = .,
         mapping = aes(x = (val) )) +
  geom_histogram() +
  facet_wrap(~var, scales = "free")

# square root transform exotic_grass_mean and exotic_herbs because it is a bit skewed


# set up models

# variables to conserve in all models (i.e. moderators or covariates)
cons_vars <- c("annual_precipitation_mean")

# set up a list of explanatory variables for the four models
exp_vars <- list(c('1'),
                 c(cons_vars),
                 c("exotic_proportion_mean", cons_vars),
                 c("total_grass_mean", cons_vars),
                 c("exotic_grass_mean", cons_vars),
                 c("native_perennial_grass_mean", cons_vars),
                 c("exotic_herbs_mean", cons_vars))

# set up a data frame to model with the transformed variables
mod_dat <- 
  mean_dat %>%
  mutate(exotic_grass_mean = sqrt(exotic_grass_mean),
         exotic_herbs_mean = sqrt(exotic_herbs_mean))


# binomial errors for prpportion data

# create output lists
lm_out_glm <- vector("list", length = length(exp_vars))
cof_out_glm <- vector("list", length = length(exp_vars))
diag_out_glm <- vector("list", length = length(exp_vars))

# run a loop to fit the different models
for (i in seq_along( 1:length(exp_vars) ) ) {
  
  lm_out_glm[[i]] <- 
    glm(reformulate(exp_vars[[i]], "seedling_y_n_mean"), 
        data = mod_dat,
        family = binomial)
  
  cof_out_glm[[i]] <- tidy(lm_out_glm[[i]])
  
  diag_out_glm[[i]] <- 
    bind_cols(glance(lm_out_glm[[i]]), piecewiseSEM::rsquared(lm_out_glm[[i]]) ) %>%
    mutate(overdisp = lm_out_glm[[i]]$deviance/lm_out_glm[[i]]$df.residual,
           AICc = AICc(lm_out_glm[[i]]))
  
}

# check model assumptions by plotting residuals against each variable
lapply(lm_out_glm[2:length(lm_out_glm)], function(x) { 
  
  mod_dat %>%
    mutate(resids = residuals(x, type = "deviance"),
           preds = predict(x, type = "response")) %>%
    gather(unique(unlist(exp_vars[2:length(lm_out_glm)])), preds,
           key = "explan", value = "var") %>%
    ggplot(data = .,
           mapping = aes(x = var, y = resids)) +
    geom_point() +
    geom_smooth(se = FALSE) +
    facet_wrap(~explan, scales = "free") +
    theme_classic()
  
  })

# check model diagnostics and overdispersion
diag_out_glm <- 
  diag_out_glm %>%
  bind_rows(.id = "model")

# calculate the aicc weights
diag_out_glm$aicc_weights <- akaike.weights(diag_out_glm$AICc)$weights

# add the explanatory variables
diag_out_glm$predictors <- 
  lapply(exp_vars, function(x) { 
  paste(x, collapse = ".")  }) %>%
  unlist()

View(diag_out_glm)

# create a table to export
names(diag_out_glm)

diag_out_glm %>%
  mutate(delta_AICc = AICc - min(AICc)) %>%
  select(model, predictors, method, R.squared, AICc, delta_AICc, aicc_weights) %>%
  arrange(desc(aicc_weights) ) %>%
  write_csv(here("figures_tables/table_1.csv"))

cof_out_glm %>%
  bind_rows(.id = "model")

drop1(lm_out_glm[[1]], test = c("Chisq"))

confint(lm_out_glm[[1]])


# plot the predictions
ggplot(data = mean_dat %>% filter(Property != "Sharrock") %>%
         mutate(pred = predict(lm_out_glm[[3]], type = "response")),
       mapping = aes(x = exotic_proportion_mean, 
                     y = seedling_y_n_mean)) +
  geom_point() +
  geom_point(mapping = aes(y = pred), colour = "red") +
  theme_classic()


mean_dat$Property %>% unique()





# paired analysis for well-sampled sites

fit_vars2 <- 
  c("MrVBF",
    "K_perc", "Th_ppm", "U_ppm",
    "annual_precipitation",
    "bare_all", "Rock_cover",
    "native_perennial_grass", "NativePerennialHerb_cover", "native_non_woody",
    "Litter_cover",
    "ExoticAnnualGrass_cover", "ExoticAnnualHerb_cover", "exotic_annual",
    "ExoticPerennialHerb_cover", "ExoticPerennialGrass_cover",
    "exotic_herbs", "exotic_grass", "exotic_non_woody",
    "total_grass", "total_non_woody",
    "Litter_cover",
    "shrub",
    "total_cover", 
    "exotic_proportion",
    "Euc_canopy_cover", "distance_to_Eucalypt_canopy_m",
    "euc_sdlgs0_50cm", "euc_sdlgs_50cm_2m", "euc_sdlgs_2m",
    "seedlings_all", "seedling_y_n", "young_seedling_y_n")

# there are many sites without any recruitment
# then, for sites where there is some recruitment, there is like one plot
# let's only use 'well-sampled' sites

well_samps <- 
  euc_ana %>%
  group_by(Property, Season, Property_season, seedling_y_n) %>%
  summarise(n = n()) %>%
  spread(key = "seedling_y_n", value = "n") %>%
  filter(`0` >= 2, `1` >= 2) %>%
  pull(Property_season)

# explore these data
euc_ana %>%
  filter(Property_season %in% well_samps) %>%
  gather("exotic_herbs", "exotic_grass", "Litter_cover", 
         "distance_to_Eucalypt_canopy_m",
         "native_perennial_grass", "total_non_woody",
         key = "variable", value = "val") %>%
  split(., .$variable) %>%
  lapply(., function(x) { 
    
    ggplot(data = x,
           mapping = aes(x = seedling_y_n, y = val, colour = Season)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) +
      facet_wrap(~Property, scales = "free") +
      theme_classic() +
      theme(legend.position = "none")
    
    } )

# get site-year combinations where seedlings were observed
recruit_sites

g_vars <- c("total_perennial_grass",
            "ExoticAnnualGrass_cover",
            "native_perennial_grass")

pair_dat <- 
  euc_ana %>%
  filter(Property_season %in% recruit_sites) %>%
  group_by(Property, Season, Property_season, seedling_y_n) %>%
  summarise_at(vars(g_vars),
               list(m = ~ mean(., na.rm = TRUE),
                    n = ~ n()) ) %>%
  ungroup() %>%
  group_by(Property, Season, Property_season) %>%
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
  facet_wrap(~Season, scales = "free") +
  theme_classic()

split(pair_dat, pair_dat$grass_variable) %>%
  lapply(., function(x) { t.test(x = x$cover, alternative = c("two.sided"), mu = 0) })









# random slopes
# - native perennial grass
# - litter cover
# - exotic herbs
# - exotic grass
  
scale_this <- 
  function(x){
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm=TRUE)
}

full_dat <- 
  euc_ana %>%
  filter(Property_season %in% well_samps) %>%
  mutate_at(vars(fit_vars2[!(fit_vars2 %in% c("seedling_y_n", "young_seedling_y_n"))]),
               ~ scale_this(.)) %>%
  ungroup()

  



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
                     y = seedlings_all,
                     colour = Property)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  #facet_wrap(~Property, scales = "free") +
  theme_classic() +
  theme(legend.position = "none")


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

recruits_prop %>%
  length()

euc_ana$Property %>%
  unique() %>%
  length()

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

