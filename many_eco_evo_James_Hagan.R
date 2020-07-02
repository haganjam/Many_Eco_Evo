
# Title: ManyEcoEvo analysis (James G. Hagan)

# load relevant libraries
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(broom)
library(viridis)
library(here)
library(corrplot)
library(piecewiseSEM)
library(quantreg)
library(ggpubr)

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
  c(chisq = Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
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

# scale and center function
scale_this <- 
  function(x){
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm=TRUE)
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

# examine the rows that were removed
rem_rows <- 
  euc_dat %>%
  filter_all( any_vars( (is.na(.)) ) ) %>%
  pull(SurveyID)

euc_dat %>%
  filter(SurveyID %in% rem_rows) %>%
  View()

# remove rows without complete data
euc_dat <- 
  euc_dat %>%
  filter_all( all_vars( (!is.na(.)) ) )

# check if the correct rows were removed
euc_dat %>%
  filter(SurveyID %in% rem_rows) %>% 
  nrow()

# check basic summary statistics
summary(euc_dat)



# create vectors with different variable combinations

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



# data exploration

# do remotely sensed sun variables represent aspect?

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


# do the remotely sensed variables represent landscape position?
ggplot(data = euc_dat,
       mapping = aes(x = landscape_position, y = MrVBF)) +
  geom_point() +
  theme_classic()


# how is Euc_canopy_cover and distance_to_Eucalypt_canopy_m correlated?
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


# check the precipitation variables
prec_vars

euc_dat %>%
  split(., .$Season) %>%
  lapply(function (x) {select(x, prec_vars) %>% pairs()} )

euc_dat %>%
  split(., .$Season) %>%
  lapply(function (x) {select(x, prec_vars) %>% cor() } )

# all precipitation variables are highly correlated


# check the soil variables
soil_vars

euc_dat %>%
  split(., .$Season) %>%
  lapply(function (x) {select(x, soil_vars) %>% pairs()} )

euc_dat %>%
  split(., .$Season) %>%
  lapply(function (x) {select(x, soil_vars) %>% cor() } )


# cover variables
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


# seedling variables
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



# analysis data

euc_ana <- 
  euc_dat %>%
  filter(SurveyID != 239) # removes data point with strange canopy cover value


# create new analysis variables

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


# reduce the euc_ana variables
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



# summary statistics

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

mean(prop_cov$prop, na.rm = TRUE)
sd(prop_cov$prop, na.rm = TRUE) 

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

euc_ana %>%
  group_by(Property) %>%
  summarise(n = unique(Property_season) %>% length())



# between property analysis

names(euc_ana)

# choose a subset of variables to use
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


# summarise data at the property scale:
# (1) calculate proportion of plots with seedlings at each property for each time point
# (2) calculate the average for each variable (fit_vars) across plots for each time point at each property
# (3) calcualte the average and standard deviation for each property across time points

mean_dat <- 
  euc_ana %>%
  group_by(Property, Season) %>%
  mutate(seedling_y_n = sum(seedling_y_n, na.rm = TRUE)/n()) %>%
  summarise_at(vars(fit_vars),
               ~ mean(., na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Property) %>%
  summarise_at(vars(all_of(fit_vars)), 
               list( ~ mean(., na.rm = TRUE),
                     ~ sd(., na.rm = TRUE) ) )


# (Results: General site characteristics)

# range of mean exotic cover
range(mean_dat$exotic_proportion_mean)

mean_dat %>%
  filter(exotic_proportion_mean == max(exotic_proportion_mean) | exotic_proportion_mean == min(exotic_proportion_mean)) %>%
  select(contains("exotic_proportion"))

# mean and standard deviation proportion of exotic annual grass cover of total non-woody exotic plant cover
mean_dat %>%
  mutate(egrass_prop = ExoticAnnualGrass_cover_mean/exotic_non_woody_mean) %>%
  summarise(m = mean(egrass_prop), sd = sd(egrass_prop))

# calculate the range in total annual precipitation
mean_dat %>%
  filter(annual_precipitation_mean == max(annual_precipitation_mean) | annual_precipitation_mean == min(annual_precipitation_mean)) %>%
  select(contains("annual_precipitation"))

# are the plant cover variables correlated?
mean_dat %>%
  select(paste(covs, c("mean"), sep = "_")) %>%
  pairs()

mean_dat %>%
  select(paste(covs, c("mean"), sep = "_")) %>%
  cor(method = "spearman")


# is litter cover correlated with total annual precipitation?
ggplot(data = mean_dat,
       mapping = aes(x = annual_precipitation_mean, y = Litter_cover_mean)) +
  geom_point()

cor.test(mean_dat$annual_precipitation_mean, mean_dat$Litter_cover_mean,
         method = "spearman")


# is total non-woody plant cover correlated with total annual precipitation?
ggplot(data = mean_dat,
       mapping = aes(x = annual_precipitation_mean, y = total_non_woody_mean)) +
  geom_point()

cor.test(mean_dat$annual_precipitation_mean, mean_dat$total_non_woody_mean,
         method = "spearman")


# is exotic plant cover correlated with total annual precipitation?
ggplot(data = mean_dat,
       mapping = aes(x = annual_precipitation_mean, y = exotic_non_woody_mean)) +
  geom_point()

cor.test(mean_dat$annual_precipitation_mean, mean_dat$exotic_non_woody_mean,
         method = "spearman")

# these increases are due to increases in both exotic grass cover and exotic herb cover
ggplot(data = mean_dat,
       mapping = aes(x = annual_precipitation_mean, y = exotic_grass_mean)) +
  geom_point()

ggplot(data = mean_dat,
       mapping = aes(x = annual_precipitation_mean, y = exotic_herbs_mean)) +
  geom_point()


# is native plant cover correlated with total annual precipitation?
ggplot(data = mean_dat,
       mapping = aes(x = annual_precipitation_mean, y = native_non_woody_mean)) +
  geom_point()

ggplot(data = mean_dat,
       mapping = aes(x = annual_precipitation_mean, y = native_perennial_grass_mean)) +
  geom_point()


# is total grass cover correlated with total annual precipitation?
ggplot(data = mean_dat,
       mapping = aes(x = annual_precipitation_mean, y = total_grass_mean)) +
  geom_point()

cor.test(mean_dat$annual_precipitation_mean, mean_dat$total_grass_mean,
         method = "spearman")


# does exotic proportion correlate with total annual precipitation?
ggplot(data = mean_dat,
       mapping = aes(x = annual_precipitation_mean, y = exotic_proportion_mean)) +
  geom_point()

cor.test(mean_dat$annual_precipitation_mean, mean_dat$exotic_proportion_mean,
         method = "spearman")


# does distance to eucalypt canopy correlate with total annual precipitation?
ggplot(data = mean_dat,
       mapping = aes(x = annual_precipitation_mean, y = distance_to_Eucalypt_canopy_m_mean)) +
  geom_point()

cor.test(mean_dat$annual_precipitation_mean, mean_dat$distance_to_Eucalypt_canopy_m_mean,
         method = "spearman")


# additional bivariate correlations for interpretation purposes (results not presented)
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


# check other potentially important explanatory variables that I did not include
mean_dat %>%
  gather("MrVBF_mean", "K_perc_mean", "Th_ppm_mean", "U_ppm_mean", "shrub_mean",
         "Euc_canopy_cover_mean",
         key = "pred", val = "value") %>%
  ggplot(data = .,
         mapping = aes(x = value, y = seedling_y_n_mean)) +
  geom_point() +
  facet_wrap(~pred, scales = "free")


# (Results: Effect of grass cover on seedling recruitment between properties)

# examine distributions of the chosen explanatory variables
mean_dat %>%
  select(exotic_proportion_mean, Litter_cover_mean, annual_precipitation_mean,
         distance_to_Eucalypt_canopy_m_mean, total_grass_mean, exotic_grass_mean,
         native_perennial_grass_mean, exotic_herbs_mean) %>%
  gather(key = "var", value = "val") %>%
  ggplot(data = .,
         mapping = aes(x = (val) )) +
  geom_histogram() +
  facet_wrap(~var, scales = "free")


# set up models to fit

# variables to conserve in all models (i.e. moderators or covariates)
cons_vars <- c("annual_precipitation_mean")

# set up a list of explanatory variables for the seven models including an intercept-only null model
exp_vars <- list(c('1'),
                 c(cons_vars),
                 c("exotic_proportion_mean", cons_vars),
                 c("total_grass_mean", cons_vars),
                 c("exotic_grass_mean", cons_vars),
                 c("native_perennial_grass_mean", cons_vars),
                 c("exotic_herbs_mean", cons_vars))


# use a loop to fit the different models using the glm() function with binomial errors

# create output lists
lm_out_glm <- vector("list", length = length(exp_vars))
cof_out_glm <- vector("list", length = length(exp_vars))
diag_out_glm <- vector("list", length = length(exp_vars))

# run a loop to fit the different models
for (i in seq_along( 1:length(exp_vars) ) ) {
  
  lm_out_glm[[i]] <- 
    glm(reformulate(exp_vars[[i]], "seedling_y_n_mean"), 
        data = mean_dat,
        family = binomial)
  
  cof_out_glm[[i]] <- tidy(lm_out_glm[[i]])
  
  diag_out_glm[[i]] <- 
    bind_cols(glance(lm_out_glm[[i]]), piecewiseSEM::rsquared(lm_out_glm[[i]]) ) %>%
    mutate(overdisp = lm_out_glm[[i]]$deviance/lm_out_glm[[i]]$df.residual,
           AICc = AICc(lm_out_glm[[i]]))
  
}

# check model assumptions by plotting residuals against each variable
lapply(lm_out_glm[2:length(lm_out_glm)], function(x) { 
  
  mean_dat %>%
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
diag_out_glm %>%
  mutate(delta_AICc = AICc - min(AICc)) %>%
  select(model, predictors, method, R.squared, AICc, delta_AICc, aicc_weights) %>%
  arrange(desc(aicc_weights) ) %>%
  write_csv(here("figures_tables/table_2.csv"))


# data exploration revealed triangular distirbutions between total_grass_cover and exotic grass cover with Eucalypt seedling proportion
ggplot(data = mean_dat,
       mapping = aes(x = total_grass_mean, y = seedling_y_n_mean)) +
  geom_point()

ggplot(data = mean_dat,
       mapping = aes(x = exotic_grass_mean, y = seedling_y_n_mean)) +
  geom_point()


# fit linear quantile regression model to the 95th percentile of Eucalypt seedling proportion and total grass cover using the rq() function
qm_1 <- rq(seedling_y_n_mean ~ total_grass_mean, tau = 0.95, data = mean_dat)
summary.rq(qm_1, se = "ker")

# output model predictions
pred_dat_qm_1 <- 
  data.frame(total_grass_mean = seq(min(mean_dat$total_grass_mean), 
                                    max(mean_dat$total_grass_mean), 0.5),
             annual_precipitation_mean = mean(mean_dat$annual_precipitation_mean))

pred_dat_qm_1 <- 
  bind_cols(pred_dat_qm_1, as_tibble(predict(qm_1, pred_dat_qm_1, interval = "confidence")) )

pred_dat_qm_1

# rename annual_precipitation_mean TAP for plotting
p_dat_1 <- 
  mean_dat %>%
  mutate(MAP = annual_precipitation_mean)

# plot the data and quantile regression model predictions
p_quant_1 <- 
  ggplot() +
  geom_ribbon(data = pred_dat_qm_1,
              mapping = aes(x = total_grass_mean, ymin = lower, ymax = higher),
              fill = "black", alpha = 0.1) +
  geom_line(data = pred_dat_qm_1,
            mapping = aes(x = total_grass_mean, y = fit), colour = "black") +
  geom_point(data = p_dat_1,
             mapping = aes(x = total_grass_mean, y = seedling_y_n_mean, colour = MAP)) +
  geom_errorbar(data = p_dat_1,
                mapping = aes(x = total_grass_mean,
                              ymin = seedling_y_n_mean - seedling_y_n_sd,
                              ymax = seedling_y_n_mean + seedling_y_n_sd,
                              colour = MAP) ) +
  scale_colour_viridis_c() +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  ylab("Eucalypt seedling proportion") +
  xlab("total grass cover (%)") +
  theme_classic() +
  theme(legend.position = "none")


# fit linear quantile regression model to the 95th percentile of Eucalypt seedling proportion and exotic grass cover using the rq() function
qm_2 <- rq(seedling_y_n_mean ~ exotic_grass_mean, tau = 0.95, 
           data = mean_dat, ci = FALSE)
summary.rq(qm_2, se = "ker")

# output model predictions
pred_dat_qm_2 <- 
  data.frame(exotic_grass_mean = seq(min(mean_dat$exotic_grass_mean), 
                                    max(mean_dat$exotic_grass_mean), 0.5),
             annual_precipitation_mean = mean(mean_dat$annual_precipitation_mean))

pred_dat_qm_2 <- 
  bind_cols(pred_dat_qm_2, 
            as_tibble(predict(qm_2, pred_dat_qm_2, interval = "confidence")) )

pred_dat_qm_2

# rename the total annual precipitation TAP for plotting
p_dat_2 <- 
  mean_dat %>%
  mutate(TAP = annual_precipitation_mean)

p_quant_2 <- 
  ggplot() +
  geom_ribbon(data = pred_dat_qm_2,
              mapping = aes(x = exotic_grass_mean, ymin = lower, ymax = higher),
              fill = "black", alpha = 0.1) +
  geom_line(data = pred_dat_qm_2,
            mapping = aes(x = exotic_grass_mean, y = fit), colour = "black") +
  geom_point(data = p_dat_2,
             mapping = aes(x = exotic_grass_mean, y = seedling_y_n_mean, colour = TAP)) +
  geom_errorbar(data = p_dat_2,
                mapping = aes(x = exotic_grass_mean,
                              ymin = seedling_y_n_mean - seedling_y_n_sd,
                              ymax = seedling_y_n_mean + seedling_y_n_sd,
                              colour = TAP) ) +
  scale_colour_viridis_c() +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  ylab("") +
  xlab("exotic grass cover (%)") +
  theme_classic()

# arrange the plots together
fig_1 <- 
  ggarrange(p_quant_1, p_quant_2, labels = c("(a)", "(b)"),
          widths = c(1, 1.3), font.label = list(face = "plain", size = 11),
          hjust = c(-0.2,-0.2))

fig_1

ggsave(here("figures_tables/fig_1.png"), fig_1, dpi = 300,
       width = 17, height = 8, units = "cm")



# (Results: Effect of grass cover on seedling recruitment within properties)

# choose variables to consider in the analysis
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

# how many properties do not have any Eucalypt seedling recruitment?
euc_ana %>%
  group_by(Property) %>%
  summarise(n = sum(seedling_y_n, na.rm = TRUE)) %>%
  filter(n == 0) %>%
  pull(Property)

# how many properties have at least two plots with recruitment?
euc_ana %>%
  group_by(Property) %>%
  summarise(n = sum(seedling_y_n, na.rm = TRUE)) %>%
  filter(n >= 2) %>%
  pull(Property) %>%
  length()

# create a vector of property-time combinations with two more plots with and without Eucalypt seedlings
well_samps <- 
  euc_ana %>%
  group_by(Property, Season, Property_season, seedling_y_n) %>%
  summarise(n = n()) %>%
  spread(key = "seedling_y_n", value = "n") %>%
  filter(`0` >= 2, `1` >= 2) %>%
  pull(Property_season)

# how many properties to these well-sampled property-time combinations come from?
euc_ana %>%
  filter(Property_season %in% well_samps) %>%
  pull(Property) %>%
  unique()

# table s1 - create a table of the property-time combinations used in the within-site analysis
euc_ana %>%
  group_by(Property, Season, Property_season, seedling_y_n) %>%
  summarise(n = n()) %>%
  spread(key = "seedling_y_n", value = "n") %>%
  filter(`0` >= 2, `1` >= 2) %>%
  ungroup() %>%
  select(-Property_season) %>%
  rename(property = Property,
         season = Season,
         'plots without Eucalypt seedling (n)' = `0`,
         'plots with Eucalypt seedling (n)' = `1`) %>%
  write_csv(., here("figures_tables/table_s1.csv"))


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


# define the variables to use in the analysis
g_vars <- c("native_perennial_grass",
            "exotic_grass",
            "total_grass",
            "total_non_woody",
            "Litter_cover",
            "exotic_herbs",
            "annual_precipitation",
            "Euc_canopy_cover", 
            "distance_to_Eucalypt_canopy_m")


# for all property-time combinations in well_samps (n = 16):
# calculate the average across plots with and without Eucalypt seedlings for each variable in g_vars
# calculate the difference between the average value with and without Eucalypt seedlings for each variable
# use gather() to make these data tidy

pair_dat <- 
  euc_ana %>%
  filter(Property_season %in% well_samps) %>%
  group_by(Property, Season, Property_season, seedling_y_n) %>%
  summarise_at(vars(g_vars),
               list(m = ~ mean(., na.rm = TRUE),
                    n = ~ n()) ) %>%
  ungroup() %>%
  group_by(Property, Season, Property_season) %>%
  summarise_at(vars(paste(g_vars, c("m"), sep = "_")),
               ~ diff(., na.rm = TRUE)) %>%
  ungroup() %>%
  gather(paste(g_vars, c("m"), sep = "_"),
         key = "grass_variable", value = "cover")


# use two-tailed, one-sample t-tests to test whether the difference in:
# total grass cover
# exotic grass cover
# native grass cover
# between plots with and without Eucalypt seedlings differs from zero

# create a vector of these grass variables
grass_vars <- c("total_grass_m",
                "native_perennial_grass_m",
                "exotic_grass_m")

t_test_grass <- 
  pair_dat %>%
  filter(grass_variable %in% grass_vars) %>%
  split(., .$grass_variable) %>%
  lapply(., function(x) { t.test(x = x$cover, alternative = c("two.sided"), mu = 0) })

t_test_grass


# fig 2

# for each variable, calculate a mean difference for plotting
pair_dat_means <- 
  pair_dat %>%
  group_by(grass_variable) %>%
  summarise(cover = mean(cover, na.rm = TRUE) ) %>%
  ungroup()

# set grass axes names
grass_axes <- c("total grass cover diff. (%)",
                "native grass cover diff. (%)",
                "exotic grass cover diff. (%)")

grass_out <- vector("list")

for(i in seq_along(1:length(grass_vars))) {
  
  grass_out[[i]] <- 
    pair_dat %>%
    filter(grass_variable == grass_vars[i]) %>%
    ggplot(data = ., mapping = aes(x = cover)) +
    geom_histogram(bins = 25, alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_vline(data = filter(pair_dat_means, grass_variable == grass_vars[i]),
               mapping = aes(xintercept = cover), colour = "red") +
    annotate("text", x = Inf, y = Inf, 
             label = paste(c("t = "), round(t_test_grass[[i]]$statistic[[1]], 2), sep = "" ), 
             hjust = 1, vjust = 2) +
    annotate("text", x = Inf, y = Inf, 
             label = paste(c("P = "), round(t_test_grass[[i]]$p.value[[1]], 2), sep = "" ), 
             hjust = 1, vjust = 4) +
    xlab(grass_axes[i]) +
    ylab("frequency") +
    theme_classic()
}

fig_2 <- ggarrange(grass_out[[1]], grass_out[[2]], grass_out[[3]],
          ncol = 3, labels = c("(a)", "(b)", "(c)"),
          font.label = list(face = "plain", size = 11),
          hjust = c(-0.2,-0.2,-0.2))

ggsave(here("figures_tables/fig_2.png"), fig_2, dpi = 300,
       width = 20, height = 9, units = "cm")



# explore relationships between potential confounding variables
pair_dat_explan <- 
  euc_ana %>%
  filter(Property_season %in% well_samps) %>%
  group_by(Property, Season, Property_season, seedling_y_n) %>%
  summarise_at(vars(g_vars),
               list(m = ~ mean(., na.rm = TRUE),
                    n = ~ n()) ) %>%
  ungroup() %>%
  group_by(Property, Season, Property_season) %>%
  summarise_at(vars(paste(g_vars, c("m"), sep = "_")),
               ~ diff(., na.rm = TRUE)) %>%
  ungroup()

pair_dat_explan %>%
  select(-Property, -Season, -Property_season) %>%
  pairs()

# test these potentially confounding variables
conf_vars <- c("Litter_cover_m",
               "distance_to_Eucalypt_canopy_m_m",
               "total_non_woody_m")

t_test_conf <- 
  pair_dat %>%
  filter(grass_variable %in% conf_vars) %>%
  split(., .$grass_variable) %>%
  lapply(., function(x) { t.test(x = x$cover, alternative = c("two.sided"), mu = 0) })

t_test_conf


# in these property-time combinations, are the plots with and without eucalypt seedlings unbalanced
euc_ana %>%
  filter(Property_season %in% well_samps) %>%
  group_by(Property, Season, Property_season, seedling_y_n) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(Property, Season, Property_season) %>%
  summarise(n_diff = diff(n, na.rm = TRUE)) %>%
  ungroup() %>%
  pull(n_diff) %>%
  t.test(x = ., alternative = c("two.sided"), mu = 0)


# (Methods and materials: Within-property analysis):

# repeat the within-property analysis but using properties as replicates and not property-time combinations as these are not strictly indepedent

pair_dat_ind <- 
  euc_ana %>%
  filter(Property_season %in% well_samps) %>%
  group_by(Property, Season, seedling_y_n) %>%
  summarise_at(vars(g_vars),
               list(m = ~ mean(., na.rm = TRUE),
                    n = ~ n()) ) %>%
  ungroup() %>%
  group_by(Property, seedling_y_n) %>%
  summarise_at(vars(paste(g_vars, c("m"), sep = "_")),
               ~ mean(., na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Property) %>%
  summarise_at(vars(paste(g_vars, c("m"), sep = "_")),
               ~ diff(., na.rm = TRUE)) %>%
  ungroup() %>%
  gather(paste(g_vars, c("m"), sep = "_"),
         key = "grass_variable", value = "cover")

# perform the one-sample t-tests
pair_dat_ind %>%
  filter(grass_variable %in% grass_vars) %>%
  split(., .$grass_variable) %>%
  lapply(., function(x) { t.test(x = x$cover, alternative = c("two.sided"), mu = 0) })

# perform one-sample t-tests on confounding variables
pair_dat_ind %>%
  filter(grass_variable %in% conf_vars) %>%
  split(., .$grass_variable) %>%
  lapply(., function(x) { t.test(x = x$cover, alternative = c("two.sided"), mu = 0) })


# are the differences in grass cover variables correlated with differences in confounding variables?

grass_cover <- c("native_perennial_grass",
                 "exotic_grass",
                 "total_grass")

diff_cor <- 
  euc_ana %>%
  filter(Property_season %in% well_samps) %>%
  group_by(Property, Season, Property_season, seedling_y_n) %>%
  summarise_at(vars(grass_cover), ~ mean(., na.rm = TRUE) ) %>%
  ungroup() %>%
  group_by(Property, Season, Property_season) %>%
  summarise_at(vars(grass_cover),
               ~ diff(., na.rm = TRUE)) %>%
  ungroup() %>%
  gather(grass_cover, 
         key = "grass_variable",
         value = "cover_difference")

diff_cor

# examine the correlations among these variables

conf <- c("total_non_woody",
          "Litter_cover",
          "distance_to_Eucalypt_canopy_m",
          "annual_precipitation")

corr_conf <- 
  euc_ana %>%
  filter(Property_season %in% well_samps) %>%
  group_by(Property, Season, Property_season) %>%
  summarise_at(vars(conf), ~ mean(., na.rm = TRUE)) %>%
  ungroup() %>%
  select(-Property, -Season) %>%
  gather(conf, 
         key = "conf_variable",
         value = "value")

# join these dataframes

diff_cor <- 
  full_join(corr_conf,
            diff_cor,
            by = "Property_season")

View(diff_cor)

# rename the grass variables for plotting

diff_cor <- 
  diff_cor %>%
  mutate(grass_variable = if_else(grass_variable == "exotic_grass",
                                  "exotic grass (%)",
                                  if_else(grass_variable == "total_grass",
                                          "total grass (%)", "native per. grass (%)"))) %>%
  rename('grass variable' = grass_variable)
  


# plot out the relationship between these confounding variables and the grass cover difference

diff_cor_plot <- 
  split(diff_cor, diff_cor$conf_variable)
  
xlabs <- c("total annual precipitation (mm)", "distance to Eucalypt canopy (m)",
           "litter cover (%)", "non-woody plant cover (%)")

ylabs <- c("cover difference (%)", "", "cover difference (%)", "")

cor_plots_out <- vector("list", length = length(diff_cor_plot))

for (i in seq_along(1:length(diff_cor_plot)) ) {
  
  cor_plots_out[[i]] <- 
    
    ggplot(data = diff_cor_plot[[i]],
         mapping = aes(x = value, y = cover_difference, colour = `grass variable`)) +
    geom_point(size = 2, shape = 16) +
    xlab(c(xlabs[i])) +
    ylab(c(ylabs[i])) +
    scale_colour_viridis_d() +
    theme_classic() +
    theme(legend.position = "none")
  
  
}

fig_3 <- 
  ggarrange(cor_plots_out[[1]],
          cor_plots_out[[2]],
          cor_plots_out[[3]],
          cor_plots_out[[4]],
          labels = c("(a)", "(b)", "(c)", "(d)"),
          common.legend = TRUE, legend = "right",
          font.label = list(face = "plain", size = 11))

ggsave(here("figures_tables/fig_3.png"), fig_3, dpi = 300,
       width = 20, height = 15, units = "cm")

spearman_out <- 
  
  lapply(diff_cor_plot, function(x) {
  
  split(x, x$`grass variable`) %>%
    lapply(., function(y) {
      cor.test(x = y$value,
               y = y$cover_difference,
               method = "spearman") %>%
        tidy() 
    }
    ) %>%
    bind_rows(., .id = "grass_variable")
  
})

# bind into a dataframe  
spearman_out <- 
  spearman_out %>%
  bind_rows(., .id = "variable")

table_3 <- 
  spearman_out %>%
  select(variable, grass_variable, estimate, statistic, p.value) %>%
  gather(estimate, statistic, p.value,
         key = "stat", value = "val") %>%
  spread(grass_variable, val) %>%
  mutate_at(vars(c("exotic grass (%)", "native per. grass (%)", "total grass (%)")),
               ~ round(., digits = 2))


write_csv(table_3, path = here("figures_tables/table_3.csv"))




# examine whether differences grass cover are explained by confounding variables
diff_cor %>%
  select(-Property, -Season, -Property_season) %>%
  pairs()

diff_cor %>%
  select(-Property, -Season, -Property_season) %>%
  cor(method = "spearman") %>%
  corrplot(method = "number", type = "lower")


# plot three correlation matrices for each grass cover variables

diff_cor %>%
  select(paste(grass_cover[1], "m", sep = "_"), conf) %>%
  cor(method = "spearman") %>%
  corrplot(method = "number", type = "lower", diag = FALSE)
  
diff_cor %>%
  select(paste(grass_cover[1], "m", sep = "_"), conf) %>%
  pairs()

diff_cor %>%
  gather(native_perennial_grass_m)








